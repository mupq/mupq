/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */


#include "bitstream.h"
#include "compat.h"
#include "io.h"
#include "kdf_shake.h"
#include "lowmc.h"
#include "mpc_lowmc.h"
#include "picnic_impl.h"
#include "randomness.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* max number of ZKB++ rounds */
#define MAX_NUM_ROUNDS 219

/* max view size per round */
#define MAX_VIEW_SIZE 75

/* max number of LowMC rounds */
#define MAX_LOWMC_R 20

typedef struct {
  uint8_t seeds[SC_PROOF][MAX_SEED_SIZE];
  uint8_t input_shares[SC_PROOF][MAX_LOWMC_BLOCK_SIZE];
  uint8_t communicated_bits[SC_PROOF][MAX_VIEW_SIZE];
  uint8_t output_shares[SC_PROOF][MAX_LOWMC_BLOCK_SIZE];
} proof_round_t;

typedef struct {
  uint8_t input_shares[SC_PROOF][MAX_LOWMC_BLOCK_SIZE];
  uint8_t communicated_bits[SC_PROOF][MAX_VIEW_SIZE];
  uint8_t output_shares[SC_PROOF][MAX_LOWMC_BLOCK_SIZE];
  const uint8_t* seeds[SC_VERIFY]; /* from signature */
  const uint8_t* commitment; /* from signature */
} verify_round_t;


static inline void clear_padding_bits(uint8_t* v, const unsigned int diff) {
  (void)v;
  (void)diff;
}

/**
 * Collapse challenge from one char per challenge to bit array.
 */
static void collapse_challenge(uint8_t* collapsed, const picnic_instance_t* pp,
                               const uint8_t* challenge) {
  bitstream_t bs;
  bs.buffer.w = collapsed;
  bs.position = 0;

  for (unsigned int i = 0; i < pp->num_rounds; ++i) {
    // flip challenge bits according to spec
    bitstream_put_bits_8(&bs, (challenge[i] >> 1) | ((challenge[i] & 1) << 1), 2);
  }
}

/**
 * Expand challenge from bit array to one char per challenge.
 */
static bool expand_challenge(uint8_t* challenge, const picnic_instance_t* pp,
                             const uint8_t* collapsed) {
  bitstream_t bs;
  bs.buffer.r = collapsed;
  bs.position = 0;

  for (unsigned int i = 0; i < pp->num_rounds; ++i) {
    const uint8_t ch = bitstream_get_bits_8(&bs, 2);
    if (ch == 3) {
      return false;
    }
    // flip challenge bits according to spec
    challenge[i] = (ch & 1) << 1 | (ch >> 1);
  }

  const size_t remaining_bits = (pp->collapsed_challenge_size << 3) - bs.position;
  if (remaining_bits && bitstream_get_bits(&bs, remaining_bits)) {
    return false;
  }

  return true;
}

static void kdf_init_from_seed(kdf_shake_t* kdf, const uint8_t* seed, const uint8_t* salt,
                               uint16_t round_number, uint16_t player_number,
                               bool include_input_size, const picnic_instance_t* pp) {
  const size_t digest_size = pp->digest_size;

  // Hash the seed with H_2.
  kdf_shake_init_prefix(kdf, digest_size, HASH_PREFIX_2);
  kdf_shake_update_key(kdf, seed, pp->seed_size);
  kdf_shake_finalize_key(kdf);

  uint8_t tmp[MAX_DIGEST_SIZE];
  kdf_shake_get_randomness(kdf, tmp, digest_size);
  kdf_shake_clear(kdf);

  // Initialize KDF with H_2(seed) || salt || round_number || player_number || output_size.
  kdf_shake_init(kdf, digest_size);
  kdf_shake_update_key(kdf, tmp, digest_size);
  kdf_shake_update_key(kdf, salt, SALT_SIZE);
  kdf_shake_update_key_uint16_le(kdf, round_number);
  kdf_shake_update_key_uint16_le(kdf, player_number);
  kdf_shake_update_key_uint16_le(kdf, pp->view_size + (include_input_size ? pp->input_size : 0));
  kdf_shake_finalize_key(kdf);
}

static void uint64_to_bitstream_10(bitstream_t* bs, const uint64_t v) {
  bitstream_put_bits(bs, v >> (64 - 30), 30);
}

static uint64_t uint64_from_bitstream_10(bitstream_t* bs) {
  return bitstream_get_bits(bs, 30) << (64 - 30);
}

static void compress_view(uint8_t* dst, const picnic_instance_t* pp, const view_t* views,
                          const unsigned int idx) {
  const size_t num_views = pp->lowmc.r;

  bitstream_t bs;
  bs.buffer.w = dst;
  bs.position = 0;

  const view_t* v = &views[0];
  if (pp->lowmc.m == 10) {
    for (size_t i = 0; i < num_views; ++i, ++v) {
      uint64_to_bitstream_10(&bs, v->t[idx]);
    }
    return;
  }
  UNREACHABLE;
}

static void decompress_view(view_t* views, const picnic_instance_t* pp, const uint8_t* src,
                            const unsigned int idx) {
  const size_t num_views = pp->lowmc.r;

  bitstream_t bs;
  bs.buffer.r = src;
  bs.position = 0;

  view_t* v = &views[0];
  if (pp->lowmc.m == 10) {
    for (size_t i = 0; i < num_views; ++i, ++v) {
      v->t[idx] = uint64_from_bitstream_10(&bs);
    }
    return;
  }
  UNREACHABLE;
}

static void decompress_random_tape(rvec_t* rvec, const picnic_instance_t* pp, const uint8_t* src,
                                   const unsigned int idx) {
  decompress_view(rvec, pp, src, idx);
}


/**
 * Compute commitment to a view.
 */
static void hash_commitment(const picnic_instance_t* pp, uint8_t* dst, const uint8_t* seed,
                            const uint8_t* input_share, const uint8_t* communicated_bits,
                            const uint8_t* output_share) {
  const size_t hashlen = pp->digest_size;

  hash_context ctx;
  // hash the seed
  hash_init_prefix(&ctx, hashlen, HASH_PREFIX_4);
  hash_update(&ctx, seed, pp->seed_size);
  hash_final(&ctx);
  uint8_t tmp[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, tmp, hashlen);
  hash_clear(&ctx);

  // compute H_0(H_4(seed), view)
  hash_init_prefix(&ctx, hashlen, HASH_PREFIX_0);
  hash_update(&ctx, tmp, hashlen);
  // hash input share
  hash_update(&ctx, input_share, pp->input_size);
  // hash communicated bits
  hash_update(&ctx, communicated_bits, pp->view_size);
  // hash output share
  hash_update(&ctx, output_share, pp->output_size);
  hash_final(&ctx);
  hash_squeeze(&ctx, dst, hashlen);
  hash_clear(&ctx);
}

static void hash_commitment_to_hash(const picnic_instance_t* pp, hash_context* dst,
                                    const uint8_t* seed, const uint8_t* input_share,
                                    const uint8_t* communicated_bits, const uint8_t* output_share) {
  const size_t hashlen = pp->digest_size;

  hash_context ctx;
  // hash the seed
  hash_init_prefix(&ctx, hashlen, HASH_PREFIX_4);
  hash_update(&ctx, seed, pp->seed_size);
  hash_final(&ctx);
  uint8_t tmp[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, tmp, hashlen);
  hash_clear(&ctx);

  // compute H_0(H_4(seed), view)
  hash_init_prefix(&ctx, hashlen, HASH_PREFIX_0);
  hash_update(&ctx, tmp, hashlen);
  // hash input share
  hash_update(&ctx, input_share, pp->input_size);
  // hash communicated bits
  hash_update(&ctx, communicated_bits, pp->view_size);
  // hash output share
  hash_update(&ctx, output_share, pp->output_size);
  hash_final(&ctx);
  hash_squeeze(&ctx, tmp, hashlen);
  hash_clear(&ctx);
  hash_update(dst, tmp, hashlen);
}

/**
 * Compute challenge from transform dependent hash - outputs {1,2 or 3}^t
 */
static void H3_compute(const picnic_instance_t* pp, hash_context* ctx, uint8_t* hash, uint8_t* ch) {
  const size_t digest_size_bits = pp->digest_size << 3;

  // Pick bits from hash
  uint8_t* eof   = ch + pp->num_rounds;
  size_t bit_idx = 0;
  while (ch < eof) {
    if (bit_idx >= digest_size_bits) {
      hash_init_prefix(ctx, pp->digest_size, HASH_PREFIX_1);
      hash_update(ctx, hash, pp->digest_size);
      hash_final(ctx);
      hash_squeeze(ctx, hash, pp->digest_size);
      hash_clear(ctx);
      bit_idx = 0;
    }

    const uint8_t twobits = (hash[bit_idx >> 3] >> ((6 - (bit_idx & 0x7)))) & 0x3;
    if (twobits != 0x3) {
      *ch++ = twobits;
    }
    bit_idx += 2;
  }
}

/**
 * Hash public key, salt and message
 */
static void H3_public_key_message(hash_context* ctx, const picnic_instance_t* pp,
                                  const uint8_t* salt, const picnic_context_t* context) {
  // hash circuit out and input (public key)
  hash_update(ctx, context->public_key, pp->output_size);
  hash_update(ctx, context->plaintext, pp->input_size);
  // hash salt
  hash_update(ctx, salt, SALT_SIZE);
  // hash message
  hash_update(ctx, context->msg, context->msglen);
}

static void H3_finalize(const picnic_instance_t* pp, hash_context* ctx, const uint8_t* salt,
                        const picnic_context_t* context, uint8_t* challenge) {
  // hash public key, salt, and message
  H3_public_key_message(ctx, pp, salt, context);
  hash_final(ctx);

  uint8_t hash[MAX_DIGEST_SIZE];
  hash_squeeze(ctx, hash, pp->digest_size);
  hash_clear(ctx);
  /* parts of this hash will be published as challenge so is public anyway */
  picnic_declassify(hash, MAX_DIGEST_SIZE);
  H3_compute(pp, ctx, hash, challenge);
}

/**
 * Re-compute challenge for verification
 */
static void H3_verify_process_round_1(const picnic_instance_t* pp, hash_context* ctx,
                                      const verify_round_t* round, uint8_t challenge) {
  const size_t output_size = pp->output_size;

  // hash output shares and commitments
  switch (challenge) {
  case 0: {
    hash_update(ctx, round->output_shares[0], output_size);
    hash_update(ctx, round->output_shares[1], output_size);
    hash_update(ctx, round->output_shares[2], output_size);
    break;
  }
  case 1: {
    hash_update(ctx, round->output_shares[2], output_size);
    hash_update(ctx, round->output_shares[0], output_size);
    hash_update(ctx, round->output_shares[1], output_size);
    break;
  }
  default: {
    hash_update(ctx, round->output_shares[1], output_size);
    hash_update(ctx, round->output_shares[2], output_size);
    hash_update(ctx, round->output_shares[0], output_size);
    break;
  }
  }
}

static void H3_verify_process_round_2(const picnic_instance_t* pp, hash_context* ctx,
                                      const verify_round_t* round, uint8_t challenge) {
  const size_t digest_size = pp->digest_size;

  // hash output shares and commitments
  switch (challenge) {
  case 0: {
    hash_commitment_to_hash(pp, ctx, round->seeds[0], round->input_shares[0],
                            round->communicated_bits[0], round->output_shares[0]);
    hash_commitment_to_hash(pp, ctx, round->seeds[1], round->input_shares[1],
                            round->communicated_bits[1], round->output_shares[1]);
    hash_update(ctx, round->commitment, digest_size);
    break;
  }
  case 1: {
    hash_update(ctx, round->commitment, digest_size);
    hash_commitment_to_hash(pp, ctx, round->seeds[0], round->input_shares[0],
                            round->communicated_bits[0], round->output_shares[0]);
    hash_commitment_to_hash(pp, ctx, round->seeds[1], round->input_shares[1],
                            round->communicated_bits[1], round->output_shares[1]);
    break;
  }
  default: {
    hash_commitment_to_hash(pp, ctx, round->seeds[1], round->input_shares[1],
                            round->communicated_bits[1], round->output_shares[1]);
    hash_update(ctx, round->commitment, digest_size);
    hash_commitment_to_hash(pp, ctx, round->seeds[0], round->input_shares[0],
                            round->communicated_bits[0], round->output_shares[0]);
    break;
  }
  }
}

/**
 * Compute challenge
 */
static void H3_process_round_1(const picnic_instance_t* pp, hash_context* ctx,
                               const proof_round_t* round) {
  // hash output shares
  // TODO process in one go
  hash_update(ctx, round->output_shares[0], pp->output_size);
  hash_update(ctx, round->output_shares[1], pp->output_size);
  hash_update(ctx, round->output_shares[2], pp->output_size);
}

static void H3_process_round_2(const picnic_instance_t* pp, hash_context* ctx,
                               const proof_round_t* round) {
  // compute and hash commitments
  hash_commitment_to_hash(pp, ctx, round->seeds[0], round->input_shares[0],
                          round->communicated_bits[0], round->output_shares[0]);
  hash_commitment_to_hash(pp, ctx, round->seeds[1], round->input_shares[1],
                          round->communicated_bits[1], round->output_shares[1]);
  hash_commitment_to_hash(pp, ctx, round->seeds[2], round->input_shares[2],
                          round->communicated_bits[2], round->output_shares[2]);
}

static uint8_t* serialize_round(const picnic_instance_t* pp, const proof_round_t* round,
                                uint8_t* tmp, uint8_t challenge) {
  // TODO: move serialization of values here to avoid work for unused values
  {
    const unsigned int c = (challenge + 2) % 3;

    // compute write commitment
    hash_commitment(pp, tmp, round->seeds[c], round->input_shares[c], round->communicated_bits[c],
                    round->output_shares[c]);
    tmp += pp->digest_size;

  }

  {
    const unsigned int b = (challenge + 1) % 3;
    // write views
    memcpy(tmp, round->communicated_bits[b], pp->view_size);
    tmp += pp->view_size;

    // write seeds
    memcpy(tmp, round->seeds[challenge], pp->seed_size);
    tmp += pp->seed_size;
    memcpy(tmp, round->seeds[b], pp->seed_size);
    tmp += pp->seed_size;
  }

  if (challenge) {
    // write input share
    memcpy(tmp, round->input_shares[SC_PROOF - 1], pp->input_size);
    tmp += pp->input_size;
  }

  return tmp;
}

static const uint8_t* deserialize_round(const picnic_instance_t* pp, verify_round_t* round,
                                        const uint8_t* data, size_t* len, uint8_t challenge) {
  const size_t digest_size            = pp->digest_size;
  const size_t seed_size              = pp->seed_size;
  const size_t input_size             = pp->input_size;
  const size_t view_size              = pp->view_size;
  const unsigned int view_diff        = pp->view_size * 8 - pp->view_round_size * pp->lowmc.r;
  const unsigned int input_share_diff = pp->input_size * 8 - pp->lowmc.k;

  const size_t base_size = pp->digest_size + pp->view_size + 2 * pp->seed_size;
  size_t requested_size  = base_size + (challenge ? pp->input_size : 0);

  if (sub_overflow_size_t(*len, requested_size, len)) {
    return NULL;
  }

  // read commitments
  round->commitment = data;
  data += digest_size;


  // read view
  memcpy(round->communicated_bits[1], data, view_size);
  if (check_padding_bits(round->communicated_bits[1][view_size - 1], view_diff)) {
    return NULL;
  }
  data += view_size;

  // read seeds
  round->seeds[0] = data;
  data += seed_size;
  round->seeds[1] = data;
  data += seed_size;

  // read input shares
  switch (challenge) {
  case 1:
    memcpy(round->input_shares[1], data, input_size);
    if (check_padding_bits(round->input_shares[1][input_size - 1], input_share_diff)) {
      return NULL;
    }
    data += input_size;
    break;
  case 2:
    memcpy(round->input_shares[0], data, input_size);
    if (check_padding_bits(round->input_shares[0][input_size - 1], input_share_diff)) {
      return NULL;
    }
    data += input_size;
    break;
  default:
    break;
  }

  return data;
}

static void generate_salt(const picnic_instance_t* pp, const picnic_context_t* context,
                          kdf_shake_t* ctx, uint8_t* salt) {
  kdf_shake_init(ctx, pp->digest_size);
  // sk || m || C || p
  kdf_shake_update_key(ctx, context->private_key, pp->input_size);
  kdf_shake_update_key(ctx, context->msg, context->msglen);
  kdf_shake_update_key(ctx, context->public_key, pp->output_size);
  kdf_shake_update_key(ctx, context->plaintext, pp->output_size);
  // N as 16 bit LE integer
  kdf_shake_update_key_uint16_le(ctx, pp->lowmc.n);
  kdf_shake_finalize_key(ctx);

  // Generate seeds but do not store them
  uint8_t buf[32];
  size_t len = pp->seed_size * SC_PROOF * pp->num_rounds;
  while (len) {
    const size_t to_process = MIN(sizeof(buf), len);
    kdf_shake_get_randomness(ctx, buf, to_process);
    len -= to_process;
  }

  // Generate salt
  kdf_shake_get_randomness(ctx, salt, SALT_SIZE);
}

static void generate_seeds(const picnic_instance_t* pp, const picnic_context_t* context,
                           kdf_shake_t* ctx) {
  kdf_shake_init(ctx, pp->digest_size);
  // sk || m || C || p
  kdf_shake_update_key(ctx, context->private_key, pp->input_size);
  kdf_shake_update_key(ctx, context->msg, context->msglen);
  kdf_shake_update_key(ctx, context->public_key, pp->output_size);
  kdf_shake_update_key(ctx, context->plaintext, pp->output_size);
  // N as 16 bit LE integer
  kdf_shake_update_key_uint16_le(ctx, pp->lowmc.n);
  kdf_shake_finalize_key(ctx);
}

int impl_sign(const picnic_instance_t* pp, const picnic_context_t* context, uint8_t* sig,
              size_t* siglen) {
  const size_t num_rounds  = pp->num_rounds;
  const unsigned int diff  = pp->input_size * 8 - pp->lowmc.n;

  const zkbpp_lowmc_implementation_f lowmc_impl       = pp->impls.zkbpp_lowmc;
  const zkbpp_share_implementation_f mzd_share        = pp->impls.mzd_share;

  // Generate salt
  uint8_t* salt = sig + pp->collapsed_challenge_size;
  kdf_shake_t seed_ctx;
  generate_salt(pp, context, &seed_ctx, salt);
  // Reset seed_ctx to produce seeds
  kdf_shake_clear(&seed_ctx);
  generate_seeds(pp, context, &seed_ctx);

  hash_context h3_ctx;
  hash_init_prefix(&h3_ctx, pp->digest_size, HASH_PREFIX_1);

  for (size_t i = 0; i < num_rounds; ++i) {
    proof_round_t round = { 0 };
    in_out_shares_t in_out_shares;
    rvec_t rvec[MAX_LOWMC_R]; // random tapes for AND-gates

    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      kdf_shake_t kdf;
      kdf_shake_get_randomness(&seed_ctx, round.seeds[j], pp->seed_size);
      kdf_init_from_seed(&kdf, round.seeds[j], salt, i, j, j != SC_PROOF - 1, pp);

      // compute sharing
      if (j < SC_PROOF - 1) {
        kdf_shake_get_randomness(&kdf, round.input_shares[j], pp->input_size);
        clear_padding_bits(&round.input_shares[j][pp->input_size - 1], diff);
        mzd_from_char_array(in_out_shares.s[j], round.input_shares[j], pp->input_size);
      } else {
        mzd_share(in_out_shares.s[2], in_out_shares.s[0], in_out_shares.s[1],
                  context->m_key);
        mzd_to_char_array(round.input_shares[SC_PROOF - 1], in_out_shares.s[SC_PROOF - 1],
                          pp->input_size);
      }

      // compute random tapes
      assert(pp->view_size <= MAX_VIEW_SIZE);
      uint8_t tape_bytes[MAX_VIEW_SIZE];
      kdf_shake_get_randomness(&kdf, tape_bytes, pp->view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);
      kdf_shake_clear(&kdf);
    }

    {
      // perform ZKB++ LowMC evaluation
      view_t views[MAX_LOWMC_R];
      lowmc_impl(context->m_plaintext, views, &in_out_shares, rvec);

      // copy output shares
      for (unsigned int j = 0; j < SC_PROOF; ++j) {
        mzd_to_char_array(round.output_shares[j], in_out_shares.s[j], pp->output_size);
      }
    }

    H3_process_round_1(pp, &h3_ctx, &round);
  }

  // reset seed_ctx to reproduce seeds
  kdf_shake_clear(&seed_ctx);
  generate_seeds(pp, context, &seed_ctx);

  for (size_t i = 0; i < num_rounds; ++i) {
    proof_round_t round = { 0 };
    in_out_shares_t in_out_shares;
    rvec_t rvec[MAX_LOWMC_R]; // random tapes for AND-gates

    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      kdf_shake_t kdf;
      kdf_shake_get_randomness(&seed_ctx, round.seeds[j], pp->seed_size);
      kdf_init_from_seed(&kdf, round.seeds[j], salt, i, j, j != SC_PROOF - 1, pp);

      // compute sharing
      if (j < SC_PROOF - 1) {
        kdf_shake_get_randomness(&kdf, round.input_shares[j], pp->input_size);
        clear_padding_bits(&round.input_shares[j][pp->input_size - 1], diff);
        mzd_from_char_array(in_out_shares.s[j], round.input_shares[j], pp->input_size);
      } else {
        mzd_share(in_out_shares.s[2], in_out_shares.s[0], in_out_shares.s[1],
                  context->m_key);
        mzd_to_char_array(round.input_shares[SC_PROOF - 1], in_out_shares.s[SC_PROOF - 1],
                          pp->input_size);
      }

      // compute random tapes
      assert(pp->view_size <= MAX_VIEW_SIZE);
      uint8_t tape_bytes[MAX_VIEW_SIZE];
      kdf_shake_get_randomness(&kdf, tape_bytes, pp->view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);
      kdf_shake_clear(&kdf);
    }

    {
      // perform ZKB++ LowMC evaluation
      view_t views[MAX_LOWMC_R];
      lowmc_impl(context->m_plaintext, views, &in_out_shares, rvec);

      // serialize view
      for (unsigned int j = 0; j < SC_PROOF; ++j) {
        mzd_to_char_array(round.output_shares[j], in_out_shares.s[j], pp->output_size);
        compress_view(round.communicated_bits[j], pp, views, j);
      }
    }

    H3_process_round_2(pp, &h3_ctx, &round);
  }
  uint8_t challenge[MAX_NUM_ROUNDS];
  H3_finalize(pp, &h3_ctx, salt, context, challenge);
  hash_clear(&h3_ctx);

  // reset seed_ctx to reproduce seeds
  kdf_shake_clear(&seed_ctx);
  generate_seeds(pp, context, &seed_ctx);

  uint8_t* tmp = sig;

  // write challenge
  collapse_challenge(tmp, pp, challenge);
  tmp += pp->collapsed_challenge_size;
  // "write salt"
  tmp += SALT_SIZE;

  for (size_t i = 0; i < num_rounds; ++i) {
    proof_round_t round = { 0 };
    in_out_shares_t in_out_shares;
    rvec_t rvec[MAX_LOWMC_R]; // random tapes for AND-gates

    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      kdf_shake_t kdf;
      kdf_shake_get_randomness(&seed_ctx, round.seeds[j], pp->seed_size);
      kdf_init_from_seed(&kdf, round.seeds[j], salt, i, j, (j != SC_PROOF - 1), pp);

      // compute sharing
      if (j < SC_PROOF - 1) {
        kdf_shake_get_randomness(&kdf, round.input_shares[j], pp->input_size);
        clear_padding_bits(&round.input_shares[j][pp->input_size - 1], diff);
        mzd_from_char_array(in_out_shares.s[j], round.input_shares[j], pp->input_size);
      } else {
        mzd_share(in_out_shares.s[2], in_out_shares.s[0], in_out_shares.s[1],
                  context->m_key);
        mzd_to_char_array(round.input_shares[SC_PROOF - 1], in_out_shares.s[SC_PROOF - 1],
                          pp->input_size);
      }

      // compute random tapes
      assert(pp->view_size <= MAX_VIEW_SIZE);
      uint8_t tape_bytes[MAX_VIEW_SIZE];
      kdf_shake_get_randomness(&kdf, tape_bytes, pp->view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);
      kdf_shake_clear(&kdf);
    }

    {
      // perform ZKB++ LowMC evaluation
      view_t views[MAX_LOWMC_R];
      lowmc_impl(context->m_plaintext, views, &in_out_shares, rvec);

      // serializes views
      for (unsigned int j = 0; j < SC_PROOF; ++j) {
        mzd_to_char_array(round.output_shares[j], in_out_shares.s[j], pp->output_size);
        compress_view(round.communicated_bits[j], pp, views, j);
      }
    }
    tmp = serialize_round(pp, &round, tmp, challenge[i]);
  }

  *siglen = tmp - sig;
  return 0;
}

int impl_verify(const picnic_instance_t* pp, const picnic_context_t* context, const uint8_t* sig,
                size_t siglen) {
  const size_t num_rounds  = pp->num_rounds;
  const size_t input_size  = pp->input_size;
  const size_t output_size = pp->output_size;
  const size_t view_size   = pp->view_size;
  const unsigned int diff  = input_size * 8 - pp->lowmc.n;

  const zkbpp_lowmc_verify_implementation_f lowmc_verify_impl = pp->impls.zkbpp_lowmc_verify;
  const zkbpp_share_implementation_f mzd_share                = pp->impls.mzd_share;

  // read and process challenge
  if (sub_overflow_size_t(siglen, pp->collapsed_challenge_size, &siglen)) {
    return -1;
  }
  uint8_t original_challenge[MAX_NUM_ROUNDS];
  if (!expand_challenge(original_challenge, pp, sig)) {
    return -1;
  }
  sig += pp->collapsed_challenge_size;

  // read salt
  if (sub_overflow_size_t(siglen, SALT_SIZE, &siglen)) {
    return -1;
  }
  const uint8_t* salt = sig;
  sig += SALT_SIZE;

  const uint8_t* old_sig = sig;
  const size_t old_siglen = siglen;

  hash_context h3_ctx;
  hash_init_prefix(&h3_ctx, pp->digest_size, HASH_PREFIX_1);

  for (size_t i = 0; i < num_rounds; ++i) {
    const unsigned int a_i = original_challenge[i];
    const unsigned int b_i = (a_i + 1) % 3;
    const unsigned int c_i = (a_i + 2) % 3;

    verify_round_t round = { 0 };
    sig = deserialize_round(pp, &round, sig, &siglen, original_challenge[i]);
    if (sig == NULL) {
      return -1;
    }

    in_out_shares_t in_out_shares;
    rvec_t rvec[MAX_LOWMC_R]; // random tapes for AND-gates
    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      kdf_shake_t kdf;
      kdf_init_from_seed(&kdf, round.seeds[j], salt, i, (j == 0) ? a_i : b_i,
                         (j == 0 && b_i) || (j == 1 && c_i), pp);

      // compute input shares if necessary
      if (j == 0 && b_i) {
        kdf_shake_get_randomness(&kdf, round.input_shares[0], input_size);
        clear_padding_bits(&round.input_shares[0][input_size - 1], diff);
      }
      if (j == 1 && c_i) {
        kdf_shake_get_randomness(&kdf, round.input_shares[1], input_size);
        clear_padding_bits(&round.input_shares[1][input_size - 1], diff);
      }

      mzd_from_char_array(in_out_shares.s[j], round.input_shares[j], input_size);

      // compute random tapes
      assert(pp->view_size <= MAX_VIEW_SIZE);
      uint8_t tape_bytes[MAX_VIEW_SIZE];
      kdf_shake_get_randomness(&kdf, tape_bytes, view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);

      kdf_shake_clear(&kdf);
    }

    {
      view_t views[MAX_LOWMC_R];
      decompress_view(views, pp, round.communicated_bits[1], 1);
      // perform ZKB++ LowMC evaluation
      lowmc_verify_impl(context->m_plaintext, views, &in_out_shares, rvec, a_i);
    }

    // recompute third output share and serialize them
    mzd_share(in_out_shares.s[2], in_out_shares.s[0], in_out_shares.s[1], context->m_key);
    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      mzd_to_char_array(round.output_shares[j], in_out_shares.s[j], output_size);
    }

    H3_verify_process_round_1(pp, &h3_ctx, &round, original_challenge[i]);
  }

  // not consumed all of the signature
  if (siglen) {
    return -1;
  }

  // reset sig and siglen to process them again
  sig = old_sig;
  siglen = old_siglen;

  for (size_t i = 0; i < num_rounds; ++i) {
    const unsigned int a_i = original_challenge[i];
    const unsigned int b_i = (a_i + 1) % 3;
    const unsigned int c_i = (a_i + 2) % 3;

    verify_round_t round = { 0 };
    sig = deserialize_round(pp, &round, sig, &siglen, original_challenge[i]);

    in_out_shares_t in_out_shares;
    rvec_t rvec[MAX_LOWMC_R]; // random tapes for AND-gates
    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      kdf_shake_t kdf;
      kdf_init_from_seed(&kdf, round.seeds[j], salt, i, (j == 0) ? a_i : b_i,
                         (j == 0 && b_i) || (j == 1 && c_i), pp);

      // compute input shares if necessary
      if (j == 0 && b_i) {
        kdf_shake_get_randomness(&kdf, round.input_shares[0], input_size);
        clear_padding_bits(&round.input_shares[0][input_size - 1], diff);
      }
      if (j == 1 && c_i) {
        kdf_shake_get_randomness(&kdf, round.input_shares[1], input_size);
        clear_padding_bits(&round.input_shares[1][input_size - 1], diff);
      }

      mzd_from_char_array(in_out_shares.s[j], round.input_shares[j], input_size);

      // compute random tapes
      assert(pp->view_size <= MAX_VIEW_SIZE);
      uint8_t tape_bytes[MAX_VIEW_SIZE];
      kdf_shake_get_randomness(&kdf, tape_bytes, view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);

      kdf_shake_clear(&kdf);
    }

    {
      view_t views[MAX_LOWMC_R];
      decompress_view(views, pp, round.communicated_bits[1], 1);
      // perform ZKB++ LowMC evaluation
      lowmc_verify_impl(context->m_plaintext, views, &in_out_shares, rvec, a_i);
      compress_view(round.communicated_bits[0], pp, views, 0);
    }

    // recompute third output share and serialize them
    mzd_share(in_out_shares.s[2], in_out_shares.s[0], in_out_shares.s[1], context->m_key);
    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      mzd_to_char_array(round.output_shares[j], in_out_shares.s[j], output_size);
    }

    H3_verify_process_round_2(pp, &h3_ctx, &round, original_challenge[i]);
  }

  assert(pp->num_rounds <= MAX_NUM_ROUNDS);
  unsigned char challenge[MAX_NUM_ROUNDS] = {0};
  H3_finalize(pp, &h3_ctx, salt, context, challenge);
  hash_clear(&h3_ctx);
  return memcmp(challenge, original_challenge, pp->num_rounds);
}
