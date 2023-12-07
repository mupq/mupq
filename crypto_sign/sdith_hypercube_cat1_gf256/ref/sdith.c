#include "sdith.h"

#include <string.h>

#include "assertions.h"
#include "gf256.h"
#include "gf2p32.h"
#include "param.h"
#include "rng.h"
#include "randombytes.h"

extern uint8_t f_poly[];
extern uint8_t leading_coefficients_of_lj_for_S[];

// If (X-alpha) divides P_in, returns P_in / (X-alpha)
void remove_one_degree_factor_from_monic(fpoly_t* out, const fpoly_t* in, fpoly_t in_degree, fpoly_t alpha) {
    out[in_degree-1] = 1;
    for(int i=in_degree-2; i>=0; i--)
        out[i] = in[i+1] ^ mul_gf256_naive(alpha, out[i+1]);
}

typedef struct isd_instance_struct {
  seed_t H_a_seed;
  ha_slice_t H_a[PAR_ha_nslice][PARAM_k] __attribute__((aligned (32))); // H_a is given by slices of 128 columns
  fsd_t y[PAR_y_size];
  fpoly_t s_A[PARAM_k];
  fpoly_t q_poly[PARAM_d][PAR_wd];
  fpoly_t p_poly[PARAM_d][PAR_wd];
} isd_instance_t;

isd_instance_t generate_inhomogeneous_sd_instance(sdith_compressed_key_t const* sk) {
  isd_instance_t reps;
  ASSERT_DRAMATICALLY(((uint64_t) (reps.H_a[0])) % 32 == 0, "bug!");
  XOF_CTX rng_ctx;
  sdith_rng_create_xof_ctx(&rng_ctx, (void*)sk->m_seed, PARAM_seed_size);
  uint8_t random_tape[8192] = {0};
  sdith_xof_next_bytes(&rng_ctx, random_tape, 8192);

#ifndef NDEBUG
  fprintf(stdout, "sk seed: %s\n", hexmem(sk->m_seed, 16));
  fprintf(stdout, "derived random_tape: %s\n", hexmem(random_tape, 64));
#endif

  uint8_t* random_tape_ptr = random_tape;

  // Generate x and s
  fsd_t x[PARAM_d][PAR_md];
  memset(x, 0, sizeof(x));
  fpoly_t s[PARAM_d][PAR_md];
  memset(s, 0, sizeof(s));

  for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
    uint8_t non_zero_pos[PAR_wd];
    fsd_t non_zero_val[PAR_wd];
    
    uint64_t i = 0;
    while (i < PAR_wd) {
      // Sample random position
      non_zero_pos[i] = *(random_tape_ptr++);
      if (non_zero_pos[i] >= PAR_md) {
        // Retry sampling if position invalid
        continue;
      }

      // Check for redundancy
      uint8_t is_redundant = 0;
      for (uint64_t j = 0; j < i; ++j) {
        if (non_zero_pos[i] == non_zero_pos[j]) {
          is_redundant = 1;
        }
      }
      if (!is_redundant) {
        i++;
      }
    }

    i = 0;
    while (i < PAR_wd) {
      // Sample random value
      non_zero_val[i] = *(random_tape_ptr++);
      if (non_zero_val[i] == 0) {
        // Retry sampling if value is zero
        continue;
      }
      i++;
    }

    // Reconstruct x
    for (uint64_t i = 0; i < PAR_md; ++i) {
      for (uint64_t j = 0; j < PAR_wd; ++j) {
        x[i_d][i] ^= non_zero_val[j] * (i == non_zero_pos[j]);
      }
    }

    // Build q poly = prod (X - j) where x_j is non zero.
    fpoly_t q[PAR_wd + 1];
    memset(q, 0, sizeof(q));
    q[0] = 1;
    for (uint64_t i = 0; i < PAR_wd; ++i) {
      for (uint64_t j = i + 1; j >= 1; --j) {
        q[j] = q[j-1] ^ mul_gf256_naive(q[j], non_zero_pos[i]);
      }
      q[0] = mul_gf256_naive(q[0], non_zero_pos[i]);
    }
    memcpy(reps.q_poly[i_d], q, sizeof(reps.q_poly[i_d]));

    // Build s poly = sum_i l_i * x_i * F(X)/(X-f_i)
    // Build p poly = sum_i l_i * x_i * Q(X)/(X-f_i)
    fpoly_t p[PAR_wd];
    memset(p, 0, sizeof(p));
    for (uint64_t i = 0; i < PAR_md; ++i) {
      fpoly_t scalar = mul_gf256_naive(leading_coefficients_of_lj_for_S[i], x[i_d][i]);

      fpoly_t temp[PAR_md];
      memset(temp, 0, sizeof(temp));
      remove_one_degree_factor_from_monic(temp, f_poly, PAR_md, i);
      for (uint64_t j = 0; j < PAR_md; j++) {
        s[i_d][j] ^= mul_gf256_naive(temp[j], scalar);
      }

      remove_one_degree_factor_from_monic(temp, q, PAR_wd, i);
      for (uint64_t j = 0; j < PAR_wd; j++) {
        p[j] ^= mul_gf256_naive(temp[j], scalar);
      }
    }
    memcpy(reps.p_poly[i_d], p, sizeof(reps.p_poly[i_d]));
  }

  // s -> s_A || s_B
  fpoly_t *s_full = (fpoly_t*) s;
  for (uint64_t i = 0; i < PARAM_k; i++) {
    reps.s_A[i] = s_full[i];
  }

  fpoly_t s_B[PAR_y_size];
  for (uint64_t i = 0; i < PAR_y_size; i++) {
    s_B[i] = s_full[PARAM_k + i];
  }

  // Sample H_a
  seed_t H_a_seed;
  memcpy(H_a_seed, random_tape_ptr, PARAM_seed_size);
  random_tape_ptr += PARAM_seed_size;
  XOF_CTX H_a_rng_ctx;
  sdith_rng_create_xof_ctx(&H_a_rng_ctx, H_a_seed, PARAM_seed_size);
  fsd_t H_a[PARAM_k][PAR_y_size];
  sdith_xof_next_bytes(&H_a_rng_ctx, H_a, PARAM_k*PAR_y_size);
  sdith_rng_free_xof_ctx(&H_a_rng_ctx);
  // copy it slice by slice to reps
  memset(reps.H_a, 0, sizeof(reps.H_a));
  for (uint64_t j=0; j<PARAM_k; ++j) {
#if (PAR_ha_nslice>1)
    for (uint64_t slice = 0; slice < PAR_ha_nslice - 1u; ++slice) {
      memcpy(reps.H_a[slice][j], H_a[j] + slice * 128, 128);
    }
#endif
    memcpy(reps.H_a[PAR_ha_nslice - 1][j], H_a[j] + (PAR_ha_nslice - 1) * 128, PAR_y_size%128);
  }

  // Compute y = s_A * H_a + s_B
  ha_slice_t y[PAR_ha_nslice];
  {
    memset(y, 0, sizeof(y));
    memcpy(y, s_B, PAR_y_size);
    for (uint64_t slice = 0; slice < PAR_ha_nslice; ++slice) {
      gf256_vec_mat128cols_muladd_ct(y[slice], reps.s_A, reps.H_a[slice], PARAM_k);
    }
    memcpy(reps.y, y, PAR_y_size);
  }

  // Validate we have a working SD instance
#ifndef NDEBUG
  fpoly_t recomputed_s[PARAM_d][PAR_md];
  fpoly_t* s_raw = (fpoly_t*) recomputed_s;
  memcpy(s_raw, reps.s_A, PARAM_k);
  memcpy(&s_raw[PARAM_k], reps.y, PAR_y_size);
  for (uint64_t k=0; k<PARAM_k; ++k) {
    for (uint64_t j = 0; j < PAR_y_size; ++j) {
      s_raw[PARAM_k + j] ^= mul_gf256_naive(H_a[k][j], reps.s_A[k]);
    }
  }

  for (uint64_t i_d=0; i_d<PARAM_d; ++i_d) {
    // Rebuild x
    fsd_t reconstructed_x[PAR_md];
    for (uint64_t k=0; k<PAR_md; ++k) {
      fsd_t c = 0;
      for (uint64_t j=PAR_md; j>=1; --j) {
        c = mul_gf256_naive(c, k) ^ recomputed_s[i_d][j-1];
      }
      reconstructed_x[k] = c;
    }
    int64_t weight = 0;
    for (uint64_t i = 0; i < PAR_md; i++) {
      if (reconstructed_x[i] != 0) weight++;
    }
    REQUIRE_DRAMATICALLY(weight == PAR_wd, "Hamming weight of solution exceeds required weight: %ld != %d", weight, PAR_wd);

    // Verify Q(i) = 0 iff x_i != 0
    for (uint64_t k=0; k<PAR_md; ++k) {
      if (reconstructed_x[k] == 0) continue;
      // Q is a monic polynomial with degree PAR_wd + 1
      fsd_t c = 1;
      for (uint64_t j=PAR_wd; j>=1; --j) {
        c = mul_gf256_naive(c, k) ^ reps.q_poly[i_d][j-1];
      }
      REQUIRE_DRAMATICALLY(c == 0, "Q poly didn't cancel out on non-zero x[%ld]: %d", k, c);
    }

    // TODO(stevenyue): verify SQ-PF=0
    fpoints_t r = 1084129906;
    fpoints_t sr = 0;
    for (uint64_t j=PAR_md; j>=1; --j) {
      sr = fpoints_mul(sr, r) ^ recomputed_s[i_d][j-1];
    }
    fpoints_t qr = 1;
    for (uint64_t j=PAR_wd; j>=1; --j) {
      qr = fpoints_mul(qr, r) ^ reps.q_poly[i_d][j-1];
    }
    fpoints_t pr = 0;
    for (uint64_t j=PAR_wd; j>=1; --j) {
      pr = fpoints_mul(pr, r) ^ reps.p_poly[i_d][j-1];
    }
    fpoints_t fr = 1;
    for (uint64_t j=PAR_md; j>=1; --j) {
      fr = fpoints_mul(fr, r) ^ f_poly[j-1];
    }
    REQUIRE_DRAMATICALLY((fpoints_mul(sr, qr) ^ fpoints_mul(pr, fr)) == 0, "SQ-PF!=0");
  }
#endif
  sdith_rng_free_xof_ctx(&rng_ctx);
  memcpy(reps.H_a_seed, H_a_seed, PARAM_seed_size);
  return reps;
}

void keygen(sdith_compressed_pubkey_t* pk, sdith_compressed_key_t* sk) {
  // Generate master seed
  randombytes(sk->m_seed, PARAM_seed_size);

  // Generate I-SD instance
  isd_instance_t instance = generate_inhomogeneous_sd_instance(sk);
  memcpy(pk->H_a_seed, instance.H_a_seed, PARAM_seed_size);
  for (uint64_t i = 0; i < PAR_y_size; i++) {
    pk->y[i] = instance.y[i];
  }
}

void uncompress_key(sdith_compressed_pubkey_t const* pk, sdith_compressed_key_t const* sk, sdith_full_pubkey_t* u_pk,
                    sdith_full_key_t* u_sk) {
  // Regenerate I-SD instance
  isd_instance_t instance = generate_inhomogeneous_sd_instance(sk);

  // copy polynomial to secret key
  memcpy(u_sk->q_poly, instance.q_poly, sizeof(instance.q_poly));
  memcpy(u_sk->p_poly, instance.p_poly, sizeof(instance.p_poly));
  
  // copy other materials
  memcpy(&u_sk->compressed_pubkey, pk, sizeof(u_sk->compressed_pubkey));
  memcpy(&u_pk->compressed_pubkey, pk, sizeof(u_pk->compressed_pubkey));

  memcpy(u_sk->s_A, instance.s_A, PARAM_k);
  memcpy(u_pk->H_a, instance.H_a, sizeof(instance.H_a));
}

void uncompress_pubkey(sdith_compressed_pubkey_t const* pk, sdith_full_pubkey_t* u_pk) {
  // Sample H_a
  seed_t H_a_seed;
  memcpy(H_a_seed, pk->H_a_seed, PARAM_seed_size);
  XOF_CTX H_a_rng_ctx;

  sdith_rng_create_xof_ctx(&H_a_rng_ctx, H_a_seed, PARAM_seed_size);
  fsd_t H_a[PARAM_k][PAR_y_size];
  sdith_xof_next_bytes(&H_a_rng_ctx, H_a, PARAM_k*PAR_y_size);
  sdith_rng_free_xof_ctx(&H_a_rng_ctx);
  // copy it slice by slice to reps
  memset(u_pk->H_a, 0, sizeof(u_pk->H_a));
  for (uint64_t j=0; j<PARAM_k; ++j) {
#if (PAR_ha_nslice>1)
    for (uint64_t slice = 0; slice < PAR_ha_nslice - 1; ++slice) {
      memcpy(u_pk->H_a[slice][j], H_a[j] + slice * 128, 128);
    }
#endif
    memcpy(u_pk->H_a[PAR_ha_nslice - 1][j], H_a[j] + (PAR_ha_nslice - 1) * 128, PAR_y_size%128);
  }
  u_pk->compressed_pubkey = *pk;
}

void serialize_pk(FILE* F, sdith_compressed_pubkey_t const* pk) {
  //TODO check for errors
  if (1!=fwrite(pk->H_a_seed, PARAM_seed_size, 1, F)) abort();
  if (1!=fwrite(pk->y, PAR_y_size, 1, F)) abort();
}

void deserialize_pk(FILE* F, sdith_compressed_pubkey_t* pk) {
  if (1!=fread(pk->H_a_seed, PARAM_seed_size, 1, F)) abort();
  if (1!=fread(pk->y, PAR_y_size, 1, F)) abort();
}

void serialize_sk(FILE* F, sdith_compressed_key_t const* sk) {
  if (1!=fwrite(sk->m_seed, PARAM_seed_size, 1, F)) abort();
}

void deserialize_sk(FILE* F, sdith_compressed_key_t* sk) {
  if (1!=fread(sk->m_seed, PARAM_seed_size, 1, F)) abort();
}

typedef struct mpc_share_struct {
  fpoly_t s_A[PARAM_k];
  fpoly_t q_poly[PARAM_d][PAR_wd];
  fpoly_t p_poly[PARAM_d][PAR_wd];
  fpoints_t a[PARAM_d][PARAM_t];
  fpoints_t b[PARAM_d][PARAM_t];
  fpoints_t c[PARAM_d][PARAM_t];
} mpc_share_t;

void xor_equals(mpc_share_t* x, const mpc_share_t* y) {
  uint8_t* x_raw = (uint8_t*) x;
  uint8_t const* y_raw = (uint8_t const*) y;
  for (uint64_t i = 0; i < sizeof(mpc_share_t); i++) {
    x_raw[i] ^= y_raw[i];
  }
}

#if 0
typedef struct {
  uint8_t share_raw[(params::s_A_size) /* s_A */ + params::w * 2 /* q_poly, p_poly */ +
                    params::t * 3 /* a, b, c */ * 3 /* 3 bytes for GF2p24 */];
} mpc_share_raw_t;

inline void operator^=(mpc_share_raw_t& x, mpc_share_raw_t y) {
  for (uint64_t i = 0; i < sizeof(x.share_raw); i++) {
    x.share_raw[i] ^= y.share_raw[i];
  }
}
#endif

/** @brief expand a leaf share from seed (any leaf except the last one) */
void expand_mpc_share_from_seed(mpc_share_t* share, const seed_t seed) {
  XOF_CTX rng_ctx;
  sdith_rng_create_xof_ctx(&rng_ctx, (void*)seed, PARAM_seed_size);
  sdith_xof_next_bytes(&rng_ctx, share, sizeof(mpc_share_t));
  for (uint64_t i_d=0; i_d < PARAM_d; ++i_d) {
    for (uint64_t i = 0; i < PARAM_t; i++) {
      share->a[i_d][i] &= PAR_fpoint_mask;
      share->b[i_d][i] &= PAR_fpoint_mask;
      share->c[i_d][i] &= PAR_fpoint_mask;
    }
  }
  sdith_rng_free_xof_ctx(&rng_ctx);
}

void expand_mpc_share_from_seed4(mpc_share_t** share, seed_t *seed) {
  
  XOF4_CTX* rng_ctx = sdith_rng_create_xof4_ctx((void**)seed, PARAM_seed_size);
  sdith_xof4_next_bytes(rng_ctx, (void**)share, sizeof(mpc_share_t));
  for (uint64_t k = 0; k < 4; ++k) {
    for (uint64_t i_d=0; i_d < PARAM_d; ++i_d) {
      for (uint64_t i = 0; i < PARAM_t; i++) {
        share[k]->a[i_d][i] &= PAR_fpoint_mask;
        share[k]->b[i_d][i] &= PAR_fpoint_mask;
        share[k]->c[i_d][i] &= PAR_fpoint_mask;
      }
    }
  }
  sdith_rng_free_xof4_ctx(rng_ctx);
  
}

/** @brief expand the random part of the last leaf share (i.e. a and b) from seed */
void expand_last_mpc_share_from_seed(mpc_share_t* share, const seed_t seed) {
  // cur_share.a, cur_share.b <- PRNG()
  XOF_CTX rng_ctx;
  sdith_rng_create_xof_ctx(&rng_ctx, (void*)seed, PARAM_seed_size);
  sdith_xof_next_bytes(&rng_ctx, &share->a, PARAM_d * PARAM_t * 2 * sizeof(fpoints_t));
  for (uint64_t i_d=0; i_d < PARAM_d; ++i_d) {
    for (uint64_t i = 0; i < PARAM_t; ++i) {
      share->a[i_d][i] &= PAR_fpoint_mask;
      share->b[i_d][i] &= PAR_fpoint_mask;
    }
  }
  sdith_rng_free_xof_ctx(&rng_ctx);
}

/**
 * @brief commit a leaf share from seed (any leaf except the last one) -- sha3 commitments
 * @param commitment : output commitment (of size params::commit_size)
 * @param leaf_seed  : the leaf seed (of size params::seed_size)
 * @param salt       : the global salt
 * @param iteration  : iteration between 0 and tau-1
 * @param leaf_idx   : leaf index between 0 and max_leaf
 */
void commit_leaf(void* commitment, void const* leaf_seed, void const* leaf_rho, void const* salt, uint16_t iteration, uint16_t leaf_idx) {
  HASH_CTX share_commit_ctx;

  sdith_hash_create_hash_ctx(&share_commit_ctx, HASH_COM);
  // H0(salt || e || i ...)
  sdith_hash_digest_update(&share_commit_ctx, salt, PARAM_salt_size);
  sdith_hash_digest_update(&share_commit_ctx, &iteration, sizeof(iteration));
  sdith_hash_digest_update(&share_commit_ctx, &leaf_idx, sizeof(leaf_idx));
  // H0(... || rho ...)
  sdith_hash_digest_update(&share_commit_ctx, leaf_rho, PARAM_rho_size);
  // H0(... || state)
  sdith_hash_digest_update(&share_commit_ctx, leaf_seed, PARAM_seed_size);
  sdith_hash_final(&share_commit_ctx, commitment);
  sdith_hash_free_hash_ctx(&share_commit_ctx);
}

void commit4_leaf(void **commitment, void **leaf_seed, void **leaf_rho, void * salt, uint16_t iteration, uint16_t *leaf_idx) {
  HASH4_CTX share_commit4_ctx;
  sdith_hash_create_hash4_ctx(&share_commit4_ctx, HASH_COM);
  // H0(salt || e || i ...)
  void *salts[4] = {salt, salt, salt, salt};
  sdith_hash4_digest_update(&share_commit4_ctx, salts, PARAM_salt_size);
  void *iters[4] = {&iteration, &iteration, &iteration, &iteration};
  sdith_hash4_digest_update(&share_commit4_ctx, iters, sizeof(iteration));
  void *leaves[4] = {&leaf_idx[0], &leaf_idx[1], &leaf_idx[2], &leaf_idx[3]};
  sdith_hash4_digest_update(&share_commit4_ctx, leaves, sizeof(uint16_t));
  // H0(... || rho ...)
  sdith_hash4_digest_update(&share_commit4_ctx, leaf_rho, PARAM_rho_size);
  // H0(... || state)
  sdith_hash4_digest_update(&share_commit4_ctx, leaf_seed, PARAM_seed_size);
  sdith_hash4_final(&share_commit4_ctx, commitment);
  sdith_hash_free_hash4_ctx(&share_commit4_ctx);
}

/**
 * @brief expand a leaf share from seed (any leaf except the last one)
 * @param commitment : output commitment (of size params::commit_size)
 * @param leaf_seed  : the leaf seed (of size params::seed_size)
 * @param aux        : the auxiliary share (only s_A, p_poly, q_poly and c are accessed)
 * @param salt       : the global salt
 * @param iteration  : iteration between 0 and tau-1
 */
void commit_last_leaf(void* commitment, void const* leaf_seed, void const* leaf_rho, mpc_share_t const* aux, void const* salt,
                      uint16_t iteration) {
  HASH_CTX share_commit_ctx;
  sdith_hash_create_hash_ctx(&share_commit_ctx, HASH_COM);
  // H0(salt || e || i ...)
  sdith_hash_digest_update(&share_commit_ctx, salt, PARAM_salt_size);
  sdith_hash_digest_update(&share_commit_ctx, &iteration, sizeof(iteration));
  uint16_t leaf_idx = (1 << PARAM_D) - 1;
  sdith_hash_digest_update(&share_commit_ctx, &leaf_idx, sizeof(leaf_idx));
  // H0(... || rho ...)
  sdith_hash_digest_update(&share_commit_ctx, leaf_rho, PARAM_rho_size);
  // H0(... || state || s_A || q_poly || p_poly || c)
  sdith_hash_digest_update(&share_commit_ctx, leaf_seed, PARAM_seed_size);
  sdith_hash_digest_update(&share_commit_ctx, aux, PARAM_k + PARAM_d * PAR_wd * 2);
  sdith_hash_digest_update(&share_commit_ctx, aux->c, PARAM_d * PARAM_t * sizeof(fpoints_t));
  sdith_hash_final(&share_commit_ctx, commitment);
  sdith_hash_free_hash_ctx(&share_commit_ctx);
}

/**
 * @brief commit the last leaf share from seed  + aux
 * @param commitment : output commitment (of size params::commit_size)
 * @param leaf_seed  : the leaf seed (of size params::seed_size)
 * @param aux        : the auxiliary share (only s_A, p_poly, q_poly and c are accessed)
 * @param salt       : the global salt
 * @param iteration  : iteration between 0 and tau-1
 */
void commit_last_leaf2(void* commitment, void const* leaf_seed, void const* leaf_rho, aux_share_t const* aux, void const* salt,
                      uint16_t iteration) {
  HASH_CTX share_commit_ctx;
  sdith_hash_create_hash_ctx(&share_commit_ctx, HASH_COM);
  // H0(salt || e || i ...)
  sdith_hash_digest_update(&share_commit_ctx, salt, PARAM_salt_size);
  sdith_hash_digest_update(&share_commit_ctx, &iteration, sizeof(iteration));
  uint16_t leaf_idx = (1 << PARAM_D) - 1;
  sdith_hash_digest_update(&share_commit_ctx, &leaf_idx, sizeof(leaf_idx));
  // H0(... || rho ...)
  sdith_hash_digest_update(&share_commit_ctx, leaf_rho, PARAM_rho_size);
  // H0(... || state || s_A || q_poly || p_poly || c)
  sdith_hash_digest_update(&share_commit_ctx, leaf_seed, PARAM_seed_size);
  sdith_hash_digest_update(&share_commit_ctx, aux, PARAM_k + PARAM_d * PAR_wd * 2);
  sdith_hash_digest_update(&share_commit_ctx, aux->c, PARAM_d * PARAM_t * sizeof(fpoints_t));
  sdith_hash_final(&share_commit_ctx, commitment);
  sdith_hash_free_hash_ctx(&share_commit_ctx);
}

/**
 * @brief expands the full treePRG from seed
 * @param sk                the SD (expanded) secret key
 * @param root_seed         the root seed (of size params::seed_size)
 * @param salt              the salt (of size params::salt_size)
 * @param msg_commit_ctx    hash context to compute h1
 * @param iteration         iteration number (between 0 and tau-1)
 * @param aux               output: auxiliary [must be memset to zero before calling the function]
 * @param sum_share         output: sum of all shares (a.k.a. plaintext) [must be memset to zero before calling the
 * function]
 * @param main_party_shares output: the D main party shares of index 0  [must be memset to zero before calling the
 * function]
 */
void expand_seed_binary_tree_bfs(sdith_full_key_t const* sk, seed_t root_seed,
                                 salt_t salt, HASH_CTX* msg_commit_ctx, uint8_t iteration,
                                 mpc_share_t* aux, mpc_share_t* sum_share,
                                 mpc_share_t (*main_party_shares)[/* N= */ 2],
                                 seed_t* seeds, commit_t* commits) {


  const uint64_t num_leafs = (1ul << PARAM_D);
  // Initialize TreePRG state
  TREE_PRG_CTX* tree_ctx = sdith_create_tree_prg_ctx(salt);
  memcpy(seeds[0], root_seed, PARAM_seed_size);
  // #ifndef NDEBUG
  //   cout << "Root seed: " << hexmem(seeds[0], params::seed_size) << std::endl;
  // #endif

  // deduce the seeds of the other levels
  seed_t *previous_level = seeds;
  seed_t *current_level = seeds + 1;
  uint64_t current_level_n = 2;
  uint64_t first_tweak = 1;

  for (uint64_t d = 1; d <= PARAM_D; ++d) {
    sdith_tree_prg_seed_expand(tree_ctx, current_level, previous_level, first_tweak, iteration, current_level_n);
    previous_level = current_level;
    current_level = previous_level + current_level_n;
    current_level_n <<= 1;
    first_tweak <<= 1;
  }
//  #ifndef NDEBUG
//    for (uint64_t i=1; i<2*num_leafs; ++i) {
//      fprintf(stdout, "Snode: %ld = %s\n", i, hexmem(seeds[i-1], PARAM_seed_size));
//    }
//  #endif

  
  // generate the last level commitments (except the last leaf)
  seed_t *const leaf_level = previous_level;
  // expand the leaf states (except the last one) from their seed and update the sum and the main shares
  mpc_share_t cur_share;
  
  for (uint64_t i = 0; i < num_leafs - 1; i += 4) {
    seed_t leaf_seed[4];
    uint8_t leaf_rho[4][PARAM_rho_size];
    void *leaf_rho_ptr[4] = {leaf_rho[0], leaf_rho[1], leaf_rho[2], leaf_rho[3]};
    sdith_tree_prg_leaf_expand4(tree_ctx, &leaf_level[i], leaf_seed, (uint8_t**)leaf_rho_ptr);
    
    void *leaf_seed_ptr[4] = {leaf_seed[0], leaf_seed[1], leaf_seed[2], leaf_seed[3]};
    uint16_t idx[4] = {i, i+1, i+2, i+3};
    void *commits_ptr[4] = {commits[i], commits[i+1], commits[i+2], commits[i+3]};
    commit4_leaf((void**)commits_ptr, (void**)leaf_seed_ptr, (void**)leaf_rho_ptr, salt, iteration, idx);

    mpc_share_t shares[4];
    memset(shares, 0, sizeof(shares));
    mpc_share_t *shares_ptr[4] = {&shares[0], &shares[1], &shares[2], &shares[3]};
    expand_mpc_share_from_seed4((mpc_share_t**)shares_ptr, (seed_t*)leaf_seed_ptr);

    for (uint64_t k = 0; k < 4; ++k) {
      // skip the last share as it's specially derived
      if (i + k >= num_leafs - 1) {
        break;
      }
      // add it to the aux
      xor_equals(aux, &shares[k]);
      // add it to the main parties
      for (uint64_t j = 0; j < PARAM_D; ++j) {
        if ((((i + k) >> (PARAM_D - 1 - j)) & 1) == 0) {
          xor_equals(&main_party_shares[j][0], &shares[k]);
        }
  #ifndef NDEBUG
        if (((i + k) >> (PARAM_D - 1 - j)) & 1) {
          xor_equals(&main_party_shares[j][1], &shares[k]);
        }
  #endif
      }
    }
  }
  
  
  // generate the last leaf share
  {
    seed_t leaf_seed;
    uint8_t leaf_rho[PARAM_rho_size];
    sdith_tree_prg_leaf_expand(tree_ctx, leaf_level[num_leafs - 1], leaf_seed, leaf_rho);

    // expand the random part from the seed
    memset(&cur_share, 0, sizeof(mpc_share_t));
    expand_last_mpc_share_from_seed(&cur_share, leaf_seed);
    xor_equals(aux, &cur_share);
    // correct the current share so that the sum is equal to the plaintext

    memcpy(sum_share->s_A, sk->s_A, PARAM_k);
    memcpy(sum_share->q_poly, sk->q_poly, PARAM_d * PAR_wd);
    memcpy(sum_share->p_poly, sk->p_poly, PARAM_d * PAR_wd);
    // cur_share.c = aux.c ^ (aux.a * aux.b)  element wise
    for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
      for (uint64_t i = 0; i < PARAM_t; i++) {
        sum_share->a[i_d][i] = aux->a[i_d][i];
        sum_share->b[i_d][i] = aux->b[i_d][i];
        sum_share->c[i_d][i] = fpoints_mul_ct(aux->a[i_d][i], aux->b[i_d][i]);
        cur_share.c[i_d][i] = sum_share->c[i_d][i] ^ aux->c[i_d][i];
      }
    }
    // cur_share.xa = aux.xa ^ sk.xa
    for (uint64_t i = 0; i < PARAM_k; i++) {
      cur_share.s_A[i] = aux->s_A[i] ^ sk->s_A[i];
    }
    // cur_share.Q = aux.Q ^ sk.Q
    // cur_share.P = aux.P ^ sk.P
    for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
      for (uint64_t i = 0; i < PAR_wd; i++) {
        cur_share.q_poly[i_d][i] = aux->q_poly[i_d][i] ^ sk->q_poly[i_d][i];
        cur_share.p_poly[i_d][i] = aux->p_poly[i_d][i] ^ sk->p_poly[i_d][i];
      }
    }
    *aux = cur_share;

#ifndef NDEBUG
    // add it to the main parties
    for (uint64_t j = 0; j < PARAM_D; ++j) {
      xor_equals(&main_party_shares[j][1], &cur_share);
    }
#endif

    // commit to the last leaf (sha3 commitment)
    commit_last_leaf(commits[num_leafs - 1], leaf_seed, leaf_rho, aux, salt, iteration);
  }
//  #ifndef NDEBUG
//    for (uint64_t i=0; i<num_leafs; ++i) {
//       fprintf(stdout, "Scommit: %ld = %s\n", i, hexmem(commits[i], PARAM_commit_size));
//    }
//  #endif

  // finally, generate the commitment hash
  sdith_hash_digest_update(msg_commit_ctx, commits, num_leafs * PARAM_commit_size);

 sdith_free_tree_prg_ctx(tree_ctx);
 
}

/*
Tweaks on the tree follow the hierarchical indexing of nodes, where
the root has number 1, left child is 2*parent and right child is 2*parent+1
root:              1
lvl1:        2           3
lvl2:    4      5      6       7
lvlD:  8  9  10  11  12  13  14  15

The hidden path index (between 0 and 2^D-1) corresponds to the leaf number:
path=6=0b110  will targets node index 2^D+path=8+6=14 in this example
*/

/**
 * @brief extracts the sibling path and challenge seed from the tree that hides position 'path'.
 * @param seeds             the seeds of the current tree
 * @param commits           the commits of the current tree
 * @param path              the leaf to hide (must be between 0 and 2^D-1)
 * @param tree_prg_seeds    output: the D sibling seeds from the root
 * @param chal_commit       output: commitment of the challenge node
 */
void walk_full_tree_prg_bfs(const seed_t* seeds, const commit_t* commits, const uint32_t path,
                       seed_t tree_prg_seeds[PARAM_D],  commit_t chal_commit) {
  uint32_t pp = path ^ 1;
  for (uint64_t j = 0; j < PARAM_D; ++j) {
    memcpy(tree_prg_seeds[PARAM_D - j - 1], seeds[(1 << (PARAM_D - j)) - 1 + pp], PARAM_seed_size);
    pp = (pp >> 1) ^ 1;
  }
  memcpy(chal_commit, commits[path], PARAM_commit_size);
  return;
}

/**
 * @brief extracts the sibling path and challenge seed from the tree that hides position 'path'.
 * @param root_seed         the root seed (of size params::seed_size)
 * @param salt              the salt (of size params::salt_size)
 * @param iteration
 * @param path              the leaf to hide (must be between 0 and 2^D-1)
 * @param aux               the auxiliary (only accessed if path == last leaf)
 * @param tree_prg_seeds    output: the D sibling seeds from the root
 * @param chal_commit       output: commitment of the challenge node
 */
void walk_tree_prg_bfs(const seed_t root_seed, const salt_t salt,
                       uint32_t iteration, uint32_t path, mpc_share_t const* aux,
                       seed_t tree_prg_seeds[PARAM_D],
                       commit_t chal_commit) {
  const uint64_t num_leafs = (1ul << PARAM_D);
  seed_t *seeds = malloc(PARAM_seed_size * 3);
  //posix_memalign((void**)&seeds, 32, PARAM_seed_size * 3);
  uint8_t *const prev_seed = seeds[0];
  seed_t *const cur_seeds = seeds + 1;

  // Initialize TreePRG state
  TREE_PRG_CTX* tree_ctx = sdith_create_tree_prg_ctx(salt);
  memcpy(prev_seed, root_seed, PARAM_seed_size);

  // deduce the seeds of the other levels
  uint64_t current_tweak = 1;                // root tweak
  for (uint64_t d = 1; d <= PARAM_D; ++d) {  // WIP!!!!
    sdith_tree_prg_seed_expand(tree_ctx, cur_seeds, prev_seed, current_tweak, iteration, 2);
    if ((path >> (PARAM_D - d)) & 1) {
      // store sibling left seed
      memcpy(tree_prg_seeds[d - 1], cur_seeds[0], PARAM_seed_size);
      // extend from right
      current_tweak = 2 * current_tweak + 1;
      memcpy(prev_seed, cur_seeds[1], PARAM_seed_size);
    } else {
      // store sibling right seed
      memcpy(tree_prg_seeds[d - 1], cur_seeds[1], PARAM_seed_size);
      // extend from left
      current_tweak = 2 * current_tweak;
      memcpy(prev_seed, cur_seeds[0], PARAM_seed_size);
    }
  }
  // finally, the challenge seed is the last live seed
  ASSERT_DRAMATICALLY(current_tweak == (1 << PARAM_D) + path, "Bug!!")
  // finally, the challenge seed is the last live seed
  // memcpy(chal_seed, prev_seed, params::seed_size);
  uint8_t const* const chal_seed = prev_seed;

  seed_t leaf_seed;
  uint8_t leaf_rho[PARAM_rho_size];
  sdith_tree_prg_leaf_expand(tree_ctx, (void*)chal_seed, leaf_seed, leaf_rho);

  if (path == num_leafs - 1) {
    commit_last_leaf(chal_commit, leaf_seed, leaf_rho, aux, salt, iteration);
  } else {
    commit_leaf(chal_commit, leaf_seed, leaf_rho, salt, iteration, path);
  }
  // #ifndef NDEBUG
  //   for (uint64_t i=0; i<params::D; ++i) {
  //     std::cout << "Wsibl: " << i << " " << hexmem(tree_prg_seeds[i], params::seed_size) << std::endl;
  //   }
  //   std::cout << "Wchall: " << path << " " << hexmem(chal_seed, params::seed_size) << std::endl;
  //   std::cout << "Wcommit: " << path << " " << hexmem(chal_commit, params::seed_size) << std::endl;
  // #endif
  free(seeds);
  sdith_free_tree_prg_ctx(tree_ctx);
}

/** @brief reconstructs the main party shares available from the sibling path
 * @param seed_hint         the sibling seeds as produced by walk_tree_prg
 * @param hint              index of the hidden leaf
 * @param salt              global salt
 * @param msg_commit_ctx    committment hash context
 * @param iteration         iteration number between 0 and tau-1
 * @param com               hidden leaf commitment
 * @param aux               auxiliary share
 * @param main_party_shares output: main party shares that are full
 */
// TODO: rename hint -> path and remove sum_share from the signature
void expand_seed_binary_tree_with_hint_bfs(seed_t *seed_hint, uint32_t hint,
                                           const salt_t salt, HASH_CTX* msg_commit_ctx,
                                           uint8_t iteration, commit_t com, aux_share_t const* aux,
                                           mpc_share_t (*main_party_shares)[/* N= */ 2]) {
  const uint64_t num_leafs = (1ul << PARAM_D);
  seed_t *seeds = malloc(num_leafs * PARAM_seed_size * 2);
  //posix_memalign((void**)&seeds, 32, num_leafs * PARAM_seed_size * 2);
  commit_t *commits = malloc(num_leafs * PARAM_commit_size);
  //posix_memalign((void**)&commits, 32, num_leafs * PARAM_commit_size);

  // Initialize TreePRG state
  TREE_PRG_CTX* tree_ctx = sdith_create_tree_prg_ctx(salt);
  // root seed is not available, we set it to zero
  memset(seeds, 0, PARAM_seed_size);
  seed_t *previous_level = seeds + 0;
  seed_t *current_level = seeds + 1;
  uint64_t current_level_n = 2;
  uint64_t current_tweak = 1;
  for (uint64_t d = 1; d <= PARAM_D; ++d) {
    // expand all the seeds from level d-1 to d
    sdith_tree_prg_seed_expand(tree_ctx, current_level, previous_level, current_tweak, iteration, current_level_n);
    // correct the seed using the sibling seed
    uint64_t index_sibling = (hint >> (PARAM_D - d)) ^ 1;
    memcpy(current_level[index_sibling], seed_hint[d - 1], PARAM_seed_size);
    previous_level = current_level;
    current_level = previous_level + current_level_n;
    current_level_n <<= 1;
    current_tweak <<= 1;
  }
//  #ifndef NDEBUG
//    fprintf(stdout, "Vpath: %d\n", hint);
//    for (uint64_t i=1; i<2*num_leafs; ++i) {
//      fprintf(stdout, "Vnode: %ld = %s\n", i, hexmem(seeds[i-1], PARAM_seed_size));
//    }
//  #endif

  // generate the last level commitments (except the last leaf)
  seed_t *const leaf_level = previous_level;
  // expand the leaf states (except the last one) from their seed and update the sum and the main shares
  mpc_share_t cur_share;
  for (uint64_t i = 0; i < num_leafs - 1; i += 4) {
    seed_t leaf_seed[4];
    uint8_t leaf_rho[4][PARAM_rho_size];
    void *leaf_rho_ptr[4] = {leaf_rho[0], leaf_rho[1], leaf_rho[2], leaf_rho[3]};
    sdith_tree_prg_leaf_expand4(tree_ctx, &leaf_level[i], leaf_seed, (uint8_t**)leaf_rho_ptr);
    
    void *leaf_seed_ptr[4] = {leaf_seed[0], leaf_seed[1], leaf_seed[2], leaf_seed[3]};
    uint16_t idx[4] = {i, i+1, i+2, i+3};
    void *commits_ptr[4] = {commits[i], commits[i+1], commits[i+2], commits[i+3]};
    commit4_leaf((void**)commits_ptr, (void**)leaf_seed_ptr, (void**)leaf_rho_ptr, (void*)salt, iteration, idx);

    mpc_share_t shares[4];
    memset(shares, 0, sizeof(shares));
    mpc_share_t *shares_ptr[4] = {&shares[0], &shares[1], &shares[2], &shares[3]};
    expand_mpc_share_from_seed4((mpc_share_t**)shares_ptr, (seed_t*)leaf_seed_ptr);

    for (uint64_t k = 0; k < 4; ++k) {
      // skip the last share as it's specially derived
      if (i + k >= num_leafs - 1) {
        break;
      }
      if (i + k == hint) continue;  // skip the hidden leaf
      // add it to the main parties
      for (uint64_t j = 0; j < PARAM_D; ++j) {
        if (((((i + k) ^ hint) >> (PARAM_D - 1 - j)) & 1) == 1) {
          // TODO get rid of leaf>>(D-1-j)&1 with the proper function signature
          xor_equals(&main_party_shares[j][((i + k) >> (PARAM_D - 1 - j)) & 1], &shares[k]);
        }
      }
    }
  }
  // replace the hidden leaf commit by the right one
  memcpy(commits[hint], com, PARAM_commit_size);
  // generate the last leaf commitment and share
  if (hint != num_leafs - 1) {
    seed_t leaf_seed;
    uint8_t leaf_rho[PARAM_rho_size];
    sdith_tree_prg_leaf_expand(tree_ctx, leaf_level[num_leafs - 1], leaf_seed, leaf_rho);

    // commitment of the last leaf is from seed + aux
    commit_last_leaf2(commits[num_leafs - 1], leaf_seed, leaf_rho, aux, salt, iteration);

    const uint32_t leaf = num_leafs - 1;
    // expand the random part from the seed
    memset(&cur_share, 0, sizeof(mpc_share_t));
    expand_last_mpc_share_from_seed(&cur_share, leaf_seed);
    // correct the current share using the aux
    memcpy(cur_share.s_A, aux->s_A, PARAM_k);
    memcpy(cur_share.q_poly, aux->q_poly, PARAM_d * PAR_wd);
    memcpy(cur_share.p_poly, aux->p_poly, PARAM_d * PAR_wd);
    memcpy(cur_share.c, aux->c, PARAM_d * PARAM_t * sizeof(fpoints_t));
    // add it to the main parties
    for (uint64_t j = 0; j < PARAM_D; ++j) {
      if ((((leaf ^ hint) >> (PARAM_D - 1 - j)) & 1) == 1) {
        // TODO get rid of leaf>>j&1 with the proper function signature
        xor_equals(&main_party_shares[j][(leaf >> (PARAM_D - 1 - j)) & 1], &cur_share);
      }
    }
  }

  // finally, generate the commitment hash
  sdith_hash_digest_update(msg_commit_ctx, commits, num_leafs * PARAM_commit_size);
  sdith_free_tree_prg_ctx(tree_ctx);
  free(seeds);
  free(commits);
}

#define compressed_helper_t 16  // TODO generalize or not?
typedef struct mpc_helper_struct {
  uint8_t compressed_pow_r[PAR_md + 1][compressed_helper_t]; __attribute__ ((aligned (32)))      // log(r^p)
  uint8_t compressed_eps_pow_r[PAR_md + 1][compressed_helper_t];  __attribute__ ((aligned (32))) // log(eps*r^p)
  uint32_t eps_f_r[PARAM_t];             // log(eps*F(r))
} mpc_helper_t;

void init(mpc_helper_t* helper) { memset(helper, 0, sizeof(*helper)); }
void helper_set_eps_pow_r(mpc_helper_t* helper, uint64_t p, uint64_t i_t, const fpoints_t value) {
  memcpy(&helper->compressed_eps_pow_r[p][PARAM_fpoint_size * i_t], &value, PARAM_fpoint_size);
}
void helper_set_pow_r(mpc_helper_t* helper, uint64_t p, uint64_t i_t, const fpoints_t value) {
  memcpy(&helper->compressed_pow_r[p][PARAM_fpoint_size * i_t], &value, PARAM_fpoint_size);
}
void helper_unpack_row(fpoints_t dst[PARAM_t], const uint8_t src[compressed_helper_t]) {
  memset(dst, 0, PARAM_t * sizeof(fpoints_t));
  for (uint64_t i_t = 0; i_t < PARAM_t; ++i_t) {
    memcpy(&dst[i_t], &src[PARAM_fpoint_size * i_t], PARAM_fpoint_size);
  }
}

void generate_mpc_helper(mpc_helper_t* mpc_helper, uint32_t* r, uint32_t* eps) {
  for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
    uint32_t log_eps = fpoints_dlog(eps[i_t]);
    uint32_t log_r = fpoints_dlog(r[i_t]);
    for (uint64_t p = 0; p < PAR_md + 1; p++) {
      uint32_t log_r_p = fpoints_dlog_pow_log(log_r, p);
      helper_set_eps_pow_r(mpc_helper, p, i_t, fpoints_dexp(fpoints_dlog_mul_log_log(log_eps, log_r_p)));
      helper_set_pow_r(mpc_helper, p, i_t, fpoints_dexp(log_r_p));
    }
  }
  gf256_vec_mat16cols_muladd(mpc_helper->eps_f_r, f_poly, mpc_helper->compressed_eps_pow_r, PAR_md + 1);
}

typedef int bool;

void mpc_compute_plain_broadcasts(mpc_share_t const *share, mpc_helper_t const helper[PARAM_d],
                            sdith_full_pubkey_t const* pk, fpoints_t alphas[PARAM_d][PARAM_t],
                            uint32_t betas [PARAM_d][PARAM_t]) {
  // x = s_A || s_B
  uint8_t x[PARAM_k + PAR_ha_nslice * 128];
  //posix_memalign((void**)&x, 32, PARAM_k + PAR_ha_nslice * 128);
  uint8_t* s_A = x;
  uint8_t* s_B = x + PARAM_k;
  memcpy(s_A, share->s_A, PARAM_k);
  memset(s_B, 0, PAR_ha_nslice * 128);
  memcpy(s_B, pk->compressed_pubkey.y, PAR_y_size);

  // s_B = s_A H_A ^ y
  // s_A H_A = [<s_A, H_A_0>, ..., <s_A, H_A_k-1>]
  for (uint64_t slice=0; slice<PAR_ha_nslice; ++slice) {
    gf256_vec_mat128cols_muladd_ct(s_B + 128*slice, s_A, pk->H_a[slice], PARAM_k);
  }

  for (uint64_t i_d=0; i_d < PARAM_d; ++i_d) {
    uint32_t s_r[PARAM_t] = {0};
    uint8_t compressed_s_r[16] __attribute__ ((aligned (32))) = {0};
    gf256_vec_mat16cols_muladd_ct(compressed_s_r, x + i_d * PAR_md, helper[i_d].compressed_pow_r, PAR_md);
    helper_unpack_row(s_r, compressed_s_r);

    memset(alphas[i_d], 0, PARAM_t * sizeof(fpoints_t));
    uint8_t compressed_alphas[16]  __attribute__ ((aligned (32))) = {0};
    memcpy(compressed_alphas, helper[i_d].compressed_eps_pow_r[PAR_wd], compressed_helper_t);
    gf256_vec_mat16cols_muladd_ct(compressed_alphas, share->q_poly[i_d], helper[i_d].compressed_eps_pow_r, PAR_wd);
    helper_unpack_row(alphas[i_d], compressed_alphas);

    for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
      alphas[i_d][i_t] ^= share->a[i_d][i_t];
      betas[i_d][i_t] = share->b[i_d][i_t] ^ s_r[i_t];
    }
  }
  return;
}

void mpc_compute_communications(mpc_share_t const* share, bool with_offsets, mpc_helper_t const helper[PARAM_d],
                                sdith_full_pubkey_t const* pk, const fpoints_t alphas[PARAM_d][PARAM_t],
                                const fpoints_t betas[PARAM_d][PARAM_t], fpoints_t alpha[PARAM_d][PARAM_t], fpoints_t beta[PARAM_d][PARAM_t],
                                fpoints_t v[PARAM_t], bool CONSTANT_TIME) {
  // x = s_A || s_B
  uint8_t x[PARAM_k + PAR_ha_nslice * 128];
  //posix_memalign((void**)&x, 32, PARAM_k + PAR_ha_nslice * 128);
  uint8_t* s_A = x;
  uint8_t* s_B = x + PARAM_k;
  memcpy(s_A, share->s_A, PARAM_k);
  memset(s_B, 0, PAR_ha_nslice * 128);
  if (with_offsets) {
    memcpy(s_B, pk->compressed_pubkey.y, PAR_y_size);
  }

  // s_B = s_A H_A ^ y
  // s_A H_A = [<s_A, H_A_0>, ..., <s_A, H_A_k-1>]
  if (CONSTANT_TIME) {  // if on a template
    for (uint64_t slice=0; slice<PAR_ha_nslice; ++slice) {
      gf256_vec_mat128cols_muladd_ct(s_B + 128*slice, s_A, pk->H_a[slice], PARAM_k);
    }
  } else {
    for (uint64_t slice=0; slice<PAR_ha_nslice; ++slice) {
      gf256_vec_mat128cols_muladd(s_B + 128*slice, s_A, pk->H_a[slice], PARAM_k);
    }
  }

  memset(v, 0, PARAM_t*sizeof(fpoints_t));
  for (uint64_t i_d=0; i_d < PARAM_d; ++i_d) {
    fpoints_t s_r[PARAM_t] = {0};
    uint8_t compressed_s_r[16]  __attribute__ ((aligned (32))) = {0};
    if (CONSTANT_TIME) {  // if on a template
      gf256_vec_mat16cols_muladd_ct(compressed_s_r, x + i_d * PAR_md, helper[i_d].compressed_pow_r, PAR_md);
    } else {
      gf256_vec_mat16cols_muladd(compressed_s_r, x + i_d * PAR_md, helper[i_d].compressed_pow_r, PAR_md);
    }
    helper_unpack_row(s_r, compressed_s_r);
#ifndef NDEBUG
    for (uint64_t i_t = 0; i_t < PARAM_t; ++i_t) {
      REQUIRE_DRAMATICALLY((s_r[i_t] & (~PAR_fpoint_mask)) == 0, "bug");
    }
#endif

    memset(alpha[i_d], 0, PARAM_t * sizeof(fpoints_t));
    uint8_t compressed_alpha[16]  __attribute__ ((aligned (32)))= {0};
    if (with_offsets) memcpy(compressed_alpha, helper[i_d].compressed_eps_pow_r[PAR_wd], compressed_helper_t);
    if (CONSTANT_TIME) {
      gf256_vec_mat16cols_muladd_ct(compressed_alpha, share->q_poly[i_d], helper[i_d].compressed_eps_pow_r, PAR_wd);
    } else {
      gf256_vec_mat16cols_muladd(compressed_alpha, share->q_poly[i_d], helper[i_d].compressed_eps_pow_r, PAR_wd);
    }
    helper_unpack_row(alpha[i_d], compressed_alpha);

    fpoints_t sh_p_r[PARAM_t] = {0};  ///< share of P(r_it)
    uint8_t compressed_sh_p_r[16] __attribute__ ((aligned (32))) = {0};
    if (CONSTANT_TIME) {
      gf256_vec_mat16cols_muladd_ct(compressed_sh_p_r, share->p_poly[i_d], helper[i_d].compressed_pow_r, PAR_wd);
    } else {
      gf256_vec_mat16cols_muladd(compressed_sh_p_r, share->p_poly[i_d], helper[i_d].compressed_pow_r, PAR_wd);
    }
    helper_unpack_row(sh_p_r, compressed_sh_p_r);

    for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
      alpha[i_d][i_t] ^= share->a[i_d][i_t];
      beta[i_d][i_t] = share->b[i_d][i_t] ^ s_r[i_t];

      if (CONSTANT_TIME) {
        v[i_t] ^= share->c[i_d][i_t] ^ fpoints_mul_ct(sh_p_r[i_t], helper[i_d].eps_f_r[i_t]) ^
                 fpoints_mul_ct(alphas[i_d][i_t], share->b[i_d][i_t]) ^ fpoints_mul_ct(betas[i_d][i_t], share->a[i_d][i_t]) ^
                 (with_offsets ? fpoints_mul_ct(alphas[i_d][i_t], betas[i_d][i_t]) : 0);
      } else {
        v[i_t] ^= share->c[i_d][i_t] ^ fpoints_mul(sh_p_r[i_t], helper[i_d].eps_f_r[i_t]) ^
                 fpoints_mul(alphas[i_d][i_t], share->b[i_d][i_t]) ^ fpoints_mul(betas[i_d][i_t], share->a[i_d][i_t]) ^
                 (with_offsets ? fpoints_mul(alphas[i_d][i_t], betas[i_d][i_t]) : 0);
      }
    }
  }
}

typedef struct sdith_ctx_struct {
  salt_t salt;
  HASH_CTX msg_commit_ctx;
  mpc_share_t main_party_shares[PARAM_tau][PARAM_D][/* N= */ 2];
  mpc_share_t aux[PARAM_tau];
  mpc_share_t sum_shares[PARAM_tau];
  seed_t root_seeds[PARAM_tau];
#ifdef FULL_TREE
  seed_t all_seeds[PARAM_tau][(1 << PARAM_D) * 2];
  commit_t all_commits[PARAM_tau][1 << PARAM_D];
#endif
} sdith_ctx_t;

sdith_ctx_t* new_sdith_ctx() { return (sdith_ctx_t*) malloc(sizeof(sdith_ctx_t)); }

void free_sdith_ctx(sdith_ctx_t* ctx) { free(ctx); }

void sign_offline(sdith_ctx_t* ctx, sdith_full_pubkey_t __attribute__((unused)) const* pk, sdith_full_key_t const* sk) {


  // Sample random salt
  randombytes(ctx->salt, PARAM_salt_size);
#ifndef NDEBUG
  fprintf(stdout, "Salt: %s\n", hexmem(ctx->salt, PARAM_salt_size));
#endif
  
  // Create message commitment context
  sdith_hash_create_hash_ctx(&ctx->msg_commit_ctx, HASH_H1);
  // We intentially skip H(m) and move it to the online phase
  // H1(pk || salt ...)
  sdith_hash_digest_update(&ctx->msg_commit_ctx , pk->compressed_pubkey.H_a_seed, PARAM_seed_size);
  sdith_hash_digest_update(&ctx->msg_commit_ctx , ctx->salt, PARAM_salt_size);

  
  // Initialize main party shares
  memset(ctx->main_party_shares, 0, sizeof(ctx->main_party_shares));

  // Initialize aux and sum shares data structure
  memset(ctx->aux, 0, sizeof(ctx->aux));
  memset(ctx->sum_shares, 0, sizeof(ctx->sum_shares));


#ifdef BENCHMARK
  sdith_bench_timer_init(&t0);
  sdith_bench_timer_start(&t0);
#endif
  // Expand the seed and obtain the MPC shares for tau times
  for (uint64_t e = 0; e < PARAM_tau; e++) {
    // Sample random TreePRG root seed
    randombytes(ctx->root_seeds[e], PARAM_seed_size);
#ifdef FULL_TREE
    seed_t* seeds = ctx->all_seeds[e];
    commit_t *commits = ctx->all_commits[e];
#else
    seed_t seeds[1<<(PARAM_D+1)];
    commit_t commits[1<<PARAM_D];
#endif
    expand_seed_binary_tree_bfs(sk, ctx->root_seeds[e], ctx->salt, &ctx->msg_commit_ctx , e, &ctx->aux[e],
                                &ctx->sum_shares[e], ctx->main_party_shares[e],
                              seeds, commits);
  }
#ifdef BENCHMARK
  sdith_bench_timer_end(&t0);
  sdith_bench_timer_count(&t0);
#endif
}

void sign_online(sdith_ctx_t* ctx, sdith_full_pubkey_t const* pk, sdith_full_key_t __attribute__((unused)) const* sk,
                 uint8_t const* msg, int msgBytes, uint8_t* sig, int* sigBytes) {
#ifndef IDS_3_ROUND
  // H1(... || m)
  sdith_hash_digest_update(&ctx->msg_commit_ctx , msg, msgBytes);
#endif
  // Obtain H1
  hash_t h1;
  sdith_hash_final(&ctx->msg_commit_ctx , h1);
  sdith_hash_free_hash_ctx(&ctx->msg_commit_ctx );
#ifndef NDEBUG
  fprintf(stdout, "H1: %s\n", hexmem(h1, sizeof(h1)));
#endif

  // Create H2 context
  sdith_hash_create_hash_ctx(&ctx->msg_commit_ctx, HASH_H2);
  // H2(m ...)
  sdith_hash_digest_update(&ctx->msg_commit_ctx , msg, msgBytes);
  // H2(... || salt ...)
  sdith_hash_digest_update(&ctx->msg_commit_ctx , ctx->salt, PARAM_salt_size);
  // H2(... || h1 ...)
  sdith_hash_digest_update(&ctx->msg_commit_ctx , h1, PARAM_hash_size);

  // Use H1 as PRNG seed to extract challenge values
  XOF_CTX chal_prg;
  sdith_rng_create_xof_ctx(&chal_prg, h1, PARAM_hash_size);

  fpoints_t r[PARAM_tau][PARAM_d][PARAM_t];
  fpoints_t eps[PARAM_tau][PARAM_d][PARAM_t];
  sdith_xof_next_bytes(&chal_prg, r, sizeof(r));
  sdith_xof_next_bytes(&chal_prg, eps, sizeof(eps));
  sdith_rng_free_xof_ctx(&chal_prg);
  for (uint64_t i = 0; i < PARAM_tau; i++) {
    for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
      for (uint64_t j = 0; j < PARAM_t; j++) {
        r[i][i_d][j] &= PAR_fpoint_mask;
        eps[i][i_d][j] &= PAR_fpoint_mask;
      }
    }
  }

  signature_t signature;
  memset(&signature, 0, sizeof(signature));

#ifdef BENCHMARK
  sdith_bench_timer_init(&t1);
  sdith_bench_timer_start(&t1);
#endif
  // precompute the mpc helpers
  mpc_helper_t mpc_helpers[PARAM_tau][PARAM_d];
  memset(mpc_helpers, 0, sizeof(mpc_helpers));
  for (uint64_t e = 0; e < PARAM_tau; e++) {
    for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
      generate_mpc_helper(&mpc_helpers[e][i_d], r[e][i_d], eps[e][i_d]);
    }
  }
//#ifndef NDEBUG
//  fprintf(stdout, "helpers hash: %s\n", hashmem(mpc_helpers, sizeof(mpc_helpers)));
//#endif

  fpoints_t alpha[PARAM_tau][PARAM_d][PARAM_t];
  fpoints_t beta[PARAM_tau][PARAM_d][PARAM_t];
  memset(alpha, 0, sizeof(alpha));
  memset(beta, 0, sizeof(beta));

  for (uint64_t e = 0; e < PARAM_tau; e++) {
    // uint32_t alphas[params::t];
    // uint32_t betas[params::t];
    mpc_compute_plain_broadcasts(&ctx->sum_shares[e], mpc_helpers[e], pk, alpha[e],
                           beta[e]);

    // H2(... || alpha[e] ...)
    sdith_hash_digest_update(&ctx->msg_commit_ctx , alpha[e], sizeof(alpha[e]));
    // H2(... || beta[e] ...)
    sdith_hash_digest_update(&ctx->msg_commit_ctx , beta[e], sizeof(beta[e]));

    for (uint64_t i = 0; i < PARAM_D; i++) {
      uint32_t sh_alpha[2][PARAM_d][PARAM_t];
      uint32_t sh_beta[2][PARAM_d][PARAM_t];
      uint32_t sh_v[2][PARAM_t];
      mpc_compute_communications(&ctx->main_party_shares[e][i][0], /* with_offsets= */ 0, mpc_helpers[e], pk,
                                       alpha[e], beta[e], sh_alpha[0], sh_beta[0], sh_v[0], 1);

#ifndef NDEBUG
      mpc_compute_communications(&ctx->main_party_shares[e][i][1], /* with_offsets= */ 1, mpc_helpers[e], pk,
                                       alpha[e], beta[e], sh_alpha[1], sh_beta[1], sh_v[1], 1);
      for (uint64_t i_d = 0; i_d < PARAM_d; i_d++) {
        for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
          REQUIRE_DRAMATICALLY((sh_alpha[0][i_d][i_t] ^ sh_alpha[1][i_d][i_t]) == alpha[e][i_d][i_t], "MPC alphas doesn't sum to alpha.");
          REQUIRE_DRAMATICALLY((sh_beta[0][i_d][i_t] ^ sh_beta[1][i_d][i_t]) == beta[e][i_d][i_t], "MPC betas doesn't sum to beta.");
        }
      }
      for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
        REQUIRE_DRAMATICALLY((sh_v[0][i_t] ^ sh_v[1][i_t]) == 0, "MPC v doesn't sum to 0.");
      }
#endif

      // we include only party 0 in the communications hash
      for (uint64_t party = 0; party < 1; party++) {
        // H2(... || sh_alpha[party][i_t] ...)
        sdith_hash_digest_update(&ctx->msg_commit_ctx , sh_alpha[party], sizeof(sh_alpha[party]));
        // H2(... || sh_beta[party] ...)
        sdith_hash_digest_update(&ctx->msg_commit_ctx , sh_beta[party], sizeof(sh_beta[party]));
        // H2(... || sh_v[party] ...)
        sdith_hash_digest_update(&ctx->msg_commit_ctx , sh_v[party], sizeof(sh_v[party]));
      }
    }
  }

  uint8_t* alpha_raw = (uint8_t*)alpha;
  uint8_t* beta_raw = (uint8_t*)beta;
  uint32_t ctr = 0;
  for (uint64_t i = 0; i < sizeof(alpha) / sizeof(fpoints_t); ++i) {
    memcpy(&((uint8_t*)signature.compressed_alpha)[ctr], &alpha_raw[i * sizeof(fpoints_t)], PARAM_fpoint_size);
    memcpy(&((uint8_t*)signature.compressed_beta)[ctr], &beta_raw[i * sizeof(fpoints_t)], PARAM_fpoint_size);
    ctr += PARAM_fpoint_size;
  }

#ifdef BENCHMARK
  sdith_bench_timer_end(&t1);
  sdith_bench_timer_count(&t1);
#ifdef IDS_3_ROUND
  sdith_bench_timer_end(&t3);
  sdith_bench_timer_count(&t3);
  sdith_bench_timer_init(&t4);
  sdith_bench_timer_start(&t4);
#endif
#endif

  // Obtain H2
  hash_t h2;
  sdith_hash_final(&ctx->msg_commit_ctx , h2);
  sdith_hash_free_hash_ctx(&ctx->msg_commit_ctx );
#ifndef NDEBUG
  fprintf(stdout, "H2: %s\n", hexmem(h2, sizeof(h2)));
#endif

#ifdef PROOF_OF_WORK
  uint8_t h2_pow_in[PARAM_hash_size * 2];
  memcpy(h2_pow_in, h2, PARAM_hash_size);
  uint8_t h2_pow[PARAM_hash_size];
  memcpy(h2_pow, h2, PARAM_hash_size);
  for (uint64_t j = 0; j < (1 << PAR_POW_iter); j++) {
    memcpy(&h2_pow_in[PARAM_hash_size], h2_pow, PARAM_hash_size);
    HASH_CTX* h2_pow_ctx = sdith_hash_create_hash_ctx(HASH_H2);
    sdith_hash_digest_update(h2_pow_ctx, h2_pow_in, sizeof(h2_pow_in));
    sdith_hash_final(h2_pow_ctx, h2_pow);
    sdith_hash_free_hash_ctx(h2_pow_ctx);
  }
#endif

  // Use H2 as PRNG seed to extract challenge values
#ifdef PROOF_OF_WORK
  sdith_rng_create_xof_ctx(&chal_prg, h2_pow, PARAM_hash_size);
#else
  sdith_rng_create_xof_ctx(&chal_prg, h2, PARAM_hash_size);
#endif
  uint64_t chal[PARAM_tau];
  sdith_xof_next_bytes(&chal_prg, chal, sizeof(chal));
  sdith_rng_free_xof_ctx(&chal_prg);

  memcpy(signature.salt, ctx->salt, sizeof(ctx->salt));
  memcpy(signature.h2, h2, sizeof(h2));

#ifdef BENCHMARK
  sdith_bench_timer_init(&t2);
  sdith_bench_timer_start(&t2);
#endif

  uint64_t chal_mask = (1 << PARAM_D) - 1; // mod leaf parties
  for (uint64_t e = 0; e < PARAM_tau; e++) {
    uint64_t challenge = chal[e] & chal_mask;

    if (challenge != chal_mask) {
      // Populate the signature with aux info
      memcpy(signature.aux[e].s_A, ctx->aux[e].s_A, sizeof(ctx->aux[e].s_A));
      memcpy(signature.aux[e].q_poly, ctx->aux[e].q_poly, sizeof(ctx->aux[e].q_poly));
      memcpy(signature.aux[e].p_poly, ctx->aux[e].p_poly, sizeof(ctx->aux[e].p_poly));
      memcpy(signature.aux[e].c, ctx->aux[e].c, sizeof(ctx->aux[e].c));
    }
#ifdef FULL_TREE
    walk_full_tree_prg_bfs(ctx->all_seeds[e], ctx->all_commits[e],challenge, signature.tree_prg_seeds[e], signature.com[e]);
#else
    walk_tree_prg_bfs(ctx->root_seeds[e], ctx->salt, e, challenge, &ctx->aux[e], signature.tree_prg_seeds[e],
                      signature.com[e]);
#endif
  }

#ifdef BENCHMARK
  sdith_bench_timer_end(&t2);
  sdith_bench_timer_count(&t2);
#endif

  memcpy(sig, &signature, sizeof(signature));
  *sigBytes = sizeof(signature);
}

void sign(sdith_full_pubkey_t const* pk, sdith_full_key_t const* sk, void const* msg, int msgBytes,
          void* sig, int* sigBytes) {
  sdith_ctx_t *ctx = new_sdith_ctx();

#ifdef BENCHMARK
  sdith_bench_timer_init(&t3);
  sdith_bench_timer_start(&t3);
#endif
  sign_offline(ctx, pk, sk);
#ifdef BENCHMARK
#ifndef IDS_3_ROUND
  sdith_bench_timer_end(&t3);
  sdith_bench_timer_count(&t3);
  sdith_bench_timer_init(&t4);
  sdith_bench_timer_start(&t4);
#endif
#endif
  sign_online(ctx, pk, sk, msg, msgBytes, sig, sigBytes);
  free(ctx);
#ifdef BENCHMARK
  sdith_bench_timer_end(&t4);
  sdith_bench_timer_count(&t4);
  fprintf(stdout, "\tcube gen and commit (ms): %lf\n", sdith_bench_timer_get(&t0));
  fprintf(stdout, "\tcompute plain broadcasts (ms): %lf\n", sdith_bench_timer_get(&t1));
  fprintf(stdout, "\tcompute sibling path (ms): %lf\n", sdith_bench_timer_get(&t2));
  fprintf(stdout, "sign_offline (ms): %lf\n", sdith_bench_timer_get(&t3));
  fprintf(stdout, "sign_online (ms): %lf\n", sdith_bench_timer_get(&t4));
#endif
}

int verify(sdith_full_pubkey_t const* pk, void const* msg, int msgBytes, void const* sig) {
  signature_t signature;
  memcpy(&signature, sig, sizeof(signature));

  XOF_CTX chal_prg;
  // Use H2 as PRNG seed to extract challenge values
#ifdef PROOF_OF_WORK
  uint8_t h2_pow_in[PARAM_hash_size * 2];
  memcpy(h2_pow_in, signature.h2, PARAM_hash_size);
  uint8_t h2_pow[PARAM_hash_size];
  memcpy(h2_pow, signature.h2, PARAM_hash_size);
  for (uint64_t j = 0; j < (1 << PAR_POW_iter); j++) {
    memcpy(&h2_pow_in[PARAM_hash_size], h2_pow, PARAM_hash_size);
    HASH_CTX *h2_pow_ctx = sdith_hash_create_hash_ctx(HASH_H2);
    sdith_hash_digest_update(h2_pow_ctx, h2_pow_in, sizeof(h2_pow_in));
    sdith_hash_final(h2_pow_ctx, h2_pow);
    sdith_hash_free_hash_ctx(h2_pow_ctx);
  }
  sdith_rng_create_xof_ctx(&chal_prg, h2_pow, PARAM_hash_size);
#else
  sdith_rng_create_xof_ctx(&chal_prg, signature.h2, PARAM_hash_size);
#endif
  uint64_t chal[PARAM_tau];
  sdith_xof_next_bytes(&chal_prg, chal, sizeof(chal));
  sdith_rng_free_xof_ctx(&chal_prg);

  // Create message commitment context
  HASH_CTX msg_commit_ctx;
  sdith_hash_create_hash_ctx(&msg_commit_ctx, HASH_H1);
  // H1(pk || salt ...)
  sdith_hash_digest_update(&msg_commit_ctx, pk->compressed_pubkey.H_a_seed, PARAM_seed_size);
  sdith_hash_digest_update(&msg_commit_ctx, signature.salt, PARAM_salt_size);

  // Create main party shares
  mpc_share_t main_party_shares[PARAM_tau][PARAM_D][/* N= */ 2];
  memset(main_party_shares, 0, sizeof(main_party_shares));

  // Create aux and sum shares data structure
  mpc_share_t sum_shares[PARAM_tau];
  memset(sum_shares, 0, sizeof(sum_shares));

  uint64_t chal_mask = (1 << PARAM_D) - 1;
  for (uint64_t e = 0; e < PARAM_tau; e++) {
    uint64_t challenge = chal[e] & chal_mask;

    if (challenge == chal_mask) {
      aux_share_t empty_share;
      memset(&empty_share, 0, sizeof(empty_share));
      if (memcmp(&signature.aux[e], &empty_share, sizeof(empty_share))) {
#ifndef NDEBUG
        fprintf(stderr, "Verification failed: Non-empty aux share.\n");
#endif
        return -1;
      }
    }

    // expand_seed_binary_tree_with_hint(signature.tree_prg_seeds[e], ~challenge & chal_mask, signature.salt,
    //                                       &msg_commit_ctx, e, signature.com[e], signature.aux[e], &sum_shares[e],
    //                                       main_party_shares[e]);
    expand_seed_binary_tree_with_hint_bfs(signature.tree_prg_seeds[e], challenge, signature.salt, &msg_commit_ctx, e,
                                          signature.com[e], &signature.aux[e], main_party_shares[e]);
//#ifndef NDEBUG
//    for (uint64_t d=0; d<PARAM_D; ++d) {
//      fprintf(stdout, "main party share: %ld %ld %d %s\n", e, d, 0,
//              hashmem(&main_party_shares[e][d][0], sizeof(mpc_share_t)));
//      fprintf(stdout, "main party share: %ld %ld %d %s\n", e, d, 1,
//              hashmem(&main_party_shares[e][d][1], sizeof(mpc_share_t)));
//    }
//#endif
  }

#ifndef IDS_3_ROUND
  // H1(... || m)
  sdith_hash_digest_update(&msg_commit_ctx, msg, msgBytes);
#endif

  // Obtain H1
  hash_t h1;
  sdith_hash_final(&msg_commit_ctx, h1);
  sdith_hash_free_hash_ctx(&msg_commit_ctx);

#ifndef NDEBUG
  fprintf(stdout, "H1: %s\n", hexmem(h1, sizeof(h1)));
#endif

  // Extract challenges r, eps, and i*
  // Use H1 as PRNG seed to extract challenge values
  sdith_rng_create_xof_ctx(&chal_prg, h1, PARAM_hash_size);

  fpoints_t r[PARAM_tau][PARAM_d][PARAM_t];
  fpoints_t eps[PARAM_tau][PARAM_d][PARAM_t];
  sdith_xof_next_bytes(&chal_prg, r, sizeof(r));
  sdith_xof_next_bytes(&chal_prg, eps, sizeof(eps));
  sdith_rng_free_xof_ctx(&chal_prg);
  for (uint64_t i = 0; i < PARAM_tau; i++) {
    for (uint64_t i_d=0; i_d<PARAM_d; ++i_d) {
      for (uint64_t j = 0; j < PARAM_t; j++) {
        r[i][i_d][j] &= PAR_fpoint_mask;
        eps[i][i_d][j] &= PAR_fpoint_mask;
      }
    }
  }

  // Create H2 context
  sdith_hash_create_hash_ctx(&msg_commit_ctx, HASH_H2);
  // H2(m ...)
  sdith_hash_digest_update(&msg_commit_ctx, msg, msgBytes);
  // H2(... || salt ...)
  sdith_hash_digest_update(&msg_commit_ctx, signature.salt, PARAM_salt_size);
  // H2(... || h1 ...)
  sdith_hash_digest_update(&msg_commit_ctx, h1, PARAM_hash_size);

  // precompute the mpc helpers
  mpc_helper_t mpc_helpers[PARAM_tau][PARAM_d];
  memset(mpc_helpers, 0, sizeof(mpc_helpers));
  for (uint64_t e = 0; e < PARAM_tau; e++) {
    for (uint64_t i_d = 0; i_d < PARAM_d; ++i_d) {
      generate_mpc_helper(&mpc_helpers[e][i_d], r[e][i_d], eps[e][i_d]);
    }
  }
//#ifndef NDEBUG
//  fprintf(stdout, "helpers hash: %s\n", hashmem(mpc_helpers, sizeof(mpc_helpers)));
//#endif

  fpoints_t alpha[PARAM_tau][PARAM_d][PARAM_t];
  fpoints_t beta[PARAM_tau][PARAM_d][PARAM_t];
  memset(alpha, 0, sizeof(alpha));
  memset(beta, 0, sizeof(beta));

  uint8_t* alpha_raw = (uint8_t*)alpha;
  uint8_t* beta_raw = (uint8_t*)beta;
  uint32_t ctr = 0;
  for (uint64_t i = 0; i < sizeof(alpha) / sizeof(fpoints_t); ++i) {
    memcpy(&alpha_raw[i * sizeof(fpoints_t)], &((uint8_t*)signature.compressed_alpha)[ctr], PARAM_fpoint_size);
    memcpy(&beta_raw[i * sizeof(fpoints_t)], &((uint8_t*)signature.compressed_beta)[ctr], PARAM_fpoint_size);
    ctr += PARAM_fpoint_size;
  }

  for (uint64_t e = 0; e < PARAM_tau; e++) {
    uint64_t challenge = chal[e] & chal_mask;

    // H2(... || alpha[e] ...)
    sdith_hash_digest_update(&msg_commit_ctx, alpha[e], sizeof(alpha[e]));
    // H2(... || beta[e] ...)
    sdith_hash_digest_update(&msg_commit_ctx, beta[e], sizeof(beta[e]));

    // // Add alpha, beta of i* to the global alpha, beta
    // for (uint64_t i_t = 0; i_t < params::t; i_t++) {
    //   alphas[i_t] ^= signature.sh_alpha[e][i_t];
    //   betas[i_t] ^= signature.sh_beta[e][i_t];
    // }

    for (uint64_t i = 0; i < PARAM_D; i++) {
      uint32_t sh_alpha[2][PARAM_d][PARAM_t];
      uint32_t sh_beta[2][PARAM_d][PARAM_t];
      uint32_t sh_v[2][PARAM_t];

      uint32_t chal_party = (challenge >> (PARAM_D - 1 - i)) & 1;
      // uint32_t chal_party = (challenge >> i) & 1;
      mpc_compute_communications(&main_party_shares[e][i][1 - chal_party], chal_party != 1, mpc_helpers[e], pk,
                                        alpha[e], beta[e], sh_alpha[1 - chal_party],
                                        sh_beta[1 - chal_party], sh_v[1 - chal_party], 0);
      if (chal_party == 0) {
        for (uint64_t i_d = 0; i_d < PARAM_d; i_d++) {
          for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
            sh_alpha[chal_party][i_d][i_t] = alpha[e][i_d][i_t] ^ sh_alpha[1 - chal_party][i_d][i_t];
            sh_beta[chal_party][i_d][i_t] = beta[e][i_d][i_t] ^ sh_beta[1 - chal_party][i_d][i_t];
          }
        }
        for (uint64_t i_t = 0; i_t < PARAM_t; i_t++) {
          sh_v[chal_party][i_t] = sh_v[1 - chal_party][i_t];
        }
      }

//#ifndef NDEBUG
//          fprintf(stdout, "VCOM %ld %ld %s\n", e, i, hexmem(&sh_alpha[0], sizeof(sh_alpha[0])));
//#endif
      for (uint64_t party = 0; party < 1; party++) {
        // H2(... || sh_alpha[party][i_t] ...)
        sdith_hash_digest_update(&msg_commit_ctx, sh_alpha[party], sizeof(sh_alpha[party]));
        // H2(... || sh_beta[party] ...)
        sdith_hash_digest_update(&msg_commit_ctx, sh_beta[party], sizeof(sh_beta[party]));
        // H2(... || sh_v[party] ...)
        sdith_hash_digest_update(&msg_commit_ctx, sh_v[party], sizeof(sh_v[party]));
      }
    }
  }

  // Obtain H2
  hash_t h2;
  sdith_hash_final(&msg_commit_ctx, h2);
  sdith_hash_free_hash_ctx(&msg_commit_ctx);
  if (memcmp(h2, signature.h2, PARAM_hash_size)) {
#ifndef NDEBUG
  fprintf(stderr, "Verification failed: H2 should be equal.\n");
#endif
    return -1;
  }

  return 0;
}

int sig_size() { return sizeof(signature_t); }

void field_init() {
  gf256_init(0);
  gf2p32_mul_init();
}
