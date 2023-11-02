#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"

#include "fips202_wrapper.h"

#include "params.h"

#include "api.h"
#include "randombytes.h"

#include "meds.h"

#include "seed.h"
#include "util.h"
#include "bitstream.h"

#include "matrixmod.h"

#define CEILING(x,y) (((x) + (y) - 1) / (y))


int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  )
{
  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);


  pmod_mat_t G_data[MEDS_k * MEDS_m * MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];

  uint8_t sigma_G0[MEDS_pub_seed_bytes];
  uint8_t sigma[MEDS_sec_seed_bytes];

  XOF((uint8_t*[]){sigma_G0, sigma},
       (size_t[]){MEDS_pub_seed_bytes, MEDS_sec_seed_bytes},
       delta, MEDS_sec_seed_bytes,
       2);

  LOG_VEC(sigma, MEDS_sec_seed_bytes);
  LOG_VEC_FMT(sigma_G0, MEDS_pub_seed_bytes, "sigma_G0");


  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, sigma_G0, MEDS_pub_seed_bytes);

  LOG_MAT(G[0], MEDS_k, MEDS_m*MEDS_n);

  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_m * MEDS_m];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  for (int i = 1; i < MEDS_s; i++)
  {
    pmod_mat_t A[MEDS_m * MEDS_m] = {0};
    pmod_mat_t B[MEDS_n * MEDS_n] = {0};

    while (1 == 1) // redo generation for this index until success
    {
      uint8_t sigma_Ti[MEDS_sec_seed_bytes];
      uint8_t sigma_a[MEDS_sec_seed_bytes];

      XOF((uint8_t*[]){sigma_a, sigma_Ti, sigma},
          (size_t[]){MEDS_sec_seed_bytes, MEDS_sec_seed_bytes, MEDS_sec_seed_bytes},
          sigma, MEDS_sec_seed_bytes,
          3);

      pmod_mat_t Ti[MEDS_k * MEDS_k];

      rnd_inv_matrix(Ti, MEDS_k, MEDS_k, sigma_Ti, MEDS_sec_seed_bytes);

      GFq_t Amm;

      {
        keccak_state Amm_shake;
        shake256_absorb_once(&Amm_shake, sigma_a, MEDS_sec_seed_bytes);

        Amm = rnd_GF(&Amm_shake);
      }


      LOG_MAT(Ti, MEDS_k, MEDS_k);
      LOG_VAL(Amm);


      pmod_mat_t G0prime[MEDS_k * MEDS_m * MEDS_n];

      pmod_mat_mul(G0prime, MEDS_k, MEDS_m * MEDS_n,
          Ti, MEDS_k, MEDS_k,
          G[0], MEDS_k, MEDS_m * MEDS_n);

      LOG_MAT(G0prime, MEDS_k, MEDS_m * MEDS_n);


      if (solve(A, B_inv[i], G0prime, Amm) < 0)
      {
        LOG("no sol");
        continue;
      }

      if (pmod_mat_inv(B, B_inv[i], MEDS_n, MEDS_n) < 0)
      {
        LOG("no inv B");
        continue;
      }

      if (pmod_mat_inv(A_inv[i], A, MEDS_m, MEDS_m) < 0)
      {
        LOG("no inv A_inv");
        continue;
      }

      LOG_MAT_FMT(A, MEDS_m, MEDS_m, "A[%i]", i);
      LOG_MAT_FMT(A_inv[i], MEDS_m, MEDS_m, "A_inv[%i]", i);
      LOG_MAT_FMT(B, MEDS_n, MEDS_n, "B[%i]", i);
      LOG_MAT_FMT(B_inv[i], MEDS_n, MEDS_n, "B_inv[%i]", i);


      pi(G[i], A, B, G[0]);


      if (pmod_mat_syst_ct(G[i], MEDS_k, MEDS_m*MEDS_n) != 0)
      {
        LOG("redo G[%i]", i);
        continue; // Not systematic; try again for index i.
      }

      LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

      // successfull generated G[s]; break out of while loop
      break;
    }
  }


  // copy pk data
  {
    uint8_t *tmp_pk = pk;

    memcpy(tmp_pk, sigma_G0, MEDS_pub_seed_bytes);
    LOG_VEC(tmp_pk, MEDS_pub_seed_bytes, "sigma_G0 (pk)");
    tmp_pk += MEDS_pub_seed_bytes;

    bitstream_t bs;

    bs_init(&bs, tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = (MEDS_m-1)*MEDS_n; j < MEDS_m*MEDS_n; j++)
        bs_write(&bs, G[si][MEDS_m*MEDS_n + j], GFq_bits);

      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
          bs_write(&bs, G[si][r*MEDS_m*MEDS_n + j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_VEC(tmp_pk, MEDS_PK_BYTES - MEDS_pub_seed_bytes, "G[1:] (pk)");
    tmp_pk += MEDS_PK_BYTES - MEDS_pub_seed_bytes;

    LOG_HEX(pk, MEDS_PK_BYTES);

    if (MEDS_PK_BYTES != MEDS_pub_seed_bytes + bs.byte_pos + (bs.bit_pos > 0 ? 1 : 0))
      return -1;
  }

  // copy sk data
  {
    memcpy(sk, delta, MEDS_sec_seed_bytes);
    memcpy(sk + MEDS_sec_seed_bytes, sigma_G0, MEDS_pub_seed_bytes);

    bitstream_t bs;

    bs_init(&bs, sk + MEDS_sec_seed_bytes + MEDS_pub_seed_bytes, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        bs_write(&bs, A_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        bs_write(&bs, B_inv[si][j], GFq_bits);

      bs_finalize(&bs);
    }

    LOG_HEX(sk, MEDS_SK_BYTES);
  }

  return 0;
}

int _crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  )
{
  uint8_t delta[MEDS_sec_seed_bytes];

  randombytes(delta, MEDS_sec_seed_bytes);


  // skip secret seed
  sk += MEDS_sec_seed_bytes;

  pmod_mat_t G_0[MEDS_k * MEDS_m * MEDS_n];


  rnd_sys_mat(G_0, MEDS_k, MEDS_m*MEDS_n, sk, MEDS_pub_seed_bytes);

  sk += MEDS_pub_seed_bytes;


  pmod_mat_t A_inv_data[MEDS_s * MEDS_m * MEDS_m];
  pmod_mat_t B_inv_data[MEDS_s * MEDS_n * MEDS_n];

  pmod_mat_t *A_inv[MEDS_s];
  pmod_mat_t *B_inv[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
  {
    A_inv[i] = &A_inv_data[i * MEDS_m * MEDS_m];
    B_inv[i] = &B_inv_data[i * MEDS_n * MEDS_n];
  }

  // Load secret key matrices.
  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t*)sk, MEDS_SK_BYTES - MEDS_sec_seed_bytes - MEDS_pub_seed_bytes);

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        A_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    for (int si = 1; si < MEDS_s; si++)
    {
      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        B_inv[si][j] = bs_read(&bs, GFq_bits);

      bs_finalize(&bs);
    }

    bs_finalize(&bs);
  }


  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT(A_inv[i], MEDS_m, MEDS_m);

  for (int i = 1; i < MEDS_s; i++)
    LOG_MAT(B_inv[i], MEDS_n, MEDS_n);

  LOG_MAT(G_0, MEDS_k, MEDS_m*MEDS_m);

  LOG_VEC(delta, MEDS_sec_seed_bytes);


  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};
  uint8_t alpha[MEDS_st_salt_bytes];

  uint8_t *rho = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(0,0)];

  XOF((uint8_t*[]){rho, alpha},
      (size_t[]){MEDS_st_seed_bytes, MEDS_st_salt_bytes},
      delta, MEDS_sec_seed_bytes,
      2);

  t_hash(stree, alpha, 0, 0);

  uint8_t *sigma = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  for (int i = 0; i < MEDS_t; i++)
  {
     LOG_HEX_FMT((&sigma[i*MEDS_st_seed_bytes]), MEDS_st_seed_bytes, "sigma[%i]", i);
  }

  pmod_mat_t A_tilde_data[MEDS_t * MEDS_m * MEDS_m];
  pmod_mat_t B_tilde_data[MEDS_t * MEDS_m * MEDS_m];

  pmod_mat_t *A_tilde[MEDS_t];
  pmod_mat_t *B_tilde[MEDS_t];

  for (int i = 0; i < MEDS_t; i++)
  {
    A_tilde[i] = &A_tilde_data[i * MEDS_m * MEDS_m];
    B_tilde[i] = &B_tilde_data[i * MEDS_n * MEDS_n];
  }


  uint8_t seed_buf[MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t)] = {0};
  memcpy(seed_buf, alpha, MEDS_st_salt_bytes);

  uint8_t *addr_pos = seed_buf + MEDS_st_salt_bytes + MEDS_st_seed_bytes;


  keccak_state h_shake;
  shake256_init(&h_shake);

  for (int i = 0; i < MEDS_t; i++)
  {
    pmod_mat_t G_tilde_ti[MEDS_k * MEDS_m * MEDS_m];

    while (1 == 1)
    {
      uint8_t sigma_A_tilde_i[MEDS_pub_seed_bytes];
      uint8_t sigma_B_tilde_i[MEDS_pub_seed_bytes];

      for (int j = 0; j < 4; j++)
        addr_pos[j] = (i >> (j*8)) & 0xff;

      memcpy(seed_buf + MEDS_st_salt_bytes, &sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes);

      LOG_HEX_FMT(seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + 4, "sigma_prime[%i]", i);

      XOF((uint8_t*[]){sigma_A_tilde_i, sigma_B_tilde_i, &sigma[i*MEDS_st_seed_bytes]},
           (size_t[]){MEDS_pub_seed_bytes, MEDS_pub_seed_bytes, MEDS_st_seed_bytes},
           seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + 4,
           3);

      LOG_HEX_FMT(sigma_A_tilde_i, MEDS_pub_seed_bytes, "sigma_A_tilde[%i]", i);
      rnd_inv_matrix(A_tilde[i], MEDS_m, MEDS_m, sigma_A_tilde_i, MEDS_pub_seed_bytes);

      LOG_HEX_FMT(sigma_B_tilde_i, MEDS_pub_seed_bytes, "sigma_B_tilde[%i]", i);
      rnd_inv_matrix(B_tilde[i], MEDS_n, MEDS_n, sigma_B_tilde_i, MEDS_pub_seed_bytes);

      LOG_MAT_FMT(A_tilde[i], MEDS_m, MEDS_m, "A_tilde[%i]", i);
      LOG_MAT_FMT(B_tilde[i], MEDS_n, MEDS_n, "B_tilde[%i]", i);


      pi(G_tilde_ti, A_tilde[i], B_tilde[i], G_0);


      LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

      if (pmod_mat_syst_ct(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n) == 0)
        break;
    }

    LOG_MAT_FMT(G_tilde_ti, MEDS_k, MEDS_m*MEDS_n, "G_tilde[%i]", i);

    bitstream_t bs;
    uint8_t bs_buf[CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8)];

    bs_init(&bs, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
        bs_write(&bs, G_tilde_ti[r * MEDS_m*MEDS_n + j], GFq_bits);

    shake256_absorb(&h_shake, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));
  }

  shake256_absorb(&h_shake, (uint8_t*)m, mlen);

  shake256_finalize(&h_shake);

  uint8_t digest[MEDS_digest_bytes];

  shake256_squeeze(digest, MEDS_digest_bytes, &h_shake);

  LOG_VEC(digest, MEDS_digest_bytes);


  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);

  LOG_VEC(h, MEDS_t);


  bitstream_t bs;

  bs_init(&bs, sm, MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8)));

  uint8_t *path = sm + MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8));

  t_hash(stree, alpha, 0, 0);

  stree_to_path(stree, h, path, alpha);

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      {
        pmod_mat_t mu[MEDS_m*MEDS_m];

        pmod_mat_mul(mu, MEDS_m, MEDS_m, A_tilde[i], MEDS_m, MEDS_m, A_inv[h[i]], MEDS_m, MEDS_m);

        LOG_MAT(mu, MEDS_m, MEDS_m);

        for (int j = 0; j < MEDS_m*MEDS_m; j++)
          bs_write(&bs, mu[j], GFq_bits);
      }

      bs_finalize(&bs);

      {
        pmod_mat_t nu[MEDS_n*MEDS_n];

        pmod_mat_mul(nu, MEDS_n, MEDS_n, B_inv[h[i]], MEDS_n, MEDS_n, B_tilde[i], MEDS_n, MEDS_n);

        LOG_MAT(nu, MEDS_n, MEDS_n);

        for (int j = 0; j < MEDS_n*MEDS_n; j++)
          bs_write(&bs, nu[j], GFq_bits);
      }

      bs_finalize(&bs);
    }
  }

  memcpy(sm + MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes, digest, MEDS_digest_bytes);
  memcpy(sm + MEDS_SIG_BYTES - MEDS_st_salt_bytes, alpha, MEDS_st_salt_bytes);
  memcpy(sm + MEDS_SIG_BYTES, m, mlen);

  *smlen = MEDS_SIG_BYTES + mlen;

  LOG_HEX(sm, MEDS_SIG_BYTES + mlen);

  return 0;
}

int _crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  )
{
  LOG_HEX(pk, MEDS_PK_BYTES);
  LOG_HEX(sm, smlen);

  pmod_mat_t G_data[MEDS_k*MEDS_m*MEDS_n * MEDS_s];
  pmod_mat_t *G[MEDS_s];

  for (int i = 0; i < MEDS_s; i++)
    G[i] = &G_data[i * MEDS_k * MEDS_m * MEDS_n];


  rnd_sys_mat(G[0], MEDS_k, MEDS_m*MEDS_n, pk, MEDS_pub_seed_bytes);

  {
    bitstream_t bs;

    bs_init(&bs, (uint8_t*)pk + MEDS_pub_seed_bytes, MEDS_PK_BYTES - MEDS_pub_seed_bytes);

    for (int i = 1; i < MEDS_s; i++)
    {
      for (int r = 0; r < MEDS_k; r++)
        for (int c = 0; c < MEDS_k; c++)
          if (r == c)
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 1);
          else
            pmod_mat_set_entry(G[i], MEDS_k, MEDS_m * MEDS_n, r, c, 0);

      for (int j = (MEDS_m-1)*MEDS_n; j < MEDS_m*MEDS_n; j++)
        G[i][MEDS_m*MEDS_n + j] = bs_read(&bs, GFq_bits);

      for (int r = 2; r < MEDS_k; r++)
        for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
          G[i][r*MEDS_m*MEDS_n + j] = bs_read(&bs, GFq_bits);

      for (int ii = 0; ii < MEDS_m; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][ii*MEDS_n + j] = ii == j ? 1 : 0;

      for (int ii = 0; ii < MEDS_m-1; ii++)
        for (int j = 0; j < MEDS_n; j++)
          G[i][MEDS_m*MEDS_n + ii*MEDS_n + j] = (ii+1) == j ? 1 : 0;

      bs_finalize(&bs);
    }
  }

  for (int i = 0; i < MEDS_s; i++)
    LOG_MAT_FMT(G[i], MEDS_k, MEDS_m*MEDS_n, "G[%i]", i);

  LOG_HEX_FMT(sm, MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8)), "munu");
  LOG_HEX_FMT(sm + MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8)),
      MEDS_max_path_len * MEDS_st_seed_bytes, "path");

  uint8_t *digest = (uint8_t*)sm + (MEDS_SIG_BYTES - MEDS_digest_bytes - MEDS_st_salt_bytes);

  uint8_t *alpha = (uint8_t*)sm + (MEDS_SIG_BYTES - MEDS_st_salt_bytes);

  LOG_HEX(digest, MEDS_digest_bytes);
  LOG_HEX(alpha, MEDS_st_salt_bytes);

  uint8_t h[MEDS_t];

  parse_hash(digest, MEDS_digest_bytes, h, MEDS_t);


  bitstream_t bs;

  bs_init(&bs, (uint8_t*)sm, MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8)));

  uint8_t *path = (uint8_t*)sm + MEDS_w * (CEILING(MEDS_m*MEDS_m * GFq_bits, 8) + CEILING(MEDS_n*MEDS_n * GFq_bits, 8));

  uint8_t stree[MEDS_st_seed_bytes * SEED_TREE_size] = {0};

  path_to_stree(stree, h, path, alpha);

  uint8_t *sigma = &stree[MEDS_st_seed_bytes * SEED_TREE_ADDR(MEDS_seed_tree_height, 0)];

  pmod_mat_t G_hat_i[MEDS_k*MEDS_m*MEDS_n];

  pmod_mat_t mu[MEDS_m*MEDS_m];
  pmod_mat_t nu[MEDS_n*MEDS_n];


  uint8_t seed_buf[MEDS_st_salt_bytes + MEDS_st_seed_bytes + sizeof(uint32_t)] = {0};
  memcpy(seed_buf, alpha, MEDS_st_salt_bytes);

  uint8_t *addr_pos = seed_buf + MEDS_st_salt_bytes + MEDS_st_seed_bytes;


  keccak_state shake;
  shake256_init(&shake);

  for (int i = 0; i < MEDS_t; i++)
  {
    if (h[i] > 0)
    {
      for (int j = 0; j < MEDS_m*MEDS_m; j++)
        mu[j] = bs_read(&bs, GFq_bits) % MEDS_p;

      bs_finalize(&bs);

      for (int j = 0; j < MEDS_n*MEDS_n; j++)
        nu[j] = bs_read(&bs, GFq_bits) % MEDS_p;

      bs_finalize(&bs);


      LOG_MAT_FMT(mu, MEDS_m, MEDS_m, "mu[%i]", i);
      LOG_MAT_FMT(nu, MEDS_n, MEDS_n, "nu[%i]", i);

      // Check if mu is invetible.
      {
        pmod_mat_t tmp_mu[MEDS_m*MEDS_m];

        memcpy(tmp_mu, mu, MEDS_m*MEDS_m*sizeof(GFq_t));

        if (pmod_mat_syst_ct(tmp_mu, MEDS_m, MEDS_m) != 0)
          return -1;
      }

      // Check if nu is invetible.
      {
        pmod_mat_t tmp_nu[MEDS_n*MEDS_n];

        memcpy(tmp_nu, nu, MEDS_n*MEDS_n*sizeof(GFq_t));

        if (pmod_mat_syst_ct(tmp_nu, MEDS_n, MEDS_n) != 0)
          return -1;
      }


      pi(G_hat_i, mu, nu, G[h[i]]);


      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);

      if (pmod_mat_syst_ct(G_hat_i, MEDS_k, MEDS_m*MEDS_n) < 0)
        return -1;

      LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
    }
    else
    {
      while (1 == 1)
      {
        LOG_VEC_FMT(&sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes, "seeds[%i]", i);

        uint8_t sigma_A_hat_i[MEDS_pub_seed_bytes];
        uint8_t sigma_B_hat_i[MEDS_pub_seed_bytes];

        for (int j = 0; j < 4; j++)
          addr_pos[j] = (i >> (j*8)) & 0xff;

        memcpy(seed_buf + MEDS_st_salt_bytes, &sigma[i*MEDS_st_seed_bytes], MEDS_st_seed_bytes);

        LOG_HEX_FMT(seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + 4, "sigma_prime[%i]", i);

        XOF((uint8_t*[]){sigma_A_hat_i, sigma_B_hat_i, &sigma[i*MEDS_st_seed_bytes]},
            (size_t[]){MEDS_pub_seed_bytes, MEDS_pub_seed_bytes, MEDS_st_seed_bytes},
            seed_buf, MEDS_st_salt_bytes + MEDS_st_seed_bytes + 4,
            3);

        pmod_mat_t A_hat_i[MEDS_m*MEDS_m];
        pmod_mat_t B_hat_i[MEDS_n*MEDS_n];

        LOG_HEX_FMT(sigma_A_hat_i, MEDS_pub_seed_bytes, "sigma_A_hat[%i]", i);
        rnd_inv_matrix(A_hat_i, MEDS_m, MEDS_m, sigma_A_hat_i, MEDS_pub_seed_bytes);

        LOG_HEX_FMT(sigma_B_hat_i, MEDS_pub_seed_bytes, "sigma_B_hat[%i]", i);
        rnd_inv_matrix(B_hat_i, MEDS_n, MEDS_n, sigma_B_hat_i, MEDS_pub_seed_bytes);

        LOG_MAT_FMT(A_hat_i, MEDS_m, MEDS_m, "A_hat[%i]", i);
        LOG_MAT_FMT(B_hat_i, MEDS_n, MEDS_n, "B_hat[%i]", i);


        pi(G_hat_i, A_hat_i, B_hat_i, G[0]);

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);


        if (pmod_mat_syst_ct(G_hat_i, MEDS_k, MEDS_m*MEDS_n) == 0)
        {
          LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
          break;
        }

        LOG_MAT_FMT(G_hat_i, MEDS_k, MEDS_m*MEDS_n, "G_hat[%i]", i);
      }
    }


    bitstream_t bs;
    uint8_t bs_buf[CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8)];

    bs_init(&bs, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));

    for (int r = 0; r < MEDS_k; r++)
      for (int j = MEDS_k; j < MEDS_m*MEDS_n; j++)
        bs_write(&bs, G_hat_i[r * MEDS_m*MEDS_n + j], GFq_bits);

    shake256_absorb(&shake, bs_buf, CEILING((MEDS_k * (MEDS_m*MEDS_n - MEDS_k)) * GFq_bits, 8));
  }

  shake256_absorb(&shake, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);

  shake256_finalize(&shake);

  uint8_t check[MEDS_digest_bytes];

  shake256_squeeze(check, MEDS_digest_bytes, &shake);

  if (memcmp(digest, check, MEDS_digest_bytes) != 0)
    return -1;

  memcpy(m, (uint8_t*)(sm + MEDS_SIG_BYTES), smlen - MEDS_SIG_BYTES);
  *mlen = smlen - MEDS_SIG_BYTES;

  return 0;
}

