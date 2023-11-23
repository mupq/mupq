#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "utils.h"

#include "biscuit.h"

#include "batch_tools.h"

#if !defined(DEGREE) || (DEGREE > 2)
#define HIGH_DEGREE
#endif

static void
extract_arg (void *output, size_t outlen, void *arg)
{
  expand_extract (output, outlen, arg);
}

#define batch_sample(dest, q, n, tape) \
  batch_generate (dest, q, n, extract_arg, tape)

#define circuit_sample(dest, q, n, m, d, tape) \
  circuit_generate (dest, q, n, m, d, extract_arg, tape)

static uint16_t
index_sample (int N, int lN, h_ctx_t h_ctx)
{
  const uint16_t mask = (1 << lN) - 1;
  uint8_t v8[2] = {0, 0};
  uint16_t v;
  int n = (lN <= 8) ? 1 : 2;

  do
    {
      expand_extract (v8, n, h_ctx);
      v = (uint16_t) (v8[0] | (v8[1] << 8)) & mask;
    }
  while (v >= N);

  return v;
}

#ifdef SEC_LEVEL
#define lambda SEC_LEVEL
#endif
#ifdef NB_ITERATIONS
#define tau NB_ITERATIONS
#endif
#ifdef NB_PARTIES
#define N NB_PARTIES
#define lN LOG2(NB_PARTIES)
#endif
#ifdef FIELD_SIZE
#define q FIELD_SIZE
#endif
#ifdef NB_VARIABLES
#define n NB_VARIABLES
#endif
#ifdef NB_EQUATIONS
#define m NB_EQUATIONS
#endif
#ifdef DEGREE
#define d DEGREE
#endif
#if defined(DEGREE) && defined(NB_EQUATIONS)
#define C ((DEGREE-1)*NB_EQUATIONS)
#define Cm ((DEGREE-1)*NB_EQUATIONS)
#endif
#if defined(FIELD_SIZE) && defined(NB_VARIABLES)
#define sklen ((LOG2(FIELD_SIZE)*NB_VARIABLES+7)>>3)
#define sklenX CONVX(sklen)
#endif
#if defined(FIELD_SIZE) && defined(NB_EQUATIONS)
#define pklen ((LOG2(FIELD_SIZE)*NB_EQUATIONS+7)>>3)
#define pklenX CONVX(pklen)
#endif
#if defined(FIELD_SIZE) && defined(DEGREE) && defined(NB_EQUATIONS)
#define Clen ((LOG2(FIELD_SIZE)*(DEGREE-1)*NB_EQUATIONS+7)>>3)
#define Cmlen ((LOG2(FIELD_SIZE)*(DEGREE-2)*NB_EQUATIONS+7)>>3)
#define ClenX CONVX(Clen)
#endif

int
keygen (uint8_t *sk, uint8_t *pk, const uint8_t *entropy,
        const params_t *params)
{
  int offset;

#ifndef SEC_LEVEL
  const int lambda = params->lambda;
#endif
#ifndef FIELD_SIZE
  const int q = params->q;
#endif
#ifndef NB_VARIABLES
  const int n = params->n;
#endif
#ifndef NB_EQUATIONS
  const int m = params->m;
#endif
#ifndef DEGREE
  const int d = params->d;
#endif
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
#ifndef COMPACT_SK
#if !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int C = (d - 1) * m;
#ifdef HIGH_DEGREE
  const int Cm = (d - 2) * m;
#endif
#endif
#if !defined(FIELD_SIZE) || !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int ClenX = batch_getlenX (q, C);
#endif
#endif

  const uint8_t *const seedF = entropy;
  const uint8_t *const seedS = entropy + (lambda >> 3);

  uintX_t s[sklenX];
  uintX_t t[pklenX];
#ifndef COMPACT_SK
  uintX_t y[ClenX];
#ifdef HIGH_DEGREE
  uintX_t z[ClenX];
#endif
#endif

  h_ctx_t expand_ctx;

  BATCH_PARAMS (q, n, m, d);

  /* Expand seedS into s */
  expand_init (expand_ctx, seedS, lambda >> 3);
  batch_sample (s, q, n, expand_ctx);

  /* Expand seedF and evaluate the circuit */
  expand_init (expand_ctx, seedF, lambda >> 3);
#ifndef COMPACT_SK
#ifdef HIGH_DEGREE
  eval_circuit_seed (y, z, t, s, q, n, m, d, extract_arg, expand_ctx);
#else
  eval_circuit_seed (y, NULL, t, s, q, n, m, d, extract_arg, expand_ctx);
#endif
#else
  eval_circuit_seed (NULL, NULL, t, s, q, n, m, d, extract_arg, expand_ctx);
#endif

#ifndef COMPACT_SK
  /* Copy seedF, s, t, y and z in sk */
  memcpy (sk, seedF, (lambda >> 3));
  sk += lambda >> 3;
  offset = 0;
  batch_export (sk, s, q, n, offset);
  offset += n;
  batch_export (sk, t, q, m, offset);
  offset += m;
  batch_export (sk, y, q, C, offset);
#ifdef HIGH_DEGREE
  offset += C;
  batch_export (sk, z, q, Cm, offset);
#endif
#else
  /* Copy seedF and seedS in sk */
  memcpy (sk, entropy, lambda >> 2);
#endif

  /* Copy seedF and t in pk */
  memcpy (pk, seedF, lambda >> 3);
  pk += lambda >> 3;
  offset = 0;
  batch_export (pk, t, q, m, offset);

  return 0;
}

int
sign (uint8_t *sig, const uint8_t *msg, uint64_t msglen, const uint8_t *sk,
      const uint8_t *entropy, const params_t *params)
{
  int e, i;
  int offset;

#ifndef SEC_LEVEL
  const int lambda = params->lambda;
#endif
#ifndef NB_ITERATIONS
  const int tau = params->tau;
#endif
#ifndef NB_PARTIES
  const int N = params->N;
  const int lN = ilog2 (N);
#endif
#ifndef FIELD_SIZE
  const int q = params->q;
#endif
#ifndef NB_VARIABLES
  const int n = params->n;
#endif
#ifndef NB_EQUATIONS
  const int m = params->m;
#endif
#ifndef DEGREE
  const int d = params->d;
#endif
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklen = batch_getlen (q, n);
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
#if !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int C = (d - 1) * m;
#ifdef HIGH_DEGREE
  const int Cm = (d - 2) * m;
#endif
#endif
#if !defined(FIELD_SIZE) || !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int Clen = batch_getlen (q, C);
  const int ClenX = batch_getlenX (q, C);
#ifdef HIGH_DEGREE
  const int Cmlen = batch_getlen (q, Cm);
#endif
#endif

  const uint8_t *const seedF = sk;
  const uint8_t *const sk_data = sk + (lambda >> 3);

  uint8_t *const salt = sig;
  uint8_t *const h1 = sig + (lambda >> 2);
  uint8_t *const h2 = sig + 2 * (lambda >> 2);
  uint8_t (*const sig_path)[lN * (lambda >> 3) + (lambda >> 2)]
    = (void *) (sig + 3 * (lambda >> 2));
  uint8_t (*const sig_com)[lN * (lambda >> 3) + (lambda >> 2)]
    = (void *) (sig + 3 * (lambda >> 2) + lN * (lambda >> 3));
  uint8_t *sigma;

  h_ctx_t hash_ctx;
  h_ctx_t expand_ctx;

  uint8_t root[tau][lambda >> 3];
  uint8_t com[tau][N][lambda >> 2];

  uintX_t a[tau][N][ClenX], c[tau][N][ClenX];
  uintX_t x[tau][N][ClenX], y[tau][N][ClenX], z[tau][N][ClenX];

  uintX_t sk_s[sklenX], sk_y[ClenX], sk_t[pklenX];
#ifdef HIGH_DEGREE
  uintX_t sk_z[ClenX];
#endif
  uintX_t f[(d + 1) * (pklenX + sklenX * m)];

  BATCH_PARAMS (q, n, m, d);

  /* Expand seedF to obtain the system */
  expand_init (expand_ctx, seedF, lambda >> 3);
  circuit_sample (f, q, n, m, d, expand_ctx);

#ifndef COMPACT_SK
  offset = 0;
  batch_import (sk_s, sk_data, q, n, offset);
  offset += n;
  batch_import (sk_t, sk_data, q, m, offset);
  offset += m;
  batch_import (sk_y, sk_data, q, C, offset);
#ifdef HIGH_DEGREE
  offset += C;
  batch_import (sk_z, sk_data, q, Cm, offset);
#endif
#else
  /* Expand seedS from sk to obtain secret s */
  expand_init (expand_ctx, sk_data, lambda >> 3);
  batch_sample (sk_s, q, n, expand_ctx);
#endif

#ifdef COMPACT_SK
  /* Evaluate the circuit for intermediate values */
#ifdef HIGH_DEGREE
  eval_circuit (sk_y, sk_z, sk_t, sk_s, q, n, m, d, f);
#else
  eval_circuit (sk_y, NULL, sk_t, sk_s, q, n, m, d, f);
#endif
#endif


  /* Phase 1: Commiting to the seeds and views of the parties */

  /* initialize offset to write to the signature */
  sigma = sig + 3 * (lambda >> 2) + tau * (lN * (lambda >> 3) + (lambda >> 2));
  offset = 0;

  /* Prepare PRF for salt and root seeds */
  prf_init (expand_ctx, entropy, lambda >> 2);
#ifndef COMPACT_SK
#ifdef HIGH_DEGREE
  prf_update (expand_ctx, sk, (lambda >> 3)
              + CDIV8 (batch_getbitlen (q, n + m + C + Cm)));
#else
  prf_update (expand_ctx, sk, (lambda >> 3)
              + CDIV8 (batch_getbitlen (q, n + m + C)));
#endif
#else
  prf_update (expand_ctx, sk, lambda >> 2);
#endif
  prf_update (expand_ctx, msg, msglen);
  prf_ready (expand_ctx);

  /* Generate salt and root seeds from PRF */
  prf_generate (salt, lambda >> 2, expand_ctx);
  prf_generate ((void *) root, tau * (lambda >> 3), expand_ctx);

  /* Initialize sigma1 with salt and msg */
  H1_init (hash_ctx, salt, msg, msglen, lambda);
  for (e = 0; e < tau; e++)
    {
      uintX_t s_e[N][sklenX];

      uintX_t delta_s_e[sklenX];
#ifdef HIGH_DEGREE
      uintX_t delta_z_e[ClenX];
#endif
      uintX_t delta_c_e[ClenX];

      uint8_t seed_e[N][lambda >> 3];

      /* use x[e][0] as temporary buffer to compute delta values */
      uintX_t *const tmp_e = x[e][0];

      /* Get N seeds from root_e */
      get_seeds ((void *) seed_e, root[e], salt, e, N, lambda);

      /* Prepare computation of delta_s_e, init with sk_s */
      batch_copy (delta_s_e, sk_s, q, n);
#ifdef HIGH_DEGREE
      /* Prepare computation of delta_z_e, init with sk_z */
      batch_copy (delta_z_e, sk_z, q, Cm);
#endif
      /* Prepare computation of delta_c_e with tmp_e, init with 0 */
      batch_clear (delta_c_e, q, C);
      batch_clear (tmp_e, q, C);
      for (i = 0; i < N; i++)
        {
          h_ctx_t tape_ei_ctx;

          /* Commit on seed_e[i] */
          commit (com[e][i], salt, e, i, seed_e[i], lambda);

          /* Expand random from tape_ei */
          expandtape_init (tape_ei_ctx, salt, e, i, seed_e[i], lambda);
          batch_sample (s_e[i], q, n, tape_ei_ctx);
#ifdef HIGH_DEGREE
          batch_sample (z[e][i], q, Cm, tape_ei_ctx);
#endif
          batch_sample (a[e][i], q, C, tape_ei_ctx);
          batch_sample (c[e][i], q, C, tape_ei_ctx);

          /* Iterated to obtain delta_s_e = sk_s - sum(s_e[i]) */
          batch_sub (delta_s_e, s_e[i], q, n);
#ifdef HIGH_DEGREE
          /* Iterated to obtain delta_z_e = sk_z - sum(z_e[i]) */
          batch_sub (delta_z_e, z[e][i], q, Cm);
#endif
          /* Iterated to obtain tmp_e = sum(a_e[i]) */
          batch_add (tmp_e, a[e][i], q, C);
          /* Iterated to obtain delta_c_e = - sum(c_e[i]) */
          batch_sub (delta_c_e, c[e][i], q, C);

          /* Add com[e][i] to sigma1 */
          H1_update (hash_ctx, com[e][i], lambda >> 2);
        }
      /* Obtain delta_c_e = sk_y * sum(a_e[i]) - sum(c_e[i]) */
      batch_mul (tmp_e, sk_y, q, C);
      batch_add (delta_c_e, tmp_e, q, C);

      /* Correct s_e[0], z[e][0] and c[e][0] */
      batch_add (s_e[0], delta_s_e, q, n);
#ifdef HIGH_DEGREE
      batch_add (z[e][0], delta_z_e, q, Cm);
#endif
      batch_add (c[e][0], delta_c_e, q, C);

      for (i = 0; i < N; i++)
        {
          /* Compute the shares x[e][i], y[e][i], z[e][i] */
          /* involved in multiplications */
          linear_circuit (x[e][i], y[e][i], z[e][i], s_e[i], sk_t, i,
                          q, n, m, d, f);
        }

      /* Add delta values to sigma1 */
      H1_update (hash_ctx, (void *) delta_s_e, sklen);
#ifdef HIGH_DEGREE
      H1_update (hash_ctx, (void *) delta_z_e, Cmlen);
#endif
      H1_update (hash_ctx, (void *) delta_c_e, Clen);

      /* Export delta values in signature */
      batch_export (sigma, delta_s_e, q, n, offset);
      offset += n;
#ifdef HIGH_DEGREE
      batch_export (sigma, delta_z_e, q, Cm, offset);
      offset += Cm;
#endif
      batch_export (sigma, delta_c_e, q, C, offset);
      offset += C;
    }
  /* Finalize computation of h1 using sigma1 */
  H1_final (h1, hash_ctx, lambda);


  /* Phase 2: Challenging the checking protocol */

  /* Prepare expansion of epsilon_e challenges */
  expand_init (expand_ctx, h1, lambda >> 2);


  /* Phase 3: Commit to simulation of the checking protocol */

  /* Initialize sigma2 with salt and h1 */
  H2_init (hash_ctx, salt, h1, lambda);
  for (e = 0; e < tau; e++)
    {
      uintX_t open_alpha_e[ClenX];
      uintX_t epsilon_e[ClenX];

      /* Obtain challenge epsilon_e */
      batch_sample (epsilon_e, q, C, expand_ctx);

      /* Use variables a[e][i] to store shares of alpha_e */
      /* Prepare computation of open_alpha_e, init with 0 */
      batch_clear (open_alpha_e, q, C);
      for (i = 0; i < N; i++)
        {
          /* Compute alpha[e][i] = x[e][i] * epsilon_e + a[e][i] */
          batch_mul (x[e][i], epsilon_e, q, C);
          batch_add (a[e][i], x[e][i], q, C);

          /* Iterated to obtain open_alpha_e = sum(alpha[e][i]) */
          batch_add (open_alpha_e, a[e][i], q, C);

          /* Add alpha[e][i] in sigma2 */
          H2_update (hash_ctx, (void *) a[e][i], Clen);
        }

      /* Use variables y[e][i] to store v[e][i] values */
      for (i = 0; i < N; i++)
        {
          /* Compute v[e][i] = y[e][i] * open_alpha_e */
          /*                 - z[e][i] * epsilon_e - c[e][i] */
          batch_mul (y[e][i], open_alpha_e, q, C);
          batch_mul (z[e][i], epsilon_e, q, C);
          batch_sub (y[e][i], z[e][i], q, C);
          batch_sub (y[e][i], c[e][i], q, C);
        }
      for (i = 0; i < N; i++)
        {
          /* Add v[e][i] in sigma2 */
          H2_update (hash_ctx, (void *) y[e][i], Clen);
        }
    }
  /* Finalize computation of h2 using sigma2 */
  H2_final (h2, hash_ctx, lambda);


  /* Phase 4: Challenging the views of the MPC protocol */

  /* Prepare expansion of ibar_e challenges */
  expand_init (expand_ctx, h2, lambda >> 2);


  /* Phase 5: Opening the views of the MPC and checking protocols */

  for (e = 0; e < tau; e++)
    {
      uint16_t ibar_e;

      /* Obtain challenge ibar_e */
      ibar_e = index_sample (N, lN, expand_ctx);
      /* Compute path to recover seed[e][i], for i != ibar_e in signature */
      get_path (sig_path[e], root[e], salt, e, ibar_e, N, lambda);
      /* Copy com[e][ibar_e] in signature */
      memcpy (sig_com[e], com[e][ibar_e], lambda >> 2);
      /* Export alpha[e][ibar_e] in signature */
      batch_export (sigma, a[e][ibar_e], q, C, offset);
      offset += C;
    }

  return 0;
}

int
verify (const uint8_t *sig, const uint8_t *msg, uint64_t msglen,
        const uint8_t *pk, const params_t *params)
{
  int e, i;
  int offset;

#ifndef SEC_LEVEL
  const int lambda = params->lambda;
#endif
#ifndef NB_ITERATIONS
  const int tau = params->tau;
#endif
#ifndef NB_PARTIES
  const int N = params->N;
  const int lN = ilog2 (N);
#endif
#ifndef FIELD_SIZE
  const int q = params->q;
#endif
#ifndef NB_VARIABLES
  const int n = params->n;
#endif
#ifndef NB_EQUATIONS
  const int m = params->m;
#endif
#ifndef DEGREE
  const int d = params->d;
#endif
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklen = batch_getlen (q, n);
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
#if !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int C = (d - 1) * m;
#ifdef HIGH_DEGREE
  const int Cm = (d - 2) * m;
#endif
#endif
#if !defined(FIELD_SIZE) || !defined(DEGREE) || !defined(NB_EQUATIONS)
  const int Clen = batch_getlen (q, C);
  const int ClenX = batch_getlenX (q, C);
#ifdef HIGH_DEGREE
  const int Cmlen = batch_getlen (q, Cm);
#endif
#endif

  const uint8_t *const seedF = pk;
  const uint8_t *const pk_data = pk + (lambda >> 3);

  const uint8_t *const salt = sig;
  const uint8_t *const h1 = sig + (lambda >> 2);
  const uint8_t *const h2 = sig + (lambda >> 2) + (lambda >> 2);
  const uint8_t (*const sig_path)[lN * (lambda >> 3) + (lambda >> 2)]
    = (void *) (sig + 3 * (lambda >> 2));
  const uint8_t (*const sig_com)[lN * (lambda >> 3) + (lambda >> 2)]
    = (void *) (sig + 3 * (lambda >> 2) + lN * (lambda >> 3));
  const uint8_t *const sigma
    = sig + 3 * (lambda >> 2) + tau * (lN * (lambda >> 3) + (lambda >> 2));

  h_ctx_t hash1_ctx, hash2_ctx;
  h_ctx_t expand_ctx, expand2_ctx;

  uint8_t h1p[lambda >> 2];
  uint8_t h2p[lambda >> 2];

  uintX_t pk_t[pklenX];
  uintX_t f[(d + 1) * (pklenX + sklenX * m)];

  BATCH_PARAMS (q, n, m, d);

  /* Recover t from pk */
  batch_import (pk_t, pk_data, q, m, 0);

  /* Expand seedF to obtain the system */
  expand_init (expand_ctx, seedF, lambda >> 3);
  circuit_sample (f, q, n, m, d, expand_ctx);

  /* Prepare expansion of epsilon_e challenges for the checking protocol */
  expand_init (expand_ctx, h1, lambda >> 2);

  /* Prepare expansion of ibar_e challenges for the views of the MPC protocol */
  expand_init (expand2_ctx, h2, lambda >> 2);

  /* initialize offset to read in the signature */
  offset = 0;

  /* Initialize sigma1 with salt and msg */
  H1_init (hash1_ctx, salt, msg, msglen, lambda);
  /* Initialize sigma2 with salt and h1 */
  H2_init (hash2_ctx, salt, h1, lambda);
  for (e = 0; e < tau; e++)
    {
      uint8_t com_e[N][lambda >> 2];

      uintX_t a_e[N][ClenX], c_e[N][ClenX];
      uintX_t x_e[N][ClenX], y_e[N][ClenX], z_e[N][ClenX];

      uintX_t epsilon_e[ClenX];
      uint16_t ibar_e;

      uintX_t s_e[N][sklenX];

      uintX_t delta_s_e[sklenX];
#ifdef HIGH_DEGREE
      uintX_t delta_z_e[ClenX];
#endif
      uintX_t delta_c_e[ClenX];

      uint8_t seed_e[N][lambda >> 3];

      uintX_t open_alpha_e[ClenX];

      /* Obtain challenge epsilon_e */
      batch_sample (epsilon_e, q, C, expand_ctx);
      /* Obtain challenge ibar_e */
      ibar_e = index_sample (N, lN, expand2_ctx);

      /* Recompute all the seeds from path (extracted from signature) */
      get_path_seeds ((void *) seed_e, sig_path[e], salt, e, ibar_e, N, lambda);

      /* Extract missing commitment com_e[ibar_e] from signature */
      memcpy (com_e[ibar_e], sig_com[e], lambda >> 2);

      /* Extract delta values from signature */
#ifdef HIGH_DEGREE
      offset = e * (n + Cm + C);
#else
      offset = e * (n + C);
#endif
      batch_import (delta_s_e, sigma, q, n, offset);
      offset += n;
#ifdef HIGH_DEGREE
      batch_import (delta_z_e, sigma, q, Cm, offset);
      offset += Cm;
#endif
      batch_import (delta_c_e, sigma, q, C, offset);
      offset += C;

      for (i = 0; i < N; i++)
        {
          if (i != ibar_e)
            {
              h_ctx_t tape_ei_ctx;

              /* Compute commitment on seed_e[i] */
              commit (com_e[i], salt, e, i, seed_e[i], lambda);

              /* Expand random from tape_ei */
              expandtape_init (tape_ei_ctx, salt, e, i, seed_e[i], lambda);
              batch_sample (s_e[i], q, n, tape_ei_ctx);
#ifdef HIGH_DEGREE
              batch_sample (z_e[i], q, Cm, tape_ei_ctx);
#endif
              batch_sample (a_e[i], q, C, tape_ei_ctx);
              batch_sample (c_e[i], q, C, tape_ei_ctx);
              if (i == 0)
                {
                  /* Correct s_e[0], z_e[0] and c_e[0] */
                  batch_add (s_e[0], delta_s_e, q, n);
#ifdef HIGH_DEGREE
                  batch_add (z_e[0], delta_z_e, q, Cm);
#endif
                  batch_add (c_e[0], delta_c_e, q, C);
                }

              /* Compute the shares x_e[i], y_e[i], z_e[i] */
              /* involved in multiplications */
              linear_circuit (x_e[i], y_e[i], z_e[i], s_e[i], pk_t, i,
                              q, n, m, d, f);
            }
          /* Add com_e[i] to sigma1 */
          H1_update (hash1_ctx, com_e[i], lambda >> 2);
        }

      /* Add delta values to sigma1 */
      H1_update (hash1_ctx, (void *) delta_s_e, sklen);
#ifdef HIGH_DEGREE
      H1_update (hash1_ctx, (void *) delta_z_e, Cmlen);
#endif
      H1_update (hash1_ctx, (void *) delta_c_e, Clen);

      /* Use variables a_e[i] to store shares of alpha_e */
      /* Extract alpha_e[ibar_e] from signature */
#ifdef HIGH_DEGREE
      offset = tau * (n + Cm + C) + e * C;
#else
      offset = tau * (n + C) + e * C;
#endif
      batch_import (a_e[ibar_e], sigma, q, C, offset);
      offset += C;

      /* Prepare computation of open_alpha_e, init with alpha_e[ibar_e] */
      batch_copy (open_alpha_e, a_e[ibar_e], q, C);
      for (i = 0; i < N; i++)
        {
          if (i != ibar_e)
            {
              /* Compute alpha_e[i] = x_e[i] * epsilon_e + a_e[i] */
              batch_mul (x_e[i], epsilon_e, q, C);
              batch_add (a_e[i], x_e[i], q, C);

              /* Iterated to obtain open_alpha_e = sum(alpha_e[i]) */
              batch_add (open_alpha_e, a_e[i], q, C);
            }

          /* Add alpha_e[i] in sigma2 */
          H2_update (hash2_ctx, (void *) a_e[i], Clen);
        }

      /* Use variables y_e[i] to store v_e[i] values */
      /* Prepare computation of v_e[ibar_e], init with 0 */
      batch_clear (y_e[ibar_e], q, C);
      for (i = 0; i < N; i++)
        {
          if (i != ibar_e)
            {
              /* Compute v_e[i] = y_e[i] * open_alpha_e */
              /*                - z_e[i] * epsilon_e - c_e[i] */
              batch_mul (y_e[i], open_alpha_e, q, C);
              batch_mul (z_e[i], epsilon_e, q, C);
              batch_sub (y_e[i], z_e[i], q, C);
              batch_sub (y_e[i], c_e[i], q, C);

              /* Iterated to obtain v_e[ibar_e] = 0 - sum(v_e[i]) */
              /* for i != ibar[e] */
              batch_sub (y_e[ibar_e], y_e[i], q, C);
            }
        }
      for (i = 0; i < N; i++)
        {
          /* Add v_e[i] in sigma2 */
          H2_update (hash2_ctx, (void *) y_e[i], Clen);
        }
    }
  /* Finalize computation of h1p using sigma1 */
  H1_final (h1p, hash1_ctx, lambda);
  /* Finalize computation of h2p using sigma2 */
  H2_final (h2p, hash2_ctx, lambda);

  /* Return comparison result of hash values */
  return memcmp (h1, h1p, lambda >> 3) | memcmp (h2, h2p, lambda >> 3);
}
