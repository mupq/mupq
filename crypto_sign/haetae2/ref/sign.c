// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Team HAETAE
// Copyright (c) 2026 YoungBeom Kim (pqm4 adaptation)

#include "api.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "polyfix.h"
#include "polymat.h"
#include "polyvec.h"
#include "randombytes.h"
#include "symmetric.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/*************************************************
 * Name:        crypto_sign_keypair_internal
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *vk: pointer to output public key (allocated
 *              array of HAETAE_CRYPTO_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key (allocated
 *              array of HAETAE_CRYPTO_SECRETKEYBYTES bytes)
 *              - uint8_t seed : input seed
 *
 * Returns 0 (success)
 *
 * Specification: Implements @[Algorithm 54, KeyGenInternal]
 **************************************************/

int crypto_sign_keypair_internal(uint8_t *vk, uint8_t *sk,
                                 uint8_t seed[HAETAE_SEEDBYTES]) {
  uint8_t seedbuf[2 * HAETAE_SEEDBYTES + HAETAE_CRHBYTES] = {0};
  uint16_t counter = 0;
  const uint8_t *rhoprime, *sigma, *key;
  polyvecm A[HAETAE_K], s1, s1hat;
  polyveck b, s2;
#if (HAETAE_MODE == HAETAE_MODE2) || (HAETAE_MODE == HAETAE_MODE3)
  polyveck a, b0;
#elif HAETAE_MODE == HAETAE_MODE5
  polyveck s2hat;
#endif
  xof256_state state;

  // Get seed_A
  memcpy(seedbuf, seed, HAETAE_SEEDBYTES);

  // Sample seeds with entropy \rho
  xof256_absorb_once(&state, seedbuf, HAETAE_SEEDBYTES);
  xof256_squeeze(seedbuf, 2 * HAETAE_SEEDBYTES + HAETAE_CRHBYTES, &state);

  rhoprime = seedbuf;
  sigma = rhoprime + HAETAE_SEEDBYTES;
  key = sigma + HAETAE_CRHBYTES;

  // Expand Matrix A0 and vector a
  polymatkm_expand_matA(A, rhoprime);

#if (HAETAE_MODE == HAETAE_MODE2) || (HAETAE_MODE == HAETAE_MODE3)
  /**********************************************
   * If there is rounding (D > 0), we need another polyveck a.
   * Then, b = a + A0 * s1 + s2 and the lower D bits are
   * rounded from b. The lower D bits are subsequently
   * subtracted from s2.
   **********************************************/
  polyveck_expand_vecA(&a, rhoprime);

reject:
  // Sample secret vectors s1 and s2
  polyvecmk_expand_S(&s1, &s2, sigma, counter);
  counter += HAETAE_M + HAETAE_K;

  // b = a + A0 * s1 + s2 mod q
  s1hat = s1;
  polyvecm_ntt(&s1hat);
  polymatkm_pointwise_montgomery(&b, A, &s1hat);
  polyveck_invntt_tomont(&b);
  polyveck_add(&b, &b, &s2);
  polyveck_add(&b, &b, &a);
  polyveck_freeze(&b);

  // round off D bits
  polyveck_decompose_vk(&b0, &b);
  polyveck_sub(&s2, &s2, &b0);

  int64_t squared_singular_value = polyvecmk_sk_singular_value(&s1, &s2);
  if (squared_singular_value > HAETAE_GAMMA * HAETAE_GAMMA * HAETAE_N) {
    goto reject;
  }
#elif HAETAE_MODE == HAETAE_MODE5
  /**********************************************
   * If there is no rounding (D == 0), we store
   * -2b directly in NTT domain into the public key.
   **********************************************/
reject:
  // Sample secret vectors s1 and s2
  polyvecmk_expand_S(&s1, &s2, sigma, counter);
  counter += HAETAE_M + HAETAE_K;
  int64_t squared_singular_value = polyvecmk_sk_singular_value(&s1, &s2);
  if (squared_singular_value > HAETAE_GAMMA * HAETAE_GAMMA * HAETAE_N) {
    goto reject;
  }

  // b = A0 * s1 + s2 mod q
  s1hat = s1;
  s2hat = s2;
  polyvecm_ntt(&s1hat);
  polyveck_ntt(&s2hat);
  polymatkm_pointwise_montgomery(&b, A, &s1hat);
  polyveck_frommont(&b);
  polyveck_add(&b, &b, &s2hat);
  polyveck_double_negate(&b);
  polyveck_caddq(&b); // directly compute and store NTT(-2b)
#endif

  pack_vk(vk, &b, rhoprime);
  pack_sk(sk, vk, &s1, &s2, key);

  return 0;
}

/*************************************************
 * Name:        crypto_sign_keypair_internal
 *
 * Description: Generates public and private key.
 *
 * Arguments:   - uint8_t *vk: pointer to output public key (allocated
 *              array of HAETAE_CRYPTO_PUBLICKEYBYTES bytes)
 *              - uint8_t *sk: pointer to output private key (allocated
 *              array of HAETAE_CRYPTO_SECRETKEYBYTES bytes)
 *
 * Returns 0 (success)
 *
 * Specification: Implements @[Algorithm 51, KeyGen]
 **************************************************/
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
  int ret = 0;
  uint8_t seed[HAETAE_SEEDBYTES] = {0};

  randombytes(seed, HAETAE_SEEDBYTES);

  ret = crypto_sign_keypair_internal(pk, sk, seed);

  return ret;
}

/*************************************************
 * Name:        crypto_sign_signature_internal
 *
 * Description: Computes signature. Internal API.
 *
 * Arguments:   - uint8_t *sig: pointer to output signature (of length
 *              HAETAE_CRYPTO_BYTES)
 *              - size_t *siglen: pointer to output length of signature
 *              - uint8_t *m: pointer to message to be signed
 *              - size_t mlen: length of message
 *              - uint8_t *pre: pointer to prefix string
 *              - size_t prelen: length of prefix string
 *              - uint8_t *rnd: pointer to random seed
 *              - uint8_t *sk: pointer to bit-packed secret key
 *
 * Returns 0 (success)
 *
 * Specification: Implements @[Algorithm 55, SignInternal]
 **************************************************/
int crypto_sign_signature_internal(uint8_t *sig, size_t *siglen,
                                   const uint8_t *m, size_t mlen,
                                   const uint8_t *pre, size_t prelen,
                                   const uint8_t rnd[HAETAE_SEEDBYTES],
                                   const uint8_t *sk) {

  uint8_t buf[HAETAE_POLYVECK_HIGHBITS_PACKEDBYTES + HAETAE_POLYC_PACKEDBYTES] =
      {0};
  uint8_t seedbuf[HAETAE_CRHBYTES] = {0}, key[HAETAE_SEEDBYTES] = {0};
  uint8_t mu[HAETAE_CRHBYTES] = {0};
  uint8_t b = 0; // one bit
  uint16_t counter = 0;
  uint64_t reject1, reject2;

  polyvecm s1;
  polyvecl A1[HAETAE_K], cs1;
  polyveck s2, cs2, highbits, Ay;
  polyfixvecl y1, z1, z1tmp;
  polyfixveck y2, z2, z2tmp;
  polyvecl z1rnd; // round of z1
  polyvecl hb_z1, lb_z1;
  polyveck z2rnd, h, htmp; // round of z2
  poly c, chat, z1rnd0, lsb;

  xof256_state state;

  unsigned int i;

  // Unpack secret key
  unpack_sk(A1, &s1, &s2, key, sk);

  xof256_absorb_thrice(&state, sk, HAETAE_CRYPTO_PUBLICKEYBYTES, pre, prelen, m,
                       mlen);
  xof256_squeeze(mu, HAETAE_CRHBYTES, &state);

  xof256_absorb_thrice(&state, key, HAETAE_SEEDBYTES, rnd, HAETAE_SEEDBYTES, mu,
                       HAETAE_CRHBYTES);
  xof256_squeeze(seedbuf, HAETAE_CRHBYTES, &state);

  polyvecm_ntt(&s1);
  polyveck_ntt(&s2);

reject:

  /*------------------ 1. Sample y1 and y2 from hyperball ------------------*/
  counter = polyfixveclk_sample_hyperball(&y1, &y2, &b, seedbuf, counter);

  /*------------------- 2. Compute a chanllenge c --------------------------*/
  // Round y1 and y2
  polyfixvecl_round(&z1rnd, &y1);
  polyfixveck_round(&z2rnd, &y2);

  // A * round(y) mod q = A1 * round(y1) + 2 * round(y2) mod q
  z1rnd0 = z1rnd.vec[0];
  polyvecl_ntt(&z1rnd);
  polymatkl_pointwise_montgomery(&Ay, A1, &z1rnd);
  polyveck_invntt_tomont(&Ay);
  polyveck_double(&z2rnd);
  polyveck_add(&Ay, &Ay, &z2rnd);

  // recover A * round(y) mod 2q
  polyveck_poly_fromcrt(&Ay, &Ay, &z1rnd0);
  polyveck_freeze2q(&Ay);

  // HighBits of (A * round(y) mod 2q)
  polyveck_highbits_hint(&highbits, &Ay);

  // LSB(round(y_0) * j)
  poly_lsb(&lsb, &z1rnd0);

  // Pack HighBits of A * round(y) mod 2q and LSB of round(y0)
  pack_vec_highbits(buf, &highbits);
  pack_poly_lsb(buf + HAETAE_POLYVECK_HIGHBITS_PACKEDBYTES, &lsb);

  // c = challenge(highbits, lsb, mu)
  poly_challenge(&c, buf, mu);

  /*------------------- 3. Compute z = y + (-1)^b c * s --------------------*/
  // cs = c * s = c * (si1 || s2)
  cs1.vec[0] = c;
  chat = c;
  poly_ntt(&chat);

  for (i = 1; i < HAETAE_L; ++i) {
    poly_pointwise_montgomery(&cs1.vec[i], &chat, &s1.vec[i - 1]);
    poly_invntt_tomont(&cs1.vec[i]);
  }
  polyveck_poly_pointwise_montgomery(&cs2, &s2, &chat);
  polyveck_invntt_tomont(&cs2);

  // z = y + (-1)^b cs = z1 + z2
  polyvecl_cneg(&cs1, b & 1);
  polyveck_cneg(&cs2, b & 1);
  polyfixvecl_add(&z1, &y1, &cs1);
  polyfixveck_add(&z2, &y2, &cs2);

  // reject if norm(z) >= B'
  reject1 = ((uint64_t)HAETAE_B1SQ * HAETAE_LN * HAETAE_LN -
             polyfixveclk_sqnorm2(&z1, &z2)) >>
            63;
  reject1 &= 1;

  polyfixvecl_double(&z1tmp, &z1);
  polyfixveck_double(&z2tmp, &z2);

  polyfixfixvecl_sub(&z1tmp, &z1tmp, &y1);
  polyfixfixveck_sub(&z2tmp, &z2tmp, &y2);

  // reject if norm(2z-y) < B and b' = 0
  reject2 = (polyfixveclk_sqnorm2(&z1tmp, &z2tmp) -
             (uint64_t)HAETAE_B0SQ * HAETAE_LN * HAETAE_LN) >>
            63;
  reject2 &= 1;
  reject2 &= (b & 0x2) >> 1;

  if (reject1 | reject2) {
    goto reject;
  }

  /*------------------- 4. Make a hint -------------------------------------*/
  // Round z1 and z2
  polyfixvecl_round(&z1rnd, &z1);
  polyfixveck_round(&z2rnd, &z2);

  // recover A1 * round(z1) - qcj mod 2q
  polyveck_double(&z2rnd);
  polyveck_sub(&htmp, &Ay, &z2rnd);
  polyveck_freeze2q(&htmp);

  // HighBits of (A * round(z) - qcj mod 2q) and (A1 * round(z1) - qcj mod 2q)
  polyveck_highbits_hint(&htmp, &htmp);
  polyveck_sub(&h, &highbits, &htmp);
  polyveck_caddDQ2ALPHA(&h);

  /*------------------ Decompose(z1) and Pack signature -------------------*/
  polyvecl_lowbits(&lb_z1, &z1rnd); // TODO do this in one function together!
  polyvecl_highbits(&hb_z1, &z1rnd);

  if (pack_sig(sig, &c, &lb_z1, &hb_z1,
               &h)) { // reject if signature is too big
    goto reject;
  }
  *siglen = HAETAE_CRYPTO_BYTES;

  return 0;
}

/*************************************************
 * Name:        crypto_sign_signature_ctx
 *
 * Description: Computes signature.
 *
 * Arguments:   - uint8_t *sig:   pointer to output signature (of length
 *                                HAETAE_CRYPTO_BYTES)
 *              - size_t *siglen: pointer to output length of signature
 *              - uint8_t *m:     pointer to message to be signed
 *              - size_t mlen:    length of message
 *              - uint8_t *ctx:   pointer of contex string
 *              - size_t ctxlen:  length of contex string
 *              - uint8_t *sk:    pointer to bit-packed secret key
 *
 * Returns 0 (success)
 *
 * Specification: Implements @[Algorithm 52, Sign]
 **************************************************/
int crypto_sign_signature_ctx(uint8_t *sig, size_t *siglen, const uint8_t *m,
                              size_t mlen, const uint8_t *ctx, size_t ctxlen,
                              const uint8_t *sk) {
  uint8_t pre[256];
  uint8_t rnd[HAETAE_SEEDBYTES];

  if (ctxlen > 255)
    return -1;

  /* Prepare pre = (ctxlen, ctx) */
  pre[0] = ctxlen;
  memcpy(pre + 1, ctx, ctxlen);

  randombytes(rnd, HAETAE_SEEDBYTES);

  return crypto_sign_signature_internal(sig, siglen, m, mlen, pre, 1 + ctxlen,
                                        rnd, sk);
}

/*************************************************
 * Name:        crypto_sign_ctx
 *
 * Description: Compute signed message.
 *
 * Arguments:   - uint8_t *sm: pointer to output signed message (allocated
 *                             array with HAETAE_CRYPTO_BYTES + mlen bytes),
 *                             can be equal to m
 *              - size_t *smlen: pointer to output length of signed
 *                               message
 *              - const uint8_t *m: pointer to message to be signed
 *              - size_t mlen: length of message
 *              - const uint8_t *ctx: pointer to context string
 *              - size_t ctxlen: length of context string
 *              - const uint8_t *sk: pointer to bit-packed secret key
 *
 * Returns 0 (success) or -1 (context string too long)
 **************************************************/
int crypto_sign_ctx(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                    const uint8_t *ctx, size_t ctxlen, const uint8_t *sk) {
  size_t i;

  for (i = 0; i < mlen; ++i)
    sm[HAETAE_CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  int ret = crypto_sign_signature_ctx(sm, smlen, sm + HAETAE_CRYPTO_BYTES, mlen,
                                      ctx, ctxlen, sk);
  *smlen += mlen;
  return ret;
}

/*************************************************
 * Name:        crypto_sign_verify_internal
 *
 * Description: Verifies signature.
 *
 * Arguments:   - uint8_t *m: pointer to input signature
 *              - size_t siglen: length of signature
 *              - const uint8_t *m: pointer to message
 *              - size_t mlen: length of message
 *              - const uint8_t *pre: pointer to prefix string
 *              - size_t prelen: length of prefix string
 *              - const uint8_t *vk: pointer to bit-packed public key
 *
 * Returns 0 if signature could be verified correctly and -1 otherwise
 *
 * Specification: Implements @[Algorithm 56, VerifyInternal]
 **************************************************/
int crypto_sign_verify_internal(const uint8_t *sig, size_t siglen,
                                const uint8_t *m, size_t mlen,
                                const uint8_t *pre, size_t prelen,
                                const uint8_t *vk) {
  unsigned int i;
  uint8_t buf[HAETAE_POLYVECK_HIGHBITS_PACKEDBYTES + HAETAE_POLYC_PACKEDBYTES] =
      {0};
  uint8_t mu[HAETAE_SEEDBYTES];
  uint64_t sqnorm2;
  polyvecl A1[HAETAE_K], highz, lowz, z1;
  polyveck highbits, h, z2, w;
  poly c, cprime, wprime;

  xof256_state state;

  // Check signature length
  if (siglen != HAETAE_CRYPTO_BYTES) {
    return -1;
  }
  /*------------------- 1. Recover A1 --------------------------------------*/
  // A1 = (-2b + qj || 2 * A0)
  // Unpack public key
  unpack_vk(A1, vk);

  // Unpack signature and Check conditions -- A1 is used only as intermediate
  // buffer for the low bits
  if (unpack_sig(&c, &lowz, &highz, &h, sig)) {
    return -1;
  }

  // Compose z1 out of HighBits(z1) and LowBits(z1)
  for (i = 0; i < HAETAE_L; ++i) {
    poly_compose(&z1.vec[i], &highz.vec[i], &lowz.vec[i]);
  }

  /*------------------- 2. Compute \tilde{z}_2 -----------------------------*/
  // compute squared norm of z1 and w' before NTT
  sqnorm2 = polyvecl_sqnorm2(&z1);
  poly_sub(&wprime, &z1.vec[0], &c);
  poly_lsb(&wprime, &wprime);

  // A1 * round(z1) - qcj mod q
  polyvecl_ntt(&z1);
  polymatkl_pointwise_montgomery(&highbits, A1, &z1);
  polyveck_invntt_tomont(&highbits);

  // recover A1 * round(z1) - qcj mod 2q
  polyveck_poly_fromcrt(&highbits, &highbits, &wprime);
  polyveck_freeze2q(&highbits);

  // recover w1
  polyveck_highbits_hint(&w, &highbits);
  polyveck_add(&w, &w, &h);
  polyveck_csubDQ2ALPHA(&w);

  // recover \tilde{z}_2 mod q
  polyveck_mul_alpha(&z2, &w);
  polyveck_sub(&z2, &z2, &highbits);
  poly_add(&z2.vec[0], &z2.vec[0], &wprime);
  polyveck_reduce2q(&z2);
  polyveck_div2(&z2);

  // check final norm of \tilde{z}
  if (sqnorm2 + polyveck_sqnorm2(&z2) > HAETAE_B2SQ) {
    return -1;
  }

  /*------------------- 3. Compute c_seed' and Compare ---------------------*/

  // Pack highBits(A * round(z) - qcj mod 2q) and h'
  pack_vec_highbits(buf, &w);
  pack_poly_lsb(buf + HAETAE_POLYVECK_HIGHBITS_PACKEDBYTES, &wprime);

  xof256_absorb_thrice(&state, vk, HAETAE_CRYPTO_PUBLICKEYBYTES, pre, prelen, m,
                       mlen);
  xof256_squeeze(mu, HAETAE_SEEDBYTES, &state);
  // c_seed = H(HighBits(A * y mod 2q), LSB(round(y0) * j), HAETAE_M)
  poly_challenge(&cprime, buf, mu);

  for (i = 0; i < HAETAE_N; ++i) {
    if (c.coeffs[i] != cprime.coeffs[i]) {
      return -1;
    }
  }
  return 0;
}

/*************************************************
 * Name:        crypto_sign_verify_ctx
 *
 * Description: Verifies signature.
 *
 * Arguments:   - uint8_t *m: pointer to input signature
 *              - size_t siglen: length of signature
 *              - const uint8_t *m: pointer to message
 *              - size_t mlen: length of message
 *              - const uint8_t *ctx: pointer to context string
 *              - size_t ctxlen: length of context string
 *              - const uint8_t *vk: pointer to bit-packed public key
 *
 * Returns 0 if signature could be verified correctly and -1 otherwise
 *
 * Specification: Implements @[Algorithm 53, Verify]
 **************************************************/
int crypto_sign_verify_ctx(const uint8_t *sig, size_t siglen, const uint8_t *m,
                           size_t mlen, const uint8_t *ctx, size_t ctxlen,
                           const uint8_t *vk) {
  size_t i;
  uint8_t pre[257];

  if (ctxlen > 255)
    return -1;

  pre[0] = ctxlen;
  memcpy(pre + 1, ctx, ctxlen);

  return crypto_sign_verify_internal(sig, siglen, m, mlen, pre, 1 + ctxlen, vk);
}

/*************************************************
 * Name:        crypto_sign_open_ctx
 *
 * Description: Verify signed message.
 *
 * Arguments:   - uint8_t *m: pointer to output message (allocated
 *                            array with smlen bytes), can be equal to sm
 *              - size_t *mlen: pointer to output length of message
 *              - const uint8_t *sm: pointer to signed message
 *              - size_t smlen: length of signed message
 *              - const uint8_t *vk: pointer to bit-packed public key
 *
 * Returns 0 if signed message could be verified correctly and -1 otherwise
 **************************************************/
int crypto_sign_open_ctx(uint8_t *m, size_t *mlen, const uint8_t *sm,
                         size_t smlen, const uint8_t *ctx, size_t ctxlen,
                         const uint8_t *vk) {
  size_t i;

  if (smlen < HAETAE_CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - HAETAE_CRYPTO_BYTES;
  if (crypto_sign_verify_ctx(sm, HAETAE_CRYPTO_BYTES, sm + HAETAE_CRYPTO_BYTES,
                             *mlen, ctx, ctxlen, vk))
    goto badsig;
  else {
    /* All good, copy msg, return 0 */
    for (i = 0; i < *mlen; ++i)
      m[i] = sm[HAETAE_CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for (i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
