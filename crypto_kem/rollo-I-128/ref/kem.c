/**
 * \file kem.c
 * \brief Implementation of api.h
 */

#include "api.h"
#include "ffi_qre.h"
#include "sha2.h"
#include "parameters.h"
#include "string.h"
#include "rsr_algorithm.h"
#include "rolloI_types.h"
#include "parsing.h"
#include "randombytes.h"

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk) {
  secretKey skTmp;
  publicKey pkTmp;

  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  unsigned char sk_seed[SEEDEXPANDER_SEED_BYTES];
  randombytes(sk_seed, SEEDEXPANDER_SEED_BYTES);

  rolloI_secret_key_from_string(&skTmp, sk_seed);

  ffi_qre invX;
  ffi_qre_init(&invX);
  ffi_qre_inv(invX, skTmp.x);

  ffi_qre_init(&(pkTmp.h));
  ffi_qre_mul(pkTmp.h, invX, skTmp.y);

  rolloI_secret_key_to_string(sk, sk_seed);
  rolloI_public_key_to_string(pk, &pkTmp);

  ffi_qre_clear(invX);
  ffi_vspace_clear(skTmp.F, PARAM_D);
  ffi_qre_clear(skTmp.x);
  ffi_qre_clear(skTmp.y);
  ffi_qre_clear(pkTmp.h);
  ffi_qre_clear_modulus();

  return 0;
}

int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk) {
  publicKey pkTmp;
  ciphertext ctTmp;

  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  rolloI_public_key_from_string(&pkTmp, pk);

  ffi_vspace E;
  ffi_vspace_init(&E, PARAM_R);

  //Support
  ffi_vspace_set_random_full_rank2(E, PARAM_R);

  ffi_qre E1, E2;

  ffi_qre_init(&E1);
  ffi_qre_init(&E2);
  ffi_qre_init(&(ctTmp.syndrom));

  //Random error vectors
  ffi_qre_set_random_from_support2(E1, E, PARAM_R);
  ffi_qre_set_random_from_support2(E2, E, PARAM_R);

  ffi_qre_mul(ctTmp.syndrom, E2, pkTmp.h);
  ffi_qre_add(ctTmp.syndrom, ctTmp.syndrom, E1);

  rolloI_ciphertext_to_string(ct, &ctTmp);

  ffi_vec_echelonize(E, PARAM_R);

  unsigned char support[FFI_VEC_R_BYTES];
  ffi_vec_to_string_compact(support, E, PARAM_R);
  sha512(ss, support, FFI_VEC_R_BYTES);

  ffi_vspace_clear(E, PARAM_R);
  ffi_qre_clear(E1);
  ffi_qre_clear(E2);
  ffi_qre_clear(pkTmp.h);
  ffi_qre_clear(ctTmp.syndrom);
  ffi_qre_clear_modulus();

  return 0;
}

int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk) {
  secretKey skTmp;
  ciphertext ctTmp;

  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  rolloI_secret_key_from_string(&skTmp, sk);
  rolloI_ciphertext_from_string(&ctTmp, ct);

  ffi_qre xc;
  ffi_qre_init(&xc);

  ffi_qre_mul(xc, skTmp.x, ctTmp.syndrom);

  ffi_vspace E;
  unsigned int dimE = 0;

  ffi_vspace_init(&E, PARAM_N);

  dimE = rank_support_recoverer(E, PARAM_R, skTmp.F, PARAM_D, xc, PARAM_N);

  if(dimE != 0) {
    unsigned char support[FFI_VEC_R_BYTES];
    ffi_vec_to_string_compact(support, E, PARAM_R);
    sha512(ss, support, FFI_VEC_R_BYTES);
  } else {
    memset(ss, 0, SHARED_SECRET_BYTES);
  }

  ffi_vspace_clear(E, PARAM_R);
  ffi_qre_clear(xc);
  ffi_vspace_clear(skTmp.F, PARAM_D);
  ffi_qre_clear(skTmp.x);
  ffi_qre_clear(skTmp.y);
  ffi_qre_clear(ctTmp.syndrom);
  ffi_qre_clear_modulus();

  if(dimE != PARAM_R) return 1;

  return 0;
}
