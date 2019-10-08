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
#include "parsing.h"
#include "randombytes.h"

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk) {
  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  unsigned char sk_seed[SEEDEXPANDER_SEED_BYTES];
  randombytes(sk_seed, SEEDEXPANDER_SEED_BYTES);
  AES_XOF_struct* sk_seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
  seedexpander_init(sk_seedexpander, sk_seed, sk_seed + 32, SEEDEXPANDER_MAX_LENGTH);

  unsigned char pk_seed[SEEDEXPANDER_SEED_BYTES];
  randombytes(pk_seed, SEEDEXPANDER_SEED_BYTES);
  AES_XOF_struct* pk_seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
  seedexpander_init(pk_seedexpander, pk_seed, pk_seed + 32, SEEDEXPANDER_MAX_LENGTH);

  //Secret key
  ffi_vspace F;
  ffi_vspace_init(&F, PARAM_W);
  ffi_vspace_set_random_full_rank_with_one(F, PARAM_W, sk_seedexpander);

  ffi_qre x, y;
  ffi_qre_init(&x);
  ffi_qre_init(&y);

  ffi_qre_set_random_from_support(x, F, PARAM_W, sk_seedexpander);
  ffi_qre_set_random_from_support(y, F, PARAM_W, sk_seedexpander);

  //Public key
  ffi_qre h, s;
  ffi_qre_init(&h);
  ffi_qre_init(&s);
  ffi_vec_set_random(h->v, PARAM_N, pk_seedexpander);
  ffi_qre_mul(s, h, y);
  ffi_qre_add(s, s, x);

  rolloIII_secret_key_to_string(sk, sk_seed);
  rolloIII_public_key_to_string(pk, s, pk_seed);

  ffi_vspace_clear(F, PARAM_W);
  ffi_qre_clear(x);
  ffi_qre_clear(y);
  ffi_qre_clear(h);
  ffi_qre_clear(s);
  ffi_qre_clear_modulus();

  free(pk_seedexpander);
  free(sk_seedexpander);

  return 0;
}

int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk) {
  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  ffi_qre h, s;
  rolloIII_public_key_from_string(&h, &s, pk);

  ffi_vspace E;
  ffi_vspace_init(&E, PARAM_W_R);

  //Support
  ffi_vspace_set_random_full_rank2(E, PARAM_W_R);

  ffi_qre r1, r2, er;

  ffi_qre_init(&r1);
  ffi_qre_init(&r2);
  ffi_qre_init(&er);

  //Random error vectors
  ffi_qre_set_random_from_support2(r1, E, PARAM_W_R);
  ffi_qre_set_random_from_support2(r2, E, PARAM_W_R);
  ffi_qre_set_random_from_support2(er, E, PARAM_W_R);

  ffi_qre sr;
  ffi_qre_init(&sr);
  ffi_qre_mul(sr, h, r2);
  ffi_qre_add(sr, sr, r1);

  ffi_qre se;
  ffi_qre_init(&se);
  ffi_qre_mul(se, s, r2);
  ffi_qre_add(se, se, er);

  rolloIII_ciphertext_to_string(ct, sr, se);

  ffi_vec_echelonize(E, PARAM_W_R);

  unsigned char support[FFI_VEC_R_BYTES];
  ffi_vec_to_string_compact(support, E, PARAM_W_R);
  sha512(ss, support, FFI_VEC_R_BYTES);

  ffi_qre_clear(h);
  ffi_qre_clear(s);
  ffi_vspace_clear(E, PARAM_W_R);
  ffi_qre_clear(r1);
  ffi_qre_clear(r2);
  ffi_qre_clear(er);
  ffi_qre_clear(sr);
  ffi_qre_clear(se);
  ffi_qre_clear_modulus();

  return 0;
}

int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk) {
  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  ffi_qre sr, se;
  rolloIII_ciphertext_from_string(&sr, &se, ct);

  ffi_vspace F;
  ffi_qre x, y;
  rolloIII_secret_key_from_string(&x, &y, &F, sk);

  ffi_qre ec;
  ffi_qre_init(&ec);
  ffi_qre_mul(ec, y, sr);
  ffi_qre_add(ec, ec, se);

  ffi_vspace E;
  unsigned int dimE = 0;

  ffi_vspace_init(&E, PARAM_N);

  dimE = rank_support_recoverer(E, PARAM_W_R, F, PARAM_W, ec, PARAM_N);

  if(dimE != 0) {
    unsigned char support[FFI_VEC_R_BYTES];
    ffi_vec_to_string_compact(support, E, PARAM_W_R);
    sha512(ss, support, FFI_VEC_R_BYTES);
  } else {
    memset(ss, 0, SHARED_SECRET_BYTES);
  }

  ffi_vspace_clear(E, PARAM_N);
  ffi_qre_clear(sr);
  ffi_qre_clear(se);
  ffi_vspace_clear(F, PARAM_W);
  ffi_qre_clear(x);
  ffi_qre_clear(y);
  ffi_qre_clear(ec);
  ffi_qre_clear_modulus();

  if(dimE != PARAM_W_R) return 1;

  return 0;
}
