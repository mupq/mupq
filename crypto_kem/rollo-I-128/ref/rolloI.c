/**
 * \file kem.c
 * \brief Implementation of api.h
 */

#include "api.h"
#include "rbc_qre.h"
#include "rbc_vec.h"
#include "sha2.h"
#include "fips202.h"
#include "parameters.h"
#include "string.h"
#include "lrpc.h"
#include "types.h"
#include "parsing.h"


int crypto_kem_keypair(uint8_t* pk, uint8_t* sk) {
  rolloI_secretKey skTmp;
  rolloI_publicKey pkTmp;

  rbc_qre invX;

  rbc_field_init();
  rbc_qre_init_modulus(ROLLOI_PARAM_N);

  uint8_t sk_seed[SEEDEXPANDER_SEED_BYTES];
  randombytes(sk_seed, SEEDEXPANDER_SEED_BYTES);

  rolloI_secret_key_from_string(&skTmp, sk_seed);

  rbc_qre_init(&invX);
  rbc_qre_inv(invX, skTmp.x);

  rbc_qre_init(&(pkTmp.h));
  rbc_qre_mul(pkTmp.h, invX, skTmp.y);

  rolloI_secret_key_to_string(sk, sk_seed);
  rolloI_public_key_to_string(pk, &pkTmp);

  rbc_qre_clear(invX);

  rbc_vspace_clear(skTmp.F);
  rbc_qre_clear(skTmp.x);
  rbc_qre_clear(skTmp.y);
  rbc_qre_clear(pkTmp.h);
  rbc_qre_clear_modulus();

  return 0;
}

int crypto_kem_enc(uint8_t* ct, uint8_t* ss, const uint8_t* pk) {
  rolloI_publicKey pkTmp;
  rolloI_ciphertext ctTmp;

  rbc_vspace E;
  rbc_qre E1, E2;

  rbc_field_init();
  rbc_qre_init_modulus(ROLLOI_PARAM_N);

  rolloI_public_key_from_string(&pkTmp, pk);

  rbc_vspace_init(&E, ROLLOI_PARAM_R);

  //Support
  rbc_vspace_set_random_full_rank2(E, ROLLOI_PARAM_R);

  rbc_qre_init(&E1);
  rbc_qre_init(&E2);
  rbc_qre_init(&(ctTmp.syndrom));

  //Random error vectors
  rbc_qre_set_random_pair_from_support2(E1, E2, E, ROLLOI_PARAM_R);

  rbc_qre_mul(ctTmp.syndrom, E2, pkTmp.h);
  rbc_qre_add(ctTmp.syndrom, ctTmp.syndrom, E1);

  rolloI_rolloI_ciphertext_to_string(ct, &ctTmp);

  rbc_vec_echelonize(E, ROLLOI_PARAM_R);

  uint8_t support[ROLLOI_RBC_VEC_R_BYTES];
  rbc_vec_to_string(support, E, ROLLOI_PARAM_R);
  sha512(ss, support, ROLLOI_RBC_VEC_R_BYTES);

  rbc_vspace_clear(E);
  rbc_qre_clear(E1);
  rbc_qre_clear(E2);
  rbc_qre_clear(pkTmp.h);
  rbc_qre_clear(ctTmp.syndrom);
  rbc_qre_clear_modulus();

  return 0;
}

int crypto_kem_dec(uint8_t* ss, const uint8_t* ct, const uint8_t* sk) {
  rolloI_secretKey skTmp;
  rolloI_ciphertext ctTmp;

  rbc_qre xc;
  rbc_vspace E;
  uint8_t support[ROLLOI_RBC_VEC_R_BYTES];

  rbc_field_init();
  rbc_qre_init_modulus(ROLLOI_PARAM_N);

  rolloI_secret_key_from_string(&skTmp, sk);
  rolloI_rolloI_ciphertext_from_string(&ctTmp, ct);

  rbc_qre_init(&xc);
  rbc_qre_mul(xc, skTmp.x, ctTmp.syndrom);

  uint32_t dimE = 0;

  rbc_vspace_init(&E, ROLLOI_PARAM_N);

  dimE = rbc_lrpc_RSR(E, ROLLOI_PARAM_R, skTmp.F, ROLLOI_PARAM_D, xc, ROLLOI_PARAM_N);

  if(dimE != 0) {
    rbc_vec_to_string(support, E, ROLLOI_PARAM_R);
    sha512(ss, support, ROLLOI_RBC_VEC_R_BYTES);
  } else {
    memset(ss, 0, ROLLOI_SHARED_SECRET_BYTES);
  }

  rbc_vspace_clear(E);
  rbc_qre_clear(xc);
  rbc_vspace_clear(skTmp.F);
  rbc_qre_clear(skTmp.x);
  rbc_qre_clear(skTmp.y);
  rbc_qre_clear(ctTmp.syndrom);
  rbc_qre_clear_modulus();

  if(dimE != ROLLOI_PARAM_R) return 1;

  return 0;
}
