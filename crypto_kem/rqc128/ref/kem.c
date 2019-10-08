/**
 * \file kem.c
 * \brief Implementation of api.h
 */

#include "api.h"
#include "ffi_qre.h"
#include "sha2.h"
#include "parameters.h"
#include "parsing.h"
#include "rqc.h"
#include "string.h"

/**
 * \fn int crypto_kem_keypair(unsigned char* pk, unsigned char* sk)
 * \brief Keygen of the RQC_KEM IND-CCA2 scheme
 *
 * The public key is composed of the syndrom <b>s</b> as well as the seed used to generate vectors <b>g</b> and <b>h</b>.
 *
 * The secret key is composed of the seed used to generate the vectors <b>x</b> and <b>y</b>.
 * As a technicality, the public key is appended to the secret key in order to respect the NIST API.
 *
 * \param[out] pk String containing the public key
 * \param[out] sk String containing the secret key
 * \return 0 if keygen is sucessfull
 */
int crypto_kem_keypair(unsigned char* pk, unsigned char* sk) {
  rqc_pke_keygen(pk, sk);
  return 0;
}



/**
 * \fn int crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk)
 * \brief Encapsulation of the RQC_KEM IND-CCA2 scheme
 *
 * \param[out] ct String containing the ciphertext
 * \param[out] ss String containing the shared secret
 * \param[in] pk String containing the public key
 * \return 0 if encapsulation is sucessfull
 */
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk) {

  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  // Computing m
  ffi_vec m;
  ffi_vec_init(&m, PARAM_K);
  ffi_vec_set_random2(m, PARAM_K);

  // Generating G function
  AES_XOF_struct G_seedexpander;
  unsigned char seed_G[VEC_K_BYTES];
  ffi_vec_to_string_compact(seed_G, m, PARAM_K);
  seedexpander_init(&G_seedexpander, seed_G, seed_G + 32, SEEDEXPANDER_MAX_LENGTH);

  // Computing theta
  unsigned char theta[SEEDEXPANDER_SEED_BYTES];
  seedexpander(&G_seedexpander, theta, SEEDEXPANDER_SEED_BYTES);

  // Encrypting m
  ffi_qre u, v;
  ffi_qre_init(&u);
  ffi_qre_init(&v);

  rqc_pke_encrypt(u, v, m, theta, pk);

  // Computing d
  unsigned char d[SHA512_BYTES];
  sha512(d, (unsigned char*) m, FFI_VEC_K_BYTES);

  // Computing ciphertext
  rqc_kem_ciphertext_to_string(ct, u, v, d);

  // Computing shared secret
  unsigned char mc[FFI_VEC_K_BYTES + 2 * FFI_VEC_N_BYTES];
  ffi_vec_to_string(mc, m, PARAM_K);
  ffi_qre_to_string(mc + FFI_VEC_K_BYTES, u);
  ffi_qre_to_string(mc + FFI_VEC_K_BYTES + FFI_VEC_N_BYTES, v);
  sha512(ss, mc, FFI_VEC_K_BYTES + 2 * FFI_VEC_N_BYTES);

  ffi_vec_clear(m, PARAM_K);
  ffi_qre_clear(u);
  ffi_qre_clear(v);
  ffi_qre_clear_modulus();

  return 0;
}



/**
 * \fn int crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk)
 * \brief Decapsulation of the RQC_KEM IND-CCA2 scheme
 *
 * \param[out] ss String containing the shared secret
 * \param[in] ct String containing the ciphertext
 * \param[in] sk String containing the secret key
 * \return 0 if decapsulation is successfull, -1 otherwise
 */
int crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk) {

  ffi_field_init();
  ffi_qre_init_modulus(PARAM_N);

  // Retrieving u, v and d from ciphertext
  ffi_qre u, v;
  ffi_qre_init(&u);
  ffi_qre_init(&v);
  unsigned char d[SHA512_BYTES];
  rqc_kem_ciphertext_from_string(u, v, d, ct);

  // Retrieving pk from sk
  unsigned char pk[PUBLIC_KEY_BYTES];
  memcpy(pk, sk + SEEDEXPANDER_SEED_BYTES, PUBLIC_KEY_BYTES);

  // Decrypting
  ffi_vec m;
  ffi_vec_init(&m, PARAM_K);
  rqc_pke_decrypt(m, u, v, sk);

  // Generating G function
  AES_XOF_struct G_seedexpander;
  unsigned char seed_G[VEC_K_BYTES];
  ffi_vec_to_string_compact(seed_G, m, PARAM_K);
  seedexpander_init(&G_seedexpander, seed_G, seed_G + 32, SEEDEXPANDER_MAX_LENGTH);

  // Computing theta
  unsigned char theta[SEEDEXPANDER_SEED_BYTES];
  seedexpander(&G_seedexpander, theta, SEEDEXPANDER_SEED_BYTES);

  // Encrypting m'
  ffi_qre u2, v2;
  ffi_qre_init(&u2);
  ffi_qre_init(&v2);

  rqc_pke_encrypt(u2, v2, m, theta, pk);

  // Checking that c = c' and abort otherwise
  int abort = 0;
  if(ffi_qre_is_equal_to(u, u2) == 0 || ffi_qre_is_equal_to(v, v2) == 0) {
    abort = 1;
  }

  // Computing d'
  unsigned char d2[SHA512_BYTES];
  sha512(d2, (unsigned char*) m, FFI_VEC_K_BYTES);

  // Checking that d = d' and abort otherwise
  if(memcmp(d, d2, SHA512_BYTES) != 0) {
    abort = 1;
  }

  if(abort == 1) {
    memset(ss, 0, SHARED_SECRET_BYTES);
    ffi_vec_clear(m, PARAM_K);
    ffi_qre_clear(u);
    ffi_qre_clear(v);
    ffi_qre_clear(u2);
    ffi_qre_clear(v2);
    ffi_qre_clear_modulus();
    return -1;
  }

  // Computing shared secret
  unsigned char mc[FFI_VEC_K_BYTES + 2 * FFI_VEC_N_BYTES];
  ffi_vec_to_string(mc, m, PARAM_K);
  ffi_qre_to_string(mc + FFI_VEC_K_BYTES, u);
  ffi_qre_to_string(mc + FFI_VEC_K_BYTES + FFI_VEC_N_BYTES, v);
  sha512(ss, mc, FFI_VEC_K_BYTES + 2 * FFI_VEC_N_BYTES);

  ffi_vec_clear(m, PARAM_K);
  ffi_qre_clear(u);
  ffi_qre_clear(v);
  ffi_qre_clear(u2);
  ffi_qre_clear(v2);
  ffi_qre_clear_modulus();

  return 0;
}

