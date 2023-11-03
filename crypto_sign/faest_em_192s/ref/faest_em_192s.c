/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest_em_192s.h"
#include "compat.h"
#include "randomness.h"
#include "owf.h"
#include "instances.h"
#include "faest.h"

#include <stdlib.h>
#include <string.h>

int FAEST_CALLING_CONVENTION faest_em_192s_keygen(uint8_t* pk, uint8_t* sk) {
  if (!pk || !sk) {
    return -1;
  }

  uint8_t randomness[48 / 2 + FAEST_EM_192S_PRIVATE_KEY_SIZE];

  uint8_t* input = pk;
  uint8_t* output = pk + 48 / 2;

  bool done = false;
  while (!done) {
    rand_bytes(randomness, 48 / 2 + FAEST_EM_192S_PRIVATE_KEY_SIZE);
    memcpy(input, randomness, 48 / 2);
    memcpy(sk, randomness + 48 / 2, FAEST_EM_192S_PRIVATE_KEY_SIZE);
    done = faest_em_192s_owf(sk, input, output);
  }

  faest_declassify(pk, 48);
  return 0;
}

int FAEST_CALLING_CONVENTION faest_em_192s_validate_keypair(const uint8_t* pk, const uint8_t* sk) {
  if (!sk || !pk) {
    return -1;
  }

  const uint8_t* input = pk;
  const uint8_t* expected_output = pk + 48 / 2;

  uint8_t output[48 / 2] = { 0 };
  if (!faest_em_192s_owf(sk, input, output)) {
    // zero bytes in SubBytes input
    return 1;
  }

  return faest_timingsafe_bcmp(output, expected_output, sizeof(output)) == 0 ? 0 : 2;
}

int FAEST_CALLING_CONVENTION faest_em_192s_sign(const uint8_t* sk, const uint8_t* pk, const uint8_t* message, size_t message_len, const uint8_t* rho, size_t rho_len, uint8_t* signature, size_t* signature_len) {
  if (!sk || !signature || !signature_len || *signature_len < FAEST_EM_192S_SIGNATURE_SIZE || (!rho && rho_len)) {
    return -1;
  }

  const faest_paramset_t params = faest_get_paramset(FAEST_EM_192S);
  signature_t sig = init_signature(&params);

  sign(message, message_len, sk, pk, rho, rho_len, &params, &sig);
  serialize_signature(signature, signature_len, &sig, &params);
  free_signature(sig, &params);
  if (*signature_len != FAEST_EM_192S_SIGNATURE_SIZE) {
    // something went wrong
    memset(signature, 0, FAEST_EM_192S_SIGNATURE_SIZE);
    return -2;
  }

  return 0;
}

int FAEST_CALLING_CONVENTION faest_em_192s_verify(const uint8_t* pk, const uint8_t* message, size_t message_len, const uint8_t* signature, size_t signature_len) {
  if (!pk || !signature || signature_len != FAEST_EM_192S_SIGNATURE_SIZE) {
    return -1;
  }

  const faest_paramset_t params = faest_get_paramset(FAEST_EM_192S);
  signature_t sig = deserialize_signature(signature, &params);

  int ret = verify(message, message_len, pk, &params, &sig);
  free_signature(sig, &params);

  return ret;
}

void FAEST_CALLING_CONVENTION faest_em_192s_clear_private_key(uint8_t* key) {
  faest_explicit_bzero(key, FAEST_EM_192S_PRIVATE_KEY_SIZE);
}

// vim: ft=c
