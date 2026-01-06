/*
 * Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLD_INTEGRATION_PQM4_API_H
#define MLD_INTEGRATION_PQM4_API_H

#define MLD_CONFIG_NO_SUPERCOP
#include "mldsa/mldsa_native.h"

#define CRYPTO_SECRETKEYBYTES MLDSA44_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES MLDSA44_PUBLICKEYBYTES
#define CRYPTO_BYTES MLDSA44_BYTES

#define crypto_sign_keypair mldsa_keypair

static inline int crypto_sign(uint8_t *sm, size_t *smlen,
                              const uint8_t *m, size_t mlen,
                              const uint8_t *sk)
{
    return mldsa_sign(sm, smlen, m, mlen, NULL, 0, sk);
}

static inline int crypto_sign_open(uint8_t *m, size_t *mlen,
                                   const uint8_t *sm, size_t smlen,
                                   const uint8_t *pk)
{
    return mldsa_open(m, mlen, sm, smlen, NULL, 0, pk);
}

#endif