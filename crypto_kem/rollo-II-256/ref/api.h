/**
 * \file api.h
 * \brief NIST KEM API used by ROLLO
 */


#ifndef API_H
#define API_H

#define CRYPTO_ALGNAME "ROLLO-II-256"

#include "rolloII_types.h"

#define CRYPTO_SECRETKEYBYTES 2533
#define CRYPTO_PUBLICKEYBYTES 2493
#define CRYPTO_BYTES 64
#define CRYPTO_CIPHERTEXTBYTES 2621

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, unsigned char* pk);
int crypto_kem_dec(unsigned char* ss, unsigned char* ct, unsigned char* sk);

#endif
