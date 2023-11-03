/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef CRYPTO_SIGN_192F_H
#define CRYPTO_SIGN_192F_H

#define CRYPTO_SECRETKEYBYTES (24 + 64 / 2)
#define CRYPTO_PUBLICKEYBYTES 64
#define CRYPTO_BYTES 16792
#define CRYPTO_ALGNAME "faest_192f"

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk);

#include <stddef.h>

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk);

int _crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m,
                unsigned long long mlen, const unsigned char* sk);
int _crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm,
                     unsigned long long smlen, const unsigned char* pk);



#endif

// vim: ft=c
