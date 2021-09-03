#ifndef PICNIC_L1_FS_API_H
#define PICNIC_L1_FS_API_H

#define CRYPTO_SECRETKEYBYTES (1 + 2 * 16 + 16)
#define CRYPTO_PUBLICKEYBYTES (1 + 2 * 16)
#define CRYPTO_BYTES (4 + 34032)
#define CRYPTO_ALGNAME "picnicl1fs"
#define CRYPTO_VERSION "3.0"
#define CRYPTO_DETERMINISTIC 1

#include <stddef.h>

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk);
int crypto_sign(unsigned char* sm, size_t* smlen, const unsigned char* m,
                size_t mlen, const unsigned char* sk);
int crypto_sign_open(unsigned char* m, size_t* mlen, const unsigned char* sm,
                     size_t smlen, const unsigned char* pk);

#endif
