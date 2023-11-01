#ifndef API_H
#define API_H

#define CRYPTO_SECRETKEYBYTES 5872
#define CRYPTO_PUBLICKEYBYTES 55604
#define CRYPTO_BYTES 54736

#define CRYPTO_ALGNAME "MEDS55604"

#ifndef PQM4
#define PQM4
#endif


int crypto_sign_keypair(
    unsigned char *pk,
    unsigned char *sk
  );

#ifndef PQM4

int crypto_sign(
    unsigned char *sm, unsigned long long *smlen,
    const unsigned char *m, unsigned long long mlen,
    const unsigned char *sk
  );

int crypto_sign_open(
    unsigned char *m, unsigned long long *mlen,
    const unsigned char *sm, unsigned long long smlen,
    const unsigned char *pk
  );
#else
#include <stddef.h>

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk);

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk);
#endif

#endif

