#ifndef api_h
#define api_h

#define CRYPTO_SECRETKEYBYTES 40
#define CRYPTO_PUBLICKEYBYTES 1194
#define CRYPTO_BYTES 32
#define CRYPTO_CIPHERTEXTBYTES 1307
#define CRYPTO_ALGNAME "MamaBear"

int crypto_kem_keypair(
    unsigned char *pk,
    unsigned char *sk
);

int crypto_kem_enc(
    unsigned char *ct,
    unsigned char *ss,
    const unsigned char *pk
);

int crypto_kem_dec(
    unsigned char *ss,
    const unsigned char *ct,
    const unsigned char *sk
);

#endif /* api_h */

