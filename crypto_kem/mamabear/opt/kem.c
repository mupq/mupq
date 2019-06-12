#include <string.h> /* for memcpy */

#include "api.h"
#include "common.h"
#include "randombytes.h"
#include "threebears.h"
#include "params.h"

int crypto_kem_keypair(
    unsigned char *pk,
    unsigned char *sk
) {
    randombytes(sk,PRIVATE_KEY_BYTES);
    get_pubkey(pk,sk);
    return 0;
}

int crypto_kem_enc(
    unsigned char *ct,
    unsigned char *ss,
    const unsigned char *pk
) {
    unsigned char seed[ENC_SEED_BYTES+IV_BYTES];
    randombytes(seed,sizeof(seed));
    encapsulate(ss,ct,pk,seed);
    secure_bzero(seed,sizeof(seed));
    return 0;
}

int crypto_kem_dec(
    unsigned char *ss,
    const unsigned char *ct,
    const unsigned char *sk
) {
    return decapsulate(ss,ct,sk);
}

