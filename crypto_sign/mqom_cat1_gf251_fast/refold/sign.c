#ifdef SUPERCOP
#include "crypto_sign.h"
#else
#include "api.h"
#endif

#include <string.h>
#include "parameters.h"
#include "keygen.h"
#include "sign-mpcith.h"
#include "rnd.h"
#include "randombytes.h"

static void serialize_uint32(uint8_t buf[4], uint32_t x) {
    buf[0] = x & 0xFF;
    buf[1] = (x >> 8) & 0xFF;
    buf[2] = (x >> 16) & 0xFF;
    buf[3] = x >> 24;
}
static uint32_t deserialize_uint32(const uint8_t buf[4]) {
    return ((uint32_t) buf[0]) | ((uint32_t) buf[1] << 8) | ((uint32_t) buf[2] << 16) | ((uint32_t) buf[3] << 24);
}

int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
    // Sample salt
    uint8_t salt[PARAM_SALT_SIZE];
    sample_salt(salt);

    // Sample seed
    uint8_t seed[PARAM_SEED_SIZE];
    sample_seed(seed);

    size_t siglen2 = 0;
    // Compute the signature
    int ret = mqom_sign(sig, &siglen2, m, mlen, sk, salt, seed);
    *siglen = (size_t) siglen2;
    return ret;
}

int crypto_sign(
        unsigned char *sm, size_t *smlen,
        const unsigned char *m, size_t mlen,
        const unsigned char *sk) {

    uint32_t signature_len = 0;
    memmove(sm + sizeof(signature_len), m, mlen);
    int ret = crypto_sign_signature(sm + sizeof(signature_len) + mlen, smlen, sm + sizeof(signature_len), mlen, sk);
    signature_len = (uint32_t) *smlen;
    *smlen += mlen + sizeof(uint32_t);
    serialize_uint32(sm, signature_len);
    return ret;
}

int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
    return mqom_sign_verify(sig, siglen, m, mlen, pk);
}

int crypto_sign_open(
        unsigned char *m, size_t *mlen,
        const unsigned char *sm, size_t smlen,
        const unsigned char *pk) {

    uint32_t signature_len = 0;
    // The signature is too short to hold the signature length.
    if (smlen < sizeof(signature_len)) {
        return -1;
    }
    
    signature_len = deserialize_uint32(sm);
    // The signature is too short to hold the signature
    if (signature_len + sizeof(uint32_t) > smlen) {
        return -1;
    }

    const size_t message_len = smlen - signature_len - sizeof(uint32_t);
    const uint8_t* message   = sm + sizeof(uint32_t);
    const uint8_t* sig       = sm + sizeof(uint32_t) + message_len;

    int ret = crypto_sign_verify(sig, signature_len, message, message_len, pk);
    if (ret) {
        return ret;
    }

    memmove(m, message, message_len);
    *mlen = message_len;
    return 0;
}

