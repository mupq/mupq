#include "api.h"

int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk) {

    unsigned long long smlen_ll;
    int rc = _crypto_sign(sm, &smlen_ll, m, mlen, sk);
    *smlen = smlen_ll;
    return rc;
}

int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk) {
    unsigned long long mlen_ll;
    int rc = _crypto_sign_open(m, &mlen_ll, sm, smlen, pk);
    *mlen = mlen_ll;
    return rc;
}