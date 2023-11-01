#include "api.h"
#include "meds.h"

#ifndef PQM4
int
crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk) {
    return _crypto_sign(sm, smlen, m, mlen, sk);
}

int
crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk) {
    return _crypto_sign_open(m, mlen, sm, smlen, pk);
}


#else
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
#endif