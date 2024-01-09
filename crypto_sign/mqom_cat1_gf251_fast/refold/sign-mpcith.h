#ifndef MQOM_SIGN_HYPERCUBE_7R_H
#define MQOM_SIGN_HYPERCUBE_7R_H

#include <stdint.h>
#include <stddef.h>

int mpcith_hypercube_7r_sign(uint8_t* sig, size_t* siglen,
                const uint8_t* m, size_t mlen,
                const uint8_t* sk,
                const uint8_t* salt, const uint8_t* seed
                );
int mpcith_hypercube_7r_sign_verify(const uint8_t* sig, size_t siglen,
                const uint8_t* m, size_t mlen,
                const uint8_t* pk
                );
#define mqom_sign(sig, siglen, m, mlen, sk, salt, seed) mpcith_hypercube_7r_sign(sig, siglen, m, mlen, sk, salt, seed)
#define mqom_sign_verify(sig, siglen, m, mlen, pk) mpcith_hypercube_7r_sign_verify(sig, siglen, m, mlen, pk)

#endif /* MQOM_SIGN_HYPERCUBE_7R_H */
