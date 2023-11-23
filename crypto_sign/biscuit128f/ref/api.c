#include <string.h>

#include "biscuit.h"
#include "api.h"
#include <randombytes.h>

#ifndef FIXED_PARAMS
static const params_t params[1] = { MAKE_PARAMS (
    SEC_LEVEL,
    NB_ITERATIONS,
    NB_PARTIES,
    FIELD_SIZE,
    NB_VARIABLES,
    NB_EQUATIONS,
    DEGREE
) };
#endif

static const int siglen = CRYPTO_BYTES;

int
crypto_sign_keypair (unsigned char *pk, unsigned char *sk)
{
  /* use pk buffer to store entropy input for key generation */
  randombytes (sk, SEC_LEVEL >> 2);
  return keygen (sk, pk, sk, params);
}

int
crypto_sign (unsigned char *sm, size_t *smlen,
             const unsigned char *m, size_t mlen,
             const unsigned char *sk)
{
#ifdef DETERMINISTIC_SIGNATURE
  sign (sm, m, mlen, sk, NULL, params);
#else
  /* use sig buffer to store entropy input for signature */
  randombytes (sm, SEC_LEVEL >> 2);
  sign (sm, m, mlen, sk, sm, params);
#endif
  memmove (sm + siglen, m, mlen);
  *smlen = siglen + mlen;
  return 0;
}

int
crypto_sign_open (unsigned char *m, size_t *mlen,
                  const unsigned char *sm, size_t smlen,
                  const unsigned char *pk)
{
  if (verify (sm, sm + siglen, smlen - siglen, pk, params) == 0)
    {
      memmove (m, sm + siglen, smlen - siglen);
      *mlen = smlen - siglen;
      return 0;
    }
  return 1;
}
