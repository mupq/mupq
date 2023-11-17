#ifndef _BISCUIT_H_
#define _BISCUIT_H_

#include <stdint.h>

#include "params_mpcith.h"
#include "params_posso.h"
#include "params_biscuit.h"

#if defined(SEC_LEVEL) &&                      \
  defined(NB_ITERATIONS) &&                    \
  defined(NB_PARTIES) &&                       \
  defined(FIELD_SIZE) &&                       \
  defined(NB_VARIABLES) &&                     \
  defined(NB_EQUATIONS) &&                     \
  defined(DEGREE)
#define FIXED_PARAMS
#endif

#ifndef FIXED_PARAMS
typedef struct
{
#ifndef SEC_LEVEL
  int lambda;
#endif
#ifndef NB_ITERATIONS
  int tau;
#endif
#ifndef NB_PARTIES
  int N;
#endif
#ifndef FIELD_SIZE
  int q;
#endif
#ifndef NB_VARIABLES
  int n;
#endif
#ifndef NB_EQUATIONS
  int m;
#endif
#ifndef DEGREE
  int d;
#endif
} params_t;
#else
typedef int params_t;
#endif

#ifdef SEC_LEVEL
#define GET_SEC_LEVEL(lambda)
#else
#define GET_SEC_LEVEL(lambda) lambda,
#endif

#ifdef NB_ITERATIONS
#define GET_NB_ITERATIONS(tau)
#else
#define GET_NB_ITERATIONS(tau) tau,
#endif

#ifdef NB_PARTIES
#define GET_NB_PARTIES(N)
#else
#define GET_NB_PARTIES(N) N,
#endif

#ifdef FIELD_SIZE
#define GET_FIELD_SIZE(q)
#else
#define GET_FIELD_SIZE(q) q,
#endif

#ifdef NB_VARIABLES
#define GET_NB_VARIABLES(n)
#else
#define GET_NB_VARIABLES(n) n,
#endif

#ifdef NB_EQUATIONS
#define GET_NB_EQUATIONS(m)
#else
#define GET_NB_EQUATIONS(m) m,
#endif

#ifdef DEGREE
#define GET_DEGREE(d)
#else
#define GET_DEGREE(d) d,
#endif

#ifndef FIXED_PARAMS

#define MAKE_PARAMS(lambda, tau, N, q, n, m, d) {\
  GET_SEC_LEVEL (lambda) \
  GET_NB_ITERATIONS (tau) \
  GET_NB_PARTIES (N) \
  GET_FIELD_SIZE (q) \
  GET_NB_VARIABLES (n) \
  GET_NB_EQUATIONS (m) \
  GET_DEGREE (d) \
}

int
keygen (uint8_t *sk, uint8_t *pk, const uint8_t *entropy,
        const params_t *params);

int
sign (uint8_t *sig, const uint8_t *msg, uint64_t msglen, const uint8_t *sk,
      const uint8_t *entropy, const params_t *params);

int
verify (const uint8_t *sig, const uint8_t *msg, uint64_t msglen,
        const uint8_t *pk, const params_t *params);

#else

#define MAKE_PARAMS(lambda, tau, N, q, n, m, d) 0

int
_keygen (uint8_t *sk, uint8_t *pk, const uint8_t *entropy);

int
_sign (uint8_t *sig, const uint8_t *msg, uint64_t msglen, const uint8_t *sk,
       const uint8_t *entropy);

int
_verify (const uint8_t *sig, const uint8_t *msg, uint64_t msglen,
         const uint8_t *pk);

#define keygen(sk, pk, entropy, params) _keygen (sk, pk, entropy)
#define sign(sig, msg, msglen, sk, entropy, params) \
  _sign (sig, msg, msglen, sk, entropy)
#define verify(sig, msg, msglen, pk, params) _verify(sig, msg, msglen, pk)

#endif

#endif
