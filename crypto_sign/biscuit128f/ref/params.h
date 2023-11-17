#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "params_mpcith.h"
#include "params_posso.h"
#include "params_biscuit.h"
#include "params_instance.h"

#if !defined(SEC_LEVEL) ||                      \
  !defined(NB_ITERATIONS) ||                    \
  !defined(NB_PARTIES) ||                       \
  !defined(FIELD_SIZE) ||                       \
  !defined(NB_VARIABLES) ||                     \
  !defined(NB_EQUATIONS) ||                     \
  !defined(DEGREE)
#error Missing parameters
#endif

#define PRIVKEY_BYTES                                                   \
  GET_PRIVKEY_BYTES ((SEC_LEVEL),                                       \
                     (NB_ITERATIONS), (NB_PARTIES),                     \
                     (FIELD_SIZE), (NB_VARIABLES), (NB_EQUATIONS), (DEGREE))
#define PUBKEY_BYTES                                                    \
  GET_PUBKEY_BYTES ((SEC_LEVEL),                                        \
                    (NB_ITERATIONS), (NB_PARTIES),                      \
                    (FIELD_SIZE), (NB_VARIABLES), (NB_EQUATIONS), (DEGREE))
#define SIGNATURE_BYTES                                                 \
  GET_SIGNATURE_BYTES ((SEC_LEVEL),                                     \
                       (NB_ITERATIONS), (NB_PARTIES),                   \
                       (FIELD_SIZE), (NB_VARIABLES), (NB_EQUATIONS), (DEGREE))


#endif
