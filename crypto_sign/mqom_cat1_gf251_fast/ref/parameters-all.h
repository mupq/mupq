#ifndef MQOM_ALL_PARAMETERS_H
#define MQOM_ALL_PARAMETERS_H

#include "parameters.h"
#include "mpc-all.h"

// Name of the variant
#define PARAM_VARIANT "hypercube-7r"

/****************************************************
 *    Check that the required parameters are set    *
 ****************************************************/

// Require PARAM_NB_EXECUTIONS
#ifndef PARAM_NB_EXECUTIONS
#error Missing parameter: PARAM_NB_EXECUTIONS
#endif

// Require PARAM_NB_PARTIES
#ifndef PARAM_NB_PARTIES
#error Missing parameter: PARAM_NB_PARTIES
#endif

// Require PARAM_LOG_NB_PARTIES
#ifndef PARAM_LOG_NB_PARTIES
#error Missing parameter: PARAM_LOG_NB_PARTIES
#endif

// Require PARAM_HYPERCUBE_DIMENSION or PARAM_LOG_NB_PARTIES
#if !defined(PARAM_HYPERCUBE_DIMENSION) && !defined(PARAM_LOG_NB_PARTIES)
#error Missing parameter: PARAM_HYPERCUBE_DIMENSION or PARAM_LOG_NB_PARTIES
#endif
#ifndef PARAM_HYPERCUBE_DIMENSION
#define PARAM_HYPERCUBE_DIMENSION (PARAM_LOG_NB_PARTIES)
#endif
#ifndef PARAM_LOG_NB_PARTIES
#define PARAM_LOG_NB_PARTIES (PARAM_HYPERCUBE_DIMENSION)
#endif

/****************************************************
 *     Definition of the maximal signature size     *
 ****************************************************/

#define PARAM_SIGNATURE_SIZEBYTES (                   \
    4 + PARAM_DIGEST_SIZE + PARAM_SALT_SIZE +         \
    + PARAM_NB_EXECUTIONS*(                             \
        PARAM_WIT_SHORT_SIZE + PARAM_HINT_SHORT_SIZE    \
            + PARAM_COMPRESSED_BR_SIZE                  \
            + PARAM_HYPERCUBE_DIMENSION*PARAM_SEED_SIZE \
            + 2*PARAM_DIGEST_SIZE                       \
    )                                                   \
)


#include "serialization.h"

#endif /* MQOM_ALL_PARAMETERS_H */
