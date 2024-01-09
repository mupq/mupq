#ifndef MQOM_MPC_UTILS_H
#define MQOM_MPC_UTILS_H

#include "mpc-struct.h"
#include <stdlib.h>

// No pointer, to have a way to have concatenate
typedef struct mpc_share_t {
    mpc_wit_t wit;
    mpc_unif_t unif;
    mpc_hint_t hint;
} mpc_share_t;

#define PARAM_WIT_SIZE (sizeof(mpc_wit_t))
#define PARAM_UNIF_SIZE (sizeof(mpc_unif_t))
#define PARAM_HINT_SIZE (sizeof(mpc_hint_t))
#define PARAM_SHARE_SIZE (PARAM_WIT_SIZE+PARAM_UNIF_SIZE+PARAM_HINT_SIZE)
#define PARAM_BR_SIZE (sizeof(mpc_broadcast_t))

static inline mpc_wit_t* get_wit(mpc_share_t* sh) {
    return (mpc_wit_t*) sh;
}

static inline mpc_unif_t* get_unif(mpc_share_t* sh) {
    return (mpc_unif_t*) (((uint8_t*)sh) + PARAM_WIT_SIZE);
}

static inline mpc_hint_t* get_hint(mpc_share_t* sh) {
    return (mpc_hint_t*) (((uint8_t*)sh) + PARAM_WIT_SIZE + PARAM_UNIF_SIZE);
}

static inline void vec_set(void* dst, const void* src, int size) {
    memcpy(dst, src, size);
}

static inline void vec_set_zero(void* x, int size) {
    memset(x, 0, size);
}

#endif /* MQOM_MPC_UTILS_H */
