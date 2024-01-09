#ifndef MPC_STRUCT_H
#define MPC_STRUCT_H

#include "parameters.h"
#include "witness.h"

#define PARAM_PRECOMPUTE_LIN_A

typedef solution_t mpc_wit_t;

typedef struct mpc_unif_t {
    uint8_t a[PARAM_n2][PARAM_eta];
} mpc_unif_t;

typedef struct mpc_challenge_1_t {
    uint8_t gamma[PARAM_m][PARAM_eta];
    // Precomputed
    #ifdef PARAM_PRECOMPUTE_LIN_A
    uint8_t linA[PARAM_MATRIX_BYTESIZE][PARAM_eta];
    #endif
} mpc_challenge_1_t;

typedef struct mpc_hint_t {
    uint8_t q_poly[2*PARAM_n1-1][PARAM_eta];
} mpc_hint_t;

typedef struct mpc_challenge_2_t {
    uint8_t r[PARAM_eta];
    // Precomputed
    uint8_t r_powers[2*PARAM_n1][PARAM_eta];
} mpc_challenge_2_t;

typedef struct mpc_broadcast_t {
    uint8_t alpha[PARAM_n2][PARAM_eta];
    uint8_t v[PARAM_eta];
} mpc_broadcast_t;

#include "field.h"
static inline void vec_add(void* x, const void* y, uint32_t size) {
    add_tab_points(x, y, size);
}
static inline void vec_sub(void* x, const void* y, uint32_t size) {
    sub_tab_points(x, y, size);
}
static inline void vec_neg(void* x, uint32_t size) {
    neg_tab_points(x, size);
}
static inline void vec_rnd(void* x, uint32_t size, samplable_t* entropy) {
    random_points(x, size, entropy);
}
static inline void vec_rnd_x4(void* const* x, uint32_t size, samplable_x4_t* entropy) {
    random_points_x4(x, size, entropy);
}
static inline void vec_normalize(void* x, uint32_t size) {
    (void) x;
    (void) size;
}

#include "mpc-compression.h"

#endif /* MPC_STRUCT_H */
