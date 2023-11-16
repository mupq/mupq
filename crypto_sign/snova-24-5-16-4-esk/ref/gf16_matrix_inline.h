#ifndef GF16M_INLINE_H
#define GF16M_INLINE_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "deriv_params.h"
#include "gf16_init.h"
#include "gf16_matrix.h"

// POD -> entry[a][b] * (entry[c][d] * entry[e][f] + entry[g][h] * entry[i][j])
#define POD(entry, a, b, c, d, e, f, g, h, i, j)                          \
    gf16_get_mul(                                                         \
        get_gf16m(entry, a, b),                                           \
        gf16_get_add(                                                     \
            gf16_get_mul(get_gf16m(entry, c, d), get_gf16m(entry, e, f)), \
            gf16_get_mul(get_gf16m(entry, g, h), get_gf16m(entry, i, j))))

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Zeroing the GF16 Matrix a.
 */
static inline void gf16m_set_zero(gf16m_t a) { memset(a, 0, sq_rank); }

/**
 * Adding GF16 Matrices. c = a + b
 */
static inline void gf16m_add(gf16m_t a, gf16m_t b, gf16m_t c) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(c, i, j,
                      gf16_get_add(get_gf16m(a, i, j), get_gf16m(b, i, j)));
        }
    }
}

/**
 * Multiplying GF16 Matrices. c = a * b
 */
static inline void gf16m_mul(gf16m_t a, gf16m_t b, gf16m_t c) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(c, i, j,
                      gf16_get_mul(get_gf16m(a, i, 0), get_gf16m(b, 0, j)));
            for (int k = 1; k < rank; ++k) {
                set_gf16m(c, i, j,
                          gf16_get_add(get_gf16m(c, i, j),
                                       gf16_get_mul(get_gf16m(a, i, k),
                                                    get_gf16m(b, k, j))));
            }
        }
    }
}

/**
 * Scaling the GF16 Matrix. c = Scaling "a" by a factor of "k"
 */
static inline void gf16m_scale(gf16m_t a, gf16_t k, gf16m_t c) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(c, i, j, gf16_get_mul(get_gf16m(a, i, j), k));
        }
    }
}

/**
 * Transposing the GF16 Matrix. ap = aT
 */
static inline void gf16m_transpose(gf16m_t a, gf16m_t ap) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(ap, i, j, get_gf16m(a, j, i));
        }
    }
}

/**
 * Cloning the GF16 Matrix target = source
 */
static inline void gf16m_clone(gf16m_t target, gf16m_t source) {
    memcpy(target, source, sq_rank);
}

/**
 * be_aI
 */
static inline void be_aI(gf16m_t target, gf16_t a) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(target, i, j, (i == j) ? a : 0);
        }
    }
}

/**
 * be_the_S
 */
static inline void be_the_S(gf16m_t target) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(target, i, j, (8 - (i + j)));
        }
    }
}

/**
 * gf16m_det
 */
static inline gf16_t gf16m_det(gf16m_t entry) {
#if rank == 2
    return gf16_get_add(
        gf16_get_mul(get_gf16m(entry, 0, 0), get_gf16m(entry, 1, 1)),
        gf16_get_mul(get_gf16m(entry, 0, 1), get_gf16m(entry, 1, 0)));
    // (entry[0][0] * entry[1][1] + entry[0][1] * entry[1][0]);
#elif rank == 3
    return gf16_get_add(
        gf16_get_add(
            gf16_get_mul(get_gf16m(entry, 0, 0),
                         gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 1),
                                                   get_gf16m(entry, 2, 2)),
                                      gf16_get_mul(get_gf16m(entry, 1, 2),
                                                   get_gf16m(entry, 2, 1)))),
            // AAAAA(entry, 0, 0, 1, 1, 2, 2, 1, 2, 2, 1),
            gf16_get_mul(get_gf16m(entry, 0, 1),
                         gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0),
                                                   get_gf16m(entry, 2, 2)),
                                      gf16_get_mul(get_gf16m(entry, 1, 2),
                                                   get_gf16m(entry, 2, 0))))),
        gf16_get_mul(
            get_gf16m(entry, 0, 2),
            gf16_get_add(
                gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2, 1)),
                gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 0)))));

    /*
    (
            (entry[0][0] * (entry[1][1] * entry[2][2] + entry[1][2] *
    entry[2][1])) + (entry[0][1] * (entry[1][0] * entry[2][2] + entry[1][2] *
    entry[2][0])) + (entry[0][2] * (entry[1][0] * entry[2][1] + entry[1][1] *
    entry[2][0]))
    )

    */

    // gf16_get_mul(gf16_get_mul(get_gf16m(entry, 0, 1),
    // gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2,
    // 2)), gf16_get_mul(get_gf16m(entry, 1, 2), get_gf16m(entry, 2, 0)))))),
    // gf16_get_mul(gf16_get_mul(get_gf16m(entry, 0, 2),
    // gf16_get_add(gf16_get_mul(get_gf16m(entry, 1, 0), get_gf16m(entry, 2,
    // 1)), gf16_get_mul(get_gf16m(entry, 1, 1), get_gf16m(entry, 2, 0))))));

#elif rank == 4

    gf16_t d0 = gf16_get_mul(
        get_gf16m(entry, 0, 0),
        gf16_get_add(gf16_get_add(POD(entry, 1, 1, 2, 2, 3, 3, 2, 3, 3, 2),
                                  POD(entry, 1, 2, 2, 1, 3, 3, 2, 3, 3, 1)),
                     POD(entry, 1, 3, 2, 1, 3, 2, 2, 2, 3, 1)));

    gf16_t d1 = gf16_get_mul(
        get_gf16m(entry, 0, 1),
        gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 2, 3, 3, 2, 3, 3, 2),
                                  POD(entry, 1, 2, 2, 0, 3, 3, 2, 3, 3, 0)),
                     POD(entry, 1, 3, 2, 0, 3, 2, 2, 2, 3, 0)));

    gf16_t d2 = gf16_get_mul(
        get_gf16m(entry, 0, 2),
        gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 3, 2, 3, 3, 1),
                                  POD(entry, 1, 1, 2, 0, 3, 3, 2, 3, 3, 0)),
                     POD(entry, 1, 3, 2, 0, 3, 1, 2, 1, 3, 0)));

    gf16_t d3 = gf16_get_mul(
        get_gf16m(entry, 0, 3),
        gf16_get_add(gf16_get_add(POD(entry, 1, 0, 2, 1, 3, 2, 2, 2, 3, 1),
                                  POD(entry, 1, 1, 2, 0, 3, 2, 2, 2, 3, 0)),
                     POD(entry, 1, 2, 2, 0, 3, 1, 2, 1, 3, 0)));

    return gf16_get_add(gf16_get_add(gf16_get_add(d0, d1), d2), d3);
    /*
    (
            entry[0][0] * (
                    (entry[1][1] * (entry[2][2] * entry[3][3] + entry[2][3] *
    entry[3][2])) + (entry[1][2] * (entry[2][1] * entry[3][3] + entry[2][3] *
    entry[3][1])) + (entry[1][3] * (entry[2][1] * entry[3][2] + entry[2][2] *
    entry[3][1]))

            ) +

            entry[0][1] * (
                    (entry[1][0] * (entry[2][2] * entry[3][3] + entry[2][3] *
    entry[3][2])) + (entry[1][2] * (entry[2][0] * entry[3][3] + entry[2][3] *
    entry[3][0])) + (entry[1][3] * (entry[2][0] * entry[3][2] + entry[2][2] *
    entry[3][0])) ) +

            entry[0][2] * (
                    (entry[1][0] * (entry[2][1] * entry[3][3] + entry[2][3] *
    entry[3][1])) + (entry[1][1] * (entry[2][0] * entry[3][3] + entry[2][3] *
    entry[3][0])) + (entry[1][3] * (entry[2][0] * entry[3][1] + entry[2][1] *
    entry[3][0])) ) +

            entry[0][3] * (
                    (entry[1][0] * (entry[2][1] * entry[3][2] + entry[2][2] *
    entry[3][1])) + (entry[1][1] * (entry[2][0] * entry[3][2] + entry[2][2] *
    entry[3][0])) + (entry[1][2] * (entry[2][0] * entry[3][1] + entry[2][1] *
    entry[3][0]))
    )
    */
#endif
    return 0;
}

/**
 * be_aS
 */
static inline void be_aS(gf16m_t target, gf16_t a) {
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            set_gf16m(target, i, j, gf16_get_mul((8 - (i + j)), a));
        }
    }
}

/**
 * be_invertible_by_add_aS
 */
static inline void be_invertible_by_add_aS(gf16m_t source) {
    gf16m_t temp;
    if (gf16m_det(source) == 0) {
        for (uint8_t a = 1; a < 16; ++a) {
            be_aS(temp, a);
            gf16m_add(temp, source, source);
            if (gf16m_det(source) != 0) {
                return;
            }
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif