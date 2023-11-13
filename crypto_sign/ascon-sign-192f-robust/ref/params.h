#ifndef ASGN_PARAMS_H
#define ASGN_PARAMS_H

// ASCON-SIGN-192f

/* Hash output length in bytes. */
#define ASGN_N 24
/* Height of the hypertree. */
#define ASGN_FULL_HEIGHT 66
/* Number of subtree layer. */
#define ASGN_D 22
/* FORS tree dimensions. */
#define ASGN_FORS_HEIGHT 8
#define ASGN_FORS_TREES 33
/* Winternitz parameter, */
#define ASGN_WOTS_W 16

/* The hash function is defined by linking a different hash.c file, as opposed
   to setting a #define constant. */

/* For clarity */
#define ASGN_ADDR_BYTES 32

/* WOTS parameters. */
#if ASGN_WOTS_W == 256
    #define ASGN_WOTS_LOGW 8
#elif ASGN_WOTS_W == 16
    #define ASGN_WOTS_LOGW 4
#else
    #error ASGN_WOTS_W assumed 16 or 256
#endif

#define ASGN_WOTS_LEN1 (8 * ASGN_N / ASGN_WOTS_LOGW)

/* ASGN_WOTS_LEN2 is floor(log(len_1 * (w - 1)) / log(w)) + 1; we precompute */
#if ASGN_WOTS_W == 256
    #if ASGN_N <= 1
        #define ASGN_WOTS_LEN2 1
    #elif ASGN_N <= 256
        #define ASGN_WOTS_LEN2 2
    #else
        #error Did not precompute ASGN_WOTS_LEN2 for n outside {2, .., 256}
    #endif
#elif ASGN_WOTS_W == 16
    #if ASGN_N <= 8
        #define ASGN_WOTS_LEN2 2
    #elif ASGN_N <= 136
        #define ASGN_WOTS_LEN2 3
    #elif ASGN_N <= 256
        #define ASGN_WOTS_LEN2 4
    #else
        #error Did not precompute ASGN_WOTS_LEN2 for n outside {2, .., 256}
    #endif
#endif

#define ASGN_WOTS_LEN (ASGN_WOTS_LEN1 + ASGN_WOTS_LEN2)
#define ASGN_WOTS_BYTES (ASGN_WOTS_LEN * ASGN_N)
#define ASGN_WOTS_PK_BYTES ASGN_WOTS_BYTES

/* Subtree size. */
#define ASGN_TREE_HEIGHT (ASGN_FULL_HEIGHT / ASGN_D)

#if ASGN_TREE_HEIGHT * ASGN_D != ASGN_FULL_HEIGHT
    #error ASGN_D should always divide ASGN_FULL_HEIGHT
#endif

/* FORS parameters. */
#define ASGN_FORS_MSG_BYTES ((ASGN_FORS_HEIGHT * ASGN_FORS_TREES + 7) / 8)
#define ASGN_FORS_BYTES ((ASGN_FORS_HEIGHT + 1) * ASGN_FORS_TREES * ASGN_N)
#define ASGN_FORS_PK_BYTES ASGN_N

/* Resulting ASGN sizes. */
#define ASGN_BYTES (ASGN_N + ASGN_FORS_BYTES + ASGN_D * ASGN_WOTS_BYTES +\
                   ASGN_FULL_HEIGHT * ASGN_N)
#define ASGN_PK_BYTES (2 * ASGN_N)
#define ASGN_SK_BYTES (2 * ASGN_N + ASGN_PK_BYTES)

#include "ascon_offsets.h"

#endif
