/*
 * Copyright (c) The mlkem-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLK_INTEGRATION_PQM4_FIPS202X4_GLUE_H
#define MLK_INTEGRATION_PQM4_FIPS202X4_GLUE_H

/* Include pqm4's own FIPS202 header */
#include "fips202.h"

typedef struct {
    shake128ctx ctx[4];
} mlk_shake128x4ctx;

#define mlk_shake128x4_absorb_once(S, IN0, IN1, IN2, IN3, INLEN) \
    do {                                                         \
        shake128_absorb(&S->ctx[0], IN0, INLEN);                 \
        shake128_absorb(&S->ctx[1], IN1, INLEN);                 \
        shake128_absorb(&S->ctx[2], IN2, INLEN);                 \
        shake128_absorb(&S->ctx[3], IN3, INLEN);                 \
    } while(0)

#define mlk_shake128x4_squeezeblocks(OUT0, OUT1, OUT2, OUT3, NBLOCKS, S) \
    do {                                                                 \
        shake128_squeezeblocks(OUT0, NBLOCKS, &S->ctx[0]);               \
        shake128_squeezeblocks(OUT1, NBLOCKS, &S->ctx[1]);               \
        shake128_squeezeblocks(OUT2, NBLOCKS, &S->ctx[2]);               \
        shake128_squeezeblocks(OUT3, NBLOCKS, &S->ctx[3]);               \
    } while(0)

#define mlk_shake128x4_init(S) do { } while(0) /* no-op */
#define mlk_shake128x4_release(S) do { } while(0) /* no-op */

#define mlk_shake256x4(OUT0, OUT1, OUT2, OUT3, OUTLEN, IN0, IN1, IN2, IN3, INLEN) \
    do {                                                                          \
        shake256(OUT0, OUTLEN, IN0, INLEN);                                       \
        shake256(OUT1, OUTLEN, IN1, INLEN);                                       \
        shake256(OUT2, OUTLEN, IN2, INLEN);                                       \
        shake256(OUT3, OUTLEN, IN3, INLEN);                                       \
    } while (0)

#endif
