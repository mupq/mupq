/*
 * Copyright (c) The mlkem-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLD_INTEGRATION_PQM4_FIPS202_GLUE_H
#define MLD_INTEGRATION_PQM4_FIPS202_GLUE_H

#include "mldsa/src/common.h"

/* Include pqm4's own FIPS202 header */
#include "fips202.h"


#define SHAKE128_RATE 168
#define SHAKE256_RATE 136

#define mld_shake128ctx shake128incctx
#define mld_shake256ctx shake256incctx

static MLD_INLINE void mld_shake128_absorb_once(mld_shake128ctx *state,
                                                const uint8_t *input,
                                                size_t inlen) {
  shake128_inc_init(state);
  shake128_inc_absorb(state, input, inlen);
  shake128_inc_finalize(state);
}

static MLD_INLINE void mld_shake128_squeezeblocks(uint8_t *output,
                                                  size_t nblocks,
                                                  mld_shake128ctx *state) {
  shake128_inc_squeeze(output, nblocks * SHAKE128_RATE, state);
}

static MLD_INLINE void mld_shake128_init(mld_shake128ctx *state) {
  shake128_inc_init(state);
}

static MLD_INLINE void mld_shake128_absorb(mld_shake128ctx *state,
                                           const uint8_t *input, size_t inlen) {
  shake128_inc_absorb(state, input, inlen);
}

static MLD_INLINE void mld_shake128_finalize(mld_shake128ctx *state) {
  shake128_inc_finalize(state);
}

static MLD_INLINE void mld_shake128_squeeze(uint8_t *output, size_t outlen,
                                            mld_shake128ctx *state) {
  shake128_inc_squeeze(output, outlen, state);
}

static MLD_INLINE void mld_shake128_release(mld_shake128ctx *state) {
  shake128_inc_ctx_release(state);
}

static MLD_INLINE void mld_shake256_absorb_once(mld_shake256ctx *state,
                                                const uint8_t *input,
                                                size_t inlen) {
  shake256_inc_init(state);
  shake256_inc_absorb(state, input, inlen);
  shake256_inc_finalize(state);
}

static MLD_INLINE void mld_shake256_squeezeblocks(uint8_t *output,
                                                  size_t nblocks,
                                                  mld_shake256ctx *state) {
  shake256_inc_squeeze(output, nblocks * SHAKE128_RATE, state);
}

static MLD_INLINE void mld_shake256_init(mld_shake256ctx *state) {
  shake256_inc_init(state);
}

static MLD_INLINE void mld_shake256_absorb(mld_shake256ctx *state,
                                           const uint8_t *input, size_t inlen) {
  shake256_inc_absorb(state, input, inlen);
}

static MLD_INLINE void mld_shake256_finalize(mld_shake256ctx *state) {
  shake256_inc_finalize(state);
}

static MLD_INLINE void mld_shake256_squeeze(uint8_t *output, size_t outlen,
                                            mld_shake256ctx *state) {
  shake256_inc_squeeze(output, outlen, state);
}

static MLD_INLINE void mld_shake256_release(mld_shake256ctx *state) {
  shake256_inc_ctx_release(state);
}

static MLD_INLINE void mld_shake256(uint8_t *output, size_t outlen,
                                    const uint8_t *input, size_t inlen) {
  shake256(output, outlen, input, inlen);
}

static MLD_INLINE void mld_sha3_256(uint8_t *output, const uint8_t *input,
                                    size_t inlen) {
  sha3_256(output, input, inlen);
}

static MLD_INLINE void mld_sha3_512(uint8_t *output, const uint8_t *input,
                                    size_t inlen) {
  sha3_512(output, input, inlen);
}

#endif  // MLD_INTEGRATION_PQM4_FIPS202_GLUE_H