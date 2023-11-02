#ifndef FIPS202_WRAPPER_H
#define FIPS202_WRAPPER_H

#include <stddef.h>
#include <stdint.h>

#include "fips202.h"

#define keccak_state shake256incctx
#define shake256_init shake256_inc_init
#define shake256_absorb shake256_inc_absorb
#define shake256_finalize shake256_inc_finalize
#define shake256_squeeze shake256_inc_squeeze

static inline void shake256_absorb_once(keccak_state *state, const uint8_t *in, size_t inlen)
{
  shake256_init(state);
  shake256_absorb(state, in, inlen);
  shake256_finalize(state);
}

#endif
