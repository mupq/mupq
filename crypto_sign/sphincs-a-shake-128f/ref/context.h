#ifndef SPX_CONTEXT_H
#define SPX_CONTEXT_H

#include <stdint.h>

#include "params.h"
#include "sha2.h"

typedef struct {
    uint8_t pub_seed[SPX_N];
    uint8_t sk_seed[SPX_N];

#ifdef SPX_SHA2
    // sha256 state that absorbed pub_seed
    sha256ctx state_seeded;

# if SPX_SHA512
    // sha512 state that absorbed pub_seed
    sha512ctx state_seeded_512;
# endif
#endif

#ifdef SPX_HARAKA
    uint64_t tweaked512_rc64[10][8];
    uint32_t tweaked256_rc32[10][8];
#endif
} spx_ctx;

#endif
