#ifndef ASGN_CONTEXT_H
#define ASGN_CONTEXT_H

#include <stdint.h>

#include "params.h"

typedef struct {
    uint8_t pub_seed[ASGN_N];
    uint8_t sk_seed[ASGN_N];
} ascon_sign_ctx;

#endif
