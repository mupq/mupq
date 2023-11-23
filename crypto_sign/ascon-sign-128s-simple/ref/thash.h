#ifndef ASGN_THASH_H
#define ASGN_THASH_H

#include "context.h"
#include "params.h"

#include <stdint.h>

void thash(unsigned char *out, const unsigned char *in, unsigned int inblocks,
           const ascon_sign_ctx *ctx, uint32_t addr[8]);

#endif
