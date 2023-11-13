#ifndef ASGN_HASH_H
#define ASGN_HASH_H

#include <stdint.h>
#include "context.h"
#include "params.h"

void initialize_hash_function(ascon_sign_ctx *ctx);

void prf_addr(unsigned char *out, const ascon_sign_ctx *ctx,
              const uint32_t addr[8]);

void gen_message_random(unsigned char *R, const unsigned char *sk_prf,
                        const unsigned char *optrand,
                        const unsigned char *m, unsigned long long mlen,
                        const ascon_sign_ctx *ctx);

void hash_message(unsigned char *digest, uint64_t *tree, uint32_t *leaf_idx,
                  const unsigned char *R, const unsigned char *pk,
                  const unsigned char *m, unsigned long long mlen,
                  const ascon_sign_ctx *ctx);

#endif
