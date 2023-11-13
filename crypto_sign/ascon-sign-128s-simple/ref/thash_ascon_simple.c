#include <stdint.h>
#include <string.h>

#include "thash.h"
#include "address.h"
#include "params.h"
#include "utils.h"

#include "ascon_api.h"

/**
 * Takes an array of inblocks concatenated arrays of ASGN_N bytes.
 */
void thash(unsigned char *out, const unsigned char *in, unsigned int inblocks,
           const ascon_sign_ctx *ctx, uint32_t addr[8])
{
    ASGN_VLA(uint8_t, buf, ASGN_N + ASGN_ADDR_BYTES + inblocks*ASGN_N);

    memcpy(buf, ctx->pub_seed, ASGN_N);
    memcpy(buf + ASGN_N, addr, ASGN_ADDR_BYTES);
    memcpy(buf + ASGN_N + ASGN_ADDR_BYTES, in, inblocks * ASGN_N);

    ascon_hash(out, ASGN_N, buf, ASGN_N + ASGN_ADDR_BYTES + inblocks*ASGN_N); 
}
