#ifndef SPX_WOTS_H
#define SPX_WOTS_H

#include <stdint.h>

#include "params.h"
#include "context.h"
#include "uintx.h"

/**
 * Takes a WOTS signature and an n-byte message, computes a WOTS public key.
 *
 * Writes the computed public key to 'pk'.
 */
#define wots_pk_from_sig SPX_NAMESPACE(wots_pk_from_sig)
void wots_pk_from_sig(unsigned char *pk,
                      const unsigned char *sig, const unsigned char *msg,
                      const spx_ctx *ctx, uint32_t addr[8], uint256_t dp[SPX_WOTS_LEN + 1][(SPX_WOTS_LEN*(SPX_WOTS_W-1)/2) + 1]);

/*
 * Compute the chain lengths needed for a given message hash
 */
#define chain_lengths SPX_NAMESPACE(chain_lengths)
void chain_lengths(unsigned int *lengths, const unsigned char *msg, uint256_t dp[SPX_WOTS_LEN + 1][(SPX_WOTS_LEN*(SPX_WOTS_W-1)/2) + 1]);


#define precompute_dp SPX_NAMESPACE(precompute_dp)
void precompute_dp(uint256_t dp[SPX_WOTS_LEN + 1][(SPX_WOTS_LEN*(SPX_WOTS_W-1)/2) + 1], int l, int w);


#endif
