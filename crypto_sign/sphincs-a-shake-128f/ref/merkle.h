#if !defined( MERKLE_H_ )
#define MERKLE_H_

#include <stdint.h>

/* Generate a Merkle signature (WOTS signature followed by the Merkle */
/* authentication path) */
#define merkle_sign SPX_NAMESPACE(merkle_sign)
void merkle_sign(uint8_t *sig, unsigned char *root,
        const spx_ctx* ctx,
        uint32_t wots_addr[8], uint32_t tree_addr[8],
        uint32_t idx_leaf, uint256_t dp[SPX_WOTS_LEN + 1][(SPX_WOTS_LEN*(SPX_WOTS_W-1)/2) + 1]);

/* Compute the root node of the top-most subtree. */
#define merkle_gen_root SPX_NAMESPACE(merkle_gen_root)
void merkle_gen_root(unsigned char *root, const spx_ctx* ctx, uint256_t dp[SPX_WOTS_LEN + 1][(SPX_WOTS_LEN*(SPX_WOTS_W-1)/2) + 1]);

#endif /* MERKLE_H_ */
