#ifndef SIMPLE_EXAMPLE_RNG_H
#define SIMPLE_EXAMPLE_RNG_H

#include "param.h"
#include "types.h"
#include <stdint.h>

#define EXPORT

#define HASH_COM 0
#define HASH_H1 1
#define HASH_H2 2
#define HASH_TREE 3

/** @brief opaque structure that represents a HASH context */
typedef struct HASH_struct HASH_CTX;
HASH_CTX *sdith_hash_create_hash_ctx(uint8_t prefix);
void sdith_hash_free_hash_ctx(HASH_CTX *ctx);
void sdith_hash_digest_update(HASH_CTX *ctx, void const *in, int inBytes);
void sdith_hash_final(HASH_CTX *ctx, void *dest);
void sdith_hash(uint8_t prefix, void *dest, void const *data, int dataBytes);

typedef struct HASH4_struct HASH4_CTX;
HASH4_CTX *sdith_hash_create_hash4_ctx(uint8_t prefix);
void sdith_hash_free_hash4_ctx(HASH4_CTX *ctx);
void sdith_hash4_digest_update(HASH4_CTX *ctx, void **in, int inBytes);
void sdith_hash4_final(HASH4_CTX *ctx, void **dest);
void sdith_hash4(uint8_t prefix, void **dest, void **data, int dataBytes);

/** @brief opaque structure that represents a XOF context */
typedef struct XOF_struct XOF_CTX; // opaque structure
/** @brief creates a xof context for randomness expansion */
EXPORT XOF_CTX *sdith_rng_create_xof_ctx(void *in, int inBytes);
/** @brief deletes a rng context instantiated with sdith_rng_create_xof_ctx */
EXPORT void sdith_rng_free_xof_ctx(XOF_CTX *ctx);
/** @brief produces the next random bytes out of this context */
EXPORT void sdith_xof_next_bytes(XOF_CTX *ctx, void *out, int outLen);
/** @brief produces the next random bytes out of this context, but each byte is
 * sampled within [0, 251) using rejection sampling. */
EXPORT void sdith_xof_next_bytes_mod251(XOF_CTX *ctx, void *out, int outLen);

/** @brief opaque structure that represents a XOF context */
typedef struct XOF4_struct XOF4_CTX; // opaque structure
/** @brief creates a xof context for randomness expansion */
EXPORT XOF4_CTX *sdith_rng_create_xof4_ctx(void **in, int inBytes);
/** @brief deletes a rng context instantiated with sdith_rng_create_xof_ctx */
EXPORT void sdith_rng_free_xof4_ctx(XOF4_CTX *ctx);
/** @brief produces the next random bytes out of this context */
EXPORT void sdith_xof4_next_bytes(XOF4_CTX *ctx, void **out, int outLen);
/** @brief produces the next random bytes out of this context, but each byte is
 * sampled within [0, 251) using rejection sampling. */
EXPORT void sdith_xof4_next_bytes_mod251(XOF4_CTX *ctx, void **out, int outLen);

/** @brief opaque structure that represents a TREE_PRG context */
typedef struct TREE_PRG_struct TREE_PRG_CTX;
/** @brief takes as input n/2 seeds and produces n seeds:
 * out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_prg_seed_expand(TREE_PRG_CTX *key, void *out,
                                       const void *in,
                                       const uint16_t first_tweak,
                                       const uint16_t iteration, uint64_t n);
/** @brief takes a tree_prg context */
EXPORT TREE_PRG_CTX *sdith_create_tree_prg_ctx(void const *const root_salt);
/** @brief free a tree_prg context */
EXPORT void sdith_free_tree_prg_ctx(TREE_PRG_CTX *key);
/** @brief expand the leaf into seed and rho */
EXPORT void sdith_tree_prg_leaf_expand(TREE_PRG_CTX *key, seed_t in_seed,
                                       seed_t out_seed,
                                       uint8_t out_rho[PARAM_rho_size]);
/** @brief expand the leaf into seed and rho 4x */
EXPORT void sdith_tree_prg_leaf_expand4(TREE_PRG_CTX *key, seed_t *in_seed,
                                        seed_t *out_seed, uint8_t **out_rho);

#endif // SIMPLE_EXAMPLE_RNG_H
