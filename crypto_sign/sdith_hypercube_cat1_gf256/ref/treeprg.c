#include <stdlib.h>
#include <string.h>

#include "assertions.h"
#include "param.h"
#include "rng.h"
#include "types.h"



typedef struct seed_expand_vec {
  salt_t salt;
  uint16_t iteration;
  uint16_t node_idx;
  seed_t parent_seed;
} seed_expand_vec_t;

/** @brief takes as input n/2 seeds and produces n seeds:
 * out[i]=F_{first_tweak+i}(in[i/2]) */
EXPORT void sdith_tree_prg_seed_expand(TREE_PRG_CTX *context, void *out,
                                       const void *in,
                                       const uint16_t first_tweak,
                                       const uint16_t iteration, uint64_t n) {


  ASSERT_DRAMATICALLY(n % 2 == 0, "n must be even for this function");
  const seed_t *in_seeds = (seed_t *)in;
  seed_t(*out_seeds)[2] = (seed_t(*)[2])out;

  // the size of the current layer is also its starting index.
  uint32_t tweak = first_tweak;

  if (n < 8) {
    seed_expand_vec_t vec;
    memcpy(vec.salt, context->salt, PARAM_salt_size);
    vec.iteration = iteration;

    for (uint64_t i = 0; i < n / 2; ++i) {
      // out[i, i+1] = sha3(salt||iteration||node_idx||in[i])
      vec.node_idx = tweak;
      memcpy(vec.parent_seed, in_seeds + i, PARAM_seed_size);
      sdith_hash(HASH_TREE, out_seeds[i], &vec, sizeof(vec));
      ++tweak;
    }
    return;
  }
  
  seed_expand_vec_t vec4[4];
  void *src4[4] = {&vec4[0], &vec4[1], &vec4[2], &vec4[3]};
  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(vec4[i].salt, context->salt, PARAM_salt_size);
    vec4[i].iteration = iteration;
  }
  for (uint64_t i = 0; i < n / 2; i += 4) {
    vec4[0].node_idx = tweak;
    memcpy(vec4[0].parent_seed, in_seeds + i, PARAM_seed_size);
    vec4[1].node_idx = tweak + 1;
    memcpy(vec4[1].parent_seed, in_seeds + i + 1, PARAM_seed_size);
    vec4[2].node_idx = tweak + 2;
    memcpy(vec4[2].parent_seed, in_seeds + i + 2, PARAM_seed_size);
    vec4[3].node_idx = tweak + 3;
    memcpy(vec4[3].parent_seed, in_seeds + i + 3, PARAM_seed_size);

    void *dst4[4] = {out_seeds[i], out_seeds[i + 1], out_seeds[i + 2],
                     out_seeds[i + 3]};
    sdith_hash4(HASH_TREE, dst4, src4, sizeof(seed_expand_vec_t));
  }
  
}

/** @brief takes a tree_prg context */
EXPORT void sdith_create_tree_prg_ctx(TREE_PRG_CTX* ctx, void const *const root_salt) {
  memcpy(ctx->salt, root_salt, PARAM_salt_size);
}

/** @brief free a tree_prg context */
EXPORT void sdith_free_tree_prg_ctx(TREE_PRG_CTX *key) {
   (void) key;
}

/** @brief expand the leaf into seed and rho */
EXPORT void sdith_tree_prg_leaf_expand(TREE_PRG_CTX *context, seed_t in_seed,
                                       seed_t out_seed,
                                       uint8_t out_rho[PARAM_rho_size]) {
  uint8_t in[PARAM_seed_size + PARAM_salt_size];
  memcpy(in, in_seed, PARAM_seed_size);
  memcpy(&in[PARAM_seed_size], context->salt, PARAM_salt_size);

  uint8_t out[PARAM_seed_size + PARAM_rho_size];
  sdith_hash(HASH_TREE, out, in, sizeof(in));
  memcpy(out_seed, out, PARAM_seed_size);
  memcpy(out_rho, &out[PARAM_seed_size], PARAM_rho_size);
}

/** @brief expand the leaf into seed and rho x4 */
EXPORT void sdith_tree_prg_leaf_expand4(TREE_PRG_CTX *context, seed_t *in_seed,
                                        seed_t *out_seed, uint8_t **out_rho) {
  uint8_t in4[4][PARAM_seed_size + PARAM_salt_size];
  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(in4[i], in_seed[i], PARAM_seed_size);
    memcpy(&in4[i][PARAM_seed_size], context->salt, PARAM_salt_size);
  }
  void *inVec[4] = {&in4[0], &in4[1], &in4[2], &in4[3]};
  uint8_t out4[4][PARAM_seed_size + PARAM_rho_size];
  void *outVec[4] = {&out4[0], &out4[1], &out4[2], &out4[3]};

  sdith_hash4(HASH_TREE, outVec, inVec, sizeof(in4[0]));

  for (uint64_t i = 0; i < 4; ++i) {
    memcpy(out_seed[i], out4[i], PARAM_seed_size);
    memcpy(out_rho[i], &out4[i][PARAM_seed_size], PARAM_rho_size);
  }
}
