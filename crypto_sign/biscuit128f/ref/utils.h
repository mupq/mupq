#ifndef _UTILS_H_
#define _UTILS_H_

#include <stdint.h>

// #include "sha3/KeccakHash.h"
#include <fips202.h>

// typedef Keccak_HashInstance h_ctx_t[1];
typedef shake256incctx h_ctx_t[1];

void
commit (uint8_t *output, const uint8_t *salt, uint16_t e, uint16_t i,
        uint8_t *seed, int seclevel);

void
expand_init (h_ctx_t h_ctx, const uint8_t *input, size_t inlen);

void
expandtape_init (h_ctx_t h_ctx, const uint8_t *salt, uint16_t e, uint16_t i,
                 uint8_t *seed, int seclevel);

void
expand_extract (uint8_t *output, size_t outlen, h_ctx_t h_ctx);

void
prf_init (h_ctx_t h_ctx, const uint8_t *input, size_t inlen);

void
prf_update (h_ctx_t h_ctx, const uint8_t *input, size_t inlen);

void
prf_ready (h_ctx_t h_ctx);

void
prf_generate (uint8_t *output, size_t outlen, h_ctx_t h_ctx);

void
H1_init (h_ctx_t h_ctx, const uint8_t *salt, const uint8_t *msg, size_t msglen,
         int seclevel);

void
H1_update (h_ctx_t h_ctx, const uint8_t *sigma1, size_t sigma1len);

void
H1_final (uint8_t *output, h_ctx_t h_ctx, int seclevel);

void
H2_init (h_ctx_t h_ctx, const uint8_t *salt, const uint8_t *h1, int seclevel);

void
H2_update (h_ctx_t h_ctx, const uint8_t *sigma2, size_t sigma2len);

void
H2_final (uint8_t *output, h_ctx_t h_ctx, int seclevel);

int
ilog2 (int N);

void
get_seeds (uint8_t *output, const uint8_t *seed, const uint8_t *salt,
           uint16_t e, uint16_t N, int seclevel);

void
get_path (uint8_t *output, const uint8_t *seed, const uint8_t *salt,
          uint16_t e, uint16_t ibar, uint16_t N, int seclevel);

void
get_path_seed (uint8_t *output, const uint8_t *path, const uint8_t *salt,
               uint16_t e, uint16_t i, uint16_t ibar, uint16_t N,
               int seclevel);

void
get_path_seeds (uint8_t *output, const uint8_t *path, const uint8_t *salt,
                uint16_t e, uint16_t ibar, uint16_t N, int seclevel);

#endif
