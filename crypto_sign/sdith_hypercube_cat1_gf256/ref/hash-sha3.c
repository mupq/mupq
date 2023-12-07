#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>

#include "param.h"
#include "rng.h"
#include "fips202.h"



void sdith_hash_create_hash_ctx(HASH_CTX *ctx, uint8_t prefix) {
  // TODO other sec levels
#if defined(CAT_1)
  //Keccak_HashInitialize_SHA3_256(inst);
  sha3_256_inc_init(ctx);
#elif defined(CAT_3)
  Keccak_HashInitialize_SHA3_384(ctx);
#else
  Keccak_HashInitialize_SHA3_512(ctx);
#endif
  sha3_256_inc_absorb(ctx, &prefix, 1);
  //Keccak_HashUpdate(inst, &prefix, sizeof(uint8_t) << 3);
}


void sdith_hash_free_hash_ctx(HASH_CTX *ctx) { sha3_256_inc_ctx_release(ctx); }

void sdith_hash_digest_update(HASH_CTX *ctx, void const *in, int inBytes) {
  //Keccak_HashUpdate((Keccak_HashInstance *)ctx, (uint8_t const *)in,
  //                  inBytes << 3);
  sha3_256_inc_absorb(ctx, in, inBytes);
}

void sdith_hash_final(HASH_CTX *ctx, void *dest) {
  //Keccak_HashFinal((Keccak_HashInstance *)ctx, dest);
  sha3_256_inc_finalize(dest, ctx);
}

void sdith_hash(uint8_t prefix, void *dest, void const *data, int dataBytes) {
  HASH_CTX ctx;
  sdith_hash_create_hash_ctx(&ctx, prefix);
  sdith_hash_digest_update(&ctx, data, dataBytes);
  sdith_hash_final(&ctx, dest);
  sdith_hash_free_hash_ctx(&ctx);
}


void sdith_hash_create_hash4_ctx(HASH4_CTX *ctx,uint8_t prefix) {
#if defined(CAT_1)
  sha3_256_inc_init(&ctx->states[0]);
  sha3_256_inc_init(&ctx->states[1]);
  sha3_256_inc_init(&ctx->states[2]);
  sha3_256_inc_init(&ctx->states[3]);
#elif defined(CAT_3)
  Keccak_HashInitializetimes4_SHA3_384(&hash4_ctx);
#else
  Keccak_HashInitializetimes4_SHA3_512(&hash4_ctx);
#endif
  sha3_256_inc_absorb(&ctx->states[0], &prefix, 1);
  sha3_256_inc_absorb(&ctx->states[1], &prefix, 1);
  sha3_256_inc_absorb(&ctx->states[2], &prefix, 1);
  sha3_256_inc_absorb(&ctx->states[3], &prefix, 1);
  //const uint8_t *prefix_ptr[4] = {&prefix, &prefix, &prefix, &prefix};
  //Keccak_HashUpdatetimes4(&hash4_ctx, prefix_ptr, sizeof(prefix) << 3);
}

void sdith_hash_free_hash4_ctx(HASH4_CTX *ctx) {
  sha3_256_inc_ctx_release(&ctx->states[0]);
  sha3_256_inc_ctx_release(&ctx->states[1]);
  sha3_256_inc_ctx_release(&ctx->states[2]);
  sha3_256_inc_ctx_release(&ctx->states[3]);
}

void sdith_hash4_digest_update(HASH4_CTX *ctx, void **in, int inBytes) {
  sha3_256_inc_absorb(&ctx->states[0], in[0], inBytes);
  sha3_256_inc_absorb(&ctx->states[1], in[1], inBytes);
  sha3_256_inc_absorb(&ctx->states[2], in[2], inBytes);
  sha3_256_inc_absorb(&ctx->states[3], in[3], inBytes);

  //Keccak_HashUpdatetimes4((Keccak_HashInstancetimes4 *)ctx,
  //                        (const uint8_t **)in, inBytes << 3);
}

void sdith_hash4_final(HASH4_CTX *ctx, void **dest) {
  //Keccak_HashFinaltimes4((Keccak_HashInstancetimes4 *)ctx, (uint8_t **)dest);
  sha3_256_inc_finalize(dest[0], &ctx->states[0]);
  sha3_256_inc_finalize(dest[1], &ctx->states[1]);
  sha3_256_inc_finalize(dest[2], &ctx->states[2]);
  sha3_256_inc_finalize(dest[3], &ctx->states[3]);
}

void sdith_hash4(uint8_t prefix, void **dest, void **data, int dataBytes) {
  HASH4_CTX ctx;
  sdith_hash_create_hash4_ctx(&ctx, prefix);
  sdith_hash4_digest_update(&ctx, data, dataBytes);
  sdith_hash4_final(&ctx, dest);
  sdith_hash_free_hash4_ctx(&ctx);
}
