/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef HASH_SHAKE_H
#define HASH_SHAKE_H

#include <stdint.h>
#include <stdio.h>

#include "macros.h"
#include "endian_compat.h"

#include <fips202.h>
#define OQS_SHA3_shake128_inc_ctx shake128incctx
#define OQS_SHA3_shake128_inc_init shake128_inc_init
#define OQS_SHA3_shake128_inc_absorb shake128_inc_absorb
#define OQS_SHA3_shake128_inc_finalize shake128_inc_finalize
#define OQS_SHA3_shake128_inc_squeeze shake128_inc_squeeze
#define OQS_SHA3_shake128_inc_ctx_release shake128_inc_ctx_release
#define OQS_SHA3_shake256_inc_ctx shake256incctx
#define OQS_SHA3_shake256_inc_init shake256_inc_init
#define OQS_SHA3_shake256_inc_absorb shake256_inc_absorb
#define OQS_SHA3_shake256_inc_finalize shake256_inc_finalize
#define OQS_SHA3_shake256_inc_squeeze shake256_inc_squeeze
#define OQS_SHA3_shake256_inc_ctx_release shake256_inc_ctx_release

typedef struct hash_context_oqs_s {
  union {
    OQS_SHA3_shake128_inc_ctx shake128_ctx;
    OQS_SHA3_shake256_inc_ctx shake256_ctx;
  };
  unsigned char shake256;
} hash_context;

/**
 * Initialize hash context based on the security parameter. If the security parameter is 128,
 * SHAKE128 is used, otherwise SHAKE256 is used.
 */
static inline void hash_init(hash_context* ctx, unsigned int security_param) {
  if (security_param == 128) {
    OQS_SHA3_shake128_inc_init(&ctx->shake128_ctx);
    ctx->shake256 = 0;
  } else {
    OQS_SHA3_shake256_inc_init(&ctx->shake256_ctx);
    ctx->shake256 = 1;
  }
}

static inline void hash_update(hash_context* ctx, const uint8_t* data, size_t size) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_absorb(&ctx->shake256_ctx, data, size);
  } else {
    OQS_SHA3_shake128_inc_absorb(&ctx->shake128_ctx, data, size);
  }
}

static inline void hash_final(hash_context* ctx) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_finalize(&ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_inc_finalize(&ctx->shake128_ctx);
  }
}

static inline void hash_squeeze(hash_context* ctx, uint8_t* buffer, size_t buflen) {
  if (ctx->shake256) {
    OQS_SHA3_shake256_inc_squeeze(buffer, buflen, &ctx->shake256_ctx);
  } else {
    OQS_SHA3_shake128_inc_squeeze(buffer, buflen, &ctx->shake128_ctx);
  }
}

#endif
