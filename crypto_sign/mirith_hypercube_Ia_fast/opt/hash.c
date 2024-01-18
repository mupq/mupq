/*
 * Copyright 2023 Carlo Sanna, Javier Verbel, and Floyd Zweydinger.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fips202.h>
#include <string.h>
#include <stdlib.h>
#include "hash.h"

void hash_init(hash_ctx_t *ctx)
{
#if HASH_SIZE == 32
	sha3_256_inc_init(ctx);
#elif HASH_SIZE == 48
	sha3_384_inc_init(ctx);
#elif HASH_SIZE == 64
	sha3_512_inc_init(ctx);
#else
#error "HASH_SIZE not implemented!"
#endif
}

void hash_update(hash_ctx_t *ctx, const void *data, size_t length)
{
#if HASH_SIZE == 32
	sha3_256_inc_absorb(ctx, data, length);
#elif HASH_SIZE == 48
	sha3_384_inc_absorb(ctx, data, length);
#elif HASH_SIZE == 64
	sha3_512_inc_absorb(ctx, data, length);
#else
#error "HASH_SIZE not implemented!"
#endif
}

void hash_finalize(hash_ctx_t *ctx, hash_t hash)
{
#if HASH_SIZE == 32
	sha3_256_inc_finalize(hash, ctx);
#elif HASH_SIZE == 48
	sha3_384_inc_finalize(hash, ctx);
#elif HASH_SIZE == 64
	sha3_512_inc_finalize(hash, ctx);
#else
#error "HASH_SIZE not implemented!"
#endif
}

int hash_equal(hash_t hash1, hash_t hash2)
{
    return memcmp(hash1, hash2, HASH_SIZE) == 0;
}

void hash_digest0(hash_t hash, const hash_t salt, int l, int i, const seed_t seed)
{
    hash_ctx_t ctx;

    hash_init(&ctx);
    hash_update(&ctx, salt, HASH_SIZE);
    hash_update(&ctx, &l, sizeof(l));
    hash_update(&ctx, &i, sizeof(i));
    hash_update(&ctx, seed, SEED_SIZE);
    hash_finalize(&ctx, hash);
}

void hash_digest0_aux(hash_t hash,
    const hash_t salt,
    int l,
    int i,
    const seed_t seed,
    const ff_t a_aux[matrix_bytes_size(PAR_K, 1)],
    const ff_t K_aux[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    const ff_t C_aux[matrix_bytes_size(PAR_S, PAR_N - PAR_R)])
{
    hash_ctx_t ctx;

    hash_init(&ctx);

    hash_update(&ctx, salt, HASH_SIZE);
    hash_update(&ctx, &l, sizeof(l));
    hash_update(&ctx, &i, sizeof(i));
    hash_update(&ctx, seed, SEED_SIZE);
    hash_update(&ctx, a_aux, matrix_bytes_size(PAR_K, 1));
    hash_update(&ctx, K_aux, matrix_bytes_size(PAR_R, PAR_N - PAR_R));
    hash_update(&ctx, C_aux, matrix_bytes_size(PAR_S, PAR_N - PAR_R));
	hash_finalize(&ctx, hash);
}
