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

#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include "hash_types.h"
#include "matrix.h"
#include "params.h"

/* Type for a hashing context. */
// typedef void *hash_ctx_t;
#if HASH_SIZE == 32
	typedef sha3_256incctx hash_ctx_t;
#elif HASH_SIZE == 48
	typedef sha3_384incctx hash_ctx_t;
#elif HASH_SIZE == 64
	typedef sha3_512incctx hash_ctx_t;
#else
#error "HASH_SIZE not implemented!"
#endif

/* Initialize the hashing context 'ctx'. */
void hash_init(hash_ctx_t *ctx);

/* Update the hashing context 'ctx' with 'data'. */
void hash_update(hash_ctx_t *ctx, const void *data, size_t length);

/* Finalize the hashing, write the digest over 'hash'. */
void hash_finalize(hash_ctx_t *ctx, hash_t hash);

/* Return 'True' if 'hash1' and 'hash2' are equal, and 'False' otherwise. */
int hash_equal(hash_t hash1, hash_t hash2);

/* Write over 'hash' the hash digest of 'salt', 'l', 'i', 'seed'. */
void hash_digest0(hash_t hash, const hash_t salt, int l, int i, const seed_t seed);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Write over 'hash' the hash digest of 'salt', 'l', 'i', 'seed',
 * 'a_aux', 'K_aux', 'C_aux'. */
void hash_digest0_aux(hash_t hash,
    const hash_t salt,
    int l,
    int i,
    const seed_t seed,
    const ff_t a_aux[matrix_bytes_size(PAR_K, 1)],
    const ff_t K_aux[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    const ff_t C_aux[matrix_bytes_size(PAR_S, PAR_N - PAR_R)]);

#endif
