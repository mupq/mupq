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

#ifndef PRNG_H
#define PRNG_H

#include <stdint.h>
#include "fips202.h"
#include "hash_types.h"
#include "params.h"

typedef shake256incctx prng_t;

/* Initialize 'prng' from 'salt' and 'seed'.
 * If 'salt == NULL' then 'salt' is ignored.
 * If 'seed == NULL' then 'seed' is ignored. */
void prng_init(prng_t *prng, const hash_t salt, const seed_t seed);

/* Write 'length' pseudorandom bytes over 'target',
 * update the internal state of 'prng'. */
void prng_bytes(prng_t *prng, void *target, size_t length);

#endif
