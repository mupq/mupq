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

#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

/* Mode * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#if MIRITH_MODE == 0 /* Ia, fast. */

#define PAR_Q 16
#define PAR_M 15
#define PAR_N 15
#define PAR_K 78
#define PAR_R 6
#define PAR_S 5
#define TAU 39
#define N_PARTIES 2
#define D_DIMENSION 4
#define N_PARTIES_ROUND 16
#define SEED_SIZE 16
#define HASH_SIZE 32
#define TREE_HEIGHT 4
#define TREE_N_NODES 31
#define CRYPTO_PUBLICKEYBYTES 129
#define CRYPTO_SECRETKEYBYTES 145
#define CRYPTO_BYTES 7877

#elif MIRITH_MODE == 1 /* Ib, fast. */

#define PAR_Q 16
#define PAR_M 16
#define PAR_N 16
#define PAR_K 142
#define PAR_R 4
#define PAR_S 5
#define TAU 39
#define N_PARTIES 2
#define D_DIMENSION 4
#define N_PARTIES_ROUND 16
#define SEED_SIZE 16
#define HASH_SIZE 32
#define TREE_HEIGHT 4
#define TREE_N_NODES 31
#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 9105

#else

#error "MIRITH_MODE not implemented!"

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif
