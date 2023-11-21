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

#ifndef API_H
#define API_H

#include "config.h"

#include "sign.h"

/* Mode * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#if MIRITH_MODE == 0 /* Ia, fast. */

#define CRYPTO_PUBLICKEYBYTES 129
#define CRYPTO_SECRETKEYBYTES 145
#define CRYPTO_BYTES 7877

#elif MIRITH_MODE == 1 /* Ia, short. */

#define CRYPTO_PUBLICKEYBYTES 129
#define CRYPTO_SECRETKEYBYTES 145
#define CRYPTO_BYTES 5673

#elif MIRITH_MODE == 2 /* Ib, fast. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 9105

#elif MIRITH_MODE == 3 /* Ib, short. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 6309

#elif MIRITH_MODE == 4 /* IIIa, fast. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 17139

#elif MIRITH_MODE == 5 /* IIIa, short. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 12440

#elif MIRITH_MODE == 6 /* IIIb, fast. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 18459

#elif MIRITH_MODE == 7 /* IIIb, short. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 13136

#elif MIRITH_MODE == 8 /* Va, fast. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 30458

#elif MIRITH_MODE == 9 /* Va, short. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 21795

#elif MIRITH_MODE == 10 /* Vb, fast. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 33048

#elif MIRITH_MODE == 11 /* Vb, short. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 23182

#else

#error "MIRITH_MODE not implemented!"

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif
