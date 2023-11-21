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

#elif MIRITH_MODE == 2 /* Ia, shorter. */

#define CRYPTO_PUBLICKEYBYTES 129
#define CRYPTO_SECRETKEYBYTES 145
#define CRYPTO_BYTES 5036

#elif MIRITH_MODE == 3 /* Ia, shortest. */

#define CRYPTO_PUBLICKEYBYTES 129
#define CRYPTO_SECRETKEYBYTES 145
#define CRYPTO_BYTES 4536

#elif MIRITH_MODE == 4 /* Ib, fast. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 9105

#elif MIRITH_MODE == 5 /* Ib, short. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 6309

#elif MIRITH_MODE == 6 /* Ib, shorter. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 5491

#elif MIRITH_MODE == 7 /* Ib, shortest. */

#define CRYPTO_PUBLICKEYBYTES 144
#define CRYPTO_SECRETKEYBYTES 160
#define CRYPTO_BYTES 4886

#elif MIRITH_MODE == 8 /* IIIa, fast. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 17139

#elif MIRITH_MODE == 9 /* IIIa, short. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 12440

#elif MIRITH_MODE == 10 /* IIIa, shorter. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 10746

#elif MIRITH_MODE == 11 /* IIIa, shortest. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 9954

#elif MIRITH_MODE == 12 /* IIIb, fast. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 18459

#elif MIRITH_MODE == 13 /* IIIb, short. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 13136

#elif MIRITH_MODE == 14 /* IIIb, shorter. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 11202

#elif MIRITH_MODE == 15 /* IIIb, shortest. */

#define CRYPTO_PUBLICKEYBYTES 205
#define CRYPTO_SECRETKEYBYTES 229
#define CRYPTO_BYTES 10314

#elif MIRITH_MODE == 16 /* Va, fast. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 31468

#elif MIRITH_MODE == 17 /* Va, short. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 21795

#elif MIRITH_MODE == 18 /* Va, shorter. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 19393

#elif MIRITH_MODE == 19 /* Va, shortest. */

#define CRYPTO_PUBLICKEYBYTES 253
#define CRYPTO_SECRETKEYBYTES 285
#define CRYPTO_BYTES 17522

#elif MIRITH_MODE == 20 /* Vb, fast. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 34059

#elif MIRITH_MODE == 21 /* Vb, short. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 23182

#elif MIRITH_MODE == 22 /* Vb, shorter. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 20394

#elif MIRITH_MODE == 23 /* Vb, shortest. */

#define CRYPTO_PUBLICKEYBYTES 274
#define CRYPTO_SECRETKEYBYTES 306
#define CRYPTO_BYTES 18292

#else

#error "MIRITH_MODE not implemented!"

#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#endif
