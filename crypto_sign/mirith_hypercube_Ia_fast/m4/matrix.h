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

#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include "prng.h"

#define matrix_bytes_per_column(n_rows) ((n_rows >> 1u) + (n_rows & 1u))
#define matrix_bytes_size(n_rows, n_cols) ((matrix_bytes_per_column(n_rows)) * (n_cols))

/* Type for an element of the finite field. */
typedef uint8_t ff_t;

/* Return the number of bytes of a 'n_rows x n_cols' matrix. */
//int matrix_bytes_size(int n_rows, int n_cols);

/* Initialized 'matrix' with zero entries. */
void matrix_init_zero(ff_t *matrix, uint32_t n_rows, uint32_t n_cols);

/* Return the (i, j) entry of 'matrix'. */
ff_t matrix_get_entry(const ff_t *matrix, uint32_t n_rows, uint32_t i, uint32_t j);

/* Set the (i, j) entry of 'matrix' to be scalar.*/
void matrix_set_entry(ff_t *matrix, uint32_t n_rows, uint32_t i, uint32_t j, ff_t scalar);

/* Initialized 'matrix' with random entries. */
void matrix_init_random(ff_t *matrix, uint32_t n_rows, uint32_t n_cols, prng_t *prng);

/* Overwrite 'matrix1' with 'matrix2'. */
void matrix_copy(ff_t *matrix1, const ff_t *matrix2, uint32_t n_rows, uint32_t n_cols);

/* Replace 'matrix' with '-matrix'. */
void matrix_negate(ff_t *matrix, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 + matrix2'. */
void matrix_add(ff_t *matrix1, const ff_t *matrix2, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 + scalar * matrix2'. */
void matrix_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 - matrix2'. */
void matrix_subtract(ff_t *matrix1, const ff_t *matrix2, uint32_t n_rows, uint32_t n_cols);

/* Overwrite 'matrix1' with 'matrix1 - scalar * matrix2'. */
void matrix_subtract_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    uint32_t n_rows, uint32_t n_cols);

/* Write 'matrix1 * matrix2' over 'result'. */
void matrix_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    uint32_t n_rows1, uint32_t n_cols1, uint32_t n_cols2);

/* Write '[matrix1 | matrix2]' over 'result'. */
void matrix_horizontal_concatenation(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    uint32_t n_rows, uint32_t n_cols1, uint32_t n_cols2);

/* Split 'matrix' as 'matrix = [matrix1 | matrix2]. */
void matrix_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
    uint32_t n_rows, uint32_t n_cols1, uint32_t n_cols2);

/* Pack 'matrix' over 'dest (bytes) + bit_offset (bits)'.
 * Update 'dest' and 'bit_offset' for the next call of 'matrix_pack'.
 * If 'bit_offset == NULL', then an offset of 0 bits is used. */
void matrix_pack(uint8_t **dest, uint32_t *bit_offset, const ff_t *matrix,
    uint32_t n_rows, uint32_t n_cols);

/* Unpack 'matrix' from 'source (bytes) + bit_offset (bits)'.
 * Update 'source' and 'bit_offset' for the next call of 'matrix_unpack'.
 * If 'bit_offset == NULL', then an offset of 0 bits is used. */
void matrix_unpack(ff_t *matrix, uint8_t **source, uint32_t *bit_offset,
    uint32_t n_rows, uint32_t n_cols);

#endif
