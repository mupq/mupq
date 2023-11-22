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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "blas_m4.h"
#include "matrix.h"
#include "matrix_constants.h"


void matrix_init_zero(ff_t *matrix, int n_rows, int n_cols)
{
    memset(matrix, 0, matrix_bytes_size(n_rows, n_cols));
}

ff_t matrix_get_entry(const ff_t *matrix, int n_rows, int i, int j)
{
    int nbytes_col;
    nbytes_col = matrix_bytes_per_column(n_rows);
    if (i & 1) // i is odd
    {
        return  matrix[nbytes_col * j + (i >> 1)] >> 4;
    }
    else
    {
        return matrix[nbytes_col * j + (i >> 1)] & 0x0f;
    }
}

void matrix_set_entry(ff_t *matrix, int n_rows, int i, int j, ff_t scalar)
{
    int nbytes_col, target_byte_id;
    nbytes_col = matrix_bytes_per_column(n_rows);
    target_byte_id = nbytes_col * j + (i >> 1);
    if (i & 1) // i is odd
    {
        matrix[target_byte_id] &= 0x0f;
        matrix[target_byte_id] |= (scalar << 4);
    }
    else
    {
        matrix[target_byte_id] &= 0xf0;
        matrix[target_byte_id] |= scalar;
    }
}

void matrix_init_random(ff_t *matrix, int n_rows, int n_cols, prng_t *prng)
{
    int matrix_height, matrix_height_x, i;

    matrix_height =  matrix_bytes_per_column(n_rows);
    matrix_height_x = matrix_height -  1;
    
    /* i = 0. */
    for (i = 0; i < n_cols; i++)
    {
        prng_bytes(prng, matrix + i * matrix_height , matrix_height);
    }
    
    if (n_rows & 1)
    {
        for (i = 0; i < n_cols; i++)
        {
            matrix[i * matrix_height + matrix_height_x ] &= 0x0f;
        }
    }
}

void matrix_copy(ff_t *matrix1, const ff_t *matrix2, int n_rows, int n_cols)
{
    int n_bytes;

    n_bytes = matrix_bytes_size(n_rows, n_cols);
    memcpy(matrix1, matrix2, n_bytes);
}

void matrix_negate(ff_t *matrix, int n_rows, int n_cols)
{
    /* Nothing to do in characteristic 2. */

    /* Suppress 'unused parameter' warnings. */
    (void)(matrix); (void)(n_rows); (void)(n_cols);
}

void matrix_add(ff_t *matrix1, const ff_t *matrix2, int n_rows, int n_cols)
{
    const unsigned n_bytes = matrix_bytes_size(n_rows, n_cols);
    gf256v_add_neon(matrix1, matrix2, n_bytes);
}

void matrix_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    int n_rows, int n_cols)
{
    const int n_bytes = matrix_bytes_size(n_rows, n_cols);

    gf16mat_intmat_mul_ref(matrix1, matrix2, n_bytes, scalar);
    return;

    for (uint32_t i = 0; i < n_bytes; ++i) {
        matrix1[i] ^= mult_table[(matrix2[i]&0xf) + 16 * scalar];
        matrix1[i] ^= mult_table[(matrix2[i]>>4)  + 16 * scalar] << 4;
    }
}


void matrix_add_product(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
                 int n_rows1, int n_cols1, int n_cols2) {

    const uint32_t n_bytes_per_column1 = matrix_bytes_per_column(n_rows1);
    gf16mat_rowmat_mul_ref(matrix1, matrix2, n_cols2, n_bytes_per_column1, matrix3, 2*((n_cols1+1)>>1));
    return;

    const uint32_t matrix_height =  matrix_bytes_per_column(n_rows1);
    const uint32_t matrix_height_x = matrix_height -  1;

    for (uint32_t  i = 0; i < n_rows1; i++) {
        for (uint32_t  j = 0; j < n_cols2; j++) {
            ff_t entry_i_j = 0;

            for (uint32_t  k = 0; k < n_cols1; k++) {
                uint32_t entry_i_k = matrix_get_entry(matrix2, n_rows1, i, k);
                uint32_t entry_k_j = matrix_get_entry(matrix3, n_cols1, k, j);
                entry_i_j ^= mult_table[entry_i_k + 16 * entry_k_j];
            }

            matrix1[matrix_height * j + (i >> 1)] ^= (entry_i_j << (4* (i&1)));
        }
    }

    if (n_rows1 & 1) {
        for (uint32_t i = 0; i < n_cols2; i++) {
            matrix1[i * matrix_height + matrix_height_x] &= 0x0f;
        }
    }
}

void matrix_subtract_product(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
                 int n_rows1, int n_cols1, int n_cols2) {

    /* This works because we are in characteristic 2. */
    matrix_add_product(matrix1, matrix2, matrix3, n_rows1, n_cols1, n_cols2);
}

void matrix_subtract(ff_t *matrix1, const ff_t *matrix2, int n_rows, int n_cols)
{
    matrix_add(matrix1, matrix2, n_rows, n_cols);
}

void matrix_subtract_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    int n_rows, int n_cols)
{
    matrix_add_multiple(matrix1, scalar, matrix2, n_rows, n_cols);
}

void matrix_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    int n_rows1, int n_cols1, int n_cols2)
{
    const uint32_t n_bytes_per_column1 = matrix_bytes_per_column(n_rows1);
    gf16mat_colmat_mul(result, matrix1, n_bytes_per_column1, n_cols1, matrix2, n_cols2);
    return;

    const uint32_t matrix_height =  matrix_bytes_per_column(n_rows1);
    const uint32_t matrix_height_x = matrix_height -  1;

    for (uint32_t  i = 0; i < n_rows1; i++) {
        for (uint32_t  j = 0; j < n_cols2; j++) {
            ff_t entry_i_j = 0;

            for (uint32_t  k = 0; k < n_cols1; k++) {
                uint32_t entry_i_k = matrix_get_entry(matrix1, n_rows1, i, k);
                uint32_t entry_k_j = matrix_get_entry(matrix2, n_cols1, k, j);
                entry_i_j ^= mult_table[entry_i_k + 16 * entry_k_j];
            }

            matrix_set_entry(result, n_rows1, i, j, entry_i_j);
        }
    }

    if (n_rows1 & 1) {
        for (uint32_t i = 0; i < n_cols2; i++) {
            result[i * matrix_height + matrix_height_x] &= 0x0f;
        }
    }
}

void matrix_horizontal_concatenation(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    int n_rows, int n_cols1, int n_cols2)
{
    int n_bytes1, n_bytes2;

    n_bytes1 = matrix_bytes_size(n_rows, n_cols1);
    n_bytes2 = matrix_bytes_size(n_rows, n_cols2);

    memcpy(result, matrix1, n_bytes1);
    memcpy(result + n_bytes1, matrix2, n_bytes2);
}

void matrix_horizontal_split(ff_t *matrix1, ff_t *matrix2, const ff_t *matrix,
    int n_rows, int n_cols1, int n_cols2)
{
    int n_bytes1, n_bytes2;

    n_bytes1 = matrix_bytes_size(n_rows, n_cols1);
    n_bytes2 = matrix_bytes_size(n_rows, n_cols2);
    
    memcpy(matrix1, matrix, n_bytes1);
    memcpy(matrix2, matrix + n_bytes1, n_bytes2);
}

void _matrix_pack_nrows_even(uint8_t **dest, int *bit_offset, const ff_t *matrix,
                             int n_rows, int n_cols)
{

    /* the packing is done row-wise */
    int bo, n_bytes;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    n_bytes = matrix_bytes_size(n_rows, n_cols);

    if (bo)
    {
        /* Pack last entry of matrix (in column-major order) in the higher bits of dest[0]. */
        ((uint8_t *)*dest)[0] |= matrix[n_bytes - 1] & 0xf0;

        /* Pack all the bytes in matrix except the last one */
        memcpy(&(((uint8_t *)*dest)[1]), matrix, n_bytes - 1); /* TODO: improve with avx2. */

        /* Pack the second-to-last entry of matrix. */
        ((uint8_t *)*dest)[n_bytes] = matrix[n_bytes - 1] & 0x0f;
    }
    else
    {
        memcpy((uint8_t *)*dest, matrix, n_bytes);
    }

    *dest = &(((uint8_t *)*dest)[n_bytes]);
}

void _matrix_unpack_nrows_even(ff_t *matrix, uint8_t **source, int *bit_offset,
                               int n_rows, int n_cols)
{
    int bo, n_bytes;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    n_bytes = matrix_bytes_size(n_rows, n_cols);

    if (bo)
    {
        /* Unpack all the bytes in matrix except the last one */
        memcpy(matrix, &(((uint8_t *)*source)[1]), n_bytes - 1); /* unpack all the bytes but the last one. */

        /* Unpack the last two entries of matrix. */
        matrix[n_bytes - 1] = (((uint8_t *)*source)[n_bytes] & 0x0f) | (((uint8_t *)*source)[0] & 0xf0);
    }
    else
    {
        memcpy(matrix, (uint8_t *)*source, n_bytes);
    }

    *source = &(((uint8_t *)*source)[n_bytes]);
}


/* Remove the last row of matrix and append it to matrix as additional column(s) */
void _matrix_pack_nrows_odd(uint8_t **dest, int *bit_offset, const ff_t *matrix,
                            int n_rows, int n_cols)
{
    assert((n_rows & 1) == 1);

    int j, bo, next_bo, matrix_height, matrix_height_x, n_bytes_not_in_last_row, n_bytes;
    int ad_bytes, jump_nbytes;
    uint8_t row_entry_j, row_entry_j_1;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }

    matrix_height = (n_rows >> 1) + 1;
    matrix_height_x =  matrix_height - 1;
    n_bytes_not_in_last_row = matrix_height_x * n_cols;
    n_bytes = matrix_bytes_size(n_rows, n_cols);

    /* Bytes that are not part of the last row. */
    for (j = 0; j < n_cols; j++)
    {
        memcpy(&(((uint8_t *)*dest)[bo + j * matrix_height_x]), &matrix[matrix_height * j], matrix_height_x);
    }
    /* When n_cols is odd the maximum value of j is j_max = n_cols - 3, hence j_max + 1 = n_cols - 2.
     * Hence in the following loop wont add the entry n_cols - 1 (the last entry) of the last row.
     * When n_cols is even the maximum value of j is j_max = n_cols - 4, hence j_max + 1 = n_cols - 3.
     * Hence in the following loop wont add the entries n_cols - 2 and n_cols - 1 (the last entry) of the last row. */
    ad_bytes = bo;
    for (j = 0; j < n_cols - 2; j+=2)
    {
        row_entry_j = matrix[matrix_height * j + matrix_height_x] & 0x0f; /* j-th entry of the last row. */
        row_entry_j_1 = matrix[matrix_height * (j + 1) + matrix_height_x] & 0x0f; /* (j + 1)-th entry of the last row. */
        ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes]  =  (row_entry_j_1 << 4) | row_entry_j;
        ad_bytes +=1;
    }
    /* When the is an odd number of columns and
     * bit_off_set = 1, we locate the last entry of matrix in higher bits
     * of the first byte of the source. Otherwise, if bit_off_set = 0, we locate ast entry of matrix
     * the next byte of dest. */
    if  (bo)
    {
        ((uint8_t *)*dest)[0]&= 0x0f;
        ((uint8_t *)*dest)[0] |= (matrix[n_bytes - 1] << 4);

        if ((n_cols & 1) == 0) /* case n_cols is even. */
        {
            /* Packing the second-last entry of the last row. */
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = matrix[n_bytes - matrix_height - 1] & 0x0f;

        }
    }
    else
    {
        /* If n_cols is even and bo = 0,
         * we pack the entries n_cols - 2 and n_cols - 1
         * in the last row in the byte
         * of the current local buffer. */
        if ((n_cols & 1) == 0)
        {
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = (matrix[n_bytes - 1] << 4) | matrix[n_bytes - matrix_height - 1];
        }
            /* Odd number of columns case. */
            /* In this case, we locate the last entry the next byte of dest. */
        else
        {
            ((uint8_t *)*dest)[n_bytes_not_in_last_row + ad_bytes] = matrix[n_bytes - 1] & 0x0f;
        }
    }

    jump_nbytes = matrix_height_x * n_cols + (n_cols >> 1);

    if (bo)
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 1;
        }
        else
        {
            next_bo = 0;
            jump_nbytes += (n_cols & 1);
        }
    }
    else
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 0;
        }
        else
        {
            next_bo = 1;
        }
    }

    if (bit_offset != NULL)
    {
        *bit_offset = next_bo;
    }

    *dest = &(((uint8_t *)*dest)[jump_nbytes]);
}

void _matrix_unpack_nrows_odd(ff_t *matrix, uint8_t **source, int *bit_offset,
                              int n_rows, int n_cols)
{
    assert((n_rows & 1) == 1);

    int j, bo, next_bo, matrix_height, matrix_height_x, n_bytes_not_in_last_row, n_bytes;
    int ad_bytes, jump_nbytes;
    uint8_t row_entries_j_and_j_1;

    if (bit_offset != NULL)
    {
        bo = *bit_offset;
    }
    else
    {
        bo = 0;
    }


    matrix_height = (n_rows >> 1) + 1;
    matrix_height_x =  matrix_height - 1;
    n_bytes_not_in_last_row = matrix_height_x * n_cols;
    n_bytes = matrix_bytes_size(n_rows, n_cols);


    /* Bytes that are not part of the last row. */
    for (j = 0; j < n_cols; j++)
    {
        memcpy(&matrix[matrix_height * j], &(((uint8_t *)*source)[bo + j * matrix_height_x]), matrix_height_x);
    }

    /* When n_cols is odd the maximum value of j is j_max = n_cols - 3, hence j_max + 1 = n_cols - 2.
     * Hence in the following loop wont add the entry n_cols - 1 (the last entry) of the last row .
     * When n_cols is even the maximum value of j is j_max = n_cols - 4, hence j_max + 1 = n_cols - 3.
     * Hence in the following loop wont add the entries n_cols - 2 and n_cols - 1 (the last entry) of the last row. */
    ad_bytes = bo;
    for (j = 0; j < n_cols - 2; j+=2)
    {
        row_entries_j_and_j_1 = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes];
        matrix[matrix_height * j + matrix_height_x] = row_entries_j_and_j_1 & 0x0f;
        matrix[matrix_height * (j + 1) + matrix_height_x] =  row_entries_j_and_j_1 >> 4;
        ad_bytes +=1;
    }
    /* When the is an odd number of columns and
     * bit_off_set = 1, we locate the last entry of matrix in higher bits
     * of the first byte of the source. Otherwise, if bit_off_set = 0, we locate ast entry of matrix
     * the next byte of dest. */
    if  (bo)
    {
        matrix[n_bytes - 1]  = ((uint8_t *)*source)[0] >> 4;

        if ((n_cols & 1) == 0) /* case n_rows is even. */
        {
            matrix[n_bytes - matrix_height - 1] = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes] & 0x0f;
        }
    }
    else
    {
        /* If n_cols is even and bo = 0,
        * we unpack the last row in the byte of the current local buffer
        * into the entries n_cols - 2 and n_cols - 1
        * of the last row of matrix .*/
        if ((n_cols & 1) == 0)
        {
            row_entries_j_and_j_1 =  ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes];
            matrix[n_bytes - 1]  = row_entries_j_and_j_1 >> 4;
            matrix[n_bytes - matrix_height - 1] = row_entries_j_and_j_1 & 0x0f;
        }
            /* Odd number of columns case.
             * In this case, we locate the last entry the next byte of dest. */
        else
        {
            matrix[n_bytes -1 ] = ((uint8_t *)*source)[n_bytes_not_in_last_row + ad_bytes] & 0x0f;
        }

    }

    jump_nbytes = matrix_height_x * n_cols + (n_cols >> 1);

    if (bo)
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 1;
        }
        else
        {
            next_bo = 0;
            jump_nbytes +=  (n_cols & 1);
        }
    }
    else
    {
        if ((n_cols & 1) == 0)
        {
            next_bo = 0;
        }
        else
        {
            next_bo = 1;
        }
    }

    if (bit_offset != NULL)
    {
        *bit_offset = next_bo;
    }

    *source = &(((uint8_t *)*source)[jump_nbytes]);
}

void matrix_pack(uint8_t **dest, int *bit_offset, const ff_t *matrix,
                 int n_rows, int n_cols)
{
    /* An even number of rows. */
    if ((n_rows & 1) == 0)
    {
        _matrix_pack_nrows_even(dest, bit_offset, matrix, n_rows, n_cols);
    }
    else
    {
        _matrix_pack_nrows_odd(dest, bit_offset, matrix, n_rows, n_cols);
    }

}

void matrix_unpack(ff_t *matrix, uint8_t **source, int *bit_offset,
                   int n_rows, int n_cols)
{
    /* An even number of rows. */
    if ((n_rows & 1) == 0)
    {
        _matrix_unpack_nrows_even(matrix, source, bit_offset, n_rows, n_cols);
    }
    else
    {
        _matrix_unpack_nrows_odd(matrix, source, bit_offset, n_rows, n_cols);
    }
}


