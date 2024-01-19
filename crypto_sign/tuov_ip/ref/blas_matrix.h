/// @file blas_matrix.h
/// @brief linear algebra functions for matrix opt.
///
#ifndef _BLAS_MATRIX_H_
#define _BLAS_MATRIX_H_

#include <stdint.h>



#ifdef  __cplusplus
extern  "C" {
#endif


///////////////// Section: multiplications  ////////////////////////////////
///// matrix-vector /////

/// @brief matrix-vector multiplication:  c = matA * b , in GF(16)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   b          - the vector b.
///
void gf16mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);


/// @brief matrix-vector multiplication:  c = matA * b , in GF(256)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   b          - the vector b.
///
void gf256mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);


/// @brief matrix-vector multiplication:  c = matA * b , in GF(16)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   multab_b     - the multiplication tables of the vector b.
///
void gf16mat_prod_multab(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *multab_b);


/// @brief matrix-vector multiplication:  c = matA * b , in GF(256)
///
/// @param[out]  c         - the output vector c
/// @param[in]   matA      - a column-major matrix A.
/// @param[in]   n_A_vec_byte  - the size of column vectors in bytes.
/// @param[in]   n_A_width   - the width of matrix A.
/// @param[in]   multab_b     - the multiplication tabls of the vector b.
///
void gf256mat_prod_multab(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *multab_b);


/// @brief F is an batched-matrix, lambda is a vector with batch/2 elements,
///         we want to compute M = lambda_0 * F_0 + ... + lambda_{batch/2-1} * F_{batch/2-1}. in GF(16)
/// @param[out]  M               - a col-major F_rows * F_cols matrix             -
/// @param[in]   F                - a batched row-major F_rows * F_cols matrix F
/// @param[in]   F_rows           - number of rows of F
/// @param[in]   F_cols           - number of cols of F
/// @param[in]   lambda           - the coefficient vector, with batch_size/2 elements
/// @param[in]   batch            - batch size of F
///
void gf16mat_linearcomb_half(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch);


/// @brief F is an batched-matrix, lambda is a vector with batch/2 elements,
///         we want to compute M = lambda_0 * F_0 + ... + lambda_{batch/2-1} * F_{batch/2-1}. in GF(256)
/// @param[out]  M               - a col-major F_rows * F_cols matrix             -
/// @param[in]   F                - a batched row-major F_rows * F_cols matrix F
/// @param[in]   F_rows           - number of rows of F
/// @param[in]   F_cols           - number of cols of F
/// @param[in]   lambda           - the coefficient vector, with batch_size/2 elements
/// @param[in]   batch            - batch size of F
///
void gf256mat_linearcomb_half(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch);




//////////////////////////  Gaussian for solving lienar equations ///////////////////////////


/// @brief Computing the row echelon form of a matrix, in GF(16)
///
/// @param[in,out]  sq_col_mat_a   - square matrix parts of a linear system. a is a column major matrix.
///                                  The returned matrix of row echelon form is a row major matrix.
/// @param[in,out]  constant       - constant parts of a linear system.
/// @param[in]           len       - the width of the matrix a, i.e., the number of column vectors.
/// @return   1(true) if success. 0(false) if the matrix is singular.
///
unsigned gf16mat_gaussian_elim(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

/// @brief Back substitution of the constant terms with a row echelon form of a matrix, in GF(16)
///
/// @param[in,out]  constant       - constant parts of a linear system.
/// @param[in]     sq_row_mat_a    - row echelon form of a linear system.
/// @param[in]           len       - the height of the matrix a, i.e., the number of row vectors.
///
void gf16mat_back_substitute( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);


unsigned gf16mat_gaussian_elim_unde(unsigned char *sqmat_a, unsigned A_rows, unsigned A_cols, unsigned char *constant);

void gf16mat_back_substitute_unde(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col);


/// @brief Computing the row echelon form of a matrix, in GF(256)
///
/// @param[in,out]  sq_col_mat_a   - square matrix parts of a linear system. a is a column major matrix.
///                                  The returned matrix of row echelon form is a row major matrix.
/// @param[in,out]  constant       - constant parts of a linear system.
/// @param[in]           len       - the width of the matrix a, i.e., the number of column vectors.
/// @return   1(true) if success. 0(false) if the matrix is singular.
///
unsigned gf256mat_gaussian_elim(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

/// @brief Back substitution of the constant terms with a row echelon form of a matrix, in GF(16)
///
/// @param[in,out]  constant       - constant parts of a linear system.
/// @param[in]     sq_row_mat_a    - row echelon form of a linear system.
/// @param[in]           len       - the height of the matrix a, i.e., the number of row vectors.
///
void gf256mat_back_substitute( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);


unsigned gf256mat_gaussian_elim_unde(unsigned char *sqmat_a, unsigned A_rows, unsigned A_cols, unsigned char *constant);

void gf256mat_back_substitute_unde(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col);

#define _TUOVSIGN_DEBUG_
#ifdef _TUOVSIGN_DEBUG_
void gfmat_transpose_vxo_change(unsigned char *M);
void gfmat_transpose_oxo_change(unsigned char *M);
#endif


#ifdef  __cplusplus
}
#endif

#endif  // _BLAS_MATRIX_H_

