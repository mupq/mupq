/// @file blas_matrix_ref.h
/// @brief linear algebra functions for matrix opt.
///
#ifndef _BLAS_MATRIX_REF_H_
#define _BLAS_MATRIX_REF_H_

#include <stdint.h>



#ifdef  __cplusplus
extern  "C" {
#endif


///////////////// Section: multiplications  ////////////////////////////////

// matrix-vector

void gf16mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);

void gf256mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b);


/// @brief calculate vector a inner-product vector b. GF(256)
///
/// @param[out]  accu_c     - the output.
/// @param[in]   a          - the vector a.
/// @param[in]   b          - the vector b.
/// @param[in]  _num_byte   - number of bytes for the vector a/b.
/// @return  1(true) if a is 0. 0(false) else.
///
void _gf256v_mul_add(unsigned char *accu_c, const unsigned char *a, const unsigned char *b, unsigned _num_byte);

/// @brief calculate vector a inner-product vector b. GF(16)
///
/// @param[out]  accu_c     - the output.
/// @param[in]   a          - the vector a.
/// @param[in]   b          - the vector b.
/// @param[in]  _num_byte   - number of bytes for the vector a/b.
/// @return  1(true) if a is 0. 0(false) else.
///
void _gf16v_mul_add(unsigned char *accu_c, const unsigned char *a, const unsigned char *b, unsigned _num_byte);

/////////////////////////////////////////////////////

void gf16mat_linearcomb_half_ref(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch);

void gf256mat_linearcomb_half_ref(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch);

unsigned gf256mat_gaussian_elim_ref(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

void gf256mat_back_substitute_ref( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);

unsigned gf16mat_gaussian_elim_ref(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len);

void gf16mat_back_substitute_ref( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len);

unsigned gf16mat_gaussian_elim_unde_ref(unsigned char *A, unsigned A_rows, unsigned A_cols, unsigned char *constant);

unsigned gf256mat_gaussian_elim_unde_ref(unsigned char *A, unsigned A_rows, unsigned A_cols, unsigned char *constant);

void gf16mat_back_substitute_unde_ref(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col);

void gf256mat_back_substitute_unde_ref(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col);

#ifdef  __cplusplus
}
#endif

#endif  // _BLAS_MATRIX_REF_H_

