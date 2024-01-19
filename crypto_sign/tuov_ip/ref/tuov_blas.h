
/** \file   tuov_blas.h
 *  \brief  all linear algebra related function we used.
 */
#ifndef _OV_BLAS_H_
#define _OV_BLAS_H_

#include "params.h"
#include "utils_randombytes.h"
#include <stdint.h>


/* -------------------------- gf operations -------------------------- */

#include "gf16.h"

/** \brief gf multiplication
 *  \param[in] a    - gf element
 *  \param[in] b    - gf element
 *  \return a * b
 */
static inline uint8_t gf_mul(uint8_t a, uint8_t b){
#ifdef _USE_GF16
    return gf16_mul(a, b);
#else
    return gf256_mul(a, b);
#endif
}

/** \brief gf inversion
 *  \param a        - gf element
 *  \return the inverse of a, return 0 if a is 0
 */
static inline uint8_t gf_inv(uint8_t a){
#ifdef _USE_GF16
    return gf16_inv(a);
#else
    return gf256_inv(a);
#endif
}

/** \brief gf square
 *  \param a        - gf element
 *  \return the square of a
 */
static inline uint8_t gf_squ(uint8_t a){
#ifdef _USE_GF16
    return gf16_squ(a);
#else
    return gf256_squ(a);
#endif
}




/* -------------------------- gf vector operations -------------------------- */

/** \brief get an element from gf vector.
 *
 *  \param[in]  a         - the input vector a.
 *  \param[in]  i         - the index in the vector a.
 *  \return  the value of the element.
 */
static inline uint8_t gfv_get_ele(const uint8_t *a, unsigned i){
#ifdef _USE_GF16
    uint8_t r = a[i >> 1];
    return (i & 1) ? (r >> 4) : (r & 0xf);
#else
    return a[i];
#endif
}

/** \brief set an element for a gf vector .
 *
 *  \param[in,out]   a   - the vector a.
 *  \param[in]  i        - the index in the vector a.
 *  \param[in]  v        - the value for the i-th element in vector a.
 *  \return  the value of the element.
 */
static inline uint8_t gfv_set_ele(uint8_t *a, unsigned i, uint8_t v) {
#ifdef _USE_GF16
    uint8_t ai = a[i >> 1];
    uint8_t i_1_or_16 = (i & 1) * 15 + 1; // 0 -> 1 , 1 -> 16
    ai &= ~(0xf * i_1_or_16); // 0 -> clear lower nibble, 1 -> clear high nibble.
    // v &= 0xf;
    a[i >> 1] = ai + v * i_1_or_16;
    return v;
#else
    a[i] = v;
    return v;
#endif
}

/** \brief return the dot product of two gf vectors 
 *  \param[in] a            - input vector
 *  \param[in] b            - input vector
 *  \param[in] _num_byte    - number of bytes of each vectors
 *  \return the dot_product
*/
static inline uint8_t gfv_dot(const uint8_t *a, const uint8_t *b, unsigned _num_byte){
    uint8_t ret = 0x0;
#ifdef _USE_GF16
    _num_byte <<= 1;
#endif
    for (unsigned i = 0; i < _num_byte; i++){
        ret ^= gf_mul(gfv_get_ele(a, i), gfv_get_ele(b, i));
    }
    return ret;
}


#include "blas.h"

/** \brief add src to dst
 *  \param[in,out] dst      - a gf vector.
 *  \param[in] src          - a gf vector.
 *  \param[in]  _num_byte   - number of bytes for the vector dst.
 */
static inline void gfv_add(uint8_t *dst, const uint8_t *src, unsigned _num_byte){
    gf256v_add(dst, src, _num_byte);
}

/** \brief set a vector to 0.
 *
 *  \param[in,out]   dst      - a gf vector.
 *  \param[in]  _num_byte   - number of bytes for the vector dst.
 */
static inline void gfv_set_zero(uint8_t *dst, unsigned _num_byte){
    gfv_add(dst, dst, _num_byte);
}

/** \brief multiple a vector by a.
 *
 *  \param[in,out]   dst      - a gf vector.
 *  \param[in] a              - a gf element
 *  \param[in]  _num_byte   - number of bytes for the vector dst.
 */
static inline void gfv_mul_scalar(uint8_t *dst, uint8_t a, unsigned _num_byte){
#ifdef _USE_GF16
    gf16v_mul_scalar(dst, a, _num_byte);
#else
    gf256v_mul_scalar(dst, a, _num_byte);
#endif
}

/** \brief multiple vector src by a, and add it to dst
 *
 *  \param[in,out]   dst      - a gf vector.
 *  \param[in] src            - a gf vector.
 *  \param[in] a              - a gf element
 *  \param[in]  _num_byte   - number of bytes for the vector dst.
 */
static inline void gfv_madd(uint8_t *dst, const uint8_t *src, uint8_t a, unsigned _num_byte){
#ifdef _USE_GF16
    gf16v_madd(dst, src, a, _num_byte);
#else
    gf256v_madd(dst, src, a, _num_byte);
#endif
}

#if defined( _BLAS_AVX2_ ) || defined( _BLAS_SSE_ )

/** \brief multiple vector src by an gf element, and add it to dst
 *
 *  \param[in,out]   dst      - a gf vector.
 *  \param[in] src            - a gf vector.
 *  \param[in] multab         - multab of the gf element
 *  \param[in]  _num_byte   - number of bytes for the vector dst.
 */
static inline void gfv_madd_multab(uint8_t *dst, const uint8_t *src, const uint8_t *multab, unsigned _num_byte){
#ifdef _USE_GF16
    gf16v_madd_multab(dst, src, multab, _num_byte);
#else
    gf256v_madd_multab(dst, src, multab, _num_byte);
#endif
}

/** \brief generate multabs for each element of vector v, the multab of v[i] is stored in 
 *          multabs + 32 * i, the input multabs should be aligned to 32 bytes.
 *  \param[out] multabs         - the multabs
 *  \param[in] v                - the input vector
 *  \param[in] n_ele            - number of element in v
*/
static inline void gfv_generate_multabs(uint8_t *multabs, const uint8_t *v, unsigned n_ele){
#ifdef _USE_GF16
    gf16v_generate_multabs(multabs, v, n_ele);
#else
    gf256v_generate_multabs_avx2(multabs, v, n_ele);
#endif
}

#endif



/* -------------------------- gf matrix operations -------------------------- */

#include "blas_matrix.h"

/** \brief matrix-vector multiplication:  c = matA * b
 *
 *  \param[out]  c         - the output vector c
 *  \param[in]   matA      - a column-major matrix A.
 *  \param[in]   n_A_vec_byte  - the size of column vectors in bytes.
 *  \param[in]   n_A_width   - the width of matrix A.
 *  \param[in]   b          - the vector b.
 */
static inline void gfmat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b){
#ifdef _USE_GF16
    gf16mat_prod(c, matA, n_A_vec_byte, n_A_width, b);
#else
    gf256mat_prod(c, matA, n_A_vec_byte, n_A_width, b);
#endif
}

/** \brief matrix-vector multiplication:  c = matA * b
 *
 *  \param[out]  c         - the output vector c
 *  \param[in]   matA      - a column-major matrix A.
 *  \param[in]   n_A_vec_byte  - the size of column vectors in bytes.
 *  \param[in]   n_A_width   - the width of matrix A.
 *  \param[in]   multab_b     - the multiplication tables of the vector b.
 */
static inline void gfmat_prod_multab(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *multab_b){
#ifdef _USE_GF16
    gf16mat_prod_multab(c, matA, n_A_vec_byte, n_A_width, multab_b);
#else
    gf256mat_prod_multab(c, matA, n_A_vec_byte, n_A_width, multab_b);
#endif
}

/** \brief Computing the row echelon form of a matrix
 *
 *  \param[in,out]  sq_col_mat_a   - square matrix parts of a linear system. a is a column major matrix.
 *                                  The returned matrix of row echelon form is a row major matrix.
 *  \param[in,out]  constant       - constant parts of a linear system.
 *  \param[in]           len       - the width of the matrix a, i.e., the number of column vectors.
 *  \return   1(true) if success. 0(false) if the matrix is singular.
 */
static inline unsigned gfmat_gaussian_elim(uint8_t *sq_col_mat_a, uint8_t *constant, unsigned len){
#ifdef _USE_GF16
    return gf16mat_gaussian_elim(sq_col_mat_a, constant, len);
#else
    return gf256mat_gaussian_elim(sq_col_mat_a, constant, len);
#endif
}

/** \brief Back substitution of the constant terms with a row echelon form of a matrix
 *
 *  \param[in,out]  constant       - constant parts of a linear system.
 *  \param[in]     sq_row_mat_a    - row echelon form of a linear system.
 *  \param[in]           len       - the height of the matrix a, i.e., the number of row vectors.
 */
static inline void gfmat_back_substitute( uint8_t *constant, const uint8_t *sq_row_mat_a, unsigned len){
#ifdef _USE_GF16
    gf16mat_back_substitute(constant, sq_row_mat_a, len);
#else
    gf256mat_back_substitute(constant, sq_row_mat_a, len);
#endif
}

/** \brief Computing the row echelon form of a undetermined linear system Ax = b
 *  \param[in, out] A           - a row-major matrix
 *  \param[in] A_rows           - number of rows of A
 *  \param[in] A_cols           - number of cols of A
 *  \param[in, out] b           - the coefficient part b, a vector with A_rows elements
 *  \return 1 if success. 0 if fail.
 */
static inline unsigned gfmat_gaussian_elim_unde(unsigned char *A, unsigned A_rows, unsigned A_cols, unsigned char *b){

#ifdef _USE_GF16
    return gf16mat_gaussian_elim_unde(A, A_rows, A_cols, b);
#else
    return gf256mat_gaussian_elim_unde(A, A_rows, A_cols, b);
#endif
}

/** \brief Computing a random solution v of linear system Mv = b
 *  \param[in, out] M           - a row-major matrix
 *  \param[in] M_rows           - number of rows of M
 *  \param[in] M_cols           - number of cols of M
 *  \param[in, out] v           - a vector with M_cols elements, the first M_rows element = b, the last elements = random element
 *  \return 1 if success. 0 if fail.
 */
static inline void gfmat_back_substitute_unde(unsigned char *v, unsigned char *M, unsigned M_rows, unsigned M_cols, unsigned char *b){
#ifdef _USE_GF16
    memcpy(v, b, M_rows >> 1);
    // randomization in the post V-M dimension can maximize the desirable range of vinegar variables as much as possible
    randombytes(v + (M_rows >> 1), (M_cols - M_rows) >> 1);
    gf16mat_back_substitute_unde(v, M, M_rows, M_cols);
#else
    memcpy(v, b, M_rows);
    randombytes(v + M_rows, M_cols - M_rows);
    gf256mat_back_substitute_unde(v, M, M_rows, M_cols);
#endif

}

/** \brief F is an batched-matrix, lambda is a vector with batch/2 elements, 
 *          we want to compute M = lambda_0 * F_0 + ... + lambda_{batch/2-1} * F_{batch/2-1}.
 *  \param[out] M               - a col-major F_rows * F_cols matrix             - 
 *  \param[in] F                - a batched row-major F_rows * F_cols matrix F
 *  \param[in] F_rows           - number of rows of F
 *  \param[in] F_cols           - number of cols of F
 *  \param[in] lambda           - the coefficient vector, with batch_size/2 elements
 *  \param[in] batch            - batch size of F
 */ 
static inline void gfmat_linearcomb_half(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch){
#ifdef _USE_GF16
    gf16mat_linearcomb_half(M, F, F_rows, F_cols, lambda, batch);
#else
    gf256mat_linearcomb_half(M, F, F_rows, F_cols, lambda, batch);
#endif
}

/** \brief transpose a _V * _O matrix
 *  \param[in, out] M       - a _V * _O matrix in row-major form, output the transpose of M (still in row major form)
*/
static inline void gfmat_transpose_vxo(unsigned char *M){
    gfmat_transpose_vxo_change(M);
}

/** \brief transpose a _O * _O matrix
 *  \param[in, out] M       - a _O * _O matrix in row-major form, output the transpose of M (still in row major form)
*/
static inline void gfmat_transpose_oxo(unsigned char *M){
    gfmat_transpose_oxo_change(M);
}


/* -------------------------- gf batched matrix operations -------------------------- */

#include "parallel_matrix_op.h"


/** \brief  bC += btriA * B
 * 
 *  \param[out]  bC         - the batched matrix C.
 *  \param[in]   btriA      - a batched UT matrix A.
 *  \param[in]   B          - a column-major matrix B.
 *  \param[in]   Bheight          - the height of B.
 *  \param[in]   size_Bcolvec     - the size of the column vector in B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_trimat_madd( unsigned char *bC, const unsigned char *btriA,
                             const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_trimat_madd_gf16(bC, btriA, B, Bheight, size_Bcolvec, Bwidth, size_batch);
#else
    batch_trimat_madd_gf256(bC, btriA, B, Bheight, size_Bcolvec, Bwidth, size_batch);
#endif
}



/** \brief  bC = A^Tr * bB  and then upper(bC)
 * 
 *  \param[out]  bC           - the batched matrix C.
 *  \param[in]   A_to_tr      - a column-major matrix A. The operand for multiplication is A^Tr.
 *  \param[in]   Aheight      - the height of A.
 *  \param[in]   size_Acolvec    - the size of a column vector in A.
 *  \param[in]   Awidth           - the width of A.
 *  \param[in]   bB          - a batched matrix B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */ 
static inline void batch_upper_matTr_x_mat( unsigned char *bC,
                                   const unsigned char *A_to_tr, unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
                                   const unsigned char *bB, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_upper_matTr_x_mat_gf16(bC, A_to_tr, Aheight, size_Acolvec, Awidth, bB, Bwidth, size_batch);
#else
    batch_upper_matTr_x_mat_gf256(bC, A_to_tr, Aheight, size_Acolvec, Awidth, bB, Bwidth, size_batch);
#endif
}



/** \brief  bC +=  (btriA + btriA^Tr) *B
 * 
 *  \param[out]  bC         - the batched matrix C.
 *  \param[in]   btriA      - a batched UT matrix A. The operand for multiplication is (btriA + btriA^Tr).
 *  \param[in]   B          - a column-major matrix B.
 *  \param[in]   Bheight          - the height of B.
 *  \param[in]   size_Bcolvec     - the size of the column vector in B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_2trimat_madd( unsigned char *bC, const unsigned char *btriA,
                              const unsigned char *B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_2trimat_madd_gf16(bC, btriA, B, Bheight, size_Bcolvec, Bwidth, size_batch);
#else
    batch_2trimat_madd_gf256(bC, btriA, B, Bheight, size_Bcolvec, Bwidth, size_batch);
#endif
}



/** \brief  y =  x^Tr * trimat * x 
 * 
 *  \param[out]  y          - the returned batched element y.
 *  \param[in]   trimat     - a batched matrix.
 *  \param[in]   x          - an input vector x.
 *  \param[in]   dim        - the dimension of matrix trimat (and x).
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_quad_trimat_eval( unsigned char *y, const unsigned char *trimat, const unsigned char *x, unsigned dim, unsigned size_batch ){
#ifdef _USE_GF16
    batch_quad_trimat_eval_gf16(y, trimat, x, dim, size_batch);
#else
    batch_quad_trimat_eval_gf256(y, trimat, x, dim, size_batch);
#endif
}
           
#if defined( _BLAS_AVX2_ ) || defined( _BLAS_SSE_)

/** \brief  bC +=  (btriA + btriA^Tr) *B
 * 
 *  \param[out]  bC         - the batched matrix C.
 *  \param[in]   btriA      - a batched UT matrix A. The operand for multiplication is (btriA + btriA^Tr).
 *  \param[in]   multab_B   - a column-major matrix B.
 *  \param[in]   Bheight          - the height of B.
 *  \param[in]   size_Bcolvec     - the size of the column vector in B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_2trimat_madd_multab( unsigned char *bC, const unsigned char *btriA,
                                     const unsigned char *multab_B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_2trimat_madd_multab_gf16(bC, btriA, multab_B, Bheight, size_Bcolvec, Bwidth, size_batch);
#else
    batch_2trimat_madd_multab_gf256(bC, btriA, multab_B, Bheight, size_Bcolvec, Bwidth, size_batch);
#endif
}



/** \brief  bC = A^Tr * bB  and then upper(bC)
 * 
 *  \param[out]  bC           - the batched matrix C.
 *  \param[in]   multab_A_to_tr      - a column-major matrix A. The operand for multiplication is A^Tr.
 *  \param[in]   Aheight      - the height of A.
 *  \param[in]   size_Acolvec    - the size of a column vector in A.
 *  \param[in]   Awidth           - the width of A.
 *  \param[in]   bB          - a batched matrix B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_upper_matTr_x_mat_multab( unsigned char *bC,
        const unsigned char *multab_A_to_tr, unsigned Aheight, unsigned size_Acolvec, unsigned Awidth,
        const unsigned char *bB, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_upper_matTr_x_mat_multab_gf16(bC, multab_A_to_tr, Aheight, size_Acolvec, Awidth, bB, Bwidth, size_batch);
#else
    batch_upper_matTr_x_mat_multab_gf256(bC, multab_A_to_tr, Aheight, size_Acolvec, Awidth, bB, Bwidth, size_batch);
#endif
}



/** \brief  bC += btriA * B
 * 
 *  \param[out]  bC         - the batched matrix C.
 *  \param[in]   btriA      - a batched UT matrix A.
 *  \param[in]   multab_B   - a column-major matrix B.
 *  \param[in]   Bheight          - the height of B.
 *  \param[in]   size_Bcolvec     - the size of the column vector in B.
 *  \param[in]   Bwidth           - the width of B.
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_trimat_madd_multab( unsigned char *bC, const unsigned char *btriA,
                                    const unsigned char *multab_B, unsigned Bheight, unsigned size_Bcolvec, unsigned Bwidth, unsigned size_batch ){
#ifdef _USE_GF16
    batch_trimat_madd_multab_gf16(bC, btriA, multab_B, Bheight, size_Bcolvec, Bwidth, size_batch);
#else
    batch_trimat_madd_multab_gf256(bC, btriA, multab_B, Bheight, size_Bcolvec, Bwidth, size_batch);
#endif
}



/** \brief  y =  x^Tr * trimat * x
 * 
 *  \param[out]  y          - the returned batched element y.
 *  \param[in]   trimat     - a batched matrix.
 *  \param[in]   multabs_x  - multabs of input vector x.
 *  \param[in]   dim        - the dimension of matrix trimat (and x).
 *  \param[in]   size_batch - number of the batched elements in the corresponding position of the matrix.
 */
static inline void batch_quad_trimat_eval_multab( unsigned char *y, const unsigned char *trimat, const unsigned char *multabs_x, unsigned dim, unsigned size_batch ){
#ifdef _USE_GF16
    batch_quad_trimat_eval_multab_gf16(y, trimat, multabs_x, dim, size_batch);
#else
    batch_quad_trimat_eval_multab_gf256(y, trimat, multabs_x, dim, size_batch);
#endif
}

#endif

#endif
