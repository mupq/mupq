//// @file blas_matrix.c
/// @brief The standard implementations for blas_matrix.h
///

//#include "blas_comm.h"
#include "blas_matrix.h"
//#include "blas.h"
#include "params.h"  // for macro _USE_GF16


#include "config.h"
// choosing the implementations depends on the macros _BLAS_AVX2_ and _BLAS_SSE_

//
// These functions of matrix operations are considered heavy funcitons.
// The cost of an extra funciton call is relatively smaller than computations.
//


#if defined( _BLAS_AVX2_ )

#include "blas_matrix_avx2.h"
#include "blas_matrix_ref.h" // need to delete it

#define gf16mat_prod_impl             gf16mat_prod_avx2
#define gf256mat_prod_impl            gf256mat_prod_avx2

#define gf16mat_prod_multab_impl      gf16mat_prod_multab_avx2
#define gf256mat_prod_multab_impl     gf256mat_prod_multab_avx2

#define gf16mat_linearcomb_half_impl  gf16mat_linearcomb_half_ref
#define gf256mat_linearcomb_half_impl gf256mat_linearcomb_half_ref

#define gf256mat_gaussian_elim_impl   gf256mat_gaussian_elim_avx2
#define gf256mat_back_substitute_impl gf256mat_back_substitute_avx2
#define gf16mat_gaussian_elim_impl   gf16mat_gaussian_elim_avx2
#define gf16mat_back_substitute_impl gf16mat_back_substitute_avx2

#define gf16mat_gaussian_elim_unde_impl gf16mat_gaussian_elim_unde_ref
#define gf256mat_gaussian_elim_unde_impl gf256mat_gaussian_elim_unde_ref

#define gf16mat_back_substitute_unde_impl gf16mat_back_substitute_unde_ref
#define gf256mat_back_substitute_unde_impl gf256mat_back_substitute_unde_ref

#else

#include "blas_matrix_ref.h"

#define gf16mat_prod_impl             gf16mat_prod_ref
#define gf256mat_prod_impl            gf256mat_prod_ref

#define gf16mat_linearcomb_half_impl  gf16mat_linearcomb_half_ref
#define gf256mat_linearcomb_half_impl gf256mat_linearcomb_half_ref

#define gf256mat_gaussian_elim_impl   gf256mat_gaussian_elim_ref
#define gf256mat_back_substitute_impl gf256mat_back_substitute_ref
#define gf16mat_gaussian_elim_impl   gf16mat_gaussian_elim_ref
#define gf16mat_back_substitute_impl gf16mat_back_substitute_ref

#define gf16mat_gaussian_elim_unde_impl gf16mat_gaussian_elim_unde_ref
#define gf256mat_gaussian_elim_unde_impl gf256mat_gaussian_elim_unde_ref

#define gf16mat_back_substitute_unde_impl gf16mat_back_substitute_unde_ref
#define gf256mat_back_substitute_unde_impl gf256mat_back_substitute_unde_ref
#endif


#ifdef _USE_GF16

void gf16mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    gf16mat_prod_impl( c, matA, n_A_vec_byte, n_A_width, b);
}

#if defined(_MUL_WITH_MULTAB_)
void gf16mat_prod_multab(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    gf16mat_prod_multab_impl( c, matA, n_A_vec_byte, n_A_width, b);
}
#endif

void gf16mat_linearcomb_half(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch) {
    gf16mat_linearcomb_half_impl(M, F, F_rows, F_cols, lambda, batch);
}

unsigned gf16mat_gaussian_elim(uint8_t *sqmat_a, uint8_t *constant, unsigned len) {
    return gf16mat_gaussian_elim_impl(sqmat_a, constant, len );
}

void gf16mat_back_substitute( uint8_t *constant, const uint8_t *sqmat_a, unsigned len) {
    gf16mat_back_substitute_impl( constant, sqmat_a, len );
}

unsigned gf16mat_gaussian_elim_unde(unsigned char *sqmat_a, unsigned A_rows, unsigned A_cols, unsigned char *constant){
    return gf16mat_gaussian_elim_unde_impl(sqmat_a, A_rows, A_cols, constant);
}

void gf16mat_back_substitute_unde(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col) {
    gf16mat_back_substitute_unde_impl(v, M, M_row, M_col);
}

#else  //  _USE_GF256

void gf256mat_prod(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    gf256mat_prod_impl( c, matA, n_A_vec_byte, n_A_width, b);
}

#if defined(_MUL_WITH_MULTAB_)
void gf256mat_prod_multab(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    gf256mat_prod_multab_impl( c, matA, n_A_vec_byte, n_A_width, b);
}
#endif

void gf256mat_linearcomb_half(unsigned char *M, const unsigned char *F, unsigned F_rows, unsigned F_cols, const unsigned char *lambda, unsigned batch) {
    gf256mat_linearcomb_half_impl(M, F, F_rows, F_cols, lambda, batch);
}

unsigned gf256mat_gaussian_elim(uint8_t *sqmat_a, uint8_t *constant, unsigned len) {
    return gf256mat_gaussian_elim_impl(sqmat_a, constant, len );
}

void gf256mat_back_substitute( uint8_t *constant, const uint8_t *sqmat_a, unsigned len) {
    gf256mat_back_substitute_impl( constant, sqmat_a, len );
}

unsigned gf256mat_gaussian_elim_unde(unsigned char *sqmat_a, unsigned A_rows, unsigned A_cols, unsigned char *constant){
    return gf256mat_gaussian_elim_unde_impl(sqmat_a, A_rows, A_cols, constant);
}

void gf256mat_back_substitute_unde(uint8_t * v, const uint8_t *M, unsigned M_row,  unsigned M_col) {
    gf256mat_back_substitute_unde_impl(v, M, M_row, M_col);
}

#endif  // #ifdef _USE_GF16


#ifdef _TUOVSIGN_DEBUG_

#include "stdlib.h"
#include "gf16.h"

static inline uint8_t __gfv_get_ele(const uint8_t *a, unsigned i){
#ifdef _USE_GF16
    uint8_t r = a[i >> 1];
    return (i & 1) ? (r >> 4) : (r & 0xf);
#else
    return a[i];
#endif
}

static inline uint8_t __gfv_set_ele(uint8_t *a, unsigned i, uint8_t v) {
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


// use rand() so every thing is deterministic
uint8_t Uniform_gf(){
#ifdef _USE_GF16
    return rand()&0xf;
#else
    return rand()&0xff;
#endif
}

void gfmat_transpose_vxo_change(unsigned char *M){
    uint8_t _M[_V * _O];
    for (long i = 0; i < _V; i++){
        for (long j = 0; j < _O; j++){
            _M[i * _O + j] = __gfv_get_ele(M + i * _O_BYTE, j);
        }
    }
    for (long i = 0; i < _V; i++){
        for (long j = 0; j < _O; j++){
             __gfv_set_ele(M + j * _V_BYTE, i, _M[i * _O + j]);
        }
    }
}

void gfmat_transpose_oxo_change(unsigned char *M){
    uint8_t _M[_O * _O];
    for (long i = 0; i < _O; i++){
        for (long j = 0; j < _O; j++){
            _M[i * _O + j] = __gfv_get_ele(M + i * _O_BYTE, j);
        }
    }
    for (long i = 0; i < _O; i++){
        for (long j = 0; j < _O; j++){
             __gfv_set_ele(M + j * _O_BYTE, i, _M[i * _O + j]);
        }
    }
}
#endif
