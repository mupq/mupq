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
#include "stdio.h"
#include "gf16.h"

static inline uint8_t __gf_mul(uint8_t a, uint8_t b){
#ifdef _USE_GF16
    return gf16_mul(a, b);
#else
    return gf256_mul(a, b);
#endif
}

static inline uint8_t __gf_inv(uint8_t a){
#ifdef _USE_GF16
    return gf16_inv(a);
#else
    return gf256_inv(a);
#endif
}

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

unsigned gfmat_gaussian_elim_unde_change(unsigned char *A, unsigned A_rows, unsigned A_cols, unsigned char *b){
    uint8_t **mat = calloc(A_rows, sizeof(uint8_t *));
    uint8_t *__ptr = malloc(A_rows * A_cols + A_cols);
    for (long i = 0; i < A_rows; i++){
        mat[i] = &__ptr[i * (A_cols + 1)];
    }

    long v = A_cols;
    long o = A_rows;

    long ind = 0;
    for (long i = 0; i < o; i++){
        for (long j = 0; j < v; j++){
            mat[i][j] = __gfv_get_ele(A, ind);
            ind++;
        }
        mat[i][v] = __gfv_get_ele(b, i);
    }

    // preparation done

    do {
        long ind = 0;
        long jnd = 0;
        while (ind < o) {
            for (long i = ind; i < o; i++){
                if (mat[i][jnd]) {
                    if (ind != i) {
                        uint8_t *tmp = mat[i];
                        mat[i] = mat[ind];
                        mat[ind] = tmp;
                    }
                    break;
                }
            }
            if (!mat[ind][jnd]) {
                jnd++;
                continue;
            }
            if (jnd > v) break;
            uint8_t inv = __gf_inv(mat[ind][jnd]);
            for (long j = jnd; j < v+1; j++) mat[ind][j] = __gf_mul(mat[ind][j], inv);
            for (long i = ind+1; i < o; i++){
                uint8_t r = mat[i][jnd];
                for (long j = jnd; j < v+1; j++) mat[i][j] = mat[i][j] ^ __gf_mul(r, mat[ind][j]);
            }
            ind++;
        }
    } while (0);

    unsigned ret = 0;
    for (long i = 0; i < v; i++){
        if (mat[o-1][i]) ret = 1;
    }
    /*printf("row elim form of MT = \n");
    for (long i = 0; i < o; i++){
        printf("[");
        for (long j = 0; j < v; j++){
            printf("%hx ", mat[i][j]);
        }
        printf("%hx]\n", mat[i][v]);
    }*/


    // computation done

    ind = 0;
    for (long i = 0; i < o; i++){
        for (long j = 0; j < v; j++){
            __gfv_set_ele(A, ind, mat[i][j]);
            ind++;
        }
        __gfv_set_ele(b, i, mat[i][v]);
    }

    free(mat);
    free(__ptr);
    return ret;
}

void gfmat_back_substitute_unde_change(unsigned char *sol, unsigned char *A, unsigned A_rows, unsigned A_cols, unsigned char *b){
    uint8_t **mat = calloc(A_rows, sizeof(uint8_t *));
    uint8_t *__ptr = malloc(A_rows * A_cols + A_cols);
    for (long i = 0; i < A_rows; i++){
        mat[i] = &__ptr[i * (A_cols + 1)];
    }

    long v = A_cols;
    long o = A_rows;

    long ind = 0;
    for (long i = 0; i < o; i++){
        for (long j = 0; j < v; j++){
            mat[i][j] = __gfv_get_ele(A, ind);
            ind++;
        }
        mat[i][v] = __gfv_get_ele(b, i);
    }

    // preparation done

    uint8_t *_sol = malloc(v);
    do {
        long ind = o - 1;
        long jnd = v - 1;
        while (jnd >= 0) {
            long knd = 0;
            while (!mat[ind][knd]) knd++;
            if (knd > jnd) {
                // no sol exist!
                fprintf(stderr, "[Error] solve: something must be wrong, aborted.\n");
                exit(1);
            }
            while (jnd > knd) {
                _sol[jnd] = Uniform_gf();
                jnd--;
            }
            if (mat[ind][knd] != 1) {
                fprintf(stderr, "[Error] solve: something must be wrong, aborted.\n");
                exit(1);
            }
            _sol[jnd] = mat[ind][v];
            for (long j = v - 1; j > jnd; j--) _sol[jnd] ^= __gf_mul(_sol[j], mat[ind][j]);
            ind--;
            jnd--;
            if (ind < 0){
                while (jnd >= 0) {
                    _sol[jnd] = Uniform_gf();
                    jnd--;
                }
            }
        }
    } while(0);

    for (long i = 0; i < v; i++) __gfv_set_ele(sol, i, _sol[i]);


    // computation done

    ind = 0;
    for (long i = 0; i < o; i++){
        for (long j = 0; j < v; j++){
            __gfv_set_ele(A, ind, mat[i][j]);
            ind++;
        }
        __gfv_set_ele(b, i, mat[i][v]);
    }
    free(_sol);
    free(mat);
    free(__ptr);
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
