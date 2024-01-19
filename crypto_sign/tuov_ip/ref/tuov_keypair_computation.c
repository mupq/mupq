/** \file   tuov_keypair_computation.h
 *  \brief  Functions for calculating pk/sk while generating keys.
 *          used in tuov_keypair.c
 */
#include "params.h"
#include "config.h"
#include "tuov_keypair_computation.h"

#include "tuov_blas.h"

#include <stdlib.h>
#include <string.h>
#include "stdint.h"


void multiple_S(unsigned char *M, const unsigned char *S, unsigned D){
    uint8_t tmp[_PUB_M1_BYTE];

    #ifdef _MUL_WITH_MULTAB_
    // multab for all field element
    uint8_t multab_all[_GFSIZE * 32] __attribute__((aligned(32)));
    // multabs for one col of M
    uint8_t multab_Mcol[_PUB_M1 * 32] __attribute__((aligned(32)));
    #ifdef _USE_GF16
    #define _GFBYTE (_GFSIZE >> 1)
    #else
    #define _GFBYTE (_GFSIZE)
    #endif
    uint8_t ele_all[_GFBYTE] = {0};
    for (unsigned i = 0; i < _GFSIZE; i++){
        gfv_set_ele(ele_all, i, i&0xff);
    }
    gfv_generate_multabs(multab_all, ele_all, _GFSIZE);
    for (unsigned i = 0; i < D; i++){
        for (unsigned j = 0; j < _PUB_M1; j++){
            memcpy(&multab_Mcol[j * 32], &multab_all[((unsigned int)gfv_get_ele(M, j)) * 32], 32);
        }
        gfmat_prod_multab(tmp, S, _PUB_M1_BYTE, _PUB_M1, multab_Mcol);
        memcpy(M, tmp, _PUB_M1_BYTE);
        M += _PUB_M1_BYTE;
    }
    #undef _GFBYTE

    #else
    for (unsigned i = 0; i < D; i++){
        gfmat_prod(tmp, S, _PUB_M1_BYTE, _PUB_M1, M);
        memcpy(M, tmp, _PUB_M1_BYTE);
        M += _PUB_M1_BYTE;
    }
    #endif
}

/** \brief combine two m1 batched matrix to a m-batched matrix
 *  \param[in,out] M        - M is an array of size M_ele * _PUB_M_BYTE
 *                          the first M1_ele * _PUB_M1_BYTE bytes of M is an m1-batched matrix,
 *                          the last M1_ele * _PUB_M1_BYTE bytes of M is another m1-batched matrix,
 *                          the output is a m-batched matrix
 *  \param[in] M_ele        - number of elements in M
 */
#define _combine_m1batched_matrix(M, M_ele) do {                            \
    unsigned char tmp[_PUB_M1_BYTE * M_ele] __attribute__((aligned(32)));   \
    memcpy(tmp, M, _PUB_M1_BYTE * M_ele);                                   \
    unsigned char *M1 = tmp;                                                \
    unsigned char *M2 = M + _PUB_M1_BYTE * M_ele;                           \
    unsigned char *dst = M;                                                 \
    for (unsigned _i = 0; _i < M_ele - 1; _i++){                            \
        memcpy(dst, M1, _PUB_M1_BYTE);                                      \
        dst += _PUB_M1_BYTE;                                                \
        memcpy(dst, M2, _PUB_M1_BYTE);                                      \
        dst += _PUB_M1_BYTE;                                                \
        M1 += _PUB_M1_BYTE;                                                 \
        M2 += _PUB_M1_BYTE;                                                 \
    }                                                                       \
    memcpy(dst, M1, _PUB_M1_BYTE);                                          \
} while (0)

void combine_F(unsigned char *F){
    _combine_m1batched_matrix(F + P11_BIAS, N_TRIANGLE_TERMS(_V));
    uint8_t *new_F23 = F + P21_BIAS;
    uint8_t F2[_PUB_M_BYTE * _V];
    uint8_t F3[_PUB_M_BYTE * _V * (_O - 1)];
    uint8_t *F21 = F2;
    uint8_t *F22 = F2 + _PUB_M1_BYTE * _V;
    uint8_t *F31 = F3;
    uint8_t *F32 = F3 + _PUB_M1_BYTE * _V * (_O - 1);
    memcpy(F21, F + P21_BIAS, _PK_P2_BYTE >> 1);
    memcpy(F22, F + P22_BIAS, _PK_P2_BYTE >> 1);
    memcpy(F31, F + P31_BIAS, _PK_P3_BYTE >> 1);
    memcpy(F32, F + P32_BIAS, _PK_P3_BYTE >> 1);
    for (long b = 0; b < _PUB_M1; b++){
        // F21
        for (long i = 0; i < _V; i++){
            gfv_set_ele(new_F23 + i * _O_BYTE, _O-1, gfv_get_ele(F21 + i * _PUB_M1_BYTE, b));
        }
        // F31
        for (long i = 0; i < _V; i++){
            for (long j = 0; j < _O - 1; j++){
                gfv_set_ele(new_F23 + i * _O_BYTE, j, gfv_get_ele(F31 + (i * (_O - 1) + j) * _PUB_M1_BYTE, b));
            }
        }
        new_F23 += _O * _V_BYTE;
    }
    for (long b = 0; b < _PUB_M1; b++){
        // F22
        for (long i = 0; i < _V; i++){
            gfv_set_ele(new_F23 + i * _O_BYTE, _O-1, gfv_get_ele(F22 + i * _PUB_M1_BYTE, b));
        }
        // F32
        for (long i = 0; i < _V; i++){
            for (long j = 0; j < _O - 1; j++){
                gfv_set_ele(new_F23 + i * _O_BYTE, j, gfv_get_ele(F32 + (i * (_O - 1) + j) * _PUB_M1_BYTE, b));
            }
        }
        new_F23 += _O * _V_BYTE;
    }

    //_combine_m1batched_matrix(F + P21_BIAS, (_V));
    //_combine_m1batched_matrix(F + P31_BIAS, (_V * (_PUB_M - 1)));
}

void combine_P(unsigned char *P){
    _combine_m1batched_matrix(P + P11_BIAS, N_TRIANGLE_TERMS(_V));
    _combine_m1batched_matrix(P + P21_BIAS, (_V));
    _combine_m1batched_matrix(P + P31_BIAS, (_V * (_PUB_M - 1)));
    _combine_m1batched_matrix(P + P51_BIAS, 1);
    _combine_m1batched_matrix(P + P61_BIAS, (_PUB_M - 1));
    _combine_m1batched_matrix(P + P91_BIAS, (N_TRIANGLE_TERMS(_PUB_M - 1)));

    // convert_123569_to_123

    // P1 --> new P1, already ok
    // P2 & P3 --> new P2
    unsigned char tmp[_PK_P2_BYTE + _PK_P3_BYTE] __attribute__((aligned(32)));
    unsigned char *P2 = tmp;
    unsigned char *P3 = tmp + _PK_P2_BYTE;
    unsigned char *new_P2 = P + _PK_P1_BYTE;
    memcpy(tmp, new_P2, _PK_P2_BYTE + _PK_P3_BYTE);
    for (unsigned i = 0; i < _V; i++){
        // the i-th row of new_P2
        memcpy(new_P2, P2, _PUB_M_BYTE);
        memcpy(new_P2 + _PUB_M_BYTE, P3, _PUB_M_BYTE * (_PUB_M - 1));
        P2 += _PUB_M_BYTE;
        P3 += _PUB_M_BYTE * (_PUB_M - 1);
        new_P2 += _PUB_M_BYTE * _PUB_M;
    }
    // P5, P6, P9 --> new P3, already ok
}

#ifdef _MUL_WITH_MULTAB_
/** \brief dist += Q2 * T3
 *
 *  \param[in,out] dst      - output
 *  \param[in] Q2           - an m1-batched row-major (n - m) * 1 matrix
 *  \param[in] multab_T3    - T3 is a col-major 1 * (m - 1) matrix, multab_T3 is the multab of T3
*/
void add_Q2_mul_T3_multab(unsigned char *dst, const unsigned char *Q2, const unsigned char *multab_T3){
    for (unsigned i = 0; i < _V; i++){
        for (unsigned j = 0; j < _PUB_M - 1; j++){
            // add the (i, j)-th pos of Q2 * T3 to (i, j)-th pos of dst
            gfv_madd_multab(dst, Q2, &multab_T3[j * 32], _PUB_M1_BYTE);
            dst += _PUB_M1_BYTE;
        }
        Q2 += _PUB_M1_BYTE;
    }
}
#endif

/** \brief dist += Q2 * T3
 *
 *  \param[in,out] dst      - output
 *  \param[in] Q2           - an m1-batched row-major (n - m) * 1 matrix
 *  \param[in] T3           - T3 is a col-major 1 * (m - 1) matrix
*/
void add_Q2_mul_T3(unsigned char *dst, const unsigned char *Q2, const unsigned char *T3){
    for (unsigned i = 0; i < _V; i++){
        for (unsigned j = 0; j < _PUB_M - 1; j++){
            // add the (i, j)-th pos of Q2 * T3 to (i, j)-th pos of dst
            gfv_madd(dst, Q2, gfv_get_ele(T3, j), _PUB_M1_BYTE);
            dst += _PUB_M1_BYTE;
        }
        Q2 += _PUB_M1_BYTE;
    }
}


void calculate_Q51(unsigned char *Q51, const unsigned char *T1, const unsigned char *Q11, const unsigned char *Q21){
    uint8_t tmp[_PUB_M1_BYTE];

    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T1[_V * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T1, T1, _V);
    // Q51 = T1_tr * Q11 * T1
    batch_quad_trimat_eval_multab(Q51, Q11, multab_T1, _V, _PUB_M1_BYTE);
    // tmp = T1_tr * Q21
    gfmat_prod_multab(tmp, Q21, _PUB_M1_BYTE, _V, multab_T1);
    #else
    batch_quad_trimat_eval(Q51, Q11, T1, _V, _PUB_M1_BYTE);
    gfmat_prod(tmp, Q21, _PUB_M1_BYTE, _V, T1);
    #endif
    // Q51 += tmp
    gfv_add(Q51, tmp, _PUB_M1_BYTE);
    return;
}


void calculate_Q61(unsigned char *Q61, const unsigned char *T1, const unsigned char *T3, const unsigned char *T4,
                    const unsigned char *Q11, const unsigned char *Q21, const unsigned char *Q31, const unsigned char *Q51){
    ((void) Q51);       // Q51 + Q51_tr = 0

    uint8_t tmp1[_PUB_M1_BYTE * _V];
    uint8_t tmp2[_PUB_M1_BYTE];

    // tmp1 = Q21
    memcpy(tmp1, Q21, _PUB_M1_BYTE * _V);

    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T1[_V * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T1, T1, _V);
    // tmp1 += T1_tr * (Q11 + Q11_tr)
    batch_2trimat_madd_multab(tmp1, Q11, multab_T1, _V, _V_BYTE, 1, _PUB_M1_BYTE);
    // Q61 = T1_tr * Q31
    gfmat_prod_multab(Q61, Q31, _PUB_M1_BYTE * (_PUB_M - 1), _V, multab_T1);
    // tmp2 = T1_tr * Q21
    gfmat_prod_multab(tmp2, Q21, _PUB_M1_BYTE, _V, multab_T1);
    #else
    batch_2trimat_madd(tmp1, Q11, T1, _V, _V_BYTE, 1, _PUB_M1_BYTE);
    gfmat_prod(Q61, Q31, _PUB_M1_BYTE * (_PUB_M - 1), _V, T1);
    gfmat_prod(tmp2, Q21, _PUB_M1_BYTE, _V, T1);
    #endif

    #ifdef _MUL_WITH_MULTAB_
    // Q61 += tmp2 * T3
    uint8_t multab_T3[(_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T3, T3, (_PUB_M - 1));
    uint8_t *_Q61  = Q61;
    uint8_t *_tab_T3 = multab_T3;
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        gfv_madd_multab(_Q61, tmp2, _tab_T3, _PUB_M1_BYTE);
        _tab_T3 += 32;
        _Q61 += _PUB_M1_BYTE;
    }

    // Q61 += tmp1 * T4
    uint8_t multab_T4[_V * (_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T4, T4, _V * (_PUB_M - 1));
    _Q61 = Q61;
    uint8_t *_tab_T4 = multab_T4;
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        for (unsigned j = 0; j < _V; j++){
            gfv_madd_multab(_Q61, tmp1 + _PUB_M1_BYTE * j, _tab_T4, _PUB_M1_BYTE);
            _tab_T4 += 32;
        }
        _Q61 += _PUB_M1_BYTE;
    }
    #else
    // Q61 += tmp2 * T3
    uint8_t *_Q61  = Q61;
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        gfv_madd(_Q61, tmp2, gfv_get_ele(T3, i), _PUB_M1_BYTE);
        _Q61 += _PUB_M1_BYTE;
    }

    // Q61 += tmp1 * T4
    _Q61 = Q61;
    const uint8_t *_T4 = T4;
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        for (unsigned j = 0; j < _V; j++){
            gfv_madd(_Q61, tmp1 + _PUB_M1_BYTE * j, gfv_get_ele(_T4, j), _PUB_M1_BYTE);
        }
        _T4 += _V_BYTE;
        _Q61 += _PUB_M1_BYTE;
    }
    #endif

    return;
}


void calculate_Q91(unsigned char *Q91, const unsigned char *T3, const unsigned char *T4, const unsigned char *Q11,
                    const unsigned char *Q21, const unsigned char *Q31, const unsigned char *Q51, const unsigned char *Q61){

    uint8_t tmp1[_PK_P3_BYTE >> 1];
    uint8_t tmp2[(_PUB_M - 1) * _PUB_M1_BYTE];
    uint8_t tmp3[_PK_P9_BYTE >> 1];

    // tmp1 = Q31, tmp2 = Q61
    memcpy(tmp1, Q31, _PK_P3_BYTE >> 1);
    memcpy(tmp2, Q61, _PK_P6_BYTE >> 1);

    #ifdef _MUL_WITH_MULTAB_
    // tmp1 += Q21 * T3
    uint8_t multab_T3[(_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T3, T3, (_PUB_M - 1));
    add_Q2_mul_T3_multab(tmp1, Q21, multab_T3);

    // tmp2 += Q51 * T3
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        gfv_madd_multab(&tmp2[i * _PUB_M1_BYTE], Q51, &multab_T3[i * 32], _PUB_M1_BYTE);
    }
    #else
    add_Q2_mul_T3(tmp1, Q21, T3);
    for (unsigned i = 0; i < _PUB_M - 1; i++){
        gfv_madd(&tmp2[i * _PUB_M1_BYTE], Q51, gfv_get_ele(T3, i), _PUB_M1_BYTE);
    }
    #endif


    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T4[_V * (_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T4, T4, _V * (_PUB_M - 1));
    // tmp1 += Q1 * T4
    batch_trimat_madd_multab(tmp1, Q11, multab_T4, _V, _V_BYTE, _PUB_M - 1, _PUB_M1_BYTE);
    // Q91 = Upper(T4_tr * tmp1)
    batch_upper_matTr_x_mat_multab(Q91, multab_T4, _V, _V_BYTE, _PUB_M - 1, tmp1, _PUB_M - 1, _PUB_M1_BYTE);
    // tmp3 = Upper(T3_tr * tmp2)
    // warning: I assumed size_Acolvec is not used in the gf16 versions of batch_upper_matTr_x_mat_multab
    batch_upper_matTr_x_mat_multab(tmp3, multab_T3, 1, 1, _PUB_M - 1, tmp2, _PUB_M - 1, _PUB_M1_BYTE);
    #else
    batch_trimat_madd(tmp1, Q11, T4, _V, _V_BYTE, _PUB_M - 1, _PUB_M1_BYTE);
    batch_upper_matTr_x_mat(Q91, T4, _V, _V_BYTE, _PUB_M - 1, tmp1, _PUB_M - 1, _PUB_M1_BYTE);
    batch_upper_matTr_x_mat(tmp3, T3, 1, 1, _PUB_M - 1, tmp2, _PUB_M - 1, _PUB_M1_BYTE);
    #endif



    gfv_add(Q91, tmp3, _PK_P9_BYTE >> 1);
    // Q91 += T3_tr * tmp2


}


void calculate_F21(unsigned char *F21, const unsigned char *T1, const unsigned char *Q11, const unsigned char *Q21){
    // F21 = Q21
    memcpy(F21, Q21, _PK_P2_BYTE >> 1);

    // F21 += (Q11 + Q11_tr) * T1
    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T1[_V * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T1, T1, _V);
    batch_2trimat_madd_multab(F21, Q11, multab_T1, _V, _V_BYTE, 1, _PUB_M1_BYTE);
    #else
    batch_2trimat_madd(F21, Q11, T1, _V, _V_BYTE, 1, _PUB_M1_BYTE);
    #endif
    return;
}


void calculate_F31(unsigned char *F31, const unsigned char *T3, const unsigned char *T4, const unsigned char *Q11,
                    const unsigned char *Q21, const unsigned char *Q31){
    // F31 = Q31
    memcpy(F31, Q31, _PK_P3_BYTE >> 1);

    // F31 += Q21 * T3
    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T3[(_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T3, T3, (_PUB_M - 1));
    add_Q2_mul_T3_multab(F31, Q21, multab_T3);
    #else
    add_Q2_mul_T3(F31, Q21, T3);
    #endif

    // F31 += (Q11 + Q11_tr) * T4
    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T4[_V * (_PUB_M - 1) * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T4, T4, _V * (_PUB_M - 1));
    batch_2trimat_madd_multab(F31, Q11, multab_T4, _V, _V_BYTE, _PUB_M - 1, _PUB_M1_BYTE);
    #else
    batch_2trimat_madd(F31, Q11, T4, _V, _V_BYTE, _PUB_M - 1, _PUB_M1_BYTE);
    #endif
    return;
}


void calculate_F52(unsigned char *F52, const unsigned char *T1, const unsigned char *Q12, const unsigned char *Q22,
                    const unsigned char *Q52){
    uint8_t tmp[_PUB_M1_BYTE];

    #ifdef _MUL_WITH_MULTAB_
    uint8_t multab_T1[_V * 32] __attribute__((aligned(32)));
    gfv_generate_multabs(multab_T1, T1, _V);
    // F52 = T1_tr * Q12 * T1
    batch_quad_trimat_eval_multab(F52, Q12, multab_T1, _V, _PUB_M1_BYTE);
    // tmp = T1_tr * Q22
    gfmat_prod_multab(tmp, Q22, _PUB_M1_BYTE, _V, multab_T1);
    #else
    batch_quad_trimat_eval(F52, Q12, T1, _V, _PUB_M1_BYTE);
    gfmat_prod(tmp, Q22, _PUB_M1_BYTE, _V, T1);
    #endif
    // F51 += tmp + Q52
    gfv_add(F52, tmp, _PUB_M1_BYTE);
    gfv_add(F52, Q52, _PUB_M1_BYTE);
    return;
}


void calculate_F62(unsigned char *F62, const unsigned char *T1, const unsigned char *T3, const unsigned char *T4,
                    const unsigned char *Q12, const unsigned char *Q22, const unsigned char *Q32, const unsigned char *Q52,
                    const unsigned char *Q62){
    calculate_Q61(F62, T1, T3, T4, Q12, Q22, Q32, Q52);
    gfv_add(F62, Q62, _PK_P6_BYTE >> 1);
    return;
}
