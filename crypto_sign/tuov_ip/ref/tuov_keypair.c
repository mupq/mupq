/** \file tuov_keypair.c 
 *  \brief The standard implementations for key related functions in tuov_keypair.h   
*/

#include "tuov_keypair.h"
#include "tuov.h"
#include "tuov_keypair_computation.h"
#include "tuov_blas.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "utils_prng.h"
#include "utils_hash.h"
#include "utils_malloc.h"

#if 96 < _V
#define _MALLOC_
#endif




int _generate_keypair( pk_t *pk, sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed){
    memcpy( sk->sk_seed, sk_seed, LEN_SKSEED );

    // prng for sk
    hash_ctx hctx;
    hash_init(&hctx);
    hash_update(&hctx, sk_seed, LEN_SKSEED);
    hash_final_digest(sk->ST, sizeof(sk->ST), &hctx);
    #ifdef _USE_GF16
    gfv_set_ele(sk->ST+T3_BIAS, _PUB_M - 1, 0x0);
    #endif

    // prng for pk
    prng_publicinputs_t prng1;
    prng_set_publicinputs(&prng1, pk_seed);
    // generate P1, P2, P3
    prng_gen_publicinputs(&prng1, pk->pk, _PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE);
    // generate P52, P62
    prng_gen_publicinputs(&prng1, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    prng_gen_publicinputs(&prng1, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);

    // Q123 stores Q11, Q21, Q31 in order
    unsigned char Q123[(_PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE) >> 1];
    unsigned char *Q11 = Q123;
    unsigned char *Q21 = Q123 + (_PK_P1_BYTE >> 1);
    unsigned char *Q31 = Q123 + ((_PK_P1_BYTE + _PK_P2_BYTE) >> 1);
    // line 5 of KeyGen
    memcpy(Q11, pk->pk + P12_BIAS, _PK_P1_BYTE >> 1);
    memcpy(Q21, pk->pk + P22_BIAS, _PK_P2_BYTE >> 1);
    memcpy(Q31, pk->pk + P32_BIAS, _PK_P3_BYTE >> 1);
    multiple_S(Q123, sk->ST, N_TRIANGLE_TERMS(_V) + _V * _PUB_M);
    gfv_add(Q11, pk->pk + P11_BIAS, _PK_P1_BYTE >> 1);
    gfv_add(Q21, pk->pk + P21_BIAS, _PK_P2_BYTE >> 1);
    gfv_add(Q31, pk->pk + P31_BIAS, _PK_P3_BYTE >> 1);

    // line 8 - 12 of KeyGen
    calculate_Q51(pk->pk + P51_BIAS, sk->ST + T1_BIAS, Q11, Q21);
    calculate_Q61(pk->pk + P61_BIAS, sk->ST + T1_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31, pk->pk + P51_BIAS);
    calculate_Q91(pk->pk + P91_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31, pk->pk + P51_BIAS, pk->pk + P61_BIAS);
    calculate_Q91(pk->pk + P92_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, 
                pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS, 
                pk->pk + P52_BIAS, pk->pk + P62_BIAS);

    // Q569 stores Q52, Q62, Q92 in order
    unsigned char Q569[(_PK_P5_BYTE + _PK_P6_BYTE + _PK_P9_BYTE) >> 1];
    unsigned char *Q52 = Q569;
    unsigned char *Q62 = Q569 + (_PK_P5_BYTE >> 1);
    unsigned char *Q92 = Q569 + ((_PK_P5_BYTE + _PK_P6_BYTE) >> 1);
    // line 13 of KeyGen
    memcpy(Q52, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    memcpy(Q62, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);
    memcpy(Q92, pk->pk + P92_BIAS, _PK_P9_BYTE >> 1);
    multiple_S(Q569, sk->ST, N_TRIANGLE_TERMS(_PUB_M));
    gfv_add(pk->pk + P51_BIAS, Q52, _PK_P5_BYTE >> 1);
    gfv_add(pk->pk + P61_BIAS, Q62, _PK_P6_BYTE >> 1);
    gfv_add(pk->pk + P91_BIAS, Q92, _PK_P9_BYTE >> 1);

    // now public key is done, we need to generate the F's in secret key

    // F1 = Q1
    memcpy(sk->F, Q123, _PK_P1_BYTE >> 1);
    memcpy(sk->F + P12_BIAS, pk->pk + P12_BIAS, _PK_P1_BYTE >> 1);
    // line 5-12 in ExpandSK
    calculate_F21(sk->F + P21_BIAS, sk->ST + T1_BIAS, Q11, Q21);
    calculate_F31(sk->F + P31_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31);
    calculate_F21(sk->F + P22_BIAS, sk->ST + T1_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS);
    calculate_F31(sk->F + P32_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS);
    calculate_F52(sk->F + P51_BIAS, sk->ST + T1_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P52_BIAS);
    calculate_F62(sk->F + P51_BIAS + (_PK_P5_BYTE >> 1), sk->ST + T1_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, 
                    pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS, pk->pk + P52_BIAS, pk->pk + P62_BIAS);
    
    return 0;
}

int generate_keypair(pk_t *pk, sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed){
    int r = _generate_keypair(pk, sk, pk_seed, sk_seed);
    combine_F(sk->F);
    combine_P(pk->pk);
    return r;
}

int generate_keypair_pkc( cpk_t *rpk, sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed ){
    #if defined(_MALLOC_)
    pk_t *pk = malloc(sizeof(pk_t));
    if (NULL == pk) {
        return -1;
    }
    #else
    pk_t _pk;
    pk_t *pk = &_pk;
    #endif
    int r = _generate_keypair(pk, sk, pk_seed, sk_seed);
    combine_F(sk->F);
    memcpy(rpk->pk_seed, pk_seed, LEN_PKSEED);
    memcpy(rpk->P51, pk->pk + P51_BIAS, _PK_P5_BYTE >> 1);
    memcpy(rpk->P61, pk->pk + P61_BIAS, _PK_P6_BYTE >> 1);
    memcpy(rpk->P9, pk->pk + P91_BIAS, _PK_P9_BYTE);

    #if defined(_MALLOC_)
    free(pk);
    #endif
    return r;
}

int generate_keypair_pkc_skc( cpk_t *rpk, csk_t *rsk, const unsigned char *pk_seed, const unsigned char *sk_seed ){
    pk_t *pk = malloc(sizeof(pk_t));
    sk_t *sk = malloc(sizeof(sk_t));
    memcpy( sk->sk_seed, sk_seed, LEN_SKSEED );

    // prng for sk
    hash_ctx hctx;
    hash_init(&hctx);
    hash_update(&hctx, sk_seed, LEN_SKSEED);
    hash_final_digest(sk->ST, sizeof(sk->ST), &hctx);
    #ifdef _USE_GF16
    gfv_set_ele(sk->ST+T3_BIAS, _PUB_M - 1, 0x0);
    #endif

    // prng for pk
    prng_publicinputs_t prng1;
    prng_set_publicinputs(&prng1, pk_seed);
    // generate P1, P2, P3
    prng_gen_publicinputs(&prng1, pk->pk, _PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE);
    // generate P52, P62
    prng_gen_publicinputs(&prng1, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    prng_gen_publicinputs(&prng1, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);

    // Q123 stores Q11, Q21, Q31 in order
    unsigned char Q123[(_PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE) >> 1];
    unsigned char *Q11 = Q123;
    unsigned char *Q21 = Q123 + (_PK_P1_BYTE >> 1);
    unsigned char *Q31 = Q123 + ((_PK_P1_BYTE + _PK_P2_BYTE) >> 1);
    // line 5 of KeyGen
    memcpy(Q11, pk->pk + P12_BIAS, _PK_P1_BYTE >> 1);
    memcpy(Q21, pk->pk + P22_BIAS, _PK_P2_BYTE >> 1);
    memcpy(Q31, pk->pk + P32_BIAS, _PK_P3_BYTE >> 1);
    multiple_S(Q123, sk->ST, N_TRIANGLE_TERMS(_V) + _V * _PUB_M);
    gfv_add(Q11, pk->pk + P11_BIAS, _PK_P1_BYTE >> 1);
    gfv_add(Q21, pk->pk + P21_BIAS, _PK_P2_BYTE >> 1);
    gfv_add(Q31, pk->pk + P31_BIAS, _PK_P3_BYTE >> 1);

    // line 8 - 12 of KeyGen
    calculate_Q51(pk->pk + P51_BIAS, sk->ST + T1_BIAS, Q11, Q21);
    calculate_Q61(pk->pk + P61_BIAS, sk->ST + T1_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31, pk->pk + P51_BIAS);
    calculate_Q91(pk->pk + P91_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31, pk->pk + P51_BIAS, pk->pk + P61_BIAS);
    calculate_Q91(pk->pk + P92_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, 
                pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS, 
                pk->pk + P52_BIAS, pk->pk + P62_BIAS);

    // Q569 stores Q52, Q62, Q92 in order
    unsigned char Q569[(_PK_P5_BYTE + _PK_P6_BYTE + _PK_P9_BYTE) >> 1];
    unsigned char *Q52 = Q569;
    unsigned char *Q62 = Q569 + (_PK_P5_BYTE >> 1);
    unsigned char *Q92 = Q569 + ((_PK_P5_BYTE + _PK_P6_BYTE) >> 1);
    // line 13 of KeyGen
    memcpy(Q52, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    memcpy(Q62, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);
    memcpy(Q92, pk->pk + P92_BIAS, _PK_P9_BYTE >> 1);
    multiple_S(Q569, sk->ST, N_TRIANGLE_TERMS(_PUB_M));
    gfv_add(pk->pk + P51_BIAS, Q52, _PK_P5_BYTE >> 1);
    gfv_add(pk->pk + P61_BIAS, Q62, _PK_P6_BYTE >> 1);
    gfv_add(pk->pk + P91_BIAS, Q92, _PK_P9_BYTE >> 1);

    // now public key is done, we don't need to generate the F's in secret key
    memcpy(rsk->pk_seed, pk_seed, LEN_PKSEED);
    memcpy(rsk->sk_seed, sk_seed, LEN_SKSEED);
    memcpy(rpk->pk_seed, pk_seed, LEN_PKSEED);
    memcpy(rpk->P51, pk->pk + P51_BIAS, _PK_P5_BYTE >> 1);
    memcpy(rpk->P61, pk->pk + P61_BIAS, _PK_P6_BYTE >> 1);
    memcpy(rpk->P9, pk->pk + P91_BIAS, _PK_P9_BYTE);
    free(sk);
    free(pk);
    return 0;
}

int expand_pk( pk_t *pk, const cpk_t *cpk ){
    prng_publicinputs_t prng1;
    prng_set_publicinputs(&prng1, cpk->pk_seed);

    // generate P1, P2, P3
    prng_gen_publicinputs(&prng1, pk->pk, _PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE);
    // generate P52, P62
    prng_gen_publicinputs(&prng1, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    prng_gen_publicinputs(&prng1, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);

    memcpy(pk->pk + P51_BIAS, cpk->P51, sizeof(cpk->P51));
    memcpy(pk->pk + P61_BIAS, cpk->P61, sizeof(cpk->P61));
    memcpy(pk->pk + P91_BIAS, cpk->P9, sizeof(cpk->P9));
    combine_P(pk->pk);
    return 0;
}

int expand_sk( sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed ){
    pk_t *pk = malloc(sizeof(pk_t));
    memcpy( sk->sk_seed, sk_seed, LEN_SKSEED );

    // prng for sk
    hash_ctx hctx;
    hash_init(&hctx);
    hash_update(&hctx, sk_seed, LEN_SKSEED);
    hash_final_digest(sk->ST, sizeof(sk->ST), &hctx);
    #ifdef _USE_GF16
    gfv_set_ele(sk->ST+T3_BIAS, _PUB_M - 1, 0x0);
    #endif

    // prng for pk
    prng_publicinputs_t prng1;
    prng_set_publicinputs(&prng1, pk_seed);
    // generate P1, P2, P3
    prng_gen_publicinputs(&prng1, pk->pk, _PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE);
    // generate P52, P62
    prng_gen_publicinputs(&prng1, pk->pk + P52_BIAS, _PK_P5_BYTE >> 1);
    prng_gen_publicinputs(&prng1, pk->pk + P62_BIAS, _PK_P6_BYTE >> 1);

    // Q123 stores Q11, Q21, Q31 in order
    unsigned char Q123[(_PK_P1_BYTE + _PK_P2_BYTE + _PK_P3_BYTE) >> 1];
    unsigned char *Q11 = Q123;
    unsigned char *Q21 = Q123 + (_PK_P1_BYTE >> 1);
    unsigned char *Q31 = Q123 + ((_PK_P1_BYTE + _PK_P2_BYTE) >> 1);
    // line 5 of KeyGen
    memcpy(Q11, pk->pk + P12_BIAS, _PK_P1_BYTE >> 1);
    memcpy(Q21, pk->pk + P22_BIAS, _PK_P2_BYTE >> 1);
    memcpy(Q31, pk->pk + P32_BIAS, _PK_P3_BYTE >> 1);
    multiple_S(Q123, sk->ST, N_TRIANGLE_TERMS(_V) + _V * _PUB_M);
    gfv_add(Q11, pk->pk + P11_BIAS, _PK_P1_BYTE >> 1);
    gfv_add(Q21, pk->pk + P21_BIAS, _PK_P2_BYTE >> 1);
    gfv_add(Q31, pk->pk + P31_BIAS, _PK_P3_BYTE >> 1);

    // F1 = Q1
    memcpy(sk->F, Q123, _PK_P1_BYTE >> 1);
    memcpy(sk->F + P12_BIAS, pk->pk + P12_BIAS, _PK_P1_BYTE >> 1);
    // line 5-12 in ExpandSK
    calculate_F21(sk->F + P21_BIAS, sk->ST + T1_BIAS, Q11, Q21);
    calculate_F31(sk->F + P31_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, Q11, Q21, Q31);
    calculate_F21(sk->F + P22_BIAS, sk->ST + T1_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS);
    calculate_F31(sk->F + P32_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS);
    calculate_F52(sk->F + P51_BIAS, sk->ST + T1_BIAS, pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P52_BIAS);
    calculate_F62(sk->F + P51_BIAS + (_PK_P5_BYTE >> 1), sk->ST + T1_BIAS, sk->ST + T3_BIAS, sk->ST + T4_BIAS, 
                    pk->pk + P12_BIAS, pk->pk + P22_BIAS, pk->pk + P32_BIAS, pk->pk + P52_BIAS, pk->pk + P62_BIAS);
    
    combine_F(sk->F);
    free(pk);
    return 0;
}
