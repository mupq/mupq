/** \file tuov.c 
 *  \brief The standard implementations for sign and verify functions in tuov.h   
*/

#include "tuov.h"
#include "tuov_blas.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "utils_prng.h"
#include "utils_hash.h"
#include "utils_malloc.h"


#define MAX_ATTEMPT_VINEGAR  256

int tuov_sign( uint8_t *signature, const sk_t *sk, const uint8_t *message, unsigned mlen ){

    uint8_t t[_PUB_M_BYTE];
    uint8_t lambda[_PUB_M1_BYTE + (_V_BYTE - _O_BYTE)];
    uint8_t M[_PUB_M * _V_BYTE];
    uint8_t v[_V_BYTE];
    uint8_t y[_PUB_M_BYTE];
    uint8_t x;
    uint8_t L[_PUB_M * _PUB_M_BYTE];
    uint8_t *random_byte_for_MTv = lambda + _PUB_M1_BYTE;

    hash_ctx h_m_secret;
    hash_ctx h_lambda_copy;

    // compute t = H(M) 
    hash_init  (&h_m_secret);
    hash_update(&h_m_secret, message, mlen);
    hash_ctx_copy(&h_lambda_copy, &h_m_secret);
    hash_final_digest(t, _PUB_M_BYTE, &h_m_secret);
    hash_update(&h_lambda_copy, sk->sk_seed, LEN_SKSEED);

    // compute t_tilde
    uint8_t t_tmp[_PUB_M1_BYTE];
    // t_tmp = S * t2
    gfmat_prod(t_tmp, sk->ST, _PUB_M1_BYTE, _PUB_M1, t + _PUB_M1_BYTE);
    // t1 += t_tmp
    gfv_add(t, t_tmp, _PUB_M1_BYTE);

    unsigned n_attempt = 0;
    while (MAX_ATTEMPT_VINEGAR > n_attempt){
        uint8_t ctr = n_attempt & 0xff;
        n_attempt++;

        // compute lambda = H(M||sk_seed||ctr)
        hash_ctx h_lambda;
        hash_ctx_copy(&h_lambda, &h_lambda_copy);
        hash_update(&h_lambda, &ctr, 1);
        hash_final_digest(lambda, _PUB_M1_BYTE + (_V_BYTE - _O_BYTE), &h_lambda);

        // compute M
        // linear comb of F21, F31
        gfv_set_zero(M, _O_BYTE * _V);
        for (long i = 0; i < _PUB_M1; i++){
            gfv_madd(M, sk->F + P21_BIAS + i * _O_BYTE * _V, gfv_get_ele(lambda, i), _O_BYTE * _V);
        }

        memcpy(v, random_byte_for_MTv, _V_BYTE - _O_BYTE);
        uint8_t e1[_PUB_M_BYTE];
        gfmat_prod(e1, M, _O_BYTE, _V - _O, v);
        gfv_set_ele(e1, _O - 1, 0x1 ^ gfv_get_ele(e1, _O - 1));
        uint8_t *_M = M + (_O_BYTE) * (_V - _O);
        unsigned succ = gfmat_gaussian_elim(_M, e1, _O);
        if (!succ) continue;
        gfmat_back_substitute(e1, _M, _O);
        memcpy(v + _V_BYTE - _O_BYTE, e1, _O_BYTE);

        uint8_t a; 
        while (1) {
            randombytes(&a, 1);
            #ifdef _USE_GF16
            a &= 0xf;
            // a cannot be 0
            if (gf16_is_nonzero(a)) {break;}
            #else 
            if (gf256_is_nonzero(a)) {break;}
            #endif
        }

        // v = a * v
        gfv_mul_scalar(v, a, _V_BYTE);

        uint8_t vtrf2[_PUB_M_BYTE];
        #if defined(_MUL_WITH_MULTAB_)
        // generate multab for v
        uint8_t multabs[(_V) * 32] __attribute__((aligned(32)));
        gfv_generate_multabs(multabs, v, _V);
        for (long i = 0; i < _O; i++){
            // compute the i-th row of L
            gfmat_prod_multab(L + i * _O_BYTE, sk->F + P21_BIAS + i * _O_BYTE * _V, _O_BYTE, _V, multabs);
            gfv_set_ele(vtrf2, i, gfv_get_ele(L + i * _O_BYTE, _O - 1));
            gfv_set_ele(L + i * _O_BYTE, _O - 1, 0x0);
        }
        batch_quad_trimat_eval_multab(y, sk->F, multabs, _V, _PUB_M_BYTE);
        #else
        for (long i = 0; i < _O; i++){
            // compute the i-th row of L
            gfmat_prod(L + i * _O_BYTE, sk->F + P21_BIAS + i * _O_BYTE * _V, _O_BYTE, _V, v);
            gfv_set_ele(vtrf2, i, gfv_get_ele(L + i * _O_BYTE, _O - 1));
            gfv_set_ele(L + i * _O_BYTE, _O - 1, 0x0);
        }
        batch_quad_trimat_eval(y, sk->F, v, _V, _PUB_M_BYTE);
        #endif
        gfmat_transpose_oxo(L);

        // compute x 
        uint8_t x_tmp[_PUB_M1_BYTE];
        memcpy(x_tmp, t, _PUB_M1_BYTE);
        gfv_add(x_tmp, y, _PUB_M1_BYTE);
        // compute the dot product of lambda and x_tmp
        x = 0x0;
        for (unsigned i = 0; i < _PUB_M1; i++){
            x ^= gf_mul(gfv_get_ele(x_tmp, i), gfv_get_ele(lambda, i));
        }
        x = gf_mul(x, gf_inv(a));
        uint8_t x2 = gf_squ(x);

        // compute x^2 * F52 and add it to y
        gfv_madd(y+_PUB_M1_BYTE, sk->F + P51_BIAS, x2, _PUB_M1_BYTE);
        // compute x * v_tr * F2 and add it to y
        gfv_madd(y, vtrf2, x, _PUB_M_BYTE);
        // y = y - t
        gfv_add(y, t, _PUB_M_BYTE);
        // compute x * F62 and add it to L (L is col-major) 
        for (unsigned i = 0; i < _PUB_M - 1; i++){
            gfv_madd(L + (i * _PUB_M_BYTE + _PUB_M1_BYTE), sk->F + P51_BIAS + (_PK_P5_BYTE >> 1) + _PUB_M1_BYTE * i, x, _PUB_M1_BYTE);
        }

        // find an index k such that lambda[k] nonzero.
        unsigned int k = 0;
        for (long i = _PUB_M1 - 1; i >= 0; i--){
            if (gfv_get_ele(lambda, i)) k = i;
        }
        uint8_t *L_last_col = L + (unsigned int)(_PUB_M_BYTE * (_PUB_M - 1));
        gfv_set_zero(L_last_col, _PUB_M_BYTE);
        gfv_set_ele(L_last_col, k, 0x1);

        // solve the equation, the result is stored in y
        // the last element of y must be zero.
        int succ2 = gfmat_gaussian_elim(L, y, _PUB_M);
        if (!succ2) continue;
        gfmat_back_substitute(y, L, _PUB_M);
        break;
    }
    hash_final_digest(lambda, _PUB_M1_BYTE, &h_lambda_copy); // free
    if (n_attempt >= MAX_ATTEMPT_VINEGAR) return -1;

    uint8_t *w = signature;
    for (long i = 0; i < _PUB_M - 1; i++){
        gfv_set_ele(w+_V_BYTE, i+1, gfv_get_ele(y, i));
    }
    gfv_set_ele(w + _V_BYTE, 0, x ^ gfv_dot(y, sk->ST + T3_BIAS, _SK_T3_BYTE));
    gfmat_prod(w, sk->ST + T4_BIAS, _V_BYTE, _PUB_M, y);

    gfv_madd(w, sk->ST + T1_BIAS, x, _V_BYTE);
    gfv_add(w, v, _V_BYTE);
    return 0;
}




/** \brief check if H( message || salt) == digest_ck
 *
 *  \param[in] message      - the message
 *  \param[in] mlen         - the length of the message
 *  \param[in] salt         - the salt used in the hash
 *  \param[in] digest_ck    - the digest to check
 *  \return 0 if H( message || salt ) == digest_ck, -1 otherwise.
 */
static
int _tuov_verify( const uint8_t *message, unsigned mlen, const unsigned char *digest_ck ) {
    /* compute Hash( message || salt ) */
    unsigned char correct[_PUB_M_BYTE];
    hash_ctx hctx;
    hash_init(&hctx);
    hash_update(&hctx, message, mlen);
    hash_final_digest(correct, _PUB_M_BYTE, &hctx);

    /* check consistency */ 
    unsigned char cc = 0;
    for (unsigned i = 0; i < _PUB_M_BYTE; i++) {
        cc |= (digest_ck[i] ^ correct[i]);
    }
    return (0 == cc) ? 0 : -1;
}


int tuov_verify( const uint8_t *message, unsigned mlen, const uint8_t *signature, const pk_t *pk ) {
    unsigned char digest_ck[_PUB_M_BYTE];
    tuov_publicmap( digest_ck, pk->pk, signature );
    return _tuov_verify( message, mlen, digest_ck );
}


#if defined(_TUOV_PKC_SKC)
int tuov_expand_and_sign( uint8_t *signature, const csk_t *csk, const uint8_t *message, unsigned mlen ) {
    sk_t _sk;
    sk_t *sk = &_sk;
    expand_sk( sk, csk->pk_seed, csk->sk_seed );    // generating classic secret key.
    int r = tuov_sign( signature, sk, message, mlen );
    return r;
}
#endif

#if defined(_TUOV_PKC) || defined(_TUOV_PKC_SKC)
int ov_expand_and_verify( const uint8_t *message, unsigned mlen, const uint8_t *signature, const cpk_t *cpk ) {
    pk_t _pk;
    pk_t *pk = &_pk;
    expand_pk( pk, cpk );
    return tuov_verify( message, mlen, signature, pk );
}
#endif
