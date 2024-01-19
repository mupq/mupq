/** \file   sign.c
 *  \brief  the implementations for functions in api.h
 */
#include <stdlib.h>
#include <string.h>

#include "params.h"
#include "tuov_keypair.h"
#include "tuov.h"

#include "api.h"

#include "utils_randombytes.h"      /* for randombytes() */


int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
    unsigned char sk_seed[LEN_SKSEED];
    unsigned char pk_seed[LEN_PKSEED];
    randombytes( sk_seed, LEN_SKSEED );
    randombytes( pk_seed, LEN_PKSEED );

    #if defined _TUOV_CLASSIC
    int r = generate_keypair( (pk_t *) pk, (sk_t *) sk, pk_seed, sk_seed );
    #elif defined _TUOV_PKC
    int r = generate_keypair_pkc( (cpk_t *) pk, (sk_t *) sk, pk_seed, sk_seed );
    #elif defined _TUOV_PKC_SKC
    int r = generate_keypair_pkc_skc( (cpk_t *) pk, (csk_t *) sk, pk_seed, sk_seed );
    #else
    #error here
    #endif

    return r;
}


int
crypto_sign(unsigned char *sm, size_t *smlen, const unsigned char *m, size_t mlen, const unsigned char *sk)
{
    int r = -1;

    #if defined _TUOV_CLASSIC
    r = tuov_sign( sm + mlen, (const sk_t *)sk, m, mlen );
    #elif defined _TUOV_PKC
    r = tuov_sign( sm + mlen, (const sk_t *)sk, m, mlen );
    #elif defined _TUOV_PKC_SKC
    r = tuov_expand_and_sign( sm + mlen, (const csk_t *)sk, m, mlen );
    #else
    #error here
    #endif

    memcpy( sm, m, mlen );
    smlen[0] = mlen + TUOV_SIGNATUREBYTES;

    return r;
}


int
crypto_sign_open(unsigned char *m, size_t *mlen, const unsigned char *sm, size_t smlen, const unsigned char *pk)
{
    if ( TUOV_SIGNATUREBYTES > smlen ) {
        return -1;
    }
    unsigned long mesg_len = smlen - TUOV_SIGNATUREBYTES;
    int r = -1;

    #if defined _TUOV_CLASSIC
    r = tuov_verify( sm, mesg_len, sm + mesg_len, (const pk_t *)pk );
    #elif defined _TUOV_PKC
    r = ov_expand_and_verify( sm, mesg_len, sm + mesg_len, (const cpk_t *)pk );
    #elif defined _TUOV_PKC_SKC
    r = ov_expand_and_verify( sm, mesg_len, sm + mesg_len, (const cpk_t *)pk );
    #else
    #error here
    #endif

    memcpy( m, sm, mesg_len );
    mlen[0] = mesg_len;

    return r;
}
