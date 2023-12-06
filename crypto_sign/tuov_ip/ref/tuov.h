/** \file tuov.h
 *  \brief functions used in sign.c 
 */

#ifndef _TUOV_H_
#define _TUOV_H_

#include "params.h"
#include "tuov_keypair.h"

#include <stdint.h>

#ifdef __cplusplus
extern  "C" {
#endif


////////////////////////  Key pair generation  ///////////////////////////////////



/** \brief Key pair generation for classic tuov.
 * 
 *  \param[out] pk        - the public key.
 *  \param[out] sk        - the secret key.
 *  \param[in]  pk_seed   - seed for generating parts of the public key.
 *  \param[in]  sk_seed   - seed for generating the secret key.
 *  \return 0 for success. -1 otherwise.
 */
int generate_keypair( pk_t *pk, sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed);


/** \brief Generate key pair for pkc
 * 
 *  \param[out] pk        - the compressed public key.
 *  \param[out] sk        - the secret key.
 *  \param[in]  pk_seed   - seed for generating parts of public key.
 *  \param[in]  sk_seed   - seed for generating secret key.
 *  \return 0 for success. -1 otherwise.
 */
int generate_keypair_pkc( cpk_t *pk, sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed );


/** \brief Generate key pair for pkc+skc
 * 
 *  \param[out] pk        - the compressed public key.
 *  \param[out] sk        - the compressed secret key.
 *  \param[in]  pk_seed   - seed for generating parts of the public key.
 *  \param[in]  sk_seed   - seed for generating the secret key.
 *  \return 0 for success. -1 otherwise.
 */
int generate_keypair_pkc_skc( cpk_t *pk, csk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed );



///////////////////////////// Sign and Verify ////////////////////////////////



/** \brief Signing function for classic secret key.
 * 
 *  \param[out] signature - the signature.
 *  \param[in]  sk        - the secret key.
 *  \param[in]  message   - the message to be signed.
 *  \param[in]  mlen      - the length of the message.
 *  \return 0 for success. -1 otherwise.
 */
int tuov_sign( uint8_t *signature, const sk_t *sk, const uint8_t *message, unsigned mlen );


/** \brief Verifying function.
 * 
 *  \param[in]  message   - the message.
 *  \param[in]  mlen      - the length of the message.
 *  \param[in]  signature - the signature.
 *  \param[in]  pk        - the public key.
 *  \return 0 for successful verified. -1 for failed verification.
 */
int tuov_verify( const uint8_t *message, unsigned mlen, const uint8_t *signature, const pk_t *pk );


/** \brief Signing function for compressed secret key.
 * 
 *  \param[out] signature - the signature.
 *  \param[in]  sk        - the compressed secret key.
 *  \param[in]  message   - the message to be signed.
 *  \param[in]  mlen      - the length of the message.
 *  \return 0 for success. -1 otherwise.
 */ 
int tuov_expand_and_sign( uint8_t *signature, const csk_t *sk, const uint8_t *message, unsigned mlen );


/** \brief Verifying function for compressed public keys.
 * 
 *  \param[in]  message   - the message to be signed.
 *  \param[in]  mlen      - the length of the message.
 *  \param[in]  signature - the signature.
 *  \param[in]  pk        - the compressed public key.
 *  \return 0 for successful verified. -1 for failed verification.
 */ 
int ov_expand_and_verify( const uint8_t *message, unsigned mlen, const uint8_t *signature, const cpk_t *pk );



/// Key conversion : compressed key pair -> classic key pair //////////////////////////////



/** \brief converting formats of public keys : from compressed to classic public key
 *
 *  \param[out] pk       - the classic public key.
 *  \param[in]  cpk      - the compressed public key.
 *  \return always return 0
 */
int expand_pk( pk_t *pk, const cpk_t *cpk );



/** \brief Generate secret key
 * 
 *  \param[out] sk        - the secret key (pre-allocated).
 *  \param[in]  pk_seed   - seed for generating parts of the public key.
 *  \param[in]  sk_seed   - seed for generating the secret key.
 *  \return 0 for success. -1 otherwise.
 */
int expand_sk( sk_t *sk, const unsigned char *pk_seed, const unsigned char *sk_seed );



////////////////////////  Public map  ///////////////////////////////////



/** \brief Public-key evaluattion
 *
 *  \param[out] z         - results of the evluation of public polynomials at the w (pre-allocated).
 *  \param[in]  pk        - the classic public key.
 *  \param[in]  w         - the input vector w.
 */
void tuov_publicmap( unsigned char *z, const unsigned char *pk, const unsigned char *w );



#ifdef __cplusplus
}
#endif


#endif
