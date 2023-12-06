#ifndef _API_H_
#define _API_H_


#include "params.h"


#define CRYPTO_SECRETKEYBYTES   TUOV_SECRETKEYBYTES
#define CRYPTO_PUBLICKEYBYTES   TUOV_PUBLICKEYBYTES
#define CRYPTO_BYTES            TUOV_SIGNATUREBYTES
#define CRYPTO_ALGNAME          TUOV_ALGNAME

#ifdef __cplusplus
extern  "C" {
#endif

#include <stddef.h>

/** Generates a pair of public key and secret key.
 *
 *  \param[out] pk      - pointer to output public key (allocated array of TUOV_PUBLICKEYBYTES bytes)
 *  \param[out] sk      - pointer to output private key (allocated array of TUOV_SECRETKEYBYTES bytes)
 *  \return 0 for success. -1 otherwise.
 */
int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk);


/** Compute signed message.
 *
 *  \param[out] sm      - pointer to output signed message (allocated array with TUOV_SIGNATUREBYTES + mlen bytes)
 *  \param[out] smlen   - pointer to output length of signed message
 *  \param[in] m        - pointer to message to be signed
 *  \param[in] mlen     - length of message
 *  \param[in] sk       - pointer to the secret key
 *  \return 0 for success. -1 otherwise.
 */
int
crypto_sign(unsigned char *sm, size_t *smlen,
            const unsigned char *m, size_t mlen,
            const unsigned char *sk);


/** Verify signed message.
 *
 *  \param[out] m       - pointer to output message, allocated array with at least (smlen - TUOV_SIGNATUREBYTES) bytes
 *  \param[out] mlen    - pointer to output length of message
 *  \param[in] sm       - pointer to signed message
 *  \param[in] smlen    - length of signed message
 *  \param[in] pk       - pointer to the public key
 *  \return 0 for successful verified. -1 for failed verification.
 */
int
crypto_sign_open(unsigned char *m, size_t *mlen,
                 const unsigned char *sm, size_t smlen,
                 const unsigned char *pk);

#ifdef __cplusplus
}
#endif


#endif
