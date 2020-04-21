/**
 * \file api.h
 * \brief NIST KEM API used by the RQC_KEM IND-CCA2 scheme
 */


#ifndef API_H
#define API_H

#define CRYPTO_ALGNAME "RQC-192"

#define CRYPTO_SECRETKEYBYTES 2893
#define CRYPTO_PUBLICKEYBYTES 2853
#define CRYPTO_BYTES 64
#define CRYPTO_CIPHERTEXTBYTES 5690

// As a technicality, the public key is appended to the secret key in order to respect the NIST API.
// Without this constraint, CRYPTO_SECRETKEYBYTES would be defined as 40


int crypto_kem_keypair(unsigned char* pk, unsigned char* sk);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk);
int crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk);

#endif
