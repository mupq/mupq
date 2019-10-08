/**
 * \file parameters.h
 * \brief Parameters of the ROLLO scheme
 */

#ifndef ROLLO_PARAMETER_H
#define ROLLO_PARAMETER_H

#include "api.h"

#define PARAM_M 101 /**< Parameter m of the scheme (finite field GF(2^m)) */
#define PARAM_N 47 /**< Parameter n of the scheme (code length) */
#define PARAM_W 5 /**< Parameter d of the scheme (weight of vectors) */
#define PARAM_W_R 6 /**< Parameter r of the scheme (weight of vectors) */
#define PARAM_DFR 30 /**< Decryption Failure Rate (2^-30) */
#define PARAM_SECURITY 128 /**< Expected security level */

#define PUBLIC_KEY_BYTES CRYPTO_PUBLICKEYBYTES
#define SECRET_KEY_BYTES CRYPTO_SECRETKEYBYTES
#define CIPHERTEXT_BYTES CRYPTO_CIPHERTEXTBYTES
#define SHARED_SECRET_BYTES CRYPTO_BYTES

#define FFI_VEC_R_BYTES 76 //Number of bytes to store R elements of GF2^m
#define FFI_VEC_N_BYTES 594 //Number of bytes to store N elements of GF2^m
#define SHA512_BYTES 64 /**< Size of SHA512 output */

#define SEEDEXPANDER_SEED_BYTES 40 /**< Seed size of the NIST seed expander */
#define SEEDEXPANDER_MAX_LENGTH 4294967295 /**< Max length of the NIST seed expander */

#endif
