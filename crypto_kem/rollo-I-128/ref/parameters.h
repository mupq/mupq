/**
 * \file parameters.h
 * \brief Parameters of the ROLLO scheme
 */

#ifndef ROLLOI_PARAMETER_H
#define ROLLOI_PARAMETER_H

#define ROLLOI_PARAM_M 67 /**< Parameter m of the scheme (finite field GF(2^m)) */
#define ROLLOI_PARAM_N 83 /**< Parameter n of the scheme (code length) */
#define ROLLOI_PARAM_D 8 /**< Parameter d of the scheme (weight of vectors) */
#define ROLLOI_PARAM_R 7 /**< Parameter r of the scheme (weight of vectors) */
#define ROLLOI_PARAM_DFR 28 /**< Decryption Failure Rate (2^-30) */
#define ROLLOI_PARAM_SECURITY 128 /**< Expected security level */

#define ROLLOI_PUBLIC_KEY_BYTES 696
#define ROLLOI_SECRET_KEY_BYTES 40
#define ROLLOI_CIPHERTEXT_BYTES 696
#define ROLLOI_SHARED_SECRET_BYTES 64

#define ROLLOI_RBC_VEC_R_BYTES 59 //Number of bytes to store R elements of GF2^m
#define ROLLOI_RBC_VEC_N_BYTES 696 //Number of bytes to store N elements of GF2^m

#define SHA512_BYTES 64 /**< Size of SHA512 output */

#define SEEDEXPANDER_SEED_BYTES 40 /**< Seed size of the NIST seed expander */
#define SEEDEXPANDER_MAX_LENGTH 4294967295 /**< Max length of the NIST seed expander */

#endif
