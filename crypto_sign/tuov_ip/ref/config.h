/** \file config.h
*  \brief This file contains all implementation options that do not affect testvectors
*   params.h contains other options.
*/

#ifndef CONFIG_H_
#define CONFIG_H_


/* The following macros choose optimizations for basic linear algebra functions. */
#if (!defined(_BLAS_AVX2_))&&(!defined(_BLAS_SSE_))&&(!defined(_BLAS_UINT64_))
//#define _BLAS_AVX2_
//#define _BLAS_SSE_
#define _BLAS_UINT64_
#endif

/* GF multiplication with multiplication tables */
#ifndef _BLAS_UINT64_
#define _MUL_WITH_MULTAB_
#endif

/* choose implementations for SHAKE256, AES128CTR, and randombytes() functions */
#define _UTILS_OPENSSL_
//#define _UTILS_SUPERCOP_
//#define _UTILS_PQM4_

/* x86 aes instrucitons */
//#define _UTILS_AESNI_


#endif

