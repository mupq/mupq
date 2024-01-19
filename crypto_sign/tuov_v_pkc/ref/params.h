/** \file param.h
*  \brief parameter file contains all options that affect testvectors
*   config.h contains implementation options.
*/

#ifndef _PARAMS_H_
#define _PARAMS_H_

#if (!defined(_TUOV16_160_64))&&(!defined(_TUOV256_112_44))&&(!defined(_TUOV256_184_72))&&(!defined(_TUOV256_244_96))
// #define _TUOV16_160_64
// #define _TUOV256_112_44
// #define _TUOV256_184_72
#define _TUOV256_244_96
#endif

#if (!defined(_TUOV_CLASSIC))&&(!defined(_TUOV_PKC))&&(!defined(_TUOV_PKC_SKC))
// #define _TUOV_CLASSIC
#define _TUOV_PKC
//#define _TUOV_PKC_SKC
#endif


#if defined _TUOV16_160_64
#define _USE_GF16
#define _GFSIZE 16
#define _PUB_N 160
#define _PUB_M 64
#define _PUB_M1 32
#define _HASH_LEN 32

#elif defined _TUOV256_112_44
#define _GFSIZE 256
#define _PUB_N 112
#define _PUB_M 44
#define _PUB_M1 22
#define _HASH_LEN 32

#elif defined _TUOV256_184_72
#define _GFSIZE 256
#define _PUB_N 184
#define _PUB_M 72
#define _PUB_M1 36
#define _HASH_LEN 48

#elif defined _TUOV256_244_96
#define _GFSIZE 256
#define _PUB_N 244
#define _PUB_M 96
#define _PUB_M1 48
#define _HASH_LEN 64
#else
#error here
#endif

/* make sure m = 2 * m1 */
#if (_PUB_M1 != (_PUB_M >> 1))
#error m shound equal to 2 * m1
#endif

#define _O (_PUB_M)
#define _V ((_PUB_N)-(_PUB_M))


#define STR1(x) #x
#define THE_NAME(gf,n,m) "TUOV(" STR1(gf) "," STR1(n) "," STR1(m) ")"
#define TUOV_PARAMNAME THE_NAME(_GFSIZE,_PUB_N,_PUB_M)


#ifdef _USE_GF16
// GF16
#define _V_BYTE (_V/2)
#define _O_BYTE (_O/2)
#define _PUB_N_BYTE  (_PUB_N/2)
#define _PUB_M_BYTE  (_PUB_M/2)
#define _PUB_M1_BYTE (_PUB_M1/2)
#else
// GF256
#define _V_BYTE (_V)
#define _O_BYTE (_O)
#define _PUB_N_BYTE  (_PUB_N)
#define _PUB_M_BYTE  (_PUB_M)
#define _PUB_M1_BYTE (_PUB_M1)
#endif

/* length of seed for public key, in # bytes */
#define LEN_PKSEED 16

/* length of seed for secret key, in # bytes */
#define LEN_SKSEED 32

/* length of a signature */
#define TUOV_SIGNATUREBYTES (_PUB_N_BYTE)

#define N_TRIANGLE_TERMS(n_var) ((n_var)*((n_var)+1)/2)

/* length of components of public key, in # bytes */
#define _PK_BYTE        (_PUB_M_BYTE * N_TRIANGLE_TERMS(_PUB_N))
#define _PK_P1_BYTE     (_PUB_M_BYTE * N_TRIANGLE_TERMS(_V))
#define _PK_P2_BYTE     (_PUB_M_BYTE * _V)
#define _PK_P3_BYTE     (_PUB_M_BYTE * _V * (_PUB_M - 1))
#define _PK_P5_BYTE     (_PUB_M_BYTE)
#define _PK_P6_BYTE     (_PUB_M_BYTE * (_PUB_M - 1))
#define _PK_P9_BYTE     (_PUB_M_BYTE * N_TRIANGLE_TERMS(_PUB_M - 1))

/* length of components of secret key, in # bytes */
#define _SK_T1_BYTE     (_V_BYTE)
#ifdef _USE_GF16
#define _SK_T3_BYTE     (_PUB_M_BYTE)
#else
#define _SK_T3_BYTE     (_PUB_M_BYTE - 1)
#endif
#define _SK_T4_BYTE     (_V_BYTE * (_PUB_M - 1))
#define _SK_F1_BYTE     (_PK_P1_BYTE)
#define _SK_F2_BYTE     (_PK_P2_BYTE)
#define _SK_F3_BYTE     (_PK_P3_BYTE)
#define _SK_F5_BYTE     (_PUB_M1_BYTE)
#define _SK_F6_BYTE     (_PUB_M1_BYTE * (_PUB_M - 1))
#define _SK_S_BYTE      (_PUB_M1_BYTE * (_PUB_M - _PUB_M1))
#define _SK_F_BYTE      (_SK_F1_BYTE + _SK_F2_BYTE + _SK_F3_BYTE + _SK_F5_BYTE + _SK_F6_BYTE)
#define _SK_T_BYTE      (_SK_T1_BYTE + _SK_T3_BYTE + _SK_T4_BYTE)

/* Pij of public key pk is stored in pk->pk + Pij_BIAS */
#define P11_BIAS    0
#define P12_BIAS    (_PK_P1_BYTE >> 1)
#define P21_BIAS    (_PK_P1_BYTE)
#define P22_BIAS    (P21_BIAS + (_PK_P2_BYTE >> 1))
#define P31_BIAS    (P21_BIAS + _PK_P2_BYTE)
#define P32_BIAS    (P31_BIAS + (_PK_P3_BYTE >> 1))
#define P51_BIAS    (P31_BIAS + _PK_P3_BYTE)
#define P52_BIAS    (P51_BIAS + (_PK_P5_BYTE >> 1))
#define P61_BIAS    (P51_BIAS + _PK_P5_BYTE)
#define P62_BIAS    (P61_BIAS + (_PK_P6_BYTE >> 1))
#define P91_BIAS    (P61_BIAS + _PK_P6_BYTE)
#define P92_BIAS    (P91_BIAS + (_PK_P9_BYTE >> 1))

/* Ti of secret key sk is stored in sk->ST + Ti_BIAS */
#define T1_BIAS (_SK_S_BYTE)
#define T3_BIAS (T1_BIAS + _SK_T1_BYTE)
#define T4_BIAS (T3_BIAS + _SK_T3_BYTE)


#define TUOV_SK_UNCOMPRESSED_BYTES  (LEN_SKSEED + _SK_S_BYTE +  \
_SK_T_BYTE + _SK_F_BYTE)
#define TUOV_PK_UNCOMPRESSED_BYTES  (_PK_BYTE)
#define TUOV_SK_COMPRESSED_BYTES    (LEN_SKSEED + LEN_PKSEED)
#define TUOV_PK_COMPRESSED_BYTES    (LEN_PKSEED + _PK_P9_BYTE + \
(_PK_P5_BYTE >> 1) + (_PK_P6_BYTE >> 1))


#if defined _TUOV_CLASSIC
#define TUOV_ALGNAME (TUOV_PARAMNAME "-classic")
#define TUOV_SECRETKEYBYTES TUOV_SK_UNCOMPRESSED_BYTES
#define TUOV_PUBLICKEYBYTES TUOV_PK_UNCOMPRESSED_BYTES

#elif defined _TUOV_PKC
#define TUOV_ALGNAME (TUOV_PARAMNAME "-cpk")
#define TUOV_SECRETKEYBYTES TUOV_SK_UNCOMPRESSED_BYTES
#define TUOV_PUBLICKEYBYTES TUOV_PK_COMPRESSED_BYTES

#elif defined _TUOV_PKC_SKC
#define TUOV_ALGNAME (TUOV_PARAMNAME "-cpk-csk")
#define TUOV_SECRETKEYBYTES TUOV_SK_COMPRESSED_BYTES
#define TUOV_PUBLICKEYBYTES TUOV_PK_COMPRESSED_BYTES
#endif



#endif
