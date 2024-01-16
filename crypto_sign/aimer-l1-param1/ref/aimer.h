// -----------------------------------------------------------------------------
// File Name   : aimer.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIMER_H
#define AIMER_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include "api.h"


#if       _AIMER_L == 1
  #define AIMER_PUBLICKEY_SIZE        32    // size of AIMer public key
  #define AIMER_PRIVATEKEY_SIZE       16    // size of AIMer private key
  #define AIMER_MAX_SIGNATURE_SIZE  5904    // maximum size of AIMer signature

#elif     _AIMER_L == 3
  #define AIMER_PUBLICKEY_SIZE        48    // size of AIMer public key
  #define AIMER_PRIVATEKEY_SIZE       24    // size of AIMer private key
  #define AIMER_MAX_SIGNATURE_SIZE 13080    // maximum size of AIMer signature

#elif     _AIMER_L == 5
  #define AIMER_PUBLICKEY_SIZE        64    // size of AIMer public key
  #define AIMER_PRIVATEKEY_SIZE       32    // size of AIMer private key
  #define AIMER_MAX_SIGNATURE_SIZE 25152    // maximum size of AIMer signature

#else
  #error  "does not support"
#endif

#if _AIMER_L == 1
  #define AIMER_SALT_SIZE 32
  #define AIMER_DIGEST_SIZE 32
  #define AIMER_SEED_SIZE 16
  #define AIMER_FIELD_SIZE 16
  #define AIMER_BLOCK_SIZE 16
  #define AIMER_NUM_INPUT_SBOXES 2
  #define AIMER_NUM_BITS 128
#elif _AIMER_L == 3
  #define AIMER_SALT_SIZE 48
  #define AIMER_DIGEST_SIZE 48
  #define AIMER_SEED_SIZE 24
  #define AIMER_FIELD_SIZE 24
  #define AIMER_BLOCK_SIZE 24
  #define AIMER_NUM_INPUT_SBOXES 2
  #define AIMER_NUM_BITS 192
#elif _AIMER_L == 5
  #define AIMER_SALT_SIZE 64
  #define AIMER_DIGEST_SIZE 64
  #define AIMER_SEED_SIZE 32
  #define AIMER_FIELD_SIZE 32
  #define AIMER_BLOCK_SIZE 32
  #define AIMER_NUM_INPUT_SBOXES 3
  #define AIMER_NUM_BITS 256
#else
  #error  "does not support"
#endif

  //{AIM_128_PARAMS, 32, 32, 16, 16, 33,   16, AIMER_L1_PARAM1},
  //{AIM_128_PARAMS, 32, 32, 16, 16, 23,   57, AIMER_L1_PARAM2},
  //{AIM_128_PARAMS, 32, 32, 16, 16, 17,  256, AIMER_L1_PARAM3},
  //{AIM_128_PARAMS, 32, 32, 16, 16, 13, 1615, AIMER_L1_PARAM4},
  //{AIM_192_PARAMS, 48, 48, 24, 24, 49,   16, AIMER_L3_PARAM1},
  //{AIM_192_PARAMS, 48, 48, 24, 24, 33,   64, AIMER_L3_PARAM2},
  //{AIM_192_PARAMS, 48, 48, 24, 24, 25,  256, AIMER_L3_PARAM3},
  //{AIM_192_PARAMS, 48, 48, 24, 24, 19, 1621, AIMER_L3_PARAM4},
  //{AIM_256_PARAMS, 64, 64, 32, 32, 65,   16, AIMER_L5_PARAM1},
  //{AIM_256_PARAMS, 64, 64, 32, 32, 44,   62, AIMER_L5_PARAM2},
  //{AIM_256_PARAMS, 64, 64, 32, 32, 33,  256, AIMER_L5_PARAM3},
  //{AIM_256_PARAMS, 64, 64, 32, 32, 25, 1623, AIMER_L5_PARAM4}
#if _AIMER_L == 1 && AIMER_PARAM == 1
  #define AIMER_T 33
  #define AIMER_N 16
  #define AIMER_LOGN 4
  #define AIMER_INSTANCE AIMER_L1_PARAM1
#elif _AIMER_L == 1 && AIMER_PARAM == 2
  #define AIMER_T 23
  #define AIMER_N 57
  #define AIMER_LOGN 6
  #define AIMER_INSTANCE AIMER_L1_PARAM2
#elif _AIMER_L == 1 && AIMER_PARAM == 3
  #define AIMER_T 17
  #define AIMER_N 256
  #define AIMER_LOGN 8
  #define AIMER_INSTANCE AIMER_L1_PARAM3
#elif _AIMER_L == 1 && AIMER_PARAM == 4
  #define AIMER_T 13
  #define AIMER_N 1615
  #define AIMER_LOGN 11
  #define AIMER_INSTANCE AIMER_L1_PARAM4
#elif _AIMER_L == 3 && AIMER_PARAM == 1
  #define AIMER_T 49
  #define AIMER_N 16
  #define AIMER_LOGN 4
  #define AIMER_INSTANCE AIMER_L3_PARAM1
#elif _AIMER_L == 3 && AIMER_PARAM == 2
  #define AIMER_T 33
  #define AIMER_N 64
  #define AIMER_LOGN 6
  #define AIMER_INSTANCE AIMER_L3_PARAM2
#elif _AIMER_L == 3 && AIMER_PARAM == 3
  #define AIMER_T 25
  #define AIMER_N 256
  #define AIMER_LOGN 8
  #define AIMER_INSTANCE AIMER_L3_PARAM3
#elif _AIMER_L == 3 && AIMER_PARAM == 4
  #define AIMER_T 19
  #define AIMER_N 1621
  #define AIMER_LOGN 11
  #define AIMER_INSTANCE AIMER_L3_PARAM4
#elif _AIMER_L == 5 && AIMER_PARAM == 1
  #define AIMER_T 65
  #define AIMER_N 16
  #define AIMER_LOGN 4
  #define AIMER_INSTANCE AIMER_L5_PARAM1
#elif _AIMER_L == 5 && AIMER_PARAM == 2
  #define AIMER_T 44
  #define AIMER_N 62
  #define AIMER_LOGN 6
  #define AIMER_INSTANCE AIMER_L5_PARAM2
#elif _AIMER_L == 5 && AIMER_PARAM == 3
  #define AIMER_T 33
  #define AIMER_N 256
  #define AIMER_LOGN 8
  #define AIMER_INSTANCE AIMER_L5_PARAM3
#elif _AIMER_L == 5 && AIMER_PARAM == 4
  #define AIMER_T 25
  #define AIMER_N 1623
  #define AIMER_LOGN 11
  #define AIMER_INSTANCE AIMER_L5_PARAM4
#else
  #error  "does not support"
#endif


#define AIMER_TREE_NUM_NODES (((1 << (AIMER_LOGN+1)) - 1) - ((1 << (AIMER_LOGN)) - AIMER_N))


// Parameter set names
typedef enum
{
  PARAMETER_SET_INVALID   =  0,
  AIMER_L1_PARAM1         =  1,
  AIMER_L1_PARAM2         =  2,
  AIMER_L1_PARAM3         =  3,
  AIMER_L1_PARAM4         =  4,
  AIMER_L3_PARAM1         =  5,
  AIMER_L3_PARAM2         =  6,
  AIMER_L3_PARAM3         =  7,
  AIMER_L3_PARAM4         =  8,
  AIMER_L5_PARAM1         =  9,
  AIMER_L5_PARAM2         = 10,
  AIMER_L5_PARAM3         = 11,
  AIMER_L5_PARAM4         = 12,
  PARAMETER_SET_MAX_INDEX = 13
} aimer_params_t;


// Public key
typedef struct
{
  uint8_t data[AIMER_PUBLICKEY_SIZE];
  aimer_params_t params;
} aimer_publickey_t;

// Private key
typedef struct
{
  uint8_t data[AIMER_PRIVATEKEY_SIZE];
} aimer_privatekey_t;

// Signature API

// Key generation function.
// Generates a public and private key pair, for the specified parameter set.
//
// @param[in]  parameters The parameter set to use when generating a key.
// @param[out] pk         The new public key.
// @param[out] sk         The new private key.
//
// @return Returns 0 for success, or a nonzero value indicating an error.

int aimer_keygen(aimer_params_t param,
                 aimer_publickey_t* public_key,
                 aimer_privatekey_t* private_key);

int aimer_sign(const aimer_publickey_t*  public_key,
               const aimer_privatekey_t* private_key,
               const uint8_t* message, const size_t message_len,
               uint8_t* signature, size_t* signature_len);

int aimer_verify(const aimer_publickey_t* public_key,
                 const uint8_t* signature, const size_t signature_len,
                 const uint8_t* message, const size_t message_len);

#endif // AIMER_H
