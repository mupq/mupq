// -----------------------------------------------------------------------------
// File Name   : field.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef FIELD_H
#define FIELD_H

#include "aimer.h"

#if       _AIMER_L == 1
  #include "field128.h"

  #define GF                      GF2_128                   // field element
  #define GF_set0                 GF2_128_set0              // set zero
  #define GF_copy                 GF2_128_copy              // copy field element
  #define GF_to_bytes             GF2_128_to_bytes          // convert field element to bytes
  #define GF_from_bytes           GF2_128_from_bytes        // convert bytes to field element
  #define GF_add                  GF2_128_add               // field addition
  #define GF_mul                  GF2_128_mul               // field multiplication
  #define GF_sqr                  GF2_128_sqr               // field squaring
  #define GF_transposed_matmul    GF2_128_transposed_matmul // transposed matrix multiplication

#elif     _AIMER_L == 3
  #include "field192.h"

  #define GF                      GF2_192                   // field element
  #define GF_set0                 GF2_192_set0              // set zero
  #define GF_copy                 GF2_192_copy              // copy field element
  #define GF_to_bytes             GF2_192_to_bytes          // convert field element to bytes
  #define GF_from_bytes           GF2_192_from_bytes        // convert bytes to field element
  #define GF_add                  GF2_192_add               // field addition
  #define GF_mul                  GF2_192_mul               // field multiplication
  #define GF_sqr                  GF2_192_sqr               // field squaring
  #define GF_transposed_matmul    GF2_192_transposed_matmul // transposed matrix multiplication

#elif     _AIMER_L == 5
  #include "field256.h"
  
  #define GF                      GF2_256                   // field element
  #define GF_set0                 GF2_256_set0              // set zero
  #define GF_copy                 GF2_256_copy              // copy field element
  #define GF_to_bytes             GF2_256_to_bytes          // convert field element to bytes
  #define GF_from_bytes           GF2_256_from_bytes        // convert bytes to field element
  #define GF_add                  GF2_256_add               // field addition
  #define GF_mul                  GF2_256_mul               // field multiplication
  #define GF_sqr                  GF2_256_sqr               // field squaring
  #define GF_transposed_matmul    GF2_256_transposed_matmul // transposed matrix multiplication

#endif

#endif // FIELD_H
