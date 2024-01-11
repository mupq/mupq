// -----------------------------------------------------------------------------
// File Name   : field128.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef FIELD128_H
#define FIELD128_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NUMWORDS_FIELD      2  // number of 64-bit words in field element
#define NUMBITS_WORD       64  // number of bits in word
#define NUMBITS_FIELD     128  // number of bits in field element

typedef uint64_t GF2_128[2];

void GF2_128_set0(GF2_128 a);
void GF2_128_copy(const GF2_128 in, GF2_128 out);
void GF2_128_to_bytes(const GF2_128 in, uint8_t* out);
void GF2_128_from_bytes(const uint8_t* in, GF2_128 out);

void GF2_128_add(const GF2_128 a, const GF2_128 b, GF2_128 c);
void GF2_128_mul(const GF2_128 a, const GF2_128 b, GF2_128 c);
void GF2_128_sqr(const GF2_128 a, GF2_128 c);
void GF2_128_transposed_matmul(const GF2_128 a, const GF2_128* b, GF2_128 c);

#endif // FIELD128_H
