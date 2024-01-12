// -----------------------------------------------------------------------------
// File Name   : field256.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef FIELD256_H
#define FIELD256_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NUMWORDS_FIELD      4  // number of 64-bit words in field element
#define NUMBITS_WORD       64  // number of bits in word
#define NUMBITS_FIELD     256  // number of bits in field element

typedef uint64_t GF2_256[4];

void GF2_256_set0(GF2_256 a);
void GF2_256_copy(const GF2_256 in, GF2_256 out);
void GF2_256_to_bytes(const GF2_256 in, uint8_t* out);
void GF2_256_from_bytes(const uint8_t* in, GF2_256 out);

void GF2_256_add(const GF2_256 a, const GF2_256 b, GF2_256 c);
void GF2_256_mul(const GF2_256 a, const GF2_256 b, GF2_256 c);
void GF2_256_sqr(const GF2_256 a, GF2_256 c);
void GF2_256_transposed_matmul(const GF2_256 a, const GF2_256* b, GF2_256 c);

#endif // FIELD256_H
