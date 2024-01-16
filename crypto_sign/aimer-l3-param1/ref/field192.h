// -----------------------------------------------------------------------------
// File Name   : field192.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef FIELD192_H
#define FIELD192_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define NUMWORDS_FIELD      3  // number of 64-bit words in field element
#define NUMBITS_WORD       64  // number of bits in word
#define NUMBITS_FIELD     192  // number of bits in field element

typedef uint64_t GF2_192[3];

void GF2_192_set0(GF2_192 a);
void GF2_192_copy(const GF2_192 in, GF2_192 out);
void GF2_192_to_bytes(const GF2_192 in, uint8_t* out);
void GF2_192_from_bytes(const uint8_t* in, GF2_192 out);

void GF2_192_add(const GF2_192 a, const GF2_192 b, GF2_192 c);
void GF2_192_mul(const GF2_192 a, const GF2_192 b, GF2_192 c);
void GF2_192_sqr(const GF2_192 a, GF2_192 c);
void GF2_192_transposed_matmul(const GF2_192 a, const GF2_192* b, GF2_192 c);

#endif // FIELD192_H
