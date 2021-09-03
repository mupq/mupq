/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */


#include "randomness.h"
#include "macros.h"

// randombytes from the PQClean
extern void randombytes(uint8_t* x, size_t xlen);

static int rand_bytes(uint8_t* dst, size_t len) {
  randombytes(dst, len);
  return 0;
}

int rand_bits(uint8_t* dst, size_t num_bits) {
  const size_t num_bytes = (num_bits + 7) / 8;
  const size_t num_extra_bits = num_bits % 8;

  if (rand_bytes(dst, num_bytes)) {
    return -1;
  }

  if (num_extra_bits) {
    dst[num_bytes - 1] &= UINT8_C(0xff) << (8 - num_extra_bits);
  }

  return 0;
}
