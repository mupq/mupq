/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "randomness.h"
#include "macros.h"

// randombytes from the MUPQ
#include "randombytes.h"

int rand_bytes(uint8_t* dst, size_t len) {
  return randombytes(dst, len);
}
