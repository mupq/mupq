/*
 * Copyright (c) The mlkem-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLK_INTEGRATION_PQM4_CONFIG_H
#define MLK_INTEGRATION_PQM4_CONFIG_H

#define MLK_CONFIG_NAMESPACE_PREFIX mlkem_native
#define MLK_CONFIG_PARAMETER_SET 512

#define MLK_CONFIG_FIPS202_CUSTOM_HEADER "fips202_glue.h"
#define MLK_CONFIG_FIPS202X4_CUSTOM_HEADER "fips202x4_glue.h"

#define MLK_CONFIG_CUSTOM_RANDOMBYTES
#if !defined(__ASSEMBLER__)
#include "randombytes.h"
#include "sys.h"

static MLK_INLINE void mlk_randombytes(uint8_t *ptr, size_t len)
{
  randombytes(ptr, len);
}
#endif /* !__ASSEMBLER__ */

#endif
