/*
 * Copyright (c) The mldsa-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLD_INTEGRATION_PQM4_CONFIG_H
#define MLD_INTEGRATION_PQM4_CONFIG_H

#define MLD_CONFIG_PARAMETER_SET 87
#define MLD_CONFIG_NAMESPACE_PREFIX mldsa
#define MLD_CONFIG_FIPS202_CUSTOM_HEADER "fips202_glue.h"
#define MLD_CONFIG_INTERNAL_API_QUALIFIER static
#define MLD_CONFIG_SERIAL_FIPS202_ONLY
//#define MLD_CONFIG_REDUCE_RAM

#define MLD_CONFIG_CUSTOM_RANDOMBYTES
#if !defined(__ASSEMBLER__)
#include "randombytes.h"
#include "mldsa/src/sys.h"

static MLD_INLINE void mld_randombytes(uint8_t *ptr, size_t len)
{
  randombytes(ptr, len);
}
#endif /* !__ASSEMBLER__ */

#endif