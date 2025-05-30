/*
 * Copyright (c) The mlkem-native project authors
 * SPDX-License-Identifier: Apache-2.0 OR ISC OR MIT
 */
#ifndef MLK_INTEGRATION_PQM4_FIPS202_GLUE_H
#define MLK_INTEGRATION_PQM4_FIPS202_GLUE_H

/* Include pqm4's own FIPS202 header */
#include "fips202.h"

#define mlk_shake128ctx shake128ctx
#define mlk_shake128_absorb_once shake128_absorb
#define mlk_shake128_squeezeblocks shake128_squeezeblocks
#define mlk_shake128_init(S) do { } while(0) /* no-op */
#define mlk_shake128_release(S) do { } while(0) /* no-op */
#define mlk_shake256 shake256
#define mlk_sha3_256 sha3_256
#define mlk_sha3_512 sha3_512


#endif
