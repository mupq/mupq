/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(LOWMC_INSTANCE)
#define N_LOWMC CONCAT(lowmc, CONCAT(IMPL, LOWMC_INSTANCE))
#define SBOX_FUNC CONCAT(sbox, CONCAT(IMPL, LOWMC_INSTANCE))
#if defined(LOWMC_PARTIAL)
#define SBOX(x) sbox_layer_10_uint64(&BLOCK(x, 0)->w64[(LOWMC_N / (sizeof(word) * 8)) - 1])
#include "lowmc_impl_partial.c.i"
#else
#define SBOX(x) SBOX_FUNC(BLOCK(x, 0))
#include "lowmc_impl.c.i"
#endif


#undef N_LOWMC
#undef SBOX
#undef SBOX_FUNC
#endif

// vim: ft=c
