/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

lowmc_partial_round_t const* round = LOWMC_INSTANCE.rounds;
  mzd_local_t nl_part[shares][(LOWMC_R * 32 + 255) / 256];
  MPC_LOOP_CONST_C(XOR, x, x, LOWMC_INSTANCE.precomputed_constant_linear, shares, ch);
  MPC_LOOP_CONST(MUL_MC, nl_part, in_out_shares->s,
                 LOWMC_INSTANCE.precomputed_non_linear_part_matrix, shares);
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, LOWMC_INSTANCE.precomputed_constant_non_linear, shares, ch);
  for (unsigned i = 0; i < (LOWMC_R-1); ++i, ++views, ++round, ++rvec) {
    SBOX_uint64(sbox, y, x, views, rvec, LOWMC_N, shares, shares);
    for (unsigned int k = 0; k < shares; ++k) {
      const word nl = CONST_BLOCK(nl_part[k], i >> 3)->w64[(i & 0x7) >> 1];
      BLOCK(y[k], 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
    }
    MPC_LOOP_CONST(MUL_Z, x, y, round->z_matrix, shares);

    for(unsigned int k = 0; k < shares; ++k) {
      SHUFFLE(y[k], round->r_mask);
    }

    MPC_LOOP_CONST(ADDMUL_R, x, y, round->r_matrix, shares);
    for(unsigned int k = 0; k < shares; ++k) {
      BLOCK(y[k], 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] &= WORD_C(0x00000003FFFFFFFF); //clear nl part
    }
    MPC_LOOP_SHARED(XOR, x, x, y, shares);
  }
  unsigned i = (LOWMC_R-1);
  SBOX_uint64(sbox, y, x, views, rvec, LOWMC_N, shares, shares);

  for (unsigned int k = 0; k < shares; ++k) {
    const word nl = CONST_BLOCK(nl_part[k], i >> 3)->w64[(i & 0x7) >> 1];
    BLOCK(y[k], 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
      (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
  }
  MPC_LOOP_CONST(MUL, x, y, LOWMC_INSTANCE.zr_matrix, shares);

// vim: ft=c
