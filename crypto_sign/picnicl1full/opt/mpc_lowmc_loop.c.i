/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

lowmc_round_t const* round = LOWMC_INSTANCE.rounds;
for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round, ++rvec) {
  SBOX(sbox, y, x, views, rvec, LOWMC_N, shares, shares);
  MPC_LOOP_CONST(MUL, x, y, round->l_matrix, shares);
  MPC_LOOP_CONST_C(XOR, x, x, round->constant, shares, ch);
  MPC_LOOP_CONST(ADDMUL, x, in_out_shares->s, round->k_matrix, shares);
}

// vim: ft=c
