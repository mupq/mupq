// -----------------------------------------------------------------------------
// File Name   : aim128.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aim128.h"
#include "hash.h"

#include "aim_common.inc"
#include "portable_endian.h"

void mersenne_exp_3(const GF in, GF out);
void mersenne_exp_5(const GF in, GF out);
void mersenne_exp_27(const GF in, GF out);

void mersenne_exp_3(const GF in, GF out)
{
  GF t1 = {0,};

  // t1 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t1);

  // out = a ^ (2^3 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, out);
}

void mersenne_exp_5(const GF in, GF out)
{
  GF t1 = {0,};
  GF t2 = {0,};

  // t2 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t2);

  // t1 = a ^ (2^4 - 1)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, t1);

  // out = a ^ (2^5 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, out);
}

void mersenne_exp_27(const GF in, GF out)
{
  int i;
  GF t1 = {0,};
  GF t2 = {0,};
  GF t3 = {0,};

  // t1 = a ^ (2^2 - 1)
  GF_sqr(in, t1);
  GF_mul(t1, in, t1);

  // t2 = a ^ (2^3 - 1)
  GF_sqr(t1, t1);
  GF_mul(t1, in, t2);

  // t3 = a ^ (2^6 - 1)
  GF_sqr(t2, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, t3);

  // t3 = a ^ (2^12 - 1)
  GF_sqr(t3, t1);
  for (i = 1; i < 6; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t3);

  // t3 = a ^ (2^24 - 1)
  GF_sqr(t3, t1);
  for (i = 1; i < 12; i++)
  {
    GF_sqr(t1, t1);
  }
  GF_mul(t1, t3, t3);

  // out = a ^ (2^27 - 1)
  GF_sqr(t3, t1);
  GF_sqr(t1, t1);
  GF_sqr(t1, t1);
  GF_mul(t1, t2, out);
}

void generate_matrices_L_and_U(const uint8_t* iv, GF matrix_A[2*AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS], GF vector_b)
{
  hash_instance ctx;
  size_t squeeze_len;
  size_t block_index;

  // Initialize Hash
  hash_init(&ctx, (BLOCK_SIZE << 1));
  hash_update(&ctx, iv, BLOCK_SIZE);
  hash_final(&ctx);

  uint8_t out[TAPE_LEN];
  uint8_t* out_ptr;

  for (size_t num = 0; num < NUM_INPUT_SBOX; num++)
  {
    hash_squeeze(&ctx, out, TAPE_LEN);

    // Lower
    squeeze_len = 8;
    block_index = 0;

    out_ptr = out;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      memcpy((uint8_t *)(&tmp) + (8 - squeeze_len), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      matrix_A[2 * num][row][block_index] =
        (le64toh(tmp) & and_mask_lower[row & 0x3f]) | or_mask[row & 0x3f];
      for (size_t i = block_index + 1; i < NUMWORDS_FIELD; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        matrix_A[2 * num][row][i] = le64toh(tmp);
      }
      if ((row & 7) == 6)
      {
        squeeze_len--;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 8;
      }
    }

    // Upper
    squeeze_len = 0;
    block_index = 0;
    for (size_t row = 0; row < NUMBITS_FIELD; row++)
    {
      uint64_t tmp = 0;
      for (size_t i = 0; i < block_index; i++)
      {
        memcpy((uint8_t *)(&tmp), out_ptr, 8);
        out_ptr += 8;
        matrix_A[2 * num + 1][row][i] = le64toh(tmp);
      }
      memcpy((uint8_t *)(&tmp), out_ptr, squeeze_len);
      out_ptr += squeeze_len;
      matrix_A[2 * num + 1][row][block_index] =
        (le64toh(tmp) & and_mask_upper[row & 0x3f]) | or_mask[row & 0x3f];

      if ((row & 7) == 0)
      {
        squeeze_len++;
      }
      if ((row & 0x3f) == 0x3f)
      {
        block_index++;
        squeeze_len = 0;
      }
    }
  }

  // generate vector B
  hash_squeeze(&ctx, (uint8_t*)vector_b, BLOCK_SIZE);
}

void generate_matrix_LU(const uint8_t* iv, GF matrix_A[AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS], GF vector_b)
{
  GF temp_matrix[2*AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS] = {0};

  generate_matrices_L_and_U(iv, temp_matrix, vector_b);

  for (size_t i = 0; i < NUM_INPUT_SBOX; i++)
  {
    for (size_t j = 0; j < NUMBITS_FIELD; j++)
    {
      GF_transposed_matmul(temp_matrix[2 * i + 1][j], temp_matrix[2 * i],
                           matrix_A[i][j]);
    }
  }
}

void compute_sbox_outputs(const uint8_t* pt, GF* sbox_outputs)
{
  GF pt_GF = {0,};

  GF_from_bytes(pt, pt_GF);

  // Mersenne layer
  mersenne_exp_3(pt_GF, sbox_outputs[0]);
  mersenne_exp_27(pt_GF, sbox_outputs[1]);
}

void aim(const uint8_t* pt, const uint8_t* iv, uint8_t* ct)
{
  GF state[NUM_INPUT_SBOX] = {0,};
  GF pt_GF = {0,};
  GF vector_b = {0,};
  GF matrix_A[2*AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS] = {0};

  GF_from_bytes(pt, pt_GF);

  // Mersenne layer
  mersenne_exp_3(pt_GF, state[0]);
  mersenne_exp_27(pt_GF, state[1]);

  // generate linear layer
  generate_matrices_L_and_U(iv, matrix_A, vector_b);

  // linear layer
  GF_transposed_matmul(state[0], matrix_A[1], state[0]);
  GF_transposed_matmul(state[0], matrix_A[0], state[0]);

  GF_transposed_matmul(state[1], matrix_A[3], state[1]);
  GF_transposed_matmul(state[1], matrix_A[2], state[1]);

  GF_add(state[0], state[1], state[0]);
  GF_add(state[0], vector_b, state[0]);

  // Mersenne layer
  mersenne_exp_5(state[0], state[0]);

  // feedback
  GF_add(state[0], pt_GF, state[0]);

  // ct conversion
  GF_to_bytes(state[0], ct);
}

void aim_mpc(const uint8_t* pt, const GF matrix_A[AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS],
             const GF vector_b, const uint8_t* ct,
             const size_t num_parties, GF *z_shares, GF *x_shares)
{
  GF state[NUM_INPUT_SBOX] = {0,};
  GF ct_GF = {0,};
  GF temp = {0,};
  size_t party_index;

  GF_from_bytes(ct, ct_GF);

  for (size_t party = 0; party < num_parties; party++)
  {
    party_index = party * (NUM_INPUT_SBOX + 1);

    GF_from_bytes(pt + party * BLOCK_SIZE, x_shares[party_index]);

    // linear layer
    GF_transposed_matmul(z_shares[party_index], matrix_A[0], state[0]);
    GF_transposed_matmul(z_shares[party_index + 1], matrix_A[1], state[1]);
    if (party == 0)
    {
      GF_add(state[1], vector_b, state[1]);
    }
    GF_add(state[0], state[1], x_shares[party_index + 2]);

    // optimization B.1.2 in eprint 2022/588
    GF_copy(z_shares[party_index], state[0]);
    GF_copy(z_shares[party_index + 1], state[1]);

    GF_sqr(x_shares[party_index], temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, z_shares[party_index]);

    GF_transposed_matmul(x_shares[party_index], e2_power_matrix,
                         z_shares[party_index + 1]);

    GF_copy(state[0], x_shares[party_index]);
    GF_copy(state[1], x_shares[party_index + 1]);

    GF_sqr(x_shares[party_index + 2], temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);
    GF_sqr(temp, temp);

    GF_mul(ct_GF, x_shares[party_index + 2], state[0]);
    GF_add(temp, state[0], z_shares[party_index + 2]);
  }
}
