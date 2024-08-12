// SPDX-License-Identifier: MIT

#ifndef SIGN_H
#define SIGN_H

#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef struct tape_t
{
  GF pt_share;
  GF t_shares[AIMER_L];
  GF a_share;
  GF c_share;
} tape_t;

typedef struct mult_chk_t
{
  GF x_shares[AIMER_L + 1];
  GF z_shares[AIMER_L + 1];
} mult_chk_t;

typedef struct proof_t
{
  uint8_t reveal_path[AIMER_LOGN][AIMER_SEED_SIZE];
  uint8_t missing_commitment[AIMER_COMMIT_SIZE];
  uint8_t delta_pt_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t delta_ts_bytes[AIMER_L][AIM2_NUM_BYTES_FIELD];
  uint8_t delta_c_bytes[AIM2_NUM_BYTES_FIELD];
  uint8_t missing_alpha_share_bytes[AIM2_NUM_BYTES_FIELD];
} proof_t;

typedef struct signature_t
{
  uint8_t salt[AIMER_SALT_SIZE];
  uint8_t h_1[AIMER_COMMIT_SIZE];
  uint8_t h_2[AIMER_COMMIT_SIZE];
  proof_t proofs[AIMER_T];
} signature_t;

#define aim2_mpc AIMER_NAMESPACE(aim2_mpc)
void aim2_mpc(mult_chk_t *mult_chk,
              const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD],
              const GF ct_GF);

#define commit_and_expand_tape AIMER_NAMESPACE(commit_and_expand_tape)
void commit_and_expand_tape(tape_t *tape, uint8_t *commit,
                            const hash_instance *ctx_precom,
                            const uint8_t seed[AIMER_SEED_SIZE],
                            size_t rep, size_t party);

#define run_phase_1 AIMER_NAMESPACE(run_phase_1)
void run_phase_1(signature_t *sign,
                 const GF pt_GF, const GF sbox_outputs[AIMER_L],
                 const uint8_t mu[AIMER_COMMIT_SIZE],
                 const uint8_t root_seeds[AIMER_T][AIMER_SEED_SIZE]);

#define run_phase_1_to_3 AIMER_NAMESPACE(run_phase_1_to_3)
void run_phase_1_to_3(signature_t *sign,
                      const uint8_t root_seeds[AIMER_T][AIMER_SEED_SIZE],
                      const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD],
                      const GF vector_b, const GF ct_GF);

#define run_phase_1_to_5 AIMER_NAMESPACE(run_phase_1_to_5)
void run_phase_1_to_5(signature_t *sign,
                      const uint8_t root_seeds[AIMER_T][AIMER_SEED_SIZE],
                      const GF matrix_A[AIMER_L][AIM2_NUM_BITS_FIELD],
                      const GF vector_b);

#endif // SIGN_H
