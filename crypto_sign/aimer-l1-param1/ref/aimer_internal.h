// -----------------------------------------------------------------------------
// File Name   : aimer_internal.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIMER_INTERNAL_H
#define AIMER_INTERNAL_H

#include "hash.h"
#include "tree.h"
#include "field.h"
#include "aimer_instances.h"

typedef struct random_tape_t 
{
  uint8_t tape[AIMER_T * AIMER_N * (AIMER_BLOCK_SIZE +
                                  AIMER_NUM_INPUT_SBOXES * AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE)];
  size_t  random_tape_size;
} random_tape_t;

typedef struct
{
  reveal_list_t reveal_list;
  uint8_t      missing_commitment[AIMER_DIGEST_SIZE];
  uint8_t      pt_delta[AIMER_BLOCK_SIZE];
  GF           z_delta[AIMER_NUM_INPUT_SBOXES];
  GF            c_delta;
  GF            missing_alpha_share;
} proof_t;

typedef struct
{
  uint8_t salt[AIMER_SALT_SIZE];
  uint8_t h_1[AIMER_DIGEST_SIZE];
  uint8_t h_2[AIMER_DIGEST_SIZE];
  proof_t proofs[AIMER_T];
} signature_t;


void commit_to_seed_and_expand_tape(const aimer_instance_t* instance, const uint8_t* seed,
                                    const uint8_t* salt, size_t repetition, size_t party,
                                    uint8_t* commitment,
                                    random_tape_t* tapes);

void commit_to_seed_and_expand_tape_x4(const aimer_instance_t* instance,
                                       const uint8_t* seed0, const uint8_t* seed1,
                                       const uint8_t* seed2, const uint8_t* seed3,
                                       const uint8_t* salt,
                                       size_t repetition, size_t party,
                                       uint8_t* commitments,
                                       random_tape_t* tapes);

void h_1_commitment(const aimer_instance_t* instance,
                    const signature_t* sig,
                    const aimer_publickey_t* public_key,
                    const uint8_t* message, size_t message_len,
                    const uint8_t* party_seed_commitments, uint8_t* h_1);

void h_1_expand(const aimer_instance_t* instance, const uint8_t* h_1,
                GF epsilons[AIMER_T][AIMER_NUM_INPUT_SBOXES+1]);

void h_2_commitment(const aimer_instance_t* instance, const uint8_t* salt,
                    const uint8_t* h_1, const GF* repetition_alpha_shares,
                    const GF v_shares[AIMER_T][AIMER_N], uint8_t* h_2);

void h_2_expand(const aimer_instance_t* instance, const uint8_t* h_2,
                uint16_t* opened_parties);

#endif // AIMER_INTERNAL_H
