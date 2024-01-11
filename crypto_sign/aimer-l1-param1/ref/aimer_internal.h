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
  uint8_t* tape;
  size_t  random_tape_size;
} random_tape_t;

typedef struct
{
  reveal_list_t reveal_list;
  uint8_t*      missing_commitment;
  uint8_t*      pt_delta;
  GF*           z_delta;
  GF            c_delta;
  GF            missing_alpha_share;
} proof_t;

typedef struct
{
  uint8_t* salt;
  uint8_t* h_1;
  uint8_t* h_2;
  proof_t* proofs;
} signature_t;

void allocate_proof(const aimer_instance_t* instance, proof_t* proof);

void allocate_signature(const aimer_instance_t* instance, signature_t* sig);

void free_proof(proof_t* proof);

void free_signature(const aimer_instance_t* instance, signature_t* sig);

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
                GF** epsilons);

void h_2_commitment(const aimer_instance_t* instance, const uint8_t* salt,
                    const uint8_t* h_1, const GF* repetition_alpha_shares,
                    const GF** v_shares, uint8_t* h_2);

void h_2_expand(const aimer_instance_t* instance, const uint8_t* h_2,
                uint16_t* opened_parties);

#endif // AIMER_INTERNAL_H
