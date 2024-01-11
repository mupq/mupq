// -----------------------------------------------------------------------------
// File Name   : aimer_internal.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aimer_internal.h"
#include "portable_endian.h"
#include <stdlib.h>

void commit_to_seed_and_expand_tape(const aimer_instance_t* instance,
                                    const uint8_t* seed, const uint8_t* salt,
                                    size_t repetition, size_t party,
                                    uint8_t* commitment,
                                    random_tape_t* tapes)
{
  hash_instance ctx;

  hash_init_prefix(&ctx, instance->digest_size, HASH_PREFIX_4);
  hash_update(&ctx, salt, instance->salt_size);
  hash_update_uint16(&ctx, (uint16_t)repetition);
  hash_update_uint16(&ctx, (uint16_t)party);
  hash_update(&ctx, seed, instance->seed_size);
  hash_final(&ctx);

  hash_squeeze(&ctx, 
               commitment + (repetition * instance->num_MPC_parties + party)
               * instance->digest_size, instance->digest_size);

  hash_squeeze(&ctx,
              tapes->tape +
              (repetition * instance->num_MPC_parties + party) *
              tapes->random_tape_size,
              tapes->random_tape_size);
}

void commit_to_seed_and_expand_tape_x4(const aimer_instance_t* instance,
                                       const uint8_t* seed0, const uint8_t* seed1,
                                       const uint8_t* seed2, const uint8_t* seed3,
                                       const uint8_t* salt,
                                       size_t repetition, size_t party,
                                       uint8_t* commitments,
                                       random_tape_t* tapes)
{
  hash_instance_x4 ctx;

  hash_init_prefix_x4(&ctx, instance->digest_size, HASH_PREFIX_4);
  hash_update_x4_1(&ctx, salt, instance->salt_size);
  hash_update_x4_uint16(&ctx, (uint16_t)repetition);

  uint16_t parties[4] = {party, party + 1, party + 2, party + 3};
  hash_update_x4_uint16s(&ctx, parties);
  hash_update_x4_4(&ctx, seed0, seed1, seed2, seed3, instance->seed_size);
  hash_final_x4(&ctx);

  size_t index = (repetition * instance->num_MPC_parties + party);
  hash_squeeze_x4_4(&ctx, commitments + index * instance->digest_size,
                    commitments + (index + 1) * instance->digest_size,
                    commitments + (index + 2) * instance->digest_size,
                    commitments + (index + 3) * instance->digest_size,
                    instance->digest_size);

  hash_squeeze_x4_4(&ctx, tapes->tape + index * tapes->random_tape_size,
                    tapes->tape + (index + 1) * tapes->random_tape_size,
                    tapes->tape + (index + 2) * tapes->random_tape_size,
                    tapes->tape + (index + 3) * tapes->random_tape_size,
                    tapes->random_tape_size);
}

void h_1_commitment(const aimer_instance_t* instance,
                    const signature_t* sig,
                    const aimer_publickey_t* public_key,
                    const uint8_t* message, size_t message_len,
                    const uint8_t* party_seed_commitments, uint8_t* h_1)
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init_prefix(&ctx, instance->digest_size, HASH_PREFIX_1);
  hash_update(&ctx, message, message_len);
  hash_update(&ctx, public_key->data, sizeof(public_key->data));
  hash_update(&ctx, sig->salt, instance->salt_size);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];
    for (size_t party = 0; party < instance->num_MPC_parties; party++)
    {
      hash_update(&ctx, party_seed_commitments +
                  ((repetition * instance->num_MPC_parties + party) *
                  instance->digest_size), instance->digest_size);
    }

    hash_update(&ctx, proof->pt_delta, instance->aim_params.block_size);

    for (size_t ell = 0; ell < instance->aim_params.num_input_sboxes; ell++)
    {
      GF_to_bytes(proof->z_delta[ell], buffer);
      hash_update(&ctx, buffer, instance->field_size);
    }

    GF_to_bytes(proof->c_delta, buffer);
    hash_update(&ctx, buffer, instance->field_size);
  }
  hash_final(&ctx);

  hash_squeeze(&ctx, h_1, instance->digest_size);
}

void h_1_expand(const aimer_instance_t* instance, const uint8_t* h_1,
                 GF epsilons[AIMER_T][AIMER_NUM_INPUT_SBOXES+1])
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init(&ctx, instance->digest_size);
  hash_update(&ctx, h_1, instance->digest_size);
  hash_final(&ctx);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    for (size_t ell = 0; ell < instance->aim_params.num_input_sboxes + 1; ell++)
    {
      hash_squeeze(&ctx, buffer, instance->field_size);
      GF_from_bytes(buffer, epsilons[repetition][ell]);
    }
  }
}

void h_2_commitment(const aimer_instance_t* instance, const uint8_t* salt,
                    const uint8_t* h_1, const GF* repetition_alpha_shares,
                    const GF v_shares[AIMER_T][AIMER_N], uint8_t* h_2)
{
  hash_instance ctx;
  uint8_t buffer[instance->field_size];

  memset(buffer, 0, sizeof(buffer));

  hash_init_prefix(&ctx, instance->digest_size, HASH_PREFIX_2);
  hash_update(&ctx, salt, instance->salt_size);
  hash_update(&ctx, h_1 , instance->digest_size);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    for (size_t party = 0; party < instance->num_MPC_parties; party++)
    {
      GF_to_bytes(repetition_alpha_shares
                  [repetition * instance->num_MPC_parties + party],
                  buffer);
      hash_update(&ctx, buffer, instance->field_size);

      GF_to_bytes(v_shares[repetition][party], buffer);
      hash_update(&ctx, buffer, instance->field_size);
    }
  }
  hash_final(&ctx);

  hash_squeeze(&ctx, h_2, instance->digest_size);
}

void h_2_expand(const aimer_instance_t* instance, const uint8_t* h_2,
                uint16_t* opened_parties)
{
  size_t num_squeeze_bytes = instance->num_MPC_parties > 256 ? 2 : 1;
  uint16_t mask = (1ULL << ceil_log2(instance->num_MPC_parties)) - 1;
  uint16_t party;

  hash_instance ctx;
  hash_init(&ctx, instance->digest_size);
  hash_update(&ctx, h_2, instance->digest_size);
  hash_final(&ctx);

  for (size_t repetition = 0; repetition < instance->num_repetitions;
       repetition++)
  {
    do
    {
      hash_squeeze(&ctx, (uint8_t *)&party, num_squeeze_bytes);
      party = le16toh(party) & mask;
    }
    while (party >= instance->num_MPC_parties);

    opened_parties[repetition] = party;
  }
}
