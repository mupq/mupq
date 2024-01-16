// -----------------------------------------------------------------------------
// File Name   : aimer.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aim.h"
#include "aimer.h"
#include "aimer_internal.h"
#include "randombytes.h"
#include "tree.h"

#include <stdlib.h>

int _aimer_sign();
int _aimer_verify();
int serialize_signature();
int deserialize_signature();

int aimer_keygen(aimer_params_t param, aimer_publickey_t* public_key,
                 aimer_privatekey_t* private_key)
{
  if (!public_key || !private_key)
  {
    return -1;
  }

  const aimer_instance_t* instance = aimer_instance_get(param);
  if (!instance)
  {
    return -1;
  }

  randombytes(public_key->data, instance->aim_params.block_size);
  randombytes(private_key->data, instance->aim_params.block_size);

  aim(private_key->data, public_key->data,
      &public_key->data[instance->aim_params.block_size]);

  public_key->params = param;
  return 0;
}

int _aimer_sign(const aimer_instance_t*   instance,
                const aimer_publickey_t*  public_key,
                const aimer_privatekey_t* private_key,
                const uint8_t* message, const size_t message_len,
                signature_t* sig)
{
  const size_t L                = instance->aim_params.num_input_sboxes;
  const size_t block_size       = instance->aim_params.block_size;
  const size_t random_tape_size = AIMER_BLOCK_SIZE +
                                  AIMER_NUM_INPUT_SBOXES * AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE;
  const size_t tau              = instance->num_repetitions;
  const size_t N                = AIMER_N;
  const size_t digest_size      = instance->digest_size;

  int ret = 0;
  GF pt_GF = {0,};
  GF vector_b = {0,};

  uint8_t iv[block_size];
  uint8_t pt[block_size];
  uint8_t ct[block_size];

  /////////////////////////////////////////////////////////////////////////////
  // Phase 1: Committing to the seeds and the execution views of the parties.
  /////////////////////////////////////////////////////////////////////////////

  memcpy(pt, private_key->data, block_size);
  memcpy(iv, public_key->data, block_size);  
  memcpy(ct, public_key->data + block_size, block_size);
  GF_from_bytes(pt, pt_GF);

  randombytes(sig->salt, instance->salt_size);

  GF sbox_outputs[(AIMER_NUM_INPUT_SBOXES + 1)];
  compute_sbox_outputs(pt, sbox_outputs);

  GF matrix_A[AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS];
  generate_matrix_LU(iv, matrix_A, vector_b);

  tree_t seed_trees[tau];
  uint8_t party_seed_commitments[AIMER_T * AIMER_N * AIMER_DIGEST_SIZE];

  random_tape_t random_tapes;
  random_tapes.random_tape_size = random_tape_size;

  uint8_t master_seed[AIMER_SEED_SIZE];
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    randombytes(master_seed, instance->seed_size);
    make_seed_tree(&seed_trees[repetition], master_seed, instance->seed_size, sig->salt,
                     instance->salt_size, N, repetition);

    size_t party = 0;
    uint8_t *seed0, *seed1, *seed2, *seed3;
    for (; party < (N / 4) * 4; party += 4)
    {
      seed0 = get_leaf(&seed_trees[repetition], party);
      seed1 = get_leaf(&seed_trees[repetition], party + 1);
      seed2 = get_leaf(&seed_trees[repetition], party + 2);
      seed3 = get_leaf(&seed_trees[repetition], party + 3);

      commit_to_seed_and_expand_tape_x4(instance, seed0, seed1, seed2, seed3,
                                        sig->salt, repetition, party,
                                        party_seed_commitments,
                                        &random_tapes);
    }

    for (; party < N; party++)
    {
      seed0 = get_leaf(&seed_trees[repetition], party);

      commit_to_seed_and_expand_tape(instance, seed0, sig->salt,
                                     repetition, party,
                                     party_seed_commitments,
                                     &random_tapes);
    }
  }

  uint8_t repetition_shared_pt[AIMER_T*AIMER_N*AIMER_BLOCK_SIZE];

  GF repetition_shared_x     [AIMER_T * AIMER_N  * (AIMER_NUM_INPUT_SBOXES + 1)];
  GF repetition_shared_z     [AIMER_T * AIMER_N  * (AIMER_NUM_INPUT_SBOXES + 1)];
  GF repetition_shared_dot_a [AIMER_T * AIMER_N];
  GF repetition_shared_dot_c [AIMER_T * AIMER_N];

  uint8_t pt_delta[block_size];
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    // Generate sharing of secret key
    memset(pt_delta, 0, block_size);
    for (size_t party = 0; party < N; party++)
    {
      uint8_t* shared_pt =
        repetition_shared_pt + (repetition * N + party) * block_size;
      uint8_t* random_share =
        random_tapes.tape + (repetition * N + party) * random_tape_size;

      memcpy(shared_pt, random_share, block_size);
      for (size_t i = 0; i < block_size; i++)
      {
        pt_delta[i] ^= random_share[i];
      }
    }

    for (size_t i = 0; i < block_size; i++)
    {
      pt_delta[i] ^= pt[i];
    }

    // Fix first share
    uint8_t* first_pt_share = 
      repetition_shared_pt + (repetition * N * block_size);
    for (size_t i = 0; i < block_size; i++)
    {
      first_pt_share[i] ^= pt_delta[i];
    }
    memcpy(proof->pt_delta, pt_delta, block_size);

    memset(proof->z_delta, 0, L * instance->field_size);
    for (size_t party = 0; party < N; party++)
    {
      GF* shared_z =
        repetition_shared_z + (repetition * N + party) * (L + 1);

      uint8_t* random_z_shares =
        random_tapes.tape + (repetition * N + party) * random_tape_size +
        block_size;

      for (size_t ell = 0; ell < L; ell++)
      {
        GF_from_bytes(random_z_shares + ell * block_size, shared_z[ell]);
        GF_add(proof->z_delta[ell], shared_z[ell], proof->z_delta[ell]);
      }
    }

    GF* first_shared_z = repetition_shared_z + (repetition * N * (L + 1));
    for (size_t ell = 0; ell < L; ell++)
    {
      GF_add(proof->z_delta[ell], sbox_outputs[ell], proof->z_delta[ell]);
      GF_add(first_shared_z[ell], proof->z_delta[ell], first_shared_z[ell]);
    }

    // Get shares of sbox inputs by simulating MPC execution
    GF* shared_x_offset = repetition_shared_x + (repetition * N * (L + 1));
    GF* shared_z_offset = first_shared_z;

    uint8_t* shared_pt_offset =
      repetition_shared_pt + (repetition * N * block_size);

    aim_mpc(shared_pt_offset, matrix_A, vector_b, ct,
            N, shared_z_offset, shared_x_offset);
  }

  GF a = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    GF_set0(a);
    GF_set0(proof->c_delta);
    for (size_t party = 0; party < N; party++)
    {
      uint8_t* random_triple = random_tapes.tape +
        ((repetition * N + party) * random_tape_size +
          block_size + L * instance->field_size);

      GF* a_share = repetition_shared_dot_a + (repetition * N + party);
      GF_from_bytes(random_triple, a_share[0]);
      GF_add(a, a_share[0], a);

      GF* c_share = repetition_shared_dot_c + (repetition * N + party);
      GF_from_bytes(random_triple + instance->field_size, c_share[0]);
      GF_add(proof->c_delta, c_share[0], proof->c_delta);
    }

    // Calculate c_delta = -c + a*b, c is negated above already
    // Calculate c_delta that fixes the dot triple
    GF_mul(a, pt_GF, a);
    GF_add(a, proof->c_delta, proof->c_delta);

    // Fix party 0's share
    GF_add(repetition_shared_dot_c[repetition * N], proof->c_delta,
           repetition_shared_dot_c[repetition * N]);
  }

  //////////////////////////////////////////////////////////////////////////////
  // Phase 2: Challenging the checking protocol.
  //////////////////////////////////////////////////////////////////////////////
  // Commit to salt, (all commitments of parties seeds,
  // pt_delta, z_delta, c_delta) for all repetitions
  GF epsilons[AIMER_T][AIMER_NUM_INPUT_SBOXES+1];
  h_1_commitment(instance, sig, public_key, message, message_len,
                 party_seed_commitments, sig->h_1);

  // Expand challenge hash to epsilon values
  h_1_expand(instance, sig->h_1, epsilons);

  //////////////////////////////////////////////////////////////////////////////
  // Phase 3: Committing to the simulation of the checking protocol.
  //////////////////////////////////////////////////////////////////////////////

  GF repetition_alpha_shares[AIMER_T*AIMER_N];
  GF v_shares[AIMER_T][AIMER_N];

  GF alpha = {0,}, pt_share = {0,}, temp = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    GF* x_shares     = repetition_shared_x     + (repetition * N) * (L + 1);
    GF* z_shares     = repetition_shared_z     + (repetition * N) * (L + 1);
    GF* dot_a_share  = repetition_shared_dot_a + (repetition * N);
    GF* dot_c_share  = repetition_shared_dot_c + (repetition * N);
    GF* alpha_share  = repetition_alpha_shares + (repetition * N);

    // Execute sacrificing check protocol
    // alpha = eps_i * x_i + a
    GF_set0(alpha);
    for (size_t party = 0; party < N; party++)
    {
      GF_mul(x_shares[party * (L + 1)], epsilons[repetition][0], temp);
      GF_add(dot_a_share[party], temp, alpha_share[party]);
      for (size_t ell = 1; ell < (L + 1); ell++)
      {
        GF_mul(x_shares[party * (L + 1) + ell], epsilons[repetition][ell], temp);
        GF_add(alpha_share[party], temp, alpha_share[party]);
      }
      GF_add(alpha, alpha_share[party], alpha);
    }

    // v^i = - c^i + dot(alpha, y^i) - dot(eps, z^i)
    for (size_t party = 0; party < N; party++)
    {
      GF_from_bytes(repetition_shared_pt +
                    (repetition * N + party) * instance->field_size,
                    pt_share);

      GF_mul(alpha, pt_share, temp);
      GF_add(dot_c_share[party], temp, v_shares[repetition][party]);
      for (size_t ell = 0; ell < (L + 1); ell++)
      {
        GF_mul(epsilons[repetition][ell], z_shares[party * (L + 1) + ell], temp);
        GF_add(v_shares[repetition][party], temp, v_shares[repetition][party]);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Phase 4: Challenging the views of the MPC protocol.
  //////////////////////////////////////////////////////////////////////////////
  h_2_commitment(instance, sig->salt, sig->h_1, repetition_alpha_shares,
                 v_shares, sig->h_2);

  uint16_t missing_parties[AIMER_T];
  h_2_expand(instance, sig->h_2, missing_parties);

  //////////////////////////////////////////////////////////////////////////////
  // Phase 5: Opening the views of the MPC and checking protocols.
  //////////////////////////////////////////////////////////////////////////////
  // Build signature
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    size_t missing_party = missing_parties[repetition];
    reveal_all_but(&proof->reveal_list, &seed_trees[repetition], missing_party);

    memcpy(proof->missing_commitment,
           party_seed_commitments +
           ((repetition * N + missing_party) * digest_size),
           digest_size);

    GF_copy(repetition_alpha_shares[repetition * N + missing_party],
            proof->missing_alpha_share);
  }

  return ret;
}

int serialize_signature(const aimer_instance_t* instance,
                        const signature_t* sig,
                        uint8_t* signature, size_t* signature_len)
{
  const size_t L   = instance->aim_params.num_input_sboxes;
  const size_t tau = instance->num_repetitions;
  const size_t digest_size = instance->digest_size;
  uint8_t* signature_offset = signature;

  memcpy(signature_offset, sig->salt, instance->salt_size);
  signature_offset += instance->salt_size;

  memcpy(signature_offset, sig->h_1, digest_size);
  signature_offset += digest_size;

  memcpy(signature_offset, sig->h_2, digest_size);
  signature_offset += digest_size;

  size_t reveal_list_size =
    ceil_log2(instance->num_MPC_parties) * instance->seed_size;

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];

    memcpy(signature_offset, proof->reveal_list.data, reveal_list_size);
    signature_offset += reveal_list_size;

    memcpy(signature_offset, proof->missing_commitment, digest_size);
    signature_offset += digest_size;

    memcpy(signature_offset, proof->pt_delta,
           instance->aim_params.block_size);
    signature_offset += instance->aim_params.block_size;

    for (size_t ell = 0; ell < L; ell++)
    {
      GF_to_bytes(proof->z_delta[ell], signature_offset);
      signature_offset += instance->field_size;
    }

    GF_to_bytes(proof->c_delta, signature_offset);
    signature_offset += instance->field_size;

    GF_to_bytes(proof->missing_alpha_share, signature_offset);
    signature_offset += instance->field_size;
  }

  *signature_len = (int)(signature_offset - signature);
  if (*signature_len > AIMER_MAX_SIGNATURE_SIZE)
  {
    return -1;
  }
  return 0;
}

int aimer_sign(const aimer_publickey_t*  public_key,
               const aimer_privatekey_t* private_key,
               const uint8_t* message, const size_t message_len,
               uint8_t* signature, size_t* signature_len)
{
  int ret = 0;
  signature_t sig;
  memset(&sig, 0, sizeof(signature_t));
  aimer_params_t params = public_key->params;

  const aimer_instance_t *instance = aimer_instance_get(params);
  if (!instance)
  {
    return -1;
  }

  ret =
    _aimer_sign(instance, public_key, private_key, message, message_len, &sig);
  if (ret == -1)
  {
    return -1;
  }

  ret = serialize_signature(instance, &sig, signature, signature_len);
  if (ret == -1)
  {
    return -1;
  }

  return ret;
}

int deserialize_signature(const aimer_instance_t* instance,
                          const uint8_t* signature, const size_t signature_len,
                          signature_t* sig, uint16_t* missing_parties)
{
  const size_t L   = instance->aim_params.num_input_sboxes;
  const size_t tau = instance->num_repetitions;
  const size_t digest_size = instance->digest_size;
  size_t offset = 0;

  memcpy(sig->salt, signature, instance->salt_size);
  offset += instance->salt_size;

  memcpy(sig->h_1, signature + offset, digest_size);
  offset += digest_size;
  memcpy(sig->h_2, signature + offset, digest_size);
  offset += digest_size;

  h_2_expand(instance, sig->h_2, missing_parties);

  size_t reveal_list_size = ceil_log2(instance->num_MPC_parties);

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    proof_t* proof = &sig->proofs[repetition];

    proof->reveal_list.missing_leaf = missing_parties[repetition];

    memcpy(proof->reveal_list.data,
           signature + offset, reveal_list_size * instance->seed_size);
    offset += reveal_list_size * instance->seed_size;

    memcpy(proof->missing_commitment,
           signature + offset, digest_size);
    offset += digest_size;

    memcpy(proof->pt_delta,
           signature + offset, instance->aim_params.block_size);
    offset += instance->aim_params.block_size;

    for (size_t ell = 0; ell < L; ell++)
    {
      GF_from_bytes(signature + offset, proof->z_delta[ell]);
      offset += instance->field_size;
    }

    GF_from_bytes(signature + offset, proof->c_delta);
    offset += instance->field_size;

    GF_from_bytes(signature + offset, proof->missing_alpha_share);
    offset += instance->field_size;
  }

  if (offset != signature_len)
  {
    return -1;
  }

  return 0;
}

int _aimer_verify(const aimer_instance_t*  instance,
                  const aimer_publickey_t* public_key,
                  const signature_t* sig,
                  const uint16_t* missing_parties,
                  const uint8_t* message, const size_t message_len)
{
  const size_t L                = instance->aim_params.num_input_sboxes;
  const size_t block_size       = instance->aim_params.block_size;
  const size_t random_tape_size = AIMER_BLOCK_SIZE +
                                  AIMER_NUM_INPUT_SBOXES * AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE +
                                  AIMER_FIELD_SIZE;
  const size_t tau              = instance->num_repetitions;
  const size_t N                = instance->num_MPC_parties;
  const size_t digest_size      = instance->digest_size;
  const size_t seed_size        = instance->seed_size;

  int ret = 0;

  // h_1 expansion,
  // h_2 expansion already happened in deserialize
  GF epsilons[AIMER_T][AIMER_NUM_INPUT_SBOXES+1];
  h_1_expand(instance, sig->h_1, epsilons);

  // Rebuild seed trees
  tree_t seed_trees[tau];
  uint8_t party_seed_commitments[AIMER_T * AIMER_N * AIMER_DIGEST_SIZE];

  random_tape_t random_tapes;
  random_tapes.random_tape_size = random_tape_size;

  uint8_t dummy[seed_size];
  memset(dummy, 0, seed_size);
  uint8_t *seed0, *seed1, *seed2, *seed3;
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];
    reconstruct_seed_tree(&seed_trees[repetition], &proof->reveal_list, sig->salt,
                            instance->salt_size, N, repetition);
    size_t party = 0;
    for (; party < (N / 4) * 4; party += 4)
    {
      seed0 = get_leaf(&seed_trees[repetition], party);
      seed1 = get_leaf(&seed_trees[repetition], party + 1);
      seed2 = get_leaf(&seed_trees[repetition], party + 2);
      seed3 = get_leaf(&seed_trees[repetition], party + 3);
      if (seed0 == NULL)
      {
        seed0 = dummy;
      }
      if (seed1 == NULL)
      {
        seed1 = dummy;
      }
      if (seed2 == NULL)
      {
        seed2 = dummy;
      }
      if (seed3 == NULL)
      {
        seed3 = dummy;
      }
      commit_to_seed_and_expand_tape_x4(instance, seed0, seed1, seed2, seed3,
                                        sig->salt, repetition, party,
                                        party_seed_commitments,
                                        &random_tapes);
    }

    for (; party < N; party++)
    {
      seed0 = get_leaf(&seed_trees[repetition], party);
      if (seed0 == NULL)
      {
        seed0 = dummy;
      }
      commit_to_seed_and_expand_tape(instance, seed0, sig->salt, repetition, party, 
                                     party_seed_commitments,
                                     &random_tapes);
    }
    memcpy(party_seed_commitments +
           (repetition * N + missing_parties[repetition]) * digest_size,
           proof->missing_commitment, digest_size);
  }

  // Recompute commitments to executions of block cipher
  uint8_t repetition_shared_pt[AIMER_T*AIMER_N*AIMER_BLOCK_SIZE];
  GF repetition_shared_x     [AIMER_T * AIMER_N  * (AIMER_NUM_INPUT_SBOXES + 1)];
  GF repetition_shared_z     [AIMER_T * AIMER_N  * (AIMER_NUM_INPUT_SBOXES + 1)];
  GF repetition_shared_dot_a [AIMER_T * AIMER_N];
  GF repetition_shared_dot_c [AIMER_T * AIMER_N];

  uint8_t iv[block_size];
  uint8_t ct[block_size];

  memcpy(iv, public_key->data, block_size);
  memcpy(ct, public_key->data + block_size, block_size);

  // generate linear layer for AIM mpc
  GF matrix_A[AIMER_NUM_INPUT_SBOXES][AIMER_NUM_BITS];
  GF vector_b = {0,};
  generate_matrix_LU(iv, matrix_A, vector_b);

  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];
    // Generate sharing of secret key
    for (size_t party = 0; party < N; party++)
    {
      uint8_t* shared_pt =
        repetition_shared_pt + (repetition * N + party) * block_size;
      uint8_t* random_share =
        random_tapes.tape + (repetition * N + party) * random_tape_size;

      memcpy(shared_pt, random_share, block_size);
    }

    // Fix first share
    uint8_t* first_pt_share = 
      repetition_shared_pt + (repetition * N * block_size);
    for (size_t j = 0; j < block_size; j++)
    {
      first_pt_share[j] ^= proof->pt_delta[j];
    }

    // Generate sharing of z values
    for (size_t party = 0; party < N; party++)
    {
      GF* shared_z = repetition_shared_z + (repetition * N + party) * (L + 1);
      uint8_t* random_z_shares =
        random_tapes.tape + (repetition * N + party) * random_tape_size +
        block_size;

      for (size_t ell = 0; ell < L; ell++)
      {
        GF_from_bytes(random_z_shares + ell * block_size, shared_z[ell]);
      }
    }

    // Fix first share
    GF* first_shared_z = repetition_shared_z + (repetition * N * (L + 1));
    for (size_t ell = 0; ell < L; ell++)
    {
      GF_add(first_shared_z[ell], proof->z_delta[ell], first_shared_z[ell]);
    }

    // Get shares of sbox inputs by executing MPC AIM
    GF* shared_x_offset = repetition_shared_x + (repetition * N * (L + 1));
    GF* shared_z_offset = first_shared_z;

    uint8_t* shared_pt_offset =
      repetition_shared_pt + (repetition * N * block_size);

    aim_mpc(shared_pt_offset, matrix_A, vector_b, ct,
            N, shared_z_offset, shared_x_offset);
  }

  // Recompute shares of dot product
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];
    // Also generate valid dot triple a,z,c and save c_delta
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_parties[repetition])
      {
        uint8_t* random_triple =
          random_tapes.tape +
          ((repetition * N + party) * random_tape_size +
          block_size + L * instance->field_size);

        GF* a_share = repetition_shared_dot_a + (repetition * N + party);
        GF_from_bytes(random_triple, a_share[0]);

        GF* c_share = repetition_shared_dot_c + (repetition * N + party);
        GF_from_bytes(random_triple + instance->field_size, c_share[0]);
      }
    }

    // Fix party 0's share
    if (0 != missing_parties[repetition])
    {
      // Adjust first share with delta from signature
      GF_add(repetition_shared_dot_c[repetition * N], proof->c_delta,
             repetition_shared_dot_c[repetition * N]);
    }
  }

  // Recompute views of sacrificing checks
  GF repetition_alpha_shares[AIMER_T*AIMER_N];
  GF v_shares[AIMER_T][AIMER_N]= {0};
  GF alpha = {0,};
  GF temp = {0,};
  for (size_t repetition = 0; repetition < tau; repetition++)
  {
    const proof_t* proof = &sig->proofs[repetition];
    size_t missing_party = missing_parties[repetition];

    GF_set0(alpha);
    for (size_t party = 0; party < N; party++)
    {
      GF* alpha_share = repetition_alpha_shares + (repetition * N + party);
      GF* a_share     = repetition_shared_dot_a + (repetition * N + party);
      GF* x_shares    = repetition_shared_x + (repetition * N + party) * (L + 1);

      if (party != missing_party)
      {
        GF_mul(x_shares[0], epsilons[repetition][0], temp);
        GF_add(a_share[0], temp, alpha_share[0]);
        for (size_t ell = 1; ell < (L + 1); ell++)
        {
          GF_mul(x_shares[ell], epsilons[repetition][ell], temp);
          GF_add(alpha_share[0], temp, alpha_share[0]);
        }
        GF_add(alpha, alpha_share[0], alpha);
      }
    }

    // Fill missing shares
    GF* missing_alpha =
      repetition_alpha_shares + (repetition * N + missing_party);

    GF_copy(proof->missing_alpha_share, missing_alpha[0]);
    GF_add(alpha, missing_alpha[0], alpha);

    // v^i = - c^i + dot(alpha, y) - dot(eps, z^i)
    GF pt_share = {0,};
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_party)
      {
        GF* z_shares = repetition_shared_z + (repetition * N + party) * (L + 1);
        GF* c_share  = repetition_shared_dot_c + (repetition * N + party);

        GF_from_bytes(repetition_shared_pt +
                      (repetition * N + party) * instance->field_size,
                      pt_share);

        GF_mul(alpha, pt_share, v_shares[repetition][party]);
        GF_add(v_shares[repetition][party], c_share[0],
               v_shares[repetition][party]);
        for (size_t ell = 0; ell < (L + 1); ell++)
        {
          GF_mul(epsilons[repetition][ell], z_shares[ell], temp);
          GF_add(v_shares[repetition][party], temp,
                 v_shares[repetition][party]);
        }
      }
    }

    // Calculate missing shares as 0 - sum_{i!=missing} v^i
    for (size_t party = 0; party < N; party++)
    {
      if (party != missing_party)
      {
        GF_add(v_shares[repetition][missing_party], v_shares[repetition][party],
               v_shares[repetition][missing_party]);
      }
    }
  }

  // Recompute h_1 and h_2
  uint8_t h_1_prime[digest_size];
  uint8_t h_2_prime[digest_size];

  h_1_commitment(instance, sig, public_key, message, message_len,
                 party_seed_commitments, h_1_prime);

  h_2_commitment(instance, sig->salt, h_1_prime, repetition_alpha_shares,
                 v_shares, h_2_prime);

  if (memcmp(h_1_prime, sig->h_1, digest_size) != 0)
  {
    ret = -1;
  }

  if (memcmp(h_2_prime, sig->h_2, digest_size) != 0)
  {
    ret = -1;
  }

  return ret;
}

int aimer_verify(const aimer_publickey_t* public_key,
                 const uint8_t* signature, const size_t signature_len,
                 const uint8_t* message, const size_t message_len)
{
  int ret = 0;
  signature_t sig;
  aimer_params_t params = public_key->params;

  const aimer_instance_t *instance = aimer_instance_get(params);
  if (!instance)
  {
    return -1;
  }

  uint16_t missing_parties[AIMER_T];

  ret = deserialize_signature(instance, signature, signature_len, &sig,
                              missing_parties);
  if (ret == -1)
  {
    return -1;
  }

  ret = _aimer_verify(instance, public_key, &sig, missing_parties, message,
                      message_len);
  if (ret == -1)
  {
    return -1;
  }
  return ret;
}
