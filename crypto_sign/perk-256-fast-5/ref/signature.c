
/**
 * @file signature.c
 * @brief Implementation of sign function
 */

#include "signature.h"
#include "parameters.h"

#include <stdio.h>
#include <string.h>
#include "common.h"
#include "crypto_memset.h"
#include "data_structures.h"
#include "parsing.h"
#include "randombytes.h"
#include "symmetric.h"
#include "theta_tree.h"
#include "verbose.h"
/**
 * @brief Generate the commitments
 *
 * @param [out] signature a pointer to a signature structure
 * @param [out,in] instances an array of PARAM_TAU instances
 * @param [in] pi a variable that containing the permutation pi
 * @param [in] H a variable containing the public matrix H
 */
static void sig_perk_gen_commitment(perk_signature_t *signature, perk_instance_t instances[PARAM_TAU], const perm_t pi,
                                    const mat_t H) {
    uint8_t rand_bytes[SEED_BYTES + SALT_BYTES] = {0};
    sig_perk_prg_state_t prg_state;
    vect1_t v = {0};
    seed_t mseed = {0};

    randombytes(rand_bytes, sizeof(rand_bytes));
    memcpy(mseed, rand_bytes, SEED_BYTES);
    memcpy(signature->salt, rand_bytes + SEED_BYTES, SALT_BYTES);
    sig_perk_prg_init(&prg_state, PRG1, signature->salt, mseed);

    SIG_PERK_VERBOSE_PRINT_uint8_t_array("mseed", mseed, SEED_BYTES);
    SIG_PERK_VERBOSE_PRINT_uint8_t_array("salt", signature->salt, SALT_BYTES);

    for (uint8_t i = 0; i < PARAM_TAU; ++i) {
        theta_t *tree = instances[i].theta_tree;
        sig_perk_prg(&prg_state, instances[i].theta_tree[0], sizeof(theta_t));
        sig_perk_expand_theta_tree(signature->salt, tree);
        sig_perk_gen_instance_commitments(&instances[i], signature->salt, i, pi);
        sig_perk_gen_vect_v(v, instances[i].pi_i, instances[i].v_i, pi);
        sig_perk_gen_instance_commitment_cmt_1(&instances[i], signature->salt, i, H, v);

        SIG_PERK_VERBOSE_PRINT_theta_cmt_1_and_v(instances[i].theta_tree[0], instances[i].cmt_1, v, i + 1);
    }
}

void sig_perk_gen_first_challenge(challenge_t challenge[PARAM_TAU], sig_perk_hash_state_t *saved_state, digest_t h1,
                                  const salt_t salt, const uint8_t *m, const uint64_t mlen, const uint8_t *pk_bytes,
                                  const perk_instance_t instances[PARAM_TAU]) {
    sig_perk_hash_state_t hash_state;
    sig_perk_prg_state_t prg_state;
    uint16_t tmp_kappa;

    sig_perk_hash_init(&hash_state, salt, NULL, NULL);
    sig_perk_hash_update(&hash_state, m, mlen);
    sig_perk_hash_update(&hash_state, pk_bytes, PUBLIC_KEY_BYTES);
    memcpy(saved_state, &hash_state, sizeof(hash_state));
    for (int i = 0; i < PARAM_TAU; ++i) {
        // absorb cmt_1_N to cmt_1_1
        // Commitments are stored in reverse order: from element (PARAM_N-1) down to element 0.
        sig_perk_hash_update(&hash_state, (uint8_t *)instances[i].cmt_1_i, sizeof(cmt_t) * PARAM_N);
        sig_perk_hash_update(&hash_state, (uint8_t *)instances[i].cmt_1, sizeof(cmt_t));
    }
    sig_perk_hash_final(&hash_state, h1, H1);

    sig_perk_prg_init(&prg_state, PRG1, NULL, h1);
    for (int i = 0; i < PARAM_TAU; ++i) {
        uint16_t nonzero_check = 0;
        do {
            for (int j = 0; j < PARAM_T; ++j) {
                do {
                    sig_perk_prg(&prg_state, (uint8_t *)&tmp_kappa, sizeof(tmp_kappa));
                    tmp_kappa = tmp_kappa & PARAM_Q_MASK;
                } while (tmp_kappa >= PARAM_Q);  // 0 <= tmp_kappa < PARAM_Q
                challenge[i].kappa[j] = tmp_kappa;
                nonzero_check |= tmp_kappa;
            }
        } while (!nonzero_check);
    }
}

/**
 * @brief Generate the first response
 *
 * @param [out] instances an array of PARAM_TAU instances
 * @param [in] challenges an array of PARAM_TAU challenges
 * @param [in] x an array of PARAM_T public values x_i
 */
static void sig_perk_gen_first_response(perk_instance_t instances[PARAM_TAU], const challenge_t challenges[PARAM_TAU],
                                        const vect1_t x[PARAM_T]) {
    for (int i = 0; i < PARAM_TAU; ++i) {
        vect1_t tmp;
        // compute s_0
        sig_perk_vect1_mult_scalar_vect(instances[i].s_i[0], challenges[i].kappa[0], x[0]);
        for (int j = 1; j < PARAM_T; j++) {
            sig_perk_vect1_mult_scalar_vect(tmp, challenges[i].kappa[j], x[j]);
            sig_perk_vect1_add(instances[i].s_i[0], instances[i].s_i[0], tmp);
        }
        for (int j = 0; j < PARAM_N; ++j) {
            sig_perk_perm_vect_permute(tmp, instances[i].pi_i[j], instances[i].s_i[j]);
            sig_perk_vect1_add(instances[i].s_i[j + 1], tmp, instances[i].v_i[j]);
        }
    }
}

void sig_perk_gen_second_challenge(digest_t h2, challenge_t challenge[PARAM_TAU], sig_perk_hash_state_t *saved_state,
                                   const digest_t h1, perk_instance_t instances[PARAM_TAU]) {
    sig_perk_prg_state_t prg_state;
    uint16_t tmp_alpha;

    sig_perk_hash_update(saved_state, h1, sizeof(digest_t));
    for (int i = 0; i < PARAM_TAU; ++i) {
        sig_perk_hash_update(saved_state, (uint8_t *)instances[i].s_i[1], sizeof(vect1_t) * PARAM_N);
    }
    sig_perk_hash_final(saved_state, h2, H2);

    sig_perk_prg_init(&prg_state, PRG1, NULL, h2);
    for (int i = 0; i < PARAM_TAU; ++i) {
        sig_perk_prg(&prg_state, (uint8_t *)&tmp_alpha, sizeof(tmp_alpha));
        tmp_alpha = (tmp_alpha & PARAM_N_MASK) + 1;  // 0 < tmp_alpha <= PARAM_N
        challenge[i].alpha = tmp_alpha;
    }
}

/**
 * @brief Generate second response
 *
 * @param [out] signature a pointer to a signature structure
 * @param [in] challenges an array of PARAM_TAU challenges
 * @param [in] instances an array of PARAM_TAU instances
 */
static void sig_perk_gen_second_response(perk_signature_t *signature, const challenge_t challenges[PARAM_TAU],
                                         const perk_instance_t instances[PARAM_TAU]) {
    for (int i = 0; i < PARAM_TAU; ++i) {
        const uint16_t alpha = challenges[i].alpha;
        memcpy(signature->responses[i].z1, instances[i].s_i[alpha], sizeof(vect1_t));
        if (alpha != 1) {
            memcpy(signature->responses[i].z2_pi, instances[i].pi_i[0], sizeof(perm_t));
        } else {
            for (int j = 0; j < PARAM_N1; j++) {
                signature->responses[i].z2_pi[j] = j;
            }
        }
        sig_perk_get_theta_partial_tree_seeds(signature->responses[i].z2_theta, instances[i].theta_tree, alpha - 1);
        // Commitments are stored in reverse order: from element (PARAM_N-1) down to element 0.
        memcpy(signature->responses[i].cmt_1_alpha, instances[i].cmt_1_i[(PARAM_N - 1) - (alpha - 1)], sizeof(cmt_t));
    }
}

uint8_t sig_perk_sign(perk_signature_t *signature, const perk_private_key_t *sk, const uint8_t *message_bytes,
                      const uint64_t message_length) {
    perk_public_key_t pk;
    perk_instance_t instances[PARAM_TAU] = {0};
    perk_instance_t_array_init(instances, PARAM_TAU);
    challenge_t challenges[PARAM_TAU] = {0};
    sig_perk_hash_state_t saved_state;  // stores a copy of the Keccak state ofter absorbing m, pk_bytes and com1

    SIG_PERK_VERBOSE_PRINT_uint8_t_array("message m", message_bytes, message_length);

    if (EXIT_SUCCESS != sig_perk_public_key_from_bytes(&pk, sk->pk_bytes)) {
        return EXIT_FAILURE;
    }
    sig_perk_gen_commitment(signature, instances, sk->pi, (const vect1_t *)pk.H);
    sig_perk_gen_first_challenge(challenges, &saved_state, signature->h1, signature->salt, message_bytes,
                                 message_length, sk->pk_bytes, instances);
    sig_perk_gen_first_response(instances, challenges, (const vect1_t *)pk.x);
    sig_perk_gen_second_challenge(signature->h2, challenges, &saved_state, signature->h1, instances);
    sig_perk_gen_second_response(signature, challenges, instances);

    SIG_PERK_VERBOSE_PRINT_challenges(challenges);
    SIG_PERK_VERBOSE_PRINT_signature(signature);

    return EXIT_SUCCESS;
}