
/**
 * @file verify.c
 * @brief Implementation of verify function
 */

#include "verify.h"
#include "common.h"
#include "parsing.h"
#include "string.h"

/**
 * @brief Compute the s_i in the verification step
 *
 * @param [out] s_i a pointer to the s_i
 * @param [in] pi_i a pointer to the permutations pi_i
 * @param [in] v_i a pointer to the vectors v_i
 * @param [in] alpha an integer that is the challenge alpha
 */
static void sig_perk_verify_compute_s_i(vect1_t *s_i, perm_t *pi_i, vect1_t *v_i, const uint16_t alpha) {
    vect1_t tmp;
    for (int i = 0; i < PARAM_N; ++i) {
        if ((i + 1) != alpha) {
            // compute pi_i(s_i - 1)
            for (int j = 0; j < PARAM_N1; ++j) {
                tmp[pi_i[i][j]] = s_i[i][j];
            }
            sig_perk_vect1_add(s_i[i + 1], tmp, v_i[i]);
        }
    }
}

/**
 * @brief Compute the commitments cmt_1_i for i in 1,...,N
 *
 * @param [out] instance a pointer to an instance
 * @param [in] response a pointer to a response
 * @param [in] challenge a variable containing the challenge
 * @param [in] x an array containing x_i for i in 1,...,t
 * @param [in] tau an integer that is an index in 0,...,tau
 * @param [in] salt a variable containing the salt
 */
static void sig_perk_verify_gen_instance_commitments(perk_instance_t *instance, const perk_response_t *response,
                                                     const challenge_t challenge, const vect1_t x[PARAM_T],
                                                     const uint8_t tau, salt_t salt) {
    vect1_t tmp;
    // generate a partial tree with N-1 seeds
    sig_perk_expand_theta_partial_tree(salt, instance->theta_tree, response->z2_theta, challenge.alpha - 1);

    sig_perk_gen_pi_i_and_v_i(instance->pi_i, instance->v_i, salt, (const theta_t *)instance->theta_tree);
    if (challenge.alpha != 1) {
        memcpy(instance->pi_i[0], response->z2_pi, sizeof(perm_t));
    }
    //  compute s_0
    sig_perk_vect1_mult_scalar_vect(instance->s_i[0], challenge.kappa[0], x[0]);
    for (int i = 1; i < PARAM_T; i++) {
        sig_perk_vect1_mult_scalar_vect(tmp, challenge.kappa[i], x[i]);
        sig_perk_vect1_add(instance->s_i[0], instance->s_i[0], tmp);
    }
    // copy s_alpha
    memcpy(instance->s_i[challenge.alpha], (uint8_t *)response->z1, sizeof(vect1_t));
    sig_perk_verify_compute_s_i(instance->s_i, instance->pi_i, instance->v_i, challenge.alpha);
    sig_perk_gen_commitment_cmt_1_i(instance, salt, tau);
    // compute cmt_1,1
    if (challenge.alpha != 1) {
        uint8_t idx = 0;
        sig_perk_hash_state_t state;
        uint8_t pi_1_bytes[PARAM_N1];
        for (int j = 0; j < PARAM_N1; ++j) {
            pi_1_bytes[j] = (uint8_t)instance->pi_i[0][j];
        }
        sig_perk_hash_init(&state, salt, &tau, &idx);
        sig_perk_hash_update(&state, pi_1_bytes, PARAM_N1);
        sig_perk_hash_update(&state, instance->theta_tree[THETA_SEEDS_OFFSET], sizeof(theta_t));
        // Commitments are stored in reverse order: from element (PARAM_N-1) down to element 0.
        sig_perk_hash_final(&state, instance->cmt_1_i[PARAM_N - 1], H0);
    }
    // copy cmt_1_alpha
    // Commitments are stored in reverse order: from element (PARAM_N-1) down to element 0.
    memcpy(instance->cmt_1_i[(PARAM_N - 1) - (challenge.alpha - 1)], (uint8_t *)response->cmt_1_alpha, sizeof(cmt_t));
}

/**
 * @brief Compute the commitment cmt_1
 *
 * @param [out] cmt_1 a variable that is the commitment
 * @param [in] public_key a pointer to the public key
 * @param [in] instance a pointer to an instance
 * @param [in] kappa a variable containing the challenges kappa_i for i in 0,...,t
 * @param [in] salt a variable containing the salt
 * @param [in] tau an integer that is an index in 0,...,tau
 */
static void sig_perk_verify_gen_instance_commitment_cmt(cmt_t cmt_1, const perk_public_key_t *public_key,
                                                        const perk_instance_t *instance, const uint16_t kappa[PARAM_T],
                                                        const salt_t salt, const uint8_t tau) {
    vect2_t tmp1, tmp2, tmp3;
    sig_perk_hash_state_t state;

    sig_perk_mat_vect_mul(tmp1, public_key->H, instance->s_i[PARAM_N]);
    sig_perk_vect2_mult_scalar_vect(tmp2, kappa[0], public_key->y[0]);
    for (int i = 1; i < PARAM_T; i++) {
        sig_perk_vect2_mult_scalar_vect(tmp3, kappa[i], public_key->y[i]);
        sig_perk_vect2_add(tmp2, tmp2, tmp3);
    }
    sig_perk_vect2_sub(tmp1, tmp1, tmp2);

    sig_perk_hash_init(&state, salt, &tau, NULL);
    sig_perk_hash_update(&state, (uint8_t *)tmp1, sizeof(tmp1));
    sig_perk_hash_final(&state, cmt_1, H0);
}

/**
 * @brief Compute the digest h1
 *
 * @param [out] h1 a variable containing the digest h1
 * @param [out,in] saved_state pointer to Keccak state after absorbing salt, m and pk_bytes
 * @param [in] salt a variable containing the salt
 * @param [in] message_bytes a pointer to a string containing the message
 * @param [in] message_length an integer that is the size of the message
 * @param [in] pk_bytes a string containing the public key
 * @param [in] instances an array containing tau instances
 */
static void sig_perk_verify_compute_h1(digest_t h1, sig_perk_hash_state_t *saved_state, const salt_t salt,
                                       const uint8_t *message_bytes, const uint64_t message_length,
                                       const uint8_t *pk_bytes, const perk_instance_t instances[PARAM_TAU]) {
    sig_perk_hash_state_t hash_state;

    sig_perk_hash_init(&hash_state, salt, NULL, NULL);
    sig_perk_hash_update(&hash_state, message_bytes, message_length);
    sig_perk_hash_update(&hash_state, pk_bytes, PUBLIC_KEY_BYTES);
    memcpy(saved_state, &hash_state, sizeof(hash_state));
    for (int i = 0; i < PARAM_TAU; ++i) {
        // absorb cmt_1_N to cmt_1_1
        // Commitments are stored in reverse order: from element (PARAM_N-1) down to element 0.
        sig_perk_hash_update(&hash_state, (uint8_t *)instances[i].cmt_1_i, sizeof(cmt_t) * PARAM_N);
        sig_perk_hash_update(&hash_state, (uint8_t *)instances[i].cmt_1, sizeof(cmt_t));
    }
    sig_perk_hash_final(&hash_state, h1, H1);
}

/**
 * @brief Compute the digest h2
 *
 * @param [out] h2 a variable containing the digest h2
 * @param [out,in] saved_state pointer to Keccak state after absorbing salt, m and pk_bytes
 * @param [in] h1 a variable containing the digest h1
 * @param [in] instances an array containing tau instances
 */
static void sig_perk_verify_compute_h2(digest_t h2, sig_perk_hash_state_t *saved_state, const digest_t h1,
                                       const perk_instance_t instances[PARAM_TAU]) {
    sig_perk_hash_update(saved_state, h1, sizeof(digest_t));
    for (int i = 0; i < PARAM_TAU; ++i) {
        sig_perk_hash_update(saved_state, (uint8_t *)instances[i].s_i[1], sizeof(vect1_t) * PARAM_N);
    }
    sig_perk_hash_final(saved_state, h2, H2);
}

int sig_perk_verify(perk_signature_t *signature, const challenge_t challenge[PARAM_TAU], const uint8_t *message_bytes,
                    const uint64_t message_length, const uint8_t *pk_bytes) {
    perk_public_key_t public_key = {0};
    perk_instance_t instances[PARAM_TAU] = {0};
    perk_instance_t_array_init(instances, PARAM_TAU);
    digest_t h1_prime, h2_prime;
    sig_perk_hash_state_t saved_state;

    for (int i = 0; i < PARAM_TAU; ++i) {
        if (challenge[i].alpha == 1) {
            for (int j = 0; j < PARAM_N1; j++) {
                if (signature->responses[i].z2_pi[j] != j) {
                    return EXIT_FAILURE;
                }
            }
        }
    }
    if (EXIT_SUCCESS != sig_perk_public_key_from_bytes(&public_key, pk_bytes)) {
        return EXIT_FAILURE;
    }
    for (int i = 0; i < PARAM_TAU; ++i) {
        sig_perk_verify_gen_instance_commitments(&instances[i], &signature->responses[i], challenge[i],
                                                 (const vect1_t *)public_key.x, i, signature->salt);
        sig_perk_verify_gen_instance_commitment_cmt(instances[i].cmt_1, &public_key, &instances[i], challenge[i].kappa,
                                                    signature->salt, i);
    }

    sig_perk_verify_compute_h1(h1_prime, &saved_state, signature->salt, message_bytes, message_length, pk_bytes,
                               instances);
    sig_perk_verify_compute_h2(h2_prime, &saved_state, h1_prime, instances);

    if (memcmp((uint8_t *)h1_prime, (uint8_t *)signature->h1, sizeof(digest_t))) {
        return EXIT_FAILURE;
    }
    if (memcmp((uint8_t *)h2_prime, (uint8_t *)signature->h2, sizeof(digest_t))) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}