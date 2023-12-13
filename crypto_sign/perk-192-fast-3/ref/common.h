
/**
 * @file common.h
 * @brief Header file for common.c
 */

#ifndef SIG_PERK_COMMON_H
#define SIG_PERK_COMMON_H

#include "data_structures.h"
#include "permutation.h"
#include "symmetric.h"

/**
 * @brief Generate a set of pi_i and v_i
 *
 * @param [out] pi_i a pointer to permutations
 * @param [out] v_i a pointer to vectors v_i
 * @param [in] salt a salt
 * @param [in] theta_tree a tree of seeds
 */
void sig_perk_gen_pi_i_and_v_i(perm_t *pi_i, vect1_t *v_i, const salt_t salt, const perk_theta_seeds_tree_t theta_tree);

/**
 * @brief Generate the vector v
 *
 * @param [out] o a variable containing the vector v
 * @param [in] pi_i a pointer to permutations pi_i
 * @param [in] v_i a pointer to vectors v_i
 * @param [in] pi a variable containing the secret permutation pi
 */
void sig_perk_gen_vect_v(vect1_t o, perm_t *pi_i, vect1_t *v_i, const perm_t pi);

/**
 * @brief Compute commitments cmt_1,i (for i in 1,...,N) for one round of the scheme
 *
 * @param [out,in] instance pointer to struct instance
 * @param [in] salt a salt
 * @param [in] tau an integer that is the index in 0,...,tau
 */
void sig_perk_gen_commitment_cmt_1_i(perk_instance_t *instance, salt_t salt, uint8_t tau);

/**
 * @brief Generate commitments cmt_1,i (for i in 1,...,N) for one round of the scheme
 *
 * @param [out,in] instance a pointer to an instance
 * @param [in] salt a salt
 * @param [in] tau an integer that is the index in 0,...,tau
 * @param [in] pi a variable containing the secret permutation pi
 */
void sig_perk_gen_instance_commitments(perk_instance_t *instance, salt_t salt, uint8_t tau, const perm_t pi);

/**
 * @brief Compute commitment cmt_1 for one round of the scheme
 *
 * @param [out,in] instance a pointer to an instance
 * @param [in] salt a salt
 * @param [in] tau an integer that is the index in 0,...,tau
 * @param [in] H a variable containing the public matrix H
 * @param [in] v a variable containing to vectors v_i
 */
void sig_perk_gen_instance_commitment_cmt_1(perk_instance_t *instance, const salt_t salt, uint8_t tau, const mat_t H,
                                            const vect1_t v);

#endif