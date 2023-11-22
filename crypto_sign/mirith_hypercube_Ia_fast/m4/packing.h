/*
 * Copyright 2023 Carlo Sanna, Javier Verbel, and Floyd Zweydinger.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef PACKING_H
#define PACKING_H

#include "hash.h"
#include "matrix.h"
#include "params.h"

/* IMPORTANT: To sign a message, this scheme requires both the public key 'M'
 * and the secret key 'alpha', 'K', 'E'. However, SUPERCOP does not allow
 * one to pass the public key to the signing function. Therefore, we had
 * to incomporate all the matrices 'M', 'alpha', 'K', 'E' in the secret key. */

/* Pack the public key into 'pk'. */
void pack_public_key(uint8_t *pk, const ff_t M0[matrix_bytes_size(PAR_M, PAR_N)], const seed_t seed_pk);

/* Unpack the public key matrices 'M[i]' from 'pk'. */
void unpack_public_key(ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)], const uint8_t *pk);

/* Pack the secret key into 'sk'. */
void pack_secret_key(uint8_t *sk,
    const ff_t M0[matrix_bytes_size(PAR_M, PAR_N)], const seed_t seed_pk,
    const seed_t seed_sk);

/* Unpack the secret key 'M', 'a', 'K', 'E' from 'sk'. */
void unpack_secret_key(
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t a[matrix_bytes_size(PAR_K, 1)],
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)],
    const uint8_t *sk);

/* Pack the signature into 'sig' and write the length of the signature
 * over 'sig_len'. The length of the signature never exceeds CRYPTO_BYTES. */
void pack_signature(
    uint8_t *sig,
    size_t *sig_len,
    const hash_t salt,
    const hash_t hash1,
    const hash_t hash2,
    const int i_star[TAU],
    seed_t tree[TAU][TREE_N_NODES],
    hash_t com_star[TAU],
    ff_t a_last_shr[TAU][matrix_bytes_size(PAR_K, 1)],
    ff_t K_last_shr[TAU][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t C_last_shr[TAU][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t S_star[TAU][matrix_bytes_size(PAR_S, PAR_R)]);

/* Unpack the signature from 'sig'. If no error occurs, then return 'True'
 * and set 'sig_len' to be equal to the length of the signature.
 * Otherwise, return 'False'. */
int unpack_signature(
    hash_t salt,
    hash_t hash1,
    hash_t hash2,
    int i_star[TAU],
    seed_t packed_tree[TAU][TREE_HEIGHT],
    hash_t com_star[TAU],
    ff_t a_rnd_aux[TAU][matrix_bytes_size(PAR_K, 1)],
    ff_t K_rnd_aux[TAU][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t C_rnd_aux[TAU][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t S_rnd_star[TAU][matrix_bytes_size(PAR_S, PAR_R)],
    const uint8_t *sig,
    size_t *sig_len);

#endif
