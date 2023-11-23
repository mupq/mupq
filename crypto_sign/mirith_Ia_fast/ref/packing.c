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

#include <string.h>
#include "challenges.h"
#include "packing.h"
#include "seed_tree.h"

void pack_public_key(uint8_t *pk, const ff_t M0[matrix_bytes_size(PAR_M, PAR_N)], const seed_t seed_pk)
{
    uint8_t *ptr = pk + SEED_SIZE;

    /* Pack the public key seed. */
    memcpy(pk, seed_pk, SEED_SIZE);

    /* Pack the matrix 'M0'. */
    matrix_pack(&ptr, NULL, M0, PAR_M, PAR_N);
}

void unpack_public_key(ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)], const uint8_t *pk)
{
    uint32_t i;
    seed_t seed_pk;
    prng_t prng_pk;
    uint8_t *ptr = (uint8_t *)pk + SEED_SIZE;

    /* Unpack the public key seed. */
    memcpy(seed_pk, pk, SEED_SIZE);
    
    /* Unpack the matrix 'M[0]'. */
    matrix_unpack(M[0], &ptr, NULL, PAR_M, PAR_N);

    /* Initialize the PRNG. */
    prng_init(&prng_pk, NULL, seed_pk);

    /* Generate the matrices 'M[1], M[2], ..., M[PAR_K]'. */
    for (i = 1; i <= PAR_K; i++)
    {
        matrix_init_random(M[i], PAR_M, PAR_N, &prng_pk);
    }
}

void pack_secret_key(uint8_t *sk,
    const ff_t M0[matrix_bytes_size(PAR_M, PAR_N)],
    const seed_t seed_pk,
    const seed_t seed_sk)
{
    /* Pack the secret key seed. */
    memcpy(sk, seed_sk, SEED_SIZE);

    /* Pack the public key. */
    pack_public_key(sk + SEED_SIZE, M0, seed_pk);
}

void unpack_secret_key(
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t a[matrix_bytes_size(PAR_K, 1)],
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)],
    const uint8_t *sk)
{
    /* Variant 1. */
    
    seed_t seed_sk;
    prng_t prng_sk;
    ff_t E_R[matrix_bytes_size(PAR_M, PAR_R)];
    ff_t T[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];

    /* Unpack the secret key seed. */
    memcpy(seed_sk, sk, SEED_SIZE);

    /* Initialize the PRNG. */
    prng_init(&prng_sk, NULL, seed_sk);

    /* Initialize matrices 'alpha', 'K', 'E_R'. */
    matrix_init_random(a, PAR_K, 1, &prng_sk);
    matrix_init_random(K, PAR_R, PAR_N - PAR_R, &prng_sk);
    matrix_init_random(E_R, PAR_M, PAR_R, &prng_sk);

    /* Compute the matrix 'E'. */
    matrix_product(T, E_R, K, PAR_M, PAR_R, PAR_N - PAR_R);
    matrix_horizontal_concatenation(E, T, E_R, PAR_M, PAR_N - PAR_R, PAR_R);

    /* Unpack the public key. */
    unpack_public_key(M, sk + SEED_SIZE);
}

void pack_signature(
    uint8_t *sig,
    size_t *sig_len,
    const hash_t salt,
    const hash_t hash1,
    const hash_t hash2,
    const uint32_t i_star[TAU],
    seed_t tree[TAU][TREE_N_NODES],
    hash_t com[TAU][N_PARTIES],
    ff_t a_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t C_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t S_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)])
{
    /* The signature is packed as follows:
     *
     *      pack(salt)
     *      pack(hash1)
     *      pack(hash2)
     *
     *      for l in {0, 1, ..., TAU}
     *          pack(com[l][i_star[l]])
     *          pack(packed_tree[l])
     *
     *      for l in {0, 1, ..., TAU}
     *          if i_star[l] != N_PARTIES - 1
     *              pack(a_rnd_shr[l][N_PARTIES - 1])
     *              pack(K_rnd_shr[l][N_PARTIES - 1])
     *              pack(C_rnd_shr[l][N_PARTIES - 1])
     *          pack(S_rnd_shr[l][i_star[l]])
     */

    uint8_t *ptr;
    uint32_t l, bo;
    seed_t packed_tree[TAU][TREE_HEIGHT];

    ptr = sig;
    bo = 0;

    /* Pack 'salt'. */
    memcpy(ptr, salt, HASH_SIZE);
    ptr += HASH_SIZE;

    /* Pack 'hash1'. */
    memcpy(ptr, hash1, HASH_SIZE);
    ptr += HASH_SIZE;
    
    /* Pack 'hash2'. */
    memcpy(ptr, hash2, HASH_SIZE);
    ptr += HASH_SIZE;

    for (l = 0; l < TAU; l++)
    {
        /* Pack 'com[l][i_star[l]]'. */
        memcpy(ptr, com[l][i_star[l]], HASH_SIZE); //
        ptr += HASH_SIZE;

        /* Compress the seed tree. */
        seed_tree_pack(packed_tree[l], tree[l], i_star[l]);

        /* Pack 'packed_tree[l]'. */
        memcpy(ptr, packed_tree[l], TREE_HEIGHT * SEED_SIZE);
        ptr += TREE_HEIGHT * SEED_SIZE;
    }

    /* NOTE: Matrices are packed last because their sizes are not
     * always exact multiples of 8 bits, and so one has to take care
     * of bit offsets. */

    for (l = 0; l < TAU; l++)
    {
        if (i_star[l] != N_PARTIES - 1)
        {
            /* Pack 'a_rnd_shr[i]', 'K_rnd_shr[i]', 'C_rnd_shr[i]'. */
            matrix_pack(&ptr, &bo, a_rnd_shr[l][N_PARTIES - 1], PAR_K, 1u);
            matrix_pack(&ptr, &bo, K_rnd_shr[l][N_PARTIES - 1], PAR_R, PAR_N - PAR_R);
            matrix_pack(&ptr, &bo, C_rnd_shr[l][N_PARTIES - 1], PAR_S, PAR_N - PAR_R);
        }

        /* Pack 'S_rnd_shr[l][i_star[l]]'. */
        matrix_pack(&ptr, &bo, S_rnd_shr[l][i_star[l]], PAR_S, PAR_R);
    }

    /* Compute the signature length. */
    if (bo == 0)
    {
        *sig_len = ptr - sig;
    }
    else
    {
        *sig_len = ptr + 1 - sig;
    }
}

int unpack_signature(
    hash_t salt,
    hash_t hash1,
    hash_t hash2,
    uint32_t i_star[TAU],
    seed_t packed_tree[TAU][TREE_HEIGHT],
    hash_t com_star[TAU],
    ff_t a_rnd_aux[TAU][matrix_bytes_size(PAR_K, 1)],
    ff_t K_rnd_aux[TAU][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t C_rnd_aux[TAU][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t S_rnd_star[TAU][matrix_bytes_size(PAR_S, PAR_R)],
    const uint8_t *sig,
    size_t *sig_len)
{
    uint8_t *ptr;
    uint32_t l, bo;
    
    ptr = (uint8_t *)sig;
    bo = 0;

    /* Unpack 'salt'. */
    memcpy(salt, ptr, HASH_SIZE);
    ptr += HASH_SIZE;

    /* Unpack 'hash1'. */
    memcpy(hash1, ptr, HASH_SIZE);
    ptr += HASH_SIZE;

    /* Unpack 'hash2'. */
    memcpy(hash2, ptr, HASH_SIZE);
    ptr += HASH_SIZE;

    /* Compute the challenges 'i_star'. */
    get_second_challenges(i_star, hash2);

    for (l = 0; l < TAU; l++)
    {
        /* Unpack 'com[l][i_star[l]]'. */
        memcpy(com_star[l], ptr, HASH_SIZE);
        ptr += HASH_SIZE;

        /* Unpack 'packed_tree[l]'. */
        memcpy(packed_tree[l], ptr, TREE_HEIGHT * SEED_SIZE);
        ptr += TREE_HEIGHT * SEED_SIZE;
    }

    for (l = 0; l < TAU; l++)
    {
        if (i_star[l] != N_PARTIES - 1)
        {
            /* Unpack 'a_rnd_aux[l]', 'K_rnd_aux[l]', 'C_rnd_aux[l]'. */
            matrix_unpack(a_rnd_aux[l], &ptr, &bo, PAR_K, 1);
            matrix_unpack(K_rnd_aux[l], &ptr, &bo, PAR_R, PAR_N - PAR_R);
            matrix_unpack(C_rnd_aux[l], &ptr, &bo, PAR_S, PAR_N - PAR_R);
        }

        /* Unpack 'S_rnd_star[l]'. */
        matrix_unpack(S_rnd_star[l], &ptr, &bo, PAR_S, PAR_R);
    }

    /* NOTE: If the number of packed matrices is odd, then the last byte
     * must have the 4 most significant bits all equal to zero. */ 
    if ((bo == 1) && ((*ptr) & 0xF0))
    {
        return 0;
    }

    /* Compute the signature length. */
    if (bo == 0)
    {
        *sig_len = ptr - sig;
    }
    else
    {
        *sig_len = ptr + 1 - sig;
    }

    return 1;
}
