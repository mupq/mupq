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

#include <stdlib.h>
#include <string.h>
#include "challenges.h"
#include "hash.h"
#include "matrix.h"
#include "packing.h"
#include "params.h"
#include "prng.h"
#include "seed_tree.h"
#include "sign.h"

#include <randombytes.h>

/* A round of Phase 1 of signing.
 *
 * Given as input:
 * - the salt 'salt';
 * - the number 'l' of the round;
 * - the seeds 'seed' of the round;
 * - the secret key 'a', 'K';
 * compute the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr',
 * the matrix 'A', and the commitments 'com'. */
void sign_phase1_round(

    ff_t A[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_shr[N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_shr[N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    hash_t com[N_PARTIES],

    hash_t salt,
    uint32_t l,
    seed_t seed[N_PARTIES],
    const ff_t a[matrix_bytes_size(PAR_K, 1)],
    const ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)])
{
    uint32_t i;

    for (i = 0; i < N_PARTIES; i++)
    {
        prng_t prng;

        /* Initialize PRNG from 'seed[i]'. */
        prng_init(&prng, salt, seed[i]);

        /* Generate the random matrix 'A_shr[i]'. */
        matrix_init_random(A_shr[i], PAR_S, PAR_R, &prng);

        if (i != N_PARTIES - 1)
        {
            /* Generate random matrices 'a_shr[i]', 'C_shr[i]', 'K_shr[i]'. */
            matrix_init_random(a_shr[i], PAR_K, 1, &prng);
            matrix_init_random(C_shr[i], PAR_S, PAR_N - PAR_R, &prng);
            matrix_init_random(K_shr[i], PAR_R, PAR_N - PAR_R, &prng);

            /* Compute the commitment. */
            hash_digest0(com[i], salt, l, i, seed[i]);
        }
        else
        {
            uint32_t j;

            /* Set a_shr[N_PARTIES - 1] = a - sum_{j < N_PARTIES - 1} a_shr[j]. */
            matrix_copy(a_shr[N_PARTIES - 1], a, PAR_K, 1);

            for (j = 0; j < N_PARTIES - 1; j++)
            {
                matrix_subtract(a_shr[N_PARTIES - 1], a_shr[j], PAR_K, 1);
            }
            /* * */

            /* Set K_shr[N_PARTIES - 1] = K - sum_{j < N_PARTIES - 1} K_shr[j]. */
            matrix_copy(K_shr[N_PARTIES - 1], K, PAR_R, PAR_N - PAR_R);

            for (j = 0; j < N_PARTIES - 1; j++)
            {
                matrix_subtract(K_shr[N_PARTIES - 1], K_shr[j], PAR_R, PAR_N - PAR_R);
            }
            /* * */

            /* Open 'A'. */
            matrix_init_zero(A, PAR_S, PAR_R);

            for (j = 0; j < N_PARTIES; j++)
            {
                matrix_add(A, A_shr[j], PAR_S, PAR_R);
            }
            /* * */

            /* Set C[N_PARTIES - 1] = A * K - sum_{j < N_PARTIES - 1} C[j]. */
            matrix_product(C_shr[N_PARTIES - 1], A, K, PAR_S, PAR_R, PAR_N - PAR_R);

            for (j = 0; j < N_PARTIES - 1; j++)
            {
                matrix_subtract(C_shr[N_PARTIES - 1], C_shr[j], PAR_S, PAR_N - PAR_R);
            }
            /* * */

            /* Compute the commitment. */
            hash_digest0_aux(com[i], salt, l, i, seed[N_PARTIES - 1],
                a_shr[N_PARTIES - 1], K_shr[N_PARTIES - 1], C_shr[N_PARTIES - 1]);
        }
    }
}

/* A round of Phase 3 of signing.
 *
 * Given as input:
 * - the hashing context 'hash_ctx';
 * - the public key 'M';
 * - the secret key 'E';
 * - the matrix 'A';
 * - the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr';
 * - the challenge matrix 'R';
 * computes and hashes the additive shares 'S_shr', 'V_shr'.
 * The matrix 'S_shr[i]' is written over 'A_shr'. */
void sign_phase3_round
(
    hash_ctx_t *hash_ctx,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)],
    ff_t A[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_shr[N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_shr[N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t R[matrix_bytes_size(PAR_S, PAR_M)]
)
{
    uint32_t i;

    ff_t S[matrix_bytes_size(PAR_S, PAR_R)];

    /* Open S. */
    {
        ff_t MaL[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];
        ff_t MaR[matrix_bytes_size(PAR_M, PAR_R)];

        matrix_horizontal_split(MaL, MaR, E, PAR_M, PAR_N - PAR_R, PAR_R);

        /* 'S = A - R * MaR'.
         * (or equivalently 'S = R * MaR + A' because we are in characteristic 2. */
        matrix_product(S, R, MaR, PAR_S, PAR_M, PAR_R);
        matrix_add(S, A, PAR_S, PAR_R);
    }
    /* * */

    for (i = 0; i < N_PARTIES; i++)
    {
        uint32_t j;

        ff_t MaL_shr_i[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];
        ff_t MaR_shr_i[matrix_bytes_size(PAR_M, PAR_R)];
        ff_t E_shr_i[matrix_bytes_size(PAR_M, PAR_N)];
        ff_t RxMaL_shr_i[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
        ff_t V_shr_i[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
        ff_t RxMaR_shr_i[matrix_bytes_size(PAR_S, PAR_R)];

        /* Set E_shr_i = (M[0] if i = 0 else 0) + sum_{j = 0}^{PAR_K-1} a_shr[i] * M[j + 1]. */
        if (i == 0)
        {
            matrix_copy(E_shr_i, M[0], PAR_M, PAR_N);
        }
        else
        {
            matrix_init_zero(E_shr_i, PAR_M, PAR_N);
        }

        for (j = 0; j < PAR_K; j++)
        {
            ff_t scalar;

            scalar = matrix_get_entry(a_shr[i], PAR_K, j, 0);
            matrix_add_multiple(E_shr_i, scalar, M[j + 1], PAR_M, PAR_N);
        }
        /* * */

        /* Compute 'MaL_shr_i' and 'MaR_shr_i' using the fact
         * that E_shr_i = [MaL_i | -MaR_i]. */
        matrix_horizontal_split(MaL_shr_i, MaR_shr_i, E_shr_i, PAR_M, PAR_N - PAR_R, PAR_R);
        matrix_negate(MaR_shr_i, PAR_M, PAR_R);

        /* Compute 'R * MaL_shr_i'. */
        matrix_product(RxMaL_shr_i, R, MaL_shr_i, PAR_S, PAR_M, PAR_N - PAR_R);

        /* Set V_shr_i = S * K_shr[i] - R * MaL_shr_i - C_shr[i]. */
        matrix_product(V_shr_i, S, K_shr[i], PAR_S, PAR_R, PAR_N - PAR_R);
        matrix_subtract(V_shr_i, RxMaL_shr_i, PAR_S, PAR_N - PAR_R);
        matrix_subtract(V_shr_i, C_shr[i], PAR_S, PAR_N - PAR_R);

        /* Overwrite 'A_shr[i]' with 'S_shr[i] = R * MaR_shr_i + A_shr[i]'. */
        matrix_product(RxMaR_shr_i, R, MaR_shr_i, PAR_S, PAR_M, PAR_R);
        matrix_add(A_shr[i], RxMaR_shr_i, PAR_S, PAR_R);

        /* Hash 'S_shr[i]', 'V_shr_i' (NOTE: 'A_shr[i]' has been overwritten with 'S_shr[i]'). */
        hash_update(hash_ctx, A_shr[i], matrix_bytes_size(PAR_S, PAR_R));
        hash_update(hash_ctx, V_shr_i, matrix_bytes_size(PAR_S, PAR_N - PAR_R));
    }
}

/* A round of Phase 1 of verifying.
 *
 * Given as input:
 * - the hashing context 'hash1_ctx';
 * - the salt 'salt';
 * - the number 'l' of the round;
 * - the seeds 'seed' of the round;
 * - the party of the provided commitment 'i_star';
 * - the provided commitment 'com_star';
 * - the auxiliary matrices 'C_aux', 'a_aux', 'K_aux';
 * compute the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr'
 * while hashing their commitments.
 */
void open_phase1_round(

    hash_ctx_t *hash1_ctx,
    hash_t salt,
    uint32_t l,
    seed_t seed[N_PARTIES],
    uint32_t i_star,
    hash_t com_star,
    ff_t C_aux[matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_aux[matrix_bytes_size(PAR_K, 1)],
    ff_t K_aux[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],

    ff_t A_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_shr[N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_shr[N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)])
{
    uint32_t i;

    for (i = 0; i < N_PARTIES; i++)
    {
        hash_t com;

        if (i == i_star)
        {
            /* Copy the commitment from the signature. */
            hash_update(hash1_ctx, com_star, HASH_SIZE);
        }
        else
        {
            prng_t prng;

            /* Initialize PRNG from 'seed[i]'. */
            prng_init(&prng, salt, seed[i]);

            /* Generate the random matrix 'A_shr[i]'. */
            matrix_init_random(A_shr[i], PAR_S, PAR_R, &prng);

            if (i != N_PARTIES - 1)
            {
                /* Generate random matrices 'a_shr[i]', 'C_shr[i]', 'K_shr[i]'. */
                matrix_init_random(a_shr[i], PAR_K, 1, &prng);
                matrix_init_random(C_shr[i], PAR_S, PAR_N - PAR_R, &prng);
                matrix_init_random(K_shr[i], PAR_R, PAR_N - PAR_R, &prng);

                /* Compute the commitment. */
                hash_digest0(com, salt, l, i, seed[i]);
            }
            else
            {
                /* Compute the commitment. */
                hash_digest0_aux(com, salt, l, i, seed[N_PARTIES - 1],
                    a_aux, K_aux, C_aux);

                /* Copy the matrices from the signature. */
                matrix_copy(a_shr[N_PARTIES - 1], a_aux, PAR_K, 1);
                matrix_copy(K_shr[N_PARTIES - 1], K_aux, PAR_R, PAR_N - PAR_R);
                matrix_copy(C_shr[N_PARTIES - 1], C_aux, PAR_S, PAR_N - PAR_R);
            }
            /* Update the hashing context. */
            hash_update(hash1_ctx, com, HASH_SIZE);
        }
    }
}

/* A round of Phase 3 of verifying.
 *
 * Given as input:
 * - the hashing context 'hash_ctx';
 * - the party id 'i_star';
 * - the public key 'M';
 * - the additive share 'S_star' of party 'i_star';
 * - the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr';
 * compute the additive shares 'V_shr' and hash them. */
void open_phase3_round
(
    hash_ctx_t *hash_ctx,
    uint32_t i_star,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t R[matrix_bytes_size(PAR_S, PAR_M)],
    ff_t S_star[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_shr[N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_shr[N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)]
)
{
    uint32_t i;

    ff_t MaL_shr[N_PARTIES][matrix_bytes_size(PAR_M, PAR_N - PAR_R)];
    ff_t S[matrix_bytes_size(PAR_S, PAR_R)];

    ff_t V_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t S_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];

    for (i = 0; i < N_PARTIES; i++)
    {
        uint32_t j;

        ff_t E_shr_i[matrix_bytes_size(PAR_M, PAR_N)];
        ff_t MaR_shr_i[matrix_bytes_size(PAR_M, PAR_R)];

        /* Skip the 'i_star' party. */
        if (i == i_star)
        {
            /* Copy the matrix provided by the signature. */
            matrix_copy(S_shr[i], S_star, PAR_S, PAR_R);
            continue;
        }

        /* Set E_shr_i = (M[0] if i = 0 else 0) + sum_{j = 0}^{PAR_K-1} a_shr[i] * M[j + 1]. */
        if (i == 0)
        {
            matrix_copy(E_shr_i, M[0], PAR_M, PAR_N);
        }
        else
        {
            matrix_init_zero(E_shr_i, PAR_M, PAR_N);
        }

        for (j = 0; j < PAR_K; j++)
        {
            ff_t scalar;

            scalar = matrix_get_entry(a_shr[i], PAR_K, j, 0);
            matrix_add_multiple(E_shr_i, scalar, M[j + 1], PAR_M, PAR_N);
        }
        /* * */

        /* Compute 'MaL_shr[i]' and 'MaR_shr_i' using the fact
         * that E_shr_i = [MaL_shr[i] | -MaR_shr_i]. */
        matrix_horizontal_split(MaL_shr[i], MaR_shr_i, E_shr_i, PAR_M, PAR_N - PAR_R, PAR_R);
        matrix_negate(MaR_shr_i, PAR_M, PAR_R);

        /* Set 'S_shr[i] = R * MaR_shr_i + A[i]'. */
        matrix_product(S_shr[i], R, MaR_shr_i, PAR_S, PAR_M, PAR_R);
        matrix_add(S_shr[i], A_shr[i], PAR_S, PAR_R);
    }

    /* Open S. */
    matrix_init_zero(S, PAR_S, PAR_R);

    for (i = 0; i < N_PARTIES; i++)
    {
        matrix_add(S, S_shr[i], PAR_S, PAR_R);
    }
    /* * */

    for (i = 0; i < N_PARTIES; i++)
    {
        ff_t RxMaL_shr_i[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];

        /* Skip the 'i_star' party. */
        if (i == i_star)
        {
            continue;
        }

        /* Compute 'R * MaL_shr[i]'. */
        matrix_product(RxMaL_shr_i, R, MaL_shr[i], PAR_S, PAR_M, PAR_N - PAR_R);

        /* Set V_shr[i] = S * K_shr[i] - R * MaL_shr_i - C_shr[i]. */
        matrix_product(V_shr[i], S, K_shr[i], PAR_S, PAR_R, PAR_N - PAR_R);
        matrix_subtract(V_shr[i], RxMaL_shr_i, PAR_S, PAR_N - PAR_R);
        matrix_subtract(V_shr[i], C_shr[i], PAR_S, PAR_N - PAR_R);
    }

    /* Set V_shr[i_star] = -sum_{i != i_star} V[i]. */
    matrix_init_zero(V_shr[i_star], PAR_S, PAR_N - PAR_R);

    for (i = 0; i < N_PARTIES; i++)
    {
        if (i != i_star)
        {
            matrix_subtract(V_shr[i_star], V_shr[i], PAR_S, PAR_N - PAR_R);
        }
    }
    /* * */

    /* Hash 'S_shr[i]', 'V_shr[i]'. */
    for (i = 0; i < N_PARTIES; i++)
    {
        hash_update(hash_ctx, S_shr[i], matrix_bytes_size(PAR_S, PAR_R));
        hash_update(hash_ctx, V_shr[i], matrix_bytes_size(PAR_S, PAR_N - PAR_R));
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
    /* PRNG. */
    prng_t prng;

    /* Public key. */
    seed_t seed_pk;
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)];

    /* Secret key. */
    seed_t seed_sk;
    ff_t a[matrix_bytes_size(PAR_K, 1)];
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)];
    ff_t E_R[matrix_bytes_size(PAR_M, PAR_R)];

    /* Other variables. */
    uint32_t i;
    ff_t T[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];

    /* Generate random seed of the secret key. */
    randombytes(seed_sk, SEED_SIZE);

    /* Initialize the PRNG. */
    prng_init(&prng, NULL, seed_sk);

    /* Initialize matrices 'a', 'K', 'E_R'. */
    matrix_init_random(a, PAR_K, 1, &prng);
    matrix_init_random(K, PAR_R, PAR_N - PAR_R, &prng);
    matrix_init_random(E_R, PAR_M, PAR_R, &prng);

    /* Compute the matrix 'E'. */
    matrix_product(T, E_R, K, PAR_M, PAR_R, PAR_N - PAR_R);
    matrix_horizontal_concatenation(E, T, E_R, PAR_M, PAR_N - PAR_R, PAR_R);

    /* Generate random seed of the public key. */
    randombytes(seed_pk, SEED_SIZE);

    /* Initialize matrices 'M[i]' (i = 1, ..., k). */
    prng_init(&prng, NULL, seed_pk);

    for (i = 1; i <= PAR_K; i++)
    {
        matrix_init_random(M[i], PAR_M, PAR_N, &prng);
    }
    /* * */

    /* Compute matrix 'M[0]'. */
    matrix_copy(M[0], E, PAR_M, PAR_N);

    for (i = 0; i < PAR_K; i++)
    {
        ff_t scalar;

        scalar = matrix_get_entry(a, PAR_K, i, 0);
        matrix_subtract_multiple(M[0], scalar, M[i + 1], PAR_M, PAR_N);
    }
    /* * */

    /* Pack the public key. */
    pack_public_key(pk, M[0], seed_pk);

    /* Pack the secret key. */
    pack_secret_key(sk, M[0], seed_pk, seed_sk);

    /* Success. */
    return 0;
}

int crypto_sign(uint8_t *sig_msg, size_t *sig_msg_len,
    uint8_t *msg, size_t msg_len,
    uint8_t *sk)
{
    /* Public key. */
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)];

    /* Secret key. */
    ff_t a[matrix_bytes_size(PAR_K, 1)];
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)];

    /* Unpack the secret key (and the public key). */
    unpack_secret_key(M, a, K, E, sk);

	/* copy the message to the end of the signature buffer.
	 * This allows for cases where 'sig' = 'msg'; */
    memcpy(sig_msg + CRYPTO_BYTES, msg, msg_len);

    /* Sign the message. */
    return crypto_sign_unpacked_keys(sig_msg, sig_msg_len,
        msg, msg_len, M, a, K, E);
}

int crypto_sign_unpacked_keys(uint8_t *sig_msg, size_t *sig_msg_len,
    uint8_t *msg, size_t msg_len,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t a[matrix_bytes_size(PAR_K, 1)],
    ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)])
{
    /* Hashes. */
    hash_t salt;
    hash_t hash1;
    hash_t hash2;
    hash_t com[TAU][N_PARTIES];

    /* Hashing context. */
    hash_ctx_t hash_ctx;

    /* Seed trees. */
    seed_t trees[TAU][TREE_N_NODES];

    /* Matrices. */
    ff_t A_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];
    ff_t C_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t a_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_K, 1)];
    ff_t K_rnd_shr[TAU][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
    ff_t A_rnd[TAU][matrix_bytes_size(PAR_S, PAR_R)];

    /* PRNG to generate the first challenges. */
    prng_t prng;

    /* Second challenges. */
    uint32_t i_star[TAU];

    /* Round counter. */
    uint32_t l;

    /* Signature length. */
    size_t sig_len;

    /* Initialize salt. */
    randombytes(salt, HASH_SIZE);

    /* Initialize the hashing context to compute 'hash1'. */
    hash_init(&hash_ctx);
    hash_update(&hash_ctx, salt, HASH_SIZE);
    hash_update(&hash_ctx, msg, msg_len);

    for (l = 0; l < TAU; l++)
    {
        seed_t seed;

        /* Generate a random seed. */
        randombytes(seed, SEED_SIZE);

        /* Generate a tree from the random seed. */
        seed_tree_init(trees[l], salt, seed);

        /* Compute the matrices
         * 'A_rnd[l][0...N_PARTIES-1]',
         * 'A_rnd_shr[l][0...N_PARTIES-1]',
         * 'C_rnd_shr[l][0...N_PARTIES-1]',
         * 'a_rnd_shr[l][0...N_PARTIES-1]',
         * 'K_rnd_shr[l][0...N_PARTIES-1]',
         * and the commitments 'com[l][0...N_PARTIES-1]'. */
        sign_phase1_round(A_rnd[l], A_rnd_shr[l], C_rnd_shr[l], a_rnd_shr[l], K_rnd_shr[l], com[l],
            salt, l, (uint8_t (*)[SEED_SIZE]) seed_tree_get_leaves(trees[l]), a, K);

        /* Hash the commitments for this round. */
        hash_update(&hash_ctx, (uint8_t *)com[l], N_PARTIES * HASH_SIZE);
    }

    /* Compute 'hash1'. */
    hash_finalize(&hash_ctx, hash1);

    /* Initialize the PRNG to generate the first challenges. */
    prng_init(&prng, hash1, NULL);

    /* Initialize the hashing context to compute 'hash2'. */
    hash_init(&hash_ctx);
    hash_update(&hash_ctx, salt, HASH_SIZE);
    hash_update(&hash_ctx, msg, msg_len);

    for (l = 0; l < TAU; l++)
    {
        ff_t R[matrix_bytes_size(PAR_S, PAR_M)];

        /* Generate the first challenge 'R' for this round. */
        matrix_init_random(R, PAR_S, PAR_M, &prng);

        /* Compute and hash the matrices
         * 'S_rnd_shr[l][0...N_PARTIES-1]' and 'V_rnd_shr[l][0...N_PARTIES-1]'. */
        sign_phase3_round(&hash_ctx, M, E, A_rnd[l], A_rnd_shr[l], C_rnd_shr[l], a_rnd_shr[l], K_rnd_shr[l], R);
    }

    /* Compute 'hash2'. */
    hash_update(&hash_ctx, hash1, HASH_SIZE);
    hash_finalize(&hash_ctx, hash2);

    /* Generate the second challenges. */
    get_second_challenges(i_star, hash2);

    /* Pack the signature. */
    pack_signature(sig_msg, &sig_len, salt, hash1, hash2, i_star, trees, com, a_rnd_shr, K_rnd_shr, C_rnd_shr, A_rnd_shr);

    /* Append the message to the signature. */
    memcpy(sig_msg + sig_len, sig_msg + CRYPTO_BYTES, msg_len);

    /* Update 'sign_msg_len'. */
    if (sig_msg_len != NULL)
    {
        *sig_msg_len = sig_len + msg_len;
    }

    /* Success. */
    return 0;
}

int crypto_sign_open(uint8_t *msg, size_t *msg_len,
    uint8_t *sig_msg, size_t sig_msg_len,
    uint8_t *pk)
{
    /* Public key. */
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)];

    /* Hashes. */
    hash_t salt;
    hash_t hash1;
    hash_t hash2;
    hash_t hash1_c;
    hash_t hash2_c;
    hash_t com_star[TAU];

    /* Hashing contexts. */
    hash_ctx_t hash1_ctx;
    hash_ctx_t hash2_ctx;

    /* Seeds. */
    seed_t packed_tree[TAU][TREE_HEIGHT];

    /* Matrices. */
    ff_t C_rnd_aux[TAU][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t a_rnd_aux[TAU][matrix_bytes_size(PAR_K, 1)];
    ff_t K_rnd_aux[TAU][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
    ff_t S_rnd_star[TAU][matrix_bytes_size(PAR_S, PAR_R)];

    /* Second challenges. */
    uint32_t i_star[TAU];

    /* PRNG. */
    prng_t prng;

    /* Other variables. */
    uint32_t l;

    /* Signature length. */
    size_t sig_len;

    /* Unpack the public key. */
    unpack_public_key(M, pk);

    /* Unpack the signature. */
    if (!unpack_signature(salt, hash1, hash2, i_star, packed_tree, com_star,
        a_rnd_aux, K_rnd_aux, C_rnd_aux, S_rnd_star, sig_msg, &sig_len))
    {
        /* Failure. */
        return -1;
    }

    /* Compute the message length. */
    *msg_len = sig_msg_len - sig_len;

    /* Initialize the PRNG to compute the first challenges. */
    prng_init(&prng, hash1, NULL);

    /* NOTE: The second challenges 'i_star' were already computed when
     * the signature was unpacked. */

    /* Initialize the hashing context to compute 'hash1_c'. */
    hash_init(&hash1_ctx);
    hash_update(&hash1_ctx, salt, HASH_SIZE);
    hash_update(&hash1_ctx, sig_msg + sig_len, *msg_len);

    /* Initialize the hashing context to compute 'hash2_c'. */
    hash_init(&hash2_ctx);
    hash_update(&hash2_ctx, salt, HASH_SIZE);
    hash_update(&hash2_ctx, sig_msg + sig_len, *msg_len);

    for (l = 0; l < TAU; l++)
    {
        /* Matrices of this round. */
        ff_t A_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];
        ff_t C_shr[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
        ff_t a_shr[N_PARTIES][matrix_bytes_size(PAR_K, 1)];
        ff_t K_shr[N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
        ff_t R[matrix_bytes_size(PAR_S, PAR_M)];

        /* Seeds of this round. */
        seed_t seed[N_PARTIES];

        /* Decompress the seed tree. */
        seed_tree_unpack(seed, salt, packed_tree[l], i_star[l]);

        /* Compute the matrices
         * 'A_shr[0...N_PARTIES-1]',
         * 'C_shr[0...N_PARTIES-1]',
         * 'a_shr[0...N_PARTIES-1]',
         * 'K_shr[0...N_PARTIES-1]',
         * while hashing the commitments. */
        open_phase1_round(&hash1_ctx, salt, l, seed, i_star[l], com_star[l],
            C_rnd_aux[l], a_rnd_aux[l], K_rnd_aux[l],
            A_shr, C_shr, a_shr, K_shr);

        /* Compute the first challenge of this round. */
        matrix_init_random(R, PAR_S, PAR_M, &prng);

        /* Compute and hash the matrices 'S1[0...N_PARTIES-1]', 'S2[0...N_PARTIES-1]', 'V[0...N_PARTIES-1]'. */
        open_phase3_round(&hash2_ctx, i_star[l], M, R, S_rnd_star[l], A_shr, C_shr, a_shr, K_shr);
    }

    /* Compute 'hash1_c'. */
    hash_finalize(&hash1_ctx, hash1_c);

    /* Compute 'hash2_c'. */
    hash_update(&hash2_ctx, hash1_c, HASH_SIZE);
    hash_finalize(&hash2_ctx, hash2_c);

    /* Check hashes. */
    if (hash_equal(hash1, hash1_c) && hash_equal(hash2, hash2_c))
    {
        /* Success. */
        memcpy(msg, sig_msg + sig_len, *msg_len);

        return 0;
    }

    /* Failure. */
    return -1;
}

