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

#if defined (OFFLINE_CC)
#include "internal_get_cycles.h"
#endif

/* Set E_shr_i = (M[0] if i = 0 else 0) + sum_{p = 0}^{PAR_K-1} a_shr[p] * M[p + 1]. */
void computeE_shr(ff_t E_shr[matrix_bytes_size(PAR_M, PAR_N)],
                  ff_t a_shr[matrix_bytes_size(PAR_K, 1)],
                  ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)])
{
    uint32_t p;
    matrix_init_zero(E_shr, PAR_M, PAR_N);

    for (p = 0; p < PAR_K; p++)
    {
        ff_t scalar;

        scalar = matrix_get_entry(a_shr, PAR_K, p, 0);
        matrix_add_multiple(E_shr, scalar, M[p + 1], PAR_M, PAR_N);
    }
}


/* A round of Phase 1 of signing.
 *
 * Given as input:
 * - the salt 'salt';
 * - the number 'l' of the round;
 * - the seeds 'seed' of the round;
 * - the secret key 'a', 'K';
 * compute the matrix 'A', the additive shares of the main parties
 * 'A_main_shr', 'C_main_shr', 'a_main_shr', 'K_main_shr', the matrices
 * of the last parties 'C_aux', 'a_aux', 'K_aux', and the
 * commitment 'com_l'. */
void sign_phase1_round(
    ff_t A[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t C_aux[matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_aux[matrix_bytes_size(PAR_K, 1)],
    ff_t K_aux[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    hash_t com_l,
    hash_t salt,
    uint32_t l,
    seed_t seed[N_PARTIES_ROUND],
    const ff_t a[matrix_bytes_size(PAR_K, 1)],
    const ff_t K[matrix_bytes_size(PAR_R, PAR_N - PAR_R)])
{
    uint32_t i, j, k;

    /* Hashing context. */
    hash_ctx_t hash_ctx;

    ff_t A_shr_i[matrix_bytes_size(PAR_S, PAR_R)];
    ff_t C_shr_i[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t a_shr_i[matrix_bytes_size(PAR_K, 1)];
    ff_t K_shr_i[matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
    hash_t com_i_l;

    /* Set 'a_aux = a'. */
    matrix_copy(a_aux, a, PAR_K, 1);

    /* Set 'K_aux = K'. */
    matrix_copy(K_aux, K, PAR_R, PAR_N - PAR_R);

    /* Initialize main parties shares to zero.*/
    for (i = 0; i < D_DIMENSION; i++)
    {
        for (j = 0; j < N_PARTIES; j++)
        {
            matrix_init_zero(A_main_shr[i][j], PAR_S, PAR_R);
            matrix_init_zero(C_main_shr[i][j], PAR_S, PAR_N - PAR_R);
            matrix_init_zero(a_main_shr[i][j], PAR_K, 1);
            matrix_init_zero(K_main_shr[i][j], PAR_R, PAR_N - PAR_R);
        }
    }
    /* * */

    /* Initialize the hashing context to compute 'com'. */
    hash_init(&hash_ctx);
    hash_update(&hash_ctx, salt, HASH_SIZE);
    hash_update(&hash_ctx, (uint8_t *)&l, sizeof(l));

#if defined (OFFLINE_CC)
begin_offline = get_cycles();
#endif

    /* Initialize 'A' and 'C_aux' to zero. */
    matrix_init_zero(A, PAR_S, PAR_R);
    matrix_init_zero(C_aux, PAR_S, PAR_N - PAR_R);

    for (i = 0; i < N_PARTIES_ROUND; i++)
    {
        prng_t prng;

        /* Initialize PRNG from 'seed[i]'. */
        prng_init(&prng, salt, seed[i]);

        /* Generate the random matrix 'A_shr_i'. */
        matrix_init_random(A_shr_i, PAR_S, PAR_R, &prng);

        /* Add 'A_shr_i' to 'A'. */
        matrix_add(A, A_shr_i, PAR_S, PAR_R);

        if (i != N_PARTIES_ROUND - 1)
        {
            /* Generate random matrices 'a_shr_i', 'C_shr_i', 'K_shr_i'. */
            matrix_init_random(a_shr_i, PAR_K, 1, &prng);
            matrix_init_random(C_shr_i, PAR_S, PAR_N - PAR_R, &prng);
            matrix_init_random(K_shr_i, PAR_R, PAR_N - PAR_R, &prng);

            /* Compute the commitment. */
            hash_digest0(com_i_l, salt, l, i, seed[i]);

            hash_update(&hash_ctx, com_i_l, HASH_SIZE);

            /* Substract 'C_shr_i' from 'C_aux'. In this way, before the last
             * party it would be 'C_aux = -sum_{i < N_PARTIES_ROUND - 1} C_shr_i'. */
            matrix_subtract(C_aux, C_shr_i, PAR_S, PAR_N - PAR_R);

            /* Substract 'a_shr_i' from 'a_aux'. In this way, at the end it would be
             * 'a_aux = a - sum_{i < N_PARTIES_ROUND - 1} a_shr_i'. */
            matrix_subtract(a_aux, a_shr_i, PAR_K, 1);

            /* Substract 'K_shr_i' from 'K_aux'. In this way, at the end it would be
             * 'K_aux = K - sum_{i < N_PARTIES_ROUND - 1} K_shr_i'. */
            matrix_subtract(K_aux, K_shr_i, PAR_R, PAR_N - PAR_R);

            /* Aggregate to the main parties. */
            for (k = 0; k < D_DIMENSION; k++)
            {
                uint32_t i_star_k;

                i_star_k = (i >> k) & 1u;

                matrix_add(A_main_shr[k][i_star_k], A_shr_i, PAR_S, PAR_R);
                matrix_add(C_main_shr[k][i_star_k], C_shr_i, PAR_S, PAR_N - PAR_R);
                matrix_add(a_main_shr[k][i_star_k], a_shr_i, PAR_K, 1);
                matrix_add(K_main_shr[k][i_star_k], K_shr_i, PAR_R, PAR_N - PAR_R);
            }
            /* * */
        }
    }

#if defined (OFFLINE_CC)
offline_cc = offline_cc + (get_cycles() - begin_offline);
#endif

    ff_t AxK[matrix_bytes_size(PAR_S, PAR_N - PAR_S)];

    /* Compute the product of 'A' and 'K'. */
    matrix_product(AxK, A, K, PAR_S, PAR_R, PAR_N - PAR_R);

    /* Add 'AxK' to 'C_aux', so that
        * 'C_aux = A * K - sum_{j < N_PARTIES_ROUND - 1} C_shr_i'. */
    matrix_add(C_aux, AxK, PAR_S, PAR_N - PAR_R);

    /* Compute the commitment. */
    hash_digest0_aux(com_i_l, salt, l, N_PARTIES_ROUND - 1, seed[N_PARTIES_ROUND - 1],
        a_aux, K_aux, C_aux);

    hash_update(&hash_ctx, com_i_l, HASH_SIZE);

    /* Aggregate to the main parties. */
    for (k = 0; k < D_DIMENSION; k++)
    {
        uint32_t i_star_k;

        i_star_k = ((N_PARTIES_ROUND - 1) >> k) & 1;

        matrix_add(A_main_shr[k][i_star_k], A_shr_i, PAR_S, PAR_R);
        matrix_add(C_main_shr[k][i_star_k], C_aux, PAR_S, PAR_N - PAR_R);
        matrix_add(a_main_shr[k][i_star_k], a_aux, PAR_K, 1);
        matrix_add(K_main_shr[k][i_star_k], K_aux, PAR_R, PAR_N - PAR_R);
    }
    /* * */

    /* Computing 'com'. */
    hash_finalize(&hash_ctx, com_l);
}

/* A round of Phase 3 of signing.
 *
 * Given as input:
 * - the hashing context 'hash_ctx';
 * - the public key 'M';
 * - the secret key 'E';
 * - the matrix 'A';
 * - the additive shares 'A_main_shr', 'C_main_shr', 'a_main_shr', 'K_main_shr';
 * - the challenge matrix 'R';
 * - the salt 'salt';
 * - the round number 'l';
 * compute and hash the additive shares 'S_shr', 'V_shr'.
 * The matrix 'S_shr[i]' is written over 'A_shr'. */
void sign_phase3_round
(
    hash_ctx_t *hash_ctx,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t E[matrix_bytes_size(PAR_M, PAR_N)],
    ff_t A[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    ff_t R[matrix_bytes_size(PAR_S, PAR_M)],
    hash_t  salt,
    uint32_t l
)
{
    uint32_t k, p;

    ff_t S[matrix_bytes_size(PAR_S, PAR_R)];

    /* Open S. */
    {
        ff_t MaL[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];
        ff_t MaR[matrix_bytes_size(PAR_M, PAR_R)];

        matrix_horizontal_split(MaL, MaR, E, PAR_M, PAR_N - PAR_R, PAR_R);

        /* 'S = R * MaR + A'.*/
        matrix_product(S, R, MaR, PAR_S, PAR_M, PAR_R);
        matrix_add(S, A, PAR_S, PAR_R);
    }
    /* * */

    for (k = 0; k < D_DIMENSION; k++)
    {
        /* Hash 'H_k_l'. */
        hash_t H_k_l;

        /* Hashing context to compute 'H_k_l'. */
        hash_ctx_t hash_H_k_l_ctx;


        /* Initialize the hashing context to compute 'H_k_l'. */
        hash_init(&hash_H_k_l_ctx);
        hash_update(&hash_H_k_l_ctx, salt, HASH_SIZE);
        hash_update(&hash_H_k_l_ctx, (uint8_t *)&l, sizeof(l));

        ff_t V_main_shr_k_j[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];

            ff_t MaL_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];
            ff_t MaR_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_R)];
            ff_t E_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_N)];
            ff_t RxMaL_main_shr_k_j[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
            ff_t RxMaR_main_shr_k_j[matrix_bytes_size(PAR_S, PAR_R)];

#if defined (OFFLINE_CC)
            // if (j < N_PARTIES - 1)
                begin_offline = get_cycles();
#endif
            /* Compute the shares of all the main party (k, 0) but the last one */

            /* Set E_main_shr_k_j = (M[0] if i = 0 else 0) + sum_{p = 0}^{PAR_K-1} a_shr[i] * M[p + 1]. */
            matrix_copy(E_main_shr_k_j, M[0], PAR_M, PAR_N);
            for (p = 0; p < PAR_K; p++)
            {
                ff_t scalar;

                scalar = matrix_get_entry(a_main_shr[k][0], PAR_K, p, 0);
                matrix_add_multiple(E_main_shr_k_j, scalar, M[p + 1], PAR_M, PAR_N);
            }
            /* * */

            /* Compute 'MaL_main_shr_k_j' and 'MaR_main_shr_k_j' using the fact
             * that E_main_shr_k_j = [MaL_main_shr_k_j | MaR_main_shr_k_j]. */
            matrix_horizontal_split(MaL_main_shr_k_j, MaR_main_shr_k_j, E_main_shr_k_j, PAR_M, PAR_N - PAR_R,
                                    PAR_R);

#if defined (OFFLINE_CC)
            // if (j < N_PARTIES - 1)
                offline_cc = offline_cc + (get_cycles() - begin_offline);
#endif

            /* Compute 'R * MaL_main_shr_k_j'. */
            matrix_product(RxMaL_main_shr_k_j, R, MaL_main_shr_k_j, PAR_S, PAR_M, PAR_N - PAR_R);

            /* Set V_main_shr_k_j = S * K_main_shr[k][0] - R * MaL_main_shr_k_j - C_main_shr_k_j. */
            matrix_product(V_main_shr_k_j, S, K_main_shr[k][0], PAR_S, PAR_R, PAR_N - PAR_R);
            matrix_subtract(V_main_shr_k_j, RxMaL_main_shr_k_j, PAR_S, PAR_N - PAR_R);
            matrix_subtract(V_main_shr_k_j, C_main_shr[k][0], PAR_S, PAR_N - PAR_R);

            /* Overwrite 'A_main_shr[k][0]' with 'S_main_shr[k][0] = R * MaR_main_shr_k_j + A_main_shr[k][0]'. */
            matrix_product(RxMaR_main_shr_k_j, R, MaR_main_shr_k_j, PAR_S, PAR_M, PAR_R);
            matrix_add(A_main_shr[k][0], RxMaR_main_shr_k_j, PAR_S, PAR_R);
            /* Hash 'S_main_shr_k_j', 'V_main_shr_k_j' (NOTE: 'A_main_shr[k][0]' has been overwritten with 'S_main_shr[k][0]'). */
            hash_update(&hash_H_k_l_ctx, A_main_shr[k][0], matrix_bytes_size(PAR_S, PAR_R));
            hash_update(&hash_H_k_l_ctx, V_main_shr_k_j, matrix_bytes_size(PAR_S, PAR_N - PAR_R));

            /* Compute the shares the last main party (k,1) */
            matrix_copy(A_main_shr[k][1], S, PAR_S, PAR_R);
            matrix_subtract(A_main_shr[k][1], A_main_shr[k][0], PAR_S, PAR_R);
            hash_update(&hash_H_k_l_ctx, A_main_shr[k][1], matrix_bytes_size(PAR_S, PAR_R));
            hash_update(&hash_H_k_l_ctx, V_main_shr_k_j, matrix_bytes_size(PAR_S, PAR_N - PAR_R));


        /* Compute 'H_k_l'. */
        hash_finalize(&hash_H_k_l_ctx, H_k_l);

        /* Hash 'H_k_l'. */
        hash_update(hash_ctx, H_k_l, HASH_SIZE);
    }
}

/* A round of Phase 1 of verifying.
 *
 * Given as input:
 * - the salt 'salt';
 * - the number 'l' of the round;
 * - the seeds 'seed' of the round;
 * - the party of the provided commitment 'i_star';
 * - the provided commitment 'com_star';
 * - the auxiliary matrices 'C_aux', 'a_aux', 'K_aux';
 * compute the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr'
 * and the commitment 'com_l'.
 */
void open_phase1_round(
    hash_t salt,
    uint32_t l,
    seed_t seed[N_PARTIES_ROUND],
    uint32_t i_star,
    hash_t com_star,
    ff_t C_aux[matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_aux[matrix_bytes_size(PAR_K, 1)],
    ff_t K_aux[matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    hash_t com_l,

    ff_t A_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)])
{
    uint32_t i, j, k, i_star_k;

    /* Hashing context for 'com_l'. */
    hash_ctx_t hash_ctx;

    /* Initialize main parties shares to zero.*/
    for (i = 0; i < D_DIMENSION; i++)
    {
        for (j = 0; j < N_PARTIES; j++)
        {
            matrix_init_zero(A_main_shr[i][j], PAR_S, PAR_R);
            matrix_init_zero(C_main_shr[i][j], PAR_S, PAR_N - PAR_R);
            matrix_init_zero(a_main_shr[i][j], PAR_K, 1);
            matrix_init_zero(K_main_shr[i][j], PAR_R, PAR_N - PAR_R);
        }
    }
    /* * */

    /* Initialize the hashing context to compute 'com_l'. */
    hash_init(&hash_ctx);
    hash_update(&hash_ctx, salt, HASH_SIZE);
    hash_update(&hash_ctx, (uint8_t *)&l, sizeof(l));

    for (i = 0; i < N_PARTIES_ROUND; i++)
    {
        if (i == i_star)
        {
            /* Copy the commitment from the signature. */
            hash_update(&hash_ctx, com_star, HASH_SIZE);
        }
        else
        {
            prng_t prng;
            hash_t com_i_l;

            ff_t A_shr_i[matrix_bytes_size(PAR_S, PAR_R)];
            ff_t C_shr_i[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
            ff_t a_shr_i[matrix_bytes_size(PAR_K, 1)];
            ff_t K_shr_i[matrix_bytes_size(PAR_R, PAR_N - PAR_R)];

            /* Initialize PRNG from 'seed[i]'. */
            prng_init(&prng, salt, seed[i]);

            /* Generate the random matrix 'A_shr_i'. */
            matrix_init_random(A_shr_i, PAR_S, PAR_R, &prng);

            if (i != N_PARTIES_ROUND - 1)
            {
                /* Generate random matrices 'a_shr_i', 'C_shr_i', 'K_shr_i'. */
                matrix_init_random(a_shr_i, PAR_K, 1, &prng);
                matrix_init_random(C_shr_i, PAR_S, PAR_N - PAR_R, &prng);
                matrix_init_random(K_shr_i, PAR_R, PAR_N - PAR_R, &prng);

                /* Compute the commitment. */
                hash_digest0(com_i_l, salt, l, i, seed[i]);

                hash_update(&hash_ctx, com_i_l, HASH_SIZE);
            }
            else
            {
                /* Compute the commitment. */
                hash_digest0_aux(com_i_l, salt, l, i, seed[N_PARTIES_ROUND - 1],
                                 a_aux, K_aux, C_aux);

                hash_update(&hash_ctx, com_i_l, HASH_SIZE);

                /* Copy the matrices from the signature. */
                matrix_copy(a_shr_i, a_aux, PAR_K, 1);
                matrix_copy(K_shr_i, K_aux, PAR_R, PAR_N - PAR_R);
                matrix_copy(C_shr_i, C_aux, PAR_S, PAR_N - PAR_R);
            }

            /* Aggregate to the main parties */
            for (k = 0; k < D_DIMENSION; k++)
            {
                i_star_k = (i >> k) & 1;
                matrix_add(A_main_shr[k][i_star_k], A_shr_i, PAR_S, PAR_R);
                matrix_add(C_main_shr[k][i_star_k], C_shr_i, PAR_S, PAR_N - PAR_R);
                matrix_add(a_main_shr[k][i_star_k], a_shr_i, PAR_K, 1);
                matrix_add(K_main_shr[k][i_star_k], K_shr_i, PAR_R, PAR_N - PAR_R);
            }
        }
    }

    /* Compute 'com_l'. */
    hash_finalize(&hash_ctx, com_l);
}

/* A round of Phase 3 of verifying.
 *
 * Given as input:
 * - the hashing context 'hash_ctx';
 * - the party id 'i_star';
 * - the public key 'M';
 * - the additive share 'S_star' of party 'i_star';
 * - the additive shares 'A_main_shr', 'C_main_shr', 'a_main_shr', 'K_main_shr';
 * - the salt 'salt';
 * - the round number 'l';
 * compute the additive shares 'V_shr' and hash them. */
void open_phase3_round
(
    hash_ctx_t *hash_ctx,
    uint32_t i_star,
    ff_t M[PAR_K + 1][matrix_bytes_size(PAR_M, PAR_N)],
    ff_t R[matrix_bytes_size(PAR_S, PAR_M)],
    ff_t S_star[matrix_bytes_size(PAR_S, PAR_R)],
    ff_t A_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)],
    ff_t C_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)],
    ff_t a_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)],
    ff_t K_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)],
    hash_t salt,
    uint32_t l
)
{
    uint32_t k, j, p, i_star_k, not_i_start_k;

    ff_t S[matrix_bytes_size(PAR_S, PAR_R)];

    /* Iteration `k=0`. In this iteration we compute the matrix `S`*/
    for (k = 0; k < D_DIMENSION; k++)
    {
        ff_t S_main_shr_k[N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];
        ff_t V_main_shr_k[N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
        ff_t MaL_main_shr_k[N_PARTIES][matrix_bytes_size(PAR_M, PAR_N - PAR_R)];

        hash_t H_k_l;

        /* Hash context to compute 'H_k_l'.*/
        hash_ctx_t hash_H_k_l_ctx;

        /* Initialize the hashing context to compute 'H_k_l'. */
        hash_init(&hash_H_k_l_ctx);
        hash_update(&hash_H_k_l_ctx, salt, HASH_SIZE);
        hash_update(&hash_H_k_l_ctx, (uint8_t * ) & l, sizeof(l));

        i_star_k = (i_star >> k) & 1;
        not_i_start_k = (i_star_k + 1) % 2;

        if (k == 0)
        {

            for (j = 0; j < N_PARTIES; j++)
            {

                ff_t E_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_N)];
                ff_t MaR_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_R)];

                /* Set E_main_shr_k_j = (M[0] if j = 0 else 0) + sum_{j = 0}^{PAR_K-1} a_main_shr[k][j] * M[j + 1]. */
                if (j == 0)
                {
                    matrix_copy(E_main_shr_k_j, M[0], PAR_M, PAR_N);
                }
                else
                {
                    matrix_init_zero(E_main_shr_k_j, PAR_M, PAR_N);
                }

                for (p = 0; p < PAR_K; p++) {
                    ff_t scalar;

                    scalar = matrix_get_entry(a_main_shr[k][j], PAR_K, p, 0);
                    matrix_add_multiple(E_main_shr_k_j, scalar, M[p + 1], PAR_M, PAR_N);
                }
                /* * */

                /* Compute 'MaL_main_shr[k][j]' and 'MaR_main_shr_k_j' using the fact
                 * that E_main_shr_k_j = [MaL_main_shr[k][j] | MaR_main_shr_k_j]. */
                matrix_horizontal_split(MaL_main_shr_k[j], MaR_main_shr_k_j, E_main_shr_k_j, PAR_M, PAR_N - PAR_R,
                                        PAR_R);


                /* Set 'S_main_shr[k][j] = R * MaR_main_shr_k_j + A_main_shr[k][j]'. */
                matrix_product(S_main_shr_k[j], R, MaR_main_shr_k_j, PAR_S, PAR_M, PAR_R);
                matrix_add(S_main_shr_k[j], A_main_shr[k][j], PAR_S, PAR_R);

                if (j == i_star_k) {
                    matrix_add(S_main_shr_k[j], S_star, PAR_S, PAR_R);
                }
            }

            /* Open S. */
            matrix_init_zero(S, PAR_S, PAR_R);

            for (j = 0; j < N_PARTIES; j++)
            {
                matrix_add(S, S_main_shr_k[j], PAR_S, PAR_R);
            }
            /* * */
            /* End of iteration `k=0`*/
        }
        /* Iterations `k > 0` */
        else
        {
            j = not_i_start_k;
            ff_t E_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_N)];
            ff_t MaR_main_shr_k_j[matrix_bytes_size(PAR_M, PAR_R)];

            /* Set E_main_shr_k_j = (M[0] if j = 0 else 0) + sum_{j = 0}^{PAR_K-1} a_main_shr[k][j] * M[j + 1]. */
            if (j == 0)
            {
                matrix_copy(E_main_shr_k_j, M[0], PAR_M, PAR_N);
            }
            else
            {
                matrix_init_zero(E_main_shr_k_j, PAR_M, PAR_N);
            }

            for (p = 0; p < PAR_K; p++)
            {
                ff_t scalar;

                scalar = matrix_get_entry(a_main_shr[k][j], PAR_K, p, 0);
                matrix_add_multiple(E_main_shr_k_j, scalar, M[p + 1], PAR_M, PAR_N);
            }
            /* * */

            /* Compute 'MaL_main_shr[k][j]' and 'MaR_main_shr_k_j' using the fact
             * that E_main_shr_k_j = [MaL_main_shr[k][j] | MaR_main_shr_k_j]. */
            matrix_horizontal_split(MaL_main_shr_k[j], MaR_main_shr_k_j, E_main_shr_k_j, PAR_M, PAR_N - PAR_R,
                                    PAR_R);


            /* Set 'S_main_shr[k][j] = R * MaR_main_shr_k_j + A_main_shr[k][j]'. */
            matrix_product(S_main_shr_k[j], R, MaR_main_shr_k_j, PAR_S, PAR_M, PAR_R);
            matrix_add(S_main_shr_k[j], A_main_shr[k][j], PAR_S, PAR_R);

            if (j == i_star_k)
            {
                matrix_add(S_main_shr_k[j], S_star, PAR_S, PAR_R);
            }
            matrix_copy(S_main_shr_k[i_star_k], S, PAR_S, PAR_R);
            matrix_add(S_main_shr_k[i_star_k], S_main_shr_k[j], PAR_S, PAR_R);
        }

        ff_t RxMaL_main_shr_k_j[matrix_bytes_size(PAR_S, PAR_N - PAR_R)];

        /* Compute 'R * MaL_main_shr[k][~ i_start_k]'. */
        matrix_product(RxMaL_main_shr_k_j, R, MaL_main_shr_k[not_i_start_k], PAR_S, PAR_M, PAR_N - PAR_R);

        /* Set V_main_shr[k][j] = S * K_main_shr[k][j] - R * MaL_main_shr_k_j - C_main_shr[k][j]. */
        matrix_product(V_main_shr_k[not_i_start_k], S, K_main_shr[k][not_i_start_k], PAR_S, PAR_R, PAR_N - PAR_R);
        matrix_subtract(V_main_shr_k[not_i_start_k], RxMaL_main_shr_k_j, PAR_S, PAR_N - PAR_R);
        matrix_subtract(V_main_shr_k[not_i_start_k], C_main_shr[k][not_i_start_k], PAR_S, PAR_N - PAR_R);

        /* Hash 'S_main_shr[k][j]', 'V_main_shr[k][j]'. */

        hash_update(&hash_H_k_l_ctx, S_main_shr_k[0], matrix_bytes_size(PAR_S, PAR_R));
        hash_update(&hash_H_k_l_ctx, V_main_shr_k[not_i_start_k], matrix_bytes_size(PAR_S, PAR_N - PAR_R));
        hash_update(&hash_H_k_l_ctx, S_main_shr_k[1], matrix_bytes_size(PAR_S, PAR_R));
        hash_update(&hash_H_k_l_ctx, V_main_shr_k[not_i_start_k], matrix_bytes_size(PAR_S, PAR_N - PAR_R));

        hash_finalize(&hash_H_k_l_ctx, H_k_l);
        hash_update(hash_ctx, H_k_l, HASH_SIZE);
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
    hash_t com[TAU];
    hash_t com_star[TAU];

    /* Hashing context. */
    hash_ctx_t hash1_ctx;
    hash_ctx_t hash2_ctx;

    /* Seed trees. */
    seed_t trees[TAU][TREE_N_NODES];

    /* Matrices. */

    /* Main parties shares of 'A', 'C', 'a', 'K'. */
    ff_t A_rnd_main_shr[TAU][D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];
    ff_t C_rnd_main_shr[TAU][D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t a_rnd_main_shr[TAU][D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)];
    ff_t K_rnd_main_shr[TAU][D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];

    /* Matrices 'A'. */
    ff_t A_rnd[TAU][matrix_bytes_size(PAR_S, PAR_R)];

    /* Shares of the last party. */
    ff_t C_aux[TAU][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
    ff_t a_aux[TAU][matrix_bytes_size(PAR_K, 1)];
    ff_t K_aux[TAU][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];

    /* Matrices of the 'i_star' party to recompute 'S_start'. */
    ff_t S_rnd_star[TAU][matrix_bytes_size(PAR_S, PAR_R)];
    ff_t A_star[matrix_bytes_size(PAR_S, PAR_R)];
    ff_t E_star[matrix_bytes_size(PAR_M, PAR_N)];
    ff_t a_star[matrix_bytes_size(PAR_K, 1)];
    ff_t MaR_star[matrix_bytes_size(PAR_M, PAR_R)];
    ff_t MaL_star[matrix_bytes_size(PAR_M, PAR_N - PAR_R)];

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
    hash_init(&hash1_ctx);
    hash_update(&hash1_ctx, salt, HASH_SIZE);
    hash_update(&hash1_ctx, msg, msg_len);

    for (l = 0; l < TAU; l++)
    {
        seed_t seed;

#if defined (OFFLINE_CC)
begin_offline = get_cycles();
#endif

        /* Generate a random seed. */
        randombytes(seed, SEED_SIZE);

        /* Generate a seed tree from the random seed. */
        seed_tree_init(trees[l], salt, seed);

#if defined (OFFLINE_CC)
offline_cc = offline_cc + (get_cycles() - begin_offline);
#endif

        /* Compute the matrices
         * 'A_rnd[l][0...N_PARTIES_ROUND-1]',
         * 'A_rnd_main_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'C_rnd_main_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'a_rnd_main_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'K_rnd_main_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'C_aux[l]',
         * 'a_aux[l]',
         * 'K_aux[l]'
         * and the commitments 'com[l]'. */
        sign_phase1_round(A_rnd[l], A_rnd_main_shr[l], C_rnd_main_shr[l],
                          a_rnd_main_shr[l], K_rnd_main_shr[l],
                          C_aux[l], a_aux[l], K_aux[l],
                          com[l],
                          salt, l, (uint8_t (*)[SEED_SIZE]) seed_tree_get_leaves(trees[l]),
                          a, K);

        /* Hash the commitment for this round. */
        hash_update(&hash1_ctx, com[l], HASH_SIZE);
    }

    /* Compute 'hash1'. */
    hash_finalize(&hash1_ctx, hash1);

    /* Initialize the PRNG to generate the first challenges. */
    prng_init(&prng, hash1, NULL);

    /* Initialize the hashing context to compute 'hash2'. */
    hash_init(&hash2_ctx);
    hash_update(&hash2_ctx, salt, HASH_SIZE);
    hash_update(&hash2_ctx, msg, msg_len);

    ff_t R[TAU][matrix_bytes_size(PAR_S, PAR_M)];

    for (l = 0; l < TAU; l++)
    {
        /* Generate the first challenge 'R' for this round. */
        matrix_init_random(R[l], PAR_S, PAR_M, &prng);

        /* Compute and hash the matrices
         * 'S_main_rnd_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]' and 'V_main_rnd_shr[l][0...D_DIMENSION-1][0...N_PARTIES-1]'. */
        sign_phase3_round(&hash2_ctx, M, E, A_rnd[l],
                          A_rnd_main_shr[l], C_rnd_main_shr[l],
                          a_rnd_main_shr[l], K_rnd_main_shr[l],
                          R[l], salt, l);
    }

    /* Compute 'hash2'. */
    hash_update(&hash2_ctx, hash1, HASH_SIZE);
    hash_finalize(&hash2_ctx, hash2);

    /* Generate the second challenges. */
    get_second_challenges(i_star, hash2);

    /* Recompute 'com_star[l]' and 'S_rnd_star[l]'. */
    for (l = 0; l < TAU; l++)
    {
        seed_t *seed_round;

        seed_round = (uint8_t (*)[SEED_SIZE]) seed_tree_get_leaves(trees[l]);

        /* Initialize PRNG from 'seed[i]'. */
        prng_init(&prng, salt, seed_round[i_star[l]]);

        /* Generate the random matrix 'A_shr[i]'. */
        matrix_init_random(A_star, PAR_S, PAR_R, &prng);

        if (i_star[l] != N_PARTIES_ROUND - 1)
        {
            hash_digest0(com_star[l], salt, l, i_star[l], seed_round[i_star[l]]);

            /* Generate random matrices 'a_shr[seed_i_star]', 'C_shr[seed_i_star]', 'K_shr[seed_i_star]'. */
            matrix_init_random(a_star, PAR_K, 1, &prng);
        }
        else
        {
            hash_digest0_aux(com_star[l], salt, l, i_star[l], seed_round[i_star[l]],
                             a_aux[l], K_aux[l], C_aux[l]);

            /* Copy 'a_aux[l]' into 'a_star'.*/
            matrix_copy(a_star, a_aux[l], PAR_K, 1);
        }

        /* Compute 'E_star'. */
        computeE_shr(E_star, a_star, M);

        /* Compute 'MaL_star' and 'MaR_star' using the fact
        * that E_star = [MaL_star| -MaR_star. */
        matrix_horizontal_split(MaL_star, MaR_star, E_star, PAR_M, PAR_N - PAR_R, PAR_R);
        matrix_negate(MaR_star, PAR_M, PAR_R);

        /* Set 'S_main_shr[k][j] = R * MaR_main_shr_k_j + A_main_shr[k][j]'. */
        matrix_product(S_rnd_star[l], R[l], MaR_star, PAR_S, PAR_M, PAR_R);
        matrix_add(S_rnd_star[l], A_star, PAR_S, PAR_R);
    }

    /* Pack the signature. */
    pack_signature(sig_msg, &sig_len, salt, hash1, hash2, i_star, trees, com_star,
                   a_aux, K_aux, C_aux, S_rnd_star);

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
    hash_ctx_t hash1_c_ctx;
    hash_ctx_t hash2_c_ctx;

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
    hash_init(&hash1_c_ctx);
    hash_update(&hash1_c_ctx, salt, HASH_SIZE);
    hash_update(&hash1_c_ctx, sig_msg + sig_len, *msg_len);

    /* Initialize the hashing context to compute 'hash2_c'. */
    hash_init(&hash2_c_ctx);
    hash_update(&hash2_c_ctx, salt, HASH_SIZE);
    hash_update(&hash2_c_ctx, sig_msg + sig_len, *msg_len);

    for (l = 0; l < TAU; l++)
    {
        /* Matrices of this round. */
        ff_t A_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_R)];
        ff_t C_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_S, PAR_N - PAR_R)];
        ff_t a_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_K, 1)];
        ff_t K_main_shr[D_DIMENSION][N_PARTIES][matrix_bytes_size(PAR_R, PAR_N - PAR_R)];
        ff_t R[matrix_bytes_size(PAR_S, PAR_M)];

        /* Commitment of Phase 1 of this round. */
        hash_t com_l;

        /* Seeds of this round. */
        seed_t seed[N_PARTIES_ROUND];

        /* Decompress the seed tree. */
        seed_tree_unpack(seed, salt, packed_tree[l], i_star[l]);

        /* Compute the shares of the main parties, without aggregation of the share i_star
         * 'A_main_shr[0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'C_main_shr[0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'a_main_shr[0...D_DIMENSION-1][0...N_PARTIES-1]',
         * 'K_main_shr[0...D_DIMENSION-1][0...N_PARTIES-1]',
         * while hashing the commitments. */
        open_phase1_round(salt, l, seed, i_star[l], com_star[l],
                            C_rnd_aux[l], a_rnd_aux[l], K_rnd_aux[l],
                          com_l, A_main_shr, C_main_shr, a_main_shr, K_main_shr);

        /* Hash the commitment for this round. */
        hash_update(&hash1_c_ctx, com_l, HASH_SIZE);

        /* Compute the first challenge of this round. */
        matrix_init_random(R, PAR_S, PAR_M, &prng);

        /* Compute and hash the matrices 'S1[0...N_PARTIES-1]', 'S2[0...N_PARTIES-1]', 'V[0...N_PARTIES-1]'. */
        open_phase3_round(&hash2_c_ctx, i_star[l],
                          M, R, S_rnd_star[l],
                          A_main_shr, C_main_shr, a_main_shr, K_main_shr,
                          salt, l);
    }

    /* Compute 'hash1_c'. */
    hash_finalize(&hash1_c_ctx, hash1_c);

    /* Compute 'hash2_c'. */
    hash_update(&hash2_c_ctx, hash1_c, HASH_SIZE);
    hash_finalize(&hash2_c_ctx, hash2_c);

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

