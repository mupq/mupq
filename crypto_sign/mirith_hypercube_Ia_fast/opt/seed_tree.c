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
#include "prng.h"
#include "seed_tree.h"

/* A seed tree has the format
 *
 * tree = [seed[0] | seed[1] | ... | seed[TREE_N_NODES - 1]]
 *
 * where the left, resp. right, child of 'seed[i]' is 'seed[2*i + 1]', resp.
 * 'seed[2*i + 2]'. The last 'N_PARTIES' seeds are the leaves of the tree.
 */

#define LEFT_CHILD(i) (2 * i + 1)
#define RIGHT_CHILD(i) (2 * i + 2)

/* Generate a subtree of given root and height. */
void generate_subtree(seed_t *subtree, const hash_t salt, const seed_t root, uint32_t height)
{
    uint32_t i;
    prng_t prng;

    /* Allocate the root as the first seed of 'tree'. */
    memcpy(&subtree[0], root, SEED_SIZE); //

    /* Generate the children nodes. */
    for (i = 0; i < ((1u << height) - 1u); i++)
    {
        prng_init(&prng, salt, subtree[i]);
        prng_bytes(&prng, subtree[LEFT_CHILD(i)], SEED_SIZE);
        prng_bytes(&prng, subtree[RIGHT_CHILD(i)], SEED_SIZE);
    }
}

/* Return a pointer to the leaves of a subtree of a given height. */
seed_t *get_leaves(seed_t subtree[], uint32_t height)
{
    return &subtree[(1 << height) - 1];
}

void seed_tree_init(seed_t tree[TREE_N_NODES], const hash_t salt, const seed_t seed)
{
    prng_t prng;
    seed_t root;

    /* Generate the root from 'salt' and 'seed'. */
    prng_init(&prng, salt, seed);
    prng_bytes(&prng, root, SEED_SIZE);

    /* Generate the tree from the root. */
    generate_subtree(tree, salt, root, TREE_HEIGHT);

}

uint8_t *seed_tree_get_leaves(seed_t tree[TREE_N_NODES])
{
    return tree[TREE_N_NODES - N_PARTIES_ROUND];
}

void seed_tree_pack(seed_t packed_tree[TREE_HEIGHT], seed_t tree[TREE_N_NODES], uint32_t i0)
{
    uint32_t j;
    uint32_t next_node;

    /* Start from the root node. */
    next_node = 0;

    for (j = 0; j < TREE_HEIGHT; j++)
    {
        if ((i0 >> (TREE_HEIGHT - 1 - j)) & 1)
        {
            /* Pack the left children and continue from the right one. */
            memcpy(packed_tree[j], tree[LEFT_CHILD(next_node)], SEED_SIZE);
            next_node = RIGHT_CHILD(next_node);
        }
        else
        {
            /* Pack the right children and continue from the left one. */
            memcpy(packed_tree[j], tree[RIGHT_CHILD(next_node)], SEED_SIZE);
            next_node = LEFT_CHILD(next_node);
        }
    }
}

void seed_tree_unpack(seed_t seeds[N_PARTIES_ROUND], const hash_t salt, seed_t packed_tree[TREE_HEIGHT], uint32_t i0)
{
    seed_t temp_tree[TREE_N_NODES];
    int j;
    uint32_t gap_covered = 0;
    uint32_t number_of_seeds_to_copy = N_PARTIES_ROUND >> 1;

    for (j = TREE_HEIGHT - 1; j > -1; j--)
    {
        uint32_t j_bit;
        uint32_t aux_initial_index;
        uint32_t initial_seed_index;

        j_bit = (i0 >> j) & 1u;

        aux_initial_index = (1u - j_bit) * (1u << j);

        initial_seed_index = aux_initial_index + gap_covered;

        gap_covered += j_bit * (1u << j);

        generate_subtree(temp_tree, salt, packed_tree[(TREE_HEIGHT - 1 - j)], j);

        memcpy(seeds[initial_seed_index], get_leaves(temp_tree, j),
            number_of_seeds_to_copy * SEED_SIZE);

        number_of_seeds_to_copy /= 2;

    }    
}
