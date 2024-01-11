// -----------------------------------------------------------------------------
// File Name   : tree.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef TREE_H
#define TREE_H

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "hash.h"
#include "aimer_instances.h"

typedef struct tree_t
{
  uint8_t data[AIMER_TREE_NUM_NODES*AIMER_SEED_SIZE];
  uint8_t exists[AIMER_TREE_NUM_NODES];
  uint8_t have_value[AIMER_TREE_NUM_NODES];
} tree_t;

typedef struct reveal_list_t
{
  uint8_t data[AIMER_SEED_SIZE*((AIMER_LOGN + 1) * 2)];
  size_t   missing_leaf;
} reveal_list_t;

uint32_t ceil_log2(uint32_t x);

void make_seed_tree(tree_t *tree, const uint8_t* seed, const size_t seed_size,
                       const uint8_t* salt, const size_t salt_size,
                       const size_t num_leaves, const size_t repetition_index);

void reconstruct_seed_tree(tree_t *tree, const reveal_list_t* reveal_list,
                              const uint8_t* salt, const size_t salt_size,
                              const size_t num_leaves,
                              const size_t repetition_index);

void reveal_all_but(reveal_list_t *reveal_list, const tree_t* tree, size_t leaf_index);

uint8_t* get_leaf(tree_t* tree, size_t leaf_index);

#endif // TREE_H
