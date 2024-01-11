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
  uint8_t* data;
  uint8_t* exists;
  uint8_t* have_value;
  size_t   seed_size;
  size_t   num_leaves;
  size_t   num_nodes;
} tree_t;

typedef struct reveal_list_t
{
  uint8_t* data;
  size_t   seed_size;
  size_t   missing_leaf;
} reveal_list_t;

uint32_t ceil_log2(uint32_t x);

tree_t* make_seed_tree(const uint8_t* seed, const size_t seed_size,
                       const uint8_t* salt, const size_t salt_size,
                       const size_t num_leaves, const size_t repetition_index);

tree_t* reconstruct_seed_tree(const reveal_list_t* reveal_list,
                              const uint8_t* salt, const size_t salt_size,
                              const size_t num_leaves,
                              const size_t repetition_index);

reveal_list_t reveal_all_but(const tree_t* tree, size_t leaf_index);

uint8_t* get_leaf(tree_t* tree, size_t leaf_index);

void free_tree(tree_t* tree);

#endif // TREE_H
