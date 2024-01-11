// -----------------------------------------------------------------------------
// File Name   : tree.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "tree.h"
#include <stdlib.h>

static uint32_t number_of_leading_zeros(uint32_t x);
int node_exists(const tree_t* tree, size_t index);
int node_has_value(tree_t* tree, size_t index);
static size_t get_parent(size_t node);
int has_sibling(const tree_t* tree, size_t index);
int get_sibling(const tree_t* tree, size_t index);
void expand_seed();
void expand_seed_x4();

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b)) // return the lesser of two values
#endif

// H.S. Warren, "Hacker's Delight", Pearson Education, 2013.
// http://www.hackersdelight.org/hdcodetxt/nlz.c.txt (not accessible, 2022.6.25)
// https://gitlab.com/hcs0/Hackers-Delight/-/blob/master/nlz.c.txt
// (accessible, 2022.6.25)
static uint32_t number_of_leading_zeros(uint32_t x)
{
  uint32_t n;

  if (x == 0)
  {
    return (32);
  }

  n = 1;
  if ((x >> 16) == 0) { n = n + 16; x = x << 16; }
  if ((x >> 24) == 0) { n = n +  8; x = x <<  8; }
  if ((x >> 28) == 0) { n = n +  4; x = x <<  4; }
  if ((x >> 30) == 0) { n = n +  2; x = x <<  2; }
  n = n - (x >> 31);

  return n;
}

uint32_t ceil_log2(uint32_t x)
{
  if (x == 0)
  {
    return 0;
  }

  return 32 - number_of_leading_zeros(x - 1);
}

int node_exists(const tree_t* tree, size_t index)
{
  if (index >= tree->num_nodes)
  {
    return 0;
  }

  if (tree->exists[index])
  {
    return 1;
  }

  return 0;
}

int node_has_value(tree_t* tree, size_t index)
{
  if (index >= tree->num_nodes)
  {
    return 0;
  }
  return tree->have_value[index];
}

static size_t get_parent(size_t node)
{
  assert(node != 0);
  return ((node + 1) >> 1) - 1;
}

int has_sibling(const tree_t* tree, size_t index)
{
  if (!node_exists(tree, index))
  {
    return 0;
  }
  if ((index % 2 == 1) && !node_exists(tree, index + 1))
  {
    return 0;
  }
  return 1;
}

int get_sibling(const tree_t* tree, size_t index)
{
  assert(index < tree->num_nodes);
  assert(index != 0);
  assert(has_sibling(tree, index));

  if ((index % 2 == 1))
  {
    if (index + 1 < tree->num_nodes)
    {
      return index + 1;
    }
    else
    {
      assert(!"get_sibling: request for node with no sibling");
      return 0;
    }
  }
  else
  {
    return index - 1;
  }
}

void expand_seed(const uint8_t* seed, const size_t seed_size,
                 const uint8_t* salt, const size_t salt_size,
                 const size_t repetition_index, const size_t node_index,
                 uint8_t* out)
{
  hash_instance ctx;

  hash_init_prefix(&ctx, seed_size << 1, HASH_PREFIX_3);
  hash_update(&ctx, seed, seed_size);
  hash_update(&ctx, salt, salt_size);
  hash_update_uint16(&ctx, (uint16_t)repetition_index);
  hash_update_uint16(&ctx, (uint16_t)node_index);
  hash_final(&ctx);
  hash_squeeze(&ctx, out, 2 * seed_size);
}

void expand_seed_x4(const uint8_t* seed0,
                    const uint8_t* seed1,
                    const uint8_t* seed2,
                    const uint8_t* seed3, const size_t seed_size,
                    const uint8_t* salt, const size_t salt_size,
                    const size_t repetition_index, const size_t node_index,
                    uint8_t* out0, uint8_t* out1, uint8_t* out2, uint8_t* out3)
{
  hash_instance_x4 ctx;

  hash_init_prefix_x4(&ctx, seed_size << 1, HASH_PREFIX_3);
  hash_update_x4_4(&ctx, seed0, seed1, seed2, seed3, seed_size);
  hash_update_x4_1(&ctx, salt, salt_size);
  hash_update_x4_uint16(&ctx, (uint16_t)repetition_index);

  uint16_t node_indices[4] = {(uint16_t)(node_index), 
                              (uint16_t)(node_index + 1),
                              (uint16_t)(node_index + 2),
                              (uint16_t)(node_index + 3)};
  hash_update_x4_uint16s(&ctx, node_indices);
  hash_final_x4(&ctx);
  hash_squeeze_x4_4(&ctx, out0, out1, out2, out3, 2 * seed_size);
}

tree_t* make_seed_tree(const uint8_t* seed, const size_t seed_size,
                       const uint8_t* salt, const size_t salt_size, 
                       const size_t num_leaves, const size_t repetition_index)
{
  size_t   i;
  tree_t*  tree = malloc(sizeof(tree_t));

  size_t   tree_depth = 1 + ceil_log2(num_leaves);
  size_t   last_non_leaf, num_nodes;
  uint8_t* buffer;

  tree->num_nodes  = ((1 << (tree_depth)) - 1) -
                     ((1 << (tree_depth - 1)) - num_leaves);
  tree->num_leaves = num_leaves;
  tree->seed_size  = seed_size;

  num_nodes = tree->num_nodes;

  tree->data       = malloc(num_nodes * seed_size);
  tree->have_value = calloc(num_nodes, 1);
  tree->exists     = calloc(num_nodes, 1);

  for (i = num_nodes - num_leaves; i < num_nodes; i++)
  {
    tree->exists[i] = 1;
  }

  for (i = num_nodes - num_leaves; i > 0; i--)
  {
    if (node_exists(tree, 2 * i + 1) || node_exists(tree, 2 * i + 2))
    {
      tree->exists[i] = 1;
    }
  }

  tree->exists[0]     = 1;
  tree->have_value[0] = 1;
  memcpy(tree->data, seed, seed_size);

  buffer        = malloc(2 * 4 * seed_size);
  last_non_leaf = get_parent(num_nodes - 1);

  for (i = 0; i <= MIN(last_non_leaf, 2); i++)
  {
    if (!node_exists(tree, i))
    {
      continue;
    }

    expand_seed(&tree->data[i * seed_size], seed_size, salt, salt_size,
                repetition_index, i, buffer);
    tree->have_value[2 * i + 1] = 1;

    if (node_exists(tree, 2 * i + 2))
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, 2 * seed_size);
      tree->have_value[2 * i + 2] = 1;
    }
    else
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, seed_size);
    }
  }

  uint8_t* seeds = malloc(4 * seed_size);
  for (; i < (last_non_leaf / 4) * 4; i += 4)
  {
    memset(seeds, 0x00, 4 * seed_size);
    for (size_t j = 0; j < 4; j++)
    {
      if (node_exists(tree, i + j))
      {
        memcpy(seeds + j * seed_size, &tree->data[(i + j) * seed_size], seed_size);
        tree->have_value[2 * (i + j) + 1] = 1;
      }
    }
    expand_seed_x4(seeds, seeds + seed_size, seeds + 2 * seed_size,
                   seeds + 3 * seed_size, seed_size,
                   salt, salt_size, repetition_index, i,
                   buffer, buffer + 2 * seed_size,
                   buffer + 4 * seed_size, buffer + 6 * seed_size);

    for (size_t j = 0; j < 4; j++)
    {
      size_t index = i + j;
      if (node_exists(tree, index))
      {
        if (node_exists(tree, 2 * index + 2))
        {
          memcpy(&tree->data[(2 * index + 1) * seed_size], 
                buffer + 2 * j * seed_size, 2 * seed_size);
          tree->have_value[2 * index + 2] = 1;
        }
        else
        {
          memcpy(&tree->data[(2 * index + 1) * seed_size], 
                buffer + 2 * j * seed_size, seed_size);
        }
      }
    } 
  }
  free(seeds);

  for (; i <= last_non_leaf; i++)
  {
    if (!node_exists(tree, i))
    {
      continue;
    }

    expand_seed(&tree->data[i * seed_size], seed_size, salt, salt_size,
                repetition_index, i, buffer);
    tree->have_value[2 * i + 1] = 1;

    if (node_exists(tree, 2 * i + 2))
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, 2 * seed_size);
      tree->have_value[2 * i + 2] = 1;
    }
    else
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, seed_size);
    }
  }
  free(buffer);
  return tree;
}

tree_t* reconstruct_seed_tree(const reveal_list_t* reveal_list,
                              const uint8_t* salt, const size_t salt_size,
                              const size_t num_leaves,
                              const size_t repetition_index)
{
  size_t   i;
  tree_t*  tree = malloc(sizeof(tree_t));

  size_t   tree_depth = 1 + ceil_log2(num_leaves);
  size_t   last_non_leaf, seed_size, num_nodes;
  uint8_t* buffer;

  tree->num_nodes  = ((1 << (tree_depth)) - 1) -
                     ((1 << (tree_depth - 1)) - num_leaves);
  tree->num_leaves = num_leaves;
  tree->seed_size  = seed_size = reveal_list->seed_size;

  num_nodes = tree->num_nodes;

  tree->data       = malloc(num_nodes * seed_size);
  tree->have_value = calloc(num_nodes, 1);
  tree->exists     = calloc(num_nodes, 1);

  for (i = num_nodes - num_leaves; i < num_nodes; i++)
  {
    tree->exists[i] = 1;
  }

  for (i = num_nodes - num_leaves; i > 0; i--)
  {
    if (node_exists(tree, 2 * i + 1) || node_exists(tree, 2 * i + 2))
    {
      tree->exists[i] = 1;
    }
  }

  tree->exists[0] = 1;

  size_t path_index = 0;
  for (size_t node = num_nodes - num_leaves + reveal_list->missing_leaf; node != 0;
       node = get_parent(node))
  {
    if (!has_sibling(tree, node))
    {
      path_index++;
      continue;
    }

    size_t sibling = get_sibling(tree, node);
    memcpy(&tree->data[sibling * seed_size],
           &reveal_list->data[path_index * seed_size], seed_size);

    tree->have_value[sibling] = 1;
    path_index++;
  }

  buffer        = malloc(2 * 4 * seed_size);
  last_non_leaf = get_parent(num_nodes - 1);

  for (i = 0; i <= MIN(last_non_leaf, 2); i++)
  {
    if (!node_exists(tree, i) || !node_has_value(tree, i))
    {
      continue;
    }

    expand_seed(&tree->data[i * seed_size], seed_size, salt, salt_size,
                repetition_index, i, buffer);
    tree->have_value[2 * i + 1] = 1;

    if (node_exists(tree, 2 * i + 2))
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, 2 * seed_size);
      tree->have_value[2 * i + 2] = 1;
    }
    else
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, seed_size);
    }
  }

  uint8_t* seeds = malloc(4 * seed_size);
  for (; i < (last_non_leaf / 4) * 4; i += 4)
  {
    memset(seeds, 0x00, 4 * seed_size);
    for (size_t j = 0; j < 4; j++)
    {
      if (node_exists(tree, i + j) && node_has_value(tree, i + j))
      {
        memcpy(seeds + j * seed_size, &tree->data[(i + j) * seed_size], seed_size);
        tree->have_value[2 * (i + j) + 1] = 1;
      }
    }
    expand_seed_x4(seeds, seeds + seed_size, seeds + 2 * seed_size,
                   seeds + 3 * seed_size, seed_size,
                   salt, salt_size, repetition_index, i,
                   buffer, buffer + 2 * seed_size,
                   buffer + 4 * seed_size, buffer + 6 * seed_size);

    for (size_t j = 0; j < 4; j++)
    {
      size_t index = i + j;
      if (node_exists(tree, index) && node_has_value(tree, index))
      {
        if (node_exists(tree, 2 * index + 2))
        {
          memcpy(&tree->data[(2 * index + 1) * seed_size], 
                buffer + 2 * j * seed_size, 2 * seed_size);
          tree->have_value[2 * index + 2] = 1;
        }
        else
        {
          memcpy(&tree->data[(2 * index + 1) * seed_size], 
                buffer + 2 * j * seed_size, seed_size);
        }
      }     
    } 
  }
  free(seeds);

  for (; i <= last_non_leaf; i++)
  {
    if (!node_exists(tree, i) || !node_has_value(tree, i))
    {
      continue;
    }

    expand_seed(&tree->data[i * seed_size], seed_size, salt, salt_size,
                repetition_index, i, buffer);
    tree->have_value[2 * i + 1] = 1;

    if (node_exists(tree, 2 * i + 2))
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, 2 * seed_size);
      tree->have_value[2 * i + 2] = 1;
    }
    else
    {
      memcpy(&tree->data[(2 * i + 1) * seed_size], buffer, seed_size);
    }
  }
  free(buffer);
  return tree;
}

reveal_list_t reveal_all_but(const tree_t* tree, size_t leaf_index)
{
  reveal_list_t reveal_list;

  // Calculate path up to root for missing leaf
  size_t   buffer_size = ((ceil_log2(tree->num_leaves) + 1) * 2);
  uint8_t* buffer      = calloc(buffer_size, tree->seed_size);

  size_t path_index = 0;
  for (size_t node = tree->num_nodes - tree->num_leaves + leaf_index; node != 0;
       node = get_parent(node))
  {
    if (!has_sibling(tree, node))
    {
      path_index++;
      continue;
    }

    size_t sibling = get_sibling(tree, node);
    memcpy(&buffer[path_index * tree->seed_size],
           &tree->data[sibling * tree->seed_size], tree->seed_size);

    path_index++;
  }

  reveal_list.data            = buffer;
  reveal_list.seed_size       = tree->seed_size;
  reveal_list.missing_leaf    = leaf_index;

  return reveal_list;
}

uint8_t* get_leaf(tree_t* tree, size_t leaf_index)
{
  assert(leaf_index < tree->num_leaves);

  size_t real_leaf_index = tree->num_nodes - tree->num_leaves + leaf_index;

  if(!tree->have_value[real_leaf_index])
  {
    return NULL;
  }

  return &tree->data[real_leaf_index * tree->seed_size];
}

void free_tree(tree_t* tree)
{
  if (tree != NULL)
  {
    free(tree->data);
    free(tree->have_value);
    free(tree->exists);
    free(tree);
  }
}
