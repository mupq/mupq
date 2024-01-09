#ifndef MQOM_SEED_TREE_H
#define MQOM_SEED_TREE_H

#include <stdint.h>

typedef struct seed_tree_t {
    uint32_t height;      /* The height of the tree */
    uint8_t** nodes;    /* The data for each node (1-based array) */
    uint32_t nb_nodes;    /* The total number of nodes in the tree */
    uint32_t nb_leaves;   /* The total number of leaves in the tree */
} seed_tree_t;

seed_tree_t* malloc_seed_tree(uint32_t height);
void free_seed_tree(seed_tree_t* tree);
void expand_seed_tree(seed_tree_t* tree, const uint8_t* root_seed, const uint8_t* salt);
uint8_t** get_leaves(seed_tree_t* tree);
void get_seed_path(uint8_t* path, const seed_tree_t* tree, uint16_t hidden_leaf);
void reconstruct_tree(seed_tree_t* tree, uint16_t hidden_leaf, const uint8_t* path, const uint8_t* salt);

#endif /* MQOM_SEED_TREE_H */
