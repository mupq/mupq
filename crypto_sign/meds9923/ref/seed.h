#ifndef SEED_H
#define SEED_H

#include <stdint.h>

#define SEED_TREE_ADDR(h, i) ((1 << (h)) - 1 + (i))

void print_tree(uint8_t *stree);

void t_hash(uint8_t *stree, uint8_t *salt, int h, int i);

#define STREE_TO_PATH 0
#define PATH_TO_STREE 1

void stree_to_path_to_stree(uint8_t *stree, uint8_t *h, uint8_t *path, uint8_t *salt, int mode);

#define stree_to_path(stree, h, path, salt) stree_to_path_to_stree(stree, h, path, salt, STREE_TO_PATH)
#define path_to_stree(stree, h, path, salt) stree_to_path_to_stree(stree, h, path, salt, PATH_TO_STREE)

#endif

