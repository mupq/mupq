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

#ifndef SEED_TREE_H
#define SEED_TREE_H

#include "hash_types.h"
#include "params.h"

/* Initialize 'tree' from 'salt' and 'seed'. */
void seed_tree_init(seed_t tree[TREE_N_NODES], const hash_t salt, const seed_t seed);

/* Return a pointer to the 'N_PARTIES_ROUND' leaves of 'tree'. */
uint8_t *seed_tree_get_leaves(seed_t tree[TREE_N_NODES]);

/* Pack 'tree' into 'packed_tree' to recover later all seeds but the 'i0'th. */
void seed_tree_pack(seed_t packed_tree[TREE_HEIGHT], seed_t tree[TREE_N_NODES], int i0);

/* Unpack 'packed_tree' to 'seeds', recovering all seeds but the 'i0'th. */
void seed_tree_unpack(seed_t seeds[N_PARTIES_ROUND], seed_t packed_tree[TREE_HEIGHT], int i0);

#endif
