#ifndef TYPES_H
#define TYPES_H

#include "param.h"
#include <stdint.h>

typedef uint32_t fpoints_t;
typedef uint8_t fpoly_t;
typedef uint8_t fsd_t;
// TODO(stevenyue): Select fpoints=2p32 for CAT_5.
typedef uint8_t fpoints_c_t[PARAM_fpoint_size]; // compressed representation of f_points
typedef uint8_t seed_t[PARAM_seed_size];
typedef uint8_t salt_t[PARAM_salt_size];
typedef uint8_t commit_t[PARAM_commit_size];
typedef uint8_t hash_t[PARAM_hash_size];
typedef fpoly_t ha_slice_t[128] __attribute__((aligned(32)));

#endif // TYPES_H
