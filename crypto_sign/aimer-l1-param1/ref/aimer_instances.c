// -----------------------------------------------------------------------------
// File Name   : aimer_instances.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "aimer_instances.h"

// {block_size, num_input_sboxes}
#define AIM_128_PARAMS {16, 2}  // AIM-I   parameters
#define AIM_192_PARAMS {24, 2}  // AIM-III parameters
#define AIM_256_PARAMS {32, 3}  // AIM-IV  parameters

static aimer_instance_t instances[PARAMETER_SET_MAX_INDEX] =
{
  {        {0, 0}, 0, 0, 0, 0, 0, 0, PARAMETER_SET_INVALID},

  // AIM_params, salt size, digest size, seed size, field size, T, N, parameter set name
  {AIM_128_PARAMS, 32, 32, 16, 16, 33,   16, AIMER_L1_PARAM1},
  {AIM_128_PARAMS, 32, 32, 16, 16, 23,   57, AIMER_L1_PARAM2},
  {AIM_128_PARAMS, 32, 32, 16, 16, 17,  256, AIMER_L1_PARAM3},
  {AIM_128_PARAMS, 32, 32, 16, 16, 13, 1615, AIMER_L1_PARAM4},
  {AIM_192_PARAMS, 48, 48, 24, 24, 49,   16, AIMER_L3_PARAM1},
  {AIM_192_PARAMS, 48, 48, 24, 24, 33,   64, AIMER_L3_PARAM2},
  {AIM_192_PARAMS, 48, 48, 24, 24, 25,  256, AIMER_L3_PARAM3},
  {AIM_192_PARAMS, 48, 48, 24, 24, 19, 1621, AIMER_L3_PARAM4},
  {AIM_256_PARAMS, 64, 64, 32, 32, 65,   16, AIMER_L5_PARAM1},
  {AIM_256_PARAMS, 64, 64, 32, 32, 44,   62, AIMER_L5_PARAM2},
  {AIM_256_PARAMS, 64, 64, 32, 32, 33,  256, AIMER_L5_PARAM3},
  {AIM_256_PARAMS, 64, 64, 32, 32, 25, 1623, AIMER_L5_PARAM4}
};

const aimer_instance_t *aimer_instance_get(aimer_params_t param)
{
  if (param <= PARAMETER_SET_INVALID || param >= PARAMETER_SET_MAX_INDEX)
  {
    return NULL;
  }

  return &instances[param];
}
