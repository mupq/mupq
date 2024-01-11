// -----------------------------------------------------------------------------
// File Name   : aimer_instances.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIMER_INSTANCES_H
#define AIMER_INSTANCES_H

#include "aim.h"
#include "aimer.h"

typedef struct aimer_instance_t
{
  aim_params_t aim_params;

  uint32_t     salt_size;
  uint32_t     digest_size;      // bytes
  uint32_t     seed_size;        // bytes
  uint32_t     field_size;

  uint32_t     num_repetitions;  // tau
  uint32_t     num_MPC_parties;  // N

  aimer_params_t params;

} aimer_instance_t;

const aimer_instance_t *aimer_instance_get(aimer_params_t param);

#endif // AIMER_INSTANCES_H
