// -----------------------------------------------------------------------------
// File Name   : aim.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef AIM_H
#define AIM_H

#include "aimer.h"

#if       _AIMER_L == 1
  #include "aim128.h"

#elif     _AIMER_L == 3
  #include "aim192.h"

#elif     _AIMER_L == 5
  #include "aim256.h"

#endif

typedef struct aim_params_t
{
  uint32_t block_size;
  uint32_t num_input_sboxes;
} aim_params_t;

#endif // AIM_H
