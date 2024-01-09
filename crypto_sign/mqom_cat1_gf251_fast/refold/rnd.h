#ifndef MQOM_RND_H
#define MQOM_RND_H

#include <stdint.h>
#include "parameters.h"

extern int randombytes(unsigned char* x, unsigned long long xlen);
#define sample_seed(seed) randombytes(seed, PARAM_SEED_SIZE)
#define sample_salt(seed) randombytes(seed, PARAM_SALT_SIZE)


#include "xof.h"
typedef xof_context prg_context;

#ifdef PRGX4
typedef xof_context_x4 prg_context_x4;
#endif



void prg_init(prg_context* ctx, const uint8_t* seed, const uint8_t* salt);
void prg_init_with_seed_material(prg_context* ctx, const uint8_t* seed_material, unsigned int seed_material_size);
void prg_sample(prg_context* ctx, uint8_t* out, uint32_t outlen);
void prg_init_x4_array(prg_context ctx[4], uint8_t const* const* seed, const uint8_t* salt);

#ifndef PRGX4
typedef struct prg_context_x4 {
  prg_context instances[4];
} prg_context_x4;
#endif
void prg_init_x4(prg_context_x4* ctx,  uint8_t const* const* seed, const uint8_t* salt);
void prg_init_with_seed_material_x4(prg_context_x4* ctx, uint8_t const* const* seed_material, unsigned int seed_material_size);
void prg_sample_x4(prg_context_x4* ctx, uint8_t* const* out, uint32_t outlen);


#endif /* MQOM_RND_H */
