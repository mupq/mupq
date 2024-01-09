#include "rnd.h"

void prg_init(prg_context* ctx, const uint8_t* seed, const uint8_t* salt) {
    xof_init(ctx);
    if(salt != NULL)
        xof_update(ctx, salt, PARAM_SALT_SIZE);
    xof_update(ctx, seed, PARAM_SEED_SIZE);
    xof_final(ctx);
}

void prg_init_with_seed_material(prg_context* ctx, const uint8_t* seed_material, unsigned int seed_material_size) {
    xof_init(ctx);
    xof_update(ctx, seed_material, seed_material_size);
    xof_final(ctx);
}

void prg_sample(prg_context* ctx, uint8_t* out, uint32_t outlen) {
    xof_squeeze(ctx, out, outlen);
}

void prg_init_x4_array(prg_context ctx[4], uint8_t const* const* seed, const uint8_t* salt) {
    for(unsigned int i = 0; i < 4; ++i)
        prg_init(&ctx[i], seed[i], salt);
}

#ifdef PRGX4
void prg_init_x4(prg_context_x4* ctx,  uint8_t const* const* seed, const uint8_t* salt) {
    xof_init_x4(ctx);
    if(salt != NULL)
        xof_update_x4_1(ctx, salt, PARAM_SALT_SIZE);
    xof_update_x4(ctx, seed, PARAM_SEED_SIZE);
    xof_final_x4(ctx);
}
void prg_init_with_seed_material_x4(prg_context_x4* ctx, uint8_t const* const* seed_material, unsigned int seed_material_size) {
    xof_init_x4(ctx);
    xof_update_x4(ctx, seed_material, seed_material_size);
    xof_final_x4(ctx);
}
void prg_sample_x4(prg_context_x4* ctx, uint8_t* const* out, uint32_t outlen) {
    xof_squeeze_x4(ctx, out, outlen);
}

#else
void prg_init_x4(prg_context_x4* ctx,  uint8_t const* const* seed, const uint8_t* salt) {
    for(unsigned int i = 0; i < 4; ++i)
        prg_init(&ctx->instances[i], seed[i], salt);
}
void prg_init_with_seed_material_x4(prg_context_x4* ctx, uint8_t const* const* seed_material, unsigned int seed_material_size) {
    for(unsigned int i = 0; i < 4; ++i)
        prg_init_with_seed_material(&ctx->instances[i], seed_material[i], seed_material_size);
}
void prg_sample_x4(prg_context_x4* ctx, uint8_t* const* out, uint32_t outlen) {
    for (unsigned int i = 0; i < 4; ++i)
        prg_sample(&ctx->instances[i], out[i], outlen);
}
#endif
