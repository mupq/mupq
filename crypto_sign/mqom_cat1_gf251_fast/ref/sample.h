#ifndef MQOM_SAMPLE_H
#define MQOM_SAMPLE_H

#include "xof.h"
#include "rnd.h"

// Common API
//  for PRG and XOF

typedef struct samplable_t {
    unsigned int type;
    void* ctx;
} samplable_t;

#define SAMPLE_FROM_XOF 1
#define SAMPLE_FROM_PRG 2

static inline samplable_t xof_to_samplable(xof_context* ctx) {
    samplable_t x;
    x.type = SAMPLE_FROM_XOF;
    x.ctx = ctx;
    return x;
}
static inline samplable_t prg_to_samplable(prg_context* ctx) {
    samplable_t x;
    x.type = SAMPLE_FROM_PRG;
    x.ctx = ctx;
    return x;
}
static inline void byte_sample(samplable_t* smpl_ctx, uint8_t* buffer, uint32_t byte_len) {
    if(smpl_ctx->type == SAMPLE_FROM_XOF)
        xof_squeeze((xof_context*) smpl_ctx->ctx, buffer, byte_len);
    else if(smpl_ctx->type == SAMPLE_FROM_PRG)
        prg_sample((prg_context*) smpl_ctx->ctx, buffer, byte_len);
}


typedef struct samplable_x4_t {
    unsigned int type;
    void* ctx;
} samplable_x4_t;

static inline samplable_x4_t xof_to_samplable_x4(xof_context_x4* ctx) {
    samplable_x4_t x;
    x.type = SAMPLE_FROM_XOF;
    x.ctx = ctx;
    return x;
}
static inline samplable_x4_t prg_to_samplable_x4(prg_context_x4* ctx) {
    samplable_x4_t x;
    x.type = SAMPLE_FROM_PRG;
    x.ctx = ctx;
    return x;
}
static inline void byte_sample_x4(samplable_x4_t* smpl_ctx, uint8_t* const* buffer, uint32_t byte_len) {
    if(smpl_ctx->type == SAMPLE_FROM_XOF)
        xof_squeeze_x4((xof_context_x4*) smpl_ctx->ctx, buffer, byte_len);
    else if(smpl_ctx->type == SAMPLE_FROM_PRG)
        prg_sample_x4((prg_context_x4*) smpl_ctx->ctx, buffer, byte_len);
}

#endif /* MQOM_SAMPLE_H */
