#include "xof.h"
#include "parameters.h"

void xof_init(xof_context* ctx) {
    // Keccak_HashInitialize_SHAKE_CUR(ctx);
    #if PARAM_SECURITY == 128
    shake128_inc_init(ctx);
    #else
    shake256_inc_init(ctx);
    #endif
}
void xof_update(xof_context* ctx, const uint8_t* data, size_t byte_len) {
    // Keccak_HashUpdate(ctx, data, byte_len << 3);
    #if PARAM_SECURITY == 128
    shake128_inc_absorb(ctx, data, byte_len);
    #else
    shake256_inc_absorb(ctx, data, byte_len);
    #endif
}
void xof_final(xof_context* ctx) {
    // Keccak_HashFinal(ctx, NULL);
    #if PARAM_SECURITY == 128
    shake128_inc_finalize(ctx);
    #else
    shake256_inc_finalize(ctx);
    #endif
}
void xof_squeeze(xof_context* ctx, uint8_t* buffer, uint32_t byte_len) {
    // Keccak_HashSqueeze(ctx, buffer, byte_len << 3);
    #if PARAM_SECURITY == 128
    shake128_inc_squeeze(buffer, byte_len, ctx);
    #else
    shake256_inc_squeeze(buffer, byte_len, ctx);
    #endif
}
void xof_update_uint16_le(xof_context* ctx, uint16_t data) {
    const uint8_t data_le[2] = {data&0xFF, data>>8};
    xof_update(ctx, data_le, sizeof(data_le));
}
void xof_init_prefix(xof_context* ctx, const uint8_t prefix) {
    xof_init(ctx);
    xof_update(ctx, &prefix, sizeof(prefix));
}

#ifdef XOFX4
void xof_init_x4(xof_context_x4* ctx) {
    Keccak_HashInitializetimes4_SHAKE_CUR(ctx);
}
void xof_update_x4(xof_context_x4* ctx, uint8_t const* const* data, size_t byte_len) {
    // We remove the "const" at the first layer of pointers,
    //   but Keccak_HashUpdatetimes4 does not change this layer.
    Keccak_HashUpdatetimes4(ctx, (uint8_t const**) data, byte_len << 3);
}
void xof_update_x4_4(xof_context_x4* ctx, const uint8_t* data0, const uint8_t* data1,
                                const uint8_t* data2, const uint8_t* data3, size_t byte_len) {
    const uint8_t* data[4] = { data0, data1, data2, data3 };
    xof_update_x4(ctx, data, byte_len);
}
void xof_update_x4_1(xof_context_x4* ctx, const uint8_t* data, size_t byte_len) {
    const uint8_t* tmp[4] = { data, data, data, data };
    xof_update_x4(ctx, tmp, byte_len);
}
void xof_final_x4(xof_context_x4* ctx) {
    Keccak_HashFinaltimes4(ctx, NULL);
}
void xof_squeeze_x4(xof_context_x4* ctx, uint8_t* const* buffer, uint32_t buflen) {
    // We remove the "const" at the first layer of pointers,
    //   but Keccak_HashSqueezetimes4 does not change this layer.
    Keccak_HashSqueezetimes4(ctx, (uint8_t**) buffer, buflen << 3);
}
void xof_init_prefix_x4(xof_context_x4* ctx, uint8_t prefix) {
    xof_init_x4(ctx);
    xof_update_x4_1(ctx, &prefix, sizeof(prefix));
}

#else
void xof_init_x4(xof_context_x4* ctx) {
    for(unsigned int i = 0; i < 4; ++i)
        xof_init(&ctx->instances[i]);
}
void xof_update_x4(xof_context_x4* ctx, uint8_t const* const* data, size_t byte_len) {
    for(unsigned int i = 0; i < 4; ++i)
        xof_update(&ctx->instances[i], data[i], byte_len);
}
void xof_update_x4_4(xof_context_x4* ctx, const uint8_t* data0, const uint8_t* data1,
                                const uint8_t* data2, const uint8_t* data3, size_t byte_len) {
    xof_update(&ctx->instances[0], data0, byte_len);
    xof_update(&ctx->instances[1], data1, byte_len);
    xof_update(&ctx->instances[2], data2, byte_len);
    xof_update(&ctx->instances[3], data3, byte_len);
}
void xof_update_x4_1(xof_context_x4* ctx, const uint8_t* data, size_t byte_len) {
    for (unsigned int i = 0; i < 4; ++i)
        xof_update(&ctx->instances[i], data, byte_len);
}
void xof_final_x4(xof_context_x4* ctx) {
    for(unsigned int i = 0; i < 4; ++i)
        xof_final(&ctx->instances[i]);
}
void xof_squeeze_x4(xof_context_x4* ctx, uint8_t* const* buffer, uint32_t buflen) {
    for (unsigned int i = 0; i < 4; ++i)
        xof_squeeze(&ctx->instances[i], buffer[i], buflen);
}
void xof_init_prefix_x4(xof_context_x4* ctx, const uint8_t prefix) {
    for (unsigned int i = 0; i < 4; ++i)
        xof_init_prefix(&ctx->instances[i], prefix);
}
#endif

void xof_update_x4_uint16_le(xof_context_x4* ctx, uint16_t data) {
    const uint8_t data_le[2] = {data&0xFF, data>>8};
    xof_update_x4_1(ctx, data_le, sizeof(data_le));
}

void xof_update_x4_uint16s_le(xof_context_x4* ctx, const uint16_t data[4]) {
    const uint8_t data0_le[2] = {data[0]&0xFF, data[0]>>8};
    const uint8_t data1_le[2] = {data[1]&0xFF, data[1]>>8};
    const uint8_t data2_le[2] = {data[2]&0xFF, data[2]>>8};
    const uint8_t data3_le[2] = {data[3]&0xFF, data[3]>>8};
    xof_update_x4_4(ctx, (const uint8_t*)&data0_le, (const uint8_t*)&data1_le,
                   (const uint8_t*)&data2_le, (const uint8_t*)&data3_le, sizeof(data[0]));
}

