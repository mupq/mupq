#include "hash.h"
#include "parameters.h"

void hash_init(hash_context* ctx) {
    // Keccak_HashInitialize_SHA3_CUR(ctx);
    #if PARAM_SECURITY == 128
    sha3_256_inc_init(ctx);
    #elif PARAM_SECURITY == 192
    sha3_384_inc_init(ctx);
    #elif PARAM_SECURITY == 256
    sha3_512_inc_init(ctx);
    #else
    #error "No hash implementation for this security level"
    #endif
}
void hash_update(hash_context* ctx, const uint8_t* data, size_t byte_len) {
    // Keccak_HashUpdate(ctx, data, byte_len << 3);
    #if PARAM_SECURITY == 128
    sha3_256_inc_absorb(ctx, data, byte_len);
    #elif PARAM_SECURITY == 192
    sha3_384_inc_absorb(ctx, data, byte_len);
    #elif PARAM_SECURITY == 256
    sha3_512_inc_absorb(ctx, data, byte_len);
    #else
    #error "No hash implementation for this security level"
    #endif
}
void hash_final(hash_context* ctx, uint8_t* digest) {
    // Keccak_HashFinal(ctx, digest);
    #if PARAM_SECURITY == 128
    sha3_256_inc_finalize(digest, ctx);
    #elif PARAM_SECURITY == 192
    sha3_384_inc_finalize(digest, ctx);
    #elif PARAM_SECURITY == 256
    sha3_512_inc_finalize(digest, ctx);
    #else
    #error "No hash implementation for this security level"
    #endif
}
void hash_update_uint16_le(hash_context* ctx, uint16_t data) {
    const uint8_t data_le[2] = {data&0xFF, data>>8};
    hash_update(ctx, data_le, sizeof(data_le));
}
void hash_init_prefix(hash_context* ctx, const uint8_t prefix) {
    hash_init(ctx);
    if(prefix != HASH_PREFIX_NONE)
        hash_update(ctx, &prefix, sizeof(prefix));
}

#ifdef HASHX4
void hash_init_x4(hash_context_x4* ctx) {
    Keccak_HashInitializetimes4_SHA3_CUR(ctx);
}
void hash_update_x4(hash_context_x4* ctx, uint8_t const* const* data, size_t byte_len) {
    // We remove the "const" at the first layer of pointers,
    //   but Keccak_HashUpdatetimes4 does not change this layer.
    Keccak_HashUpdatetimes4(ctx, (uint8_t const**) data, byte_len << 3);
}
void hash_update_x4_4(hash_context_x4* ctx, const uint8_t* data0, const uint8_t* data1,
                                const uint8_t* data2, const uint8_t* data3, size_t byte_len) {
    const uint8_t* data[4] = { data0, data1, data2, data3 };
    hash_update_x4(ctx, data, byte_len);
}
void hash_update_x4_1(hash_context_x4* ctx, const uint8_t* data, size_t byte_len) {
    const uint8_t* tmp[4] = { data, data, data, data };
    hash_update_x4(ctx, tmp, byte_len);
}
void hash_final_x4(hash_context_x4* ctx, uint8_t* const* digest) {
    // We remove the "const" at the first layer of pointers,
    //   but Keccak_HashSqueezetimes4 does not change this layer.
    Keccak_HashFinaltimes4(ctx, (uint8_t**) digest);
}
void hash_init_prefix_x4(hash_context_x4* ctx, uint8_t prefix) {
    hash_init_x4(ctx);
    if(prefix != HASH_PREFIX_NONE)
        hash_update_x4_1(ctx, &prefix, sizeof(prefix));
}

#else
void hash_init_x4(hash_context_x4* ctx) {
    for(unsigned int i = 0; i < 4; ++i)
        hash_init(&ctx->instances[i]);
}
void hash_update_x4(hash_context_x4* ctx, uint8_t const* const* data, size_t byte_len) {
    for(unsigned int i = 0; i < 4; ++i)
        hash_update(&ctx->instances[i], data[i], byte_len);
}
void hash_update_x4_4(hash_context_x4* ctx, const uint8_t* data0, const uint8_t* data1,
                                const uint8_t* data2, const uint8_t* data3, size_t byte_len) {
    hash_update(&ctx->instances[0], data0, byte_len);
    hash_update(&ctx->instances[1], data1, byte_len);
    hash_update(&ctx->instances[2], data2, byte_len);
    hash_update(&ctx->instances[3], data3, byte_len);
}
void hash_update_x4_1(hash_context_x4* ctx, const uint8_t* data, size_t byte_len) {
    for (unsigned int i = 0; i < 4; ++i)
        hash_update(&ctx->instances[i], data, byte_len);
}
void hash_final_x4(hash_context_x4* ctx, uint8_t* const* digest) {
    for(unsigned int i = 0; i < 4; ++i)
        hash_final(&ctx->instances[i], digest[i]);
}
void hash_init_prefix_x4(hash_context_x4* ctx, const uint8_t prefix) {
    for (unsigned int i = 0; i < 4; ++i)
        hash_init_prefix(&ctx->instances[i], prefix);
}
#endif

void hash_update_x4_uint16_le(hash_context_x4* ctx, uint16_t data) {
    const uint8_t data_le[2] = {data&0xFF, data>>8};
    hash_update_x4_1(ctx, data_le, sizeof(data_le));
}

void hash_update_x4_uint16s_le(hash_context_x4* ctx, const uint16_t data[4]) {
    const uint8_t data0_le[2] = {data[0]&0xFF, data[0]>>8};
    const uint8_t data1_le[2] = {data[1]&0xFF, data[1]>>8};
    const uint8_t data2_le[2] = {data[2]&0xFF, data[2]>>8};
    const uint8_t data3_le[2] = {data[3]&0xFF, data[3]>>8};
    hash_update_x4_4(ctx, (const uint8_t*)&data0_le, (const uint8_t*)&data1_le,
                   (const uint8_t*)&data2_le, (const uint8_t*)&data3_le, sizeof(data[0]));
}
