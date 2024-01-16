// -----------------------------------------------------------------------------
// File Name   : hash.h
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#ifndef HASH_H
#define HASH_H

#include <stdint.h>
#include "fips202.h"
#include "aimer.h"

static const uint8_t HASH_PREFIX_1    =  1;
static const uint8_t HASH_PREFIX_2    =  2;
static const uint8_t HASH_PREFIX_3    =  3;
static const uint8_t HASH_PREFIX_4    =  4;

#if _AIMER_L == 1
typedef shake128incctx hash_instance;
#elif _AIMER_L == 3 || _AIMER_L == 5
typedef shake256incctx hash_instance;
#else
#error "hash.[ch] only supports L in {1,3,5}."
#endif


void hash_init(hash_instance* ctx, size_t digest_size);
void hash_update(hash_instance* ctx, const uint8_t* data, size_t data_byte_len);
void hash_final(hash_instance* ctx);

void hash_squeeze(hash_instance* ctx, uint8_t* buffer, size_t buffer_len);
void hash_update_uint16(hash_instance* ctx, uint16_t data);
void hash_init_prefix(hash_instance* ctx, size_t digest_size,
                      const uint8_t prefix);

// x4 parallel hashing
typedef struct {
    hash_instance v[4];
} hash_instance_x4;

void hash_init_x4(hash_instance_x4* ctx, size_t digest_size);
void hash_update_x4(hash_instance_x4* ctx, const uint8_t** data, size_t data_byte_len);
void hash_update_x4_4(hash_instance_x4* ctx, const uint8_t* data0,
                      const uint8_t* data1, const uint8_t* data2, 
                      const uint8_t* data3, size_t data_byte_len);
void hash_update_x4_1(hash_instance_x4* ctx, const uint8_t *data,
                      size_t data_byte_len);
void hash_final_x4(hash_instance_x4* ctx);

void hash_init_prefix_x4(hash_instance_x4* ctx, size_t digest_size,
                         const uint8_t prefix);
void hash_squeeze_x4(hash_instance_x4* ctx, uint8_t** buffer, size_t buffer_len);
void hash_squeeze_x4_4(hash_instance_x4* ctx, uint8_t* buffer0,
                       uint8_t* buffer1, uint8_t* buffer2,
                       uint8_t* buffer3, size_t buffer_len);

void hash_update_x4_uint16(hash_instance_x4* ctx, uint16_t data);
void hash_update_x4_uint16s(hash_instance_x4* ctx, uint16_t* data);

#endif
