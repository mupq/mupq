// -----------------------------------------------------------------------------
// File Name   : hash.c
// Description : 
// SPDX-License-Identifier: MIT
// -----------------------------------------------------------------------------

#include "hash.h"
#include "portable_endian.h"

void hash_init(hash_instance* ctx, size_t digest_size)
{
  // if (digest_size <= 32)
  // {
  //   Keccak_HashInitialize_SHAKE128(ctx);
  // }
  // else
  // {
  //   Keccak_HashInitialize_SHAKE256(ctx);
  // }
  (void)digest_size;
  #if _AIMER_L == 1
  shake128_inc_init(ctx);
  #elif _AIMER_L == 3 || _AIMER_L == 5
  shake256_inc_init(ctx);
  #else
  #error "hash.[ch] only supports L in {1,3,5}."
  #endif
}

void hash_update(hash_instance* ctx, const uint8_t* data, size_t data_byte_len)
{
  //Keccak_HashUpdate(ctx, data, data_byte_len << 3);
  #if _AIMER_L == 1
  shake128_inc_absorb(ctx, data, data_byte_len);
  #elif _AIMER_L == 3 || _AIMER_L == 5
  shake256_inc_absorb(ctx, data, data_byte_len);
  #else
  #error "hash.[ch] only supports L in {1,3,5}."
  #endif
}

void hash_final(hash_instance* ctx)
{
  //Keccak_HashFinal(ctx, NULL);
  #if _AIMER_L == 1
  shake128_inc_finalize(ctx);
  #elif _AIMER_L == 3 || _AIMER_L == 5
  shake256_inc_finalize(ctx);
  #else
  #error "hash.[ch] only supports L in {1,3,5}."
  #endif
}

void hash_squeeze(hash_instance* ctx, uint8_t* buffer, size_t buffer_len)
{
  //Keccak_HashSqueeze(ctx, buffer, buffer_len << 3);
  #if _AIMER_L == 1
  shake128_inc_squeeze(buffer, buffer_len, ctx);
  #elif _AIMER_L == 3 || _AIMER_L == 5
  shake256_inc_squeeze(buffer, buffer_len, ctx);
  #else
  #error "hash.[ch] only supports L in {1,3,5}."
  #endif
}

void hash_update_uint16(hash_instance* ctx, uint16_t data)
{
  const uint16_t data_little_endian = htole16(data);
  hash_update(ctx, (const uint8_t*)&data_little_endian,
              sizeof(data_little_endian));
}

void hash_init_prefix(hash_instance* ctx, size_t digest_size,
                      const uint8_t prefix)
{
  hash_init(ctx, digest_size);
  hash_update(ctx, &prefix, sizeof(prefix));
}

// x4 parallel hashing
void hash_init_x4(hash_instance_x4* ctx, size_t digest_size)
{
  //if (digest_size <= 32)
  //{
  //  Keccak_HashInitializetimes4_SHAKE128(ctx);
  //}
  //else
  //{
  //  Keccak_HashInitializetimes4_SHAKE256(ctx);
  //}
  hash_init(&ctx->v[0], digest_size);
  hash_init(&ctx->v[1], digest_size);
  hash_init(&ctx->v[2], digest_size);
  hash_init(&ctx->v[3], digest_size);
}

void hash_update_x4(hash_instance_x4* ctx, const uint8_t** data, size_t data_byte_len)
{
  //Keccak_HashUpdatetimes4(ctx, data, data_byte_len << 3);
  hash_update(&ctx->v[0], data[0], data_byte_len);
  hash_update(&ctx->v[1], data[1], data_byte_len);
  hash_update(&ctx->v[2], data[2], data_byte_len);
  hash_update(&ctx->v[3], data[3], data_byte_len);
}

void hash_update_x4_4(hash_instance_x4* ctx, const uint8_t* data0,
                      const uint8_t* data1, const uint8_t* data2, 
                      const uint8_t* data3, size_t data_byte_len)
{
  //const uint8_t* data[4] = {data0, data1, data2, data3};
  //Keccak_HashUpdatetimes4(ctx, data, data_byte_len << 3);
  hash_update(&ctx->v[0], data0, data_byte_len);
  hash_update(&ctx->v[1], data1, data_byte_len);
  hash_update(&ctx->v[2], data2, data_byte_len);
  hash_update(&ctx->v[3], data3, data_byte_len);
}

void hash_update_x4_1(hash_instance_x4* ctx, const uint8_t *data,
                      size_t data_byte_len)
{
  //const uint8_t* temp[4] = {data, data, data, data};
  //Keccak_HashUpdatetimes4(ctx, temp, data_byte_len << 3);
  hash_update(&ctx->v[0], data, data_byte_len);
  hash_update(&ctx->v[1], data, data_byte_len);
  hash_update(&ctx->v[1], data, data_byte_len);
  hash_update(&ctx->v[1], data, data_byte_len);
}

void hash_init_prefix_x4(hash_instance_x4* ctx, size_t digest_size,
                         const uint8_t prefix)
{
  hash_init_x4(ctx, digest_size);
  hash_update_x4_1(ctx, &prefix, sizeof(prefix));
}

void hash_final_x4(hash_instance_x4* ctx)
{
  // Keccak_HashFinaltimes4(ctx, NULL);
  hash_final(&ctx->v[0]);
  hash_final(&ctx->v[1]);
  hash_final(&ctx->v[2]);
  hash_final(&ctx->v[3]);
}

void hash_squeeze_x4(hash_instance_x4* ctx, uint8_t** buffer, size_t buffer_len)
{
  //Keccak_HashSqueezetimes4(ctx, buffer, buffer_len << 3);
  hash_squeeze(&ctx->v[0], buffer[0], buffer_len);
  hash_squeeze(&ctx->v[1], buffer[1], buffer_len);
  hash_squeeze(&ctx->v[2], buffer[2], buffer_len);
  hash_squeeze(&ctx->v[3], buffer[3], buffer_len);
}

void hash_squeeze_x4_4(hash_instance_x4* ctx, uint8_t* buffer0,
                       uint8_t* buffer1, uint8_t* buffer2,
                       uint8_t* buffer3, size_t buffer_len)
{
  //uint8_t *buffer[4];
  //buffer[0] = buffer0; buffer[1] = buffer1;
  //buffer[2] = buffer2; buffer[3] = buffer3;
  //Keccak_HashSqueezetimes4(ctx, buffer, buffer_len << 3);
  hash_squeeze(&ctx->v[0], buffer0, buffer_len);
  hash_squeeze(&ctx->v[1], buffer1, buffer_len);
  hash_squeeze(&ctx->v[2], buffer2, buffer_len);
  hash_squeeze(&ctx->v[3], buffer3, buffer_len);
}

void hash_update_x4_uint16(hash_instance_x4* ctx, uint16_t data)
{
  const uint16_t data_little_endian = htole16(data);
  hash_update_x4_1(ctx, (const uint8_t*)&data_little_endian,
                   sizeof(data_little_endian));
}

void hash_update_x4_uint16s(hash_instance_x4* ctx, uint16_t* data)
{
  uint16_t data0_little_endian = htole16(data[0]);
  uint16_t data1_little_endian = htole16(data[1]);
  uint16_t data2_little_endian = htole16(data[2]);
  uint16_t data3_little_endian = htole16(data[3]);
  hash_update_x4_4(ctx, (const uint8_t*)&data0_little_endian,
                  (const uint8_t*)&data1_little_endian, (const uint8_t*)&data2_little_endian,
                  (const uint8_t*)&data3_little_endian, sizeof(data[0]));
}
