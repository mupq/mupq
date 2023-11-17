#include <string.h>

#include "batch_tools.h"

static void
memcpyX (uintX_t *dest, const uintX_t *src, size_t n)
{
  size_t i;
  for (i = 0; i < n; i++)
    {
      dest[i] = src[i];
    }
}

#if !defined(FIELD_SIZE) || IS_POW2 (FIELD_SIZE)

static void
memxorX (uintX_t *dest, const uintX_t *src, size_t n)
{
  size_t i;
  for (i = 0; i < n; i++)
    {
      dest[i] ^= src[i];
    }
}

#endif

#if !defined(FIELD_SIZE) || (FIELD_SIZE == 2)

static void
memandX (uintX_t *dest, const uintX_t *src, size_t n)
{
  size_t i;
  for (i = 0; i < n; i++)
    {
      dest[i] &= src[i];
    }
}

#ifdef GF2K_ARITH_USE_TABLES
#define P2(n) n, n^1, n^1, n
#define P4(n) P2(n), P2(n^1), P2(n^1), P2(n)
#define P6(n) P4(n), P4(n^1), P4(n^1), P4(n)
static const uint8_t parity_table[256] =
  {
    P6(0), P6(1), P6(1), P6(0)
  };
#undef P6
#undef P4
#undef P2
#endif

static uint8_t
xorbitsX (uintX_t *buf, size_t n)
{
  union vecX_t u;
  while (n > 1)
    {
      memxorX (buf, buf + (n + 1) / 2, n / 2);
      n = (n + 1) / 2;
    }
  u.tX = buf[0];
#if UINTX_BITSIZE >= 256
  u.t128[0] ^= u.t128[1];
#endif
#if UINTX_BITSIZE >= 128
  u.t64[0] ^= u.t64[1];
#endif
#if UINTX_BITSIZE >= 64
  u.t32[0] ^= u.t32[1];
#endif
#if UINTX_BITSIZE >= 32
  u.t16[0] ^= u.t16[1];
#endif
#if UINTX_BITSIZE >= 16
  u.t8[0] ^= u.t8[1];
#endif
#ifdef GF2K_ARITH_USE_TABLES
  return parity_table[u.t8[0]];
#else
  u.t8[0] ^= u.t8[0] >> 4;
  u.t8[0] ^= u.t8[0] >> 2;
  u.t8[0] ^= u.t8[0] >> 1;
  return u.t8[0] & 0x1;
#endif
}

#endif

#if !defined(FIELD_SIZE) || (FIELD_SIZE == 4)

#ifdef GF2K_ARITH_USE_TABLES
static const uint8_t parity2_table[256] = {
0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0, 1, 0, 3, 2, 0, 1, 2, 3, 3, 2, 1, 0, 2, 3, 0, 1, 2, 3, 0, 1, 3, 2, 1, 0, 0, 1, 2, 3, 1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2, 0, 1, 2, 3, 1, 0, 3, 2, 0, 1, 2, 3, 3, 2, 1, 0, 2, 3, 0, 1, 0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0, 3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2, 0, 1, 2, 3, 2, 3, 0, 1, 3, 2, 1, 0, 0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0, 0, 1, 2, 3, 1, 0, 3, 2, 3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0, 1, 0, 3, 2, 0, 1, 2, 3, 3, 2, 1, 0, 2, 3, 0, 1, 3, 2, 1, 0, 2, 3, 0, 1, 1, 0, 3, 2, 0, 1, 2, 3, 2, 3, 0, 1, 3, 2, 1, 0, 0, 1, 2, 3, 1, 0, 3, 2, 1, 0, 3, 2, 0, 1, 2, 3, 3, 2, 1, 0, 2, 3, 0, 1, 0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0
};
#endif

static uint8_t
xorbits2X (uintX_t *buf, size_t n)
{
  union vecX_t u;
  while (n > 1)
    {
      memxorX (buf, buf + (n + 1) / 2, n / 2);
      n = (n + 1) / 2;
    }
  u.tX = buf[0];
#if UINTX_BITSIZE >= 256
  u.t128[0] ^= u.t128[1];
#endif
#if UINTX_BITSIZE >= 128
  u.t64[0] ^= u.t64[1];
#endif
#if UINTX_BITSIZE >= 64
  u.t32[0] ^= u.t32[1];
#endif
#if UINTX_BITSIZE >= 32
  u.t16[0] ^= u.t16[1];
#endif
#if UINTX_BITSIZE >= 16
  u.t8[0] ^= u.t8[1];
#endif
#ifdef GF2K_ARITH_USE_TABLES
  return parity2_table[u.t8[0]];
#else
  u.t8[0] ^= u.t8[0] >> 4;
  u.t8[0] ^= u.t8[0] >> 2;
  return u.t8[0] & 0x3;
#endif
}

static uint8_t
mul_gf2_2 (uint8_t x, uint8_t y)
{
#ifdef GF2K_ARITH_USE_TABLES
#include "mul_gf2_2.inc"
  return mul_tab[x][y];
#else
  int i;
  uint8_t b;
  uint8_t z = 0;
  for (i = 0; i < 2; i++)
    {
      z = z ^ (x & ((y & 0x55) ^ (y & 0x55) << 1));
      b = (x >> 1);
      x = x << 1;
      x = x ^ (0x7 & -(b & 0x1)) ^ (0x1c & -(b & 0x4))
        ^ (0x70 & -(b & 0x10)) ^ (0x1c0 & -(b & 0x40));
      y = y >> 1;
    }
  return z;
#endif
}

#endif

#if !defined(FIELD_SIZE) || (FIELD_SIZE == 16)

#ifdef GF2K_ARITH_USE_TABLES
static const uint8_t parity4_table[256] = {
0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13, 3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8, 15, 14, 13, 12, 4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 5, 4, 7, 6, 1, 0, 3, 2, 13, 12, 15, 14, 9, 8, 11, 10, 6, 7, 4, 5, 2, 3, 0, 1, 14, 15, 12, 13, 10, 11, 8, 9, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 13, 12, 15, 14, 1, 0, 3, 2, 5, 4, 7, 6, 10, 11, 8, 9, 14, 15, 12, 13, 2, 3, 0, 1, 6, 7, 4, 5, 11, 10, 9, 8, 15, 14, 13, 12, 3, 2, 1, 0, 7, 6, 5, 4, 12, 13, 14, 15, 8, 9, 10, 11, 4, 5, 6, 7, 0, 1, 2, 3, 13, 12, 15, 14, 9, 8, 11, 10, 5, 4, 7, 6, 1, 0, 3, 2, 14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
};
#endif

static uint8_t
xornibblesX (uintX_t *buf, size_t n)
{
  union vecX_t u;
  while (n > 1)
    {
      memxorX (buf, buf + (n + 1) / 2, n / 2);
      n = (n + 1) / 2;
    }
  u.tX = buf[0];
#if UINTX_BITSIZE >= 256
  u.t128[0] ^= u.t128[1];
#endif
#if UINTX_BITSIZE >= 128
  u.t64[0] ^= u.t64[1];
#endif
#if UINTX_BITSIZE >= 64
  u.t32[0] ^= u.t32[1];
#endif
#if UINTX_BITSIZE >= 32
  u.t16[0] ^= u.t16[1];
#endif
#if UINTX_BITSIZE >= 16
  u.t8[0] ^= u.t8[1];
#endif
#ifdef GF2K_ARITH_USE_TABLES
  return parity4_table[u.t8[0]];
#else
  u.t8[0] ^= u.t8[0] >> 4;
  return u.t8[0] & 0xf;
#endif
}

static uint8_t
mul_gf2_4 (uint8_t x, uint8_t y)
{
#ifdef GF2K_ARITH_USE_TABLES
#include "mul_gf2_4.inc"
  return mul_tab[x][y];
#else
  int i;
  uint8_t b;
  uint8_t z = 0;
  for (i = 0; i < 4; i++)
    {
      z = z ^ (x & ((-(y & 0x1) & 0xf) | -(y & 0x10)));
      b = (x >> 3);
      x = x << 1;
      x = x ^ (0x13 & -(b & 0x1)) ^ (0x130 & -(b & 0x10));
      y = y >> 1;
    }
  return z;
#endif
}

#endif

#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256)

static uint8_t
xorbytesX (uintX_t *buf, size_t n)
{
  union vecX_t u;
  while (n > 1)
    {
      memxorX (buf, buf + (n + 1) / 2, n / 2);
      n = (n + 1) / 2;
    }
  u.tX = buf[0];
#if UINTX_BITSIZE >= 256
  u.t128[0] ^= u.t128[1];
#endif
#if UINTX_BITSIZE >= 128
  u.t64[0] ^= u.t64[1];
#endif
#if UINTX_BITSIZE >= 64
  u.t32[0] ^= u.t32[1];
#endif
#if UINTX_BITSIZE >= 32
  u.t16[0] ^= u.t16[1];
#endif
#if UINTX_BITSIZE >= 16
  u.t8[0] ^= u.t8[1];
#endif
  return u.t8[0];
}

static uint8_t
mul_gf2_8 (uint8_t x, uint8_t y)
{
#ifdef GF2K_ARITH_USE_TABLES
#include "mul_gf2_8.inc"
  return mul_tab[x][y];
#else
  int i;
  uint8_t b;
  uint8_t z = 0;
  for (i = 0; i < 8; i++)
    {
      z = z ^ (x & -(y & 0x1));
      b = x >> 7;
      x = x << 1;
      x = x ^ (0x1d & -b);
      y = y >> 1;
    }
  return z;
#endif
}

#endif

int
batch_getbitlen (int q, int n)
{
#ifndef FIELD_SIZE
  if (q == 2)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 2)
    {
      return n;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 4)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 4)
    {
      return n << 1;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 16)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 16)
    {
      return n << 2;
    }
#endif
#ifndef FIELD_SIZE
  else if ((q == 256) || (q == 251))
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256) || (FIELD_SIZE == 251)
    {
      return n << 3;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      return n << 4;
    }
#endif
  return 0;
}

int
batch_getlen (int q, int n)
{
  int n8 = 0;
#ifndef FIELD_SIZE
  if (q == 2)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 2)
    {
      n8 = (n + 7) >> 3;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 4)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 4)
    {
      n8 = (n + 3) >> 2;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 16)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 16)
    {
      n8 = (n + 1) >> 1;
    }
#endif
#ifndef FIELD_SIZE
  else if ((q == 256) || (q == 251))
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256) || (FIELD_SIZE == 251)
    {
      n8 = n;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      n8 = n << 1;
    }
#endif
  return n8;
}

int
batch_getlenX (int q, int n)
{
  return CONVX ((batch_getbitlen (q, n) + 7) >> 3);
}

void
batch_import (uintX_t *dest, const uint8_t *src, int q, int n, int k)
{
#ifndef FIELD_SIZE
  if ((q == 2) || (q == 4) || (q == 16))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 2) || (FIELD_SIZE == 4) || (FIELD_SIZE == 16)
    {
      int i;
      const int n1 = batch_getbitlen (q, n);
      const int k1 = batch_getbitlen (q, k);
      const int shift = (8 - k1) & 0x7;
      const int offset = (k1 + 7) >> 3;
      const int n8 = (n1 + 7) >> 3;
      uint8_t v;

      v = shift ? src[offset - 1] : 0;
      for (i = 0; i < n8; i++)
        {
          ((uint8_t *) dest)[i] = v >> (8 - shift);
          v = src[offset + i];
          ((uint8_t *) dest)[i] |= v << shift;
        }
      if (n1 & 0x7)
        {
          ((uint8_t *) dest)[i - 1] &= (1 << (n1 & 0x7)) - 1;
        }
      memset (((uint8_t *) dest) + n8, 0, CONVX8 (n8) - n8);
    }
#endif
#ifndef FIELD_SIZE
  else if ((q == 256) || (q == 251) || (q == 65521))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 256) || (FIELD_SIZE == 251) || (FIELD_SIZE == 65521)
    {
      const int o8 = batch_getlen (q, k);
      const int n8 = batch_getlen (q, n);
      memcpy (dest, src + o8, n8);
      memset (((uint8_t *) dest) + n8, 0, CONVX8 (n8) - n8);
    }
#endif
}

void
batch_export (uint8_t *dest, const uintX_t *src, int q, int n, int k)
{
#ifndef FIELD_SIZE
  if ((q == 2) || (q == 4) || (q == 16))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 2) || (FIELD_SIZE == 4) || (FIELD_SIZE == 16)
    {
      int i;
      const int n1 = batch_getbitlen (q, n);
      const int k1 = batch_getbitlen (q, k);
      const int shift = k1 & 0x7;
      const int offset = k1 >> 3;
      const int n8 = batch_getlen (q, n);
      uint8_t v, w;

      v = shift ? dest[offset] << (8 - shift) : 0;
      for (i = 0; i < n8; i++)
        {
          w = v >> (8 - shift);
          v = ((uint8_t *) src)[i];
          w |= v << shift;
          dest[offset + i] = w;
        }
      if (((n1 + shift + 7) >> 3) > n8)
        {
          dest[offset + i] = v >> (8 - shift);
        }
    }
#endif
#ifndef FIELD_SIZE
  else if ((q == 256) || (q == 251) || (q == 65521))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 256) || (FIELD_SIZE == 251) || (FIELD_SIZE == 65521)
    {
      const int o8 = batch_getlen (q, k);
      const int n8 = batch_getlen (q, n);
      memcpy (dest + o8, src, n8);
    }
#endif
}

void
batch_clear (uintX_t *dest, int q, int n)
{
  memset (dest, 0, CONVX8 (batch_getlen (q, n)));
}

void
batch_copy (uintX_t *dest, const uintX_t *src, int q, int n)
{
  memcpyX (dest, src, batch_getlenX (q, n));
}

void
batch_add (uintX_t *dest, const uintX_t *src, int q, int n)
{
#ifndef FIELD_SIZE
  if ((q == 2) || (q == 4) || (q == 16) || (q == 256))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 2) || (FIELD_SIZE == 4) || (FIELD_SIZE == 16) \
  || (FIELD_SIZE == 256)
    {
      memxorX (dest, src, batch_getlenX (q, n));
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 251)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 251)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint16_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = (v0 + v1) % 251;
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint32_t v0, v1;
          v0 = ((uint16_t *) dest)[i];
          v1 = ((uint16_t *) src)[i];
          v0 = (v0 + v1) % 65521;
          ((uint16_t *) dest)[i] = (uint16_t) v0;
        }
    }
#endif
}

void
batch_sub (uintX_t *dest, const uintX_t *src, int q, int n)
{
#ifndef FIELD_SIZE
  if ((q == 2) || (q == 4) || (q == 16) || (q == 256))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 2) || (FIELD_SIZE == 4) || (FIELD_SIZE == 16) \
  || (FIELD_SIZE == 256)
    {
      memxorX (dest, src, batch_getlenX (q, n));
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 251)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 251)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint16_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = (251 + v0 - v1) % 251;
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint32_t v0, v1;
          v0 = ((uint16_t *) dest)[i];
          v1 = ((uint16_t *) src)[i];
          v0 = (65521 + v0 - v1) % 65521;
          ((uint16_t *) dest)[i] = (uint16_t) v0;
        }
    }
#endif
}

void
batch_mul (uintX_t *dest, const uintX_t *src, int q, int n)
{
#ifndef FIELD_SIZE
  if (q == 2)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 2)
    {
      memandX (dest, src, batch_getlenX (q, n));
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 4)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 4)
    {
      const int n8 = batch_getlen (q, n);
      int i;
      for (i = 0; i < n8; i++)
        {
          uint8_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = mul_gf2_2 (v0, v1);
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 16)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 16)
    {
      const int n8 = batch_getlen (q, n);
      int i;
      for (i = 0; i < n8; i++)
        {
          uint8_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = mul_gf2_4 (v0, v1);
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 256)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint8_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = mul_gf2_8 (v0, v1);
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 251)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 251)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint16_t v0, v1;
          v0 = ((uint8_t *) dest)[i];
          v1 = ((uint8_t *) src)[i];
          v0 = (v0 * v1) % 251;
          ((uint8_t *) dest)[i] = (uint8_t) v0;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      int i;
      for (i = 0; i < n; i++)
        {
          uint32_t v0, v1;
          v0 = ((uint16_t *) dest)[i];
          v1 = ((uint16_t *) src)[i];
          v0 = (v0 * v1) % 65521;
          ((uint16_t *) dest)[i] = (uint16_t) v0;
        }
    }
#endif
}

void
batch_sum (uintX_t *dest, uintX_t *buf, int j, int q, int n)
{
#ifndef FIELD_SIZE
  if (q == 2)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 2)
    {
      const int nX = batch_getlenX (q, n);
      ((uint8_t *) dest)[j >> 3] ^= xorbitsX (buf, nX) << (j & 0x7);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 4)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 4)
    {
      const int nX = batch_getlenX (q, n);
      ((uint8_t *) dest)[j >> 2] ^= xorbits2X (buf, nX) << ((j & 0x3) << 1);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 16)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 16)
    {
      const int nX = batch_getlenX (q, n);
      ((uint8_t *) dest)[j >> 1] ^= xornibblesX (buf, nX) << ((j & 0x1) << 2);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 256)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256)
    {
      const int nX = batch_getlenX (q, n);
      ((uint8_t *) dest)[j] ^= xorbytesX (buf, nX);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 251)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 251)
    {
      int i;
      uint16_t v0, v1;
      v0 = 0;
      for (i = 0; i < n; i++)
        {
          v1 = ((uint8_t *) buf)[i];
          v0 = (v0 + v1) % 251;
        }
      v1 = ((uint8_t *) dest)[j];
      v0 = (v0 + v1) % 251;
      ((uint8_t *) dest)[j] = (uint8_t) v0;
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      int i;
      uint32_t v0, v1;
      v0 = 0;
      for (i = 0; i < n; i++)
        {
          v1 = ((uint16_t *) buf)[i];
          v0 = (v0 + v1) % 65521;
        }
      v1 = ((uint16_t *) dest)[j];
      v0 = (v0 + v1) % 65521;
      ((uint16_t *) dest)[j] = (uint16_t) v0;
    }
#endif
}

void
batch_generate (uintX_t *dest, int q, int n,
                void (*extract) (void *, size_t, void *), void *arg)
{
  uint8_t *dest8 = (uint8_t *) dest;
#ifndef FIELD_SIZE
  if ((q == 2) || (q == 4) || (q == 16))
#endif
#if !defined(FIELD_SIZE) \
  || (FIELD_SIZE == 2) || (FIELD_SIZE == 4) || (FIELD_SIZE == 16)
    {
      const int n8 = batch_getlen (q, n);
      const int rem_bits = batch_getbitlen (q, n) & 0x7;
      extract (dest8, n8, arg);
      memset (dest8 + n8, 0, CONVX8 (n8) - n8);
      if (rem_bits)
        {
          dest8[n8 - 1] &= (1 << rem_bits) - 1;
        }
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 256)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 256)
    {
      const int n8 = batch_getlen (q, n);
      extract (dest8, n8, arg);
      memset (dest8 + n8, 0, CONVX8 (n8) - n8);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 251)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 251)
    {
      int i = 0;
      const int n8 = batch_getlen (q, n);
      while (i < n8)
        {
          uint8_t v0;
          extract (dest8, 1, arg);
          v0 = dest8[0];
          if (v0 < 251)
            {
              i++;
              dest8++;
            }
        }
      memset (dest8, 0, CONVX8 (n8) - n8);
    }
#endif
#ifndef FIELD_SIZE
  else if (q == 65521)
#endif
#if !defined(FIELD_SIZE) || (FIELD_SIZE == 65521)
    {
      int i = 0;
      const int n8 = batch_getlen (q, n);
      while (i < n)
        {
          uint16_t v0;
          extract (dest8, 2, arg);
          v0 = ((uint16_t *) dest8)[0];
          if (v0 < 65521)
            {
              i++;
              dest8 += 2;
            }
        }
      memset (dest8, 0, CONVX8 (n8) - n8);
    }
#endif
}


#if defined(FIELD_SIZE) && defined(NB_VARIABLES)
#define sklen ((LOG2(FIELD_SIZE)*NB_VARIABLES+7)>>3)
#define sklenX CONVX(sklen)
#endif
#if defined(FIELD_SIZE) && defined(NB_EQUATIONS)
#define pklen ((LOG2(FIELD_SIZE)*NB_EQUATIONS+7)>>3)
#define pklenX CONVX(pklen)
#endif

/* use the extract function to obtain the whole system */
void
circuit_generate (uintX_t *f, int q, int n, int m, int d,
                  void (*extract) (void *, size_t, void *), void *arg)
{
#ifdef DEGREE
  const int d = DEGREE;
#endif
#if defined(FIELD_SIZE) && defined(DEGREE) \
  && defined(NB_VARIABLES) && defined(NB_EQUATIONS)
  const int n = NB_VARIABLES;
  const int m = NB_EQUATIONS;
#endif
  int j, k;
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif

  BATCH_PARAMS (q, n, m, d);

  for (k = 0; k < d + 1; k++)
    {
      batch_generate (f, q, m, extract, arg);
      f += pklenX;
      for (j = 0; j < m; j++)
        {
          batch_generate (f, q, n, extract, arg);
          f += sklenX;
        }
    }
}

/* compute all the multiplications input pairs shares (x, y) and */
/* the m last shares in z according to the input key share sk and ct, */
/* and the previous multiplications' output shares in z */
void
linear_circuit (uintX_t *x, uintX_t *y, uintX_t *z, const uintX_t *sk,
                const uintX_t *ct, int i, int q, int n, int m, int d,
                const uintX_t *f)
{
#ifdef DEGREE
  const int d = DEGREE;
#endif
#if defined(FIELD_SIZE) && defined(DEGREE) \
  && defined(NB_VARIABLES) && defined(NB_EQUATIONS)
#if FIELD_SIZE != 2
  const int n = NB_VARIABLES;
#endif
  const int m = NB_EQUATIONS;
#endif
  int j, k;
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
  uintX_t tmp[sklenX > pklenX ? sklenX : pklenX];
  uintX_t w[pklenX];
  uint8_t *const x_ = (void *) x;
  uint8_t *const y_ = (void *) y;
  uint8_t *const z_ = (void *) z;

  BATCH_PARAMS (q, n, m, d);

  for (k = 0; k < d + 1; k++)
    {
      if (i == 0)
        {
          batch_copy (w, f, q, m);
        }
      else
        {
          batch_clear (w, q, m);
        }
      f += pklenX;
      for (j = 0; j < m; j++)
        {
          batch_copy (tmp, f, q, n);
          batch_mul (tmp, sk, q, n);
          batch_sum (w, tmp, j, q, n);
          f += sklenX;
        }
      if (k == 0)
        {
          batch_export (x_, w, q, m, k * m);
        }
      else if (k < d)
        {
          batch_export (y_, w, q, m, (k - 1) * m);
#if !defined(DEGREE) || (DEGREE > 2)
          if (k < d - 1)
            {
              batch_import (tmp, z_, q, m, (k - 1) * m);
              batch_export (x_, tmp, q, m, k * m);
            }
#endif
        }
    }
  if (i == 0)
    {
      /* subtract from ct to get the last multiplication result */
      batch_copy (tmp, ct, q, m);
    }
  else
    {
      batch_clear (tmp, q, m);
    }
  batch_sub (tmp, w, q, m);
  batch_export (z_, tmp, q, m, (d - 2) * m);
}

/* compute the resulting ciphertext ct as well as intermediate value y and z */
/* for the multiplications according to the input key sk */
void
eval_circuit (uintX_t *y, uintX_t *z, uintX_t *ct,
              const uintX_t *sk, int q, int n, int m, int d,
              const uintX_t *f)
{
#ifdef DEGREE
  const int d = DEGREE;
#endif
#if defined(FIELD_SIZE) && defined(DEGREE) \
  && defined(NB_VARIABLES) && defined(NB_EQUATIONS)
#if FIELD_SIZE != 2
  const int n = NB_VARIABLES;
#endif
  const int m = NB_EQUATIONS;
#endif
  int j, k;
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
  uintX_t x[pklenX];
  uintX_t tmp[sklenX];
  uint8_t *const y_ = (void *) y;
#if !defined(DEGREE) || (DEGREE > 2)
  uint8_t *const z_ = (void *) z;
#endif

  BATCH_PARAMS (q, n, m, d);

  for (k = 0; k < d + 1; k++)
    {
      batch_copy (ct, f, q, m);
      f += pklenX;
      for (j = 0; j < m; j++)
        {
          batch_copy (tmp, f, q, n);
          batch_mul (tmp, sk, q, n);
          batch_sum (ct, tmp, j, q, n);
          f += sklenX;
        }
      if (k == 0)
        {
          batch_copy (x, ct, q, m);
        }
      else if (k < d)
        {
          if (y)
            {
              batch_export (y_, ct, q, m, (k - 1) * m);
            }
          batch_mul (x, ct, q, m);
#if !defined(DEGREE) || (DEGREE > 2)
          if (z && (k < d - 1))
            {
              batch_export (z_, x, q, m, (k - 1) * m);
            }
#endif
        }
    }
  batch_add (ct, x, q, m);
}

/* compute the resulting ciphertext ct as well as intermediate value y and z */
/* for the multiplications according to the input key sk */
/* use the extract function to obtain the whole system on the fly */
void
eval_circuit_seed (uintX_t *y, uintX_t *z, uintX_t *ct,
                   const uintX_t *sk, int q, int n, int m, int d,
                   void (*extract) (void *, size_t, void *), void *arg)
{
#ifdef DEGREE
  const int d = DEGREE;
#endif
#if defined(FIELD_SIZE) && defined(DEGREE) \
  && defined(NB_VARIABLES) && defined(NB_EQUATIONS)
  const int n = NB_VARIABLES;
  const int m = NB_EQUATIONS;
#endif
  int j, k;
#if !defined(FIELD_SIZE) || !defined(NB_VARIABLES)
  const int sklenX = batch_getlenX (q, n);
#endif
#if !defined(FIELD_SIZE) || !defined(NB_EQUATIONS)
  const int pklenX = batch_getlenX (q, m);
#endif
  uintX_t x[pklenX];
  uintX_t tmp[sklenX];
  uint8_t *const y_ = (void *) y;
#if !defined(DEGREE) || (DEGREE > 2)
  uint8_t *const z_ = (void *) z;
#endif

  BATCH_PARAMS (q, n, m, d);

  for (k = 0; k < d + 1; k++)
    {
      batch_generate (ct, q, m, extract, arg);
      for (j = 0; j < m; j++)
        {
          batch_generate (tmp, q, n, extract, arg);
          batch_mul (tmp, sk, q, n);
          batch_sum (ct, tmp, j, q, n);
        }
      if (k == 0)
        {
          batch_copy (x, ct, q, m);
        }
      else if (k < d)
        {
          if (y)
            {
              batch_export (y_, ct, q, m, (k - 1) * m);
            }
          batch_mul (x, ct, q, m);
#if !defined(DEGREE) || (DEGREE > 2)
          if (z && (k < d - 1))
            {
              batch_export (z_, x, q, m, (k - 1) * m);
            }
#endif
        }
    }
  batch_add (ct, x, q, m);
}
