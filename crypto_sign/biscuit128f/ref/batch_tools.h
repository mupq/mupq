#ifndef _BATCH_TOOLS_H_
#define _BATCH_TOOLS_H_

#ifndef UINTX_BITSIZE
#define UINTX_BITSIZE 64
#endif

#if UINTX_BITSIZE > 64
#include <immintrin.h>
#endif
#include <stdint.h>

#include "params_posso.h"

#define NBITS1(x) ((x) >= (1 << 0) ? 1 : 0)
#define NBITS2(x) ((x) >= (1 << 1) ? 1 + NBITS1 ((x) >> 1) : NBITS1 (x))
#define NBITS4(x) ((x) >= (1 << 2) ? 2 + NBITS2 ((x) >> 2) : NBITS2 (x))
#define NBITS8(x) ((x) >= (1 << 4) ? 4 + NBITS4 ((x) >> 4) : NBITS4 (x))
#define NBITS16(x) ((x) >= (1 << 8) ? 8 + NBITS8 ((x) >> 8) : NBITS8 (x))
#define NBITS32(x) ((x) >= (1 << 16) ? 16 + NBITS16 ((x) >> 16) : NBITS16 (x))
#define LOG2(x) NBITS32 ((x) - 1)
#define IS_POW2(x) (!((x) & ((x) - 1)))

#if UINTX_BITSIZE == 256
#define _CONVX(x) (((x) + 31) >> 5)
#define _CONVX8(x) (((x) + 31) & ~0x1f)
typedef __m256i uintX_t;
union vecX_t {
  uintX_t tX;
  __m256i t256[1];
  __m128i t128[2];
  uint64_t t64[4];
  uint32_t t32[8];
  uint16_t t16[16];
  uint8_t t8[32];
};
#elif UINTX_BITSIZE == 128
#define _CONVX(x) (((x) + 15) >> 4)
#define _CONVX8(x) (((x) + 15) & ~0xf)
typedef __m128i uintX_t;
union vecX_t {
  uintX_t tX;
  __m128i t128[1];
  uint64_t t64[2];
  uint32_t t32[4];
  uint16_t t16[8];
  uint8_t t8[16];
};
#elif UINTX_BITSIZE == 64
#define _CONVX(x) (((x) + 7) >> 3)
#define _CONVX8(x) (((x) + 7) & ~0x7)
typedef uint64_t uintX_t;
union vecX_t {
  uintX_t tX;
  uint64_t t64[1];
  uint32_t t32[2];
  uint16_t t16[4];
  uint8_t t8[8];
};
#elif UINTX_BITSIZE == 32
#define _CONVX(x) (((x) + 3) >> 2)
#define _CONVX8(x) (((x) + 3) & ~0x3)
typedef uint32_t uintX_t;
union vecX_t {
  uintX_t tX;
  uint32_t t32[1];
  uint16_t t16[2];
  uint8_t t8[4];
};
#elif UINTX_BITSIZE == 16
#define _CONVX(x) (((x) + 1) >> 1)
#define _CONVX8(x) (((x) + 1) & ~0x1)
typedef uint16_t uintX_t;
union vecX_t {
  uintX_t tX;
  uint16_t t16[1];
  uint8_t t8[2];
};
#elif UINTX_BITSIZE == 8
#define _CONVX(x) (x)
#define _CONVX8(x) (x)
typedef uint8_t uintX_t;
union vecX_t {
  uintX_t tX;
  uint8_t t8[1];
};
#endif

#define BATCH_PARAMS(q, n, m, d)

#ifdef FIELD_SIZE
#if defined(DEGREE) && defined(NB_VARIABLES) && defined(NB_EQUATIONS)

#define NB_MAX ((NB_VARIABLES) > ((DEGREE - 1) * NB_EQUATIONS) ?   \
                (NB_VARIABLES) : ((DEGREE - 1) * NB_EQUATIONS))
#define CONVX(x) _CONVX(((LOG2 (FIELD_SIZE) * NB_MAX) + 7) >> 3)
#define CONVX8(x) _CONVX8(((LOG2 (FIELD_SIZE) * NB_MAX) + 7) >> 3)

#define batch_getbitlen(q, n) batch_getbitlen_impl (n)
#define batch_getlen(q, n) batch_getlen_impl (n)
#define batch_getlenX(q, n) batch_getlenX_impl ()
#define batch_import(dest, src, q, n, k) batch_import_impl (dest, src, n, k)
#define batch_export(dest, src, q, n, k) batch_export_impl (dest, src, n, k)
#define batch_clear(dest, q, n) batch_clear_impl (dest)
#define batch_copy(dest, src, q, n) batch_copy_impl (dest, src)
#if IS_POW2(FIELD_SIZE)
#define batch_add(dest, src, q, n) batch_add_impl (dest, src)
#define batch_sub(dest, src, q, n) batch_sub_impl (dest, src)
#define batch_sum(dest, buf, j, q, n) batch_sum_impl (dest, buf, j)
#else
#define batch_add(dest, src, q, n) batch_add_impl (dest, src, n)
#define batch_sub(dest, src, q, n) batch_sub_impl (dest, src, n)
#define batch_sum(dest, buf, j, q, n) batch_sum_impl (dest, buf, j, n)
#endif
#if FIELD_SIZE == 2
#define batch_mul(dest, src, q, n) batch_mul_impl (dest, src)
#else
#define batch_mul(dest, src, q, n) batch_mul_impl (dest, src, n)
#endif
#define batch_generate(output, q, n, extract, arg) \
  batch_generate_impl (output, n, extract, arg)

#define circuit_generate(f, q, n, m, d, extract, arg)             \
  circuit_generate_impl (f, extract, arg)
#define linear_circuit(x, y, z, sk, ct, i, q, n, m, d, f)         \
  linear_circuit_impl (x, y, z, sk, ct, i, f)
#if DEGREE == 2
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, ct, sk, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, ct, sk, extract, arg)
#else
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, z, ct, sk, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, z, ct, sk, extract, arg)
#endif

#else

#define batch_getbitlen(q, n) batch_getbitlen_impl (n)
#define batch_getlen(q, n) batch_getlen_impl (n)
#define batch_getlenX(q, n) batch_getlenX_impl (n)
#define batch_import(dest, src, q, n, k) batch_import_impl (dest, src, n, k)
#define batch_export(dest, src, q, n, k) batch_export_impl (dest, src, n, k)
#define batch_clear(dest, q, n) batch_clear_impl (dest, n)
#define batch_copy(dest, src, q, n) batch_copy_impl (dest, src, n)
#define batch_add(dest, src, q, n) batch_add_impl (dest, src, n)
#define batch_sub(dest, src, q, n) batch_sub_impl (dest, src, n)
#define batch_mul(dest, src, q, n) batch_mul_impl (dest, src, n)
#define batch_sum(dest, buf, j, q, n) batch_sum_impl (dest, buf, j, n)
#define batch_generate(output, q, n, extract, arg) \
  batch_generate_impl (output, n, extract, arg)

#ifdef DEGREE

#define circuit_generate(f, q, n, m, d, extract, arg)             \
  circuit_generate_impl (f, n, m, extract, arg)
#define linear_circuit(x, y, z, sk, ct, i, q, n, m, d, f)         \
  linear_circuit_impl (x, y, z, sk, ct, i, n, m, f)
#if DEGREE == 2
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, ct, sk, n, m, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, ct, sk, n, m, extract, arg)
#else
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, z, ct, sk, n, m, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, z, ct, sk, n, m, extract, arg)
#endif

#else

#define circuit_generate(f, q, n, m, d, extract, arg)             \
  circuit_generate_impl (f, n, m, d, extract, arg)
#define linear_circuit(x, y, z, sk, ct, i, q, n, m, d, f)         \
  linear_circuit_impl (x, y, z, sk, ct, i, n, m, d, f)
#if DEGREE == 2
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, ct, sk, n, m, d, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, ct, sk, n, m, d, extract, arg)
#else
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, z, ct, sk, n, m, d, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, z, ct, sk, n, m, d, extract, arg)
#endif
#endif

#endif

#else

#ifdef DEGREE

#define circuit_generate(f, q, n, m, d, extract, arg)             \
  circuit_generate_impl (f, q, n, m, extract, arg)
#define linear_circuit(x, y, z, sk, ct, i, q, n, m, d, f)         \
  linear_circuit_impl (x, y, z, sk, ct, i, q, n, m, f)
#if DEGREE == 2
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, ct, sk, q, n, m, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, ct, sk, q, n, m, extract, arg)
#else
#define eval_circuit(y, z, ct, sk, q, n, m, d, f)                 \
  eval_circuit_impl (y, z, ct, sk, q, n, m, f)
#define eval_circuit_seed(y, z, ct, sk, q, n, m, d, extract, arg) \
  eval_circuit_seed_impl (y, z, ct, sk, q, n, m, extract, arg)
#endif

#endif

#endif

#ifndef CONVX
#define CONVX(x) _CONVX(x)
#endif
#ifndef CONVX8
#define CONVX8(x) _CONVX8(x)
#endif

int
batch_getbitlen (int q, int n);

int
batch_getlen (int q, int n);

int
batch_getlenX (int q, int n);

void
batch_import (uintX_t *dest, const uint8_t *src, int q, int n, int k);

void
batch_export (uint8_t *dest, const uintX_t *src, int q, int n, int k);

void
batch_clear (uintX_t *dest, int q, int n);

void
batch_copy (uintX_t *dest, const uintX_t *src, int q, int n);

void
batch_add (uintX_t *dest, const uintX_t *src, int q, int n);

void
batch_sub (uintX_t *dest, const uintX_t *src, int q, int n);

void
batch_mul (uintX_t *dest, const uintX_t *src, int q, int n);

void
batch_sum (uintX_t *dest, uintX_t *buf, int j, int q, int n);

void
batch_generate (uintX_t *dest, int q, int n,
                void (*extract) (void *, size_t, void *), void *arg);

void
circuit_generate (uintX_t *f, int q, int n, int m, int d,
                  void (*extract) (void *, size_t, void *), void *arg);

void
linear_circuit (uintX_t *x, uintX_t *y, uintX_t *z, const uintX_t *sk,
                const uintX_t *ct, int i, int q, int n, int m, int d,
                const uintX_t *f);

void
eval_circuit (uintX_t *y, uintX_t *z, uintX_t *ct,
              const uintX_t *sk, int q, int n, int m, int d,
              const uintX_t *f);

void
eval_circuit_seed (uintX_t *y, uintX_t *z, uintX_t *ct,
                   const uintX_t *sk, int q, int n, int m, int d,
                   void (*extract) (void *, size_t, void *), void *arg);


#endif
