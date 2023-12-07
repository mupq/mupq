#include "gf256.h"

#include <stdlib.h>
#include "assertions.h"
#define nullptr 0x0

void (*gf256_vec_mat16cols_muladd)(void* vz, void const* vx, void const* my, uint64_t m) = nullptr;
void (*gf256_vec_mat16cols_muladd_ct)(void* vz, void const* vx, void const* my, uint64_t m) = nullptr;
void (*gf256_vec_mat128cols_muladd)(void* vz, void const* vx, void const* my, uint64_t m) = nullptr;
void (*gf256_vec_mat128cols_muladd_ct)(void* vz, void const* vx, void const* my, uint64_t m) = nullptr;

EXPORT void gf256_init_avx2_polytable_ct();
EXPORT int gf256_init_avx2();
EXPORT void gf256_init_gfni();
EXPORT void gf256_init_pclmul_ct();
EXPORT void gf256_init(__attribute__((unused)) uint8_t verbose) {
#ifdef AVX2
  gf256_init_avx2_polytable_ct();
  gf256_init_avx2();
  gf256_init_pclmul_ct();
#ifdef __GFNI__
  gf256_init_gfni();
  if (__builtin_cpu_supports("gfni")) {
    gf256_vec_mat16cols_muladd = gf256_vec_mat16cols_muladd_gfni_ct;
    gf256_vec_mat16cols_muladd_ct = gf256_vec_mat16cols_muladd_gfni_ct;
    gf256_vec_mat128cols_muladd = gf256_vec_mat128cols_muladd_gfni_ct;
    gf256_vec_mat128cols_muladd_ct = gf256_vec_mat128cols_muladd_gfni_ct;
    return;
  }
#endif // __GFNI
  // TODO: use the processor caps to assign the best candidate,
  // TODO: currently assuming avx2 is available
  gf256_vec_mat16cols_muladd = gf256_vec_mat16cols_muladd_avx2;
  gf256_vec_mat16cols_muladd_ct = gf256_vec_mat16cols_muladd_pclmul_ct;
  gf256_vec_mat128cols_muladd = gf256_vec_mat128cols_muladd_avx2;
  gf256_vec_mat128cols_muladd_ct = gf256_vec_mat128cols_muladd_polytable_avx2_ct;
#else
  gf256_vec_mat16cols_muladd = gf256_vec_mat16cols_muladd_ref_ct;
  gf256_vec_mat16cols_muladd_ct = gf256_vec_mat16cols_muladd_ref_ct;
  gf256_vec_mat128cols_muladd = gf256_vec_mat128cols_muladd_ref_ct;
  gf256_vec_mat128cols_muladd_ct = gf256_vec_mat128cols_muladd_ref_ct;
#endif // AVX2
}

/** @brief naive gf256 multiplication */
uint8_t mul_gf256_naive(uint8_t x, uint8_t y) {
  static const uint16_t B = 0x11b;
  uint16_t xx = x;
  uint16_t r = 0;
  for (uint8_t i=0; i<8; ++i)
    r ^= (-(uint16_t)((y>>i)&1)) & (xx << i);
  for (uint8_t i=15; i>=8; --i)
    r ^= (-(uint16_t)((r>>i)&1)) & (B << (i-8));
  return r;
}

/** @brief naive multiplication in gf2^24 */
uint32_t mul_gf2p24_naive(uint32_t x, uint32_t y) {
  uint32_t res = 0;
  const uint8_t* const xx = (const uint8_t*) (void*) &x;
  const uint8_t* const yy = (const uint8_t*) (void*) &y;
  uint8_t* const zz = (uint8_t*) (void*) &res;
  zz[0] = mul_gf256_naive(xx[0],yy[0]) ^ mul_gf256_naive(2,(
                                                     mul_gf256_naive(xx[1],yy[2]) ^
                                                     mul_gf256_naive(xx[2],yy[1])));
  zz[1] = mul_gf256_naive(xx[0],yy[1]) ^ mul_gf256_naive(xx[1],yy[0]) ^
          mul_gf256_naive(2,mul_gf256_naive(xx[2],yy[2]));
  zz[2] = mul_gf256_naive(xx[0],yy[2]) ^ mul_gf256_naive(xx[1],yy[1]) ^ mul_gf256_naive(xx[2],yy[0]);
  return res;
}

/** @brief pure C constant time implementation */
EXPORT void gf256_vec_mat16cols_muladd_ref_ct(void* vz, void const* vx, void const* my, uint64_t m) {
  static const uint64_t N = 16;
  const uint8_t (*const y)[N] = (uint8_t (*)[N]) my;
  const uint8_t *const x = (uint8_t *) vx;
  uint8_t *const z = (uint8_t*) vz;
  for (uint64_t i=0; i<m; ++i) {
    for (uint64_t j = 0; j < N; ++j) {
      z[j] ^= mul_gf256_naive(x[i], y[i][j]);
    }
  }
}

/** @brief avx2 non-constant time (sbox) implementation) */
EXPORT void gf256_vec_mat128cols_muladd_ref_ct(void* vz, void const* vx, void const* my, uint64_t m) {
  static const uint64_t N = 128;
  const uint8_t (*const y)[N] = (uint8_t (*)[N]) my;
  const uint8_t *const x = (uint8_t *) vx;
  uint8_t *const z = (uint8_t*) vz;
  for (uint64_t i=0; i<m; ++i) {
    for (uint64_t j = 0; j < N; ++j) {
      z[j] ^= mul_gf256_naive(x[i], y[i][j]);
    }
  }
}

uint8_t const* sdith_gf256_dexp_table = nullptr;
uint8_t const* sdith_gf256_dlog_table = nullptr;

uint8_t gf256_dlog(uint8_t x) {
  ASSERT_DRAMATICALLY(sdith_gf256_dlog_table, "dlog tables not loaded!");
  return sdith_gf256_dlog_table[x];
}

uint8_t gf256_dexp(uint8_t x) {
  ASSERT_DRAMATICALLY(sdith_gf256_dexp_table, "dexp tables not loaded!");
  return sdith_gf256_dexp_table[x];
}

uint8_t gf256_log_pow_log(uint8_t logx, uint8_t p) {
  if (logx == 0xff) return 0xff;
  return (p * logx) % 0xff;
}

uint8_t gf256_pow_log(uint8_t logx, uint8_t p) {
  return sdith_gf256_dexp_table[gf256_log_pow_log(logx, p)];
}

uint8_t gf256_pow(uint8_t x, uint8_t p) {
  if (p == 0) return 1;
  return gf256_pow_log(gf256_dlog(x), p);
}

uint8_t mul_gf256_table(uint8_t x, uint8_t y) {
  ASSERT_DRAMATICALLY(sdith_gf256_dlog_table, "dlog tables not loaded!");
  if (x == 0 || y == 0) return 0;
  uint16_t l = (uint16_t)sdith_gf256_dlog_table[x] + (uint16_t)sdith_gf256_dlog_table[y];
  if (l >= 255) l -= 255;
  return sdith_gf256_dexp_table[l];
}

extern uint8_t dexp8_table_precomputed[];
extern uint8_t dlog8_table_precomputed[];

void gf256_create_log_tables() {
  // don't re-create the tables if they are already initialized
  if (sdith_gf256_dlog_table) return;
#ifdef NO_PRECOMPUTE_TABLES
  uint8_t* dlog8_table = malloc(sizeof(uint8_t)*(1ul << 8));
  uint8_t* dexp8_table = malloc(sizeof(uint8_t)*(1ul << 8));
  // create the dlog table over gf2p8
  dexp8_table[0] = 1;
  for (uint64_t i = 1; i < 255; ++i) {
    dexp8_table[i] = mul_gf256_naive(dexp8_table[i-1], SDITH_GEN_GF256);
  }
  dexp8_table[255] = 0;
  for (uint64_t i = 0; i < 256; ++i) {
    dlog8_table[dexp8_table[i]] = i;
  }
#ifndef NDEBUG
  for (uint64_t i = 0; i < (1ul << 8); ++i) {
    REQUIRE_DRAMATICALLY(dlog8_table[dexp8_table[i]] == i, "bug: %ld %ld %ld\n", i, (uint64_t)(dexp8_table[i]),
                         (uint64_t)(dlog8_table[dexp8_table[i]]));
  }
#endif
  sdith_gf256_dexp_table = dexp8_table;
  sdith_gf256_dlog_table = dlog8_table;
#else
  sdith_gf256_dexp_table = dexp8_table_precomputed;
  sdith_gf256_dlog_table = dlog8_table_precomputed;
#endif
}
