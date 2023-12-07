#include "gf2p32.h"

static int sdith_gf2p16_tables_initialized = 0;
static int sdith_gf2p32_tables_initialized = 0;
#ifdef NO_PRECOMPUTE_TABLES
uint16_t sdith_gf2p16_logtable[65536];
uint16_t sdith_gf2p16_exptable[65536];
uint32_t sdith_gf2p32_log1table[65536];     // for x=u+Y: table of log16(u) -> log32(u+Y)
uint16_t sdith_gf2p32_exp1table[65537][2];  // table of p -> (log16(u),log16(v))  where gen32^p = u+vY
#else
extern uint16_t sdith_gf2p16_logtable[];
extern uint16_t sdith_gf2p16_exptable[];
extern uint32_t sdith_gf2p32_log1table[];
extern uint16_t sdith_gf2p32_exp1table[65537][2];
#endif

/** @brief naive multiplication in gf2p16 */
uint16_t gf2p16_mul_naive(uint16_t x, uint16_t y) {
  uint8_t* xx = (uint8_t*)&x;
  uint8_t* yy = (uint8_t*)&y;
  uint16_t a0 = mul_gf256_naive(xx[0], yy[0]);
  uint16_t a1 = mul_gf256_naive(xx[0], yy[1]) ^ mul_gf256_naive(xx[1], yy[0]);
  uint16_t a2 = mul_gf256_naive(xx[1], yy[1]);
  a1 ^= a2;
  a0 ^= mul_gf256_naive(SDITH_IRRED_CST_GF2P16, a2);
  return (a1 << 8) | a0;
}

/** @brief naive multiplication in gf2p32 */
uint32_t gf2p32_mul_naive(uint32_t x, uint32_t y) {
  uint16_t* xx = (uint16_t*)&x;
  uint16_t* yy = (uint16_t*)&y;
  uint32_t a0 = gf2p16_mul_naive(xx[0], yy[0]);
  uint32_t a1 = gf2p16_mul_naive(xx[0], yy[1]) ^ gf2p16_mul_naive(xx[1], yy[0]);
  uint32_t a2 = gf2p16_mul_naive(xx[1], yy[1]);
  a1 ^= a2;
  a0 ^= gf2p16_mul_naive(SDITH_IRRED_CST_GF2P32, a2);
  return (a1 << 16u) | a0;
}

/** @brief naive (divide and conquer) pow in gf2p32 */
uint32_t gf2p32_pow_naive(uint32_t x, uint64_t n) {
  uint32_t y = x;
  uint32_t z = 1;
  // keep zy^n constant and decrease n->0
  while (n != 0) {
    if (n % 2 == 1) {
      z = gf2p32_mul_naive(z, y);
    }
    y = gf2p32_mul_naive(y, y);
    n /= 2;
  }
  return z;
}

/** @brief lookup-table based discrete log in gf2p32 (non constant time)
 * Please initialize the gf2p32 log tables once before using this function */
uint32_t gf2p32_dlog(const uint32_t x) {
  const uint32_t lu = sdith_gf2p16_logtable[x & 0xFFFF];
  const uint32_t lv = sdith_gf2p16_logtable[x >> 16u];
  if (lv == 0xFFFF) return 0x10001 * lu;  // x is in gf2p16
  const uint32_t ldiff = (lu == 0xFFFFU ? 0xFFFF : (lu >= lv ? lu - lv : 0xFFFF + lu - lv));
  return ((uint64_t)(sdith_gf2p32_log1table[ldiff]) + 0x10001UL * lv) % 0xFFFFFFFFUL;
}

/** @brief lookup-table based discrete exp in gf2p32 (non constant time)
 * Please initialize the gf2p32 log tables once before using this function */
uint32_t gf2p32_dexp(const uint32_t x) {
  if (x == 0xFFFFFFFFu) return 0;
  // decompose x = K.lv + lu
  const uint32_t lv = x / 0x10001u;  // lv in [0,0xFFFE]
  const uint32_t lu = x % 0x10001u;  // lu in [0,65536]
  if (lu == 0) return sdith_gf2p16_exptable[lv];
  uint16_t const* r = sdith_gf2p32_exp1table[lu];
  const uint32_t lru = r[0] == 0xFFFF ? 0xFFFF : (r[0] + lv) % 0xFFFF;
  const uint32_t lrv = r[1] == 0xFFFF ? 0xFFFF : (r[1] + lv) % 0xFFFF;
  return ((uint32_t)(sdith_gf2p16_exptable[lrv]) << 16) | (sdith_gf2p16_exptable[lru]);
}

#include <stdlib.h>
#include <stdio.h>

extern uint16_t gf2p16_logtable_precomputed[];
extern uint16_t gf2p16_exptable_precomputed[];

void create_gf2p16_log_tables() {
  if (sdith_gf2p16_tables_initialized) return;
#ifdef NO_PRECOMPUTE_TABLES
  gf256_create_log_tables();
  sdith_gf2p16_exptable[0xFFFF] = 0;
  sdith_gf2p16_exptable[0] = 1;
  sdith_gf2p16_exptable[1] = SDITH_GEN_GF2P16;
  for (uint64_t i = 2; i < 0x101; ++i) {
    sdith_gf2p16_exptable[i] = gf2p16_mul_naive(sdith_gf2p16_exptable[i - 1], SDITH_GEN_GF2P16);
  }
  for (uint64_t i = 0; i < 0x101; ++i) {
    uint16_t l = sdith_gf2p16_exptable[i];
    uint8_t l0 = sdith_gf256_dlog_table[l & 0xFF];
    uint8_t l1 = sdith_gf256_dlog_table[l >> 8];
    for (uint64_t j = 1; j < 255; ++j) {
      l0 = l0 == 255 ? 255 : l0 == 254 ? 0 : l0 + 1;
      l1 = l1 == 255 ? 255 : l1 == 254 ? 0 : l1 + 1;
      sdith_gf2p16_exptable[i + 0x101 * j] =
          (uint16_t)(sdith_gf256_dexp_table[l0]) | ((uint16_t)(sdith_gf256_dexp_table[l1]) << 8);
    }
  }
  for (uint64_t i = 0; i < 65536; ++i) {
    sdith_gf2p16_logtable[sdith_gf2p16_exptable[i]] = i;
  }
#endif
  sdith_gf2p16_tables_initialized = 1;
}

void create_gf2p32_log_tables() {
  if (sdith_gf2p32_tables_initialized) return;
#ifdef NO_PRECOMPUTE_TABLES
  create_gf2p16_log_tables();
  uint32_t z = 1;
  uint16_t* zz = (uint16_t*)&z;
  // g^0 = 1 + 0y
  sdith_gf2p32_exp1table[0][0] = 0x0000;
  sdith_gf2p32_exp1table[0][1] = 0xFFFF;
  for (uint64_t i = 1; i < 0x10001; ++i) {
    // g^i = u + Y v.
    z = gf2p32_mul_naive(z, SDITH_GEN_GF2P32);
    uint16_t lu = sdith_gf2p16_logtable[zz[0]];
    uint16_t lv = sdith_gf2p16_logtable[zz[1]];
    sdith_gf2p32_exp1table[i][0] = lu;
    sdith_gf2p32_exp1table[i][1] = lv;
    // ASSERT_DRAMATICALLY(lv!=0xFFFFu, "bug!! v cannot be zero");
    uint16_t ldiff = (lu == 0xFFFF) ? 0xFFFF : lu >= lv ? lu - lv : 0xFFFF + lu - lv;
    // g32^(i-K.lv) = (u/v) + y = g16^ldiff + Y
    sdith_gf2p32_log1table[ldiff] = (lv == 0) ? i : 0xFFFFFFFFu + i - 0x10001 * lv;
    // REQUIRE_DRAMATICALLY(gf2p32_pow_naive(SDITH_GEN_GF2P32, sdith_gf2p32_log1table[ldiff]) == 0x10000U +
    // sdith_gf2p16_exptable[ldiff], "XXXX")
  }
#endif
  sdith_gf2p32_tables_initialized = 1;
}

uint32_t gf2p32_dlog_pow_l32(uint32_t logx, uint32_t p) {
  if (logx == 0xffffffff) return 0xffffffff;
  return ((uint64_t)p * (uint64_t)logx) % 0xffffffff;
}

uint32_t gf2p32_dlog_mul_l32_l32(uint32_t logx, uint32_t logy) {
  if (logx == 0xffffffff || logy == 0xffffffff) return 0xffffffff;
  static const uint32_t order = 0xffffffff;
  uint64_t l = (uint64_t)logx + (uint64_t)logy;
  if (l >= order) l -= order;
  return l;
}

uint32_t gf2p32_mul_table(uint32_t x, uint32_t y) {
  return gf2p32_dexp(gf2p32_dlog_mul_l32_l32(gf2p32_dlog(x), gf2p32_dlog(y)));
}

uint32_t (*gf2p32_mul_ct)(uint32_t x, uint32_t y) = 0;
uint32_t (*gf2p32_mul)(uint32_t x, uint32_t y) = 0;

EXPORT void gf2p32_mul_init() {
  create_gf2p32_log_tables();
  gf2p32_mul_ct = gf2p32_mul_naive;
  gf2p32_mul = gf2p32_mul_table;
}
