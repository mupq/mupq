#ifndef SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P32_H_
#define SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P32_H_

#include "gf256.h"

#ifdef __cplusplus
#include <cstdint>
#define EXPORT extern "C"
#define EXPORT_DECL extern "C"
#else
#include <stdint.h>
#define EXPORT
#define EXPORT_DECL extern
#endif

// We are going to represent GF2p16 as P(X) mod X^2-0x20 mod over GF256.
// A field element a+bX is represented in binary by the concatenation of two bytes (a,b)
// We are going to represent GF2p32 as Q(Y) mod Y^2-0x2000 over GF2p16
// A field element c+dY is represented in binary by the concatenation of two 16-bytes (c,d)

#define SDITH_GEN_GF2P32 0x8d1f8c20
#define SDITH_GEN_GF2P16 0x118a
#define SDITH_IRRED_CST_GF2P16 0x20    // gf2p16 = gf2p8[X]/X^2+X+0x20
#define SDITH_IRRED_CST_GF2P32 0x2000  // gf2p32 = gf2p16[Y]/Y^2+Y+0x2000

/** @brief multiplication in gf2p32 constant time */
EXPORT_DECL uint32_t (*gf2p32_mul_ct)(uint32_t x, uint32_t y);

/** @brief multiplication in gf2p32 non constant time */
EXPORT_DECL uint32_t (*gf2p32_mul)(uint32_t x, uint32_t y);

/** @brief initialize material for multiplication in gf2p32 */
EXPORT void gf2p32_mul_init();

// naive constant time implementation

/** @brief naive multiplication in gf2p16 (constant time) */
EXPORT uint16_t gf2p16_mul_naive(uint16_t x, uint16_t y);

/** @brief naive multiplication in gf2p32 (constant time) */
EXPORT uint32_t gf2p32_mul_naive(uint32_t x, uint32_t y);

/** @brief naive (divide and conquer) pow in gf2p32 (constant time in x, not in n) */
EXPORT uint32_t gf2p32_pow_naive(uint32_t x, uint64_t n);

EXPORT_DECL uint16_t sdith_gf2p16_logtable[65536];
EXPORT_DECL uint16_t sdith_gf2p16_exptable[65536];
EXPORT_DECL uint32_t sdith_gf2p32_log1table[65536];     // for x=u+Y: table of log16(u) -> log32(u+Y)
EXPORT_DECL uint16_t sdith_gf2p32_exp1table[65537][2];  // table of p -> (log16(u),log16(v))  where gen32^p = u+vY

/** @brief lookup-table based discrete log in gf2p32 (non constant time)
 * Please initialize the gf2p32 log tables once before using this function */
EXPORT uint32_t gf2p32_dlog(uint32_t x);

/** @brief lookup-table based discrete exp in gf2p32 (non constant time)
 * Please initialize the gf2p32 log tables once before using this function */
EXPORT uint32_t gf2p32_dexp(uint32_t x);

EXPORT uint32_t gf2p32_dlog_pow_l32(uint32_t logx, uint32_t p);
EXPORT uint32_t gf2p32_dlog_mul_l32_l32(uint32_t logx, uint32_t logy);

/** @brief initialize the gf2p32 log tables, needed for gf2p32 dlog and dexp. */
EXPORT void create_gf2p32_log_tables();

#endif  // SIMPLE_EXAMPLE_SDITH_ARITHMETIC_GF2P32_H_
