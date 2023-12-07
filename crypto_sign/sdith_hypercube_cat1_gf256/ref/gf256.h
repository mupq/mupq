#ifndef GF256_SDITH_ARITHMETIC_GF256_H_
#define GF256_SDITH_ARITHMETIC_GF256_H_

#ifdef __cplusplus
#include <cstdint>
#define EXPORT extern "C"
#define EXPORT_DECL extern "C"
#else
#define EXPORT
#define EXPORT_DECL extern
#include <stdint.h>
#endif

EXPORT void gf256_init(__attribute__((unused)) __attribute__((unused)) uint8_t verbose);

/** @brief naive multiplication in gf256, the generator polynomial is 0x11b = x^8+x^4+x^3+x+1 */
EXPORT uint8_t mul_gf256_naive(uint8_t x, uint8_t y);
/** @brief naive multiplication in gf256^3, the generator polynomial is x^3+2 */
EXPORT uint32_t mul_gf2p24_naive(uint32_t x, uint32_t y);
/** @brief multiplication in GF2p8. (lookup table based) */
EXPORT uint8_t mul_gf256_table(uint8_t x, uint8_t y);

/**
 *  @brief Performs vz[16] += vx[m] * my[m][16] using the best available implementation (may be non constant time)
 *  @param vz accumulated result (read before write)
 *  @param vx  accumulated result (read only)
 *  @param my  matrix (read only, must be 32 bytes aligned)
*  */
EXPORT_DECL void (*gf256_vec_mat16cols_muladd)(void* vz, void const* vx, void const* my, uint64_t m);
/**
 *  @brief Performs vz[16] += vx[m] * my[m][16] using the best available constant time implementation
 *  @param vz accumulated result (read before write)
 *  @param vx  accumulated result (read only)
 *  @param my  matrix (read only, must be 32 bytes aligned)
*  */
EXPORT_DECL void (*gf256_vec_mat16cols_muladd_ct)(void* vz, void const* vx, void const* my, uint64_t m);

/** @brief pure C constant time implementation */
EXPORT void gf256_vec_mat16cols_muladd_ref_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief avx2 non-constant time (sbox) implementation) */
EXPORT void gf256_vec_mat16cols_muladd_avx2(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief avx2 constant time (polytable) implementation) */
EXPORT void gf256_vec_mat16cols_muladd_polytable_avx2_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief gfni constant time implementation) */
EXPORT void gf256_vec_mat16cols_muladd_gfni_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief pclmul+avx2 constant time implementation) */
EXPORT void gf256_vec_mat16cols_muladd_pclmul_ct(void* vz, void const* vx, void const* my, uint64_t m);

/**
 *  @brief Performs vz[128] += vx[m] * my[m][128] using the best available implementation (may be non constant time)
 *  @param vz accumulated result (read before write)
 *  @param vx  accumulated result (read only)
 *  @param my  matrix (read only, must be 32 bytes aligned)
*  */
EXPORT_DECL void (*gf256_vec_mat128cols_muladd)(void* vz, void const* vx, void const* my, uint64_t m);
/**
 *  @brief Performs vz[128] += vx[m] * my[m][128] using the best available implementation (may be non constant time)
 *  @param vz accumulated result (read before write)
 *  @param vx  accumulated result (read only)
 *  @param my  matrix (read only, must be 32 bytes aligned)
*  */
EXPORT_DECL void (*gf256_vec_mat128cols_muladd_ct)(void* vz, void const* vx, void const* my, uint64_t m);

/** @brief pure C constant time implementation */
EXPORT void gf256_vec_mat128cols_muladd_ref_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief avx2 non-constant time (sbox) implementation) */
EXPORT void gf256_vec_mat128cols_muladd_avx2(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief avx2 non-constant time (polytable) implementation) */
EXPORT void gf256_vec_mat128cols_muladd_polytable_avx2_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief gfni constant time implementation) */
EXPORT void gf256_vec_mat128cols_muladd_gfni_ct(void* vz, void const* vx, void const* my, uint64_t m);
/** @brief pclmul+avx2 constant time implementation) */
EXPORT void gf256_vec_mat128cols_muladd_pclmul_ct(void* vz, void const* vx, void const* my, uint64_t m);

#define SDITH_GEN_GF256 0x41 // multiplicative generator of gf256

EXPORT_DECL uint8_t const* sdith_gf256_dexp_table;
EXPORT_DECL uint8_t const* sdith_gf256_dlog_table;

/** @brief extended discrete log in GF2p8. log(0)=255 by convention, otherwise,
 * log in base u */
EXPORT uint8_t gf256_dlog(uint8_t x);
/** @brief extended power in GF2p8. */
EXPORT uint8_t gf256_pow(uint8_t x, uint8_t p);

EXPORT void gf256_create_log_tables();


#endif // GF256_SDITH_ARITHMETIC_GF256_H_
