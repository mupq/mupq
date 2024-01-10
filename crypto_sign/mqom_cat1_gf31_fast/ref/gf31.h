#ifndef MQOM_GF31_H
#define MQOM_GF31_H

#include <stdint.h>

/*************************************/
/***          Compression          ***/
/*************************************/

void gf31_compress_tab(uint8_t* buf, const uint8_t* v, uint32_t size);
void gf31_decompress_tab(uint8_t* v, const uint8_t* buf, uint32_t size);

/*************************************/
/***        Basic Operations       ***/
/*************************************/

// Return "a+b" over GF(31)
uint8_t gf31_add(uint8_t a, uint8_t b);

// Return "a-b" over GF(31)
uint8_t gf31_sub(uint8_t a, uint8_t b);

// Return "a*b" over GF(31)
uint8_t gf31_mul(uint8_t a, uint8_t b);

// Return "-a" over GF(31)
uint8_t gf31_neg(uint8_t a);

// Return x[0]*y[0]+...+x[size-1]*y[size-1] over GF(31)
uint8_t gf31_innerproduct(const uint8_t* x, const uint8_t* y, uint32_t size);

/*************************************/
/***       Batched Operations      ***/
/*************************************/

// vz[] += vx[] over GF(31)
void gf31_add_tab(uint8_t* vz, const uint8_t* vx, uint32_t size);

// vz[] -= vx[] over GF(31)
void gf31_sub_tab(uint8_t* vz, const uint8_t* vx, uint32_t size);

// vz[] = vx[] * y over GF(31)
void gf31_mul_tab(uint8_t* vz, const uint8_t* vx, uint8_t y, uint32_t size);

// vz[] += vx[] * y over GF(31)
void gf31_muladd_tab(uint8_t* vz, uint8_t y, const uint8_t* vx, uint32_t size);

// vz[] = -vz[] over GF(31)
void gf31_neg_tab(uint8_t* vz, uint32_t size);

/*************************************/
/***        Matrix Operations      ***/
/*************************************/

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb]
void gf31_matcols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

// vz[] += vx[1][] * y[1] + ... + vx[nb][] * y[nb]
void gf31_matrows_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb], size must be a multiple of 128
void gf31_mat128cols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

void gf31_matcols_muladd_triangular(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t step);

/*************************************/
/***        Field Extensions       ***/
/*************************************/

/// Extension Degree 6

// Check if a == b over GF(31^6)
int gf31to6_eq(const uint8_t a[6], const uint8_t b[6]);
// Set res as a+b over GF(31^6)
void gf31to6_add(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]);
// Set res as a-b over GF(31^6)
void gf31to6_sub(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]);
// Set res as a*b over GF(31^6)
void gf31to6_mul(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]);
// Set res as a*b over GF(31^6), a belongs to GF(31)
void gf31to6_mul_gf31(uint8_t res[6], uint8_t a, const uint8_t b[6]);

/// Extension Degree 7

// Check if a == b over GF(31^7)
int gf31to7_eq(const uint8_t a[7], const uint8_t b[7]);
// Set res as a+b over GF(31^7)
void gf31to7_add(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a-b over GF(31^7)
void gf31to7_sub(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a*b over GF(31^7)
void gf31to7_mul(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a*b over GF(31^7), a belongs to GF(31)
void gf31to7_mul_gf31(uint8_t res[7], uint8_t a, const uint8_t b[7]);

/// Extension Degree 8

// Check if a == b over GF(31^8)
int gf31to8_eq(const uint8_t a[8], const uint8_t b[8]);
// Set res as a+b over GF(31^8)
void gf31to8_add(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]);
// Set res as a-b over GF(31^8)
void gf31to8_sub(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]);
// Set res as a*b over GF(31^8)
void gf31to8_mul(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]);
// Set res as a*b over GF(31^8), a belongs to GF(31)
void gf31to8_mul_gf31(uint8_t res[8], uint8_t a, const uint8_t b[8]);

/// Extension Degree 10

// Check if a == b over GF(31^10)
int gf31to10_eq(const uint8_t a[10], const uint8_t b[10]);
// Set res as a+b over GF(31^10)
void gf31to10_add(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]);
// Set res as a-b over GF(31^10)
void gf31to10_sub(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]);
// Set res as a*b over GF(31^10)
void gf31to10_mul(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]);
// Set res as a*b over GF(31^10), a belongs to GF(31)
void gf31to10_mul_gf31(uint8_t res[10], uint8_t a, const uint8_t b[10]);

/// Extension Degree 11

// Check if a == b over GF(31^11)
int gf31to11_eq(const uint8_t a[11], const uint8_t b[11]);
// Set res as a+b over GF(31^11)
void gf31to11_add(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]);
// Set res as a-b over GF(31^11)
void gf31to11_sub(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]);
// Set res as a*b over GF(31^11)
void gf31to11_mul(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]);
// Set res as a*b over GF(31^11), a belongs to GF(31)
void gf31to11_mul_gf31(uint8_t res[11], uint8_t a, const uint8_t b[11]);

/*************************************/
/***           Randomness          ***/
/*************************************/

#ifndef NO_RND
#include "sample.h"
// Fill the array "points" with random GF(31) elements
void gf31_random_elements(uint8_t points[], uint32_t nb_points, samplable_t* entropy);
void gf31_random_elements_x4(uint8_t* const* points, uint32_t nb_points, samplable_x4_t* entropy);
#endif

#endif /* MQOM_GF31_H */
