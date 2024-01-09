#ifndef MQOM_GF251_H
#define MQOM_GF251_H

#include <stdint.h>

/*************************************/
/***        Basic Operations       ***/
/*************************************/

// Return "a+b" over GF(251)
uint8_t gf251_add(uint8_t a, uint8_t b);

// Return "a-b" over GF(251)
uint8_t gf251_sub(uint8_t a, uint8_t b);

// Return "a*b" over GF(251)
uint8_t gf251_mul(uint8_t a, uint8_t b);

// Return "-a" over GF(251)
uint8_t gf251_neg(uint8_t a);

// Return x[0]*y[0]+...+x[size-1]*y[size-1] over GF(251)
uint8_t gf251_innerproduct(const uint8_t* x, const uint8_t* y, uint32_t size);

/*************************************/
/***       Batched Operations      ***/
/*************************************/

// vz[] += vx[] over GF(251)
void gf251_add_tab(uint8_t* vz, const uint8_t* vx, uint32_t size);

// vz[] -= vx[] over GF(251)
void gf251_sub_tab(uint8_t* vz, const uint8_t* vx, uint32_t size);

// vz[] = vx[] * y over GF(251)
void gf251_mul_tab(uint8_t* vz, const uint8_t* vx, uint8_t y, uint32_t size);

// vz[] += vx[] * y over GF(251)
void gf251_muladd_tab(uint8_t* vz, uint8_t y, const uint8_t* vx, uint32_t size);

// vz[] = -vz[] over GF(251)
void gf251_neg_tab(uint8_t* vz, uint32_t size);

/*************************************/
/***        Matrix Operations      ***/
/*************************************/

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb]
void gf251_matcols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

// vz[] += vx[1][] * y[1] + ... + vx[nb][] * y[nb]
void gf251_matrows_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb], size must be a multiple of 128
void gf251_mat128cols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size);

void gf251_matcols_muladd_triangular(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t step);

/*************************************/
/***        Field Extensions       ***/
/*************************************/

/// Extension Degree 4

// Check if a == b over GF(251^4)
int gf251to4_eq(const uint8_t a[4], const uint8_t b[4]);
// Set res as a+b over GF(251^4)
void gf251to4_add(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]);
// Set res as a-b over GF(251^4)
void gf251to4_sub(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]);
// Set res as a*b over GF(251^4)
void gf251to4_mul(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]);
// Set res as a*b over GF(251^4), a belongs to GF(251)
void gf251to4_mul_gf251(uint8_t res[4], uint8_t a, const uint8_t b[4]);

/// Extension Degree 5

// Check if a == b over GF(251^5)
int gf251to5_eq(const uint8_t a[5], const uint8_t b[5]);
// Set res as a+b over GF(251^5)
void gf251to5_add(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]);
// Set res as a-b over GF(251^5)
void gf251to5_sub(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]);
// Set res as a*b over GF(251^5)
void gf251to5_mul(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]);
// Set res as a*b over GF(251^5), a belongs to GF(251)
void gf251to5_mul_gf251(uint8_t res[5], uint8_t a, const uint8_t b[5]);

/// Extension Degree 7

// Check if a == b over GF(251^7)
int gf251to7_eq(const uint8_t a[7], const uint8_t b[7]);
// Set res as a+b over GF(251^7)
void gf251to7_add(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a-b over GF(251^7)
void gf251to7_sub(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a*b over GF(251^7)
void gf251to7_mul(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]);
// Set res as a*b over GF(251^7), a belongs to GF(251)
void gf251to7_mul_gf251(uint8_t res[7], uint8_t a, const uint8_t b[7]);

/*************************************/
/***           Randomness          ***/
/*************************************/

#ifndef NO_RND
#include "sample.h"
// Fill the array "points" with random GF(251) elements
void gf251_random_elements(uint8_t points[], uint32_t nb_points, samplable_t* entropy);
void gf251_random_elements_x4(uint8_t* const* points, uint32_t nb_points, samplable_x4_t* entropy);
#endif

#endif /* MQOM_GF251_H */
