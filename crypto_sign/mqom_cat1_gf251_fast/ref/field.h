#ifndef FIELD_H
#define FIELD_H

#include <stdint.h>
#include <string.h>
#include "parameters.h"

#include "gf251.h"

// Unitary operations
#define add_points(a,b)      gf251_add(a,b)
#define sub_points(a,b)      gf251_sub(a,b)
#define mul_points(a,b)      gf251_mul(a,b)
#define neg_point(a)         gf251_neg(a)
#define innerproduct_points(x,y,size) gf251_innerproduct(x,y,size)

// Batched operations
static inline void set_zero_tab_points(void* vz, uint32_t size) { memset(vz, 0, size); }
static inline void set_tab_points(void* vz, const void* vx, uint32_t size) { memcpy(vz, vx, size); }
static inline void add_tab_points(void* vz, const void* vx, uint32_t size) { gf251_add_tab(vz, vx, size); }
static inline void sub_tab_points(void* vz, const void* vx, uint32_t size) { gf251_sub_tab(vz, vx, size); }
static inline void mul_tab_points(void* vz, const void* vx, uint8_t y, uint32_t size) { gf251_mul_tab(vz, vx, y, size); }
static inline void muladd_tab_points(void* vz, uint8_t y, const void* vx, uint32_t size) { gf251_muladd_tab(vz, y, vx, size); }
static inline void neg_tab_points(void* vz, uint32_t size) { gf251_neg_tab(vz, size); }
static inline void matcols_muladd(void* vz, const void* y, const void* vx, uint32_t nb_columns, uint32_t nb_rows) {
    gf251_matcols_muladd(vz, y, vx, nb_columns, nb_rows);
}
static inline void matrows_muladd(void* vz, const void* y, const void* vx, uint32_t nb_columns, uint32_t nb_rows) {
    gf251_matrows_muladd(vz, y, vx, nb_columns, nb_rows);
}
static inline void mat128cols_muladd(void* vz, const void* y, const void* vx, uint32_t nb_columns, uint32_t nb_rows) {
    gf251_mat128cols_muladd(vz, y, vx, nb_columns, nb_rows);
}
static inline void matcols_muladd_triangular(void* vz, const void* y, const void* vx, uint32_t nb_columns, uint32_t step) {
    gf251_matcols_muladd_triangular(vz, y, vx, nb_columns, step);
}

// Operations in field extension
#define JOIN3_0(a,b,c)              a ## b ## c
#define JOIN3(a,b,c)                JOIN3_0(a, b, c)

#define mul_points_ext(r,a,b)       JOIN3(gf251to,PARAM_EXT_DEGREE,_mul)(r,a,b)
#define add_points_ext(r,a,b)       JOIN3(gf251to,PARAM_EXT_DEGREE,_add)(r,a,b)
#define sub_points_ext(r,a,b)       JOIN3(gf251to,PARAM_EXT_DEGREE,_sub)(r,a,b)
#define mul_points_mixed(r,a,b)     JOIN3(gf251to,PARAM_EXT_DEGREE,_mul_gf251)(r,a,b)
#define eq_points_ext(a,b)          JOIN3(gf251to,PARAM_EXT_DEGREE,_eq)(a,b)

// Randomness
static inline void random_points(void* points, uint32_t size, samplable_t* entropy) { gf251_random_elements(points, size, entropy); }
static inline void random_points_x4(void* const* points, uint32_t size, samplable_x4_t* entropy) { gf251_random_elements_x4((uint8_t *const *) points, size, entropy); }

#endif /* FIELD_H */
