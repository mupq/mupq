#include "gf251.h"
#include <stdlib.h>
#include <string.h>

#include "gf251-internal.h"

static int memeq_ct(const uint8_t* a, const uint8_t* b, uint16_t n) {
    int res = 1;
    for(uint16_t i=0; i<n; i++)
        res &= (a[i] == b[i]);
    return res;
}

/*************************************************/
/***********      FIELD OPERATIONS     ***********/
/*************************************************/

uint8_t gf251_add(uint8_t a, uint8_t b) {
    uint16_t res = a+b;
    res -= PRIME*(PRIME <= res);
    return (uint8_t) res;
}

uint8_t gf251_sub(uint8_t a, uint8_t b) {
    uint16_t res = a+PRIME-b;
    res -= PRIME*(PRIME <= res);
    return (uint8_t) res;
}

uint8_t gf251_mul(uint8_t a, uint8_t b) {
    uint16_t res = a*b;
    return (uint8_t)_gf251_reduce16(res);
}

uint8_t gf251_neg(uint8_t a) {
    return (a != 0)*(PRIME - a);
}

uint8_t gf251_innerproduct(const uint8_t* x, const uint8_t* y, uint32_t size) {
    uint32_t z = 0;
    for(uint32_t i=0; i<size; i++)
        z += x[i]*y[i];
    return (uint8_t)_gf251_reduce32(z);
}

/*************************************************/
/***********     BATCHED OPERATIONS    ***********/
/*************************************************/

// vz[] += vx[]
void gf251_add_tab(uint8_t* vz, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i]+vx[i];
        res -= PRIME*(PRIME <= res);
        vz[i] = (uint8_t) res;
    }
}

// vz[] -= vx[]
void gf251_sub_tab(uint8_t* vz, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i]+PRIME-vx[i];
        res -= PRIME*(PRIME <= res);
        vz[i] = (uint8_t) res;
    }
}

// vz[] = vx[] * y
void gf251_mul_tab(uint8_t* vz, const uint8_t* vx, uint8_t y, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vx[i]*y;
        vz[i] = (uint8_t)_gf251_reduce16(res);
    }
}

// vz[] += vx[] * y
void gf251_muladd_tab(uint8_t* vz, uint8_t y, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i] + vx[i]*y;
        vz[i] = (uint8_t)_gf251_reduce16(res);
    }
}

// vz[] = -vz[]
void gf251_neg_tab(uint8_t* vz, uint32_t size) {
    for(uint32_t i=0; i<size; i++)
        vz[i] = (vz[i] != 0)*(PRIME - vz[i]);
}

/*************************************************/
/***********      FIELD EXTENSION      ***********/
/*************************************************/

#define SET_RES_REDUCE8(i,exp) {v=exp; res[i]=(uint8_t)(v-PRIME*(PRIME<=v));}
#define SET_RES_REDUCE16(i,exp) {v=exp; res[i]=(uint8_t)_gf251_reduce16(v);}

// X^2 - 2 = 0
static void gf251to2_mul(uint8_t res[2], const uint8_t a[2], const uint8_t b[2]) {
    uint16_t leading = a[1]*b[1];
    uint16_t cnst = a[0]*b[0];
    uint32_t v;
    v = (uint32_t) (a[0]+a[1])*(b[0]+b[1])+2*PRIME*PRIME-leading-cnst;
    res[1] = (uint8_t)_gf251_reduce32(v);
    v = (uint32_t) 2*leading+cnst;
    res[0] = (uint8_t)_gf251_reduce32(v);
}

static void gf251to2_add(uint8_t res[2], const uint8_t a[2], const uint8_t b[2]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
}

static void gf251to2_sub(uint8_t res[2], const uint8_t a[2], const uint8_t b[2]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
}

static void gf251to2_mul_0x0101(uint8_t res[2], const uint8_t a[2]) {
    uint16_t v;
    SET_RES_REDUCE16(0, 2*a[1]+a[0]);
    SET_RES_REDUCE16(1, a[0]+a[1]);
}

// X^2 - (256*1+1) = 0
void gf251to4_mul(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]) {
    uint8_t leading[2], cnst[2], sum_a[2], sum_b[2];
    gf251to2_mul(leading, &a[2], &b[2]);
    gf251to2_mul(cnst, &a[0], &b[0]);
    gf251to2_add(sum_a, &a[0], &a[2]);
    gf251to2_add(sum_b, &b[0], &b[2]);
    // Compute &res[0]
    gf251to2_mul_0x0101(&res[0],leading);
    gf251to2_add(&res[0], &res[0], cnst);
    // Compute &res[2]
    gf251to2_mul(&res[2], sum_a, sum_b);
    gf251to2_sub(&res[2], &res[2], cnst);
    gf251to2_sub(&res[2], &res[2], leading);
}

void gf251to4_add(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
}

void gf251to4_sub(uint8_t res[4], const uint8_t a[4], const uint8_t b[4]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
}

void gf251to4_mul_gf251(uint8_t res[4], uint8_t a, const uint8_t b[4]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
}

int gf251to4_eq(const uint8_t a[4], const uint8_t b[4]) {
    return memeq_ct(a, b, 4);
}

// X^5 - 3 = 0
void gf251to5_mul(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]) {
    uint32_t inter[9] = {0};
    uint32_t i, j;
    for(i=0; i<5; i++)
        for(j=0; j<5; j++)
            inter[i+j] += a[i]*b[j];
    for(i=5; i<9; i++)
        inter[i-5] += 3*inter[i];
    res[0] = (uint8_t)_gf251_reduce32(inter[0]);
    res[1] = (uint8_t)_gf251_reduce32(inter[1]);
    res[2] = (uint8_t)_gf251_reduce32(inter[2]);
    res[3] = (uint8_t)_gf251_reduce32(inter[3]);
    res[4] = (uint8_t)_gf251_reduce32(inter[4]);
}

void gf251to5_add(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
}

void gf251to5_sub(uint8_t res[5], const uint8_t a[5], const uint8_t b[5]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
}

void gf251to5_mul_gf251(uint8_t res[5], uint8_t a, const uint8_t b[5]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
}

int gf251to5_eq(const uint8_t a[5], const uint8_t b[5]) {
    return memeq_ct(a, b, 5);
}

// X^7 - X - 1 = 0
void gf251to7_mul(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint32_t inter[13] = {0};
    uint32_t i, j;
    for(i=0; i<7; i++)
        for(j=0; j<7; j++)
            inter[i+j] += a[i]*b[j];
    for(i=12; i>=7; i--) {
        inter[i-7] += inter[i];
        inter[i-6] += inter[i];
    }
    res[0] = (uint8_t)_gf251_reduce32(inter[0]);
    res[1] = (uint8_t)_gf251_reduce32(inter[1]);
    res[2] = (uint8_t)_gf251_reduce32(inter[2]);
    res[3] = (uint8_t)_gf251_reduce32(inter[3]);
    res[4] = (uint8_t)_gf251_reduce32(inter[4]);
    res[5] = (uint8_t)_gf251_reduce32(inter[5]);
    res[6] = (uint8_t)_gf251_reduce32(inter[6]);
}

void gf251to7_add(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
    SET_RES_REDUCE8(6, a[6]+b[6]);
}

void gf251to7_sub(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
    SET_RES_REDUCE8(6, a[6]+PRIME-b[6]);
}

void gf251to7_mul_gf251(uint8_t res[7], uint8_t a, const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
    SET_RES_REDUCE16(6, a * b[6]);
}

int gf251to7_eq(const uint8_t a[7], const uint8_t b[7]) {
    return memeq_ct(a, b, 7);
}

#ifndef NO_RND
/*************************************************/
/************        RANDOMNESS       ************/
/*************************************************/

#include "sample.h"

void gf251_random_elements(uint8_t points[], uint32_t nb_points, samplable_t* entropy) {
    uint32_t pos = 0;
    do {
        byte_sample(entropy, points+pos, nb_points-pos);
        for(uint32_t i=pos; i<nb_points; i++) {
            if(points[i] < PRIME) {
                points[pos] = points[i];
                pos++;
            }
        }
    } while(pos < nb_points);
}

void gf251_random_elements_x4(uint8_t* const* points, uint32_t nb_points, samplable_x4_t* entropy) {
    uint32_t buffer_size = nb_points+(nb_points>>4);
    uint8_t* buffer_mem = malloc(4*buffer_size);
    uint8_t* const buffer[4] = {
        &buffer_mem[0],             &buffer_mem[  buffer_size],
        &buffer_mem[2*buffer_size], &buffer_mem[3*buffer_size]
    };
    uint32_t pos[4] = {0};
    unsigned int ok[4] = {0};
    unsigned int nb_ok = 0;
    do {
        byte_sample_x4(entropy, buffer, buffer_size);
        for(uint32_t i=0; i<buffer_size; i++) {
            if(!ok[0] && buffer[0][i] < PRIME) {
                points[0][pos[0]] = buffer[0][i];
                pos[0]++;
                if(pos[0]==nb_points) {
                    ok[0] = 1;
                    nb_ok++;
                }
            }
            if(!ok[1] && buffer[1][i] < PRIME) {
                points[1][pos[1]] = buffer[1][i];
                pos[1]++;
                if(pos[1]==nb_points) {
                    ok[1] = 1;
                    nb_ok++;
                }
            }
            if(!ok[2] && buffer[2][i] < PRIME) {
                points[2][pos[2]] = buffer[2][i];
                pos[2]++;
                if(pos[2]==nb_points) {
                    ok[2] = 1;
                    nb_ok++;
                }
            }
            if(!ok[3] && buffer[3][i] < PRIME) {
                points[3][pos[3]] = buffer[3][i];
                pos[3]++;
                if(pos[3]==nb_points) {
                    ok[3] = 1;
                    nb_ok++;
                }
            }
        }
    } while(nb_ok < 4);
    free(buffer_mem);
}
#endif
