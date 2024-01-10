#include "gf31.h"
#include <stdlib.h>
#include <string.h>

#include "gf31-internal.h"

static int memeq_ct(const uint8_t* a, const uint8_t* b, uint16_t n) {
    int res = 1;
    for(uint16_t i=0; i<n; i++)
        res &= (a[i] == b[i]);
    return res;
}

/*************************************************/
/***********        COMPRESSION        ***********/
/*************************************************/

void gf31_compress_tab(uint8_t* buf, const uint8_t* v, uint32_t size) {
    unsigned int ind = 0;
    unsigned int pos = 0;
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<5; j++) {
            uint8_t bit = (v[i] >> j) & 0x01;
            if(pos == 0)
                buf[ind] = bit;
            else
                buf[ind] |= bit << pos;
            pos++;
            if(pos == 8) {
                pos = 0;
                ind++;
            } 
        }
    }
}

void gf31_decompress_tab(uint8_t* v, const uint8_t* buf, uint32_t size) {
    unsigned int ind = 0;
    unsigned int pos = 0;
    for(unsigned int i=0; i<size; i++) {
        v[i] = 0;
        for(unsigned int j=0; j<5; j++) {
            uint8_t bit = (buf[ind] >> pos) & 0x01;
            v[i] |= bit << j;
            pos++;
            if(pos == 8) {
                pos = 0;
                ind++;
            } 
        }
    }
}

/*************************************************/
/***********      FIELD OPERATIONS     ***********/
/*************************************************/

uint8_t gf31_add(uint8_t a, uint8_t b) {
    uint16_t res = a+b;
    res -= PRIME*(PRIME <= res);
    return (uint8_t) res;
}

uint8_t gf31_sub(uint8_t a, uint8_t b) {
    uint16_t res = a+PRIME-b;
    res -= PRIME*(PRIME <= res);
    return (uint8_t) res;
}

uint8_t gf31_mul(uint8_t a, uint8_t b) {
    uint16_t res = a*b;
    return (uint8_t) _gf31_reduce16(res);
}

uint8_t gf31_neg(uint8_t a) {
    return (a != 0)*(PRIME - a);
}

uint8_t gf31_innerproduct(const uint8_t* x, const uint8_t* y, uint32_t size) {
    uint32_t z = 0;
    for(uint32_t i=0; i<size; i++)
        z += x[i]*y[i];
    return (uint8_t) _gf31_reduce32(z);
}

/*************************************************/
/***********     BATCHED OPERATIONS    ***********/
/*************************************************/

// vz[] += vx[]
void gf31_add_tab(uint8_t* vz, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i]+vx[i];
        res -= PRIME*(PRIME <= res);
        vz[i] = (uint8_t) res;
    }
}

// vz[] -= vx[]
void gf31_sub_tab(uint8_t* vz, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i]+PRIME-vx[i];
        res -= PRIME*(PRIME <= res);
        vz[i] = (uint8_t) res;
    }
}

// vz[] = vx[] * y
void gf31_mul_tab(uint8_t* vz, const uint8_t* vx, uint8_t y, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vx[i]*y;
        vz[i] = (uint8_t)_gf31_reduce16(res);
    }
}

// vz[] += vx[] * y
void gf31_muladd_tab(uint8_t* vz, uint8_t y, const uint8_t* vx, uint32_t size) {
    uint16_t res;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i] + vx[i]*y;
        vz[i] = (uint8_t)_gf31_reduce16(res);
    }
}

// vz[] = -vz[]
void gf31_neg_tab(uint8_t* vz, uint32_t size) {
    for(uint32_t i=0; i<size; i++)
        vz[i] = (vz[i] != 0)*(PRIME - vz[i]);
}

/*************************************************/
/***********      FIELD EXTENSION      ***********/
/*************************************************/

#define SET_RES_REDUCE8(i,exp) {v=exp; res[i]=(uint8_t)(v-PRIME*(PRIME<=v));}
#define SET_RES_REDUCE16(i,exp) {v=exp; res[i]=(uint8_t)_gf31_reduce16(v);}

// X^6 - 3 = 0
void gf31to6_mul(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]) {
    uint32_t inter[11] = {0};
    unsigned int i, j;
    for(i=0; i<6; i++)
        for(j=0; j<6; j++)
            inter[i+j] += a[i]*b[j];
    for(i=6; i<11; i++)
        inter[i-6] += 3*inter[i];
    res[0] = (uint8_t)_gf31_reduce16(inter[0]);
    res[1] = (uint8_t)_gf31_reduce16(inter[1]);
    res[2] = (uint8_t)_gf31_reduce16(inter[2]);
    res[3] = (uint8_t)_gf31_reduce16(inter[3]);
    res[4] = (uint8_t)_gf31_reduce16(inter[4]);
    res[5] = (uint8_t)_gf31_reduce16(inter[5]);
}

void gf31to6_add(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
}

void gf31to6_sub(uint8_t res[6], const uint8_t a[6], const uint8_t b[6]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
}

void gf31to6_mul_gf31(uint8_t res[6], uint8_t a, const uint8_t b[6]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
}

int gf31to6_eq(const uint8_t a[6], const uint8_t b[6]) {
    return memeq_ct(a, b, 6);
}


// X^7 - X - 3 = 0
void gf31to7_mul(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint32_t inter[13] = {0};
    unsigned int i, j;
    for(i=0; i<7; i++)
        for(j=0; j<7; j++)
            inter[i+j] += a[i]*b[j];
    for(i=12; i>=7; i--) {
        inter[i-7] += 3*inter[i];
        inter[i-6] += inter[i];
    }
    res[0] = (uint8_t)_gf31_reduce16(inter[0]);
    res[1] = (uint8_t)_gf31_reduce16(inter[1]);
    res[2] = (uint8_t)_gf31_reduce16(inter[2]);
    res[3] = (uint8_t)_gf31_reduce16(inter[3]);
    res[4] = (uint8_t)_gf31_reduce16(inter[4]);
    res[5] = (uint8_t)_gf31_reduce16(inter[5]);
    res[6] = (uint8_t)_gf31_reduce16(inter[6]);
}

void gf31to7_add(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
    SET_RES_REDUCE8(6, a[6]+b[6]);
}

void gf31to7_sub(uint8_t res[7], const uint8_t a[7], const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
    SET_RES_REDUCE8(6, a[6]+PRIME-b[6]);
}

void gf31to7_mul_gf31(uint8_t res[7], uint8_t a, const uint8_t b[7]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
    SET_RES_REDUCE16(6, a * b[6]);
}

int gf31to7_eq(const uint8_t a[7], const uint8_t b[7]) {
    return memeq_ct(a, b, 7);
}

// X^8 - X^2 - 1 = 0
void gf31to8_mul(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]) {
    uint32_t inter[15] = {0};
    unsigned int i, j;
    for(i=0; i<8; i++)
        for(j=0; j<8; j++)
            inter[i+j] += a[i]*b[j];
    for(i=14; i>=8; i--) {
        inter[i-8] += inter[i];
        inter[i-6] += inter[i];
    }
    res[0] = (uint8_t)_gf31_reduce16(inter[0]);
    res[1] = (uint8_t)_gf31_reduce16(inter[1]);
    res[2] = (uint8_t)_gf31_reduce16(inter[2]);
    res[3] = (uint8_t)_gf31_reduce16(inter[3]);
    res[4] = (uint8_t)_gf31_reduce16(inter[4]);
    res[5] = (uint8_t)_gf31_reduce16(inter[5]);
    res[6] = (uint8_t)_gf31_reduce16(inter[6]);
    res[7] = (uint8_t)_gf31_reduce16(inter[7]);
}

void gf31to8_add(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
    SET_RES_REDUCE8(6, a[6]+b[6]);
    SET_RES_REDUCE8(7, a[7]+b[7]);
}

void gf31to8_sub(uint8_t res[8], const uint8_t a[8], const uint8_t b[8]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
    SET_RES_REDUCE8(6, a[6]+PRIME-b[6]);
    SET_RES_REDUCE8(7, a[7]+PRIME-b[7]);
}

void gf31to8_mul_gf31(uint8_t res[8], uint8_t a, const uint8_t b[8]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
    SET_RES_REDUCE16(6, a * b[6]);
    SET_RES_REDUCE16(7, a * b[7]);
}

int gf31to8_eq(const uint8_t a[8], const uint8_t b[8]) {
    return memeq_ct(a, b, 8);
}

// X^10 - 3 = 0
void gf31to10_mul(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]) {
    uint32_t inter[19] = {0};
    unsigned int i, j;
    for(i=0; i<10; i++)
        for(j=0; j<10; j++)
            inter[i+j] += a[i]*b[j];
    for(i=10; i<19; i++)
        inter[i-10] += 3*inter[i];
    res[0] = (uint8_t)_gf31_reduce16(inter[0]);
    res[1] = (uint8_t)_gf31_reduce16(inter[1]);
    res[2] = (uint8_t)_gf31_reduce16(inter[2]);
    res[3] = (uint8_t)_gf31_reduce16(inter[3]);
    res[4] = (uint8_t)_gf31_reduce16(inter[4]);
    res[5] = (uint8_t)_gf31_reduce16(inter[5]);
    res[6] = (uint8_t)_gf31_reduce16(inter[6]);
    res[7] = (uint8_t)_gf31_reduce16(inter[7]);
    res[8] = (uint8_t)_gf31_reduce16(inter[8]);
    res[9] = (uint8_t)_gf31_reduce16(inter[9]);
}

void gf31to10_add(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
    SET_RES_REDUCE8(6, a[6]+b[6]);
    SET_RES_REDUCE8(7, a[7]+b[7]);
    SET_RES_REDUCE8(8, a[8]+b[8]);
    SET_RES_REDUCE8(9, a[9]+b[9]);
}

void gf31to10_sub(uint8_t res[10], const uint8_t a[10], const uint8_t b[10]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
    SET_RES_REDUCE8(6, a[6]+PRIME-b[6]);
    SET_RES_REDUCE8(7, a[7]+PRIME-b[7]);
    SET_RES_REDUCE8(8, a[8]+PRIME-b[8]);
    SET_RES_REDUCE8(9, a[9]+PRIME-b[9]);
}

void gf31to10_mul_gf31(uint8_t res[10], uint8_t a, const uint8_t b[10]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
    SET_RES_REDUCE16(6, a * b[6]);
    SET_RES_REDUCE16(7, a * b[7]);
    SET_RES_REDUCE16(8, a * b[8]);
    SET_RES_REDUCE16(9, a * b[9]);
}

int gf31to10_eq(const uint8_t a[10], const uint8_t b[10]) {
    return memeq_ct(a, b, 10);
}

// X^11 - X^3 - 1 = 0
void gf31to11_mul(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]) {
    uint32_t inter[21] = {0};
    unsigned int i, j;
    for(i=0; i<11; i++)
        for(j=0; j<11; j++)
            inter[i+j] += a[i]*b[j];
    for(i=20; i>=11; i--) {
        inter[i-11] += inter[i];
        inter[i-8] += inter[i];
    }
    res[0] = (uint8_t)_gf31_reduce16(inter[0]);
    res[1] = (uint8_t)_gf31_reduce16(inter[1]);
    res[2] = (uint8_t)_gf31_reduce16(inter[2]);
    res[3] = (uint8_t)_gf31_reduce16(inter[3]);
    res[4] = (uint8_t)_gf31_reduce16(inter[4]);
    res[5] = (uint8_t)_gf31_reduce16(inter[5]);
    res[6] = (uint8_t)_gf31_reduce16(inter[6]);
    res[7] = (uint8_t)_gf31_reduce16(inter[7]);
    res[8] = (uint8_t)_gf31_reduce16(inter[8]);
    res[9] = (uint8_t)_gf31_reduce16(inter[9]);
    res[10] = (uint8_t)_gf31_reduce16(inter[10]);
}

void gf31to11_add(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+b[0]);
    SET_RES_REDUCE8(1, a[1]+b[1]);
    SET_RES_REDUCE8(2, a[2]+b[2]);
    SET_RES_REDUCE8(3, a[3]+b[3]);
    SET_RES_REDUCE8(4, a[4]+b[4]);
    SET_RES_REDUCE8(5, a[5]+b[5]);
    SET_RES_REDUCE8(6, a[6]+b[6]);
    SET_RES_REDUCE8(7, a[7]+b[7]);
    SET_RES_REDUCE8(8, a[8]+b[8]);
    SET_RES_REDUCE8(9, a[9]+b[9]);
    SET_RES_REDUCE8(10, a[10]+b[10]);
}

void gf31to11_sub(uint8_t res[11], const uint8_t a[11], const uint8_t b[11]) {
    uint16_t v;
    SET_RES_REDUCE8(0, a[0]+PRIME-b[0]);
    SET_RES_REDUCE8(1, a[1]+PRIME-b[1]);
    SET_RES_REDUCE8(2, a[2]+PRIME-b[2]);
    SET_RES_REDUCE8(3, a[3]+PRIME-b[3]);
    SET_RES_REDUCE8(4, a[4]+PRIME-b[4]);
    SET_RES_REDUCE8(5, a[5]+PRIME-b[5]);
    SET_RES_REDUCE8(6, a[6]+PRIME-b[6]);
    SET_RES_REDUCE8(7, a[7]+PRIME-b[7]);
    SET_RES_REDUCE8(8, a[8]+PRIME-b[8]);
    SET_RES_REDUCE8(9, a[9]+PRIME-b[9]);
    SET_RES_REDUCE8(10, a[10]+PRIME-b[10]);
}

void gf31to11_mul_gf31(uint8_t res[11], uint8_t a, const uint8_t b[11]) {
    uint16_t v;
    SET_RES_REDUCE16(0, a * b[0]);
    SET_RES_REDUCE16(1, a * b[1]);
    SET_RES_REDUCE16(2, a * b[2]);
    SET_RES_REDUCE16(3, a * b[3]);
    SET_RES_REDUCE16(4, a * b[4]);
    SET_RES_REDUCE16(5, a * b[5]);
    SET_RES_REDUCE16(6, a * b[6]);
    SET_RES_REDUCE16(7, a * b[7]);
    SET_RES_REDUCE16(8, a * b[8]);
    SET_RES_REDUCE16(9, a * b[9]);
    SET_RES_REDUCE16(10, a * b[10]);
}

int gf31to11_eq(const uint8_t a[11], const uint8_t b[11]) {
    return memeq_ct(a, b, 11);
}


#ifndef NO_RND
/*************************************************/
/************        RANDOMNESS       ************/
/*************************************************/

#include "sample.h"

void gf31_random_elements(uint8_t points[], uint32_t nb_points, samplable_t* entropy) {
    uint32_t buffer_size = 5*((nb_points+(nb_points>>4)+7)>>3);
    uint8_t* buffer = malloc(buffer_size);
    uint32_t pos = 0;
    uint8_t buf[8];
    uint32_t i,j;
    uint8_t* cur;
    unsigned int ok=0;
    do {
        byte_sample(entropy, buffer, buffer_size);
        cur = buffer;
        for(i=0; i<buffer_size && !ok; i+=5) {
            buf[0] = cur[0] & 0x1F;
            buf[1] = (uint8_t) ((cur[0]>>5) | (cur[1] & 0x03)<<3);
            buf[2] = (cur[1]>>2) & 0x1F;
            buf[3] = (uint8_t) ((cur[1]>>7) | (cur[2] & 0x0F)<<1);
            buf[4] = (uint8_t) ((cur[2]>>4) | (cur[3] & 0x01)<<4);
            buf[5] = (cur[3]>>1) & 0x1F;
            buf[6] = (uint8_t) ((cur[3]>>6) | (cur[4] & 0x07)<<2);
            buf[7] = cur[4]>>3;
            for(j=0;j<8 && !ok;j++) {
                if(buf[j] < PRIME) {
                    points[pos] = buf[j];
                    pos++;
                    if(pos==nb_points)
                        ok = 1;
                }
            }
            cur+=5;
        }
    } while(!ok);
    free(buffer);
}

void gf31_random_elements_x4(uint8_t* const* points, uint32_t nb_points, samplable_x4_t* entropy) {
    uint32_t buffer_size = 5*((nb_points+(nb_points>>4)+7)>>3);
    uint8_t* buffer_mem = malloc(4*buffer_size);
    uint8_t* buffer[4] = {
        &buffer_mem[0],             &buffer_mem[  buffer_size],
        &buffer_mem[2*buffer_size], &buffer_mem[3*buffer_size]
    };
    uint32_t pos[4] = {0};
    unsigned int ok[4] = {0};
    unsigned int nb_ok = 0;
    uint8_t buf[8];
    uint32_t i,j,k;
    uint8_t* cur;
    do {
        byte_sample_x4(entropy, buffer, buffer_size);
        for(k=0; k<4; k++) {
            cur = buffer[k];
            for(i=0; i<buffer_size && !ok[k]; i+=5) {
                buf[0] = cur[0] & 0x1F;
                buf[1] = (uint8_t) ((cur[0]>>5) | (cur[1] & 0x03)<<3);
                buf[2] = (cur[1]>>2) & 0x1F;
                buf[3] = (uint8_t) ((cur[1]>>7) | (cur[2] & 0x0F)<<1);
                buf[4] = (uint8_t) ((cur[2]>>4) | (cur[3] & 0x01)<<4);
                buf[5] = (cur[3]>>1) & 0x1F;
                buf[6] = (uint8_t) ((cur[3]>>6) | (cur[4] & 0x07)<<2);
                buf[7] = cur[4]>>3;
                for(j=0;j<8 && !ok[k];j++) {
                    if(buf[j] < PRIME) {
                        points[k][pos[k]] = buf[j];
                        pos[k]++;
                        if(pos[k]==nb_points) {
                            ok[k] = 1;
                            nb_ok++;
                        }
                    }
                }
                cur+=5;
            }
        }
    } while(nb_ok < 4);
    free(buffer_mem);
}
#endif
