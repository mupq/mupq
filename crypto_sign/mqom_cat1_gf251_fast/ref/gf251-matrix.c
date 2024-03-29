#include "gf251.h"
#include <stdlib.h>
#include <string.h>

#include "gf251-internal.h"

/*************************************************/
/***********     MATRIX OPERATIONS    ************/
/*************************************************/

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb]
void gf251_matcols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size) {
    uint32_t unreduced_vz[size];
    memset(unreduced_vz, 0, sizeof(uint32_t)*size);
    size_t ind=0;
    for(uint32_t j=0; j<nb; j++)
        for(uint32_t i=0; i<size; i++)
            unreduced_vz[i] += vx[ind++]*y[j];
    for(uint32_t i=0; i<size; i++)
        vz[i] = (uint8_t)_gf251_reduce32(unreduced_vz[i] + vz[i]);
}

// vz[] += vx[1][] * y[1] + ... + vx[nb][] * y[nb]
void gf251_matrows_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size) {
    uint32_t res;
    size_t ind=0;
    for(uint32_t i=0; i<size; i++) {
        res = vz[i];
        for(uint32_t j=0; j<nb; j++)
            res += vx[ind++]*y[j];
        vz[i] = (uint8_t)_gf251_reduce32(res);
    }
}

void gf251_mat128cols_muladd(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t size) {
    gf251_matcols_muladd(vz, y, vx, nb, size);
}

// vz[] += vx[][1] * y[1] + ... + vx[][nb] * y[nb]
void gf251_matcols_muladd_triangular(uint8_t* vz, const uint8_t* y, const uint8_t* vx, uint32_t nb, uint32_t step) {
    uint32_t size = nb*step;
    uint32_t unreduced_vz[size];
    memset(unreduced_vz, 0, sizeof(uint32_t)*size);
    size_t ind=0;
    for(uint32_t j=0; j<nb; j++)
        for(uint32_t i=j*step; i<size; i++)
            unreduced_vz[i] += vx[ind++]*y[j];
    for(uint32_t i=0; i<size; i++)
        vz[i] = (uint8_t)_gf251_reduce32(unreduced_vz[i] + vz[i]);
}
