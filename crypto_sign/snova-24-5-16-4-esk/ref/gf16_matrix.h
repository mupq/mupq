#ifndef GF16M_H
#define GF16M_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "deriv_params.h"
#include "gf16.h"

#if rank == 2
#define get_gf16m(gf16m, x, y) (gf16m[(((x) << 1) ^ (y))])
#define set_gf16m(gf16m, x, y, value) (gf16m[(((x) << 1) ^ (y))] = value)
#elif rank == 3
#define get_gf16m(gf16m, x, y) (gf16m[(((x)*3) + (y))])
#define set_gf16m(gf16m, x, y, value) (gf16m[(((x)*3) + (y))] = value)
#elif rank == 4
#define get_gf16m(gf16m, x, y) (gf16m[(((x) << 2) ^ (y))])
#define set_gf16m(gf16m, x, y, value) (gf16m[(((x) << 2) ^ (y))] = value)
#endif

typedef gf16_t gf16m_t[sq_rank];

#endif