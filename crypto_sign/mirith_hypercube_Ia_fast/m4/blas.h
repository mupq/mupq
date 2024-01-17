#ifndef BLAS_H
#define BLAS_H

#include "matrix.h"

#define BITSLICED
static inline void gf256v_set_zero(uint8_t *b, unsigned _num_byte) {
    for (uint32_t i = 0; i < _num_byte; ++i) {
        b[i] = 0;
    }
}
static inline
void gf256v_add(uint8_t *accu_b, const uint8_t *a , unsigned _num_byte) {
    while(_num_byte > 4 ) {
        *(uint32_t *)accu_b ^= *(uint32_t *)a;
        accu_b += 4;
        a += 4;
        _num_byte -= 4;
    }

    while( _num_byte-- ) { *accu_b ^= *a; accu_b++; a++; }
}


/// clobbert: acc0-acc3,t0-t3
/// only 4bits of b are used
/// acc += in*b
#define gv16v_madd_bitsclice(acc0, acc1, acc2, acc3, in0, in1, in2, in3, b, t0, t1, t2) \
	t0 = in0 ^ in3;    \
	t1 = in2 ^ in3;    \
	t2 = in1 ^ in2;    \
	if (b&1u){		   \
	acc0 = acc0 ^ in0; \
	acc1 = acc1 ^ in1; \
	acc2 = acc2 ^ in2; \
	acc3 = acc3 ^ in3; \
	}				   \
	if (b&2u) {		   \
	acc0 = acc0 ^ in3; \
	acc1 = acc1 ^ t0;  \
	acc2 = acc2 ^ in1; \
	acc3 = acc3 ^ in2; \
	}				   \
 	if (b&4u){		   \
	acc0 = acc0 ^ in2; \
	acc1 = acc1 ^ t1;  \
	acc2 = acc2 ^ t0;  \
	acc3 = acc3 ^ in1; \
	}				   \
 	if (b&8u){		   \
	acc0 = acc0 ^ in1; \
	acc1 = acc1 ^ t2;  \
	acc2 = acc2 ^ t1;  \
	acc3 = acc3 ^ t0;  \
	}


/// dont remove tmp
/// clobbert: out0-3 and tmp. in0-3 remain untouched
#define gv16v_bitslice(out0, out1, out2, out3, in0, in1, in2, in3, tmp) \
	out0 = in0 & 0x11111111; \
	out3 = in1 & 0x11111111; \
	out0 = out0 | (out3<<1u);\
	out3 = in2 & 0x11111111; \
	out0 = out0 | (out3<<2u);\
	out3 = in3 & 0x11111111; \
	out0 = out0 | (out3<<3u);\
							 \
	out1 = in1 & 0x22222222; \
	out3 = in0 & 0x22222222; \
	out1 = out1 | (out3>>1u);\
	out3 = in2 & 0x22222222; \
	out1 = out1 | (out3<<1u);\
	out3 = in3 & 0x22222222; \
	out1 = out1 | (out3<<2u);\
							 \
	out2 = in2 & 0x44444444; \
	out3 = in0 & 0x44444444; \
	out2 = out2 | (out3>>2u);\
	out3 = in1 & 0x44444444; \
	out2 = out2 | (out3>>1u);\
	out3 = in3 & 0x44444444; \
	out2 = out2 | (out3<<1u);\
							 \
	out3 = in3 & 0x88888888; \
	tmp  = in0 & 0x88888888; \
	out3 = out3 | (tmp>>3u); \
	tmp  = in1 & 0x88888888; \
	out3 = out3 | (tmp>>2u); \
	tmp  = in2 & 0x88888888; \
	out3 = out3 | (tmp>>1u);



/// small helper struct,
union tmp_32 {
    uint8_t u8[4];
    uint32_t u32;
};

union tmp_64 {
    uint8_t u8[8];
    uint32_t u32[2];
};

union tmp_128 {
    uint8_t u8[16];
    uint32_t u32[4];
};


/// taken from: https://github.com/pqov/pqov-paper/blob/main/src/gf16.h
/// \param a
/// \param b
/// \return a*b
static uint32_t gf16v_mul_u32(const uint32_t a, const uint8_t b) {
    uint32_t a_msb;
    uint32_t a32 = a;
    uint32_t b32 = b;
    uint32_t r32 = a32*(b32&1);

    a_msb = a32&0x88888888;  // MSB, 3rd bits
    a32 ^= a_msb;   // clear MSB
    a32 = (a32<<1)^((a_msb>>3)*3);
    r32 ^= (a32)*((b32>>1)&1);

    a_msb = a32&0x88888888;  // MSB, 3rd bits
    a32 ^= a_msb;   // clear MSB
    a32 = (a32<<1)^((a_msb>>3)*3);
    r32 ^= (a32)*((b32>>2)&1);

    a_msb = a32&0x88888888;  // MSB, 3rd bits
    a32 ^= a_msb;   // clear MSB
    a32 = (a32<<1)^((a_msb>>3)*3);
    r32 ^= (a32)*((b32>>3)&1);

    return r32;
}

/// \param accu_c
/// \param a
/// \param gf16_b
/// \param _num_byte  SHOULD <= 8
static inline void _gf16v_madd_u32_aligned(uint8_t *accu_c, const uint8_t *a, uint8_t gf16_b, unsigned _num_byte) {
    union tmp_32 t = {0};
    while ( _num_byte >= 4 ) {
        const uint32_t *ax = (const uint32_t *)a;
        uint32_t *cx = (uint32_t *)accu_c;
        cx[0] ^= gf16v_mul_u32( ax[0], gf16_b );
        a += 4;
        accu_c += 4;
        _num_byte -= 4;
    }

    // early exit
    if ( 0 == _num_byte ) {
        return;
    }

    for (unsigned i = 0; i < _num_byte; i++) {
        t.u8[i] = a[i];
    }

    t.u32 = gf16v_mul_u32(t.u32, gf16_b);
    for (unsigned i = 0; i < _num_byte; i++) {
        accu_c[i] ^= t.u8[i];
    }
}

static inline void _gf16v_madd_u32(uint8_t *accu_c, 
		const uint8_t *a,
		uint8_t gf16_b, 
		unsigned _num_byte) {
    uintptr_t ap = (uintptr_t)(const void *)a;
    uintptr_t cp = (uintptr_t)(const void *)accu_c;
    if ( !((cp & 3) || (ap & 3) || (_num_byte < 8)) ) {
        _gf16v_madd_u32_aligned(accu_c, a, gf16_b, _num_byte);
        return;
    }

    union tmp_32 t;

    while ( _num_byte >= 4 ) {
        t.u8[0] = a[0];
        t.u8[1] = a[1];
        t.u8[2] = a[2];
        t.u8[3] = a[3];
        t.u32 = gf16v_mul_u32(t.u32, gf16_b);
        accu_c[0] ^= t.u8[0];
        accu_c[1] ^= t.u8[1];
        accu_c[2] ^= t.u8[2];
        accu_c[3] ^= t.u8[3];
        a += 4;
        accu_c += 4;
        _num_byte -= 4;
    }
    if ( 0 == _num_byte ) {
        return;
    }

    for (unsigned i = 0; i < _num_byte; i++) {
        t.u8[i] = a[i];
    }
    t.u32 = gf16v_mul_u32(t.u32, gf16_b);
    for (unsigned i = 0; i < _num_byte; i++) {
        accu_c[i] ^= t.u8[i];
    }

}

/// taken from: https://github.com/pqov/pqov-paper/blob/main/src/blas_comm.h#L18
/// @brief get an element from GF(16) vector .
///
/// @param[in]  a         - the input vector a.
/// @param[in]  i         - the index in the vector a.
/// @return  the value of the element.
///
static inline uint8_t gf16v_get_ele(const uint8_t *a, const unsigned i) {
    uint8_t r = a[i >> 1];
    return (i & 1) ? (r >> 4) : (r & 0xf);
}

static void gf16mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    for (unsigned i = 0; i < n_A_width; i++) {
        uint8_t bb = gf16v_get_ele(b, i);
        _gf16v_madd_u32(c, matA, bb, n_A_vec_byte);
        matA += n_A_vec_byte;
    }
}

/// \param matC
/// \param matA
/// \param height_A
/// \param width_A_byte
/// \param matB
/// \param width_B_byte
static void gf16mat_intmat_mul_ref(uint8_t *matC,
                                   const uint8_t *matA,
                                   unsigned nr_bytes,
                                   const uint8_t b)
{
    uint32_t *matA_32 = (uint32_t *)matA;
    uint32_t *matC_32 = (uint32_t *)matC;
    uint32_t k=0;

#ifdef BITSLICED
    uint32_t out0,out1,out2,out3,tmp0,tmp1,tmp2;
    union tmp_128 t128={0};

    for( ;(k+16)<nr_bytes; k+=16) {
        for (uint16_t i = 0; i < 4; ++i) {
            t128.u32[i] = matA_32[i];
        }
        uint32_t int0=0,int1=0,int2=0,int3=0;
        gv16v_bitslice(out0,out1,out2,out3,t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],tmp0)
        gv16v_madd_bitsclice(int0,int1,int2,int3,out0,out1,out2,out3,b,tmp0,tmp1,tmp2)
        gv16v_bitslice(t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],int0,int1,int2,int3,tmp0)
        for (uint16_t i = 0; i < 4; ++i) {
            matC_32[i] ^= t128.u32[i];
        }

        matC_32 += 4;
        matA_32 += 4;
    }
#endif
    for(; (k+4)<nr_bytes; k+=4) {
        *matC_32 ^= gf16v_mul_u32(*matA_32, b);

        matC_32 += 1;
        matA_32 += 1;
    }

    union tmp_32 t32 = {0};
    for (unsigned i = 0; i < nr_bytes - k; i++) {
        t32.u8[i] = matA[k + i];
    }

    t32.u32 = gf16v_mul_u32(t32.u32, b);
    for (unsigned i = 0; i < (nr_bytes-k); i++) {
        matC[k + i] ^= t32.u8[i];
    }
}

static void gf16mat_rowmat_mul_ref(uint8_t *matC,
                                   const uint8_t *matA,
                                   unsigned height_A,
                                   unsigned width_A_byte,
                                   const uint8_t *matB,
                                   unsigned width_B_byte) {
//#ifdef BITSLICED
//    uint32_t out0,out1,out2,out3,tmp0,tmp1,tmp2;
//
//    for(uint32_t k=0; k<height_A; k++) {
//        union tmp_128 t128= {0};
//        for (uint32_t i = 0; i < width_B_byte; ++i) {
//            t128.u8[i] = matB[i];
//        }
//        uint32_t int0=0,int1=0,int2=0,int3=0;
//        gv16v_bitslice(out0,out1,out2,out3,t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],tmp0);
//        for (uint32_t i = 0; i < width_A_byte*2; ++i) {
//            uint8_t bb = gf16v_get_ele(matA, i);
//	        gv16v_madd_bitsclice(int0,int1,int2,int3,out0,out1,out2,out3,bb,tmp0,tmp1,tmp2);
//        }
//
//        gv16v_bitslice(t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],int0,int1,int2,int3,tmp0);
//        for (uint32_t i = 0; i < width_B_byte; ++i) {
//            matC[i] ^= t128.u8[i];
//        }
//        matC += width_B_byte;
//        matA += width_A_byte;
//    }
//#else
    for( unsigned i=0; i<height_A; i++) {
        gf16mat_prod_ref(matC, matB , width_B_byte , width_A_byte*2 , matA );
        matC += width_B_byte;
        matA += width_A_byte;
    }
//#endif
}



static void gf16mat_colmat_mul_ref(uint8_t *mat_c, 
		const uint8_t *mat_a, 
		unsigned a_veclen_byte, 
		unsigned a_n_vec, 
		const uint8_t *mat_b, 
		unsigned b_n_vec) {
    gf16mat_rowmat_mul_ref( mat_c , mat_b , b_n_vec , (a_n_vec+1)>>1 , mat_a , a_veclen_byte );
}

void _matrix_product(ff_t *result, const ff_t *matrix1, const ff_t *matrix2,
    const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    const uint32_t n_bytes_per_column1 = matrix_bytes_per_column(n_rows1);
    gf256v_set_zero(result , n_bytes_per_column1*n_cols2);
    gf16mat_colmat_mul_ref(result, matrix1, n_bytes_per_column1, n_cols1, matrix2, n_cols2);
    return;
}
void _matrix_add(ff_t *matrix1, const ff_t *matrix2, 
		const uint32_t n_rows, const uint32_t n_cols) {
    const unsigned n_bytes = matrix_bytes_size(n_rows, n_cols);
    gf256v_add(matrix1, matrix2, n_bytes);
}

void _matrix_add_multiple(ff_t *matrix1, ff_t scalar, const ff_t *matrix2,
    const uint32_t n_rows, const uint32_t n_cols) {
    const int n_bytes = matrix_bytes_size(n_rows, n_cols);
    gf16mat_intmat_mul_ref(matrix1, matrix2, n_bytes, scalar);
    return;
}

void _matrix_add_product(ff_t *matrix1, const ff_t *matrix2, const ff_t *matrix3,
                        const uint32_t n_rows1, const uint32_t n_cols1, const uint32_t n_cols2) {
    const uint32_t n_bytes_per_column1 = matrix_bytes_per_column(n_rows1);
    gf16mat_rowmat_mul_ref(matrix1, matrix2, n_cols2, n_bytes_per_column1, matrix3, 2*((n_cols1+1)>>1));
    return;
}
#endif
