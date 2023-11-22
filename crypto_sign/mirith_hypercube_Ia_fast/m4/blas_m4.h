/* Based on the public domain implementation in pqov-paper/src/avx2/ from
 * https://github.com/pqov/pqov-paper by Ward Beullens, Ming-Shing Chen,
 * Shih-Hao Hung, Matthias J. Kannwischer, Bo-Yuan Peng, Cheng-Jhih Shih,
 * and Bo-Yin Yang */

#ifndef MIRITH_DSS_BLAS_AVX2_H
#define MIRITH_DSS_BLAS_AVX2_H

#include <stdint.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include "matrix.h"
#include "matrix_constants.h"

/// small helper struct,
union tmp_32 {
    uint8_t u8[4];
    uint32_t u32;
} t;

union tmp_64 {
    uint8_t u8[8];
    uint32_t u32[2];
};

union tmp_128 {
    uint8_t u8[16];
    uint32_t u32[4];
};

static ff_t mult_table[256] = {
        0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        0x00,0x02,0x04,0x06,0x08,0x0a,0x0c,0x0e,0x03,0x01,0x07,0x05,0x0b,0x09,0x0f,0x0d,
        0x00,0x03,0x06,0x05,0x0c,0x0f,0x0a,0x09,0x0b,0x08,0x0d,0x0e,0x07,0x04,0x01,0x02,
        0x00,0x04,0x08,0x0c,0x03,0x07,0x0b,0x0f,0x06,0x02,0x0e,0x0a,0x05,0x01,0x0d,0x09,
        0x00,0x05,0x0a,0x0f,0x07,0x02,0x0d,0x08,0x0e,0x0b,0x04,0x01,0x09,0x0c,0x03,0x06,
        0x00,0x06,0x0c,0x0a,0x0b,0x0d,0x07,0x01,0x05,0x03,0x09,0x0f,0x0e,0x08,0x02,0x04,
        0x00,0x07,0x0e,0x09,0x0f,0x08,0x01,0x06,0x0d,0x0a,0x03,0x04,0x02,0x05,0x0c,0x0b,
        0x00,0x08,0x03,0x0b,0x06,0x0e,0x05,0x0d,0x0c,0x04,0x0f,0x07,0x0a,0x02,0x09,0x01,
        0x00,0x09,0x01,0x08,0x02,0x0b,0x03,0x0a,0x04,0x0d,0x05,0x0c,0x06,0x0f,0x07,0x0e,
        0x00,0x0a,0x07,0x0d,0x0e,0x04,0x09,0x03,0x0f,0x05,0x08,0x02,0x01,0x0b,0x06,0x0c,
        0x00,0x0b,0x05,0x0e,0x0a,0x01,0x0f,0x04,0x07,0x0c,0x02,0x09,0x0d,0x06,0x08,0x03,
        0x00,0x0c,0x0b,0x07,0x05,0x09,0x0e,0x02,0x0a,0x06,0x01,0x0d,0x0f,0x03,0x04,0x08,
        0x00,0x0d,0x09,0x04,0x01,0x0c,0x08,0x05,0x02,0x0f,0x0b,0x06,0x03,0x0e,0x0a,0x07,
        0x00,0x0e,0x0f,0x01,0x0d,0x03,0x02,0x0c,0x09,0x07,0x06,0x08,0x04,0x0a,0x0b,0x05,
        0x00,0x0f,0x0d,0x02,0x09,0x06,0x04,0x0b,0x01,0x0e,0x0c,0x03,0x08,0x07,0x05,0x0a,
};

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

/// taken from: https://github.com/pqov/pqov-paper/blob/main/src/gf16.h
///  return v[0]^v[1]^v[2]^v[3]
static uint8_t gf256v_reduce_u32(uint32_t a) {
// https://godbolt.org/z/7hirMb
    uint16_t *aa = (uint16_t *) (&a);
    uint16_t r = aa[0] ^ aa[1];
    uint8_t *rr = (uint8_t *) (&r);
    return rr[0] ^ rr[1];
}

/// taken from: https://github.com/pqov/pqov-paper/blob/main/src/gf16.h
/// \param a
/// \return
static inline uint8_t gf16v_reduce_u32(uint32_t a) {
    uint8_t r256 = gf256v_reduce_u32(a);
    return (r256 & 0xf) ^ (r256 >> 4);
}


/// taken from: https://github.com/pqov/pqov-paper/blob/main/src/gf16.h
/// \param a
/// \return a**2
static inline uint32_t gf16v_squ_u32(const uint32_t a) {
    uint32_t a01 = (a & 0x11111111) + ((a<<1) & 0x44444444);
    uint32_t a23 = (((a>>2)&0x11111111)+((a>>1)&0x44444444))*3;
    return a01^a23;
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
	acc1 = acc1 ^ tmp0;\
	acc2 = acc2 ^ in1; \
	acc3 = acc3 ^ in2; \
	}				   \
 	if (b&4u){		   \
	acc0 = acc0 ^ in2; \
	acc1 = acc1 ^ tmp1;\
	acc2 = acc2 ^ tmp0;\
	acc3 = acc3 ^ in1; \
	}				   \
 	if (b&8u){		   \
	acc0 = acc0 ^ in1; \
	acc1 = acc1 ^ tmp2;\
	acc2 = acc2 ^ tmp1;\
	acc3 = acc3 ^ tmp0;\
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


static inline
void gf256v_add_neon( uint8_t * accu_b, const uint8_t * a , unsigned _num_byte ) {
    while(_num_byte > 4 ) {
        *(uint32_t *)accu_b ^= *(uint32_t *)a;
        accu_b += 4;
        a += 4;
        _num_byte -= 4;
    }
    while( _num_byte-- ) { *accu_b ^= *a; accu_b++; a++; }
}

/////////////////////////////from blas_common.c//////////////////////////////
#define gf256v_add          gf256v_add_neon
static void gf256v_set_zero(uint8_t *b, unsigned _num_byte)
{
    //gf256v_add(b, b, _num_byte);
    for (uint32_t i = 0; i < _num_byte; ++i) {
        b[i] = 0;
    }
}
//////////////////////////////////////////////////////////////////////////////

/////////////////////////  matrix-matrix multiplication  //////////////////////

/// \param accu_c
/// \param a
/// \param gf16_b
/// \param _num_byte  SHOULD <= 8
static inline void _gf16v_madd_u32_aligned(uint8_t *accu_c, const uint8_t *a, uint8_t gf16_b, unsigned _num_byte) {
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



static inline void _gf16v_madd_u32(uint8_t *accu_c, const uint8_t *a, uint8_t gf16_b, unsigned _num_byte) {

    uintptr_t ap = (uintptr_t)(const void *)a;
    uintptr_t cp = (uintptr_t)(const void *)accu_c;
    if ( !((cp & 3) || (ap & 3) || (_num_byte < 8)) ) {
        _gf16v_madd_u32_aligned(accu_c, a, gf16_b, _num_byte);
        return;
    }

    union tmp_32 {
        uint8_t u8[4];
        uint32_t u32;
    } t;

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

static void gf16mat_prod_ref(uint8_t *c, const uint8_t *matA, unsigned n_A_vec_byte, unsigned n_A_width, const uint8_t *b) {
    for (unsigned i = 0; i < n_A_width; i++) {
        uint8_t bb = gf16v_get_ele(b, i);
        _gf16v_madd_u32(c, matA, bb, n_A_vec_byte);
        matA += n_A_vec_byte;
    }
}

//////////// specialized functions  /////////////////////

#define BITSLICED

///
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
    // TODO in both cases memory overread in A and C
#ifdef BITSLICED
    uint32_t out0,out1,out2,out3,tmp0,tmp1,tmp2;
    union tmp_128 t128={0};

    uint32_t k=0;
    for( ;k+16<nr_bytes; k+=16) {
        for (uint16_t i = 0; i < 4; ++i) {
            t128.u32[i] = matA_32[i];
        }
        uint32_t int0=0,int1=0,int2=0,int3=0;
        gv16v_bitslice(out0,out1,out2,out3,t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],tmp0)
        gv16v_madd_bitsclice(int0,int1,int2,int3,out0,out1,out2,out3,b,tmp0,tmp1,tmp2)
        gv16v_bitslice(t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],int0,int1,int2,int3,tmp0)
        for (uint16_t i = 0; i < 4; ++i) {
            matC_32[i] ^= t128.u8[i];
        }
        matC_32 += 4;
        matA_32 += 4;
    }

    for(; k<nr_bytes; k+=4) {
        *matC_32 = gf16v_mul_u32(*matA_32, b);

        matC_32 += 1;
        matA_32 += 1;
    }

#else
    union tmp_32 t32 = {0};
    for(uint32_t k=0; k<nr_bytes; k+=4) {
        t32.u32 = *matA_32;
        t32.u32 = gf16v_mul_u32(t.u32, b);
        *matC_32 ^= t32.u32;

        matC_32 += 1;
        matA_32 += 1;
    }
#endif
}


///
/// \param matC
/// \param matA
/// \param height_A
/// \param width_A_byte
/// \param matB
/// \param width_B_byte
static void gf16mat_rowmat_mul_ref(uint8_t *matC,
                                   const uint8_t *matA,
                                   unsigned height_A,
                                   unsigned width_A_byte,
                                   const uint8_t *matB,
                                   unsigned width_B_byte)
{
    gf256v_set_zero(matC , height_A*width_B_byte);
#ifdef BITSLICED
    uint32_t out0,out1,out2,out3,tmp0,tmp1,tmp2;
    union tmp_128 t128={0};

    for(uint32_t k=0; k<height_A; k++) {
        for (uint16_t i = 0; i < width_B_byte; ++i) {
            t128.u8[i] = matB[i];
        }
        uint32_t int0=0,int1=0,int2=0,int3=0;
        gv16v_bitslice(out0,out1,out2,out3,t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],tmp0)
        for (uint32_t i = 0; i < width_A_byte*2; ++i) {
            uint8_t bb = gf16v_get_ele(matA, i);
            gv16v_madd_bitsclice(int0,int1,int2,int3,out0,out1,out2,out3,bb,tmp0,tmp1,tmp2)
        }

        gv16v_bitslice(t128.u32[0],t128.u32[1],t128.u32[2],t128.u32[3],int0,int1,int2,int3,tmp0)
        for (uint16_t i = 0; i < width_B_byte; ++i) {
            matC[i] ^= t128.u8[i];
        }
        matC += width_B_byte;
        matA += width_A_byte;
    }


#else

    for( unsigned i=0; i<height_A; i++) {
        gf16mat_prod_ref(matC, matB , width_B_byte , width_A_byte*2 , matA );
        matC += width_B_byte;
        matA += width_A_byte;
    }
#endif
}



static void gf16mat_colmat_mul_ref(uint8_t *mat_c, const uint8_t *mat_a, unsigned a_veclen_byte, unsigned a_n_vec, const uint8_t *mat_b, unsigned b_n_vec)
{
    gf16mat_rowmat_mul_ref( mat_c , mat_b , b_n_vec , (a_n_vec+1)>>1 , mat_a , a_veclen_byte );
}

///////////////////////////////////////////////////////////////////////////////

#define gf16mat_rowmat_mul_impl             gf16mat_rowmat_mul_ref
//////////////////////////////from blas_matrix.c ///////////////////////////
/// @brief matrix multiplication:  matC = matA * matB , in GF(16)
///
/// @param[out]  matC         - the output row-major matrix C
/// @param[in]   matA         - a row-major matrix A.
/// @param[in]   height_A     - the number of row vectors in the matrix A.
/// @param[in]   width_A_byte  - the size of row vectors of A in bytes.
/// @param[in]   matB            - a row-major matrix B.
/// @param[in]   width_B_byte  - the size of row vectors of B in bytes.
static void gf16mat_rowmat_mul(uint8_t *matC, const uint8_t *matA, unsigned height_A, unsigned width_A_byte, const uint8_t *matB, unsigned width_B_byte)
{
    gf16mat_rowmat_mul_impl( matC, matA, height_A, width_A_byte, matB, width_B_byte);
}

#define gf16mat_colmat_mul_impl             gf16mat_colmat_mul_ref
/// @brief (column-major) matrix multiplication:  matC = matA * matB , in GF(16)
///
/// @param[out]  matC           - the output column-major matrix C
/// @param[in]   mat_a          - a column-major matrix A.
/// @param[in]   a_veclen_byte  - the vector length (height) of the matrix A in byte.
/// @param[in]   a_n_vec        - the number of vectors (width) in the matrix A.
/// @param[in]   mat_b          - a column-major matrix B.
/// @param[in]   b_n_vec        - the number of vectors in the matrix B.
///
static void gf16mat_colmat_mul(uint8_t *mat_c, const uint8_t *mat_a, unsigned a_veclen_byte, unsigned a_n_vec, const uint8_t *mat_b, unsigned b_n_vec)
{
    gf16mat_colmat_mul_impl( mat_c , mat_a , a_veclen_byte , a_n_vec , mat_b , b_n_vec );
}
///////////////////////////////////////////////////////////////////////////////



#endif /* MIRITH_DSS_BLAS_AVX2_H */
