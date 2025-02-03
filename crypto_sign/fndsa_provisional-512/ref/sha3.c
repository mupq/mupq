/*
 * SHAKE implementation.
 */

#include "inner.h"

/* Internal alias for the Keccak-f function. */
#define process_block   fndsa_sha3_process_block

/* Process the provided state.
     A   pointer to state
     r   number of data words (r = rate/8)
   The implementation may use a different internal representation for the
   non-data words, provided that the all-zero initial state is still all-zero
   in that internal representation. */
#if FNDSA_ASM_CORTEXM4
void process_block(uint64_t *A, unsigned r);
#else
/* Round constants. */
static const uint64_t RC[] = {
	0x0000000000000001, 0x0000000000008082,
	0x800000000000808A, 0x8000000080008000,
	0x000000000000808B, 0x0000000080000001,
	0x8000000080008081, 0x8000000000008009,
	0x000000000000008A, 0x0000000000000088,
	0x0000000080008009, 0x000000008000000A,
	0x000000008000808B, 0x800000000000008B,
	0x8000000000008089, 0x8000000000008003,
	0x8000000000008002, 0x8000000000000080,
	0x000000000000800A, 0x800000008000000A,
	0x8000000080008081, 0x8000000000008080,
	0x0000000080000001, 0x8000000080008008
};
void
process_block(uint64_t *A, unsigned r)
{
	uint64_t t0, t1, t2, t3, t4;
	uint64_t tt0, tt1, tt2, tt3;
	uint64_t t, kt;
	uint64_t c0, c1, c2, c3, c4, bnn;
	int j;

	(void)r;

	/* Invert some words (alternate internal representation, which
	   saves some operations). */
	A[ 1] = ~A[ 1];
	A[ 2] = ~A[ 2];
	A[ 8] = ~A[ 8];
	A[12] = ~A[12];
	A[17] = ~A[17];
	A[20] = ~A[20];

	/* Compute the 24 rounds. This loop is partially unrolled (each
	   iteration computes two rounds). */
	for (j = 0; j < 24; j += 2) {
		tt0 = A[ 1] ^ A[ 6];
		tt1 = A[11] ^ A[16];
		tt0 ^= A[21] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 4] ^ A[ 9];
		tt3 = A[14] ^ A[19];
		tt0 ^= A[24];
		tt2 ^= tt3;
		t0 = tt0 ^ tt2;

		tt0 = A[ 2] ^ A[ 7];
		tt1 = A[12] ^ A[17];
		tt0 ^= A[22] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 0] ^ A[ 5];
		tt3 = A[10] ^ A[15];
		tt0 ^= A[20];
		tt2 ^= tt3;
		t1 = tt0 ^ tt2;

		tt0 = A[ 3] ^ A[ 8];
		tt1 = A[13] ^ A[18];
		tt0 ^= A[23] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 1] ^ A[ 6];
		tt3 = A[11] ^ A[16];
		tt0 ^= A[21];
		tt2 ^= tt3;
		t2 = tt0 ^ tt2;

		tt0 = A[ 4] ^ A[ 9];
		tt1 = A[14] ^ A[19];
		tt0 ^= A[24] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 2] ^ A[ 7];
		tt3 = A[12] ^ A[17];
		tt0 ^= A[22];
		tt2 ^= tt3;
		t3 = tt0 ^ tt2;

		tt0 = A[ 0] ^ A[ 5];
		tt1 = A[10] ^ A[15];
		tt0 ^= A[20] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 3] ^ A[ 8];
		tt3 = A[13] ^ A[18];
		tt0 ^= A[23];
		tt2 ^= tt3;
		t4 = tt0 ^ tt2;

		A[ 0] = A[ 0] ^ t0;
		A[ 5] = A[ 5] ^ t0;
		A[10] = A[10] ^ t0;
		A[15] = A[15] ^ t0;
		A[20] = A[20] ^ t0;
		A[ 1] = A[ 1] ^ t1;
		A[ 6] = A[ 6] ^ t1;
		A[11] = A[11] ^ t1;
		A[16] = A[16] ^ t1;
		A[21] = A[21] ^ t1;
		A[ 2] = A[ 2] ^ t2;
		A[ 7] = A[ 7] ^ t2;
		A[12] = A[12] ^ t2;
		A[17] = A[17] ^ t2;
		A[22] = A[22] ^ t2;
		A[ 3] = A[ 3] ^ t3;
		A[ 8] = A[ 8] ^ t3;
		A[13] = A[13] ^ t3;
		A[18] = A[18] ^ t3;
		A[23] = A[23] ^ t3;
		A[ 4] = A[ 4] ^ t4;
		A[ 9] = A[ 9] ^ t4;
		A[14] = A[14] ^ t4;
		A[19] = A[19] ^ t4;
		A[24] = A[24] ^ t4;
		A[ 5] = (A[ 5] << 36) | (A[ 5] >> (64 - 36));
		A[10] = (A[10] <<  3) | (A[10] >> (64 -  3));
		A[15] = (A[15] << 41) | (A[15] >> (64 - 41));
		A[20] = (A[20] << 18) | (A[20] >> (64 - 18));
		A[ 1] = (A[ 1] <<  1) | (A[ 1] >> (64 -  1));
		A[ 6] = (A[ 6] << 44) | (A[ 6] >> (64 - 44));
		A[11] = (A[11] << 10) | (A[11] >> (64 - 10));
		A[16] = (A[16] << 45) | (A[16] >> (64 - 45));
		A[21] = (A[21] <<  2) | (A[21] >> (64 - 2));
		A[ 2] = (A[ 2] << 62) | (A[ 2] >> (64 - 62));
		A[ 7] = (A[ 7] <<  6) | (A[ 7] >> (64 -  6));
		A[12] = (A[12] << 43) | (A[12] >> (64 - 43));
		A[17] = (A[17] << 15) | (A[17] >> (64 - 15));
		A[22] = (A[22] << 61) | (A[22] >> (64 - 61));
		A[ 3] = (A[ 3] << 28) | (A[ 3] >> (64 - 28));
		A[ 8] = (A[ 8] << 55) | (A[ 8] >> (64 - 55));
		A[13] = (A[13] << 25) | (A[13] >> (64 - 25));
		A[18] = (A[18] << 21) | (A[18] >> (64 - 21));
		A[23] = (A[23] << 56) | (A[23] >> (64 - 56));
		A[ 4] = (A[ 4] << 27) | (A[ 4] >> (64 - 27));
		A[ 9] = (A[ 9] << 20) | (A[ 9] >> (64 - 20));
		A[14] = (A[14] << 39) | (A[14] >> (64 - 39));
		A[19] = (A[19] <<  8) | (A[19] >> (64 -  8));
		A[24] = (A[24] << 14) | (A[24] >> (64 - 14));

		bnn = ~A[12];
		kt = A[ 6] | A[12];
		c0 = A[ 0] ^ kt;
		kt = bnn | A[18];
		c1 = A[ 6] ^ kt;
		kt = A[18] & A[24];
		c2 = A[12] ^ kt;
		kt = A[24] | A[ 0];
		c3 = A[18] ^ kt;
		kt = A[ 0] & A[ 6];
		c4 = A[24] ^ kt;
		A[ 0] = c0;
		A[ 6] = c1;
		A[12] = c2;
		A[18] = c3;
		A[24] = c4;
		bnn = ~A[22];
		kt = A[ 9] | A[10];
		c0 = A[ 3] ^ kt;
		kt = A[10] & A[16];
		c1 = A[ 9] ^ kt;
		kt = A[16] | bnn;
		c2 = A[10] ^ kt;
		kt = A[22] | A[ 3];
		c3 = A[16] ^ kt;
		kt = A[ 3] & A[ 9];
		c4 = A[22] ^ kt;
		A[ 3] = c0;
		A[ 9] = c1;
		A[10] = c2;
		A[16] = c3;
		A[22] = c4;
		bnn = ~A[19];
		kt = A[ 7] | A[13];
		c0 = A[ 1] ^ kt;
		kt = A[13] & A[19];
		c1 = A[ 7] ^ kt;
		kt = bnn & A[20];
		c2 = A[13] ^ kt;
		kt = A[20] | A[ 1];
		c3 = bnn ^ kt;
		kt = A[ 1] & A[ 7];
		c4 = A[20] ^ kt;
		A[ 1] = c0;
		A[ 7] = c1;
		A[13] = c2;
		A[19] = c3;
		A[20] = c4;
		bnn = ~A[17];
		kt = A[ 5] & A[11];
		c0 = A[ 4] ^ kt;
		kt = A[11] | A[17];
		c1 = A[ 5] ^ kt;
		kt = bnn | A[23];
		c2 = A[11] ^ kt;
		kt = A[23] & A[ 4];
		c3 = bnn ^ kt;
		kt = A[ 4] | A[ 5];
		c4 = A[23] ^ kt;
		A[ 4] = c0;
		A[ 5] = c1;
		A[11] = c2;
		A[17] = c3;
		A[23] = c4;
		bnn = ~A[ 8];
		kt = bnn & A[14];
		c0 = A[ 2] ^ kt;
		kt = A[14] | A[15];
		c1 = bnn ^ kt;
		kt = A[15] & A[21];
		c2 = A[14] ^ kt;
		kt = A[21] | A[ 2];
		c3 = A[15] ^ kt;
		kt = A[ 2] & A[ 8];
		c4 = A[21] ^ kt;
		A[ 2] = c0;
		A[ 8] = c1;
		A[14] = c2;
		A[15] = c3;
		A[21] = c4;
		A[ 0] = A[ 0] ^ RC[j + 0];

		tt0 = A[ 6] ^ A[ 9];
		tt1 = A[ 7] ^ A[ 5];
		tt0 ^= A[ 8] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[24] ^ A[22];
		tt3 = A[20] ^ A[23];
		tt0 ^= A[21];
		tt2 ^= tt3;
		t0 = tt0 ^ tt2;

		tt0 = A[12] ^ A[10];
		tt1 = A[13] ^ A[11];
		tt0 ^= A[14] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 0] ^ A[ 3];
		tt3 = A[ 1] ^ A[ 4];
		tt0 ^= A[ 2];
		tt2 ^= tt3;
		t1 = tt0 ^ tt2;

		tt0 = A[18] ^ A[16];
		tt1 = A[19] ^ A[17];
		tt0 ^= A[15] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[ 6] ^ A[ 9];
		tt3 = A[ 7] ^ A[ 5];
		tt0 ^= A[ 8];
		tt2 ^= tt3;
		t2 = tt0 ^ tt2;

		tt0 = A[24] ^ A[22];
		tt1 = A[20] ^ A[23];
		tt0 ^= A[21] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[12] ^ A[10];
		tt3 = A[13] ^ A[11];
		tt0 ^= A[14];
		tt2 ^= tt3;
		t3 = tt0 ^ tt2;

		tt0 = A[ 0] ^ A[ 3];
		tt1 = A[ 1] ^ A[ 4];
		tt0 ^= A[ 2] ^ tt1;
		tt0 = (tt0 << 1) | (tt0 >> 63);
		tt2 = A[18] ^ A[16];
		tt3 = A[19] ^ A[17];
		tt0 ^= A[15];
		tt2 ^= tt3;
		t4 = tt0 ^ tt2;

		A[ 0] = A[ 0] ^ t0;
		A[ 3] = A[ 3] ^ t0;
		A[ 1] = A[ 1] ^ t0;
		A[ 4] = A[ 4] ^ t0;
		A[ 2] = A[ 2] ^ t0;
		A[ 6] = A[ 6] ^ t1;
		A[ 9] = A[ 9] ^ t1;
		A[ 7] = A[ 7] ^ t1;
		A[ 5] = A[ 5] ^ t1;
		A[ 8] = A[ 8] ^ t1;
		A[12] = A[12] ^ t2;
		A[10] = A[10] ^ t2;
		A[13] = A[13] ^ t2;
		A[11] = A[11] ^ t2;
		A[14] = A[14] ^ t2;
		A[18] = A[18] ^ t3;
		A[16] = A[16] ^ t3;
		A[19] = A[19] ^ t3;
		A[17] = A[17] ^ t3;
		A[15] = A[15] ^ t3;
		A[24] = A[24] ^ t4;
		A[22] = A[22] ^ t4;
		A[20] = A[20] ^ t4;
		A[23] = A[23] ^ t4;
		A[21] = A[21] ^ t4;
		A[ 3] = (A[ 3] << 36) | (A[ 3] >> (64 - 36));
		A[ 1] = (A[ 1] <<  3) | (A[ 1] >> (64 -  3));
		A[ 4] = (A[ 4] << 41) | (A[ 4] >> (64 - 41));
		A[ 2] = (A[ 2] << 18) | (A[ 2] >> (64 - 18));
		A[ 6] = (A[ 6] <<  1) | (A[ 6] >> (64 -  1));
		A[ 9] = (A[ 9] << 44) | (A[ 9] >> (64 - 44));
		A[ 7] = (A[ 7] << 10) | (A[ 7] >> (64 - 10));
		A[ 5] = (A[ 5] << 45) | (A[ 5] >> (64 - 45));
		A[ 8] = (A[ 8] <<  2) | (A[ 8] >> (64 - 2));
		A[12] = (A[12] << 62) | (A[12] >> (64 - 62));
		A[10] = (A[10] <<  6) | (A[10] >> (64 -  6));
		A[13] = (A[13] << 43) | (A[13] >> (64 - 43));
		A[11] = (A[11] << 15) | (A[11] >> (64 - 15));
		A[14] = (A[14] << 61) | (A[14] >> (64 - 61));
		A[18] = (A[18] << 28) | (A[18] >> (64 - 28));
		A[16] = (A[16] << 55) | (A[16] >> (64 - 55));
		A[19] = (A[19] << 25) | (A[19] >> (64 - 25));
		A[17] = (A[17] << 21) | (A[17] >> (64 - 21));
		A[15] = (A[15] << 56) | (A[15] >> (64 - 56));
		A[24] = (A[24] << 27) | (A[24] >> (64 - 27));
		A[22] = (A[22] << 20) | (A[22] >> (64 - 20));
		A[20] = (A[20] << 39) | (A[20] >> (64 - 39));
		A[23] = (A[23] <<  8) | (A[23] >> (64 -  8));
		A[21] = (A[21] << 14) | (A[21] >> (64 - 14));

		bnn = ~A[13];
		kt = A[ 9] | A[13];
		c0 = A[ 0] ^ kt;
		kt = bnn | A[17];
		c1 = A[ 9] ^ kt;
		kt = A[17] & A[21];
		c2 = A[13] ^ kt;
		kt = A[21] | A[ 0];
		c3 = A[17] ^ kt;
		kt = A[ 0] & A[ 9];
		c4 = A[21] ^ kt;
		A[ 0] = c0;
		A[ 9] = c1;
		A[13] = c2;
		A[17] = c3;
		A[21] = c4;
		bnn = ~A[14];
		kt = A[22] | A[ 1];
		c0 = A[18] ^ kt;
		kt = A[ 1] & A[ 5];
		c1 = A[22] ^ kt;
		kt = A[ 5] | bnn;
		c2 = A[ 1] ^ kt;
		kt = A[14] | A[18];
		c3 = A[ 5] ^ kt;
		kt = A[18] & A[22];
		c4 = A[14] ^ kt;
		A[18] = c0;
		A[22] = c1;
		A[ 1] = c2;
		A[ 5] = c3;
		A[14] = c4;
		bnn = ~A[23];
		kt = A[10] | A[19];
		c0 = A[ 6] ^ kt;
		kt = A[19] & A[23];
		c1 = A[10] ^ kt;
		kt = bnn & A[ 2];
		c2 = A[19] ^ kt;
		kt = A[ 2] | A[ 6];
		c3 = bnn ^ kt;
		kt = A[ 6] & A[10];
		c4 = A[ 2] ^ kt;
		A[ 6] = c0;
		A[10] = c1;
		A[19] = c2;
		A[23] = c3;
		A[ 2] = c4;
		bnn = ~A[11];
		kt = A[ 3] & A[ 7];
		c0 = A[24] ^ kt;
		kt = A[ 7] | A[11];
		c1 = A[ 3] ^ kt;
		kt = bnn | A[15];
		c2 = A[ 7] ^ kt;
		kt = A[15] & A[24];
		c3 = bnn ^ kt;
		kt = A[24] | A[ 3];
		c4 = A[15] ^ kt;
		A[24] = c0;
		A[ 3] = c1;
		A[ 7] = c2;
		A[11] = c3;
		A[15] = c4;
		bnn = ~A[16];
		kt = bnn & A[20];
		c0 = A[12] ^ kt;
		kt = A[20] | A[ 4];
		c1 = bnn ^ kt;
		kt = A[ 4] & A[ 8];
		c2 = A[20] ^ kt;
		kt = A[ 8] | A[12];
		c3 = A[ 4] ^ kt;
		kt = A[12] & A[16];
		c4 = A[ 8] ^ kt;
		A[12] = c0;
		A[16] = c1;
		A[20] = c2;
		A[ 4] = c3;
		A[ 8] = c4;
		A[ 0] = A[ 0] ^ RC[j + 1];

		t = A[ 5];
		A[ 5] = A[18];
		A[18] = A[11];
		A[11] = A[10];
		A[10] = A[ 6];
		A[ 6] = A[22];
		A[22] = A[20];
		A[20] = A[12];
		A[12] = A[19];
		A[19] = A[15];
		A[15] = A[24];
		A[24] = A[ 8];
		A[ 8] = t;
		t = A[ 1];
		A[ 1] = A[ 9];
		A[ 9] = A[14];
		A[14] = A[ 2];
		A[ 2] = A[13];
		A[13] = A[23];
		A[23] = A[ 4];
		A[ 4] = A[21];
		A[21] = A[16];
		A[16] = A[ 3];
		A[ 3] = A[17];
		A[17] = A[ 7];
		A[ 7] = t;
	}

	/* Invert some words back to normal representation. */
	A[ 1] = ~A[ 1];
	A[ 2] = ~A[ 2];
	A[ 8] = ~A[ 8];
	A[12] = ~A[12];
	A[17] = ~A[17];
	A[20] = ~A[20];
}
#endif

#if FNDSA_AVX2
/* Four SHAKE256 instances in parallel. The provided array contains the
   four states, which are interleaved (this is not the same layout as
   in the plain implementation). */
TARGET_AVX2
static void
process_block_x4(uint64_t *A)
{
	__m256i ya[25];

	for (int i = 0; i < 25; i ++) {
		ya[i] = _mm256_loadu_si256((const __m256i *)A + i);
	}

	/*
	 * Compute the 24 rounds. This loop is partially unrolled (each
	 * iteration computes two rounds).
	 */
	for (int j = 0; j < 24; j += 2) {
		__m256i yt0, yt1, yt2, yt3, yt4;

#define yy_rotl(yv, nn)   _mm256_or_si256( \
	_mm256_slli_epi64(yv, nn), _mm256_srli_epi64(yv, 64 - (nn)))
#define yy_andnotL(a, b)   _mm256_andnot_si256(a, b)
#define yy_xor(a, b)       _mm256_xor_si256(a, b)

#define yCOMB1(yd, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9)   do { \
		__m256i ytt0, ytt1, ytt2, ytt3; \
		ytt0 = yy_xor(ya[i0], ya[i1]); \
		ytt1 = yy_xor(ya[i2], ya[i3]); \
		ytt0 = yy_xor(ytt0, yy_xor(ya[i4], ytt1)); \
		ytt0 = yy_rotl(ytt0, 1); \
		ytt2 = yy_xor(ya[i5], ya[i6]); \
		ytt3 = yy_xor(ya[i7], ya[i8]); \
		ytt0 = yy_xor(ytt0, ya[i9]); \
		ytt2 = yy_xor(ytt2, ytt3); \
		yd = yy_xor(ytt0, ytt2); \
	} while (0)

#define yCOMB2(i0, i1, i2, i3, i4, op0, op1, op2, op3, op4)   do { \
		__m256i yc0, yc1, yc2, yc3, yc4, ykt; \
		ykt = yy_andnotL(ya[i1], ya[i2]); \
		yc0 = yy_xor(ykt, ya[i0]); \
		ykt = yy_andnotL(ya[i2], ya[i3]); \
		yc1 = yy_xor(ykt, ya[i1]); \
		ykt = yy_andnotL(ya[i3], ya[i4]); \
		yc2 = yy_xor(ykt, ya[i2]); \
		ykt = yy_andnotL(ya[i4], ya[i0]); \
		yc3 = yy_xor(ykt, ya[i3]); \
		ykt = yy_andnotL(ya[i0], ya[i1]); \
		yc4 = yy_xor(ykt, ya[i4]); \
		ya[i0] = yc0; \
		ya[i1] = yc1; \
		ya[i2] = yc2; \
		ya[i3] = yc3; \
		ya[i4] = yc4; \
	} while (0)

		/* Round j */

		yCOMB1(yt0, 1, 6, 11, 16, 21, 4, 9, 14, 19, 24);
		yCOMB1(yt1, 2, 7, 12, 17, 22, 0, 5, 10, 15, 20);
		yCOMB1(yt2, 3, 8, 13, 18, 23, 1, 6, 11, 16, 21);
		yCOMB1(yt3, 4, 9, 14, 19, 24, 2, 7, 12, 17, 22);
		yCOMB1(yt4, 0, 5, 10, 15, 20, 3, 8, 13, 18, 23);

		ya[ 0] = yy_xor(ya[ 0], yt0);
		ya[ 5] = yy_xor(ya[ 5], yt0);
		ya[10] = yy_xor(ya[10], yt0);
		ya[15] = yy_xor(ya[15], yt0);
		ya[20] = yy_xor(ya[20], yt0);
		ya[ 1] = yy_xor(ya[ 1], yt1);
		ya[ 6] = yy_xor(ya[ 6], yt1);
		ya[11] = yy_xor(ya[11], yt1);
		ya[16] = yy_xor(ya[16], yt1);
		ya[21] = yy_xor(ya[21], yt1);
		ya[ 2] = yy_xor(ya[ 2], yt2);
		ya[ 7] = yy_xor(ya[ 7], yt2);
		ya[12] = yy_xor(ya[12], yt2);
		ya[17] = yy_xor(ya[17], yt2);
		ya[22] = yy_xor(ya[22], yt2);
		ya[ 3] = yy_xor(ya[ 3], yt3);
		ya[ 8] = yy_xor(ya[ 8], yt3);
		ya[13] = yy_xor(ya[13], yt3);
		ya[18] = yy_xor(ya[18], yt3);
		ya[23] = yy_xor(ya[23], yt3);
		ya[ 4] = yy_xor(ya[ 4], yt4);
		ya[ 9] = yy_xor(ya[ 9], yt4);
		ya[14] = yy_xor(ya[14], yt4);
		ya[19] = yy_xor(ya[19], yt4);
		ya[24] = yy_xor(ya[24], yt4);
		ya[ 5] = yy_rotl(ya[ 5], 36);
		ya[10] = yy_rotl(ya[10],  3);
		ya[15] = yy_rotl(ya[15], 41);
		ya[20] = yy_rotl(ya[20], 18);
		ya[ 1] = yy_rotl(ya[ 1],  1);
		ya[ 6] = yy_rotl(ya[ 6], 44);
		ya[11] = yy_rotl(ya[11], 10);
		ya[16] = yy_rotl(ya[16], 45);
		ya[21] = yy_rotl(ya[21],  2);
		ya[ 2] = yy_rotl(ya[ 2], 62);
		ya[ 7] = yy_rotl(ya[ 7],  6);
		ya[12] = yy_rotl(ya[12], 43);
		ya[17] = yy_rotl(ya[17], 15);
		ya[22] = yy_rotl(ya[22], 61);
		ya[ 3] = yy_rotl(ya[ 3], 28);
		ya[ 8] = yy_rotl(ya[ 8], 55);
		ya[13] = yy_rotl(ya[13], 25);
		ya[18] = yy_rotl(ya[18], 21);
		ya[23] = yy_rotl(ya[23], 56);
		ya[ 4] = yy_rotl(ya[ 4], 27);
		ya[ 9] = yy_rotl(ya[ 9], 20);
		ya[14] = yy_rotl(ya[14], 39);
		ya[19] = yy_rotl(ya[19],  8);
		ya[24] = yy_rotl(ya[24], 14);

		yCOMB2(0, 6, 12, 18, 24, or, ornotL, and, or, and);
		yCOMB2(3, 9, 10, 16, 22, or, and, ornotR, or, and);
		yCOMB2(1, 7, 13, 19, 20, or, andnotR, and, or, and);
		yCOMB2(4, 5, 11, 17, 23, and, ornotR, or, and, or);
		yCOMB2(2, 8, 14, 15, 21, and, or, and, or, andnotR);

		ya[0] = yy_xor(ya[0], _mm256_set1_epi64x(RC[j + 0]));

		/* Round j + 1 */

		yCOMB1(yt0, 6, 9, 7, 5, 8, 24, 22, 20, 23, 21);
		yCOMB1(yt1, 12, 10, 13, 11, 14, 0, 3, 1, 4, 2);
		yCOMB1(yt2, 18, 16, 19, 17, 15, 6, 9, 7, 5, 8);
		yCOMB1(yt3, 24, 22, 20, 23, 21, 12, 10, 13, 11, 14);
		yCOMB1(yt4, 0, 3, 1, 4, 2, 18, 16, 19, 17, 15);

		ya[ 0] = yy_xor(ya[ 0], yt0);
		ya[ 3] = yy_xor(ya[ 3], yt0);
		ya[ 1] = yy_xor(ya[ 1], yt0);
		ya[ 4] = yy_xor(ya[ 4], yt0);
		ya[ 2] = yy_xor(ya[ 2], yt0);
		ya[ 6] = yy_xor(ya[ 6], yt1);
		ya[ 9] = yy_xor(ya[ 9], yt1);
		ya[ 7] = yy_xor(ya[ 7], yt1);
		ya[ 5] = yy_xor(ya[ 5], yt1);
		ya[ 8] = yy_xor(ya[ 8], yt1);
		ya[12] = yy_xor(ya[12], yt2);
		ya[10] = yy_xor(ya[10], yt2);
		ya[13] = yy_xor(ya[13], yt2);
		ya[11] = yy_xor(ya[11], yt2);
		ya[14] = yy_xor(ya[14], yt2);
		ya[18] = yy_xor(ya[18], yt3);
		ya[16] = yy_xor(ya[16], yt3);
		ya[19] = yy_xor(ya[19], yt3);
		ya[17] = yy_xor(ya[17], yt3);
		ya[15] = yy_xor(ya[15], yt3);
		ya[24] = yy_xor(ya[24], yt4);
		ya[22] = yy_xor(ya[22], yt4);
		ya[20] = yy_xor(ya[20], yt4);
		ya[23] = yy_xor(ya[23], yt4);
		ya[21] = yy_xor(ya[21], yt4);
		ya[ 3] = yy_rotl(ya[ 3], 36);
		ya[ 1] = yy_rotl(ya[ 1],  3);
		ya[ 4] = yy_rotl(ya[ 4], 41);
		ya[ 2] = yy_rotl(ya[ 2], 18);
		ya[ 6] = yy_rotl(ya[ 6],  1);
		ya[ 9] = yy_rotl(ya[ 9], 44);
		ya[ 7] = yy_rotl(ya[ 7], 10);
		ya[ 5] = yy_rotl(ya[ 5], 45);
		ya[ 8] = yy_rotl(ya[ 8],  2);
		ya[12] = yy_rotl(ya[12], 62);
		ya[10] = yy_rotl(ya[10],  6);
		ya[13] = yy_rotl(ya[13], 43);
		ya[11] = yy_rotl(ya[11], 15);
		ya[14] = yy_rotl(ya[14], 61);
		ya[18] = yy_rotl(ya[18], 28);
		ya[16] = yy_rotl(ya[16], 55);
		ya[19] = yy_rotl(ya[19], 25);
		ya[17] = yy_rotl(ya[17], 21);
		ya[15] = yy_rotl(ya[15], 56);
		ya[24] = yy_rotl(ya[24], 27);
		ya[22] = yy_rotl(ya[22], 20);
		ya[20] = yy_rotl(ya[20], 39);
		ya[23] = yy_rotl(ya[23],  8);
		ya[21] = yy_rotl(ya[21], 14);

		yCOMB2(0, 9, 13, 17, 21, or, ornotL, and, or, and);
		yCOMB2(18, 22, 1, 5, 14, or, and, ornotR, or, and);
		yCOMB2(6, 10, 19, 23, 2, or, andnotR, and, or, and);
		yCOMB2(24, 3, 7, 11, 15, and, ornotR, or, and, or);
		yCOMB2(12, 16, 20, 4, 8, and, or, and, or, andnotR);

		ya[0] = yy_xor(ya[0], _mm256_set1_epi64x(RC[j + 1]));

		/* Apply combined permutation for next round */

		__m256i yt = ya[ 5];
		ya[ 5] = ya[18];
		ya[18] = ya[11];
		ya[11] = ya[10];
		ya[10] = ya[ 6];
		ya[ 6] = ya[22];
		ya[22] = ya[20];
		ya[20] = ya[12];
		ya[12] = ya[19];
		ya[19] = ya[15];
		ya[15] = ya[24];
		ya[24] = ya[ 8];
		ya[ 8] = yt;
		yt = ya[ 1];
		ya[ 1] = ya[ 9];
		ya[ 9] = ya[14];
		ya[14] = ya[ 2];
		ya[ 2] = ya[13];
		ya[13] = ya[23];
		ya[23] = ya[ 4];
		ya[ 4] = ya[21];
		ya[21] = ya[16];
		ya[16] = ya[ 3];
                ya[ 3] = ya[17];
                ya[17] = ya[ 7];
                ya[ 7] = yt;

#undef yy_rotl
#undef yy_andnotL
#undef yy_xor
#undef yCOMB1
#undef yCOMB2
	}

	/*
	 * Write back state words.
	 */
	for (int i = 0; i < 25; i ++) {
		_mm256_storeu_si256((__m256i *)A + i, ya[i]);
	}
}
#endif

#if FNDSA_SSE2
/* This is a variant of the AVX2 process_block_x4() function, but with
   only SSE2 opcodes, it runs only two blocks in parallel. It still uses
   the same memory layout as the AVX2 code (the four states are
   interleaved). */
TARGET_SSE2
static void
process_block_x2(uint64_t *A)
{
	__m128i xa[25];

	for (int i = 0; i < 25; i ++) {
		xa[i] = _mm_loadu_si128((const __m128i *)A + (i << 1));
	}

	/*
	 * Compute the 24 rounds. This loop is partially unrolled (each
	 * iteration computes two rounds).
	 */
	for (int j = 0; j < 24; j += 2) {
		__m128i xt0, xt1, xt2, xt3, xt4;

#define xx_rotl(xv, nn)   _mm_or_si128( \
	_mm_slli_epi64(xv, nn), _mm_srli_epi64(xv, 64 - (nn)))
#define xx_andnotL(a, b)   _mm_andnot_si128(a, b)
#define xx_xor(a, b)       _mm_xor_si128(a, b)

#define xCOMB1(xd, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9)   do { \
		__m128i xtt0, xtt1, xtt2, xtt3; \
		xtt0 = xx_xor(xa[i0], xa[i1]); \
		xtt1 = xx_xor(xa[i2], xa[i3]); \
		xtt0 = xx_xor(xtt0, xx_xor(xa[i4], xtt1)); \
		xtt0 = xx_rotl(xtt0, 1); \
		xtt2 = xx_xor(xa[i5], xa[i6]); \
		xtt3 = xx_xor(xa[i7], xa[i8]); \
		xtt0 = xx_xor(xtt0, xa[i9]); \
		xtt2 = xx_xor(xtt2, xtt3); \
		xd = xx_xor(xtt0, xtt2); \
	} while (0)

#define xCOMB2(i0, i1, i2, i3, i4, op0, op1, op2, op3, op4)   do { \
		__m128i xc0, xc1, xc2, xc3, xc4, xkt; \
		xkt = xx_andnotL(xa[i1], xa[i2]); \
		xc0 = xx_xor(xkt, xa[i0]); \
		xkt = xx_andnotL(xa[i2], xa[i3]); \
		xc1 = xx_xor(xkt, xa[i1]); \
		xkt = xx_andnotL(xa[i3], xa[i4]); \
		xc2 = xx_xor(xkt, xa[i2]); \
		xkt = xx_andnotL(xa[i4], xa[i0]); \
		xc3 = xx_xor(xkt, xa[i3]); \
		xkt = xx_andnotL(xa[i0], xa[i1]); \
		xc4 = xx_xor(xkt, xa[i4]); \
		xa[i0] = xc0; \
		xa[i1] = xc1; \
		xa[i2] = xc2; \
		xa[i3] = xc3; \
		xa[i4] = xc4; \
	} while (0)

		/* Round j */

		xCOMB1(xt0, 1, 6, 11, 16, 21, 4, 9, 14, 19, 24);
		xCOMB1(xt1, 2, 7, 12, 17, 22, 0, 5, 10, 15, 20);
		xCOMB1(xt2, 3, 8, 13, 18, 23, 1, 6, 11, 16, 21);
		xCOMB1(xt3, 4, 9, 14, 19, 24, 2, 7, 12, 17, 22);
		xCOMB1(xt4, 0, 5, 10, 15, 20, 3, 8, 13, 18, 23);

		xa[ 0] = xx_xor(xa[ 0], xt0);
		xa[ 5] = xx_xor(xa[ 5], xt0);
		xa[10] = xx_xor(xa[10], xt0);
		xa[15] = xx_xor(xa[15], xt0);
		xa[20] = xx_xor(xa[20], xt0);
		xa[ 1] = xx_xor(xa[ 1], xt1);
		xa[ 6] = xx_xor(xa[ 6], xt1);
		xa[11] = xx_xor(xa[11], xt1);
		xa[16] = xx_xor(xa[16], xt1);
		xa[21] = xx_xor(xa[21], xt1);
		xa[ 2] = xx_xor(xa[ 2], xt2);
		xa[ 7] = xx_xor(xa[ 7], xt2);
		xa[12] = xx_xor(xa[12], xt2);
		xa[17] = xx_xor(xa[17], xt2);
		xa[22] = xx_xor(xa[22], xt2);
		xa[ 3] = xx_xor(xa[ 3], xt3);
		xa[ 8] = xx_xor(xa[ 8], xt3);
		xa[13] = xx_xor(xa[13], xt3);
		xa[18] = xx_xor(xa[18], xt3);
		xa[23] = xx_xor(xa[23], xt3);
		xa[ 4] = xx_xor(xa[ 4], xt4);
		xa[ 9] = xx_xor(xa[ 9], xt4);
		xa[14] = xx_xor(xa[14], xt4);
		xa[19] = xx_xor(xa[19], xt4);
		xa[24] = xx_xor(xa[24], xt4);
		xa[ 5] = xx_rotl(xa[ 5], 36);
		xa[10] = xx_rotl(xa[10],  3);
		xa[15] = xx_rotl(xa[15], 41);
		xa[20] = xx_rotl(xa[20], 18);
		xa[ 1] = xx_rotl(xa[ 1],  1);
		xa[ 6] = xx_rotl(xa[ 6], 44);
		xa[11] = xx_rotl(xa[11], 10);
		xa[16] = xx_rotl(xa[16], 45);
		xa[21] = xx_rotl(xa[21],  2);
		xa[ 2] = xx_rotl(xa[ 2], 62);
		xa[ 7] = xx_rotl(xa[ 7],  6);
		xa[12] = xx_rotl(xa[12], 43);
		xa[17] = xx_rotl(xa[17], 15);
		xa[22] = xx_rotl(xa[22], 61);
		xa[ 3] = xx_rotl(xa[ 3], 28);
		xa[ 8] = xx_rotl(xa[ 8], 55);
		xa[13] = xx_rotl(xa[13], 25);
		xa[18] = xx_rotl(xa[18], 21);
		xa[23] = xx_rotl(xa[23], 56);
		xa[ 4] = xx_rotl(xa[ 4], 27);
		xa[ 9] = xx_rotl(xa[ 9], 20);
		xa[14] = xx_rotl(xa[14], 39);
		xa[19] = xx_rotl(xa[19],  8);
		xa[24] = xx_rotl(xa[24], 14);

		xCOMB2(0, 6, 12, 18, 24, or, ornotL, and, or, and);
		xCOMB2(3, 9, 10, 16, 22, or, and, ornotR, or, and);
		xCOMB2(1, 7, 13, 19, 20, or, andnotR, and, or, and);
		xCOMB2(4, 5, 11, 17, 23, and, ornotR, or, and, or);
		xCOMB2(2, 8, 14, 15, 21, and, or, and, or, andnotR);

		xa[0] = xx_xor(xa[0], _mm_set1_epi64x(RC[j + 0]));

		/* Round j + 1 */

		xCOMB1(xt0, 6, 9, 7, 5, 8, 24, 22, 20, 23, 21);
		xCOMB1(xt1, 12, 10, 13, 11, 14, 0, 3, 1, 4, 2);
		xCOMB1(xt2, 18, 16, 19, 17, 15, 6, 9, 7, 5, 8);
		xCOMB1(xt3, 24, 22, 20, 23, 21, 12, 10, 13, 11, 14);
		xCOMB1(xt4, 0, 3, 1, 4, 2, 18, 16, 19, 17, 15);

		xa[ 0] = xx_xor(xa[ 0], xt0);
		xa[ 3] = xx_xor(xa[ 3], xt0);
		xa[ 1] = xx_xor(xa[ 1], xt0);
		xa[ 4] = xx_xor(xa[ 4], xt0);
		xa[ 2] = xx_xor(xa[ 2], xt0);
		xa[ 6] = xx_xor(xa[ 6], xt1);
		xa[ 9] = xx_xor(xa[ 9], xt1);
		xa[ 7] = xx_xor(xa[ 7], xt1);
		xa[ 5] = xx_xor(xa[ 5], xt1);
		xa[ 8] = xx_xor(xa[ 8], xt1);
		xa[12] = xx_xor(xa[12], xt2);
		xa[10] = xx_xor(xa[10], xt2);
		xa[13] = xx_xor(xa[13], xt2);
		xa[11] = xx_xor(xa[11], xt2);
		xa[14] = xx_xor(xa[14], xt2);
		xa[18] = xx_xor(xa[18], xt3);
		xa[16] = xx_xor(xa[16], xt3);
		xa[19] = xx_xor(xa[19], xt3);
		xa[17] = xx_xor(xa[17], xt3);
		xa[15] = xx_xor(xa[15], xt3);
		xa[24] = xx_xor(xa[24], xt4);
		xa[22] = xx_xor(xa[22], xt4);
		xa[20] = xx_xor(xa[20], xt4);
		xa[23] = xx_xor(xa[23], xt4);
		xa[21] = xx_xor(xa[21], xt4);
		xa[ 3] = xx_rotl(xa[ 3], 36);
		xa[ 1] = xx_rotl(xa[ 1],  3);
		xa[ 4] = xx_rotl(xa[ 4], 41);
		xa[ 2] = xx_rotl(xa[ 2], 18);
		xa[ 6] = xx_rotl(xa[ 6],  1);
		xa[ 9] = xx_rotl(xa[ 9], 44);
		xa[ 7] = xx_rotl(xa[ 7], 10);
		xa[ 5] = xx_rotl(xa[ 5], 45);
		xa[ 8] = xx_rotl(xa[ 8],  2);
		xa[12] = xx_rotl(xa[12], 62);
		xa[10] = xx_rotl(xa[10],  6);
		xa[13] = xx_rotl(xa[13], 43);
		xa[11] = xx_rotl(xa[11], 15);
		xa[14] = xx_rotl(xa[14], 61);
		xa[18] = xx_rotl(xa[18], 28);
		xa[16] = xx_rotl(xa[16], 55);
		xa[19] = xx_rotl(xa[19], 25);
		xa[17] = xx_rotl(xa[17], 21);
		xa[15] = xx_rotl(xa[15], 56);
		xa[24] = xx_rotl(xa[24], 27);
		xa[22] = xx_rotl(xa[22], 20);
		xa[20] = xx_rotl(xa[20], 39);
		xa[23] = xx_rotl(xa[23],  8);
		xa[21] = xx_rotl(xa[21], 14);

		xCOMB2(0, 9, 13, 17, 21, or, ornotL, and, or, and);
		xCOMB2(18, 22, 1, 5, 14, or, and, ornotR, or, and);
		xCOMB2(6, 10, 19, 23, 2, or, andnotR, and, or, and);
		xCOMB2(24, 3, 7, 11, 15, and, ornotR, or, and, or);
		xCOMB2(12, 16, 20, 4, 8, and, or, and, or, andnotR);

		xa[0] = xx_xor(xa[0], _mm_set1_epi64x(RC[j + 1]));

		/* Apply combined permutation for next round */

		__m128i xt = xa[ 5];
		xa[ 5] = xa[18];
		xa[18] = xa[11];
		xa[11] = xa[10];
		xa[10] = xa[ 6];
		xa[ 6] = xa[22];
		xa[22] = xa[20];
		xa[20] = xa[12];
		xa[12] = xa[19];
		xa[19] = xa[15];
		xa[15] = xa[24];
		xa[24] = xa[ 8];
		xa[ 8] = xt;
		xt = xa[ 1];
		xa[ 1] = xa[ 9];
		xa[ 9] = xa[14];
		xa[14] = xa[ 2];
		xa[ 2] = xa[13];
		xa[13] = xa[23];
		xa[23] = xa[ 4];
		xa[ 4] = xa[21];
		xa[21] = xa[16];
		xa[16] = xa[ 3];
                xa[ 3] = xa[17];
                xa[17] = xa[ 7];
                xa[ 7] = xt;

#undef xx_rotl
#undef xx_or
#undef xx_ornotL
#undef xx_ornotR
#undef xx_and
#undef xx_andnotL
#undef xx_andnotR
#undef xx_xor
#undef xCOMB1
#undef xCOMB2
	}

	/*
	 * Write back state words.
	 */
	for (int i = 0; i < 25; i ++) {
		_mm_storeu_si128((__m128i *)A + (i << 1), xa[i]);
	}
}
#elif FNDSA_NEON_SHA3
/* Similar to the SSE2 implementation of process_block_x2(), but using
   NEON opcodes (for aarch64). */
TARGET_NEON
static void
process_block_x2(uint64_t *A)
{
	uint64x2_t xa[25];

	for (int i = 0; i < 25; i ++) {
		xa[i] = vld1q_u64(A + (i << 2));
	}

	/* Invert some words (alternate internal representation, which
	   saves some operations). */
	uint64x2_t xones = vdupq_n_u64(~(uint64_t)0);
	xa[ 1] = veorq_u64(xa[ 1], xones);
	xa[ 2] = veorq_u64(xa[ 2], xones);
	xa[ 8] = veorq_u64(xa[ 8], xones);
	xa[12] = veorq_u64(xa[12], xones);
	xa[17] = veorq_u64(xa[17], xones);
	xa[20] = veorq_u64(xa[20], xones);

	/*
	 * Compute the 24 rounds. This loop is partially unrolled (each
	 * iteration computes two rounds).
	 */
	for (int j = 0; j < 24; j += 2) {
		uint64x2_t xt0, xt1, xt2, xt3, xt4;

#define xx_rotl(xv, nn)    vsliq_n_u64(vshrq_n_u64(xv, 64 - (nn)), xv, nn)
#define xx_or(a, b)        vorrq_u64(a, b)
#define xx_ornotL(a, b)    vornq_u64(b, a)
#define xx_ornotR(a, b)    vornq_u64(a, b)
#define xx_and(a, b)       vandq_u64(a, b)
#define xx_andnotL(a, b)   vandq_u64(veorq_u64(a, xones), b)
#define xx_andnotR(a, b)   vandq_u64(a, veorq_u64(b, xones))
#define xx_xor(a, b)       veorq_u64(a, b)

#define xCOMB1(xd, i0, i1, i2, i3, i4, i5, i6, i7, i8, i9)   do { \
		uint64x2_t xtt0, xtt1, xtt2, xtt3; \
		xtt0 = xx_xor(xa[i0], xa[i1]); \
		xtt1 = xx_xor(xa[i2], xa[i3]); \
		xtt0 = xx_xor(xtt0, xx_xor(xa[i4], xtt1)); \
		xtt0 = xx_rotl(xtt0, 1); \
		xtt2 = xx_xor(xa[i5], xa[i6]); \
		xtt3 = xx_xor(xa[i7], xa[i8]); \
		xtt0 = xx_xor(xtt0, xa[i9]); \
		xtt2 = xx_xor(xtt2, xtt3); \
		xd = xx_xor(xtt0, xtt2); \
	} while (0)

#define xCOMB2(i0, i1, i2, i3, i4, op0, op1, op2, op3, op4)   do { \
		uint64x2_t xc0, xc1, xc2, xc3, xc4, xkt; \
		xkt = xx_ ## op0(xa[i1], xa[i2]); \
		xc0 = xx_xor(xkt, xa[i0]); \
		xkt = xx_ ## op1(xa[i2], xa[i3]); \
		xc1 = xx_xor(xkt, xa[i1]); \
		xkt = xx_ ## op2(xa[i3], xa[i4]); \
		xc2 = xx_xor(xkt, xa[i2]); \
		xkt = xx_ ## op3(xa[i4], xa[i0]); \
		xc3 = xx_xor(xkt, xa[i3]); \
		xkt = xx_ ## op4(xa[i0], xa[i1]); \
		xc4 = xx_xor(xkt, xa[i4]); \
		xa[i0] = xc0; \
		xa[i1] = xc1; \
		xa[i2] = xc2; \
		xa[i3] = xc3; \
		xa[i4] = xc4; \
	} while (0)

		/* Round j */

		xCOMB1(xt0, 1, 6, 11, 16, 21, 4, 9, 14, 19, 24);
		xCOMB1(xt1, 2, 7, 12, 17, 22, 0, 5, 10, 15, 20);
		xCOMB1(xt2, 3, 8, 13, 18, 23, 1, 6, 11, 16, 21);
		xCOMB1(xt3, 4, 9, 14, 19, 24, 2, 7, 12, 17, 22);
		xCOMB1(xt4, 0, 5, 10, 15, 20, 3, 8, 13, 18, 23);

		xa[ 0] = xx_xor(xa[ 0], xt0);
		xa[ 5] = xx_xor(xa[ 5], xt0);
		xa[10] = xx_xor(xa[10], xt0);
		xa[15] = xx_xor(xa[15], xt0);
		xa[20] = xx_xor(xa[20], xt0);
		xa[ 1] = xx_xor(xa[ 1], xt1);
		xa[ 6] = xx_xor(xa[ 6], xt1);
		xa[11] = xx_xor(xa[11], xt1);
		xa[16] = xx_xor(xa[16], xt1);
		xa[21] = xx_xor(xa[21], xt1);
		xa[ 2] = xx_xor(xa[ 2], xt2);
		xa[ 7] = xx_xor(xa[ 7], xt2);
		xa[12] = xx_xor(xa[12], xt2);
		xa[17] = xx_xor(xa[17], xt2);
		xa[22] = xx_xor(xa[22], xt2);
		xa[ 3] = xx_xor(xa[ 3], xt3);
		xa[ 8] = xx_xor(xa[ 8], xt3);
		xa[13] = xx_xor(xa[13], xt3);
		xa[18] = xx_xor(xa[18], xt3);
		xa[23] = xx_xor(xa[23], xt3);
		xa[ 4] = xx_xor(xa[ 4], xt4);
		xa[ 9] = xx_xor(xa[ 9], xt4);
		xa[14] = xx_xor(xa[14], xt4);
		xa[19] = xx_xor(xa[19], xt4);
		xa[24] = xx_xor(xa[24], xt4);
		xa[ 5] = xx_rotl(xa[ 5], 36);
		xa[10] = xx_rotl(xa[10],  3);
		xa[15] = xx_rotl(xa[15], 41);
		xa[20] = xx_rotl(xa[20], 18);
		xa[ 1] = xx_rotl(xa[ 1],  1);
		xa[ 6] = xx_rotl(xa[ 6], 44);
		xa[11] = xx_rotl(xa[11], 10);
		xa[16] = xx_rotl(xa[16], 45);
		xa[21] = xx_rotl(xa[21],  2);
		xa[ 2] = xx_rotl(xa[ 2], 62);
		xa[ 7] = xx_rotl(xa[ 7],  6);
		xa[12] = xx_rotl(xa[12], 43);
		xa[17] = xx_rotl(xa[17], 15);
		xa[22] = xx_rotl(xa[22], 61);
		xa[ 3] = xx_rotl(xa[ 3], 28);
		xa[ 8] = xx_rotl(xa[ 8], 55);
		xa[13] = xx_rotl(xa[13], 25);
		xa[18] = xx_rotl(xa[18], 21);
		xa[23] = xx_rotl(xa[23], 56);
		xa[ 4] = xx_rotl(xa[ 4], 27);
		xa[ 9] = xx_rotl(xa[ 9], 20);
		xa[14] = xx_rotl(xa[14], 39);
		xa[19] = xx_rotl(xa[19],  8);
		xa[24] = xx_rotl(xa[24], 14);

		xCOMB2(0, 6, 12, 18, 24, or, ornotL, and, or, and);
		xCOMB2(3, 9, 10, 16, 22, or, and, ornotR, or, and);
		xa[19] = xx_xor(xa[19], xones);
		xCOMB2(1, 7, 13, 19, 20, or, andnotR, and, or, and);
		xa[17] = xx_xor(xa[17], xones);
		xCOMB2(4, 5, 11, 17, 23, and, ornotR, or, and, or);
		xa[8] = xx_xor(xa[8], xones);
		xCOMB2(2, 8, 14, 15, 21, and, or, and, or, andnotR);

		xa[0] = xx_xor(xa[0], vdupq_n_u64(RC[j + 0]));

		/* Round j + 1 */

		xCOMB1(xt0, 6, 9, 7, 5, 8, 24, 22, 20, 23, 21);
		xCOMB1(xt1, 12, 10, 13, 11, 14, 0, 3, 1, 4, 2);
		xCOMB1(xt2, 18, 16, 19, 17, 15, 6, 9, 7, 5, 8);
		xCOMB1(xt3, 24, 22, 20, 23, 21, 12, 10, 13, 11, 14);
		xCOMB1(xt4, 0, 3, 1, 4, 2, 18, 16, 19, 17, 15);

		xa[ 0] = xx_xor(xa[ 0], xt0);
		xa[ 3] = xx_xor(xa[ 3], xt0);
		xa[ 1] = xx_xor(xa[ 1], xt0);
		xa[ 4] = xx_xor(xa[ 4], xt0);
		xa[ 2] = xx_xor(xa[ 2], xt0);
		xa[ 6] = xx_xor(xa[ 6], xt1);
		xa[ 9] = xx_xor(xa[ 9], xt1);
		xa[ 7] = xx_xor(xa[ 7], xt1);
		xa[ 5] = xx_xor(xa[ 5], xt1);
		xa[ 8] = xx_xor(xa[ 8], xt1);
		xa[12] = xx_xor(xa[12], xt2);
		xa[10] = xx_xor(xa[10], xt2);
		xa[13] = xx_xor(xa[13], xt2);
		xa[11] = xx_xor(xa[11], xt2);
		xa[14] = xx_xor(xa[14], xt2);
		xa[18] = xx_xor(xa[18], xt3);
		xa[16] = xx_xor(xa[16], xt3);
		xa[19] = xx_xor(xa[19], xt3);
		xa[17] = xx_xor(xa[17], xt3);
		xa[15] = xx_xor(xa[15], xt3);
		xa[24] = xx_xor(xa[24], xt4);
		xa[22] = xx_xor(xa[22], xt4);
		xa[20] = xx_xor(xa[20], xt4);
		xa[23] = xx_xor(xa[23], xt4);
		xa[21] = xx_xor(xa[21], xt4);
		xa[ 3] = xx_rotl(xa[ 3], 36);
		xa[ 1] = xx_rotl(xa[ 1],  3);
		xa[ 4] = xx_rotl(xa[ 4], 41);
		xa[ 2] = xx_rotl(xa[ 2], 18);
		xa[ 6] = xx_rotl(xa[ 6],  1);
		xa[ 9] = xx_rotl(xa[ 9], 44);
		xa[ 7] = xx_rotl(xa[ 7], 10);
		xa[ 5] = xx_rotl(xa[ 5], 45);
		xa[ 8] = xx_rotl(xa[ 8],  2);
		xa[12] = xx_rotl(xa[12], 62);
		xa[10] = xx_rotl(xa[10],  6);
		xa[13] = xx_rotl(xa[13], 43);
		xa[11] = xx_rotl(xa[11], 15);
		xa[14] = xx_rotl(xa[14], 61);
		xa[18] = xx_rotl(xa[18], 28);
		xa[16] = xx_rotl(xa[16], 55);
		xa[19] = xx_rotl(xa[19], 25);
		xa[17] = xx_rotl(xa[17], 21);
		xa[15] = xx_rotl(xa[15], 56);
		xa[24] = xx_rotl(xa[24], 27);
		xa[22] = xx_rotl(xa[22], 20);
		xa[20] = xx_rotl(xa[20], 39);
		xa[23] = xx_rotl(xa[23],  8);
		xa[21] = xx_rotl(xa[21], 14);

		xCOMB2(0, 9, 13, 17, 21, or, ornotL, and, or, and);
		xCOMB2(18, 22, 1, 5, 14, or, and, ornotR, or, and);
		xa[23] = xx_xor(xa[23], xones);
		xCOMB2(6, 10, 19, 23, 2, or, andnotR, and, or, and);
		xa[11] = xx_xor(xa[11], xones);
		xCOMB2(24, 3, 7, 11, 15, and, ornotR, or, and, or);
		xa[16] = xx_xor(xa[16], xones);
		xCOMB2(12, 16, 20, 4, 8, and, or, and, or, andnotR);

		xa[0] = xx_xor(xa[0], vdupq_n_u64(RC[j + 1]));

		/* Apply combined permutation for next round */

		uint64x2_t xt = xa[ 5];
		xa[ 5] = xa[18];
		xa[18] = xa[11];
		xa[11] = xa[10];
		xa[10] = xa[ 6];
		xa[ 6] = xa[22];
		xa[22] = xa[20];
		xa[20] = xa[12];
		xa[12] = xa[19];
		xa[19] = xa[15];
		xa[15] = xa[24];
		xa[24] = xa[ 8];
		xa[ 8] = xt;
		xt = xa[ 1];
		xa[ 1] = xa[ 9];
		xa[ 9] = xa[14];
		xa[14] = xa[ 2];
		xa[ 2] = xa[13];
		xa[13] = xa[23];
		xa[23] = xa[ 4];
		xa[ 4] = xa[21];
		xa[21] = xa[16];
		xa[16] = xa[ 3];
                xa[ 3] = xa[17];
                xa[17] = xa[ 7];
                xa[ 7] = xt;

#undef xx_rotl
#undef xx_or
#undef xx_ornotL
#undef xx_ornotR
#undef xx_and
#undef xx_andnotL
#undef xx_andnotR
#undef xx_xor
#undef xCOMB1
#undef xCOMB2
	}

	/* Invert some words back to normal representation. */
	xa[ 1] = veorq_u64(xa[ 1], xones);
	xa[ 2] = veorq_u64(xa[ 2], xones);
	xa[ 8] = veorq_u64(xa[ 8], xones);
	xa[12] = veorq_u64(xa[12], xones);
	xa[17] = veorq_u64(xa[17], xones);
	xa[20] = veorq_u64(xa[20], xones);

	/*
	 * Write back state words.
	 */
	for (int i = 0; i < 25; i ++) {
		vst1q_u64(A + (i << 2), xa[i]);
	}
}
#endif

/* see inner.h */
void
shake_init(shake_context *sc, unsigned size)
{
	sc->rate = 200 - (size_t)(size >> 2);
	sc->dptr = 0;
	memset(sc->A, 0, sizeof sc->A);
}

#if FNDSA_ASM_CORTEXM4
/* Inject the specified chunk by XORing the source bytes into the
   destination. Length MUST NOT be zero. */
void fndsa_sha3_inject_chunk(void *dst, const void *src, size_t len);
#endif

/* see inner.h */
void
shake_inject(shake_context *sc, const void *in, size_t len)
{
	size_t dptr, rate;
	const uint8_t *buf;

	dptr = sc->dptr;
	rate = sc->rate;
	buf = in;
	while (len > 0) {
		size_t clen = rate - dptr;
		if (clen > len) {
			clen = len;
		}
#if FNDSA_ASM_CORTEXM4
		fndsa_sha3_inject_chunk((uint8_t *)sc->A + dptr, buf, clen);
#else
		for (size_t u = 0; u < clen; u ++) {
			size_t v;

			v = u + dptr;
			sc->A[v >> 3] ^= (uint64_t)buf[u] << ((v & 7) << 3);
		}
#endif
		dptr += clen;
		buf += clen;
		len -= clen;
		if (dptr == rate) {
			process_block(sc->A, rate >> 3);
			dptr = 0;
		}
	}
	sc->dptr = (unsigned)dptr;
}

/* see inner.h */
void
shake_flip(shake_context *sc)
{
	/* We apply padding and pre-XOR the value into the state. We
	   set dptr to the end of the buffer, so that first call to
	   shake_extract() will process the block. */
	unsigned v;

	v = (unsigned)sc->dptr;
	sc->A[v >> 3] ^= (uint64_t)0x1F << ((v & 7) << 3);
	v = (unsigned)sc->rate - 1;
	sc->A[v >> 3] ^= (uint64_t)0x80 << ((v & 7) << 3);
	sc->dptr = sc->rate;
}

/* see inner.h */
void
shake_extract(shake_context *sc, void *out, size_t len)
{
	size_t dptr, rate;
	uint8_t *buf;

	dptr = sc->dptr;
	rate = sc->rate;
	buf = out;
	while (len > 0) {
		size_t clen;

		if (dptr == rate) {
			process_block(sc->A, rate >> 3);
			dptr = 0;
		}
		clen = rate - dptr;
		if (clen > len) {
			clen = len;
		}
		len -= clen;
#if FNDSA_ASM_CORTEXM4
		/* On the ARM Cortex M4, we can assume little-endian
		   encoding. */
		if (buf != NULL) {
			memcpy(buf, (uint8_t *)sc->A + dptr, clen);
			buf += clen;
		}
		dptr += clen;
#else
		if (buf != NULL) {
			while (clen -- > 0) {
				*buf ++ = (uint8_t)(sc->A[dptr >> 3]
					>> ((dptr & 7) << 3));
				dptr ++;
			}
		} else {
			dptr += clen;
		}
#endif
	}
	sc->dptr = (unsigned)dptr;
}

#if FNDSA_SHAKE256X4
/* Little-endian 64-bit decoding. */
static inline uint64_t
dec64le(const void *src)
{
	const uint8_t *buf = (const uint8_t *)src;
	return (uint64_t)buf[0]
		| ((uint64_t)buf[1] << 8)
		| ((uint64_t)buf[2] << 16)
		| ((uint64_t)buf[3] << 24)
		| ((uint64_t)buf[4] << 32)
		| ((uint64_t)buf[5] << 40)
		| ((uint64_t)buf[6] << 48)
		| ((uint64_t)buf[7] << 56);
}

/* Little-endian 64-bit encoding. */
static inline void
enc64le(void *dst, uint64_t x)
{
	uint8_t *buf = (uint8_t *)dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
	buf[4] = (uint8_t)(x >> 32);
	buf[5] = (uint8_t)(x >> 40);
	buf[6] = (uint8_t)(x >> 48);
	buf[7] = (uint8_t)(x >> 56);
}

/* see inner.h */
void
shake256x4_init(shake256x4_context *sc, const void *seed, size_t seed_len)
{
	memset(sc->state, 0, sizeof sc->state);
	const uint8_t *sbuf = (const uint8_t *)seed;
	size_t rlen = seed_len & ~(size_t)7;
	size_t elen = seed_len - rlen;
	size_t k = rlen >> 3;
	sc->ptr = sizeof sc->buf;

#if FNDSA_AVX2
	if ((sc->use_avx2 = has_avx2()) != 0) {
		/* The AVX2 implementation interleaves the four SHAKE256
		   states. */
		for (size_t i = 0; i < rlen; i += 8) {
			uint64_t x = dec64le(sbuf + i);
			sc->state[((i >> 3) << 2) + 0] = x;
			sc->state[((i >> 3) << 2) + 1] = x;
			sc->state[((i >> 3) << 2) + 2] = x;
			sc->state[((i >> 3) << 2) + 3] = x;
		}
		uint64_t x = 0;
		for (size_t j = 0; j < elen; j ++) {
			x |= (uint64_t)sbuf[rlen + j] << (j << 3);
		}
		if (elen < 7) {
			x |= (uint64_t)0x1F << ((elen + 1) << 3);
			sc->state[(k << 2) + 0] = x;
			sc->state[(k << 2) + 1] =
				x | ((uint64_t)0x01 << (elen << 3));
			sc->state[(k << 2) + 2] =
				x | ((uint64_t)0x02 << (elen << 3));
			sc->state[(k << 2) + 3] =
				x | ((uint64_t)0x03 << (elen << 3));
		} else {
			sc->state[(k << 2) + 0] = x;
			sc->state[(k << 2) + 1] = x | ((uint64_t)0x01 << 56);
			sc->state[(k << 2) + 2] = x | ((uint64_t)0x02 << 56);
			sc->state[(k << 2) + 3] = x | ((uint64_t)0x03 << 56);
			sc->state[(k << 2) + 4] = 0x1F;
			sc->state[(k << 2) + 5] = 0x1F;
			sc->state[(k << 2) + 6] = 0x1F;
			sc->state[(k << 2) + 7] = 0x1F;
		}
		sc->state[64] ^= (uint64_t)0x80 << 56;
		sc->state[65] ^= (uint64_t)0x80 << 56;
		sc->state[66] ^= (uint64_t)0x80 << 56;
		sc->state[67] ^= (uint64_t)0x80 << 56;
		return;
	}
#endif

#if FNDSA_SSE2 || FNDSA_NEON_SHA3
	/* The SSE2 and NEON implementations interleave the four
	   SHAKE256 states. */
	for (size_t i = 0; i < rlen; i += 8) {
		uint64_t x = dec64le(sbuf + i);
		sc->state[((i >> 3) << 2) + 0] = x;
		sc->state[((i >> 3) << 2) + 1] = x;
		sc->state[((i >> 3) << 2) + 2] = x;
		sc->state[((i >> 3) << 2) + 3] = x;
	}
	uint64_t x = 0;
	for (size_t j = 0; j < elen; j ++) {
		x |= (uint64_t)sbuf[rlen + j] << (j << 3);
	}
	if (elen < 7) {
		x |= (uint64_t)0x1F << ((elen + 1) << 3);
		sc->state[(k << 2) + 0] = x;
		sc->state[(k << 2) + 1] =
			x | ((uint64_t)0x01 << (elen << 3));
		sc->state[(k << 2) + 2] =
			x | ((uint64_t)0x02 << (elen << 3));
		sc->state[(k << 2) + 3] =
			x | ((uint64_t)0x03 << (elen << 3));
	} else {
		sc->state[(k << 2) + 0] = x;
		sc->state[(k << 2) + 1] = x | ((uint64_t)0x01 << 56);
		sc->state[(k << 2) + 2] = x | ((uint64_t)0x02 << 56);
		sc->state[(k << 2) + 3] = x | ((uint64_t)0x03 << 56);
		sc->state[(k << 2) + 4] = 0x1F;
		sc->state[(k << 2) + 5] = 0x1F;
		sc->state[(k << 2) + 6] = 0x1F;
		sc->state[(k << 2) + 7] = 0x1F;
	}
	sc->state[64] ^= (uint64_t)0x80 << 56;
	sc->state[65] ^= (uint64_t)0x80 << 56;
	sc->state[66] ^= (uint64_t)0x80 << 56;
	sc->state[67] ^= (uint64_t)0x80 << 56;

#else
	/* In the plain implementation, the four SHAKE256 states are
	   successive (i.e. not interleaved). */
	for (size_t i = 0; i < rlen; i += 8) {
		sc->state[i >> 3] = dec64le(sbuf + i);
	}
	memcpy(sc->state + 25, sc->state, rlen);
	memcpy(sc->state + 50, sc->state, rlen);
	memcpy(sc->state + 75, sc->state, rlen);
	uint64_t x = 0;
	for (size_t j = 0; j < elen; j ++) {
		x |= (uint64_t)sbuf[rlen + j] << (j << 3);
	}
	if (elen < 7) {
		/* We have room for both the instance ID and the padding
		   byte. */
		x |= (uint64_t)0x1F << ((elen + 1) << 3);
		sc->state[k] = x;
		sc->state[25 + k] = x | ((uint64_t)0x01 << (elen << 3));
		sc->state[50 + k] = x | ((uint64_t)0x02 << (elen << 3));
		sc->state[75 + k] = x | ((uint64_t)0x03 << (elen << 3));
	} else {
		/* ID byte and padding byte fall on different words. */
		sc->state[k] = x;
		sc->state[25 + k] = x | ((uint64_t)0x01 << 56);
		sc->state[50 + k] = x | ((uint64_t)0x02 << 56);
		sc->state[75 + k] = x | ((uint64_t)0x03 << 56);
		sc->state[k +  1] = 0x1F;
		sc->state[k + 26] = 0x1F;
		sc->state[k + 51] = 0x1F;
		sc->state[k + 76] = 0x1F;
	}
	sc->state[16] ^= (uint64_t)0x80 << 56;
	sc->state[41] ^= (uint64_t)0x80 << 56;
	sc->state[66] ^= (uint64_t)0x80 << 56;
	sc->state[91] ^= (uint64_t)0x80 << 56;
#endif
}

/* see inner.h */
void
shake256x4_refill(shake256x4_context *sc)
{
	sc->ptr = 0;
#if FNDSA_AVX2
	if (sc->use_avx2) {
		process_block_x4(sc->state);
		memcpy(sc->buf, sc->state, sizeof sc->buf);
		return;
	}
#endif
#if FNDSA_SSE2 || FNDSA_NEON_SHA3
	process_block_x2(sc->state);
	process_block_x2(sc->state + 2);
	memcpy(sc->buf, sc->state, sizeof sc->buf);
#else
	process_block(sc->state, 17);
	process_block(sc->state + 25, 17);
	process_block(sc->state + 50, 17);
	process_block(sc->state + 75, 17);

	/* Interleave the outputs into the buffer. */
	for (int i = 0; i < 17; i ++) {
		enc64le(&sc->buf[(i << 5) +  0], sc->state[ 0 + i]);
		enc64le(&sc->buf[(i << 5) +  8], sc->state[25 + i]);
		enc64le(&sc->buf[(i << 5) + 16], sc->state[50 + i]);
		enc64le(&sc->buf[(i << 5) + 24], sc->state[75 + i]);
	}
#endif
}
#endif

/* SHA-3 is mostly the same as SHAKE, except for the padding, and the
   fact that the output size is fixed. */

/* see inner.h */
void
sha3_init(sha3_context *sc, unsigned size)
{
	shake_init(sc, size);
}

/* see inner.h */
void
sha3_update(sha3_context *sc, const void *in, size_t len)
{
	shake_inject(sc, in, len);
}

/* see inner.h */
void
sha3_close(sha3_context *sc, void *out)
{
	unsigned v;

	/* Padding starts with '01' instead of '1111'. */
	v = (unsigned)sc->dptr;
	sc->A[v >> 3] ^= (uint64_t)0x06 << ((v & 7) << 3);
	v = (unsigned)sc->rate - 1;
	sc->A[v >> 3] ^= (uint64_t)0x80 << ((v & 7) << 3);

	/* Process the padding block, then write the output. */
	process_block(sc->A, sc->rate >> 3);
	size_t len = (200 - sc->rate) >> 1;
	uint8_t *buf = out;
	for (size_t i = 0; i < len; i ++) {
		buf[i] = (uint8_t)(sc->A[i >> 3] >> ((i & 7) << 3));
	}

	/* Reinitialize the context. */
	sha3_init(sc, len << 3);
}
