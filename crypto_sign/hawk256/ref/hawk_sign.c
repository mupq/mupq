#include "hawk_inner.h"

/*
 * We use computations with polynomials modulo X^n+1 and modulo 18433.
 */
#define Q   18433
#include "modq.h"

/*
 * Binary polynomials (GF(2)[X]); inputs a and b have size N bits;
 * output d has size N or 2*N bits, depending on the operation:
 *
 *    bp_xor_N():     d <- a XOR b
 *    bp_mul_N():     d <- a*b
 *    bp_muladd_N():  d <- d + a*b
 *    bp_mulmod_N():  d <- a*b mod X^N+1
 *
 * Multiplications use Karatsuba:
 * If:
 *    a = a0 + a1*(X^(N/2))
 *    b = b0 + b1*(X^(N/2))
 * Then:
 *    d0 <- a0*b0
 *    d1 <- a1*b1
 *    d2 <- (a0 + a1)*(b0 + b1) - d0 - d1
 *    d <- d0 + d1*X^(N/2) + d2*X^N
 * Addition and subtraction are both XOR for binary polynomials; there is
 * no carry, so no need for an extra bit.
 *
 * bp_muladd_N() may use up to 4*N bits in the provided 'tmp' buffer.
 * bp_modadd_N() may use up to 3*N bits in the provided 'tmp' buffer.
 */

/*
 * Multiplication of two binary polynomials ("carryless multiplication")
 * of degree less than 32, result has degree less than 64.
 */
static uint64_t
bp_mul_32(uint32_t x, uint32_t y)
{
	/*
	 * Classic technique (rediscovered many times): integer
	 * multiplications with "holes" to absorb unwanted carries.
	 */
	uint32_t x0 = x & 0x11111111;
	uint32_t x1 = x & 0x22222222;
	uint32_t x2 = x & 0x44444444;
	uint32_t x3 = x & 0x88888888;

	uint32_t y0 = y & 0x11111111;
	uint32_t y1 = y & 0x22222222;
	uint32_t y2 = y & 0x44444444;
	uint32_t y3 = y & 0x88888888;

#define M(a, b)   ((uint64_t)(a) * (uint64_t)(b))

	uint64_t z0 = M(x0, y0) ^ M(x1, y3) ^ M(x2, y2) ^ M(x3, y1);
	uint64_t z1 = M(x0, y1) ^ M(x1, y0) ^ M(x2, y3) ^ M(x3, y2);
	uint64_t z2 = M(x0, y2) ^ M(x1, y1) ^ M(x2, y0) ^ M(x3, y3);
	uint64_t z3 = M(x0, y3) ^ M(x1, y2) ^ M(x2, y1) ^ M(x3, y0);

#undef M

	z0 &= (uint64_t)0x1111111111111111;
	z1 &= (uint64_t)0x2222222222222222;
	z2 &= (uint64_t)0x4444444444444444;
	z3 &= (uint64_t)0x8888888888888888;
	return z0 | z1 | z2 | z3;
}

static inline void
bp_xor_64(uint8_t *d, const uint8_t *a, const uint8_t *b)
{
#if HAWK_UNALIGNED
	*(uint64_t *)d = *(const uint64_t *)a ^ *(const uint64_t *)b;
#else
	for (size_t u = 0; u < 8; u ++) {
		d[u] = a[u] ^ b[u];
	}
#endif
}

static inline void
bp_xor_128(uint8_t *d, const uint8_t *a, const uint8_t *b)
{
#if HAWK_UNALIGNED
	const uint64_t *aq = (uint64_t *)a;
	const uint64_t *bq = (uint64_t *)b;
	uint64_t *dq = (uint64_t *)d;
	for (size_t u = 0; u < 2; u ++) {
		dq[u] = aq[u] ^ bq[u];
	}
#else
	for (size_t u = 0; u < 16; u ++) {
		d[u] = a[u] ^ b[u];
	}
#endif
}

static inline void
bp_xor_256(uint8_t *d, const uint8_t *a, const uint8_t *b)
{
#if HAWK_UNALIGNED
	const uint64_t *aq = (uint64_t *)a;
	const uint64_t *bq = (uint64_t *)b;
	uint64_t *dq = (uint64_t *)d;
	for (size_t u = 0; u < 4; u ++) {
		dq[u] = aq[u] ^ bq[u];
	}
#else
	for (size_t u = 0; u < 32; u ++) {
		d[u] = a[u] ^ b[u];
	}
#endif
}

static inline void
bp_xor_512(uint8_t *d, const uint8_t *a, const uint8_t *b)
{
#if HAWK_UNALIGNED
	const uint64_t *aq = (uint64_t *)a;
	const uint64_t *bq = (uint64_t *)b;
	uint64_t *dq = (uint64_t *)d;
	for (size_t u = 0; u < 8; u ++) {
		dq[u] = aq[u] ^ bq[u];
	}
#else
	for (size_t u = 0; u < 64; u ++) {
		d[u] = a[u] ^ b[u];
	}
#endif
}

static inline void
bp_xor_1024(uint8_t *d, const uint8_t *a, const uint8_t *b)
{
#if HAWK_UNALIGNED
	const uint64_t *aq = (uint64_t *)a;
	const uint64_t *bq = (uint64_t *)b;
	uint64_t *dq = (uint64_t *)d;
	for (size_t u = 0; u < 16; u ++) {
		dq[u] = aq[u] ^ bq[u];
	}
#else
	for (size_t u = 0; u < 128; u ++) {
		d[u] = a[u] ^ b[u];
	}
#endif
}

static void
bp_muladd_64(uint8_t *restrict d, const uint8_t *restrict a,
	const uint8_t *restrict b, uint8_t *restrict tmp)
{
	(void)tmp;

	uint32_t a0 = dec32le(a);
	uint32_t a1 = dec32le(a + 4);
	uint32_t b0 = dec32le(b);
	uint32_t b1 = dec32le(b + 4);

	uint64_t c0 = bp_mul_32(a0, b0);
	uint64_t c1 = bp_mul_32(a1, b1);
	uint64_t c2 = bp_mul_32(a0 ^ a1, b0 ^ b1) ^ c0 ^ c1;

	enc64le(d, dec64le(d) ^ c0 ^ (c2 << 32));
	enc64le(d + 8, dec64le(d + 8) ^ c1 ^ (c2 >> 32));
}

#define MKBP_MULADD(n, hn)    MKBP_MULADD_(n, hn)
#define MKBP_MULADD_(n, hn) \
static void \
bp_muladd_ ## n(uint8_t *restrict d, const uint8_t *restrict a, \
	const uint8_t *restrict b, uint8_t *restrict tmp) \
{ \
	uint8_t *t1 = tmp; \
	uint8_t *t2 = t1 + (n / 8); \
	uint8_t *t3 = t2 + (n / 8); \
	/* t1 <- (a0 + a1)*(b0 + b1) + d0 + d1 */ \
	bp_xor_ ## hn(t2, a, a + (hn / 8)); \
	bp_xor_ ## hn(t2 + (hn / 8), b, b + (hn / 8)); \
	bp_xor_ ## n(t1, d, d + (n / 8)); \
	bp_muladd_ ## hn(t1, t2, t2 + (hn / 8), t3); \
	/* d0 <- d0 + a0*b0 */ \
	bp_muladd_ ## hn(d, a, b, t2); \
	/* d1 <- d1 + a1*b1 */ \
	bp_muladd_ ## hn(d + (n / 8), a + (hn / 8), b + (hn / 8), t2); \
	/* t1 <- t1 + d0 + d1 = a0*b1 + a1*b0 */ \
	bp_xor_ ## n(t1, t1, d); \
	bp_xor_ ## n(t1, t1, d + (n / 8)); \
	/* d <- d + (x^{n/2})*t1 */ \
	bp_xor_ ## n(d + (hn / 8), d + (hn / 8), t1); \
}

MKBP_MULADD(128, 64)
MKBP_MULADD(256, 128)
MKBP_MULADD(512, 256)

#define MKBP_MULMOD(n, hn)    MKBP_MULMOD_(n, hn)
#define MKBP_MULMOD_(n, hn) \
static void \
bp_mulmod_ ## n(uint8_t *restrict d, const uint8_t *restrict a, \
	const uint8_t *restrict b, uint8_t *restrict tmp) \
{ \
	uint8_t *t1 = tmp; \
	uint8_t *t2 = t1 + (n / 8); \
	/* t1 <- (a0 + a1)*(b0 + b1) */ \
	bp_xor_ ## hn(d, a, a + (hn / 8)); \
	bp_xor_ ## hn(d + (hn / 8), b, b + (hn / 8)); \
	memset(t1, 0, (n / 8)); \
	bp_muladd_ ## hn(t1, d, d + (hn / 8), t2); \
	/* d <- a0*b0 + a1*b1 */ \
	memset(d, 0, (n / 8)); \
	bp_muladd_ ## hn(d, a, b, t2); \
	bp_muladd_ ## hn(d, a + (hn / 8), b + (hn / 8), t2); \
	/* t1 <- t1 + d = a0*b1 + a1*b0 */ \
	bp_xor_ ## n(t1, t1, d); \
	/* d <- d + rotate_{n/2}(t1) */ \
	bp_xor_ ## hn(d, d, t1 + (hn / 8)); \
	bp_xor_ ## hn(d + (hn / 8), d + (hn / 8), t1); \
}

MKBP_MULMOD(256, 128)
MKBP_MULMOD(512, 256)
MKBP_MULMOD(1024, 512)

/*
 * Reduce a polynomial modulo 2 (output f2 is n/8 bytes = n bits).
 */
static void
extract_lowbit(unsigned logn, uint8_t *f2, const int8_t *f)
{
	size_t n = (size_t)1 << logn;
	const uint8_t *fu = (const uint8_t *)f;
	for (size_t u = 0; u < n; u += 8) {
		f2[u >> 3] = (fu[u + 0] & 1)
			| ((fu[u + 1] & 1) << 1)
			| ((fu[u + 2] & 1) << 2)
			| ((fu[u + 3] & 1) << 3)
			| ((fu[u + 4] & 1) << 4)
			| ((fu[u + 5] & 1) << 5)
			| ((fu[u + 6] & 1) << 6)
			| ((fu[u + 7] & 1) << 7);
	}
}

/*
 * (t0, t1) <- B*h  (mod 2)
 * B = [[f F],[g G]]
 * Only the low bit of each coefficient of f, g, F and G is provided;
 * parameters f2, g2, F2 and G2 are arrays of n/8 bytes (n bits) each.
 * Inputs h0, h1 and outputs t0, t1 are also arrays of n/8 bytes.
 *
 * tmp usage: at most n/2 bytes (4*n bits)
 */
static void
basis_m2_mul(unsigned logn, uint8_t *restrict t0, uint8_t *restrict t1,
	const uint8_t *restrict h0, const uint8_t *restrict h1,
	const uint8_t *restrict f2, const uint8_t *restrict g2,
	const uint8_t *restrict F2, const uint8_t *restrict G2,
	uint8_t *tmp)
{
	/*
	 * t0 = f*h0 + F*h1
	 * t1 = g*h0 + G*h1
	 */
	size_t n = (size_t)1 << logn;
	uint8_t *w1 = tmp;
	uint8_t *w2 = w1 + (n >> 3);
	switch (logn) {
	case 8:
		bp_mulmod_256(t0, h0, f2, w2);
		bp_mulmod_256(w1, h1, F2, w2);
		bp_xor_256(t0, t0, w1);
		bp_mulmod_256(t1, h0, g2, w2);
		bp_mulmod_256(w1, h1, G2, w2);
		bp_xor_256(t1, t1, w1);
		break;
	case 9:
		bp_mulmod_512(t0, h0, f2, w2);
		bp_mulmod_512(w1, h1, F2, w2);
		bp_xor_512(t0, t0, w1);
		bp_mulmod_512(t1, h0, g2, w2);
		bp_mulmod_512(w1, h1, G2, w2);
		bp_xor_512(t1, t1, w1);
		break;
	case 10:
		bp_mulmod_1024(t0, h0, f2, w2);
		bp_mulmod_1024(w1, h1, F2, w2);
		bp_xor_1024(t0, t0, w1);
		bp_mulmod_1024(t1, h0, g2, w2);
		bp_mulmod_1024(w1, h1, G2, w2);
		bp_xor_1024(t1, t1, w1);
		break;
	}
}

/*
 * Tables for the Gaussian sampler: we have two distributions over 2*Z,
 * with the same standard deviation 2*sigma. Given:
 *    p(x) = exp(-(x^2) / 2*(2*sigma)^2)
 * Then, for integers k:
 *    D0(2*k) = p(2*k) / \sum_j p(2*j)
 *    D1(2*k) = p(2*k-1) / \sum_j p(2*j-1)
 * D0 is centred on 0, while D1 is centred on 1. Both distributions only
 * return even integers.
 *
 * Let P0(x) = P(|X0| >= x)  (with X0 selected with distribution D0).
 * Let P1(x) = P(|X1| >= x)  (with X1 selected with distribution D1).
 * For integers k >= 0, we define the table T:
 *    T[2*k]   = floor(P0(2*(k+1)) * 2^78)
 *    T[2*k+1] = floor(P1(2*(k+1)+1) * 2^78)
 * Each 78-bit value is split into a high part ("hi") of 15 bits, and
 * a low part ("lo") of 63 bits.
 */

static const uint16_t sig_gauss_hi_Hawk_256[] = {
	0x4D70, 0x268B,
	0x0F80, 0x04FA,
	0x0144, 0x0041,
	0x000A, 0x0001
};
#define SG_MAX_HI_Hawk_256  ((sizeof sig_gauss_hi_Hawk_256) / sizeof(uint16_t))

static const uint64_t sig_gauss_lo_Hawk_256[] = {
	0x71FBD58485D45050, 0x1408A4B181C718B1,
	0x54114F1DC2FA7AC9, 0x614569CC54722DC9,
	0x42F74ADDA0B5AE61, 0x151C5CDCBAFF49A3,
	0x252E2152AB5D758B, 0x23460C30AC398322,
	0x0FDE62196C1718FC, 0x01355A8330C44097,
	0x00127325DDF8CEBA, 0x0000DC8DE401FD12,
	0x000008100822C548, 0x0000003B0FFB28F0,
	0x0000000152A6E9AE, 0x0000000005EFCD99,
	0x000000000014DA4A, 0x0000000000003953,
	0x000000000000007B, 0x0000000000000000
};
#define SG_MAX_LO_Hawk_256  ((sizeof sig_gauss_lo_Hawk_256) / sizeof(uint64_t))

static const uint16_t sig_gauss_hi_Hawk_512[] = {
	0x580B, 0x35F9,
	0x1D34, 0x0DD7,
	0x05B7, 0x020C,
	0x00A2, 0x002B,
	0x000A, 0x0001
};
#define SG_MAX_HI_Hawk_512  ((sizeof sig_gauss_hi_Hawk_512) / sizeof(uint16_t))

static const uint64_t sig_gauss_lo_Hawk_512[] = {
	0x0C27920A04F8F267, 0x3C689D9213449DC9,
	0x1C4FF17C204AA058, 0x7B908C81FCE3524F,
	0x5E63263BE0098FFD, 0x4EBEFD8FF4F07378,
	0x56AEDFB0876A3BD8, 0x4628BC6B23887196,
	0x061E21D588CC61CC, 0x7F769211F07B326F,
	0x2BA568D92EEC18E7, 0x0668F461693DFF8F,
	0x00CF0F8687D3B009, 0x001670DB65964485,
	0x000216A0C344EB45, 0x00002AB6E11C2552,
	0x000002EDF0B98A84, 0x0000002C253C7E81,
	0x000000023AF3B2E7, 0x0000000018C14ABF,
	0x0000000000EBCC6A, 0x000000000007876E,
	0x00000000000034CF, 0x000000000000013D,
	0x0000000000000006, 0x0000000000000000
};
#define SG_MAX_LO_Hawk_512  ((sizeof sig_gauss_lo_Hawk_512) / sizeof(uint64_t))

static const uint16_t sig_gauss_hi_Hawk_1024[] = {
	0x58B0, 0x36FE,
	0x1E3A, 0x0EA0,
	0x0632, 0x024A,
	0x00BC, 0x0034,
	0x000C, 0x0002
};
#define SG_MAX_HI_Hawk_1024 ((sizeof sig_gauss_hi_Hawk_1024) / sizeof(uint16_t))

static const uint64_t sig_gauss_lo_Hawk_1024[] = {
	0x3AAA2EB76504E560, 0x01AE2B17728DF2DE,
	0x70E1C03E49BB683E, 0x6A00B82C69624C93,
	0x55CDA662EF2D1C48, 0x2685DB30348656A4,
	0x31E874B355421BB7, 0x430192770E205503,
	0x57C0676C029895A7, 0x5353BD4091AA96DB,
	0x3D4D67696E51F820, 0x09915A53D8667BEE,
	0x014A1A8A93F20738, 0x0026670030160D5F,
	0x0003DAF47E8DFB21, 0x0000557CD1C5F797,
	0x000006634617B3FF, 0x0000006965E15B13,
	0x00000005DBEFB646, 0x0000000047E9AB38,
	0x0000000002F93038, 0x00000000001B2445,
	0x000000000000D5A7, 0x00000000000005AA,
	0x0000000000000021, 0x0000000000000000
};
#define SG_MAX_LO_Hawk_1024 ((sizeof sig_gauss_lo_Hawk_1024) / sizeof(uint64_t))

/*
 * Generate x with the right Gaussian, for the specified parity bits.
 * x is formally generated with center t/2 and standard deviation sigma_sign
 * (with sigma_sign = 1.010, 1.278 or 1.299, depending on degree); this
 * function generates 2*x.
 *
 * Returned value is the squared norm of x.
 *
 * Internally, a random 40-byte seed is obtained from the RNG, to initialize
 * (along with a counter) four internal SHAKE instances. If sc_extra is not
 * NULL, then it should be a non-flipped SHAKE context which will be used
 * as basis for the four instances; otherwise, empty contexts will be used.
 */
static uint32_t
sig_gauss(unsigned logn,
	void (*rng)(void *ctx, void *dst, size_t len), void *rng_context,
	const shake_context *sc_extra, int8_t *x, const uint8_t *t)
{
	const uint16_t *tab_hi;
	const uint64_t *tab_lo;
	size_t hi_len, lo_len;

	switch (logn) {
	case 8:
		tab_hi = sig_gauss_hi_Hawk_256;
		tab_lo = sig_gauss_lo_Hawk_256;
		hi_len = SG_MAX_HI_Hawk_256;
		lo_len = SG_MAX_LO_Hawk_256;
		break;
	case 9:
		tab_hi = sig_gauss_hi_Hawk_512;
		tab_lo = sig_gauss_lo_Hawk_512;
		hi_len = SG_MAX_HI_Hawk_512;
		lo_len = SG_MAX_LO_Hawk_512;
		break;
	default: /* 10 */
		tab_hi = sig_gauss_hi_Hawk_1024;
		tab_lo = sig_gauss_lo_Hawk_1024;
		hi_len = SG_MAX_HI_Hawk_1024;
		lo_len = SG_MAX_LO_Hawk_1024;
		break;
	}

	/*
	 * We produce 2*n samples. Each sample requires two elements
	 * from the RNG: a low part (64 bits) and a high part (16 bits).
	 * A SHAKE instance produces 64-bit words in groups of five
	 * (little-endian encoding is assumed); the first four words are
	 * the low parts for four samples (lo0, lo1, lo2 and lo3), and the
	 * fifth word contains the high parts for the four samples
	 * (hi0 + (hi1 << 16) + (hi2 << 32) + (hi3 << 48)).
	 *
	 * We formally use four parallel SHAKE instances, instantiated
	 * over seeds obtained from the source RNG. For j = 0 to 3,
	 * instance j produces samples of index 16*i + 4*j + k, with
	 * k = 0, 1, 2 or 3. A small, RAM-constrained implementation
	 * will use the four instances sequentially, reusing the space
	 * for the SHAKE context; a large system with vector instructions
	 * (e.g. AVX2) might run the four SHAKE instances at the same time.
	 */
	size_t n = (size_t)1 << logn;
	uint8_t seed[41];
	rng(rng_context, seed, 40);

	uint32_t sn = 0;
	for (size_t j = 0; j < 4; j ++) {
		shake_context sc;
		if (sc_extra != NULL) {
			sc = *sc_extra;
		} else {
			shake_init(&sc, 256);
		}
		seed[40] = (uint8_t)j;
		shake_inject(&sc, seed, 41);
		shake_flip(&sc);
		for (size_t u = 0; u < (n << 1); u += 16) {
			union {
				uint8_t b[40];
				uint16_t w[20];
				uint64_t q[5];
			} buf;
			shake_extract(&sc, buf.b, 40);
			for (size_t k = 0; k < 4; k ++) {
				size_t v = u + (j << 2) + k;
				uint64_t lo = dec64le(buf.b + (k << 3));
				uint32_t hi = dec16le(buf.b + 32 + (k << 1));

				/* Extract sign bit. */
				uint32_t neg = -(uint32_t)(lo >> 63);
				lo &= 0x7FFFFFFFFFFFFFFF;
				hi &= 0x7FFF;

				/* Use even or odd column depending on
				   parity of t. */
				uint32_t pbit = (t[v >> 3] >> (v & 7)) & 1;
				uint64_t p_odd = -(uint64_t)pbit;
				uint32_t p_oddw = (uint32_t)p_odd;

				uint32_t r = 0;
				for (size_t i = 0; i < hi_len; i += 2) {
					uint64_t tlo0 = tab_lo[i + 0];
					uint64_t tlo1 = tab_lo[i + 1];
					uint64_t tlo = tlo0
						^ (p_odd & (tlo0 ^ tlo1));
					uint32_t cc =
						(uint32_t)((lo - tlo) >> 63);
					uint32_t thi0 = tab_hi[i + 0];
					uint32_t thi1 = tab_hi[i + 1];
					uint32_t thi = thi0
						^ (p_oddw & (thi0 ^ thi1));
					r += (hi - thi - cc) >> 31;
				}
				uint32_t hinz = (hi - 1) >> 31;
				for (size_t i = hi_len; i < lo_len; i += 2) {
					uint64_t tlo0 = tab_lo[i + 0];
					uint64_t tlo1 = tab_lo[i + 1];
					uint64_t tlo = tlo0
						^ (p_odd & (tlo0 ^ tlo1));
					uint32_t cc =
						(uint32_t)((lo - tlo) >> 63);
					r += hinz & cc;
				}

				/* Multiply by 2 and apply parity. */
				r = (r << 1) - p_oddw;

				/* Apply sign bit. */
				r = (r ^ neg) - neg;

				x[v] = (int8_t)*(int32_t *)&r;
				sn += r * r;
			}
		}
	}
	return sn;

}

/*
 * Alternate function for sampling x; the same mechanism is used, but the
 * provided RNG is used directly instead of instantiating four SHAKE
 * instances in parallel.
 *
 * Returned value is the squared norm of x.
 */
static uint32_t
sig_gauss_alt(unsigned logn,
	void (*rng)(void *ctx, void *dst, size_t len), void *rng_context,
	int8_t *x, const uint8_t *t)
{
	const uint16_t *tab_hi;
	const uint64_t *tab_lo;
	size_t hi_len, lo_len;

	switch (logn) {
	case 8:
		tab_hi = sig_gauss_hi_Hawk_256;
		tab_lo = sig_gauss_lo_Hawk_256;
		hi_len = SG_MAX_HI_Hawk_256;
		lo_len = SG_MAX_LO_Hawk_256;
		break;
	case 9:
		tab_hi = sig_gauss_hi_Hawk_512;
		tab_lo = sig_gauss_lo_Hawk_512;
		hi_len = SG_MAX_HI_Hawk_512;
		lo_len = SG_MAX_LO_Hawk_512;
		break;
	default: /* 10 */
		tab_hi = sig_gauss_hi_Hawk_1024;
		tab_lo = sig_gauss_lo_Hawk_1024;
		hi_len = SG_MAX_HI_Hawk_1024;
		lo_len = SG_MAX_LO_Hawk_1024;
		break;
	}

	/*
	 * We produce 2*n samples. Each sample requires two elements
	 * from the RNG: a low part (64 bits) and a high part (16 bits).
	 * A SHAKE instance produces 64-bit words in groups of five
	 * (little-endian encoding is assumed); the first four words are
	 * the low parts for four samples (lo0, lo1, lo2 and lo3), and the
	 * fifth word contains the high parts for the four samples
	 * (hi0 + (hi1 << 16) + (hi2 << 32) + (hi3 << 48)).
	 *
	 * We formally use four parallel SHAKE instances, instantiated
	 * over seeds obtained from the source RNG. For j = 0 to 3,
	 * instance j produces samples of index 16*i + 4*j + k, with
	 * k = 0, 1, 2 or 3. A small, RAM-constrained implementation
	 * will use the four instances sequentially, reusing the space
	 * for the SHAKE context; a large system with vector instructions
	 * (e.g. AVX2) might run the four SHAKE instances at the same time.
	 */
	size_t n = (size_t)1 << logn;

	uint32_t sn = 0;
	for (size_t u = 0; u < (n << 1); u += 16) {
		union {
			uint8_t b[160];
			uint16_t w[80];
			uint64_t q[20];
		} buf;
		rng(rng_context, buf.b, sizeof buf.b);
		for (size_t j = 0; j < 4; j ++) {
			for (size_t k = 0; k < 4; k ++) {
				size_t v = u + (j << 2) + k;
				uint64_t lo = dec64le(
					buf.b + (j << 3) + (k << 5));
				uint32_t hi = dec16le(
					buf.b + (j << 3) + 128 + (k << 1));

				/* Extract sign bit. */
				uint32_t neg = -(uint32_t)(lo >> 63);
				lo &= 0x7FFFFFFFFFFFFFFF;
				hi &= 0x7FFF;

				/* Use even or odd column depending on
				   parity of t. */
				uint32_t pbit = (t[v >> 3] >> (v & 7)) & 1;
				uint64_t p_odd = -(uint64_t)pbit;
				uint32_t p_oddw = (uint32_t)p_odd;

				uint32_t r = 0;
				for (size_t i = 0; i < hi_len; i += 2) {
					uint64_t tlo0 = tab_lo[i + 0];
					uint64_t tlo1 = tab_lo[i + 1];
					uint64_t tlo = tlo0
						^ (p_odd & (tlo0 ^ tlo1));
					uint32_t cc =
						(uint32_t)((lo - tlo) >> 63);
					uint32_t thi0 = tab_hi[i + 0];
					uint32_t thi1 = tab_hi[i + 1];
					uint32_t thi = thi0
						^ (p_oddw & (thi0 ^ thi1));
					r += (hi - thi - cc) >> 31;
				}
				uint32_t hinz = (hi - 1) >> 31;
				for (size_t i = hi_len; i < lo_len; i += 2) {
					uint64_t tlo0 = tab_lo[i + 0];
					uint64_t tlo1 = tab_lo[i + 1];
					uint64_t tlo = tlo0
						^ (p_odd & (tlo0 ^ tlo1));
					uint32_t cc =
						(uint32_t)((lo - tlo) >> 63);
					r += hinz & cc;
				}

				/* Multiply by 2 and apply parity. */
				r = (r << 1) - p_oddw;

				/* Apply sign bit. */
				r = (r ^ neg) - neg;

				x[v] = (int8_t)*(int32_t *)&r;
				sn += r * r;
			}
		}
	}
	return sn;

}

/*
 * Returned value:
 *   1   first non-zero coefficient of s is positive
 *  -1   first non-zero coefficient of s is negative
 *   0   s is entirely zero
 */
static int32_t
poly_symbreak(unsigned logn, const int16_t *s)
{
	size_t n = (size_t)1 << logn;
	uint32_t r = 0;
	uint32_t c = 0xFFFFFFFF;
	for (size_t u = 0; u < n; u ++) {
		uint32_t x = (uint32_t)s[u];
		uint32_t nz = c & tbmask(x | -x);
		c &= ~nz;
		r |= nz & (tbmask(x) | 1);
	}
	return r;
}

/*
 * Encode the signature, with output length exactly sig_len bytes.
 * Padding is applied if necessary. Returned value is 1 on success, 0
 * on error; an error is reported if the signature does not fit in the
 * provided buffer.
 */
static int
encode_sig(unsigned logn, void *sig, size_t sig_len,
	const uint8_t *salt, size_t salt_len, const int16_t *s1)
{
	size_t n = (size_t)1 << logn;
	int low = (logn == 10) ? 6 : 5;
	uint8_t *buf = (uint8_t *)sig;
	size_t buf_len = sig_len;

	/* We check the minimal size, including at least n bits for the
	   variable part. */
	if (buf_len < salt_len + ((uint32_t)(low + 2) << (logn - 3))) {
		return 0;
	}

	/* Salt */
	memcpy(buf, salt, salt_len);
	buf += salt_len;
	buf_len -= salt_len;

	/* Sign bits */
	for (size_t u = 0; u < n; u += 8) {
		unsigned x = 0;
		for (size_t v = 0; v < 8; v ++) {
			x |= (*(const uint16_t *)(s1 + u + v) >> 15) << v;
		}
		buf[u >> 3] = x;
	}
	buf += (n >> 3);
	buf_len -= (n >> 3);

	/* Fixed-size parts */
	uint32_t low_mask = ((uint32_t)1 << low) - 1;
	for (size_t u = 0; u < n; u += 8) {
		uint64_t x = 0;
		for (size_t v = 0, vv = 0; v < 8; v ++, vv += low) {
			uint32_t w = (uint32_t)s1[u + v];
			w ^= tbmask(w);
			x |= (uint64_t)(w & low_mask) << vv;
		}
		for (int i = 0; i < (low << 3); i += 8) {
			*buf ++ = (uint8_t)(x >> i);
		}
	}
	buf_len -= (size_t)low << (logn - 3);

	/* Variable-size parts */
	uint32_t acc = 0;
	int acc_len = 0;
	size_t j = 0;
	for (size_t u = j; u < n; u ++) {
		uint32_t w = (uint32_t)s1[u];
		int k = (int)((w ^ tbmask(w)) >> low);
		acc |= (uint32_t)1 << (acc_len + k);
		acc_len += 1 + k;
		while (acc_len >= 8) {
			if (buf_len == 0) {
				return 0;
			}
			*buf ++ = (uint8_t)acc;
			buf_len --;
			acc >>= 8;
			acc_len -= 8;
		}
	}
	if (acc_len > 0) {
		if (buf_len == 0) {
			return 0;
		}
		*buf ++ = (uint8_t)acc;
		buf_len --;
	}

	/* Padding. */
	memset(buf, 0, buf_len);
	return 1;
}

/* see hawk.h */
void
hawk_sign_start(shake_context *sc_data)
{
	shake_init(sc_data, 256);
}

static int
sign_finish_inner(unsigned logn, int use_shake,
	void (*rng)(void *ctx, void *dst, size_t len), void *rng_context,
	void *sig, const shake_context *sc_data,
	const void *priv, size_t priv_len, void *tmp, size_t tmp_len)
{
	/*
	 * Ensure that the tmp[] buffer has proper alignment for 64-bit
	 * access.
	 */
	if (tmp_len < 7) {
		return 0;
	}
	if (logn < 8 || logn > 10) {
		return 0;
	}
	uintptr_t utmp1 = (uintptr_t)tmp;
	uintptr_t utmp2 = (utmp1 + 7) & ~(uintptr_t)7;
	tmp_len -= (size_t)(utmp2 - utmp1);
	uint32_t *tt32 = (uint32_t *)utmp2;
	if (tmp_len < ((size_t)6 << logn)) {
		return 0;
	}

	/*
	 * Check whether the private key is decoded or encoded.
	 */
	int priv_decoded;
	if (priv_len == HAWK_PRIVKEY_SIZE(logn)) {
		priv_decoded = 0;
	} else if (priv_len == HAWK_PRIVKEY_DECODED_SIZE(logn)) {
		priv_decoded = 1;
	} else {
		return 0;
	}

	/*
	 * Hawk parameters from: https://eprint.iacr.org/2022/1155
	 *
	 *   logn                   8         9       10
	 *   salt_len (bits)      112       192      320
	 *   sigma_sign         1.010     1.278    1.299
	 *   sigma_ver          1.042     1.425    1.571
	 *   sigma_sec          1.042     1.425    1.974
	 *   sigma_pk             1.1       1.5      2.0
	 *
	 * Note: sigma_ver is changed to 1.042 for logn=8 and 1.571 for
	 * logn=10 (compared to the eprint).
	 *
	 * Maximum squared norm of x = 2*n*(sigma_ver^2).
	 * max_xnorm is used for the squared norm of 2*x (as returned
	 * by sig_gauss()).
	 */
	size_t n = (size_t)1 << logn;
	size_t salt_len;
	uint32_t max_xnorm;
	switch (logn) {
	case 8:
		salt_len = 14;
		max_xnorm = 2223;
		break;
	case 9:
		salt_len = 24;
		max_xnorm = 8317;
		break;
	case 10:
		salt_len = 40;
		max_xnorm = 20218;
		break;
	default:
		return 0;
	}
	size_t seed_len = 8 + ((size_t)1 << (logn - 5));
	size_t hpub_len = (size_t)1 << (logn - 4);

	/*
	 * Memory layout:
	 *    g    n bytes
	 *    ww   2*n bytes
	 *    x0   n bytes
	 *    x1   n bytes
	 *    f    n bytes
	 *
	 * ww is the area where we will perform mod 2 computations.
	 */
	int8_t *g = (int8_t *)tt32;
	uint8_t *ww = (uint8_t *)(g + n);
	int8_t *x0 = (int8_t *)(ww + 2 * n);
	int8_t *x1 = x0 + n;
	int8_t *f = x1 + n;

	/*
	 * Re-expand the private key.
	 */
	const uint8_t *F2;
	const uint8_t *G2;
	const void *hpub;
	if (priv_decoded) {
		f = (int8_t *)priv;
		g = f + n;
		F2 = (const uint8_t *)(g + n);
		G2 = F2 + (n >> 3);
		hpub = G2 + (n >> 3);
	} else {
		const uint8_t *secbuf = (const uint8_t *)priv;
		Hawk_regen_fg(logn, f, g, secbuf);
		F2 = secbuf + seed_len;
		G2 = F2 + (n >> 3);
		hpub = G2 + (n >> 3);
	}

#if HAWK_DEBUG
	printf("#### Sign (n=%u):\n", 1u << logn);
	print_i8(logn, "f", f);
	print_i8(logn, "g", g);
	print_u1(logn, "F2", F2);
	print_u1(logn, "G2", G2);
	print_blob("hpub", hpub, hpub_len);
#endif

	/*
	 * SHAKE256 computations will use this context structure. We do
	 * not modify the input context value.
	 */
	shake_context scd;

	/*
	 * hm <- SHAKE256(message || hpub)
	 */
	uint8_t hm[64];
	scd = *sc_data;
	shake_inject(&scd, hpub, hpub_len);
	shake_flip(&scd);
	shake_extract(&scd, hm, sizeof hm);

#if HAWK_DEBUG
	printf("# hm = SHAKE256(message || hpub) (64 bytes)\n");
	print_blob("hm", hm, sizeof hm);
#endif

	for (uint32_t attempt = 0;; attempt += 2) {
		/*
		 * We get temporary values (of n/8 bytes each) inside ww.
		 *    t0    n/8 bytes
		 *    t1    n/8 bytes
		 *    h0    n/8 bytes
		 *    h1    n/8 bytes
		 *    f2    n/8 bytes
		 *    g2    n/8 bytes
		 *    xx    10*n/8 bytes
		 */
		uint8_t *t0 = ww;
		uint8_t *t1 = t0 + (n >> 3);
		uint8_t *h0 = t1 + (n >> 3);
		uint8_t *h1 = h0 + (n >> 3);
		uint8_t *f2 = h1 + (n >> 3);
		uint8_t *g2 = f2 + (n >> 3);
		uint8_t *xx = g2 + (n >> 3);

		/*
		 * Generate the salt.
		 *
		 * In normal mode (with SHAKE-based sampling), we use
		 * SHAKE256(hm || kgseed || counter || rand), i.e. the
		 * concatenation of the message, the private key seed, an
		 * attempt counter, and some extra randomness (from the
		 * provided RNG); this ensures security even if the RNG
		 * does not provide good randomness. If (f,g) were provided
		 * already decoded, we use f and g directly instead of
		 * kgseed.
		 *
		 * In alternate mode, as per the API, we use the RNG directly.
		 */
		uint8_t salt[40];
		rng(rng_context, salt, salt_len);
		if (use_shake) {
			uint8_t tbuf[4];
			enc32le(tbuf, attempt);
			shake_init(&scd, 256);
			shake_inject(&scd, hm, sizeof hm);
			if (priv_decoded) {
				shake_inject(&scd, priv, n << 1);
			} else {
				shake_inject(&scd, priv, seed_len);
			}
			shake_inject(&scd, tbuf, sizeof tbuf);
			shake_inject(&scd, salt, salt_len);
			shake_flip(&scd);
			shake_extract(&scd, salt, salt_len);
		}

#if HAWK_DEBUG
		print_blob("salt", salt, salt_len);
#endif

		/*
		 * h <- SHAKE256(hm || salt)  (out: 2*n bits)
		 */
		shake_init(&scd, 256);
		shake_inject(&scd, hm, sizeof hm);
		shake_inject(&scd, salt, salt_len);
		shake_flip(&scd);
		shake_extract(&scd, h0, n >> 2);

#if HAWK_DEBUG
		printf("# h = SHAKE256(hm || salt) (2*n bits)\n");
		print_u1(logn, "h0", h0);
		print_u1(logn, "h1", h1);
#endif

		/*
		 * t <- B*h  (mod 2)
		 */
		extract_lowbit(logn, f2, f);
		extract_lowbit(logn, g2, g);
		basis_m2_mul(logn, t0, t1, h0, h1, f2, g2, F2, G2, xx);

#if HAWK_DEBUG
		printf("# t = B*h (mod 2)\n");
		print_u1(logn, "t0", t0);
		print_u1(logn, "t1", t1);
#endif

		/*
		 * x <- D_{t/2}(sigma_sign)
		 * Reject if the squared norm of x is larger than
		 * 2*n*(sigma_ver^2) (this is quite improbable).
		 *
		 * When using SHAKE for the sampling, we also inject
		 * the message, the private key (seed) and the attempt
		 * counter into the SHAKE instances, so that security is
		 * maintained even if the RNG source is poor.
		 */
		uint32_t xsn;
		if (use_shake) {
			uint8_t tbuf[4];
			enc32le(tbuf, attempt + 1);
			shake_init(&scd, 256);
			shake_inject(&scd, hm, sizeof hm);
			if (priv_decoded) {
				shake_inject(&scd, priv, n << 1);
			} else {
				shake_inject(&scd, priv, seed_len);
			}
			shake_inject(&scd, tbuf, sizeof tbuf);
			xsn = sig_gauss(logn, rng, rng_context, &scd, x0, t0);
		} else {
			xsn = sig_gauss_alt(logn, rng, rng_context, x0, t0);
		}
#if HAWK_DEBUG
		printf("# (dx0, dx1) = (2*x0, 2*x1)\n");
		print_i8(logn, "dx0", x0);
		print_i8(logn, "dx1", x1);
		printf("l2norm(2*x)^2 = %lu\n", (unsigned long)xsn);
#endif
		if (xsn > max_xnorm) {
			continue;
		}

		/*
		 * (x0*adj(f) + x1*adj(g))/(f*adj(f) + g*adj(g)) should
		 * have all its coefficients in [-1/2, +1/2]. For
		 * degrees n=512 and n=1024, the key has been generated
		 * so that it is extremely improbable (less than
		 * 2^(-105) and 2^(-315), respectively) that this
		 * property is not met, and we can omit the test (the
		 * risk of miscomputation through a stray cosmic ray is
		 * many orders of magnitude higher).
		 *
		 * For degree n=256, the property may be false with
		 * probability about 2^(-39.5), though the actual
		 * verification failure rate may be lower. Since n=256
		 * is a "challenge" variant, we can tolerate that
		 * failure rate.
		 */

		/*
		 * s1 = (1/2)*h1 + g*x0 - f*x1
		 * We compute 2*(f*x1 - g*x0) = h1 - 2*s1
		 * (we already have 2*x in (x0,x1)).
		 *
		 * We no longer need the mod 2 values, except h1, which we
		 * move to a stack buffer (it has size at most 128 bytes,
		 * we can afford to put it on the stack).
		 *
		 * Memory layout:
		 *    g    n bytes   w1
		 *    --   n bytes   w1
		 *    --   n bytes   w2
		 *    x0   n bytes   w2
		 *    x1   n bytes   w3
		 *    f    n bytes   w3
		 */
		uint8_t h1buf[128];
		memcpy(h1buf, h1, n >> 3);

		/* w1 <- 2*g*x0 (NTT) */
		uint16_t *w1 = (uint16_t *)tt32;
		uint16_t *w2 = w1 + n;
		uint16_t *w3 = w2 + n;
		if (priv_decoded) {
			mq18433_poly_set_small(logn, w1, g);
		} else {
			mq18433_poly_set_small_inplace_low(logn, w1);
		}
		mq18433_poly_set_small_inplace_high(logn, w2);
		mq18433_NTT(logn, w1);
		mq18433_NTT(logn, w2);
		for (size_t u = 0; u < n; u ++) {
			w1[u] = mq18433_montymul(w1[u], w2[u]);
		}

		/* w3 <- 2*(f*x1 - g*x0) = h1 - 2*s1 */
		mq18433_poly_set_small(logn, w2, x1);
		if (priv_decoded) {
			mq18433_poly_set_small(logn, w3, f);
		} else {
			mq18433_poly_set_small_inplace_high(logn, w3);
		}
		mq18433_NTT(logn, w2);
		mq18433_NTT(logn, w3);
		for (size_t u = 0; u < n; u ++) {
			w3[u] = mq18433_tomonty(mq18433_sub(
				mq18433_montymul(w2[u], w3[u]), w1[u]));
		}
		mq18433_iNTT(logn, w3);
		mq18433_poly_snorm(logn, w3);

#if HAWK_DEBUG
		printf("# w = h1 - 2*s1\n");
		print_i16(logn, "w", (int16_t *)w3);
#endif

		/*
		 * sym-break(): if h1 - 2*s1 (currently in w3) is
		 * positive, then we return s1 = (h1 - w3)/2. Otherwise,
		 * we return s1 = h1 - (h1 - w3)/2 = (h1 + w3)/2.
		 * Thus, this uses a conditional negation of w3. We
		 * set nm = -1 if the negation must happen, 0 otherwise.
		 *
		 * We also enforce a limit on the individual elements of s1.
		 */
		int16_t *s1 = (int16_t *)w3;
		uint32_t ps = poly_symbreak(logn, s1);
#if HAWK_DEBUG
		printf("symbreak = %d\n", (int)ps);
#endif
		int lim = 1 << ((logn == 10) ? 10 : 9);
		uint32_t nm = ~tbmask(ps - 1);
		for (size_t u = 0; u < n; u ++) {
			uint32_t z = (uint32_t)s1[u];
			z = ((z ^ nm) - nm) + ((h1buf[u >> 3] >> (u & 7)) & 1);
			int32_t y = *(int32_t *)&z >> 1;
			if (y < -lim || y >= lim) {
				lim = 0;
				break;
			}
			s1[u] = y;
		}
		if (lim == 0) {
#if HAWK_DEBUG
			printf("# limit on s1 exceeded, restarting\n");
#endif
			if (!priv_decoded) {
				Hawk_regen_fg(logn, f, g, priv);
			}
			continue;
		}

#if HAWK_DEBUG
		print_i16(logn, "s1", (int16_t *)s1);
#endif

		/*
		 * We have the signature (in s1). We encode it into
		 * the temporary buffer, and (optionally) into the provided
		 * output buffer.
		 */
		size_t sig_len = HAWK_SIG_SIZE(logn);
		if (encode_sig(logn, tmp, sig_len, salt, salt_len, s1)) {
#if HAWK_DEBUG
			print_blob("sig", tmp, sig_len);
#endif
			if (sig != NULL) {
				memcpy(sig, tmp, sig_len);
			}
			return 1;
		}
#if HAWK_DEBUG
		printf("# signature too large, restarting");
#endif
		if (!priv_decoded) {
			Hawk_regen_fg(logn, f, g, priv);
		}
	}
}

/* see hawk.h */
int
hawk_sign_finish(unsigned logn,
	void (*rng)(void *ctx, void *dst, size_t len), void *rng_context,
	void *sig, const shake_context *sc_data,
	const void *priv, void *tmp, size_t tmp_len)
{
	return sign_finish_inner(logn, 1, rng, rng_context, sig, sc_data,
		priv, HAWK_PRIVKEY_SIZE(logn), tmp, tmp_len);
}

/* see hawk.h */
int
hawk_sign_finish_alt(unsigned logn,
	void (*rng)(void *ctx, void *dst, size_t len), void *rng_context,
	void *sig, const shake_context *sc_data,
	const void *priv, size_t priv_len, void *tmp, size_t tmp_len)
{
	return sign_finish_inner(logn, 0, rng, rng_context, sig, sc_data,
		priv, priv_len, tmp, tmp_len);
}

/* see hawk.h */
void
hawk_decode_private_key(unsigned logn, void *priv_dec, const void *priv)
{
	size_t n = (size_t)1 << logn;
	int8_t *f = (int8_t *)priv_dec;
	int8_t *g = f + n;
	size_t seed_len = 8 + ((size_t)1 << (logn - 5));
	size_t hpub_len = (size_t)1 << (logn - 4);
	uint8_t seed[40];
	memcpy(seed, priv, seed_len);
	memmove(g + n, (const uint8_t *)priv + seed_len, (n >> 2) + hpub_len);
	Hawk_regen_fg(logn, f, g, seed);
}
