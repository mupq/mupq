/*
 * Gaussian sampling.
 */

#include "sign_inner.h"

/* Union type to get easier access to values with SIMD intrinsics. */
typedef union {
	fpr f;
#if FNDSA_NEON
	float64x1_t v;
#endif
#if FNDSA_RV64D
	f64 v;
#endif
} fpr_u;

/* 1/(2*(1.8205^2)) */
#define INV_2SQRSIGMA0   FPR(5435486223186882, -55)

/* For logn = 1 to 10, n = 2^logn:
      q = 12289
      gs_norm = (117/100)*sqrt(q)
      bitsec = max(2, n/4)
      eps = 1/sqrt(bitsec*2^64)
      smoothz2n = sqrt(log(4*n*(1 + 1/eps))/pi)/sqrt(2*pi)
      sigma = smoothz2n*gs_norm
      sigma_min = sigma/gs_norm = smoothz2n
   We store precomputed values for 1/sigma and for sigma_min, indexed
   by logn. */
static const fpr_u INV_SIGMA[] = {
	{ FPR_ZERO },                     /* unused */
	{ FPR(7961475618707097, -60) },   /* 0.0069054793295940881528 */
	{ FPR(7851656902127320, -60) },   /* 0.0068102267767177965681 */
	{ FPR(7746260754658859, -60) },   /* 0.0067188101910722700565 */
	{ FPR(7595833604889141, -60) },   /* 0.0065883354370073655600 */
	{ FPR(7453842886538220, -60) },   /* 0.0064651781207602890978 */
	{ FPR(7319528409832599, -60) },   /* 0.0063486788828078985744 */
	{ FPR(7192222552237877, -60) },   /* 0.0062382586529084365056 */
	{ FPR(7071336252758509, -60) },   /* 0.0061334065020930252290 */
	{ FPR(6956347512113097, -60) },   /* 0.0060336696681577231923 */
	{ FPR(6846791885593314, -60) }    /* 0.0059386453095331150985 */
};
static const fpr_u SIGMA_MIN[] = {
	{ FPR_ZERO },                     /* unused */
	{ FPR(5028307297130123, -52) },   /* 1.1165085072329102589 */
	{ FPR(5098636688852518, -52) },   /* 1.1321247692325272406 */
	{ FPR(5168009084304506, -52) },   /* 1.1475285353733668685 */
	{ FPR(5270355833453349, -52) },   /* 1.1702540788534828940 */
	{ FPR(5370752584786614, -52) },   /* 1.1925466358390344011 */
	{ FPR(5469306724145091, -52) },   /* 1.2144300507766139921 */
	{ FPR(5566116128735780, -52) },   /* 1.2359260567719808790 */
	{ FPR(5661270305715104, -52) },   /* 1.2570545284063214163 */
	{ FPR(5754851361258101, -52) },   /* 1.2778336969128335860 */
	{ FPR(5846934829975396, -52) }    /* 1.2982803343442918540 */
};

#if !FNDSA_ASM_CORTEXM4
/* Distribution for gaussian0() (this is the RCDT table from the
   specification, expressed in base 2^24). */
static const uint32_t GAUSS0[][3] = {
	{ 10745844,  3068844,  3741698 },
	{  5559083,  1580863,  8248194 },
	{  2260429, 13669192,  2736639 },
	{   708981,  4421575, 10046180 },
	{   169348,  7122675,  4136815 },
	{    30538, 13063405,  7650655 },
	{     4132, 14505003,  7826148 },
	{      417, 16768101, 11363290 },
	{       31,  8444042,  8086568 },
	{        1, 12844466,   265321 },
	{        0,  1232676, 13644283 },
	{        0,    38047,  9111839 },
	{        0,      870,  6138264 },
	{        0,       14, 12545723 },
	{        0,        0,  3104126 },
	{        0,        0,    28824 },
	{        0,        0,      198 },
	{        0,        0,        1 }
};
#endif

/* log(2) */
#define LOG2   FPR(6243314768165359, -53)

/* 1/log(2) */
#define INV_LOG2   FPR(6497320848556798, -52)

/* We access the PRNG through macros so that they can be overridden by some
   compatiblity tests with the original Falcon implementation. */
#ifndef prng_init
#if FNDSA_SHAKE256X4
#define prng_init       shake256x4_init
#define prng_next_u8    shake256x4_next_u8
#define prng_next_u64   shake256x4_next_u64
#else
#define prng_init(pc, seed, seed_len)   do { \
		shake_init(pc, 256); \
		shake_inject(pc, seed, seed_len); \
		shake_flip(pc); \
	} while (0)
#define prng_next_u8    shake_next_u8
#define prng_next_u64   shake_next_u64
#endif
#endif

/* see sign_inner.h */
void
sampler_init(sampler_state *ss, unsigned logn,
	const void *seed, size_t seed_len)
{
	prng_init(&ss->pc, seed, seed_len);
	ss->logn = logn;
}

#if FNDSA_ASM_CORTEXM4
int32_t fndsa_gaussian0_helper(uint64_t lo, uint32_t hi);
#endif

static inline int32_t
gaussian0(sampler_state *ss)
{
	/* Get a random 72-bit value, into three 24-bit limbs (v0..v2). */
	uint64_t lo = prng_next_u64(&ss->pc);
	uint32_t hi = prng_next_u8(&ss->pc);
#if FNDSA_ASM_CORTEXM4
	return fndsa_gaussian0_helper(lo, hi);
#else
	uint32_t v0 = (uint32_t)lo & 0xFFFFFF;
	uint32_t v1 = (uint32_t)(lo >> 24) & 0xFFFFFF;
	uint32_t v2 = (uint32_t)(lo >> 48) | (hi << 16);

	/* Sampled value is z such that v0..v2 is lower than the first
	   z elements of the table. */
	int32_t z = 0;
	for (size_t i = 0; i < (sizeof GAUSS0) / sizeof(GAUSS0[0]); i ++) {
		uint32_t cc;
		cc = (v0 - GAUSS0[i][2]) >> 31;
		cc = (v1 - GAUSS0[i][1] - cc) >> 31;
		cc = (v2 - GAUSS0[i][0] - cc) >> 31;
		z += (int32_t)cc;
	}
	return z;
#endif
}

#if FNDSA_SSE2
/* ========================= SSE2 IMPLEMENTATION ========================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_SSE2
static inline int64_t
mtwop63(__m128d x)
{
#if FNDSA_64
	static const union {
		fpr f[2];
		__m128d x;
	} twop63 = { {
		FPR(4503599627370496, 11),
		FPR(4503599627370496, 11),
	} };

	return _mm_cvttsd_si64(_mm_mul_sd(x, twop63.x));
#else
	/* 32-bit x86 does not have an SSE2 opcode to convert floating-point
	   values to 64-bit integers, only 32-bit signed integers. We must
	   do the conversion in three steps with factor 2^21. */
	static const union {
		fpr f[2];
		__m128d x;
	} twop21 = { {
		FPR(4503599627370496, -31),
		FPR(4503599627370496, -31),
	} };
	x = _mm_mul_sd(x, twop21.x);
	int32_t z2 = _mm_cvttsd_si32(x);
	x = _mm_sub_sd(x, _mm_cvtsi32_sd(_mm_setzero_pd(), z2));
	x = _mm_mul_sd(x, twop21.x);
	int32_t z1 = _mm_cvttsd_si32(x);
	x = _mm_sub_sd(x, _mm_cvtsi32_sd(_mm_setzero_pd(), z1));
	x = _mm_mul_sd(x, twop21.x);
	int32_t z0 = _mm_cvttsd_si32(x);
	return ((int64_t)z2 << 42) + ((int64_t)z1 << 21) + (int64_t)z0;
#endif
}

/* Compute ccs*exp(-x)*2^63, rounded to an integer. This function assumes
   that 0 <= x < log(2), and 0 <= ccs <= 1. It returns a value in [0,2^63]. */
TARGET_SSE2
static inline uint64_t
expm_p63(__m128d x, __m128d ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;
#if FNDSA_64
	/* On 64-bit x86, we have 64x64->128 multiplication, then we can use
	   it, it's normally constant-time.
	   MSVC uses a different syntax for this operation. */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
#if defined _MSC_VER
		y = EXPM_COEFFS[i] - __umulh(z, y);
#else
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#endif
	}
#if defined _MSC_VER
	y = __umulh(w, y);
#else
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
#endif
#else
	/* On 32-bit x86, no 64x64->128 multiplication, we must use
	   four 32x32->64 multiplications. */
	uint32_t z0 = (uint32_t)z, z1 = (uint32_t)(z >> 32);
	uint32_t w0 = (uint32_t)w, w1 = (uint32_t)(w >> 32);

	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
		uint64_t f = (uint64_t)z0 * (uint64_t)y0;
		uint64_t a = (uint64_t)z0 * (uint64_t)y1 + (f >> 32);
		uint64_t b = (uint64_t)z1 * (uint64_t)y0;
		uint64_t c = (a >> 32) + (b >> 32)
			+ (((uint64_t)(uint32_t)a
			  + (uint64_t)(uint32_t)b) >> 32)
			+ (uint64_t)z1 * (uint64_t)y1;
		y = EXPM_COEFFS[i] - c;
	}
	uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
	uint64_t f = (uint64_t)w0 * (uint64_t)y0;
	uint64_t a = (uint64_t)w0 * (uint64_t)y1 + (f >> 32);
	uint64_t b = (uint64_t)w1 * (uint64_t)y0;
	y = (a >> 32) + (b >> 32)
		+ (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32)
		+ (uint64_t)w1 * (uint64_t)y1;
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
TARGET_SSE2
static inline int
ber_exp(sampler_state *ss, __m128d x, __m128d ccs)
{
	static union { fpr f[2]; __m128d x; }
		LOG2_u = { { LOG2, LOG2 } },
		INV_LOG2_u = { { INV_LOG2, INV_LOG2 } };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use a truncating conversion because
	   x >= 0. Moreover, x is small, so we can stick to 32-bit values. */
	int32_t si = _mm_cvttsd_si32(_mm_mul_sd(x, INV_LOG2_u.x));
	__m128d r = _mm_sub_sd(x,
		_mm_mul_sd(_mm_cvtsi32_sd(_mm_setzero_pd(), si), LOG2_u.x));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

TARGET_SSE2
static int32_t
sampler_next_sse2(sampler_state *ss, __m128d mu, __m128d isigma)
{
	static union { fpr f[2]; __m128d x; }
		HALF_u = { {
			FPR(4503599627370496, -53),
			FPR(4503599627370496, -53)
		} },
		INV_2SQRSIGMA0_u = { {
			INV_2SQRSIGMA0, INV_2SQRSIGMA0
		} };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int32_t s = _mm_cvttsd_si32(mu);
	s -= _mm_comilt_sd(mu, _mm_cvtsi32_sd(_mm_setzero_pd(), s));
	__m128d r = _mm_sub_sd(mu, _mm_cvtsi32_sd(_mm_setzero_pd(), s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	__m128d dss = _mm_mul_sd(_mm_mul_sd(isigma, isigma), HALF_u.x);

	/* css = sigma_min / sigma = sigma_min * isigma  */
	__m128d ccs = _mm_mul_sd(isigma,
		_mm_load_sd((const double *)SIGMA_MIN + ss->logn));

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		__m128d x = _mm_sub_sd(_mm_cvtsi32_sd(_mm_setzero_pd(), z), r);
		x = _mm_mul_sd(_mm_mul_sd(x, x), dss);
		x = _mm_sub_sd(x, _mm_mul_sd(
			_mm_cvtsi32_sd(_mm_setzero_pd(), z0 * z0),
			INV_2SQRSIGMA0_u.x));
		if (ber_exp(ss, x, ccs)) {
			return s + z;
		}
	}
}

/* see sign_inner.h */
TARGET_SSE2
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_sse2(ss,
		_mm_load_sd((double *)&mu),
		_mm_load_sd((double *)&isigma));
}

#elif FNDSA_NEON
/* ========================= NEON IMPLEMENTATION ========================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_NEON
static inline int64_t
mtwop63(float64x1_t x)
{
	static const fpr_u twop63 = { FPR(4503599627370496, 11) };
	return vget_lane_s64(vcvt_s64_f64(vmul_f64(x, twop63.v)), 0);
}

/* Compute ccs*exp(-x)*2^63, rounded to an integer. This function assumes
   that 0 <= x < log(2), and 0 <= ccs <= 1. It returns a value in [0,2^63]. */
TARGET_NEON
static inline uint64_t
expm_p63(float64x1_t x, float64x1_t ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;

	/* We assume here that 64x64->128 multiplications are constant-time,
	   which is not exactly true on some aarch64 systems (e.g. ARM
	   Cortex A53 and A55 return the result one cycle earlier when
	   the operands fit on 32 bits).
	   ARM compiler calls the 128-bit type '__uint128_t' while GCC
	   and Clang use 'unsigned __int128'. */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
#if defined __GNUC__ || defined __clang__
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#else
		__uint128_t c = (__uint128_t)z * (__uint128_t)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
#endif
	}
#if defined __GNUC__ || defined __clang__
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
#else
	y = (uint64_t)(((__uint128_t)w * (__uint128_t)y) >> 64);
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
TARGET_NEON
static inline int
ber_exp(sampler_state *ss, float64x1_t x, float64x1_t ccs)
{
	static const fpr_u LOG2_u = { LOG2 };
	static const fpr_u INV_LOG2_u = { INV_LOG2 };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use a truncating conversion because
	   x >= 0. Moreover, x is small, so we can stick to 32-bit values. */
	int32_t si = (int32_t)vget_lane_s64(
		vcvt_s64_f64(vmul_f64(x, INV_LOG2_u.v)), 0);
	float64x1_t r = vsub_f64(x,
		vmul_f64(vcvt_f64_s64(vcreate_s64(si)), LOG2_u.v));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

TARGET_NEON
static int32_t
sampler_next_neon(sampler_state *ss, float64x1_t mu, float64x1_t isigma)
{
	static const fpr_u HALF_u = { FPR(4503599627370496, -53) };
	static const fpr_u INV_2SQRSIGMA0_u = { INV_2SQRSIGMA0 };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int32_t s = (int32_t)vcvtmd_s64_f64(mu);
	float64x1_t r = vsub_f64(mu, vcvt_f64_s64(vcreate_s64(s)));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	float64x1_t dss = vmul_f64(vmul_f64(isigma, isigma), HALF_u.v);

	/* css = sigma_min / sigma = sigma_min * isigma  */
	float64x1_t ccs = vmul_f64(isigma, SIGMA_MIN[ss->logn].v);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		float64x1_t x = vsub_f64(vcvt_f64_s64(vcreate_s64(z)), r);
		x = vmul_f64(vmul_f64(x, x), dss);
		x = vsub_f64(x, vmul_f64(
			vcvt_f64_s64(vcreate_s64(z0 * z0)),
			INV_2SQRSIGMA0_u.v));
		if (ber_exp(ss, x, ccs)) {
			return s + z;
		}
	}
}

/* see sign_inner.h */
TARGET_NEON
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_neon(ss,
		vld1_f64((const float64_t *)&mu),
		vld1_f64((const float64_t *)&isigma));
}

#elif FNDSA_RV64D
/* ========================= RISC-V IMPLEMENTATION ======================= */

/* Input: 0 <= x < log(2)
   Output: trunc(x*2^63) */
TARGET_NEON
static inline int64_t
mtwop63(f64 x)
{
	static const fpr_u twop63 = { FPR(4503599627370496, 11) };
	return f64_trunc(f64_mul(x, twop63.v));
}

static inline uint64_t
expm_p63(f64 x, f64 ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)mtwop63(x) << 1;
	uint64_t w = (uint64_t)mtwop63(ccs) << 1;

	/* We assume here that 64x64->128 multiplications are constant-time
	   and that the compiler is GCC/Clang compatible (i.e. supports
	   the 'unsigned __int128' type). */
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		unsigned __int128 c =
			(unsigned __int128)z * (unsigned __int128)y;
		y = EXPM_COEFFS[i] - (uint64_t)(c >> 64);
	}
	y = (uint64_t)(((unsigned __int128)w * (unsigned __int128)y) >> 64);
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
static inline int
ber_exp(sampler_state *ss, f64 x, f64 ccs)
{
	static const fpr_u LOG2_u = { LOG2 };
	static const fpr_u INV_LOG2_u = { INV_LOG2 };

	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use f64_trunc() because x >= 0. */
	int32_t si = (int32_t)f64_trunc(f64_mul(x, INV_LOG2_u.v));
	f64 r = f64_sub(x, f64_mul(f64_of(si), LOG2_u.v));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

static int32_t
sampler_next_rv64d(sampler_state *ss, f64 mu, f64 isigma)
{
	static const fpr_u INV_2SQRSIGMA0_u = { INV_2SQRSIGMA0 };

	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int64_t s = f64_floor(mu);
	f64 r = f64_sub(mu, f64_of(s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	f64 dss = f64_half(f64_sqr(isigma));

	/* css = sigma_min / sigma = sigma_min * isigma  */
	f64 ccs = f64_mul(isigma, SIGMA_MIN[ss->logn].v);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		f64 x = f64_mul(f64_sqr(f64_sub(f64_of(z), r)), dss);
		x = f64_sub(x, f64_mul(f64_of(z0 * z0), INV_2SQRSIGMA0_u.v));
		if (ber_exp(ss, x, ccs)) {
			return (int32_t)s + z;
		}
	}
}

/* see sign_inner.h */
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	return sampler_next_rv64d(ss, f64_from_raw(mu), f64_from_raw(isigma));
}

#else
/* ========================= PLAIN IMPLEMENTATION ======================== */

static inline uint64_t
expm_p63(fpr x, fpr ccs)
{
	/* The polynomial approximation of exp(-x) is from FACCT:
	      https://eprint.iacr.org/2018/1234
	   Specifically, the values are extracted from the implementation
	   referenced by the FACCT paper, available at:
	      https://github.com/raykzhao/gaussian  */
	static const uint64_t EXPM_COEFFS[] = {
		0x00000004741183A3,
		0x00000036548CFC06,
		0x0000024FDCBF140A,
		0x0000171D939DE045,
		0x0000D00CF58F6F84,
		0x000680681CF796E3,
		0x002D82D8305B0FEA,
		0x011111110E066FD0,
		0x0555555555070F00,
		0x155555555581FF00,
		0x400000000002B400,
		0x7FFFFFFFFFFF4800,
		0x8000000000000000
	};

	/* TODO: maybe use 64x64->128 multiplications if available? It
	   is a bit tricky to decide, because the plain code is used for
	   unknown architectures, and we do not know if the larger
	   multiplication is constant-time (it often happens that it
	   is not). */

	uint64_t y = EXPM_COEFFS[0];
	uint64_t z = (uint64_t)fpr_trunc(fpr_mul2e(x, 63)) << 1;
	uint32_t z0 = (uint32_t)z, z1 = (uint32_t)(z >> 32);
	for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i ++) {
		uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
#if FNDSA_ASM_CORTEXM4
		uint32_t tt, r0, r1;
		__asm__(
			"umull	%0, %2, %3, %5\n\t"
			"umull	%0, %1, %3, %6\n\t"
			"umaal	%2, %0, %4, %5\n\t"
			"umaal	%0, %1, %4, %6\n\t"
			: "=&r" (r0), "=&r" (r1), "=&r" (tt)
			: "r" (y0), "r" (y1), "r" (z0), "r" (z1));
		y = EXPM_COEFFS[i] - ((uint64_t)r0 | ((uint64_t)r1 << 32));
#else
		uint64_t f = (uint64_t)z0 * (uint64_t)y0;
		uint64_t a = (uint64_t)z0 * (uint64_t)y1 + (f >> 32);
		uint64_t b = (uint64_t)z1 * (uint64_t)y0;
		uint64_t c = (a >> 32) + (b >> 32)
			+ (((uint64_t)(uint32_t)a
			  + (uint64_t)(uint32_t)b) >> 32)
			+ (uint64_t)z1 * (uint64_t)y1;
		y = EXPM_COEFFS[i] - c;
#endif
	}

	/* The scaling factor must be applied at the end. Since y is now
	   in fixed-point notation, we have to convert the factor to the
	   same format, and we do an extra integer multiplication. */
	uint64_t w = (uint64_t)fpr_trunc(fpr_mul2e(ccs, 63)) << 1;
	uint32_t w0 = (uint32_t)w, w1 = (uint32_t)(w >> 32);
	uint32_t y0 = (uint32_t)y, y1 = (uint32_t)(y >> 32);
#if FNDSA_ASM_CORTEXM4
	uint32_t tt, r0, r1;
	__asm__(
		"umull	%0, %2, %3, %5\n\t"
		"umull	%0, %1, %3, %6\n\t"
		"umaal	%2, %0, %4, %5\n\t"
		"umaal	%0, %1, %4, %6\n\t"
		: "=&r" (r0), "=&r" (r1), "=&r" (tt)
		: "r" (y0), "r" (y1), "r" (w0), "r" (w1));
	y = (uint64_t)r0 | ((uint64_t)r1 << 32);
#else
	uint64_t f = (uint64_t)w0 * (uint64_t)y0;
	uint64_t a = (uint64_t)w0 * (uint64_t)y1 + (f >> 32);
	uint64_t b = (uint64_t)w1 * (uint64_t)y0;
	y = (a >> 32) + (b >> 32)
		+ (((uint64_t)(uint32_t)a + (uint64_t)(uint32_t)b) >> 32)
		+ (uint64_t)w1 * (uint64_t)y1;
#endif
	return y;
}

/* Sample a bit with probability ccs*exp(-x) (for x >= 0). */
static inline int
ber_exp(sampler_state *ss, fpr x, fpr ccs)
{
	/* Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	   and 0 <= r < log(2). We can use fpr_trunc() because x >= 0
	   (fpr_trunc() is presumably a bit faster than fpr_floor()). */
	int32_t si = (int32_t)fpr_trunc(fpr_mul(x, INV_LOG2));
	fpr r = fpr_sub(x, fpr_mul(fpr_of(si), LOG2));

	/* If s >= 64, sigma = 1.2, r = 0 and b = 1, then we get s >= 64
	   if the half-Gaussian produced z >= 13, which happens with
	   probability about 2^(-32). When s >= 64, ber_exp() will return
	   true with probability less than 2^(-64), so we can simply
	   saturate s at 63 (the bias introduced here is lower than 2^(-96),
	   and would require about 2^192 samplings to be detectable, which
	   is way beyond the formal bound of 2^64 signatures with the
	   same key. */
	uint32_t s = (uint32_t)si;
	s |= (uint32_t)(63 - s) >> 26;

	/* Compute ccs*exp(-x). Since x = s*log(2) + r, we compute
	   ccs*exp(-r)/2^s. We know that 0 <= r < log(2), so we can
	   use expm_p63(), which yields a result scaled by 63 bits. We
	   scale it up 1 bit further, then right-shift by s bits.

	   We subtract 1 to make sure that the value fits on 64 bits
	   (i.e. if r = 0 then we may get 2^64 and we prefer 2^64-1
	   in that case, to avoid the overflow). The bias is negligible
	   since expm_p63() has precision only 51 bits or so. */
	uint64_t z = fpr_ursh((expm_p63(r, ccs) << 1) - 1, s);

	/* Sample a bit. We lazily compare the value z with a uniform 64-bit
	   integer, consuming only as many bytes as necessary. Since the PRNG
	   is cryptographically strong, we leak no information from the
	   conditional jumps below. */
	for (int i = 56; i >= 0; i -= 8) {
		unsigned w = prng_next_u8(&ss->pc);
		unsigned bz = (unsigned)(z >> i) & 0xFF;
		if (w != bz) {
			return w < bz;
		}
	}
	return 0;
}

/* see sign_inner.h */
int32_t
sampler_next(sampler_state *ss, fpr mu, fpr isigma)
{
	/* Split center mu into s + r, for an integer s, and 0 <= r < 1. */
	int64_t s = fpr_floor(mu);
	fpr r = fpr_sub(mu, fpr_of(s));

	/* dss = 1/(2*sigma^2) = 0.5*(isigma^2)  */
	fpr dss = fpr_half(fpr_sqr(isigma));

	/* css = sigma_min / sigma = sigma_min * isigma  */
	fpr ccs = fpr_mul(isigma, SIGMA_MIN[ss->logn].f);

	/* We sample on centre r. */
	for (;;) {
		/* Sample z for a Gaussian distribution (non-negative only),
		   then get a random bit b to turn the sampling into a
		   bimodal distribution (we use z+1 if b = 1, or -z
		   otherwise). */
		int32_t z0 = gaussian0(ss);
		int32_t b = prng_next_u8(&ss->pc) & 1;
		int32_t z = b + ((b << 1) - 1) * z0;

		/* Rejection sampling. We want a Gaussian centred on r,
		   but we sampled against a bimodal distribution (with
		   "centres" at 0 and 1). However, we know that z is
		   always in the range where our sampling distribution is
		   greater than the Gaussian distribution, so rejection works.

		   We got z from distribution:
		      G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		   We target distribution:
		      S(z) = exp(-((z-r)^2)/(2*signa^2))
		   Rejection sampling works by keeping the value z with
		   probability S(z)/G(z), and starting again otherwise.
		   This requires S(z) <= G(z), which is the case here.
		   Thus, we simply need to keep our z with probability:
		      P = exp(-x)
		   where:
		      x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		   Here, we scale up the Bernouilli distribution, which
		   makes rejection more probable, but also makes the
		   rejection rate sufficiently decorrelated from the Gaussian
		   centre and standard deviation, so that measurement of the
		   rejection rate does not leak enough usable information
		   to attackers (which is how the implementation can claim
		   to be "constant-time").  */
		fpr x = fpr_mul(fpr_sqr(fpr_sub(fpr_of(z), r)), dss);
		x = fpr_sub(x, fpr_mul(fpr_of(z0 * z0), INV_2SQRSIGMA0));
		if (ber_exp(ss, x, ccs)) {
			return (int32_t)s + z;
		}
	}
}
#endif

TARGET_SSE2 TARGET_NEON
static void
ffsamp_fft_inner(sampler_state *ss, unsigned logn,
	fpr *t0, fpr *t1, fpr *g00, fpr *g01, fpr *g11, fpr *tmp)
{
	/* When logn = 1, arrays have length 2; we unroll the last steps. */
	if (logn == 1) {
#if FNDSA_SSE2
		static const union {
			fpr f[2];
			__m128d x;
		} one_u = { { FPR_ONE, FPR_ONE } };
		__m128d cz = _mm_castsi128_pd(
			_mm_setr_epi32(0, 0, 0, -0x80000000));

		/* Decompose G into LDL. g00 and g11 are self-adjoint,
		   thus only one (real) coefficient each. */
		__m128d g00_re = _mm_load_sd((double *)g00);
		__m128d g01_cc = _mm_loadu_pd((double *)g01);
		__m128d g11_re = _mm_load_sd((double *)g11);
		__m128d inv_g00_re = _mm_div_sd(one_u.x, g00_re);
		__m128d inv_g00 = _mm_shuffle_pd(inv_g00_re, inv_g00_re, 0);
		__m128d mu = _mm_mul_pd(g01_cc, inv_g00);
		__m128d zo = _mm_mul_pd(mu, g01_cc);
		__m128d zo_re = _mm_add_sd(zo, _mm_shuffle_pd(zo, zo, 1));
		__m128d d00_re = g00_re;
		__m128d l01 = _mm_xor_pd(cz, mu);
		__m128d d11_re = _mm_sub_sd(g11_re, zo_re);

		/* No split on d00 and d11, since they are one-coeff each. */

		/* The half-size Gram matrices for the recursive LDL tree
		   exploration are now:
		     - left sub-tree:   d00_re, zero, d00_re
		     - right sub-tree:  d11_re, zero, d11_re
		   t1 split is trivial. */
		__m128d w = _mm_loadu_pd((double *)t1);
		__m128d w0 = w;
		__m128d w1 = _mm_shuffle_pd(w, w, 3);
		__m128d leaf = _mm_mul_sd(
			_mm_sqrt_sd(_mm_setzero_pd(), d11_re),
			_mm_load_sd((const double *)INV_SIGMA + ss->logn));
		__m128d y0 = _mm_cvtsi32_sd(_mm_setzero_pd(),
			sampler_next_sse2(ss, w0, leaf));
		__m128d y1 = _mm_cvtsi32_sd(_mm_setzero_pd(),
			sampler_next_sse2(ss, w1, leaf));

		/* Merge is trivial, since logn = 1. */

		/* At this point:
		     t0 and t1 are unmodified; t1 is also [w0, w1]
		     l10 is in [l10_re, l10_im]
		     z1 is [y0, y1]
		   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
		   z1 is moved into t1. */
		__m128d y = _mm_shuffle_pd(y0, y1, 0);
		__m128d a = _mm_sub_pd(w, y);
		__m128d b1 = _mm_mul_pd(a, _mm_xor_pd(cz, l01));
		__m128d b2 = _mm_mul_pd(a, _mm_shuffle_pd(l01, l01, 1));
		__m128d b = _mm_add_pd(
			_mm_shuffle_pd(b1, b2, 2),
			_mm_shuffle_pd(b1, b2, 1));
		__m128d x = _mm_add_pd(b, _mm_loadu_pd((double *)t0));
		_mm_storeu_pd((double *)t1, y);

		/* Second recursive invocation, on the split tb0, using
		   the left sub-tree. tb0 is [x0, x1], and the split is
		   trivial since logn = 1. */
		__m128d x0 = x;
		__m128d x1 = _mm_shuffle_pd(x, x, 3);
		leaf = _mm_mul_sd(
			_mm_sqrt_sd(_mm_setzero_pd(), d00_re),
			_mm_load_sd((const double *)INV_SIGMA + ss->logn));
		x0 = _mm_cvtsi32_sd(_mm_setzero_pd(),
			sampler_next_sse2(ss, x0, leaf));
		x1 = _mm_cvtsi32_sd(_mm_setzero_pd(),
			sampler_next_sse2(ss, x1, leaf));
		_mm_store_sd((double *)t0, x0);
		_mm_store_sd((double *)t0 + 1, x1);
#elif FNDSA_NEON
		static const fpr_u one_u = { FPR_ONE };
		static const union { fpr f[2]; float64x2_t x; }
			cz = { { FPR_ZERO, FPR_NZERO } };

		/* Decompose G into LDL. g00 and g11 are self-adjoint,
		   thus only one (real) coefficient each. */
		float64x1_t g00_re = vld1_f64((float64_t *)g00);
		float64x2_t g01_cc = vld1q_f64((float64_t *)g01);
		float64x1_t g11_re = vld1_f64((float64_t *)g11);
		float64x1_t inv_g00_re = vdiv_f64(one_u.v, g00_re);
		float64x2_t inv_g00 = vdupq_lane_f64(inv_g00_re, 0);
		float64x2_t mu = vmulq_f64(g01_cc, inv_g00);
		float64x2_t zo = vmulq_f64(mu, g01_cc);
		float64x1_t zo_re = vget_low_f64(vpaddq_f64(zo, zo));
		float64x1_t d00_re = g00_re;
		float64x2_t l01 = vreinterpretq_f64_u64(
			veorq_u64(cz.x, vreinterpretq_u64_f64(mu)));
		float64x1_t d11_re = vsub_f64(g11_re, zo_re);

		/* No split on d00 and d11, since they are one-coeff each. */

		/* The half-size Gram matrices for the recursive LDL tree
		   exploration are now:
		     - left sub-tree:   d00_re, zero, d00_re
		     - right sub-tree:  d11_re, zero, d11_re
		   t1 split is trivial. */
		float64x2_t w = vld1q_f64((const float64_t *)t1);
		float64x1_t w0 = vget_low_f64(w);
		float64x1_t w1 = vget_high_f64(w);
		float64x1_t leaf = vmul_f64(
			vsqrt_f64(d11_re),
			INV_SIGMA[ss->logn].v);
		float64x1_t y0 = vcvt_f64_s64(vcreate_s64(
			sampler_next_neon(ss, w0, leaf)));
		float64x1_t y1 = vcvt_f64_s64(vcreate_s64(
			sampler_next_neon(ss, w1, leaf)));

		/* Merge is trivial, since logn = 1. */

		/* At this point:
		     t0 and t1 are unmodified; t1 is also [w0, w1]
		     l10 is in [l10_re, l10_im]
		     z1 is [y0, y1]
		   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
		   z1 is moved into t1. */
		float64x2_t y = vcombine_f64(y0, y1);
		float64x2_t a = vsubq_f64(w, y);
		float64x2_t b1 = vmulq_f64(a,
			vreinterpretq_f64_u64(veorq_u64(cz.x,
				vreinterpretq_u64_f64(l01))));
		float64x2_t b2 = vmulq_f64(a, vextq_f64(l01, l01, 1));
		float64x2_t b = vpaddq_f64(b1, b2);
		float64x2_t x = vaddq_f64(b, vld1q_f64((const float64_t *)t0));
		vst1q_f64((float64_t *)t1, y);

		/* Second recursive invocation, on the split tb0, using
		   the left sub-tree. tb0 is [x0, x1], and the split is
		   trivial since logn = 1. */
		float64x1_t x0 = vget_low_f64(x);
		float64x1_t x1 = vget_high_f64(x);
		leaf = vmul_f64(
			vsqrt_f64(d00_re),
			INV_SIGMA[ss->logn].v);
		x0 = vcvt_f64_s64(vcreate_s64(
			sampler_next_neon(ss, x0, leaf)));
		x1 = vcvt_f64_s64(vcreate_s64(
			sampler_next_neon(ss, x1, leaf)));
		vst1_f64((float64_t *)t0, x0);
		vst1_f64((float64_t *)t0 + 1, x1);
#elif FNDSA_RV64D
		/* Decompose G into LDL. g00 and g11 are self-adjoint,
		   thus only one (real) coefficient each. */
		f64 g00_re = ((const f64 *)g00)[0];
		f64 g01_re = ((const f64 *)g01)[0];
		f64 g01_im = ((const f64 *)g01)[1];
		f64 g11_re = ((const f64 *)g11)[0];
		f64 inv_g00_re = f64_inv(g00_re);
		f64 mu_re = f64_mul(g01_re, inv_g00_re);
		f64 mu_im = f64_mul(g01_im, inv_g00_re);
		f64 zo_re = f64_add(
			f64_mul(mu_re, g01_re),
			f64_mul(mu_im, g01_im));
		f64 d00_re = g00_re;
		f64 l01_re = mu_re;
		f64 l01_im = f64_neg(mu_im);
		f64 d11_re = f64_sub(g11_re, zo_re);

		/* No split on d00 and d11, since they are one-coeff each. */

		/* The half-size Gram matrices for the recursive LDL tree
		   exploration are now:
		     - left sub-tree:   d00_re, zero, d00_re
		     - right sub-tree:  d11_re, zero, d11_re
		   t1 split is trivial. */
		f64 w0 = ((const f64 *)t1)[0];
		f64 w1 = ((const f64 *)t1)[1];
		f64 leaf = f64_mul(f64_sqrt(d11_re), INV_SIGMA[ss->logn].v);
		f64 y0 = f64_of(sampler_next_rv64d(ss, w0, leaf));
		f64 y1 = f64_of(sampler_next_rv64d(ss, w1, leaf));

		/* Merge is trivial, since logn = 1. */

		/* At this point:
		     t0 and t1 are unmodified; t1 is also [w0, w1]
		     l10 is in [l10_re, l10_im]
		     z1 is [y0, y1]
		   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
		   z1 is moved into t1. */
		f64 a_re = f64_sub(w0, y0);
		f64 a_im = f64_sub(w1, y1);
		f64 b_re = f64_sub(
			f64_mul(a_re, l01_re),
			f64_mul(a_im, l01_im));
		f64 b_im = f64_add(
			f64_mul(a_im, l01_re),
			f64_mul(a_re, l01_im));
		f64 x0 = f64_add(((const f64 *)t0)[0], b_re);
		f64 x1 = f64_add(((const f64 *)t0)[1], b_im);
		((f64 *)t1)[0] = y0;
		((f64 *)t1)[1] = y1;

		/* Second recursive invocation, on the split tb0, using
		   the left sub-tree. tb0 is [x0, x1], and the split is
		   trivial since logn = 1. */
		leaf = f64_mul(f64_sqrt(d00_re), INV_SIGMA[ss->logn].v);
		((f64 *)t0)[0] = f64_of(sampler_next_rv64d(ss, x0, leaf));
		((f64 *)t0)[1] = f64_of(sampler_next_rv64d(ss, x1, leaf));
#else
		/* Decompose G into LDL. g00 and g11 are self-adjoint,
		   thus only one (real) coefficient each. */
		fpr g00_re = g00[0];
		fpr g01_re = g01[0], g01_im = g01[1];
		fpr g11_re = g11[0];
		fpr inv_g00_re = fpr_inv(g00_re);
		fpr mu_re = fpr_mul(g01_re, inv_g00_re);
		fpr mu_im = fpr_mul(g01_im, inv_g00_re);
		fpr zo_re = fpr_add(
			fpr_mul(mu_re, g01_re),
			fpr_mul(mu_im, g01_im));
		fpr d00_re = g00_re;
		fpr l01_re = mu_re;
		fpr l01_im = fpr_neg(mu_im);
		fpr d11_re = fpr_sub(g11_re, zo_re);

		/* No split on d00 and d11, since they are one-coeff each. */

		/* The half-size Gram matrices for the recursive LDL tree
		   exploration are now:
		     - left sub-tree:   d00_re, zero, d00_re
		     - right sub-tree:  d11_re, zero, d11_re
		   t1 split is trivial. */
		fpr w0 = t1[0];
		fpr w1 = t1[1];
		fpr leaf = fpr_mul(fpr_sqrt(d11_re), INV_SIGMA[ss->logn].f);
		fpr y0 = fpr_of(sampler_next(ss, w0, leaf));
		fpr y1 = fpr_of(sampler_next(ss, w1, leaf));

		/* Merge is trivial, since logn = 1. */

		/* At this point:
		     t0 and t1 are unmodified; t1 is also [w0, w1]
		     l10 is in [l10_re, l10_im]
		     z1 is [y0, y1]
		   Compute tb0 = t0 + (t1 - z1)*l10  (into [x0, x1]).
		   z1 is moved into t1. */
		fpr a_re = fpr_sub(w0, y0);
		fpr a_im = fpr_sub(w1, y1);
		fpr b_re, b_im;
		FPC_MUL(b_re, b_im, a_re, a_im, l01_re, l01_im);
		fpr x0 = fpr_add(t0[0], b_re);
		fpr x1 = fpr_add(t0[1], b_im);
		t1[0] = y0;
		t1[1] = y1;

		/* Second recursive invocation, on the split tb0, using
		   the left sub-tree. tb0 is [x0, x1], and the split is
		   trivial since logn = 1. */
		leaf = fpr_mul(fpr_sqrt(d00_re), INV_SIGMA[ss->logn].f);
		t0[0] = fpr_of(sampler_next(ss, x0, leaf));
		t0[1] = fpr_of(sampler_next(ss, x1, leaf));
#endif
		return;
	}

	/* General case: logn >= 2 */
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* Decompose G into LDL; the decomposed matrix replaces G. */
	fpoly_LDL_fft(logn, g00, g01, g11);

	/* Split d00 and d11 (currently in g00 and g11) and expand them
	   into half-size quasi-cyclic Gram matrices. We also save l10
	   (currently in g01) into tmp. */
	fpr *w0 = tmp;
	fpr *w1 = w0 + hn;
	fpoly_split_selfadj_fft(logn, w0, w1, g00);
	memcpy(g00, w0, n * sizeof(fpr));
	fpoly_split_selfadj_fft(logn, w0, w1, g11);
	memcpy(g11, w0, n * sizeof(fpr));
	memcpy(tmp, g01, n * sizeof(fpr));
	memcpy(g01, g00, hn * sizeof(fpr));
	memcpy(g01 + hn, g11, hn * sizeof(fpr));

	/* The half-size Gram matrices for the recursive LDL tree
	   exploration are now:
	     - left sub-tree:   g00[0..hn], g00[hn..n], g01[0..hn]
	     - right sub-tree:  g11[0..hn], g11[hn..n], g01[hn..n]
	   l10 is in tmp[0..n]. */
	fpr *left_00 = g00;
	fpr *left_01 = g00 + hn;
	fpr *right_00 = g11;
	fpr *right_01 = g11 + hn;
	fpr *left_11 = g01;
	fpr *right_11 = g01 + hn;

	/* We split t1 and use the first recursive call on the two
	   halves, using the right sub-tree. The result is merged
	   back into tmp[2*n..3*n]. */
	w0 = tmp + n;
	w1 = w0 + hn;
	fpr *w2 = w1 + hn;
	fpoly_split_fft(logn, w0, w1, t1);
	ffsamp_fft_inner(ss, logn - 1, w0, w1,
		right_00, right_01, right_11, w2);
	fpoly_merge_fft(logn, w2, w0, w1);

	/* At this point:
	     t0 and t1 are unmodified
	     l10 is in tmp[0..n]
	     z1 is in tmp[2*n..3*n]
	   We compute tb0 = t0 + (t1 - z1)*l10.
	   tb0 is written over t0.
	   z1 is moved into t1.
	   l10 is scratched. */
	fpr *l10 = tmp;
	fpr *w = l10 + n;
	fpr *z1 = w + n;
	memcpy(w, t1, n * sizeof(fpr));
	fpoly_sub(logn, w, z1);
	memcpy(t1, z1, n * sizeof(fpr));
	fpoly_mul_fft(logn, l10, w);
	fpoly_add(logn, t0, l10);

	/* Second recursive invocation, on the split tb0 (currently in t0),
	   using the left sub-tree.
	   tmp is free. */
	w0 = tmp;
	w1 = w0 + hn;
	w2 = w1 + hn;
	fpoly_split_fft(logn, w0, w1, t0);
	ffsamp_fft_inner(ss, logn - 1,
		w0, w1, left_00, left_01, left_11, w2);
	fpoly_merge_fft(logn, t0, w0, w1);
}

/* see sign_inner.h */
void
ffsamp_fft(sampler_state *ss,
	fpr *t0, fpr *t1, fpr *g00, fpr *g01, fpr *g11, fpr *tmp)
{
	ffsamp_fft_inner(ss, ss->logn, t0, t1, g00, g01, g11, tmp);
}
