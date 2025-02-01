/*
 * NTRU equation.
 */

#include "kgen_inner.h"

#define Q   12289

static const uint16_t MAX_BL_SMALL[11] = {
	1, 1, 2, 3, 4, 8, 14, 27, 53, 104, 207
};
static const uint16_t MAX_BL_LARGE[10] = {
	1, 2, 3, 6, 11, 21, 40, 78, 155, 308
};
static const uint16_t WORD_WIN[10] = {
	1, 1, 2, 2, 2, 3, 3, 4, 5, 7
};
static const uint16_t MIN_SAVE_FG[11] = {
	0, 0, 1, 2, 2, 2, 2, 2, 2, 3, 3
};

/* Convert source f and g into RNS+NTT, at the start of the provided tmp[]
   (one word per coefficient). */
static void
make_fg_zero(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	uint32_t *ft = tmp;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	poly_mp_set_small(logn, ft, f, p);
	poly_mp_set_small(logn, gt, g, p);
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	mp_NTT(logn, ft, gm, p, p0i);
	mp_NTT(logn, gt, gm, p, p0i);
}

#if FNDSA_AVX2
TARGET_AVX2
static void
avx2_make_fg_zero(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	uint32_t *ft = tmp;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	avx2_poly_mp_set_small(logn, ft, f, p);
	avx2_poly_mp_set_small(logn, gt, g, p);
	avx2_mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	avx2_mp_NTT(logn, ft, gm, p, p0i);
	avx2_mp_NTT(logn, gt, gm, p, p0i);
}
#endif

/* One step of computing (f,g) at a given depth.
     Input: (f,g) of degree 2^(logn_top-depth)
     Output: (f',g') of degree 2^(logn_top-(depth+1))
   Input and output values are at the start of tmp[], in RNS+NTT notation.
  
   RAM USAGE: 3*(2^logn_top) (at most)
   (assumptions: max_bl_small[0] = max_bl_small[1] = 1, max_bl_small[2] = 2) */
static void
make_fg_step(unsigned logn_top, unsigned depth, uint32_t *tmp)
{
	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t slen = MAX_BL_SMALL[depth];
	size_t tlen = MAX_BL_SMALL[depth + 1];

	/* Layout:
	     fd    output f' (hn*tlen)
	     gd    output g' (hn*tlen)
	     fs    source (n*slen)
	     gs    source (n*slen)
	     t1    NTT support (n)
	     t2    extra (max(n, slen - n))  */
	uint32_t *fd = tmp;
	uint32_t *gd = fd + hn * tlen;
	uint32_t *fs = gd + hn * tlen;
	uint32_t *gs = fs + n * slen;
	uint32_t *t1 = gs + n * slen;
	uint32_t *t2 = t1 + n;
	memmove(fs, tmp, 2 * n * slen * sizeof *tmp);

	/* First slen words: we use the input values directly, and apply
	   inverse NTT as we go, so that we get the sources in RNS (non-NTT). */
	uint32_t *xf = fs;
	uint32_t *xg = gs;
	uint32_t *yf = fd;
	uint32_t *yg = gd;
	for (size_t i = 0; i < slen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		for (size_t j = 0; j < hn; j ++) {
			yf[j] = mp_mmul(
				mp_mmul(xf[2 * j], xf[2 * j + 1], p, p0i),
				R2, p, p0i);
			yg[j] = mp_mmul(
				mp_mmul(xg[2 * j], xg[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		mp_mkigm(logn, t1, PRIMES[i].ig, p, p0i);
		mp_iNTT(logn, xf, t1, p, p0i);
		mp_iNTT(logn, xg, t1, p, p0i);
		xf += n;
		xg += n;
		yf += hn;
		yg += hn;
	}

	/* Now that fs and gs are in RNS, rebuild their plain integer
	   coefficients. */
	zint_rebuild_CRT(fs, slen, n, 2, 1, t1);

	/* Remaining output words. */
	for (size_t i = slen; i < tlen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
		mp_mkgm(logn, t1, PRIMES[i].g, p, p0i);
		for (size_t j = 0; j < n; j ++) {
			t2[j] = zint_mod_small_signed(
				fs + j, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t1, p, p0i);
		for (size_t j = 0; j < hn; j ++) {
			yf[j] = mp_mmul(
				mp_mmul(t2[2 * j], t2[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		yf += hn;
		for (size_t j = 0; j < n; j ++) {
			t2[j] = zint_mod_small_signed(
				gs + j, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t1, p, p0i);
		for (size_t j = 0; j < hn; j ++) {
			yg[j] = mp_mmul(
				mp_mmul(t2[2 * j], t2[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		yg += hn;
	}
}

#if FNDSA_AVX2
TARGET_AVX2
static void
avx2_make_fg_step(unsigned logn_top, unsigned depth, uint32_t *tmp)
{
	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t slen = MAX_BL_SMALL[depth];
	size_t tlen = MAX_BL_SMALL[depth + 1];

	/* Layout:
	     fd    output f' (hn*tlen)
	     gd    output g' (hn*tlen)
	     fs    source (n*slen)
	     gs    source (n*slen)
	     t1    NTT support (n)
	     t2    extra (max(n, slen - n))  */
	uint32_t *fd = tmp;
	uint32_t *gd = fd + hn * tlen;
	uint32_t *fs = gd + hn * tlen;
	uint32_t *gs = fs + n * slen;
	uint32_t *t1 = gs + n * slen;
	uint32_t *t2 = t1 + n;
	memmove(fs, tmp, 2 * n * slen * sizeof *tmp);

	/* First slen words: we use the input values directly, and apply
	   inverse NTT as we go, so that we get the sources in RNS (non-NTT). */
	uint32_t *xf = fs;
	uint32_t *xg = gs;
	uint32_t *yf = fd;
	uint32_t *yg = gd;
	for (size_t i = 0; i < slen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		for (size_t j = 0; j < hn; j ++) {
			yf[j] = mp_mmul(
				mp_mmul(xf[2 * j], xf[2 * j + 1], p, p0i),
				R2, p, p0i);
			yg[j] = mp_mmul(
				mp_mmul(xg[2 * j], xg[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		mp_mkigm(logn, t1, PRIMES[i].ig, p, p0i);
		mp_iNTT(logn, xf, t1, p, p0i);
		mp_iNTT(logn, xg, t1, p, p0i);
		xf += n;
		xg += n;
		yf += hn;
		yg += hn;
	}

	/* Now that fs and gs are in RNS, rebuild their plain integer
	   coefficients. */
	avx2_zint_rebuild_CRT(fs, slen, n, 2, 1, t1);

	/* Remaining output words. */
	for (size_t i = slen; i < tlen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
		avx2_mp_mkgm(logn, t1, PRIMES[i].g, p, p0i);
		if (logn >= 3) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			__m256i yRx = _mm256_set1_epi32(Rx);
			for (size_t j = 0; j < n; j += 8) {
				__m256i yt = zint_mod_small_signed_x8(
					fs + j, slen, n, yp, yp0i, yR2, yRx);
				_mm256_storeu_si256((__m256i *)(t2 + j), yt);
			}
			avx2_mp_NTT(logn, t2, t1, p, p0i);
			for (size_t j = 0; j < hn; j += 4) {
				__m256i yt = _mm256_loadu_si256(
					(__m256i *)(t2 + (2 * j)));
				yt = mp_mmul_x4(yt,
					_mm256_srli_epi64(yt, 32), yp, yp0i);
				yt = mp_mmul_x4(yt, yR2, yp, yp0i);
				yt = _mm256_shuffle_epi32(yt, 0xD8);
				yt = _mm256_permute4x64_epi64(yt, 0xD8);
				_mm_storeu_si128((__m128i *)(yf + j),
					_mm256_castsi256_si128(yt));
			}
			yf += hn;
			for (size_t j = 0; j < n; j += 8) {
				__m256i yt = zint_mod_small_signed_x8(
					gs + j, slen, n, yp, yp0i, yR2, yRx);
				_mm256_storeu_si256((__m256i *)(t2 + j), yt);
			}
			avx2_mp_NTT(logn, t2, t1, p, p0i);
			for (size_t j = 0; j < hn; j += 4) {
				__m256i yt = _mm256_loadu_si256(
					(__m256i *)(t2 + (2 * j)));
				yt = mp_mmul_x4(yt,
					_mm256_srli_epi64(yt, 32), yp, yp0i);
				yt = mp_mmul_x4(yt, yR2, yp, yp0i);
				yt = _mm256_shuffle_epi32(yt, 0xD8);
				yt = _mm256_permute4x64_epi64(yt, 0xD8);
				_mm_storeu_si128((__m128i *)(yg + j),
					_mm256_castsi256_si128(yt));
			}
			yg += hn;
			continue;
		}
		for (size_t j = 0; j < n; j ++) {
			t2[j] = zint_mod_small_signed(
				fs + j, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t2, t1, p, p0i);
		for (size_t j = 0; j < hn; j ++) {
			yf[j] = mp_mmul(
				mp_mmul(t2[2 * j], t2[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		yf += hn;
		for (size_t j = 0; j < n; j ++) {
			t2[j] = zint_mod_small_signed(
				gs + j, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t2, t1, p, p0i);
		for (size_t j = 0; j < hn; j ++) {
			yg[j] = mp_mmul(
				mp_mmul(t2[2 * j], t2[2 * j + 1], p, p0i),
				R2, p, p0i);
		}
		yg += hn;
	}
}
#endif

/* Compute (f,g) at a specified depth, in RNS+NTT notation.
   Computed values are stored at the start of the provided tmp[] (slen
   words per coefficient).
  
   This function is for depth < logn_top. For the deepest layer, use
   make_fg_deepest().
  
   RAM USAGE: 3*(2^logn_top) */
static void
make_fg_intermediate(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *tmp)
{
	make_fg_zero(logn_top, f, g, tmp);
	for (unsigned d = 0; d < depth; d ++) {
		make_fg_step(logn_top, d, tmp);
	}
}

#if FNDSA_AVX2
TARGET_AVX2
static void
avx2_make_fg_intermediate(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *tmp)
{
	avx2_make_fg_zero(logn_top, f, g, tmp);
	for (unsigned d = 0; d < depth; d ++) {
		avx2_make_fg_step(logn_top, d, tmp);
	}
}
#endif

/* Compute (f,g) at the deepest level (i.e. get Res(f,X^n+1) and
   Res(g,X^n+1)). Intermediate (f,g) values (below the save threshold)
   are copied at the end of tmp (of size save_off words).
  
   If f is not invertible modulo X^n+1 and modulo p = 2147473409, then
   this function returns 0 (the resultants and intermediate values are
   still computed). Otherwise it returns 1. */
static int
make_fg_deepest(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *tmp, size_t sav_off)
{
	make_fg_zero(logn_top, f, g, tmp);
	int r = 1;

	/* f is now in RNS+NTT, so we can test its invertibility (mod p)
	   by simply checking that all NTT coefficients are non-zero.
	   (This invertibility allows recovering of G from f, g and F
	   by working modulo p = 2147473409. It is not actually necessary
	   since the recovery of G is normally done modulo q = 12289, but
	   it was done in the original ntrugen code, so it is maintained
	   here for full compatibility of test vectors.) */
	size_t n = (size_t)1 << logn_top;
	uint32_t b = 0;
	for (size_t i = 0; i < n; i ++) {
		b |= tmp[i] - 1;
	}
	r = (int)(1 - (b >> 31));

	for (unsigned d = 0; d < logn_top; d ++) {
		make_fg_step(logn_top, d, tmp);

		/* make_fg_step() computes the (f,g) for depth d+1; we
		   save that value if d+1 is at least at the save
		   threshold, but is not the deepest level. */
		unsigned d2 = d + 1;
		if (d2 < logn_top && d2 >= MIN_SAVE_FG[logn_top]) {
			size_t slen = MAX_BL_SMALL[d2];
			size_t fglen = slen << (logn_top + 1 - d2);
			sav_off -= fglen;
			memmove(tmp + sav_off, tmp, fglen * sizeof *tmp);
		}
	}
	return r;
}

#if FNDSA_AVX2
TARGET_AVX2
static int
avx2_make_fg_deepest(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *tmp, size_t sav_off)
{
	avx2_make_fg_zero(logn_top, f, g, tmp);
	int r = 1;

	/* f is now in RNS+NTT, so we can test its invertibility (mod p)
	   by simply checking that all NTT coefficients are non-zero.
	   (This invertibility allows recovering of G from f, g and F
	   by working modulo p = 2147473409. It is not actually necessary
	   since the recovery of G is normally done modulo q = 12289, but
	   it was done in the original ntrugen code, so it is maintained
	   here for full compatibility of test vectors.) */
	size_t n = (size_t)1 << logn_top;
	uint32_t b = 0;
	for (size_t i = 0; i < n; i ++) {
		b |= tmp[i] - 1;
	}
	r = (int)(1 - (b >> 31));

	for (unsigned d = 0; d < logn_top; d ++) {
		avx2_make_fg_step(logn_top, d, tmp);

		/* make_fg_step() computes the (f,g) for depth d+1; we
		   save that value if d+1 is at least at the save
		   threshold, but is not the deepest level. */
		unsigned d2 = d + 1;
		if (d2 < logn_top && d2 >= MIN_SAVE_FG[logn_top]) {
			size_t slen = MAX_BL_SMALL[d2];
			size_t fglen = slen << (logn_top + 1 - d2);
			sav_off -= fglen;
			memmove(tmp + sav_off, tmp, fglen * sizeof *tmp);
		}
	}
	return r;
}
#endif

/* Error code: no error (so far) */
#define SOLVE_OK           0

/* Error code: GCD(Res(f,X^n+1), Res(g,X^n+1)) != 1 */
#define SOLVE_ERR_GCD      -1

/* Error code: reduction error (NTRU equation no longer fulfilled) */
#define SOLVE_ERR_REDUCE   -2

/* Error code: output (F,G) coefficients are off-limits */
#define SOLVE_ERR_LIMIT    -3

/* Solve the NTRU equation at the deepest level. This computes the
   integers F and G such that Res(f,X^n+1)*G - Res(g,X^n+1)*F = q.
   The two integers are written into tmp[].
  
   Returned value: 0 on success, a negative error code otherwise.
  
   RAM USAGE: max(3*(2^logn_top), 8*max_bl_small[depth]) */
static int
solve_NTRU_deepest(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp)
{
	/* Get (f,g) at the deepest level (i.e. Res(f,X^n+1) and Res(g,X^n+1)).
	   Obtained (f,g) are in RNS+NTT (since degree n = 1, this is
	   equivalent to RNS). */
	if (!make_fg_deepest(logn_top, f, g, tmp, (size_t)6 << logn_top)) {
		return SOLVE_ERR_GCD;
	}

	/* Reorganize memory:
	      Fp   output F (len)
	      Gp   output G (len)
	      fp   Res(f,X^n+1) (len)
	      gp   Res(g,X^n+1) (len)
	      t1   rest of temporary */
	size_t len = MAX_BL_SMALL[logn_top];
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + len;
	uint32_t *fp = Gp + len;
	uint32_t *gp = fp + len;
	uint32_t *t1 = gp + len;
	memmove(fp, tmp, 2 * len * sizeof *tmp);

	/* Convert back the resultants into plain integers. */
	zint_rebuild_CRT(fp, len, 1, 2, 0, t1);

	/* Apply the binary GCD to get a solution (F,G) such that:
	     f*G - g*F = 1  */
	if (!zint_bezout(Gp, Fp, fp, gp, len, t1)) {
		return SOLVE_ERR_GCD;
	}

	/* Multiply the obtained (F,G) by q to get a proper solution:
	     f*G - g*F = q  */
	if (zint_mul_small(Fp, len, Q) != 0
		|| zint_mul_small(Gp, len, Q) != 0)
	{
		return SOLVE_ERR_REDUCE;
	}
	return SOLVE_OK;
}
#if FNDSA_AVX2
TARGET_AVX2
static int
avx2_solve_NTRU_deepest(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp)
{
	/* Get (f,g) at the deepest level (i.e. Res(f,X^n+1) and Res(g,X^n+1)).
	   Obtained (f,g) are in RNS+NTT (since degree n = 1, this is
	   equivalent to RNS). */
	if (!avx2_make_fg_deepest(logn_top, f, g, tmp, (size_t)6 << logn_top)) {
		return SOLVE_ERR_GCD;
	}

	/* Reorganize memory:
	      Fp   output F (len)
	      Gp   output G (len)
	      fp   Res(f,X^n+1) (len)
	      gp   Res(g,X^n+1) (len)
	      t1   rest of temporary */
	size_t len = MAX_BL_SMALL[logn_top];
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + len;
	uint32_t *fp = Gp + len;
	uint32_t *gp = fp + len;
	uint32_t *t1 = gp + len;
	memmove(fp, tmp, 2 * len * sizeof *tmp);

	/* Convert back the resultants into plain integers. */
	avx2_zint_rebuild_CRT(fp, len, 1, 2, 0, t1);

	/* Apply the binary GCD to get a solution (F,G) such that:
	     f*G - g*F = 1  */
	if (!zint_bezout(Gp, Fp, fp, gp, len, t1)) {
		return SOLVE_ERR_GCD;
	}

	/* Multiply the obtained (F,G) by q to get a proper solution:
	     f*G - g*F = q  */
	if (zint_mul_small(Fp, len, Q) != 0
		|| zint_mul_small(Gp, len, Q) != 0)
	{
		return SOLVE_ERR_REDUCE;
	}
	return SOLVE_OK;
}
#endif

/* We use poly_sub_scaled() when log(n) < MIN_LOGN_FGNTT, and
   poly_sub_scaled_ntt() when log(n) >= MIN_LOGN_FGNTT. The NTT variant
   is faster at large degrees, but not at small degrees. */
#define MIN_LOGN_FGNTT   4

/* Solving the NTRU equation, intermediate level.
   Input is (F,G) from one level deeper (half-degree), in plain
   representation, at the start of tmp[]; output is (F,G) from this
   level, written at the start of tmp[].
  
   Returned value: 0 on success, a negative error code otherwise. */
static int
solve_NTRU_intermediate(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *restrict tmp)
{
	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* slen   size for (f,g) at this level (also output (F,G))
	   llen   size for unreduced (F,G) at this level
	   dlen   size for input (F,G) from deeper level
	   Note: we always have llen >= dlen */
	size_t slen = MAX_BL_SMALL[depth];
	size_t llen = MAX_BL_LARGE[depth];
	size_t dlen = MAX_BL_SMALL[depth + 1];

	/* Fd   F from deeper level (dlen*hn)
	   Gd   G from deeper level (dlen*hn)
	   ft   f from this level (slen*n)
	   gt   g from this level (slen*n) */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + dlen * hn;
	uint32_t *fgt = Gd + dlen * hn;

	/* Get (f,g) for this level (in RNS+NTT). */
	if (depth < MIN_SAVE_FG[logn_top]) {
		make_fg_intermediate(logn_top, f, g, depth, fgt);
	} else {
		uint32_t *sav_fg = tmp + ((size_t)6 << logn_top);
		for (unsigned d = MIN_SAVE_FG[logn_top];
			d <= depth; d ++)
		{
			sav_fg -= MAX_BL_SMALL[d] << (logn_top + 1 - d);
		}
		memmove(fgt, sav_fg, 2 * slen * n * sizeof *fgt);
	}

	/* Move buffers so that we have room for the unreduced (F,G) at
	   this level.
	     Ft   F from this level (unreduced) (llen*n)
	     Gt   G from this level (unreduced) (llen*n)
	     ft   f from this level (slen*n)
	     gt   g from this level (slen*n)
	     Fd   F from deeper level (dlen*hn)
	     Gd   G from deeper level (dlen*hn)  */
	uint32_t *Ft = tmp;
	uint32_t *Gt = Ft + llen * n;
	uint32_t *ft = Gt + llen * n;
	uint32_t *gt = ft + slen * n;
	Fd = gt + slen * n;
	Gd = Fd + dlen * hn;
	uint32_t *t1 = Gd + dlen * hn;
	memmove(ft, fgt, 2 * n * slen * sizeof *ft);
	memmove(Fd, tmp, 2 * hn * dlen * sizeof *tmp);

	/* Convert Fd and Gd to RNS, with output temporarily stored
	   in (Ft, Gt). Fd and Gd have degree hn only; we store the
	   values for each modulus p in the _last_ hn slots of the
	   n-word line for that modulus. */
	for (size_t i = 0; i < llen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		uint32_t Rx = mp_Rx31((unsigned)dlen, p, p0i, R2);
		uint32_t *xt = Ft + i * n + hn;
		uint32_t *yt = Gt + i * n + hn;
		for (size_t j = 0; j < hn; j ++) {
			xt[j] = zint_mod_small_signed(Fd + j, dlen, hn,
				p, p0i, R2, Rx);
			yt[j] = zint_mod_small_signed(Gd + j, dlen, hn,
				p, p0i, R2, Rx);
		}
	}

	/* Fd and Gd are no longer needed. */
	t1 = Fd;

	/* Compute (F,G) (unreduced) modulo sufficiently many small primes.
	   We also un-NTT (f,g) as we go; when slen primes have been
	   processed, we obtain (f,g) in RNS, and we apply the CRT to
	   get (f,g) in plain representation. */
	for (size_t i = 0; i < llen; i ++) {
		/* If we have processed exactly slen primes, then (f,g)
		   are in RNS, and we can rebuild them. */
		if (i == slen) {
			zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
		}

		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;

		/* Memory layout: we keep Ft, Gt, ft and gt; we append:
		     gm    NTT support (n)
		     igm   iNTT support (n)
		     fx    temporary f mod p (NTT) (n)
		     gx    temporary g mod p (NTT) (n)  */
		uint32_t *gm = t1;
		uint32_t *igm = gm + n;
		uint32_t *fx = igm + n;
		uint32_t *gx = fx + n;
		mp_mkgmigm(logn, gm, igm, PRIMES[i].g, PRIMES[i].ig, p, p0i);
		if (i < slen) {
			memcpy(fx, ft + i * n, n * sizeof *fx);
			memcpy(gx, gt + i * n, n * sizeof *gx);
			mp_iNTT(logn, ft + i * n, igm, p, p0i);
			mp_iNTT(logn, gt + i * n, igm, p, p0i);
		} else {
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			for (size_t j = 0; j < n; j ++) {
				fx[j] = zint_mod_small_signed(ft + j, slen, n,
					p, p0i, R2, Rx);
				gx[j] = zint_mod_small_signed(gt + j, slen, n,
					p, p0i, R2, Rx);
			}
			mp_NTT(logn, fx, gm, p, p0i);
			mp_NTT(logn, gx, gm, p, p0i);
		}

		/* We have (F,G) from deeper level in Ft and Gt, in
		   RNS. We apply the NTT modulo p. */
		uint32_t *Fe = Ft + i * n;
		uint32_t *Ge = Gt + i * n;
		mp_NTT(logn - 1, Fe + hn, gm, p, p0i);
		mp_NTT(logn - 1, Ge + hn, gm, p, p0i);

		/* Compute F and G (unreduced) modulo p. */
		for (size_t j = 0; j < hn; j ++) {
			uint32_t fa = fx[(j << 1) + 0];
			uint32_t fb = fx[(j << 1) + 1];
			uint32_t ga = gx[(j << 1) + 0];
			uint32_t gb = gx[(j << 1) + 1];
			uint32_t mFp = mp_mmul(Fe[j + hn], R2, p, p0i);
			uint32_t mGp = mp_mmul(Ge[j + hn], R2, p, p0i);
			Fe[(j << 1) + 0] = mp_mmul(gb, mFp, p, p0i);
			Fe[(j << 1) + 1] = mp_mmul(ga, mFp, p, p0i);
			Ge[(j << 1) + 0] = mp_mmul(fb, mGp, p, p0i);
			Ge[(j << 1) + 1] = mp_mmul(fa, mGp, p, p0i);
		}

		/* We want the new (F,G) in RNS only (no NTT). */
		mp_iNTT(logn, Fe, igm, p, p0i);
		mp_iNTT(logn, Ge, igm, p, p0i);
	}

	/* Edge case: if slen == llen, then we have not rebuilt (f,g)
	   into plain representation yet, so we do it now. */
	if (slen == llen) {
		zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
	}

	/* We now have the unreduced (F,G) in RNS. We rebuild their
	   plain representation. */
	zint_rebuild_CRT(Ft, llen, n, 2, 1, t1);

	/* We now reduce these (F,G) with Babai's nearest plane
	   algorithm. The reduction conceptually goes as follows:
	     k <- round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
	     (F, G) <- (F - k*f, G - k*g)
	   We use fixed-point approximations of (f,g) and (F, G) to get
	   a value k as a small polynomial with scaling; we then apply
	   k on the full-width polynomial. Each iteration "shaves" a
	   a few bits off F and G.
	  
	   We apply the process sufficiently many times to reduce (F, G)
	   to the size of (f, g) with a reasonable probability of success.
	   Since we want full constant-time processing, the number of
	   iterations and the accessed slots work on some assumptions on
	   the sizes of values (sizes have been measured over many samples,
	   and a margin of 5 times the standard deviation). */

	/* If depth is at least 2, and we will use the NTT to subtract
	   k*(f,g) from (F,G), then we will need to convert (f,g) to NTT over
	   slen+1 words, which requires an extra word to ft and gt. */
	int use_sub_ntt = (depth > 1 && logn >= MIN_LOGN_FGNTT);
	if (use_sub_ntt) {
		memmove(gt + n, gt, n * slen * sizeof *gt);
		gt += n;
		t1 += 2 * n;
	}

	/* New layout:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     ft    f from this level (slen*n) (+n if use_sub_ntt)
	     gt    g from this level (slen*n) (+n if use_sub_ntt)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (hn fxr = n)  */

	fxr *rt3 = (fxr *)t1;
	fxr *rt4 = rt3 + n;
	fxr *rt1 = rt4 + n;

	/* We consider only the top rlen words of (f,g). */
	size_t rlen = WORD_WIN[depth];
	if (rlen > slen) {
		rlen = slen;
	}
	size_t blen = slen - rlen;
	uint32_t *ftb = ft + blen * n;
	uint32_t *gtb = gt + blen * n;
	uint32_t scale_fg = 31 * (uint32_t)blen;
	uint32_t scale_FG = 31 * (uint32_t)llen;

	/* Convert f and g into fixed-point approximations, in rt3 and rt4,
	   respectively. They are scaled down by 2^(scale_fg + scale_x).
	   scale_fg is public (it depends only on the recursion depth), but
	   scale_x comes from a measurement on the actual values of (f,g) and
	   is thus secret.
	  
	   The value scale_x is adjusted so that the largest coefficient is
	   close to, but lower than, some limit t (in absolute value). The
	   limit t is chosen so that f*adj(f) + g*adj(g) does not overflow,
	   i.e. all coefficients must remain below 2^31.
	  
	   Let n be the degree; we know that n <= 2^10. The squared norm
	   of a polynomial is the sum of the squared norms of the
	   coefficients, with the squared norm of a complex number being
	   the product of that number with its complex conjugate. If all
	   coefficients of f are less than t (in absolute value), then
	   the squared norm of f is less than n*t^2. The squared norm of
	   FFT(f) (f in FFT representation) is exactly n times the
	   squared norm of f, so this leads to n^2*t^2 as a maximum
	   bound. adj(f) has the same norm as f. This implies that each
	   complex coefficient of FFT(f) has a maximum squared norm of
	   n^2*t^2 (with a maximally imbalanced polynomial with all
	   coefficient but one being zero). The computation of f*adj(f)
	   exactly is, in FFT representation, the product of each
	   coefficient with its conjugate; thus, the coefficients of
	   f*adj(f), in FFT representation, are at most n^2*t^2.
	  
	   Since we want the coefficients of f*adj(f)+g*adj(g) not to exceed
	   2^31, we need n^2*t^2 <= 2^30, i.e. n*t <= 2^15. We can adjust t
	   accordingly (called scale_t in the code below). We also need to
	   take care that t must not exceed scale_x. Approximation of f and
	   g are extracted with scale scale_fg + scale_x - scale_t, and
	   later fixed by dividing them by 2^scale_t. */
	uint32_t scale_xf = poly_max_bitlength(logn, ftb, rlen);
	uint32_t scale_xg = poly_max_bitlength(logn, gtb, rlen);
	uint32_t scale_x = scale_xf;
	scale_x ^= (scale_xf ^ scale_xg) & tbmask(scale_xf - scale_xg);
	uint32_t scale_t = 15 - logn;
	scale_t ^= (scale_t ^ scale_x) & tbmask(scale_x - scale_t);
	uint32_t scdiff = scale_x - scale_t;

	poly_big_to_fixed(logn, rt3, ftb, rlen, scdiff);
	poly_big_to_fixed(logn, rt4, gtb, rlen, scdiff);

	/* rt3 <- adj(f)/(f*adj(f) + g*adj(g))  (FFT)
	   rt4 <- adj(g)/(f*adj(f) + g*adj(g))  (FFT) */
	vect_FFT(logn, rt3);
	vect_FFT(logn, rt4);
	vect_norm_fft(logn, rt1, rt3, rt4);
	vect_mul2e(logn, rt3, scale_t);
	vect_mul2e(logn, rt4, scale_t);
	for (size_t i = 0; i < hn; i ++) {
		rt3[i] = fxr_div(rt3[i], rt1[i]);
		rt3[i + hn] = fxr_div(fxr_neg(rt3[i + hn]), rt1[i]);
		rt4[i] = fxr_div(rt4[i], rt1[i]);
		rt4[i + hn] = fxr_div(fxr_neg(rt4[i + hn]), rt1[i]);
	}

	/* New layout:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     ft    f from this level (slen*n) (+n if use_sub_ntt)
	     gt    g from this level (slen*n) (+n if use_sub_ntt)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (n fxr = 2*n)     |   k    (n)
	     rt2   (n fxr = 2*n)     |   t2   (3*n)
	   Exception: at depth == 1, we omit ft and gt:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (n fxr = 2*n)     |   k    (n)
	     rt2   (n fxr = 2*n)     |   t2   (3*n)  */
	if (depth == 1) {
		t1 = ft;
		fxr *nrt3 = (fxr *)t1;
		memmove(nrt3, rt3, 2 * n * sizeof *rt3);
		rt3 = nrt3;
		rt4 = rt3 + n;
		rt1 = rt4 + n;
	}
	int32_t *k = (int32_t *)rt1;
	uint32_t *t2 = (uint32_t *)(k + n);
	fxr *rt2 = (fxr *)t2;
	if (rt2 < (rt1 + n)) {
		rt2 = rt1 + n;
	}

	/* If we are going to use poly_sub_scaled_ntt(), then we convert
	   f and g to the NTT representation. Since poly_sub_scaled_ntt()
	   itself will use more than n*(slen+2) words in t2[], we can do
	   the same here. */
	if (use_sub_ntt) {
		uint32_t *gm = t2;
		uint32_t *tn = gm + n;
		for (size_t i = 0; i <= slen; i ++) {
			uint32_t p = PRIMES[i].p;
			uint32_t p0i = PRIMES[i].p0i;
			uint32_t R2 = PRIMES[i].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			mp_mkgm(logn, gm, PRIMES[i].g, p, p0i);
			for (size_t j = 0; j < n; j ++) {
				tn[j] = zint_mod_small_signed(
					ft + j, slen, n, p, p0i, R2, Rx);
			}
			mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(ft, tn, (slen + 1) * n * sizeof *tn);
		for (size_t i = 0; i <= slen; i ++) {
			uint32_t p = PRIMES[i].p;
			uint32_t p0i = PRIMES[i].p0i;
			uint32_t R2 = PRIMES[i].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			mp_mkgm(logn, gm, PRIMES[i].g, p, p0i);
			for (size_t j = 0; j < n; j ++) {
				tn[j] = zint_mod_small_signed(
					gt + j, slen, n, p, p0i, R2, Rx);
			}
			mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(gt, tn, (slen + 1) * n * sizeof *tn);
	}

	/* Reduce F and G repeatedly. */
	size_t FGlen = llen;
	uint32_t reduce_bits;
	switch (logn_top) {
	case 9:  reduce_bits = 13; break;
	case 10: reduce_bits = 11; break;
	default: reduce_bits = 16; break;
	}
	for (;;) {
		/* Convert the current F and G into fixed-point. We want
		   to apply scaling scale_FG + scale_x. */
		uint32_t tlen, toff;
		DIVREM31(tlen, toff, scale_FG);
		poly_big_to_fixed(logn, rt1,
			Ft + tlen * n, FGlen - tlen, scale_x + toff);
		poly_big_to_fixed(logn, rt2,
			Gt + tlen * n, FGlen - tlen, scale_x + toff);

		/* rt2 <- (F*adj(f) + G*adj(g)) / (f*adj(f) + g*adj(g)) */
		vect_FFT(logn, rt1);
		vect_FFT(logn, rt2);
		vect_mul_fft(logn, rt1, rt3);
		vect_mul_fft(logn, rt2, rt4);
		vect_add(logn, rt2, rt1);
		vect_iFFT(logn, rt2);

		/* k <- round(rt2) */
		for (size_t i = 0; i < n; i ++) {
			k[i] = fxr_round(rt2[i]);
		}

		/* (f,g) are scaled by scale_fg + scale_x
		   (F,G) are scaled by scale_FG + scale_x
		   Thus, k is scaled by scale_FG - scale_fg, which is public. */
		uint32_t scale_k = scale_FG - scale_fg;
		if (depth == 1) {
			poly_sub_kfg_scaled_depth1(logn_top, Ft, Gt, FGlen,
				(uint32_t *)k, scale_k, f, g, t2);
		} else if (use_sub_ntt) {
			poly_sub_scaled_ntt(logn, Ft, FGlen, ft, slen,
				k, scale_k, t2);
			poly_sub_scaled_ntt(logn, Gt, FGlen, gt, slen,
				k, scale_k, t2);
		} else {
			poly_sub_scaled(logn, Ft, FGlen, ft, slen, k, scale_k);
			poly_sub_scaled(logn, Gt, FGlen, gt, slen, k, scale_k);
		}

		/* We now assume that F and G have shrunk by at least
		   reduce_bits. We adjust FGlen accordinly. */
		if (scale_FG <= scale_fg) {
			break;
		}
		if (scale_FG <= (scale_fg + reduce_bits)) {
			scale_FG = scale_fg;
		} else {
			scale_FG -= reduce_bits;
		}
		while (FGlen > slen
			&& 31 * (FGlen - slen) > scale_FG - scale_fg + 30)
		{
			FGlen --;
		}
	}

	/* Output F is already in the right place; G is in Gt, and must be
	   moved back a bit. */
	memmove(tmp + slen * n, Gt, slen * n * sizeof *tmp);
	Gt = tmp + slen * n;

	/* Reduction is done. We test the current solution modulo a single
	   prime.
	   Exception: we cannot do that if depth == 1, since in that case
	   we did not keep (ft,gt). Reduction errors rarely occur at this
	   stage, so we can omit that test (depth-0 test will cover it).
	  
	   If use_sub_ntt != 0, then ft and gt are already in NTT
	   representation. */
	if (depth == 1) {
		return SOLVE_OK;
	}

	t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
	mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	if (use_sub_ntt) {
		t1 = ft;
		for (size_t i = 0; i < n; i ++) {
			t2[i] = zint_mod_small_signed(
				Gt + i, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t i = 0; i < n; i ++) {
			t1[i] = zint_mod_small_signed(
				ft + i, slen, n, p, p0i, R2, Rx);
			t2[i] = zint_mod_small_signed(
				Gt + i, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t1, t4, p, p0i);
		mp_NTT(logn, t2, t4, p, p0i);
	}
	for (size_t i = 0; i < n; i ++) {
		t3[i] = mp_mmul(t1[i], t2[i], p, p0i);
	}
	if (use_sub_ntt) {
		t1 = gt;
		for (size_t i = 0; i < n; i ++) {
			t2[i] = zint_mod_small_signed(
				Ft + i, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t i = 0; i < n; i ++) {
			t1[i] = zint_mod_small_signed(
				gt + i, slen, n, p, p0i, R2, Rx);
			t2[i] = zint_mod_small_signed(
				Ft + i, slen, n, p, p0i, R2, Rx);
		}
		mp_NTT(logn, t1, t4, p, p0i);
		mp_NTT(logn, t2, t4, p, p0i);
	}
	uint32_t rv = mp_mmul(Q, 1, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		uint32_t x = mp_mmul(t1[i], t2[i], p, p0i);
		if (mp_sub(t3[i], x, p) != rv) {
			return SOLVE_ERR_REDUCE;
		}
	}

	return SOLVE_OK;
}

#if FNDSA_AVX2
TARGET_AVX2
static int
avx2_solve_NTRU_intermediate(unsigned logn_top,
	const int8_t *restrict f, const int8_t *restrict g,
	unsigned depth, uint32_t *restrict tmp)
{
	unsigned logn = logn_top - depth;
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* slen   size for (f,g) at this level (also output (F,G))
	   llen   size for unreduced (F,G) at this level
	   dlen   size for input (F,G) from deeper level
	   Note: we always have llen >= dlen */
	size_t slen = MAX_BL_SMALL[depth];
	size_t llen = MAX_BL_LARGE[depth];
	size_t dlen = MAX_BL_SMALL[depth + 1];

	/* Fd   F from deeper level (dlen*hn)
	   Gd   G from deeper level (dlen*hn)
	   ft   f from this level (slen*n)
	   gt   g from this level (slen*n) */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + dlen * hn;
	uint32_t *fgt = Gd + dlen * hn;

	/* Get (f,g) for this level (in RNS+NTT). */
	if (depth < MIN_SAVE_FG[logn_top]) {
		avx2_make_fg_intermediate(logn_top, f, g, depth, fgt);
	} else {
		uint32_t *sav_fg = tmp + ((size_t)6 << logn_top);
		for (unsigned d = MIN_SAVE_FG[logn_top];
			d <= depth; d ++)
		{
			sav_fg -= MAX_BL_SMALL[d] << (logn_top + 1 - d);
		}
		memmove(fgt, sav_fg, 2 * slen * n * sizeof *fgt);
	}

	/* Move buffers so that we have room for the unreduced (F,G) at
	   this level.
	     Ft   F from this level (unreduced) (llen*n)
	     Gt   G from this level (unreduced) (llen*n)
	     ft   f from this level (slen*n)
	     gt   g from this level (slen*n)
	     Fd   F from deeper level (dlen*hn)
	     Gd   G from deeper level (dlen*hn)  */
	uint32_t *Ft = tmp;
	uint32_t *Gt = Ft + llen * n;
	uint32_t *ft = Gt + llen * n;
	uint32_t *gt = ft + slen * n;
	Fd = gt + slen * n;
	Gd = Fd + dlen * hn;
	uint32_t *t1 = Gd + dlen * hn;
	memmove(ft, fgt, 2 * n * slen * sizeof *ft);
	memmove(Fd, tmp, 2 * hn * dlen * sizeof *tmp);

	/* Convert Fd and Gd to RNS, with output temporarily stored
	   in (Ft, Gt). Fd and Gd have degree hn only; we store the
	   values for each modulus p in the _last_ hn slots of the
	   n-word line for that modulus. */
	for (size_t i = 0; i < llen; i ++) {
		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;
		uint32_t Rx = mp_Rx31((unsigned)dlen, p, p0i, R2);
		uint32_t *xt = Ft + i * n + hn;
		uint32_t *yt = Gt + i * n + hn;
		if (logn >= 4) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			__m256i yRx = _mm256_set1_epi32(Rx);
			for (size_t j = 0; j < hn; j += 8) {
				_mm256_storeu_si256((__m256i *)(xt + j),
					zint_mod_small_signed_x8(Fd + j, dlen,
						hn, yp, yp0i, yR2, yRx));
				_mm256_storeu_si256((__m256i *)(yt + j),
					zint_mod_small_signed_x8(Gd + j, dlen,
						hn, yp, yp0i, yR2, yRx));
			}
		} else {
			for (size_t j = 0; j < hn; j ++) {
				xt[j] = zint_mod_small_signed(Fd + j, dlen, hn,
					p, p0i, R2, Rx);
				yt[j] = zint_mod_small_signed(Gd + j, dlen, hn,
					p, p0i, R2, Rx);
			}
		}
	}

	/* Fd and Gd are no longer needed. */
	t1 = Fd;

	/* Compute (F,G) (unreduced) modulo sufficiently many small primes.
	   We also un-NTT (f,g) as we go; when slen primes have been
	   processed, we obtain (f,g) in RNS, and we apply the CRT to
	   get (f,g) in plain representation. */
	for (size_t i = 0; i < llen; i ++) {
		/* If we have processed exactly slen primes, then (f,g)
		   are in RNS, and we can rebuild them. */
		if (i == slen) {
			avx2_zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
		}

		uint32_t p = PRIMES[i].p;
		uint32_t p0i = PRIMES[i].p0i;
		uint32_t R2 = PRIMES[i].R2;

		/* Memory layout: we keep Ft, Gt, ft and gt; we append:
		     gm    NTT support (n)
		     igm   iNTT support (n)
		     fx    temporary f mod p (NTT) (n)
		     gx    temporary g mod p (NTT) (n)  */
		uint32_t *gm = t1;
		uint32_t *igm = gm + n;
		uint32_t *fx = igm + n;
		uint32_t *gx = fx + n;
		avx2_mp_mkgmigm(logn, gm, igm,
			PRIMES[i].g, PRIMES[i].ig, p, p0i);
		if (i < slen) {
			memcpy(fx, ft + i * n, n * sizeof *fx);
			memcpy(gx, gt + i * n, n * sizeof *gx);
			avx2_mp_iNTT(logn, ft + i * n, igm, p, p0i);
			avx2_mp_iNTT(logn, gt + i * n, igm, p, p0i);
		} else {
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			for (size_t j = 0; j < n; j ++) {
				fx[j] = zint_mod_small_signed(ft + j, slen, n,
					p, p0i, R2, Rx);
				gx[j] = zint_mod_small_signed(gt + j, slen, n,
					p, p0i, R2, Rx);
			}
			avx2_mp_NTT(logn, fx, gm, p, p0i);
			avx2_mp_NTT(logn, gx, gm, p, p0i);
		}

		/* We have (F,G) from deeper level in Ft and Gt, in
		   RNS. We apply the NTT modulo p. */
		uint32_t *Fe = Ft + i * n;
		uint32_t *Ge = Gt + i * n;
		avx2_mp_NTT(logn - 1, Fe + hn, gm, p, p0i);
		avx2_mp_NTT(logn - 1, Ge + hn, gm, p, p0i);

		/* Compute F and G (unreduced) modulo p. */
		if (hn >= 4) {
			__m256i yp = _mm256_set1_epi32(p);
			__m256i yp0i = _mm256_set1_epi32(p0i);
			__m256i yR2 = _mm256_set1_epi32(R2);
			for (size_t j = 0; j < hn; j += 4) {
				__m256i yfa = _mm256_loadu_si256(
					(__m256i *)(fx + (j << 1)));
				__m256i yga = _mm256_loadu_si256(
					(__m256i *)(gx + (j << 1)));
				__m256i yfb = _mm256_srli_epi64(yfa, 32);
				__m256i ygb = _mm256_srli_epi64(yga, 32);
				__m128i xFe = _mm_loadu_si128(
					(__m128i *)(Fe + j + hn));
				__m128i xGe = _mm_loadu_si128(
					(__m128i *)(Ge + j + hn));
				__m256i yFp = _mm256_permute4x64_epi64(
					_mm256_castsi128_si256(xFe), 0x50);
				__m256i yGp = _mm256_permute4x64_epi64(
					_mm256_castsi128_si256(xGe), 0x50);
				yFp = _mm256_shuffle_epi32(yFp, 0x30);
				yGp = _mm256_shuffle_epi32(yGp, 0x30);
				yFp = mp_mmul_x4(yFp, yR2, yp, yp0i);
				yGp = mp_mmul_x4(yGp, yR2, yp, yp0i);
				__m256i yFe0 = mp_mmul_x4(
					ygb, yFp, yp, yp0i);
				__m256i yFe1 = mp_mmul_x4(
					yga, yFp, yp, yp0i);
				__m256i yGe0 = mp_mmul_x4(
					yfb, yGp, yp, yp0i);
				__m256i yGe1 = mp_mmul_x4(
					yfa, yGp, yp, yp0i);
				_mm256_storeu_si256((__m256i *)(Fe + (j << 1)),
					_mm256_or_si256(yFe0,
						_mm256_slli_epi64(yFe1, 32)));
				_mm256_storeu_si256((__m256i *)(Ge + (j << 1)),
					_mm256_or_si256(yGe0,
						_mm256_slli_epi64(yGe1, 32)));
			}
		} else {
			for (size_t j = 0; j < hn; j ++) {
				uint32_t fa = fx[(j << 1) + 0];
				uint32_t fb = fx[(j << 1) + 1];
				uint32_t ga = gx[(j << 1) + 0];
				uint32_t gb = gx[(j << 1) + 1];
				uint32_t mFp = mp_mmul(
					Fe[j + hn], R2, p, p0i);
				uint32_t mGp = mp_mmul(
					Ge[j + hn], R2, p, p0i);
				Fe[(j << 1) + 0] = mp_mmul(gb, mFp, p, p0i);
				Fe[(j << 1) + 1] = mp_mmul(ga, mFp, p, p0i);
				Ge[(j << 1) + 0] = mp_mmul(fb, mGp, p, p0i);
				Ge[(j << 1) + 1] = mp_mmul(fa, mGp, p, p0i);
			}
		}

		/* We want the new (F,G) in RNS only (no NTT). */
		avx2_mp_iNTT(logn, Fe, igm, p, p0i);
		avx2_mp_iNTT(logn, Ge, igm, p, p0i);
	}

	/* Edge case: if slen == llen, then we have not rebuilt (f,g)
	   into plain representation yet, so we do it now. */
	if (slen == llen) {
		avx2_zint_rebuild_CRT(ft, slen, n, 2, 1, t1);
	}

	/* We now have the unreduced (F,G) in RNS. We rebuild their
	   plain representation. */
	avx2_zint_rebuild_CRT(Ft, llen, n, 2, 1, t1);

	/* We now reduce these (F,G) with Babai's nearest plane
	   algorithm. The reduction conceptually goes as follows:
	     k <- round((F*adj(f) + G*adj(g))/(f*adj(f) + g*adj(g)))
	     (F, G) <- (F - k*f, G - k*g)
	   We use fixed-point approximations of (f,g) and (F, G) to get
	   a value k as a small polynomial with scaling; we then apply
	   k on the full-width polynomial. Each iteration "shaves" a
	   a few bits off F and G.
	  
	   We apply the process sufficiently many times to reduce (F, G)
	   to the size of (f, g) with a reasonable probability of success.
	   Since we want full constant-time processing, the number of
	   iterations and the accessed slots work on some assumptions on
	   the sizes of values (sizes have been measured over many samples,
	   and a margin of 5 times the standard deviation). */

	/* If depth is at least 2, and we will use the NTT to subtract
	   k*(f,g) from (F,G), then we will need to convert (f,g) to NTT over
	   slen+1 words, which requires an extra word to ft and gt. */
	int use_sub_ntt = (depth > 1 && logn >= MIN_LOGN_FGNTT);
	if (use_sub_ntt) {
		memmove(gt + n, gt, n * slen * sizeof *gt);
		gt += n;
		t1 += 2 * n;
	}

	/* New layout:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     ft    f from this level (slen*n) (+n if use_sub_ntt)
	     gt    g from this level (slen*n) (+n if use_sub_ntt)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (hn fxr = n)  */

	fxr *rt3 = (fxr *)t1;
	fxr *rt4 = rt3 + n;
	fxr *rt1 = rt4 + n;

	/* We consider only the top rlen words of (f,g). */
	size_t rlen = WORD_WIN[depth];
	if (rlen > slen) {
		rlen = slen;
	}
	size_t blen = slen - rlen;
	uint32_t *ftb = ft + blen * n;
	uint32_t *gtb = gt + blen * n;
	uint32_t scale_fg = 31 * (uint32_t)blen;
	uint32_t scale_FG = 31 * (uint32_t)llen;

	/* Convert f and g into fixed-point approximations, in rt3 and rt4,
	   respectively. They are scaled down by 2^(scale_fg + scale_x).
	   scale_fg is public (it depends only on the recursion depth), but
	   scale_x comes from a measurement on the actual values of (f,g) and
	   is thus secret.
	  
	   The value scale_x is adjusted so that the largest coefficient is
	   close to, but lower than, some limit t (in absolute value). The
	   limit t is chosen so that f*adj(f) + g*adj(g) does not overflow,
	   i.e. all coefficients must remain below 2^31.
	  
	   Let n be the degree; we know that n <= 2^10. The squared norm
	   of a polynomial is the sum of the squared norms of the
	   coefficients, with the squared norm of a complex number being
	   the product of that number with its complex conjugate. If all
	   coefficients of f are less than t (in absolute value), then
	   the squared norm of f is less than n*t^2. The squared norm of
	   FFT(f) (f in FFT representation) is exactly n times the
	   squared norm of f, so this leads to n^2*t^2 as a maximum
	   bound. adj(f) has the same norm as f. This implies that each
	   complex coefficient of FFT(f) has a maximum squared norm of
	   n^2*t^2 (with a maximally imbalanced polynomial with all
	   coefficient but one being zero). The computation of f*adj(f)
	   exactly is, in FFT representation, the product of each
	   coefficient with its conjugate; thus, the coefficients of
	   f*adj(f), in FFT representation, are at most n^2*t^2.
	  
	   Since we want the coefficients of f*adj(f)+g*adj(g) not to exceed
	   2^31, we need n^2*t^2 <= 2^30, i.e. n*t <= 2^15. We can adjust t
	   accordingly (called scale_t in the code below). We also need to
	   take care that t must not exceed scale_x. Approximation of f and
	   g are extracted with scale scale_fg + scale_x - scale_t, and
	   later fixed by dividing them by 2^scale_t. */
	uint32_t scale_xf = poly_max_bitlength(logn, ftb, rlen);
	uint32_t scale_xg = poly_max_bitlength(logn, gtb, rlen);
	uint32_t scale_x = scale_xf;
	scale_x ^= (scale_xf ^ scale_xg) & tbmask(scale_xf - scale_xg);
	uint32_t scale_t = 15 - logn;
	scale_t ^= (scale_t ^ scale_x) & tbmask(scale_x - scale_t);
	uint32_t scdiff = scale_x - scale_t;

	poly_big_to_fixed(logn, rt3, ftb, rlen, scdiff);
	poly_big_to_fixed(logn, rt4, gtb, rlen, scdiff);

	/* rt3 <- adj(f)/(f*adj(f) + g*adj(g))  (FFT)
	   rt4 <- adj(g)/(f*adj(f) + g*adj(g))  (FFT) */
	avx2_vect_FFT(logn, rt3);
	avx2_vect_FFT(logn, rt4);
	avx2_vect_norm_fft(logn, rt1, rt3, rt4);
	avx2_vect_mul2e(logn, rt3, scale_t);
	avx2_vect_mul2e(logn, rt4, scale_t);
	for (size_t i = 0; i < hn; i ++) {
		fxr ni3 = fxr_neg(rt3[i + hn]);
		fxr ni4 = fxr_neg(rt4[i + hn]);
		fxr_div_x4_1(&rt3[i], &ni3, &rt4[i], &ni4, rt1[i]);
		rt3[i + hn] = ni3;
		rt4[i + hn] = ni4;
	}

	/* New layout:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     ft    f from this level (slen*n) (+n if use_sub_ntt)
	     gt    g from this level (slen*n) (+n if use_sub_ntt)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (n fxr = 2*n)     |   k    (n)
	     rt2   (n fxr = 2*n)     |   t2   (3*n)
	   Exception: at depth == 1, we omit ft and gt:
	     Ft    F from this level (unreduced) (llen*n)
	     Gt    G from this level (unreduced) (llen*n)
	     rt3   (n fxr = 2*n)
	     rt4   (n fxr = 2*n)
	     rt1   (n fxr = 2*n)     |   k    (n)
	     rt2   (n fxr = 2*n)     |   t2   (3*n)  */
	if (depth == 1) {
		t1 = ft;
		fxr *nrt3 = (fxr *)t1;
		memmove(nrt3, rt3, 2 * n * sizeof *rt3);
		rt3 = nrt3;
		rt4 = rt3 + n;
		rt1 = rt4 + n;
	}
	int32_t *k = (int32_t *)rt1;
	uint32_t *t2 = (uint32_t *)(k + n);
	fxr *rt2 = (fxr *)t2;
	if (rt2 < (rt1 + n)) {
		rt2 = rt1 + n;
	}

	/* If we are going to use poly_sub_scaled_ntt(), then we convert
	   f and g to the NTT representation. Since poly_sub_scaled_ntt()
	   itself will use more than n*(slen+2) words in t2[], we can do
	   the same here. */
	if (use_sub_ntt) {
		uint32_t *gm = t2;
		uint32_t *tn = gm + n;
		for (size_t i = 0; i <= slen; i ++) {
			uint32_t p = PRIMES[i].p;
			uint32_t p0i = PRIMES[i].p0i;
			uint32_t R2 = PRIMES[i].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			avx2_mp_mkgm(logn, gm, PRIMES[i].g, p, p0i);
			for (size_t j = 0; j < n; j ++) {
				tn[j] = zint_mod_small_signed(
					ft + j, slen, n, p, p0i, R2, Rx);
			}
			avx2_mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(ft, tn, (slen + 1) * n * sizeof *tn);
		for (size_t i = 0; i <= slen; i ++) {
			uint32_t p = PRIMES[i].p;
			uint32_t p0i = PRIMES[i].p0i;
			uint32_t R2 = PRIMES[i].R2;
			uint32_t Rx = mp_Rx31((unsigned)slen, p, p0i, R2);
			avx2_mp_mkgm(logn, gm, PRIMES[i].g, p, p0i);
			for (size_t j = 0; j < n; j ++) {
				tn[j] = zint_mod_small_signed(
					gt + j, slen, n, p, p0i, R2, Rx);
			}
			avx2_mp_NTT(logn, tn, gm, p, p0i);
			tn += n;
		}
		tn = gm + n;
		memmove(gt, tn, (slen + 1) * n * sizeof *tn);
	}

	/* Reduce F and G repeatedly. */
	size_t FGlen = llen;
	uint32_t reduce_bits;
	switch (logn_top) {
	case 9:  reduce_bits = 13; break;
	case 10: reduce_bits = 11; break;
	default: reduce_bits = 16; break;
	}
	for (;;) {
		/* Convert the current F and G into fixed-point. We want
		   to apply scaling scale_FG + scale_x. */
		uint32_t tlen, toff;
		DIVREM31(tlen, toff, scale_FG);
		poly_big_to_fixed(logn, rt1,
			Ft + tlen * n, FGlen - tlen, scale_x + toff);
		poly_big_to_fixed(logn, rt2,
			Gt + tlen * n, FGlen - tlen, scale_x + toff);

		/* rt2 <- (F*adj(f) + G*adj(g)) / (f*adj(f) + g*adj(g)) */
		avx2_vect_FFT(logn, rt1);
		avx2_vect_FFT(logn, rt2);
		avx2_vect_mul_fft(logn, rt1, rt3);
		avx2_vect_mul_fft(logn, rt2, rt4);
		avx2_vect_add(logn, rt2, rt1);
		avx2_vect_iFFT(logn, rt2);

		/* k <- round(rt2) */
		for (size_t i = 0; i < n; i ++) {
			k[i] = fxr_round(rt2[i]);
		}

		/* (f,g) are scaled by scale_fg + scale_x
		   (F,G) are scaled by scale_FG + scale_x
		   Thus, k is scaled by scale_FG - scale_fg, which is public. */
		uint32_t scale_k = scale_FG - scale_fg;
		if (depth == 1) {
			avx2_poly_sub_kfg_scaled_depth1(logn_top,
				Ft, Gt, FGlen,
				(uint32_t *)k, scale_k, f, g, t2);
		} else if (use_sub_ntt) {
			avx2_poly_sub_scaled_ntt(logn,
				Ft, FGlen, ft, slen, k, scale_k, t2);
			avx2_poly_sub_scaled_ntt(logn,
				Gt, FGlen, gt, slen, k, scale_k, t2);
		} else {
			avx2_poly_sub_scaled(logn,
				Ft, FGlen, ft, slen, k, scale_k);
			avx2_poly_sub_scaled(logn,
				Gt, FGlen, gt, slen, k, scale_k);
		}

		/* We now assume that F and G have shrunk by at least
		   reduce_bits. We adjust FGlen accordinly. */
		if (scale_FG <= scale_fg) {
			break;
		}
		if (scale_FG <= (scale_fg + reduce_bits)) {
			scale_FG = scale_fg;
		} else {
			scale_FG -= reduce_bits;
		}
		while (FGlen > slen
			&& 31 * (FGlen - slen) > scale_FG - scale_fg + 30)
		{
			FGlen --;
		}
	}

	/* Output F is already in the right place; G is in Gt, and must be
	   moved back a bit. */
	memmove(tmp + slen * n, Gt, slen * n * sizeof *tmp);
	Gt = tmp + slen * n;

	/* Reduction is done. We test the current solution modulo a single
	   prime.
	   Exception: we cannot do that if depth == 1, since in that case
	   we did not keep (ft,gt). Reduction errors rarely occur at this
	   stage, so we can omit that test (depth-0 test will cover it).
	  
	   If use_sub_ntt != 0, then ft and gt are already in NTT
	   representation. */
	if (depth == 1) {
		return SOLVE_OK;
	}

	t2 = t1 + n;
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;
	uint32_t Rx = mp_Rx31(slen, p, p0i, R2);
	avx2_mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	if (use_sub_ntt) {
		t1 = ft;
		for (size_t i = 0; i < n; i ++) {
			t2[i] = zint_mod_small_signed(
				Gt + i, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t i = 0; i < n; i ++) {
			t1[i] = zint_mod_small_signed(
				ft + i, slen, n, p, p0i, R2, Rx);
			t2[i] = zint_mod_small_signed(
				Gt + i, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t1, t4, p, p0i);
		avx2_mp_NTT(logn, t2, t4, p, p0i);
	}
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		for (size_t i = 0; i < n; i += 8) {
			__m256i y1 = _mm256_loadu_si256((__m256i *)(t1 + i));
			__m256i y2 = _mm256_loadu_si256((__m256i *)(t2 + i));
			__m256i y3 = mp_mmul_x8(y1, y2, yp, yp0i);
			_mm256_storeu_si256((__m256i *)(t3 + i), y3);
		}
	} else {
		for (size_t i = 0; i < n; i ++) {
			t3[i] = mp_mmul(t1[i], t2[i], p, p0i);
		}
	}
	if (use_sub_ntt) {
		t1 = gt;
		for (size_t i = 0; i < n; i ++) {
			t2[i] = zint_mod_small_signed(
				Ft + i, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t2, t4, p, p0i);
	} else {
		for (size_t i = 0; i < n; i ++) {
			t1[i] = zint_mod_small_signed(
				gt + i, slen, n, p, p0i, R2, Rx);
			t2[i] = zint_mod_small_signed(
				Ft + i, slen, n, p, p0i, R2, Rx);
		}
		avx2_mp_NTT(logn, t1, t4, p, p0i);
		avx2_mp_NTT(logn, t2, t4, p, p0i);
	}
	uint32_t rv = mp_mmul(Q, 1, p, p0i);
	if (n >= 8) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yrv = _mm256_set1_epi32(rv);
		for (size_t i = 0; i < n; i += 8) {
			__m256i y1 = _mm256_loadu_si256((__m256i *)(t1 + i));
			__m256i y2 = _mm256_loadu_si256((__m256i *)(t2 + i));
			__m256i y3 = _mm256_loadu_si256((__m256i *)(t3 + i));
			__m256i yx = mp_sub_x8(y3,
				mp_mmul_x8(y1, y2, yp, yp0i), yp);
			if ((uint32_t)_mm256_movemask_epi8(
				_mm256_cmpeq_epi32(yx, yrv)) != 0xFFFFFFFF)
			{
				return SOLVE_ERR_REDUCE;
			}
		}
	} else {
		for (size_t i = 0; i < n; i ++) {
			uint32_t x = mp_mmul(t1[i], t2[i], p, p0i);
			if (mp_sub(t3[i], x, p) != rv) {
				return SOLVE_ERR_REDUCE;
			}
		}
	}

	return SOLVE_OK;
}
#endif

/* Solving the NTRU equation, top recursion level. This is a specialized
   variant for solve_NTRU_intermediate() with depth == 0, for lower RAM
   usage and faster operation.
  
   Returned value: 0 on success, a negative error code otherwise. */
static int
solve_NTRU_depth0(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* At depth 0, all values fit on 30 bits, so we work with a
	   single modulus p. */
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;

	/* Buffer layout:
	     Fd   F from upper level (hn)
	     Gd   G from upper level (hn)
	     ft   f (n)
	     gt   g (n)
	     gm   helper for NTT  */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + hn;
	uint32_t *ft = Gd + hn;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;

	/* Load f and g, and convert to RNS+NTT. */
	mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	poly_mp_set_small(logn, ft, f, p);
	poly_mp_set_small(logn, gt, g, p);
	mp_NTT(logn, ft, gm, p, p0i);
	mp_NTT(logn, gt, gm, p, p0i);

	/* Convert Fd and Gd to RNS+NTT. */
	poly_mp_set(logn - 1, Fd, p);
	poly_mp_set(logn - 1, Gd, p);
	mp_NTT(logn - 1, Fd, gm, p, p0i);
	mp_NTT(logn - 1, Gd, gm, p, p0i);

	/* Build the unreduced (F,G) into ft and gt. */
	for (size_t i = 0; i < hn; i ++) {
		uint32_t fa = ft[(i << 1) + 0];
		uint32_t fb = ft[(i << 1) + 1];
		uint32_t ga = gt[(i << 1) + 0];
		uint32_t gb = gt[(i << 1) + 1];
		uint32_t mFd = mp_mmul(Fd[i], R2, p, p0i);
		uint32_t mGd = mp_mmul(Gd[i], R2, p, p0i);
		ft[(i << 1) + 0] = mp_mmul(gb, mFd, p, p0i);
		ft[(i << 1) + 1] = mp_mmul(ga, mFd, p, p0i);
		gt[(i << 1) + 0] = mp_mmul(fb, mGd, p, p0i);
		gt[(i << 1) + 1] = mp_mmul(fa, mGd, p, p0i);
	}

	/* Reorganize buffers:
	     Fp   unreduced F (RNS+NTT) (n)
	     Gp   unreduced G (RNS+NTT) (n)
	     t1   free (n)
	     t2   NTT support (gm) (n)
	     t3   free (n)
	     t4   free (n)  */
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + n;
	uint32_t *t1 = Gp + n;
	uint32_t *t2 = t1 + n;  /* alias on gm */
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	memmove(Fp, ft, 2 * n * sizeof *ft);

	/* Working modulo p (using the NTT), we compute:
	      t1 <- F*adj(f) + G*adj(g)
	      t2 <- f*adj(f) + g*adj(g)  */

	/* t4 <- f (RNS+NTT) */
	poly_mp_set_small(logn, t4, f, p);
	mp_NTT(logn, t4, gm, p, p0i);

	/* t1 <- F*adj(f) (RNS+NTT)
	   t3 <- f*adj(f) (RNS+NTT)  */
	for (size_t i = 0; i < n; i ++) {
		uint32_t w = mp_mmul(t4[(n - 1) - i], R2, p, p0i);
		t1[i] = mp_mmul(w, Fp[i], p, p0i);
		t3[i] = mp_mmul(w, t4[i], p, p0i);
	}

	/* t4 <- g (RNS+NTT) */
	poly_mp_set_small(logn, t4, g, p);
	mp_NTT(logn, t4, gm, p, p0i);

	/* t1 <- t1 + G*adj(g)
	   t3 <- t3 + g*adj(g)  */
	for (size_t i = 0; i < n; i ++) {
		uint32_t w = mp_mmul(t4[(n - 1) - i], R2, p, p0i);
		t1[i] = mp_add(t1[i], mp_mmul(w, Gp[i], p, p0i), p);
		t3[i] = mp_add(t3[i], mp_mmul(w, t4[i], p, p0i), p);
	}

	/* Convert back F*adj(f) + G*adj(g) and f*adj(f) + g*adj(g) to
	   plain representation, and move f*adj(f) + g*adj(g) to t2. */
	mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, t1, t4, p, p0i);
	mp_iNTT(logn, t3, t4, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		/* NOTE: no truncature to 31 bits. */
		t1[i] = (uint32_t)mp_norm(t1[i], p);
		t2[i] = (uint32_t)mp_norm(t3[i], p);
	}

	/* Buffer contents:
	     Fp   unreduced F (RNS+NTT) (n)
	     Gp   unreduced G (RNS+NTT) (n)
	     t1   F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	     t2   f*adj(f) + g*adj(g) (plain, 32-bit) (n)  */

	/* We need to divide t1 by t2, and round the result. We convert
	   them to FFT representation, downscaled by 2^10 (to avoid overflows).
	   We first convert f*adj(f) + g*adj(g), which is self-adjoint;
	   thus, its FFT representation only has half-size. */
	fxr *rt3 = (fxr *)t3;
	for (size_t i = 0; i < n; i ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t2[i] << 22;
		rt3[i] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt3);
	fxr *rt2 = (fxr *)t2;
	memmove(rt2, rt3, hn * sizeof *rt3);
	rt3 = rt2 + hn;

	/* Buffer contents:
	     Fp    unreduced F (RNS+NTT) (n)
	     Gp    unreduced G (RNS+NTT) (n)
	     t1    F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	     rt2   f*adj(f) + g*adj(g) (FFT, self-ajdoint) (hn fxr values = n)
	     rt3   free (n fxr values = 2*n)  */

	/* Convert F*adj(f) + G*adj(g) to FFT (scaled by 2^10) (into rt3). */
	for (size_t i = 0; i < n; i ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t1[i] << 22;
		rt3[i] = fxr_of_scaled32(x);
	}
	vect_FFT(logn, rt3);

	/* Divide F*adj(f) + G*adj(g) by f*adj(f) + g*adj(g) and round
	   the result into t1, with conversion to RNS. */
	vect_div_selfadj_fft(logn, rt3, rt2);
	vect_iFFT(logn, rt3);
	for (size_t i = 0; i < n; i ++) {
		t1[i] = mp_set(fxr_round(rt3[i]), p);
	}

	/* Buffer contents:
	     Fp    unreduced F (RNS+NTT) (n)
	     Gp    unreduced G (RNS+NTT) (n)
	     t1    k (RNS) (n)
	     t2    free (n)
	     t3    free (n)
	     t4    free (n)  */

	/* Convert k to RNS+NTT+Montgomery. */
	mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	mp_NTT(logn, t1, t4, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		t1[i] = mp_mmul(t1[i], R2, p, p0i);
	}

	/* Subtract k*f from F and k*g from G.
	   We also compute f*G - g*F (in RNS+NTT) to check that the solution
	   is correct. */
	for (size_t i = 0; i < n; i ++) {
		t2[i] = mp_set(f[i], p);
		t3[i] = mp_set(g[i], p);
	}
	mp_NTT(logn, t2, t4, p, p0i);
	mp_NTT(logn, t3, t4, p, p0i);
	uint32_t rv = mp_mmul(Q, 1, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		Fp[i] = mp_sub(Fp[i], mp_mmul(t1[i], t2[i], p, p0i), p);
		Gp[i] = mp_sub(Gp[i], mp_mmul(t1[i], t3[i], p, p0i), p);
		uint32_t x = mp_sub(
			mp_mmul(t2[i], Gp[i], p, p0i),
			mp_mmul(t3[i], Fp[i], p, p0i), p);
		if (x != rv) {
			return SOLVE_ERR_REDUCE;
		}
	}

	/* Convert back F and G into normal representation. */
	mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	mp_iNTT(logn, Fp, t4, p, p0i);
	mp_iNTT(logn, Gp, t4, p, p0i);
	poly_mp_norm(logn, Fp, p);
	poly_mp_norm(logn, Gp, p);

	return SOLVE_OK;
}

#if FNDSA_AVX2
TARGET_AVX2
static int
avx2_solve_NTRU_depth0(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g,
	uint32_t *restrict tmp)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;

	/* At depth 0, all values fit on 30 bits, so we work with a
	   single modulus p. */
	uint32_t p = PRIMES[0].p;
	uint32_t p0i = PRIMES[0].p0i;
	uint32_t R2 = PRIMES[0].R2;

	/* Buffer layout:
	     Fd   F from upper level (hn)
	     Gd   G from upper level (hn)
	     ft   f (n)
	     gt   g (n)
	     gm   helper for NTT  */
	uint32_t *Fd = tmp;
	uint32_t *Gd = Fd + hn;
	uint32_t *ft = Gd + hn;
	uint32_t *gt = ft + n;
	uint32_t *gm = gt + n;

	/* Load f and g, and convert to RNS+NTT. */
	avx2_mp_mkgm(logn, gm, PRIMES[0].g, p, p0i);
	avx2_poly_mp_set_small(logn, ft, f, p);
	avx2_poly_mp_set_small(logn, gt, g, p);
	avx2_mp_NTT(logn, ft, gm, p, p0i);
	avx2_mp_NTT(logn, gt, gm, p, p0i);

	/* Convert Fd and Gd to RNS+NTT. */
	avx2_poly_mp_set(logn - 1, Fd, p);
	avx2_poly_mp_set(logn - 1, Gd, p);
	avx2_mp_NTT(logn - 1, Fd, gm, p, p0i);
	avx2_mp_NTT(logn - 1, Gd, gm, p, p0i);

	/* Build the unreduced (F,G) into ft and gt. */
	if (hn >= 4) {
		__m256i yp = _mm256_set1_epi32(p);
		__m256i yp0i = _mm256_set1_epi32(p0i);
		__m256i yR2 = _mm256_set1_epi32(R2);
		for (size_t j = 0; j < hn; j += 4) {
			__m256i yfa = _mm256_loadu_si256(
				(__m256i *)(ft + (j << 1)));
			__m256i yga = _mm256_loadu_si256(
				(__m256i *)(gt + (j << 1)));
			__m256i yfb = _mm256_srli_epi64(yfa, 32);
			__m256i ygb = _mm256_srli_epi64(yga, 32);
			__m128i xFd = _mm_loadu_si128((__m128i *)(Fd + j));
			__m128i xGd = _mm_loadu_si128((__m128i *)(Gd + j));
			__m256i yFd = _mm256_permute4x64_epi64(
				_mm256_castsi128_si256(xFd), 0x50);
			__m256i yGd = _mm256_permute4x64_epi64(
				_mm256_castsi128_si256(xGd), 0x50);
			yFd = _mm256_shuffle_epi32(yFd, 0x30);
			yGd = _mm256_shuffle_epi32(yGd, 0x30);
			yFd = mp_mmul_x4(yFd, yR2, yp, yp0i);
			yGd = mp_mmul_x4(yGd, yR2, yp, yp0i);
			__m256i yFe0 = mp_mmul_x4(ygb, yFd, yp, yp0i);
			__m256i yFe1 = mp_mmul_x4(yga, yFd, yp, yp0i);
			__m256i yGe0 = mp_mmul_x4(yfb, yGd, yp, yp0i);
			__m256i yGe1 = mp_mmul_x4(yfa, yGd, yp, yp0i);
			_mm256_storeu_si256((__m256i *)(ft + (j << 1)),
				_mm256_or_si256(yFe0,
					_mm256_slli_epi64(yFe1, 32)));
			_mm256_storeu_si256((__m256i *)(gt + (j << 1)),
				_mm256_or_si256(yGe0,
					_mm256_slli_epi64(yGe1, 32)));
		}
	} else {
		for (size_t j = 0; j < hn; j ++) {
			uint32_t fa = ft[(j << 1) + 0];
			uint32_t fb = ft[(j << 1) + 1];
			uint32_t ga = gt[(j << 1) + 0];
			uint32_t gb = gt[(j << 1) + 1];
			uint32_t mFd = mp_mmul(Fd[j], R2, p, p0i);
			uint32_t mGd = mp_mmul(Gd[j], R2, p, p0i);
			ft[(j << 1) + 0] = mp_mmul(gb, mFd, p, p0i);
			ft[(j << 1) + 1] = mp_mmul(ga, mFd, p, p0i);
			gt[(j << 1) + 0] = mp_mmul(fb, mGd, p, p0i);
			gt[(j << 1) + 1] = mp_mmul(fa, mGd, p, p0i);
		}
	}

	/* Reorganize buffers:
	     Fp   unreduced F (RNS+NTT) (n)
	     Gp   unreduced G (RNS+NTT) (n)
	     t1   free (n)
	     t2   NTT support (gm) (n)
	     t3   free (n)
	     t4   free (n)  */
	uint32_t *Fp = tmp;
	uint32_t *Gp = Fp + n;
	uint32_t *t1 = Gp + n;
	uint32_t *t2 = t1 + n;  /* alias on gm */
	uint32_t *t3 = t2 + n;
	uint32_t *t4 = t3 + n;
	memmove(Fp, ft, 2 * n * sizeof *ft);

	/* Working modulo p (using the NTT), we compute:
	      t1 <- F*adj(f) + G*adj(g)
	      t2 <- f*adj(f) + g*adj(g)  */

	/* t4 <- f (RNS+NTT) */
	avx2_poly_mp_set_small(logn, t4, f, p);
	avx2_mp_NTT(logn, t4, gm, p, p0i);

	/* t1 <- F*adj(f) (RNS+NTT)
	   t3 <- f*adj(f) (RNS+NTT)  */
	for (size_t i = 0; i < n; i ++) {
		uint32_t w = mp_mmul(t4[(n - 1) - i], R2, p, p0i);
		t1[i] = mp_mmul(w, Fp[i], p, p0i);
		t3[i] = mp_mmul(w, t4[i], p, p0i);
	}

	/* t4 <- g (RNS+NTT) */
	avx2_poly_mp_set_small(logn, t4, g, p);
	avx2_mp_NTT(logn, t4, gm, p, p0i);

	/* t1 <- t1 + G*adj(g)
	   t3 <- t3 + g*adj(g)  */
	for (size_t i = 0; i < n; i ++) {
		uint32_t w = mp_mmul(t4[(n - 1) - i], R2, p, p0i);
		t1[i] = mp_add(t1[i], mp_mmul(w, Gp[i], p, p0i), p);
		t3[i] = mp_add(t3[i], mp_mmul(w, t4[i], p, p0i), p);
	}

	/* Convert back F*adj(f) + G*adj(g) and f*adj(f) + g*adj(g) to
	   plain representation, and move f*adj(f) + g*adj(g) to t2. */
	avx2_mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	avx2_mp_iNTT(logn, t1, t4, p, p0i);
	avx2_mp_iNTT(logn, t3, t4, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		/* NOTE: no truncature to 31 bits. */
		t1[i] = (uint32_t)mp_norm(t1[i], p);
		t2[i] = (uint32_t)mp_norm(t3[i], p);
	}

	/* Buffer contents:
	     Fp   unreduced F (RNS+NTT) (n)
	     Gp   unreduced G (RNS+NTT) (n)
	     t1   F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	     t2   f*adj(f) + g*adj(g) (plain, 32-bit) (n)  */

	/* We need to divide t1 by t2, and round the result. We convert
	   them to FFT representation, downscaled by 2^10 (to avoid overflows).
	   We first convert f*adj(f) + g*adj(g), which is self-adjoint;
	   thus, its FFT representation only has half-size. */
	fxr *rt3 = (fxr *)t3;
	for (size_t i = 0; i < n; i ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t2[i] << 22;
		rt3[i] = fxr_of_scaled32(x);
	}
	avx2_vect_FFT(logn, rt3);
	fxr *rt2 = (fxr *)t2;
	memmove(rt2, rt3, hn * sizeof *rt3);
	rt3 = rt2 + hn;

	/* Buffer contents:
	     Fp    unreduced F (RNS+NTT) (n)
	     Gp    unreduced G (RNS+NTT) (n)
	     t1    F*adj(f) + G*adj(g) (plain, 32-bit) (n)
	     rt2   f*adj(f) + g*adj(g) (FFT, self-ajdoint) (hn fxr values = n)
	     rt3   free (n fxr values = 2*n)  */

	/* Convert F*adj(f) + G*adj(g) to FFT (scaled by 2^10) (into rt3). */
	for (size_t i = 0; i < n; i ++) {
		uint64_t x = (uint64_t)*(int32_t *)&t1[i] << 22;
		rt3[i] = fxr_of_scaled32(x);
	}
	avx2_vect_FFT(logn, rt3);

	/* Divide F*adj(f) + G*adj(g) by f*adj(f) + g*adj(g) and round
	   the result into t1, with conversion to RNS. */
	avx2_vect_div_selfadj_fft(logn, rt3, rt2);
	avx2_vect_iFFT(logn, rt3);
	for (size_t i = 0; i < n; i ++) {
		t1[i] = mp_set(fxr_round(rt3[i]), p);
	}

	/* Buffer contents:
	     Fp    unreduced F (RNS+NTT) (n)
	     Gp    unreduced G (RNS+NTT) (n)
	     t1    k (RNS) (n)
	     t2    free (n)
	     t3    free (n)
	     t4    free (n)  */

	/* Convert k to RNS+NTT+Montgomery. */
	avx2_mp_mkgm(logn, t4, PRIMES[0].g, p, p0i);
	avx2_mp_NTT(logn, t1, t4, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		t1[i] = mp_mmul(t1[i], R2, p, p0i);
	}

	/* Subtract k*f from F and k*g from G.
	   We also compute f*G - g*F (in RNS+NTT) to check that the solution
	   is correct. */
	for (size_t i = 0; i < n; i ++) {
		t2[i] = mp_set(f[i], p);
		t3[i] = mp_set(g[i], p);
	}
	avx2_mp_NTT(logn, t2, t4, p, p0i);
	avx2_mp_NTT(logn, t3, t4, p, p0i);
	uint32_t rv = mp_mmul(Q, 1, p, p0i);
	for (size_t i = 0; i < n; i ++) {
		Fp[i] = mp_sub(Fp[i], mp_mmul(t1[i], t2[i], p, p0i), p);
		Gp[i] = mp_sub(Gp[i], mp_mmul(t1[i], t3[i], p, p0i), p);
		uint32_t x = mp_sub(
			mp_mmul(t2[i], Gp[i], p, p0i),
			mp_mmul(t3[i], Fp[i], p, p0i), p);
		if (x != rv) {
			return SOLVE_ERR_REDUCE;
		}
	}

	/* Convert back F and G into normal representation. */
	avx2_mp_mkigm(logn, t4, PRIMES[0].ig, p, p0i);
	avx2_mp_iNTT(logn, Fp, t4, p, p0i);
	avx2_mp_iNTT(logn, Gp, t4, p, p0i);
	avx2_poly_mp_norm(logn, Fp, p);
	avx2_poly_mp_norm(logn, Gp, p);

	return SOLVE_OK;
}
#endif

/* see kgen_inner.h */
int
solve_NTRU(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp)
{
	size_t n = (size_t)1 << logn;

	int err = solve_NTRU_deepest(logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return 0;
	}
	unsigned depth = logn;
	while (depth -- > 1) {
		err = solve_NTRU_intermediate(logn, f, g, depth, tmp);
		if (err != SOLVE_OK) {
			return 0;
		}
	}
	err = solve_NTRU_depth0(logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return 0;
	}

	/* F and G are at the start of tmp[] (plain, 31 bits per value).
	   We need to convert them to 8-bit representation, and check
	   that they are within the expected range. */
	int8_t *F = (int8_t *)(tmp + 2 * n);
	int8_t *G = F + n;
	int lim = 127;
	if (!poly_big_to_small(logn, F, tmp, lim)) {
		/* SOLVE_ERR_LIMIT */
		return 0;
	}
	if (!poly_big_to_small(logn, G, tmp + n, lim)) {
		/* SOLVE_ERR_LIMIT */
		return 0;
	}
	memmove(tmp, F, 2 * n);

	return 1;
}

#if FNDSA_AVX2
TARGET_AVX2
int
avx2_solve_NTRU(unsigned logn,
	const int8_t *restrict f, const int8_t *restrict g, uint32_t *tmp)
{
	size_t n = (size_t)1 << logn;

	int err = avx2_solve_NTRU_deepest(logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return 0;
	}
	unsigned depth = logn;
	while (depth -- > 1) {
		err = avx2_solve_NTRU_intermediate(logn, f, g, depth, tmp);
		if (err != SOLVE_OK) {
			return 0;
		}
	}
	err = avx2_solve_NTRU_depth0(logn, f, g, tmp);
	if (err != SOLVE_OK) {
		return 0;
	}

	/* F and G are at the start of tmp[] (plain, 31 bits per value).
	   We need to convert them to 8-bit representation, and check
	   that they are within the expected range. */
	int8_t *F = (int8_t *)(tmp + 2 * n);
	int8_t *G = F + n;
	int lim = 127;
	if (!poly_big_to_small(logn, F, tmp, lim)) {
		/* SOLVE_ERR_LIMIT */
		return 0;
	}
	if (!poly_big_to_small(logn, G, tmp + n, lim)) {
		/* SOLVE_ERR_LIMIT */
		return 0;
	}
	memmove(tmp, F, 2 * n);

	return 1;
}
#endif

/* see kgen_inner.h */
int
check_ortho_norm(unsigned logn, const int8_t *f, const int8_t *g, fxr *tmp)
{
	size_t n = (size_t)1 << logn;
	fxr *rt1 = tmp;
	fxr *rt2 = rt1 + n;
	fxr *rt3 = rt2 + n;
	vect_set(logn, rt1, f);
	vect_set(logn, rt2, g);
	vect_FFT(logn, rt1);
	vect_FFT(logn, rt2);
	vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
	vect_adj_fft(logn, rt1);
	vect_adj_fft(logn, rt2);
	vect_mul_realconst(logn, rt1, fxr_of(12289));
	vect_mul_realconst(logn, rt2, fxr_of(12289));
	vect_mul_selfadj_fft(logn, rt1, rt3);
	vect_mul_selfadj_fft(logn, rt2, rt3);
	vect_iFFT(logn, rt1);
	vect_iFFT(logn, rt2);
	fxr sn = fxr_zero;
	for (size_t i = 0; i < n; i ++) {
		sn = fxr_add(sn, fxr_add(fxr_sqr(rt1[i]), fxr_sqr(rt2[i])));
	}
	return fxr_lt(sn, fxr_of_scaled32(72251709809335));
}

#if FNDSA_AVX2
TARGET_AVX2
int
avx2_check_ortho_norm(unsigned logn, const int8_t *f, const int8_t *g, fxr *tmp)
{
	size_t n = (size_t)1 << logn;
	fxr *rt1 = tmp;
	fxr *rt2 = rt1 + n;
	fxr *rt3 = rt2 + n;
	avx2_vect_set(logn, rt1, f);
	avx2_vect_set(logn, rt2, g);
	avx2_vect_FFT(logn, rt1);
	avx2_vect_FFT(logn, rt2);
	avx2_vect_invnorm_fft(logn, rt3, rt1, rt2, 0);
	avx2_vect_adj_fft(logn, rt1);
	avx2_vect_adj_fft(logn, rt2);
	avx2_vect_mul_realconst(logn, rt1, fxr_of(12289));
	avx2_vect_mul_realconst(logn, rt2, fxr_of(12289));
	avx2_vect_mul_selfadj_fft(logn, rt1, rt3);
	avx2_vect_mul_selfadj_fft(logn, rt2, rt3);
	avx2_vect_iFFT(logn, rt1);
	avx2_vect_iFFT(logn, rt2);
	fxr sn;
	if (logn >= 2) {
		__m256i ysn = _mm256_setzero_si256();
		const __m256i *rp1 = (const __m256i *)rt1;
		const __m256i *rp2 = (const __m256i *)rt2;
		for (size_t i = 0; i < (n >> 2); i ++) {
			__m256i y1 = _mm256_loadu_si256(rp1 + i);
			__m256i y2 = _mm256_loadu_si256(rp2 + i);
			y1 = fxr_sqr_x4(y1);
			y2 = fxr_sqr_x4(y2);
			ysn = _mm256_add_epi64(ysn, _mm256_add_epi64(y1, y2));
		}
		__m128i xsn = _mm_add_epi64(
			_mm256_castsi256_si128(ysn),
			_mm256_extracti128_si256(ysn, 1));
		xsn = _mm_add_epi64(xsn, _mm_bsrli_si128(xsn, 8));
		uint32_t lo = (uint32_t)_mm_cvtsi128_si32(xsn);
		int32_t hi = _mm_cvtsi128_si32(_mm_bsrli_si128(xsn, 4));
		sn = fxr_of_scaled32(((int64_t)hi << 32) | (int64_t)lo);
	} else {
		/* Unused, kept for reference only: logn is the top-level
		   degree, and we enforce logn >= 2 throughout the
		   implementation. */
		sn = fxr_zero;
		for (size_t i = 0; i < n; i ++) {
			sn = fxr_add(sn,
				fxr_add(fxr_sqr(rt1[i]), fxr_sqr(rt2[i])));
		}
	}
	return fxr_lt(sn, fxr_of_scaled32(72251709809335));
}
#endif
