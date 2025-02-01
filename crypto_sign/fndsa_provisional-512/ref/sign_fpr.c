/*
 * Operations on floating-point values.
 */

#include "sign_inner.h"

#define M52   (((uint64_t)1 << 52) - 1)
#define M63   (((uint64_t)1 << 63) - 1)

/* Make a value out of the sign bit s, exponent e, and mantissa.
   Rules for a non-zero value:
      Only the low bit of s is used (0 or 1), other bits are ignored.
      2^54 <= m < 2^55, and the low bit is sticky.
      Value is (-1)^s * 2^e * m; this function applies proper rounding.
      No exponent overflow occurs.
   If the value is a zero, then it must be specified as:
      m = 0
      e = -1076  */
static inline fpr
make(uint64_t s, int32_t e, uint64_t m)
{
	uint64_t cc = (0xC8u >> ((unsigned)m & 7)) & 1;
	return (s << 63)
		+ ((uint64_t)(uint32_t)(e + 1076) << 52)
		+ (m >> 2) + cc;
}

/* Like make(), but it also sets e to -1076 in case m = 0. m MUST be either
   0 or in the [2^54,2^55-1] range. */
static inline fpr
make_z(uint64_t s, int32_t e, uint64_t m)
{
	uint32_t eu = (uint32_t)(e + 1076) & -(uint32_t)(m >> 54);
	uint64_t cc = (0xC8u >> ((unsigned)m & 7)) & 1;
	return (s << 63) + ((uint64_t)(eu) << 52) + (m >> 2) + cc;
}

/* Count of leading zeros in a 64-bit non-zero word. */
static inline uint32_t
lzcnt64_nonzero(uint64_t x)
{
	uint32_t x0 = (uint32_t)x;
	uint32_t x1 = (uint32_t)(x >> 32);
	uint32_t m = ~tbmask(x1 | -x1);
	x1 |= x0 & m;
	return lzcnt_nonzero(x1) + (m & 32);
}

/* Adjust m and e such that m*2^e is preserved, and m is in [2^63,2^64-1].
   If m is 0 then it remains equal to 0, and e is replaced with e-63. */
#define NORM64(m, e)   do { \
		uint64_t norm64_m = (uint64_t)(m); \
		uint32_t norm64_c = lzcnt64_nonzero(norm64_m | 1); \
		(m) = fpr_ulsh(norm64_m, norm64_c); \
		(e) -= (int32_t)norm64_c; \
	} while (0)

#if !FNDSA_ASM_CORTEXM4
/* see sign_inner.h */
fpr
fpr_scaled(int64_t i, int sc)
{
	/* Get sign and absolute value. */
	uint64_t s = (uint64_t)(i >> 63);
	uint64_t m = ((uint64_t)i ^ s) - s;

	/* For now, suppose that m != 0. We normalize m to [2^63,2^64-1]. */
	NORM64(m, sc);

	/* Divide m by 2^9 to get it in [2^54,2^55-1]. The least significant
	   bit must be sticky. */
	sc += 9;
	m = (m | ((m & 0x1FF) + 0x1FF)) >> 9;

	/* If input was zero then m = 0 at this point; otherwise, m is
	   in [2^54,2^55-1]. We can use make_z() to finish the processing. */
	return make_z(s, sc, m);
}
#endif

#if !FNDSA_SSE2 && !FNDSA_NEON

#if !FNDSA_ASM_CORTEXM4
/* see sign_inner.h */
fpr
fpr_add(fpr x, fpr y)
{
	/* Get both operands as x and y, and such that x has the greater
	   absolute value of the two. If x and y have the same absolute
	   value and different signs, when we want x to be the positive
	   value. This guarantees the following:
	     - Exponent of y is not greater than exponent of x.
	     - Result has the sign of x.
	   The special case for identical absolute values is for adding
	   z with -z for some value z. Indeed, if abs(x) = abs(y), then
	   the following situations may happen:
	      x > 0, y = x    -> result is positive
	      x < 0, y = x    -> result is negative
	      x > 0, y = -x   -> result is +0
	      x < 0, y = -x   -> result is +0   (*)
	      x = +0, y = +0  -> result is +0
	      x = +0, y = -0  -> result is +0
	      x = -0, y = +0  -> result is +0   (*)
	      x = -0, y = -0  -> result is -0
	   Enforcing a swap when absolute values are equal but the sign of
	   x is 1 (negative) avoids the two situations tagged '(*)' above.
	   For all other situations, the result indeed has the sign of x.

	   Note that for positive values, the numerical order of encoded
	   exponent||mantissa values matches the order of the encoded
	   values. */
	uint64_t za = (x & M63) - (y & M63);
	za |= ((za - 1) & x);
	uint64_t sw = (x ^ y) & (uint64_t)(*(int64_t *)&za >> 63);
	x ^= sw;
	y ^= sw;

	/* Extract sign bits, exponents and mantissas. The mantissas are
	   scaled up to [2^55,2^56-1] and the exponent is unbiased. If
	   an operand is 0, then its mantissa is set to 0 at this step,
	   and its unbiased exponent is -1078. */
	uint32_t ex = (uint32_t)(x >> 52);
	uint32_t sx = ex >> 11;
	ex &= 0x7FF;
	uint64_t xu = ((x & M52) << 3) | ((uint64_t)((ex + 0x7FF) >> 11) << 55);
	ex -= 1078;

	uint32_t ey = (uint32_t)(y >> 52);
	uint32_t sy = ey >> 11;
	ey &= 0x7FF;
	uint64_t yu = ((y & M52) << 3) | ((uint64_t)((ey + 0x7FF) >> 11) << 55);
	ey -= 1078;

	/* x has the larger exponent, hence we only need to right-shift y.
	   If the shift count is larger than 59 then we clamp the value
	   to 0. */
	uint32_t n = ex - ey;
	yu &= (uint64_t)((int64_t)(*(int32_t *)&n - 60) >> 16);
	n &= 63;

	/* Right-shift yu by n bits; the lowest bit of yu is sticky. */
	uint64_t m = fpr_ulsh(1, n) - 1;
	yu = fpr_ursh(yu | ((yu & m) + m), n);

	/* Add of subtract the mantissas, depending on the sign bits. */
	uint64_t dm = -(uint64_t)(sx ^ sy);
	uint64_t zu = xu + yu - (dm & (yu << 1));

	/* The result may be smaller than abs(x), or slightly larger, though
	   no more than twice larger. We normalize to [2^63, 2^64-1], then
	   shrink back to [2^54,2^55-1] (with a sticky bit). */
	NORM64(zu, ex);
	zu = (zu | ((zu & 0x1FF) + 0x1FF)) >> 9;
	ex += 9;

	/* Result uses the sign of x. */
	return make_z((uint64_t)sx, *(int32_t *)&ex, zu);
}
#endif

#if !FNDSA_ASM_CORTEXM4
/* see sign_inner.h */
fpr
fpr_mul(fpr x, fpr y)
{
	/* Extract absolute values of mantissas, assuming non-zero
	   operands, and multiply them together. */
	uint64_t xu = (x & M52) | ((uint64_t)1 << 52);
	uint64_t yu = (y & M52) | ((uint64_t)1 << 52);

	/* Compute the product into z0:z1:zu (z0 and z1 have size 25 bits
	   each, zu contains the upper bits. */
	uint32_t x0 = (uint32_t)xu & 0x01FFFFFF;
	uint32_t x1 = (uint32_t)(xu >> 25);
	uint32_t y0 = (uint32_t)yu & 0x01FFFFFF;
	uint32_t y1 = (uint32_t)(yu >> 25);
	uint64_t w = (uint64_t)x0 * (uint64_t)y0;
	uint32_t z0 = (uint32_t)w & 0x01FFFFFF;
	uint32_t z1 = (uint32_t)(w >> 25);
	w = (uint64_t)x0 * (uint64_t)y1;
	z1 += (uint32_t)w & 0x01FFFFFF;
	uint32_t z2 = (uint32_t)(w >> 25);
	w = (uint64_t)x1 * (uint64_t)y0;
	z1 += (uint32_t)w & 0x01FFFFFF;
	z2 += (uint32_t)(w >> 25);
	uint64_t zu = (uint64_t)x1 * (uint64_t)y1;
	z2 += (z1 >> 25);
	z1 &= 0x01FFFFFF;
	zu += z2;

	/* Since both xu and yu are in [2^52,2^53-1], product is in
	   [2^104,2^106-1]. We scale down the value into [2^54,2^56-1]
	   with a sticky bit; in practice, this just means dropping
	   z0 and z1, keeping only zu, but setting the lsb of zu to 1
	   if either of z0 or z1 is non-zero. */
	zu |= (uint64_t)(((z0 | z1) + 0x01FFFFFF) >> 25);

	/* If zu is in [2^55,2^56-1] then right-shift it by 1 bit.
	   lsb is sticky and must be preserved if non-zero. */
	uint64_t es = zu >> 55;
	zu = (zu >> es) | (zu & 1);

	/* Aggregate scaling factor:
	    - Each source exponent is biased by 1023.
	    - Integral mantiassas are scaled by 2^52, hence an extra 52
	      bias for each exponent.
	    - However, we right-shifted z by 50 + es.
	   In total: we add exponents, then subtract 2*(1023 + 52),
	   then add 50 + es. */
	uint32_t ex = (uint32_t)(x >> 52) & 0x7FF;
	uint32_t ey = (uint32_t)(y >> 52) & 0x7FF;
	uint32_t e = ex + ey - 2100 + (uint32_t)es;

	/* Sign bit is the XOR of the operand sign bits. */
	uint64_t s = (x ^ y) >> 63;

	/* Corrective action for zeros: if either of the operands is zero,
	   then the computations above are wrong and we must clear the
	   mantissa and adjust the exponent. */
	uint32_t dzu = tbmask((ex - 1) | (ey - 1));
	e ^= dzu & (e ^ (uint32_t)-1076);
	zu &= (uint64_t)(dzu & 1) - 1;
	return make(s, *(int32_t *)&e, zu);
}
#endif

#if !FNDSA_ASM_CORTEXM4
/* see sign_inner.h */
fpr
fpr_div(fpr x, fpr y)
{
	/* Extract mantissas (unsigned). */
	uint64_t xu = (x & M52) | ((uint64_t)1 << 52);
	uint64_t yu = (y & M52) | ((uint64_t)1 << 52);

	/* Perform bit-by-bit division of xu by yu; we run it for 55 bits. */
	uint64_t q = 0;
	for (int i = 0; i < 55; i ++) {
		uint64_t b = ((xu - yu) >> 63) - 1;
		xu -= b & yu;
		q |= b & 1;
		xu <<= 1;
		q <<= 1;
	}

	/* 55-bit quotient is in q, with an extra multiplication by 2.
	   Set the lsb to 1 if xu is non-zero at this point (sticky bit). */
	q |= (xu | -xu) >> 63;

	/* Quotient is at most 2^56-1, but cannot be lower than 2^54 since
	   both operands to the loop were in [2^52,2^53-1]. This is
	   similar to the situation in fpr_mul(). */
	uint64_t es = q >> 55;
	q = (q >> es) | (q & 1);

	/* Aggregate scaling factor. */
	uint32_t ex = (uint32_t)(x >> 52) & 0x7FF;
	uint32_t ey = (uint32_t)(y >> 52) & 0x7FF;
	uint32_t e = ex - ey - 55 + (uint32_t)es;

	/* Sign bit is the XOR of the operand sign bits. */
	uint64_t s = (x ^ y) >> 63;

	/* Corrective action for zeros: if x was zero, then the
	   computations above are wrong and we must clear the mantissa
	   and adjust the exponent. Since we do not support infinites,
	   we assume that y (divisor) was not zero. */
	uint32_t dzu = tbmask(ex - 1);
	e ^= dzu & (e ^ (uint32_t)-1076);
	uint64_t dm = (uint64_t)(dzu & 1) - 1;
	s &= dm;
	q &= dm;
	return make(s, *(int32_t *)&e, q);
}
#endif

#if !FNDSA_ASM_CORTEXM4
/* see sign_inner.h */
fpr
fpr_sqrt(fpr x)
{
	/* Extract exponent and mantissa. By assumption, the operand is
	   non-negative, hence we can ignore the sign bit (we must still
	   mask it out because sqrt() should work on -0.0). We want the
	   "true" exponent corresponding to a mantissa between 1 (inclusive)
	   and 2 (exclusive). */
	uint64_t xu = (x & M52) | ((uint64_t)1 << 52);
	uint32_t ex = (uint32_t)(x >> 52) & 0x7FF;
	int32_t e = *(int32_t *)&ex - 1023;

	/* If the exponent is odd, then we double the mantissa, and subtract
	   1 from the exponent. We can then halve the exponent. */
	xu += xu & -(uint64_t)((uint32_t)e & 1);
	e >>= 1;

	/* Double the mantissa to make it an integer in [2^53,2^55-1]. */
	xu <<= 1;

	/* xu represents an integer between 1 (inclusive) and 4
	   (exclusive) in a fixed-point notation (53 fractional bits).
	   We compute the square root bit by bit. */
	uint64_t q = 0;
	uint64_t s = 0;
	uint64_t r = (uint64_t)1 << 53;
	for (int i = 0; i < 54; i ++) {
		uint64_t t = s + r;
		uint64_t b = ((xu - t) >> 63) - 1;
		s += b & (r << 1);
		xu -= t & b;
		q += r & b;
		xu <<= 1;
		r >>= 1;
	}

	/* Now q is a rounded-low 54-bit value, with a leading 1, then
	   52 fractional digits, and an additional guard bit. We add an
	   extra sticky bit to account for what remains of the operand. */
	q <<= 1;
	q |= (xu | -xu) >> 63;

	/* Result q is in [2^54,2^55-1]; we bias the exponent by 54 bits
	   (since the computed e is, at this point, the "true" exponent). */
	e -= 54;

	/* If the source value was zero, then we computed the square root
	   of 2^53 and set the exponent to -512, both of which are wrong
	   and must be corrected. */
	q &= -(uint64_t)((ex + 0x7FF) >> 11);
	return make_z(0, e, q);
}
#endif

#endif
