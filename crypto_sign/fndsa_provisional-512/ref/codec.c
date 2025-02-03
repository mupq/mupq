/*
 * Encoding/decoding primitives.
 */

#include "inner.h"

/* see inner.h */
size_t
trim_i8_encode(unsigned logn, const int8_t *f, unsigned nbits, uint8_t *d)
{
	size_t n = (size_t)1 << logn;
	if (nbits == 8) {
		memmove(d, f, n);
		return n;
	}
	size_t j = 0;
	uint32_t acc = 0;
	unsigned acc_len = 0;
	uint32_t mask = ((uint32_t)1 << nbits) - 1;
	for (size_t i = 0; i < n; i ++) {
		acc = (acc << nbits) | ((uint32_t)f[i] & mask);
		acc_len += nbits;
		if (acc_len >= 8) {
			acc_len -= 8;
			d[j ++] = (uint8_t)(acc >> acc_len);
		}
	}
	return j;
}

/* see inner.h */
size_t
trim_i8_decode(unsigned logn, const uint8_t *d, int8_t *f, unsigned nbits)
{
	size_t needed = ((size_t)nbits << logn) >> 3;
	size_t j = 0;
	uint32_t acc = 0;
	unsigned acc_len = 0;
	uint32_t mask1 = (1 << nbits) - 1;
	uint32_t mask2 = 1 << (nbits - 1);
	for (size_t i = 0; i < needed; i ++) {
		acc = (acc << 8) | d[i];
		acc_len += 8;
		while (acc_len >= nbits) {
			acc_len -= nbits;
			uint32_t w = (acc >> acc_len) & mask1;
			w |= -(w & mask2);
			if (w == -mask2) {
				return 0;
			}
			f[j ++] = (int8_t)*(int32_t *)&w;
		}
	}
	return needed;
}

/* see inner.h */
size_t
mqpoly_encode(unsigned logn, const uint16_t *h, uint8_t *d)
{
	size_t n = (size_t)1 << logn;
	size_t j = 0;
	for (size_t i = 0; i < n; i += 4) {
		uint32_t h0 = h[i + 0];
		uint32_t h1 = h[i + 1];
		uint32_t h2 = h[i + 2];
		uint32_t h3 = h[i + 3];
		d[j + 0] = (uint8_t)(h0 >> 6);
		d[j + 1] = (uint8_t)((h0 << 2) | (h1 >> 12));
		d[j + 2] = (uint8_t)(h1 >> 4);
		d[j + 3] = (uint8_t)((h1 << 4) | (h2 >> 10));
		d[j + 4] = (uint8_t)(h2 >> 2);
		d[j + 5] = (uint8_t)((h2 << 6) | (h3 >> 8));
		d[j + 6] = (uint8_t)h3;
		j += 7;
	}
	return j;
}

#if !FNDSA_ASM_CORTEXM4
/* see inner.h */
size_t
mqpoly_decode(unsigned logn, const uint8_t *d, uint16_t *h)
{
	size_t n = (size_t)1 << logn;
	size_t j = 0;
	uint32_t ov = 0xFFFFFFFF;
	for (size_t i = 0; i < n; i += 4) {
		uint32_t d0 = d[j + 0];
		uint32_t d1 = d[j + 1];
		uint32_t d2 = d[j + 2];
		uint32_t d3 = d[j + 3];
		uint32_t d4 = d[j + 4];
		uint32_t d5 = d[j + 5];
		uint32_t d6 = d[j + 6];
		j += 7;
		uint32_t h0 = (d0 << 6) | (d1 >> 2);
		uint32_t h1 = ((d1 << 12) | (d2 << 4) | (d3 >> 4)) & 0x3FFF;
		uint32_t h2 = ((d3 << 10) | (d4 << 2) | (d5 >> 6)) & 0x3FFF;
		uint32_t h3 = ((d5 << 8) | d6) & 0x3FFF;
		h[i + 0] = h0;
		h[i + 1] = h1;
		h[i + 2] = h2;
		h[i + 3] = h3;
		ov &= h0 - 12289;
		ov &= h1 - 12289;
		ov &= h2 - 12289;
		ov &= h3 - 12289;
	}
	if ((ov >> 16) == 0) {
		return 0;
	} else {
		return j;
	}
}
#endif

/* see inner.h */
int
comp_encode(unsigned logn, const int16_t *s, uint8_t *d, size_t dlen)
{
	size_t n = (size_t)1 << logn;
	uint32_t acc = 0;
	unsigned acc_len = 0;
	size_t j = 0;
	for (size_t i = 0; i < n; i ++) {
		/* Invariant: acc_len <= 7 */
		int32_t x = s[i];
		if (x < -2047 || x > +2047) {
			return 0;
		}

		/* Get sign mask and absolute value. */
		uint32_t sw = (uint32_t)(x >> 16);
		uint32_t w = ((uint32_t)x ^ sw) - sw;

		/* Encode sign bit and low 7 bits of the absolute value. */
		acc = (acc << 8) | (sw & 0x80) | (w & 0x7F);
		acc_len += 8;

		/* Encode the high bits. Since |x| <= 2047, the high bits
		   have a value on [0,15], hence we add at most 16 bits. */
		unsigned wh = (w >> 7) + 1;
		acc = (acc << wh) | 1;
		acc_len += wh;

		/* We appended at most 8 + 15 + 1 = 24 bits, so the total
		   accumulated length in acc is at most 31 bits. */
		while (acc_len >= 8) {
			acc_len -= 8;
			if (j >= dlen) {
				return 0;
			}
			d[j ++] = (uint8_t)(acc >> acc_len);
		}
	}

	/* Flush remaining bits (if any). */
	if (acc_len > 0) {
		if (j >= dlen) {
			return 0;
		}
		d[j ++] = (uint8_t)(acc << (8 - acc_len));
	}

	/* Pad with zeros. */
	while (j < dlen) {
		d[j ++] = 0;
	}
	return 1;
}

#if !FNDSA_ASM_CORTEXM4
/* see inner.h */
int
comp_decode(unsigned logn, const uint8_t *d, size_t dlen, int16_t *s)
{
	size_t n = (size_t)1 << logn;
	uint32_t acc = 0;
	unsigned acc_len = 0;
	size_t j = 0;
	for (size_t i = 0; i < n; i ++) {
		/* Invariant: acc_len <= 7 */

		/* Get next 8 bits and split into sign bit (t) and low
		   bits of the absolute value (m). */
		if (j >= dlen) {
			return 0;
		}
		acc = (acc << 8) | d[j ++];
		uint32_t m = acc >> acc_len;
		uint32_t t = (m >> 7) & 1;
		m &= 0x7F;

		/* Get next bits until a 1 is reached.
		   We'll need up to two extra bytes. */
		for (;;) {
			if (acc_len == 0) {
				if (j >= dlen) {
					return 0;
				}
				acc = (acc << 8) | d[j ++];
				acc_len = 8;
			}
			acc_len --;
			if (((acc >> acc_len) & 1) != 0) {
				break;
			}
			m += 0x80;
			if (m > 2047) {
				return 0;
			}
		}

		/* Reject "-0" (which is an invalid encoding). */
		if (m == 0 && t != 0) {
			return 0;
		}

		m = (m ^ -t) + t;
		s[i] = (int16_t)*(int32_t *)&m;
	}

	/* Check that the unused bits are all zero. */
	if (acc_len > 0) {
		if ((acc & ((1 << acc_len) - 1)) != 0) {
			return 0;
		}
	}
	while (j < dlen) {
		if (d[j ++] != 0) {
			return 0;
		}
	}
	return 1;
}
#endif
