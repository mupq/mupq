/*
 * Gaussian generation of (f,g)
 */

#include "kgen_inner.h"

/* q = 12289, n = 256 -> kmax = 24 */
static const uint16_t gauss_256[] = {
	    1,     3,     6,    11,    22,    40,    73,   129,
	  222,   371,   602,   950,  1460,  2183,  3179,  4509,
	 6231,  8395, 11032, 14150, 17726, 21703, 25995, 30487,
	35048, 39540, 43832, 47809, 51385, 54503, 57140, 59304,
	61026, 62356, 63352, 64075, 64585, 64933, 65164, 65313,
	65406, 65462, 65495, 65513, 65524, 65529, 65532, 65534
};

/* q = 12289, n = 512 -> kmax = 17 */
static const uint16_t gauss_512[] = {
	    1,     4,    11,    28,    65,   146,   308,   615,
	 1164,  2083,  3535,  5692,  8706, 12669, 17574, 23285,
	29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000,
	63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
	65531, 65534
};

/* q = 12289, n = 1024 -> kmax = 12 */
static const uint16_t gauss_1024[] = {
	    2,     8,    28,    94,   280,   742,  1761,  3753,
	 7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
	61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533
};

/* see kgen_inner.h */
#if FNDSA_SHAKE256X4
void
sample_f(unsigned logn, shake256x4_context *pc, int8_t *f)
#else
void
sample_f(unsigned logn, shake_context *pc, int8_t *f)
#endif
{
	const uint16_t *tab;
	size_t tab_len;
	unsigned zz;
	switch (logn) {
	case 9:
		tab = gauss_512;
		tab_len = (sizeof gauss_512) / sizeof(uint16_t);
		zz = 1;
		break;
	case 10:
		tab = gauss_1024;
		tab_len = (sizeof gauss_1024) / sizeof(uint16_t);
		zz = 1;
		break;
	default:
		tab = gauss_256;
		tab_len = (sizeof gauss_256) / sizeof(uint16_t);
		zz = 1u << (8 - logn);
		break;
	}
	uint32_t kmax = (uint32_t)(tab_len >> 1) << 16;
	size_t n = (size_t)1 << logn;

	/* We loop until we sample an odd-parity polynomial. */
	for (;;) {
		unsigned parity = 0;
		size_t i = 0;
		while (i < n) {
			/* Sampling: we choose a random 16-bit y value.
			   We then start with -kmax, and add 1 for each
			   table entry which is lower than y. We accumulate
			   these values in the upper 16 bits of v.

			   For logn < 8, we use the table for degree 256
			   but add multiple samples together. */
			uint32_t v = 0;
			for (unsigned t = 0; t < zz; t ++) {
#if FNDSA_SHAKE256X4
				uint32_t y = shake256x4_next_u16(pc);
#else
				uint32_t y = shake_next_u16(pc);
#endif
				v -= kmax;
				for (size_t k = 0; k < tab_len; k ++) {
					v -= ((uint32_t)tab[k] - y)
						& ~(uint32_t)0xFFFF;
				}
			}
			int s = *(int32_t *)&v >> 16;
			/* If logn <= 4 then it may happen that s does
			   not fit in [-127,+127], requiring some extra
			   sampling. */
			if (s < -127 || s > +127) {
				continue;
			}
			f[i ++] = (int8_t)s;
			parity ^= v;
		}
		/* Parity has been computed in bit 16 of 'parity'. */
		if (((parity >> 16) & 1) != 0) {
			break;
		}
	}
}
