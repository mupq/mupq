/*
 * Utility functions.
 */

#include "inner.h"

/* see inner.h */
void
hash_to_point(unsigned logn,
        const uint8_t *nonce, const uint8_t *hashed_vrfy_key,
        const void *ctx, size_t ctx_len,
        const char *hash_id, const void *hv, size_t hv_len,
        uint16_t *c)
{
	/*
	 * If hash_id starts with a single byte of value 0xFF then we
	 * use the original Falcon mode (TODO: remove this feature, it
	 * is obsolescent).
	 *
	 * Raw message:
	 *   nonce || hashed_vrfy_key || 0x00 || len(ctx) || ctx || message
	 * Pre-hashed message:
	 *   nonce || hashed_vrfy_key || 0x01 || len(ctx) || ctx || id || hv
	 * Original Falcon:
	 *   nonce || message
	 */
	shake_context sc;
	shake_init(&sc, 256);
	shake_inject(&sc, nonce, 40);
	if (*(const uint8_t *)hash_id == 0xFF) {
		/* Original Falcon mode */
		shake_inject(&sc, hv, hv_len);
	} else {
		shake_inject(&sc, hashed_vrfy_key, 64);
		uint8_t hb[2];
		size_t id_len;
		if (hash_id[0] == 0x00) {
			/* Raw message, no pre-hashing */
			hb[0] = 0x00;
			id_len = 0;
		} else {
			/* Pre-hashing.
			   ID is ASN.1 OID: 0x06 len value */
			hb[0] = 0x01;
			id_len = hash_id[1] + 2;
		}
		hb[1] = ctx_len;
		shake_inject(&sc, &hb, 2);
		shake_inject(&sc, ctx, ctx_len);
		shake_inject(&sc, hash_id, id_len);
		shake_inject(&sc, hv, hv_len);
	}
	shake_flip(&sc);

	size_t n = (size_t)1 << logn;
	size_t i = 0;
#if FNDSA_ASM_CORTEXM4
	uint8_t *sbuf = (uint8_t *)(void *)&sc;
	size_t j = 136;
#else
	uint8_t sbuf[136];
	size_t j = sizeof sbuf;
#endif
	while (i < n) {
		if (j == 136) {
#if FNDSA_ASM_CORTEXM4
			shake_extract(&sc, NULL, j);
#else
			shake_extract(&sc, sbuf, j);
#endif
			j = 0;
		}
		unsigned w = ((unsigned)sbuf[j] << 8) | sbuf[j + 1];
		j += 2;
		if (w < 61445) {
			while (w >= 12289) {
				w -= 12289;
			}
			c[i ++] = w;
		}
	}
}

#if FNDSA_AVX2
#if defined __GNUC__ || defined __clang__
#include <cpuid.h>
__attribute__((target("xsave")))
int
has_avx2(void)
{
	/* __get_cpuid_count() includes a check that CPUID is callable,
	   and that the requested leaf number is available. */
	unsigned eax, ebx, ecx, edx;
	if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx)) {
		/* Check AVX2 support by the hardware. */
		if ((ebx & (1 << 5)) != 0) {
			/* Also check that YMM registers have not been
			   disabled by the OS. */
			return (_xgetbv(0) & 0x06) == 0x06;
		}
	}
	return 0;
}
#elif _MSC_VER
int
has_avx2(void)
{
	int rr[4];
	/* Check that CPUID leaf 7 is accessible. */
	__cpuid(rr, 0);
	if (rr[0] < 7) {
		return 0;
	}
	/* Check that the hardware supports AVX2. */
	__cpuidex(rr, 7, 0);
	if ((rr[1] & (1 << 5)) == 0) {
		return 0;
	}
	/* Check that the YMM registers have not been disabled by the OS. */
	return (_xgetbv(0) & 0x06) == 0x06;
}
#else
#error Missing has_avx2() implementation (not GCC/Clang/MSVC)
#endif
#endif
