/*
 * Top-level key pair generator functions.
 */

#include "kgen_inner.h"

static void
keygen_inner(unsigned logn, const void *seed, size_t seed_len,
	void *sign_key, void *vrfy_key, void *tmp)
{
	/* Ensure that tmp is 32-byte aligned. */
	tmp = (void *)(((uintptr_t)tmp + 31) & ~(uintptr_t)31);

	/* We allocate f and g at the end of the temporary area. */
	size_t n = (size_t)1 << logn;
	int8_t *f = (int8_t *)tmp + ((size_t)24 << logn);
	int8_t *g = f + n;

	/* Make a PRNG with the provided seed. */
#if FNDSA_SHAKE256X4
	shake256x4_context pc;
	shake256x4_init(&pc, seed, seed_len);
#else
	shake_context pc;
	shake_init(&pc, 256);
	shake_inject(&pc, seed, seed_len);
	shake_flip(&pc);
#endif

	for (;;) {
		/* Sample f and g, both with odd parity. */
		sample_f(logn, &pc, f);
		sample_f(logn, &pc, g);

		/* Ensure that ||(g, -f)|| < 1.17*sqrt(q),
		   i.e. that ||(g, -f)||^2 < (1.17^2)*q = 16822.4121  */
		int32_t sn = 0;
		for (size_t i = 0; i < n; i ++) {
			int32_t xf = f[i];
			int32_t xg = g[i];
			sn += xf * xf + xg * xg;
		}
		if (sn >= 16823) {
			continue;
		}

		/* f must be invertible modulo X^n+1 modulo q. */
		if (!mqpoly_is_invertible(logn, f, tmp)) {
			continue;
		}

		/* (f,g) must have an acceptable orthogonalized norm. */
		if (!check_ortho_norm(logn, f, g, tmp)) {
			continue;
		}

		/* Try to solve the NTRU equation. */
		if (!solve_NTRU(logn, f, g, tmp)) {
			continue;
		}

		/* We have F and G at the start of tmp. We encode the
		   private and public keys. */
		if (sign_key != NULL) {
			int8_t *F = tmp;
			uint8_t *buf = sign_key;
			buf[0] = 0x50 + logn;
			unsigned nbits;
			switch (logn) {
			case 2: case 3: case 4: case 5:
				nbits = 8;
				break;
			case 6: case 7:
				nbits = 7;
				break;
			case 8: case 9:
				nbits = 6;
				break;
			default:
				nbits = 5;
			}
			size_t j = 1;
			j += trim_i8_encode(logn, f, nbits, buf + j);
			j += trim_i8_encode(logn, g, nbits, buf + j);
			(void)trim_i8_encode(logn, F, 8, buf + j);
		}
		if (vrfy_key != NULL) {
			uint8_t *buf = vrfy_key;
			uint16_t *h = tmp;
			mqpoly_div_small(logn, g, f, h, h + n);
			buf[0] = 0x00 + logn;
			(void)mqpoly_encode(logn, h, buf + 1);
		}
		break;
	}
}

#if FNDSA_AVX2
TARGET_AVX2
static void
avx2_keygen_inner(unsigned logn, const void *seed, size_t seed_len,
	void *sign_key, void *vrfy_key, void *tmp)
{
	/* Ensure that tmp is 32-byte aligned. */
	tmp = (void *)(((uintptr_t)tmp + 31) & ~(uintptr_t)31);

	/* We allocate f and g at the end of the temporary area. */
	size_t n = (size_t)1 << logn;
	int8_t *f = (int8_t *)tmp + ((size_t)24 << logn);
	int8_t *g = f + n;

	/* Make a PRNG with the provided seed. */
#if FNDSA_SHAKE256X4
	shake256x4_context pc;
	shake256x4_init(&pc, seed, seed_len);
#else
	shake_context pc;
	shake_init(&pc, 256);
	shake_inject(&pc, seed, seed_len);
	shake_flip(&pc);
#endif

	for (;;) {
		/* Sample f and g, both with odd parity. */
		sample_f(logn, &pc, f);
		sample_f(logn, &pc, g);

		/* Ensure that ||(g, -f)|| < 1.17*sqrt(q),
		   i.e. that ||(g, -f)||^2 < (1.17^2)*q = 16822.4121  */
		int32_t sn = 0;
		for (size_t i = 0; i < n; i ++) {
			int32_t xf = f[i];
			int32_t xg = g[i];
			sn += xf * xf + xg * xg;
		}
		if (sn >= 16823) {
			continue;
		}

		/* f must be invertible modulo X^n+1 modulo q. */
		if (!avx2_mqpoly_is_invertible(logn, f, tmp)) {
			continue;
		}

		/* (f,g) must have an acceptable orthogonolized norm. */
		if (!avx2_check_ortho_norm(logn, f, g, tmp)) {
			continue;
		}

		/* Try to solve the NTRU equation. */
		if (!avx2_solve_NTRU(logn, f, g, tmp)) {
			continue;
		}

		/* We have F and G at the start of tmp. We encode the
		   private and public keys. */
		if (sign_key != NULL) {
			int8_t *F = tmp;
			uint8_t *buf = sign_key;
			buf[0] = 0x50 + logn;
			unsigned nbits;
			switch (logn) {
			case 2: case 3: case 4: case 5:
				nbits = 8;
				break;
			case 6: case 7:
				nbits = 7;
				break;
			case 8: case 9:
				nbits = 6;
				break;
			default:
				nbits = 5;
			}
			size_t j = 1;
			j += trim_i8_encode(logn, f, nbits, buf + j);
			j += trim_i8_encode(logn, g, nbits, buf + j);
			(void)trim_i8_encode(logn, F, 8, buf + j);
		}
		if (vrfy_key != NULL) {
			uint8_t *buf = vrfy_key;
			uint16_t *h = tmp;
			avx2_mqpoly_div_small(logn, g, f, h, h + n);
			buf[0] = 0x00 + logn;
			(void)mqpoly_encode(logn, h, buf + 1);
		}
		break;
	}
}
#endif

/* Custom wrappers to allocate the temporary buffers on the stack. Several
   wrappers are defined so that stack allocation is not always worst-case. */
#if FNDSA_AVX2
#define KEYGEN_WRAP(sz)   \
	static void keygen_ ## sz(unsigned logn, \
		const void *seed, size_t seed_len, \
		void *sign_key, void *vrfy_key) \
	{ \
		uint8_t tmp[(sz) * 26 + 31]; \
		if (has_avx2()) { \
			avx2_keygen_inner(logn, \
				seed, seed_len, sign_key, vrfy_key, tmp); \
		} else { \
			keygen_inner(logn, \
				seed, seed_len, sign_key, vrfy_key, tmp); \
		} \
	}
#else
#define KEYGEN_WRAP(sz)   \
	static void keygen_ ## sz(unsigned logn, \
		const void *seed, size_t seed_len, \
		void *sign_key, void *vrfy_key) \
	{ \
		uint8_t tmp[(sz) * 26 + 31]; \
		keygen_inner(logn, seed, seed_len, sign_key, vrfy_key, tmp); \
	}
#endif

KEYGEN_WRAP(32)
KEYGEN_WRAP(64)
KEYGEN_WRAP(128)
KEYGEN_WRAP(256)
KEYGEN_WRAP(512)
KEYGEN_WRAP(1024)

static int
keygen(unsigned logn, const void *seed, size_t seed_len,
	void *sign_key, void *vrfy_key, void *tmp, size_t tmp_len)
{
	/* If no seed is provided, uses the system RNG to get a
	   32-byte seed. */
	uint8_t seedbuf[32];
	if (seed == NULL) {
		if (!sysrng(seedbuf, sizeof seedbuf)) {
			goto fail;
		}
		seed = seedbuf;
		seed_len = sizeof seedbuf;
	}

	if (tmp == NULL) {
		/* If no temporary area is provided, call the relevant
		   wrapper to allocate it on the stack. */
		switch (logn) {
		case 6:
			keygen_64(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		case 7:
			keygen_128(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		case 8:
			keygen_256(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		case 9:
			keygen_512(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		case 10:
			keygen_1024(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		default:
			keygen_32(logn, seed, seed_len, sign_key, vrfy_key);
			break;
		}
	} else {
		/* Check that the provided temporary area is large enough.
		   We want 24*n bytes + enough room to ensure 32-byte
		   alilgnment. */
		if (tmp_len < (31 + ((size_t)24 << logn))) {
			goto fail;
		}
		keygen_inner(logn, seed, seed_len, sign_key, vrfy_key, tmp);
	}
	return 1;

fail:
	if (sign_key != NULL) {
		memset(sign_key, 0, FNDSA_SIGN_KEY_SIZE(logn));
	}
	if (vrfy_key != NULL) {
		memset(vrfy_key, 0, FNDSA_VRFY_KEY_SIZE(logn));
	}
	return 0;
}

/* see fndsa.h */
int
fndsa_keygen(unsigned logn, void *sign_key, void *vrfk_key)
{
	return keygen(logn, NULL, 0, sign_key, vrfk_key, NULL, 0);
}

/* see fndsa.h */
int
fndsa_keygen_temp(unsigned logn, void *sign_key, void *vrfk_key,
	void *tmp, size_t tmp_len)
{
	return keygen(logn, NULL, 0, sign_key, vrfk_key, tmp, tmp_len);
}

/* see fndsa.h */
void
fndsa_keygen_seeded(unsigned logn, const void *seed, size_t seed_len,
	void *sign_key, void *vrfk_key)
{
	(void)keygen(logn, seed, seed_len, sign_key, vrfk_key, NULL, 0);
}

/* see fndsa.h */
int
fndsa_keygen_seeded_temp(unsigned logn, const void *seed, size_t seed_len,
	void *sign_key, void *vrfk_key, void *tmp, size_t tmp_len)
{
	return keygen(logn, seed, seed_len, sign_key, vrfk_key, tmp, tmp_len);
}
