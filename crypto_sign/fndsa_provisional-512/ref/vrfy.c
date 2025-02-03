/*
 * Signature verification.
 */

#include "inner.h"

/*
 * Inner verification function; it enforces a specific degree range.
 */
static int
inner_verify(unsigned logn_min, unsigned logn_max,
	const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len,
	void *tmp, size_t tmp_len)
{
	/* Get header bytes for key and signature, check that they relate
	   to the same degree, and that it is acceptable. */
	if (sig_len == 0 || vrfy_key_len == 0) {
		return 0;
	}
	const uint8_t *sigbuf = (const uint8_t *)sig;
	const uint8_t *vkbuf = (const uint8_t *)vrfy_key;
	unsigned logn = vkbuf[0];
	if (logn < logn_min || logn > logn_max || sigbuf[0] != 0x30 + logn) {
		return 0;
	}

	/* Keys and signatures have known fixed sizes. */
	if (sig_len != FNDSA_SIGNATURE_SIZE(logn)
		|| vrfy_key_len != FNDSA_VRFY_KEY_SIZE(logn))
	{
		return 0;
	}

	/* Check that temporary area is large enough. */
	size_t n = (size_t)1 << logn;
	if (tmp_len < (n * 4 + 31)) {
		return 0;
	}
	tmp = (void *)(((uintptr_t)tmp + 31) & ~(uintptr_t)31);

	/* Hash verifying key (SHAKE256, 64-byte output).
	   It is omitted in the "original Falcon" mode. */
	uint8_t hk[64];
	if (*(const uint8_t *)id == 0xFF && id[1] == 0) {
		/* TODO: remove original Falcon mode? */
		memset(hk, 0, sizeof hk);
	} else {
		/* TODO: maybe try to use tmp for the SHAKE context?
		   This could save some stack space. */
		shake_context sc;
		shake_init(&sc, 256);
		shake_inject(&sc, vrfy_key, vrfy_key_len);
		shake_flip(&sc);
		shake_extract(&sc, hk, sizeof hk);
	}

	uint16_t *t1 = (uint16_t *)tmp;
	uint16_t *t2 = t1 + n;

	/* t1 <- h (verifying key, decoded, converted to ntt) */
	if (mqpoly_decode(logn, vkbuf + 1, t1) != vrfy_key_len - 1) {
		return 0;
	}
	mqpoly_ext_to_int(logn, t1);
	mqpoly_int_to_ntt(logn, t1);

	/* t2 <- s2 (signature, decoded, converted to ntt)
	   Also get the squared norm of s2. */
	if (!comp_decode(logn, sigbuf + 41, sig_len - 41, (int16_t *)t2)) {
		return 0;
	}
	uint32_t norm2 = mqpoly_sqnorm_signed(logn, t2);
	mqpoly_signed_to_int(logn, t2);
	mqpoly_int_to_ntt(logn, t2);

	/* t2 <- s2*h (converted to int) */
	mqpoly_mul_ntt(logn, t2, t1);
	mqpoly_ntt_to_int(logn, t2);

	/* Hash message into polynomial c (into t1, converted to int) */
	hash_to_point(logn, sigbuf + 1, hk,
		ctx, ctx_len, id, hv, hv_len, t1);
	mqpoly_ext_to_int(logn, t1);

	/* t1 <- s1 = c - s2*h (converted to ext), and compute its norm. */
	mqpoly_sub(logn, t1, t2);
	mqpoly_int_to_ext(logn, t1);
	uint32_t norm1 = mqpoly_sqnorm_ext(logn, t1);

	/* Signature is valid if the total squared norm of (s1,s2) is
	   small enough. Beware overflows. */
	if (norm1 >= -norm2) {
		return 0;
	}
	return mqpoly_sqnorm_is_acceptable(logn, norm1 + norm2);
}

#if FNDSA_AVX2
TARGET_AVX2
static int
avx2_inner_verify(unsigned logn_min, unsigned logn_max,
	const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len)
{
	/* Get header bytes for key and signature, check that they relate
	   to the same degree, and that it is acceptable. */
	if (sig_len == 0 || vrfy_key_len == 0) {
		return 0;
	}
	const uint8_t *sigbuf = (const uint8_t *)sig;
	const uint8_t *vkbuf = (const uint8_t *)vrfy_key;
	unsigned logn = vkbuf[0];
	if (logn < logn_min || logn > logn_max || sigbuf[0] != 0x30 + logn) {
		return 0;
	}

	/* Keys and signatures have known fixed sizes. */
	if (sig_len != FNDSA_SIGNATURE_SIZE(logn)
		|| vrfy_key_len != FNDSA_VRFY_KEY_SIZE(logn))
	{
		return 0;
	}

	uint16_t t1[1024], t2[1024];

	/* t1 <- h (verifying key, decoded, converted to ntt) */
	if (mqpoly_decode(logn, vkbuf + 1, t1) != vrfy_key_len - 1) {
		return 0;
	}
	avx2_mqpoly_ext_to_int(logn, t1);
	avx2_mqpoly_int_to_ntt(logn, t1);

	/* t2 <- s2 (signature, decoded, converted to ntt)
	   Also get the squared norm of s2. */
	if (!comp_decode(logn, sigbuf + 41, sig_len - 41, (int16_t *)t2)) {
		return 0;
	}
	uint32_t norm2 = avx2_mqpoly_sqnorm_signed(logn, t2);
	avx2_mqpoly_signed_to_int(logn, t2);
	avx2_mqpoly_int_to_ntt(logn, t2);

	/* t2 <- s2*h (converted to int) */
	avx2_mqpoly_mul_ntt(logn, t2, t1);
	avx2_mqpoly_ntt_to_int(logn, t2);

	/* Hash verifying key (SHAKE256, 64-byte output). */
	uint8_t hk[64];
	shake_context sc;
	shake_init(&sc, 256);
	shake_inject(&sc, vrfy_key, vrfy_key_len);
	shake_flip(&sc);
	shake_extract(&sc, hk, sizeof hk);

	/* Hash message into polynomial c (into t1, converted to int) */
	hash_to_point(logn, sigbuf + 1, hk,
		ctx, ctx_len, id, hv, hv_len, t1);
	avx2_mqpoly_ext_to_int(logn, t1);

	/* t1 <- s1 = c - s2*h (converted to ext), and compute its norm. */
	avx2_mqpoly_sub(logn, t1, t2);
	avx2_mqpoly_int_to_ext(logn, t1);
	uint32_t norm1 = avx2_mqpoly_sqnorm_ext(logn, t1);

	/* Signature is valid if the total squared norm of (s1,s2) is
	   small enough. Beware overflows. */
	if (norm1 >= -norm2) {
		return 0;
	}
	return mqpoly_sqnorm_is_acceptable(logn, norm1 + norm2);
}
#endif

/* see fndsa.h */
int
fndsa_verify(const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len)
{
#if FNDSA_AVX2
	if (has_avx2()) {
		return avx2_inner_verify(9, 10,
			sig, sig_len, vrfy_key, vrfy_key_len,
			ctx, ctx_len, id, hv, hv_len);
	}
#endif
	uint8_t tmp[4 * 1024 + 31];
	return inner_verify(9, 10,
		sig, sig_len, vrfy_key, vrfy_key_len,
		ctx, ctx_len, id, hv, hv_len, tmp, sizeof tmp);
}

/* see fndsa.h */
int
fndsa_verify_weak(const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len)
{
#if FNDSA_AVX2
	if (has_avx2()) {
		return avx2_inner_verify(2, 8,
			sig, sig_len, vrfy_key, vrfy_key_len,
			ctx, ctx_len, id, hv, hv_len);
	}
#endif
	uint8_t tmp[4 * 256 + 31];
	return inner_verify(2, 8,
		sig, sig_len, vrfy_key, vrfy_key_len,
		ctx, ctx_len, id, hv, hv_len, tmp, sizeof tmp);
}

/* see fndsa.h */
int
fndsa_verify_temp(const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len,
	void *tmp, size_t tmp_len)
{
	return inner_verify(9, 10,
		sig, sig_len, vrfy_key, vrfy_key_len,
		ctx, ctx_len, id, hv, hv_len, tmp, tmp_len);
}

/* see fndsa.h */
int
fndsa_verify_weak_temp(const void *sig, size_t sig_len,
        const void *vrfy_key, size_t vrfy_key_len,
        const void *ctx, size_t ctx_len,
        const char *id, const void *hv, size_t hv_len,
	void *tmp, size_t tmp_len)
{
	return inner_verify(2, 8,
		sig, sig_len, vrfy_key, vrfy_key_len,
		ctx, ctx_len, id, hv, hv_len, tmp, tmp_len);
}
