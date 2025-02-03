/*
 * Top-level signature generation functions.
 */

#include "sign_inner.h"

/* Verified properties at this point:
      degree is acceptable
      encoded signing key has the proper size
      signature buffer is large enough to receive the result
      tmp is large enough (but not necessarily aligned)  */
static size_t
sign_step1(unsigned logn, const uint8_t *sign_key,
	const uint8_t *ctx, size_t ctx_len,
	const char *id, const uint8_t *hv, size_t hv_len,
	const uint8_t *seed, size_t seed_len,
	uint8_t *sig, void *tmp)
{
	size_t n = (size_t)1 << logn;

	/* Align tmp to a 32-byte boundary. */
	tmp = (void *)(((uintptr_t)tmp + 31) & ~(uintptr_t)31);

	int8_t *f = (int8_t *)tmp + ((size_t)74 << logn);
	int8_t *g = f + n;
	int8_t *F = g + n;
	int8_t *G = F + n;

	/* Decode the private key. Header byte and length have already
	   been verified. */
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
		break;
	}
	size_t k, j = 1;
	k = trim_i8_decode(logn, sign_key + j, f, nbits);
	if (k == 0) {
		return 0;
	}
	j += k;
	k = trim_i8_decode(logn, sign_key + j, g, nbits);
	if (k == 0) {
		return 0;
	}
	j += k;
	k = trim_i8_decode(logn, sign_key + j, F, 8);
	if (k == 0) {
		return 0;
	}

	/* Rebuild G and the public polynomial h:
	      h = g/f mod X^n+1 mod q
	      G = h*F mod X^n+1 mod q
	   We also compute the SHAKE256 hash of the verifying key. */
	uint16_t *t0 = tmp;
	uint16_t *t1 = t0 + n;
	/* t0 <- h = g/f */
	mqpoly_small_to_int(logn, g, t0);
	mqpoly_small_to_int(logn, f, t1);
	mqpoly_int_to_ntt(logn, t0);
	mqpoly_int_to_ntt(logn, t1);
	if (!mqpoly_div_ntt(logn, t0, t1)) {
		/* f is not invertible; the key is not valid */
		return 0;
	}
	/* t1 <- G = h*F */
	mqpoly_small_to_int(logn, F, t1);
	mqpoly_int_to_ntt(logn, t1);
	mqpoly_mul_ntt(logn, t1, t0);
	mqpoly_ntt_to_int(logn, t1);
	if (!mqpoly_int_to_small(logn, t1, G)) {
		/* coefficients of G are out-of-range */
		return 0;
	}
	/* t0 contains h (in ntt representation), we encode and hash
	   the verifying key.
	   TODO: if the original Falcon mode is retained, then we can
	   skip both encoding and hashing. */
	mqpoly_ntt_to_int(logn, t0);
	mqpoly_int_to_ext(logn, t0);
	uint8_t *vrfy_key = (uint8_t *)(t0 + n);
	vrfy_key[0] = 0x00 + logn;
	mqpoly_encode(logn, t0, vrfy_key + 1);

	/* We can use tmp (skipping t0 and t1) for the SHAKE256 context.
	   The buffer currently starts with t0 (2*n bytes) then the
	   encoded public key (less than 2*n bytes), leaving 70*n bytes,
	   i.e. at least 280 bytes since n >= 4; the shake_context
	   structure size is normally 208 bytes, so it fits well. */
	uint8_t hashed_key[64];
	shake_context *sc = (shake_context *)(t0 + 2 * n);
	shake_init(sc, 256);
	shake_inject(sc, vrfy_key, FNDSA_VRFY_KEY_SIZE(logn));
	shake_flip(sc);
	shake_extract(sc, hashed_key, sizeof hashed_key);

	/* We now have f, g, F and G decoded and verified. Hashed public
	   key is in hashed_key[]. We can proceed to the main signing loop. */
	return sign_core(logn, f, g, F, G, hashed_key,
		ctx, ctx_len, id, hv, hv_len,
		seed, seed_len, sig, tmp);

	/* TODO: maybe explicitly overwrite the whole temporary area with
	   zeros? Arguably this is mostly wasted time if the area is
	   allocated on the stack, and if tmp is provided explicitly then
	   it is the responsibility of the caller to do any appropriate
	   zeroizing. */
}

/* Custom wrappers to allocate the temporary buffers on the stack. Several
   wrappers are defined so that stack allocation is not always worst-case. */
#define SIGN_WRAP(sz)   \
	static size_t sign_ ## sz(unsigned logn, \
		const uint8_t *sign_key, \
		const uint8_t *ctx, size_t ctx_len, \
		const char *id, const uint8_t *hv, size_t hv_len, \
		const uint8_t *seed, size_t seed_len, \
		uint8_t *sig) \
	{ \
		uint8_t tmp[(sz) * 78 + 31]; \
		return sign_step1(logn, \
			sign_key, ctx, ctx_len, id, hv, hv_len, \
			seed, seed_len, sig, tmp); \
	}

SIGN_WRAP(32)
SIGN_WRAP(64)
SIGN_WRAP(128)
SIGN_WRAP(256)
SIGN_WRAP(512)
SIGN_WRAP(1024)

static size_t
sign_wrapper(int weak,
	const uint8_t *sign_key, size_t sign_key_len,
	const uint8_t *ctx, size_t ctx_len,
	const char *id, const uint8_t *hv, size_t hv_len,
	const uint8_t *seed, size_t seed_len,
	uint8_t *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	/* Signing key defines the degree to use. */
	if (sign_key_len == 0) {
		return 0;
	}
	unsigned head = sign_key[0];
	if ((head & 0xF0) != 0x50) {
		return 0;
	}
	unsigned logn = head & 0x0F;
	if (weak) {
		if (logn < 2 || logn > 8) {
			return 0;
		}
	} else {
		if (logn < 9 || logn > 10) {
			return 0;
		}
	}
	if (sign_key_len != FNDSA_SIGN_KEY_SIZE(logn)) {
		return 0;
	}
	if (sig == NULL) {
		return FNDSA_SIGNATURE_SIZE(logn);
	}
	if (max_sig_len < FNDSA_SIGNATURE_SIZE(logn)) {
		return 0;
	}

	/* We have checked that the degree is acceptable, the signing key
	   size is correct, and the signature will fit in the output buffer. */
	if (tmp == NULL) {
		switch (logn) {
		case 6:
			return sign_64(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		case 7:
			return sign_128(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		case 8:
			return sign_256(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		case 9:
			return sign_512(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		case 10:
			return sign_1024(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		default:
			return sign_32(logn,
				sign_key, ctx, ctx_len, id, hv, hv_len,
				seed, seed_len, sig);
		}
	} else {
		if (tmp_len < (((size_t)78 << logn) + 31)) {
			return 0;
		}
		return sign_step1(logn,
			sign_key, ctx, ctx_len, id, hv, hv_len,
			seed, seed_len, sig, tmp);
	}
}

/* see fndsa.h */
size_t
fndsa_sign(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len)
{
	return sign_wrapper(0, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		NULL, 0, sig, max_sig_len, NULL, 0);
}

/* see fndsa.h */
size_t
fndsa_sign_seeded(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len)
{
	return sign_wrapper(0, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		seed, seed_len, sig, max_sig_len, NULL, 0);
}

/* see fndsa.h */
size_t
fndsa_sign_temp(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_wrapper(0, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		NULL, 0, sig, max_sig_len, tmp, tmp_len);
}

/* see fndsa.h */
size_t
fndsa_sign_seeded_temp(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_wrapper(0, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		seed, seed_len, sig, max_sig_len, tmp, tmp_len);
}

/* see fndsa.h */
size_t
fndsa_sign_weak(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len)
{
	return sign_wrapper(1, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		NULL, 0, sig, max_sig_len, NULL, 0);
}

/* see fndsa.h */
size_t
fndsa_sign_weak_seeded(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len)
{
	return sign_wrapper(1, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		seed, seed_len, sig, max_sig_len, NULL, 0);
}

/* see fndsa.h */
size_t
fndsa_sign_weak_temp(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_wrapper(1, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		NULL, 0, sig, max_sig_len, tmp, tmp_len);
}

/* see fndsa.h */
size_t
fndsa_sign_weak_seeded_temp(const void *sign_key, size_t sign_key_len,
	const void *ctx, size_t ctx_len,
	const char *id, const void *hv, size_t hv_len,
	const void *seed, size_t seed_len,
	void *sig, size_t max_sig_len,
	void *tmp, size_t tmp_len)
{
	return sign_wrapper(1, sign_key, sign_key_len,
		ctx, ctx_len, id, hv, hv_len,
		seed, seed_len, sig, max_sig_len, tmp, tmp_len);
}
