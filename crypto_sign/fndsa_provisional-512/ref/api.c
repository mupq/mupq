#include <stddef.h>
#include <string.h>

#include "api.h"
#include "fndsa.h"
#include "randombytes.h"

/* We use the signature size to autodetect the degree. */
#define LOGN   (CRYPTO_BYTES == FNDSA_SIGNATURE_SIZE(9) ? 9 : 10)

int
crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
	uint8_t seed[40];
	randombytes(seed, sizeof seed);
	fndsa_keygen_seeded(LOGN, seed, sizeof seed, sk, pk);
	return 0;
}

int
crypto_sign(unsigned char *sm, size_t *smlen,
	const unsigned char *m, size_t mlen,
	const unsigned char *sk)
{
	if (m != sm) {
		memmove(sm, m, mlen);
	}
	uint8_t seed[40];
	randombytes(seed, sizeof seed);
	size_t j = fndsa_sign_seeded(sk, FNDSA_SIGN_KEY_SIZE(LOGN),
		NULL, 0, FNDSA_HASH_ID_RAW, m, mlen,
		seed, sizeof seed, sm + mlen, FNDSA_SIGNATURE_SIZE(LOGN));
	if (j == 0) {
		return -1;
	}
	*smlen = mlen + FNDSA_SIGNATURE_SIZE(LOGN);
	return 0;
}

int
crypto_sign_open(unsigned char *m, size_t *mlen,
	const unsigned char *sm, size_t smlen,
	const unsigned char *pk)
{
	if (smlen < FNDSA_SIGNATURE_SIZE(LOGN)) {
		return -1;
	}
	size_t dlen = smlen - FNDSA_SIGNATURE_SIZE(LOGN);
	uint8_t tmp[((size_t)4 << LOGN) + 31];
	if (!fndsa_verify_temp(sm + dlen, FNDSA_SIGNATURE_SIZE(LOGN),
		pk, FNDSA_VRFY_KEY_SIZE(LOGN),
		NULL, 0, FNDSA_HASH_ID_RAW, sm, dlen, tmp, sizeof tmp))
	{
		return -1;
	}
	if (m != sm) {
		memmove(m, sm, dlen);
	}
	*mlen = dlen;
	return 0;
}
