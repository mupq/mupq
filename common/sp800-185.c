#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "keccakf1600.h"
#include "sp800-185.h"

#ifdef PROFILE_HASHING
#include "hal.h"
extern unsigned long long hash_cycles;
#endif

static size_t left_encode(uint8_t *encbuf, uint64_t value) {
    uint64_t v;
    size_t n, i;

    for (v = value, n = 0; v && (n < sizeof(uint64_t)); n++, v >>= 8) {
        ; /* empty */
    }
    if (n == 0) {
        n = 1;
    }
    for (i = 1; i <= n; i++) {
        encbuf[i] = (uint8_t)(value >> (8 * (n-i)));
    }
    encbuf[0] = (uint8_t)n;
    return n + 1;
}

void cshake128_inc_init(shake128incctx *state, const uint8_t *name, size_t namelen, const uint8_t *cstm, size_t cstmlen) {
    uint8_t encbuf[sizeof(uint64_t)+1];

    shake128_inc_init(state);

    shake128_inc_absorb(state, encbuf, left_encode(encbuf, SHAKE128_RATE));

    shake128_inc_absorb(state, encbuf, left_encode(encbuf, namelen * 8));
    shake128_inc_absorb(state, name, namelen);

    shake128_inc_absorb(state, encbuf, left_encode(encbuf, cstmlen * 8));
    shake128_inc_absorb(state, cstm, cstmlen);

    if (state->ctx[25] != 0) {
        state->ctx[25] = SHAKE128_RATE - 1;
        encbuf[0] = 0;
        shake128_inc_absorb(state, encbuf, 1);
    }
}

void cshake128_inc_absorb(shake128incctx *state, const uint8_t *input, size_t inlen) {
    shake128_inc_absorb(state, input, inlen);
}

void cshake128_inc_finalize(shake128incctx *state) {
#ifdef PROFILE_HASHING
    uint64_t t0 = hal_get_time();
#endif
    size_t i;
    uint8_t t[200];
    for (i = 0; i < SHAKE128_RATE; ++i) {
        t[i] = 0;
    }
    t[state->ctx[25]] = 0x04;
    t[SHAKE128_RATE - 1] |= 128;

    KeccakF1600_StateXORBytes(state->ctx, t, 0, SHAKE128_RATE);
    state->ctx[25] = 0;
#ifdef PROFILE_HASHING
    uint64_t t1 = hal_get_time();
    hash_cycles += (t1-t0);
#endif
}

void cshake128_inc_squeeze(uint8_t *output, size_t outlen, shake128incctx *state) {
    shake128_inc_squeeze(output, outlen, state);
}

void cshake128_inc_ctx_release(shake128incctx *state) {
    // no-op for mupq
    // this is required for compatibility with code from PQClean
    // see https://github.com/PQClean/PQClean/pull/265
    (void) state;
}

void cshake128_inc_ctx_clone(shake128incctx *dest, const shake128incctx *src) {
    memcpy(dest, src, sizeof(shake128incctx));
}


void cshake256_inc_init(shake256incctx *state, const uint8_t *name, size_t namelen, const uint8_t *cstm, size_t cstmlen) {
    uint8_t encbuf[sizeof(uint64_t)+1];

    shake256_inc_init(state);

    shake256_inc_absorb(state, encbuf, left_encode(encbuf, SHAKE256_RATE));

    shake256_inc_absorb(state, encbuf, left_encode(encbuf, namelen * 8));
    shake256_inc_absorb(state, name, namelen);

    shake256_inc_absorb(state, encbuf, left_encode(encbuf, cstmlen * 8));
    shake256_inc_absorb(state, cstm, cstmlen);

    if (state->ctx[25] != 0) {
        state->ctx[25] = SHAKE256_RATE - 1;
        encbuf[0] = 0;
        shake256_inc_absorb(state, encbuf, 1);
    }
}

void cshake256_inc_absorb(shake256incctx *state, const uint8_t *input, size_t inlen) {
    shake256_inc_absorb(state, input, inlen);
}

void cshake256_inc_finalize(shake256incctx *state) {
#ifdef PROFILE_HASHING
    uint64_t t0 = hal_get_time();
#endif
    size_t i;
    uint8_t t[200];
    for (i = 0; i < SHAKE256_RATE; ++i) {
        t[i] = 0;
    }
    t[state->ctx[25]] = 0x04;
    t[SHAKE256_RATE - 1] |= 128;

    KeccakF1600_StateXORBytes(state->ctx, t, 0, SHAKE256_RATE);
    state->ctx[25] = 0;
#ifdef PROFILE_HASHING
    uint64_t t1 = hal_get_time();
    hash_cycles += (t1-t0);
#endif
}

void cshake256_inc_squeeze(uint8_t *output, size_t outlen, shake256incctx *state) {
    shake256_inc_squeeze(output, outlen, state);
}

/*************************************************
 * Name:        cshake128
 *
 * Description: cSHAKE128 XOF with non-incremental API
 *
 * Arguments:   - uint8_t *output: pointer to output
 *              - size_t outlen: requested output length in bytes
 *              - const uint8_t *name: pointer to function-name string
 *              - size_t namelen: length of function-name string in bytes
 *              - const uint8_t *cstm: pointer to non-empty customization string
 *              - size_t cstmlen: length of customization string in bytes
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen: length of input in bytes
 **************************************************/
void cshake128(uint8_t *output, size_t outlen,
               const uint8_t *name, size_t namelen,
               const uint8_t *cstm, size_t cstmlen,
               const uint8_t *input, size_t inlen) {
    shake128incctx state;
    cshake128_inc_init(&state, name, namelen, cstm, cstmlen);
    cshake128_inc_absorb(&state, input, inlen);
    cshake128_inc_finalize(&state);
    cshake128_inc_squeeze(output, outlen, &state);
}

/*************************************************
 * Name:        cshake256
 *
 * Description: cSHAKE256 XOF with non-incremental API
 *
 * Arguments:   - uint8_t *output: pointer to output
 *              - size_t outlen: requested output length in bytes
 *              - const uint8_t *name: pointer to function-name string
 *              - size_t namelen: length of function-name string in bytes
 *              - const uint8_t *cstm: pointer to non-empty customization string
 *              - size_t cstmlen: length of customization string in bytes
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen: length of input in bytes
 **************************************************/
void cshake256(uint8_t *output, size_t outlen,
               const uint8_t *name, size_t namelen,
               const uint8_t *cstm, size_t cstmlen,
               const uint8_t *input, size_t inlen) {
    shake256incctx state;
    cshake256_inc_init(&state, name, namelen, cstm, cstmlen);
    cshake256_inc_absorb(&state, input, inlen);
    cshake256_inc_finalize(&state);
    cshake256_inc_squeeze(output, outlen, &state);
}

void cshake256_inc_ctx_release(shake256incctx *state) {
    // no-op for mupq
    // this is required for compatibility with code from PQClean
    // see https://github.com/PQClean/PQClean/pull/265
    (void) state;
}

void cshake256_inc_ctx_clone(shake256incctx *dest, const shake256incctx *src) {
    memcpy(dest, src, sizeof(shake256incctx));
}
