#ifndef RNG_H
#define RNG_H

#include "parameters.h"
#include <stdlib.h>
#include "buffer.h"
#include "LinearAlgebra.h"
#include "fips202.h"

#define __Sponge(NUM) shake##NUM##incctx
#define _Sponge(NUM) __Sponge(NUM)
#define Sponge _Sponge(SHAKENUM)

#define __shake_inc_init(NUM) shake##NUM##_inc_init
#define _shake_inc_init(NUM) __shake_inc_init(NUM)
#define shake_inc_init _shake_inc_init(SHAKENUM)
#define __shake_inc_absorb(NUM) shake##NUM##_inc_absorb
#define _shake_inc_absorb(NUM) __shake_inc_absorb(NUM)
#define shake_inc_absorb _shake_inc_absorb(SHAKENUM)
#define __shake_inc_finalize(NUM) shake##NUM##_inc_finalize
#define _shake_inc_finalize(NUM) __shake_inc_finalize(NUM)
#define shake_inc_finalize _shake_inc_finalize(SHAKENUM)
#define __shake_inc_squeeze(NUM) shake##NUM##_inc_squeeze
#define _shake_inc_squeeze(NUM) __shake_inc_squeeze(NUM)
#define shake_inc_squeeze _shake_inc_squeeze(SHAKENUM)

#define squeezeBytes(S,D,L) shake_inc_squeeze(D,L,S)

void initializeAndAbsorb(Sponge *sponge, const unsigned char * seed, int len);
uint64_t squeezeuint64_t(Sponge *sponge, int bytes);

#ifdef PRNG_KECCAK

#define BLOCK_SIZE 64
#define PRNG_STATE Sponge
#define PRNG_INIT(S,K,I) shake_inc_init(S); shake_inc_absorb(S, K, 32); shake_inc_absorb(S, I, 1); shake_inc_finalize(S);
#define PRNG_GET_BLOCK(S,O) shake_inc_squeeze(O, BLOCK_SIZE, S);

#endif

#ifdef PRNG_CHACHA
#include "chacha.h"

#define BLOCK_SIZE 64
#define PRNG_STATE chacha_ctx
#define PRNG_INIT(S,K,I) chacha_keysetup(S, K, 256); chacha_ivsetup(S, I, 0);
#define PRNG_GET_BLOCK(S,O) { unsigned char In[BLOCK_SIZE] = {0}; \
							  chacha_encrypt_bytes(S, In, O, BLOCK_SIZE); };
#endif

#endif

