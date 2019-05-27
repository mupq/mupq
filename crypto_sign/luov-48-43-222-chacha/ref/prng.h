#ifndef RNG_H
#define RNG_H

#include "parameters.h"
#include <stdlib.h>
#include "buffer.h"
#include "LinearAlgebra.h"
#include "fips202.h"

//#define Sponge Keccak_HashInstance
#define Sponge shake128incctx

//#define squeezeBytes(S,D,L) Keccak_HashSqueeze (S,D,L * 8)
#define squeezeBytes(S,D,L) shake128_inc_squeeze(D,L,S);

void initializeAndAbsorb(Sponge *sponge, const unsigned char * seed, int len);
uint64_t squeezeuint64_t(Sponge *sponge, int bytes);

#ifdef PRNG_KECCAK

#define BLOCK_SIZE 64
#define PRNG_STATE Sponge
#define PRNG_INIT(S,K,I) Keccak_HashInitialize_SHAKE128(S); Keccak_HashUpdate(S, K, 32*8 ); Keccak_HashUpdate(S, I, 8); Keccak_HashFinal(S,0);
#define PRNG_GET_BLOCK(S,O) Keccak_HashSqueeze(S, O, BLOCK_SIZE*8);

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

