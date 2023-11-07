
/// @file utils_prng.c
/// @brief The implementation of PRNG related functions.
///

#include "utils_prng.h"

#include <string.h>  // for memcpy()



//
// Defining prng_set_publicinputs() and aes128ctr_publicinputs()
//
//

static inline uint32_t br_swap32(uint32_t x) {
    x = ((x & (uint32_t)0x00FF00FF) << 8)
        | ((x >> 8) & (uint32_t)0x00FF00FF);
    return (x << 16) | (x >> 16);
}

static inline void inc1_be(uint32_t *x) {
    uint32_t t = br_swap32(*x) + 1;
    *x = br_swap32(t);
}


int prng_set_publicinputs(prng_publicinputs_t *ctx, const unsigned char prng_seed[16]) {
    #ifdef _4ROUND_AES_
    aes128_4rounds_ctr_keyexp_publicinputs(&ctx->ctx, prng_seed);
    #else
    aes128_ecb_keyexp_publicinputs(&ctx->ctx, prng_seed);
    #endif
    ctx->used = RNG_OUTPUTLEN;
    ctx->ctr = 0;
    return 0;

}

static void ov_aes128_ctr_publicinputs(unsigned char *out, size_t outlen, const unsigned char *iv, uint32_t ctr, const aes128ctx_publicinputs *ctx) {
    uint32_t ivw[4];
    size_t i;
    unsigned char tmp[16];
    memcpy(ivw, iv, 12);
    ivw[ 3] = br_swap32(ctr);
    while (outlen > 16) {
        aes128_ecb_publicinputs(out, (unsigned char *)ivw, 1, ctx);
        out += 16;
        outlen -= 16;
        inc1_be(ivw + 3);
    }
    
    if(outlen > 0){
        aes128_ecb_publicinputs(tmp, (unsigned char *)ivw, 1, ctx);
        for (i = 0; i < outlen; i++) {
            out[i] = tmp[i];
        }
    }
}



static
int aes128ctr_publicinputs( unsigned char *out, unsigned nblocks, const unsigned char *n, uint32_t ctr, const prng_publicinputs_t *ctx ) {
    ov_aes128_ctr_publicinputs(out, nblocks * RNG_OUTPUTLEN, n, ctr, &ctx->ctx);
    return 0;
}



int prng_gen_publicinputs(prng_publicinputs_t *ctx, unsigned char *out, unsigned long outlen) {
    unsigned long long xlen = outlen;
    unsigned long long ready;
    uint8_t nonce[AES128CTR_NONCELEN] = {0};

    if (ctx->used < RNG_OUTPUTLEN) {
        ready = RNG_OUTPUTLEN - ctx->used;
        if (xlen <= ready) {
            ready = xlen;
        }
        memcpy(out, &ctx->buf[ctx->used], ready);
        ctx->used += ready;
        xlen -= ready;
        out += ready;
    }


    if (xlen >= RNG_OUTPUTLEN) {
        uint32_t nblocks = xlen / RNG_OUTPUTLEN;
        aes128ctr_publicinputs(out, nblocks, nonce, ctx->ctr, ctx);
        ctx->ctr += (RNG_OUTPUTLEN / AES128_BLOCKSIZE) * nblocks;
        xlen -= nblocks * RNG_OUTPUTLEN;
        out += nblocks * RNG_OUTPUTLEN;
    }

    if (xlen > 0) {
        aes128ctr_publicinputs(ctx->buf, 1, nonce, ctx->ctr, ctx);
        ctx->ctr += (RNG_OUTPUTLEN / AES128_BLOCKSIZE);
        memcpy(out, ctx->buf, xlen);
        ctx->used = xlen;
    }
    return outlen;
}


void prng_skip_publicinputs(prng_publicinputs_t *ctx, unsigned long outlen) {
    if (ctx->used + outlen <= RNG_OUTPUTLEN ) {
        ctx->used += outlen;
        return;
    }
    outlen -= (RNG_OUTPUTLEN - ctx->used);

    unsigned long n_blocks_skip = outlen / RNG_OUTPUTLEN;
    unsigned long rem = outlen - n_blocks_skip * RNG_OUTPUTLEN;
    uint8_t nonce[AES128CTR_NONCELEN] = {0};
    if (rem) {
        ctx->ctr += n_blocks_skip * (RNG_OUTPUTLEN / AES128_BLOCKSIZE);
        ctx->used = rem;
        aes128ctr_publicinputs(ctx->buf, 1, nonce, ctx->ctr, ctx);
        ctx->ctr += (RNG_OUTPUTLEN / AES128_BLOCKSIZE);
    } else {  // 0 == rem
        ctx->ctr += n_blocks_skip * (RNG_OUTPUTLEN / AES128_BLOCKSIZE);
        ctx->used = RNG_OUTPUTLEN;
    }
}

