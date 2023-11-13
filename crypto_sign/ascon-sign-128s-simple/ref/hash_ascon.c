#include <stdint.h>
#include <string.h>

#include "address.h"
#include "utils.h"
#include "params.h"
#include "hash.h"
#include "ascon_api.h"

/* For SHAKE256, there is no immediate reason to initialize at the start,
   so this function is an empty operation. */
void initialize_hash_function(ascon_sign_ctx* ctx)
{
    (void)ctx; /* Suppress an 'unused parameter' warning. */
}

/*
 * Computes PRF(pk_seed, sk_seed, addr)
 */
void prf_addr(unsigned char *out, const ascon_sign_ctx *ctx,
              const uint32_t addr[8])
{
    ASGN_VLA(uint8_t, buf, 2*ASGN_N + ASGN_ADDR_BYTES);

    memcpy(buf, ctx->pub_seed, ASGN_N);
    memcpy(buf + ASGN_N, addr, ASGN_ADDR_BYTES);
    memcpy(buf + ASGN_N + ASGN_ADDR_BYTES, ctx->sk_seed, ASGN_N);

    ascon_hash (out, ASGN_N, buf, 2*ASGN_N + ASGN_ADDR_BYTES);
}

/**
 * Computes the message-dependent randomness R, using a secret seed and an
 * optional randomization value as well as the message.
 */
void gen_message_random(unsigned char *R, const unsigned char *sk_prf,
                        const unsigned char *optrand,
                        const unsigned char *m, unsigned long long mlen,
                        const ascon_sign_ctx *ctx)
{
    (void)ctx;
    
    ASGN_VLA(uint8_t, temp, ASGN_N + ASGN_N + mlen);
 
    memcpy(temp, sk_prf, ASGN_N);
    memcpy(temp + ASGN_N, optrand, ASGN_N);
    memcpy(temp + (ASGN_N + ASGN_N), m, mlen);

    ascon_hash (R, ASGN_N, temp, ASGN_N + ASGN_N + mlen); 
}

/*
void printbytes(const char* desc, uint8_t* content, int count)
{
    printf("%s", desc);
    for(int i=0; i < count; i++)
    {
        printf("%02X ", content[i]);        
        if((i+1) % 16 == 0) 
        printf("\r\n");
    }
}
*/

/**
 * Computes the message hash using R, the public key, and the message.
 * Outputs the message digest and the index of the leaf. The index is split in
 * the tree index and the leaf index, for convenient copying to an address.
 */
void hash_message(unsigned char *digest, uint64_t *tree, uint32_t *leaf_idx,
                  const unsigned char *R, const unsigned char *pk,
                  const unsigned char *m, unsigned long long mlen,
                  const ascon_sign_ctx *ctx)
{
    (void)ctx;
#define ASGN_TREE_BITS (ASGN_TREE_HEIGHT * (ASGN_D - 1))
#define ASGN_TREE_BYTES ((ASGN_TREE_BITS + 7) / 8)
#define ASGN_LEAF_BITS ASGN_TREE_HEIGHT
#define ASGN_LEAF_BYTES ((ASGN_LEAF_BITS + 7) / 8)
#define ASGN_DGST_BYTES (ASGN_FORS_MSG_BYTES + ASGN_TREE_BYTES + ASGN_LEAF_BYTES)

    ASGN_VLA(uint8_t, buf, ASGN_DGST_BYTES);
    
    unsigned char *bufp = buf;
   
    ASGN_VLA(uint8_t, temp, ASGN_N + ASGN_PK_BYTES + mlen);
   
    unsigned char *temp_p = temp;

    memcpy(temp_p, R, ASGN_N); temp_p += ASGN_N;
    memcpy(temp_p, pk, ASGN_PK_BYTES); temp_p += ASGN_PK_BYTES;
    memcpy(temp_p, m, mlen);

    ascon_hash (buf, ASGN_DGST_BYTES, temp, ASGN_N + ASGN_PK_BYTES + mlen);

    // printbytes("hash_message() buf : ", buf , ASGN_DGST_BYTES);

    memcpy(digest, bufp, ASGN_FORS_MSG_BYTES);
    bufp += ASGN_FORS_MSG_BYTES;

#if ASGN_TREE_BITS > 64
    #error For given height and depth, 64 bits cannot represent all subtrees
#endif

    *tree = bytes_to_ull(bufp, ASGN_TREE_BYTES);
    *tree &= (~(uint64_t)0) >> (64 - ASGN_TREE_BITS);
    bufp += ASGN_TREE_BYTES;

    *leaf_idx = (uint32_t)bytes_to_ull(bufp, ASGN_LEAF_BYTES);
    *leaf_idx &= (~(uint32_t)0) >> (32 - ASGN_LEAF_BITS);
}
