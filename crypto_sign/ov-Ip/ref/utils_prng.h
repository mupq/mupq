/// @file utils_prng.h
/// @brief the interface for adapting PRNG functions.
///
///

#ifndef _UTILS_PRNG_H_
#define _UTILS_PRNG_H_

#include <stdint.h>

#include "config.h"
#include "params.h"  // for macro _4ROUND_AES_

#include "utils_randombytes.h"   // declare randombytes(), including here for backward compatibility

//
// This is not for traditional PRNG functions. It is PRNG with PUBLIC INPUTS.
// It uses AES128CTR or CTR with only 4 rounds AES128.
//

#define AES128CTR_KEYLEN   16
#define AES128CTR_NONCELEN 16
#define AES128_BLOCKSIZE 16


#ifdef  __cplusplus
extern  "C" {
#endif

/////////////////  defination of prng_publicinputs_t  /////////////////////////////////

// fixsliced implementation processes two blocks in parallel
#define RNG_OUTPUTLEN 32

#include "aes.h"
#include "aes-publicinputs.h"


typedef struct {
    unsigned used;
    uint32_t ctr;
    unsigned char buf[RNG_OUTPUTLEN];
    aes128ctx_publicinputs ctx;
} prng_publicinputs_t;


///////////////// end of defination of prng_publicinputs_t  /////////////////////////////////

int prng_set_publicinputs(prng_publicinputs_t *ctx, const unsigned char prng_seed[16]);

int prng_gen_publicinputs(prng_publicinputs_t *ctx, unsigned char *out, unsigned long outlen);

void prng_skip_publicinputs(prng_publicinputs_t *ctx, unsigned long outlen);




#ifdef  __cplusplus
}
#endif



#endif // _UTILS_PRNG_H_


