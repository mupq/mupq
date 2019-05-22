#include "common.h"

#if defined(VECLEN)
#define VECLEN 1
#endif

typedef struct hash_ctx_s {
    union { uint8_t b[25*8]; uint64_t w[25]; } state;
    uint8_t position, squeezing;
} hash_ctx_s, hash_ctx_t[1];

/* Init cSHAKE(N="",S="ThreeBears") */
void INTERNAL threebears_cshake_init(hash_ctx_t ctx);

/* Update hash function by adding len bytes from in */
void INTERNAL hash_update (
    struct hash_ctx_s * __restrict__ ctx,
    const uint8_t *in,
    size_t len
);

/* Output len bytes to out */
void INTERNAL hash_output (
    hash_ctx_t ctx,
    uint8_t * __restrict__ out,
    size_t len
);

/* Vectorized hash */
void INTERNAL hash_times_n(
    uint8_t *output,
    unsigned outlen,
    unsigned stride,
    const hash_ctx_t ctx0,
    uint8_t iv,
    unsigned n
);

/* Destroy hash */
static inline void hash_destroy(hash_ctx_t ctx) {
    secure_bzero(ctx,sizeof(hash_ctx_t));
}
