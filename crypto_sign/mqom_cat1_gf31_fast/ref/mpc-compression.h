#ifndef MPC_COMPRESSION_H
#define MPC_COMPRESSION_H

#include "mpc-struct.h"
#include "field.h"

static inline void vec_compress(uint8_t* buf, const void* x, uint32_t size) {
    compress_points(buf, x, size);
}
static inline void vec_decompress(void* x, const uint8_t* buf, uint32_t size) {
    decompress_points(x, buf, size);
}
#define wit_compress(buf,wit) vec_compress(buf,wit,sizeof(mpc_wit_t))
#define unif_compress(buf,unif) vec_compress(buf,unif,sizeof(mpc_unif_t))
#define hint_compress(buf,hint) vec_compress(buf,hint,sizeof(mpc_hint_t))
#define br_compress(buf,br) vec_compress(buf,br,sizeof(mpc_broadcast_t))
#define wit_decompress(wit,buf) vec_decompress(wit,buf,sizeof(mpc_wit_t))
#define unif_decompress(unif,buf) vec_decompress(unif,buf,sizeof(mpc_unif_t))
#define hint_decompress(hint,buf) vec_decompress(hint,buf,sizeof(mpc_hint_t))
#define br_decompress(br,buf) vec_decompress(br,buf,sizeof(mpc_broadcast_t))

#define PARAM_WIT_SHORT_SIZE ((sizeof(mpc_wit_t)*5+7)>>3)
#define PARAM_UNIF_SHORT_SIZE ((sizeof(mpc_unif_t)*5+7)>>3)
#define PARAM_HINT_SHORT_SIZE ((sizeof(mpc_hint_t)*5+7)>>3)
#define PARAM_BR_SHORT_SIZE ((sizeof(mpc_broadcast_t)*5+7)>>3)

//#define PARAM_COMPRESSED_BR_SIZE (sizeof(mpc_unif_t))
#define PARAM_COMPRESSED_BR_SIZE ((sizeof(mpc_unif_t)*5+7)>>3)

#endif /* MPC_COMPRESSION_H */
