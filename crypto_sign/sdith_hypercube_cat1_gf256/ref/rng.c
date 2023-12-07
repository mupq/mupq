#include "rng.h"

#include "param.h"
#include "types.h"
#include "fips202.h"
#include <stdlib.h>
#include <string.h>

EXPORT void sdith_rng_create_xof_ctx(XOF_CTX *ctx, void *in, int inBytes) {
  // TODO: other security levels
#if defined(CAT_1)
  //Keccak_HashInitialize_SHAKE128(inst);
  shake128_inc_init(ctx);
#else
  Keccak_HashInitialize_SHAKE256(inst);
#endif
  //Keccak_HashUpdate(inst, in, inBytes << 3);
  shake128_inc_absorb(ctx, in, inBytes);
  //Keccak_HashFinal(inst, NULL);
  shake128_inc_finalize(ctx);
}

EXPORT void sdith_rng_free_xof_ctx(XOF_CTX *ctx) {
   shake128_inc_ctx_release(ctx);
}

EXPORT void sdith_xof_next_bytes(XOF_CTX *ctx, void *out, int outLen) {
  //Keccak_HashSqueeze((Keccak_HashInstance *)ctx, out, outLen << 3);
  shake128_inc_squeeze(out, outLen, ctx);
}

EXPORT void sdith_xof_next_bytes_mod251(XOF_CTX *ctx, void *out, int outLen) {
  // Roughly sample 1.03x of original length, which is greater than 256/251.
  int len = outLen + (outLen >> 5);
  uint8_t buf[len];
  sdith_xof_next_bytes(ctx, buf, len);
  int bytes_remaining = len;
  uint8_t *buf_ptr = buf;
  uint8_t *out_buf = (uint8_t *)out;
  int counter = 0;
  while (counter < outLen) {
    if (*buf_ptr < 251) {
      out_buf[counter++] = *buf_ptr;
    }
    buf_ptr++;
    bytes_remaining--;
    if (bytes_remaining == 0) {
      sdith_xof_next_bytes(ctx, buf, len);
      bytes_remaining = len;
      buf_ptr = buf;
    }
  }
}

// TODO: move this outside
shake128incctx states[4];

EXPORT XOF4_CTX *sdith_rng_create_xof4_ctx(void **in, int inBytes) {
  //shake128incctx *xof4_ctx =
  //    (shake128incctx *)malloc(4*sizeof(shake128incctx));
  shake128incctx *xof4_ctx = states;

#if defined(CAT_1)
  //Keccak_HashInitializetimes4_SHAKE128(&xof4_ctx);
  shake128_inc_init(&xof4_ctx[0]);
  shake128_inc_init(&xof4_ctx[1]);
  shake128_inc_init(&xof4_ctx[2]);
  shake128_inc_init(&xof4_ctx[3]);
#else
  Keccak_HashInitializetimes4_SHAKE256(&xof4_ctx);
#endif
  //Keccak_HashUpdatetimes4(&xof4_ctx, (const uint8_t **)in, inBytes << 3);
  shake128_inc_absorb(&xof4_ctx[0], in[0], inBytes);
  shake128_inc_absorb(&xof4_ctx[1], in[1], inBytes);
  shake128_inc_absorb(&xof4_ctx[2], in[2], inBytes);
  shake128_inc_absorb(&xof4_ctx[3], in[3], inBytes);

  //Keccak_HashFinaltimes4(&xof4_ctx, NULL);
  shake128_inc_finalize(&xof4_ctx[0]);
  shake128_inc_finalize(&xof4_ctx[1]);
  shake128_inc_finalize(&xof4_ctx[2]);
  shake128_inc_finalize(&xof4_ctx[3]);

  return (XOF4_CTX *)&xof4_ctx;
}

EXPORT void sdith_rng_free_xof4_ctx(XOF4_CTX *ctx) { (void) ctx; }

EXPORT void sdith_xof4_next_bytes(XOF4_CTX *ctx, void **out, int outLen) {
  // Keccak_HashSqueezetimes4((Keccak_HashInstancetimes4 *)ctx, (uint8_t **)out,
  //                          outLen << 3);
  shake128incctx *xof4_ctx = states;
  shake128_inc_squeeze(out[0], outLen, &xof4_ctx[0]);
  shake128_inc_squeeze(out[1], outLen, &xof4_ctx[1]);
  shake128_inc_squeeze(out[2], outLen, &xof4_ctx[2]);
  shake128_inc_squeeze(out[3], outLen, &xof4_ctx[3]);
}

EXPORT void sdith_xof4_next_bytes_mod251(XOF4_CTX *ctx, void **out,
                                         int outLen) {
  // Roughly sample 1.03x of original length, which is greater than 256/251.
  int len = outLen + (outLen >> 5);
  uint8_t buf[len * 4];
  uint8_t *bufs[4] = {&buf[0], &buf[len], &buf[len * 2], &buf[len * 3]};
  int bytes_counter[4] = {0};
  while (bytes_counter[0] + bytes_counter[1] + bytes_counter[2] +
             bytes_counter[3] <
         outLen * 4) {
    sdith_xof4_next_bytes(ctx, (void **)bufs, len);
    for (uint64_t i = 0; i < 4; ++i) {
      uint8_t *buf_ptr = bufs[i];
      uint8_t *out_buf = (uint8_t *)out[i];
      int bytes_remaining = len;
      while (bytes_counter[i] < outLen) {
        if (*buf_ptr < 251) {
          out_buf[bytes_counter[i]++] = *buf_ptr;
        }
        buf_ptr++;
        bytes_remaining--;
        if (bytes_remaining == 0) {
          break;
        }
      }
    }
  }
}