#ifndef SHA3_H__
#define SHA3_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "fips202.h"
#define shake_context shake256incctx
#define shake_init(sc,size) shake256_inc_init(sc);
#define shake_inject(sc, data, len) shake256_inc_absorb(sc, data, len)
#define shake_flip(sc) shake256_inc_finalize(sc);
#define shake_extract(sc, out, len) shake256_inc_squeeze(out, len, sc)

#ifdef __cplusplus
}
#endif

#endif
