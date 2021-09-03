/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */


#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "compat.h"
#include "mzd_additional.h"

#if !defined(_MSC_VER) && !defined(static_assert)
#define static_assert _Static_assert
#endif

static_assert(((sizeof(mzd_local_t) + 0x1f) & ~0x1f) == 32, "sizeof mzd_local_t not supported");



/* implementation of copy */

void mzd_copy_uint64_128(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 2; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}

void mzd_copy_uint64_192(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 3; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}

void mzd_copy_uint64_256(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 4; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}


/* implementation of mzd_xor and variants */


static void mzd_xor_uint64_block(block_t* rblock, const block_t* fblock, const block_t* sblock,
                                 const unsigned int len) {
  for (unsigned int i = 0; i < len; ++i) {
    rblock->w64[i] = fblock->w64[i] ^ sblock->w64[i];
  }
}

static void mzd_xor_uint64_blocks(block_t* rblock, const block_t* fblock, const block_t* sblock,
                                  const unsigned int len) {
  for (unsigned int i = len; i; --i, ++rblock, ++fblock, ++sblock) {
    mzd_xor_uint64_block(rblock, fblock, sblock, 4);
  }
}

void mzd_xor_uint64_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
}

void mzd_xor_uint64_192(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
}

void mzd_xor_uint64_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}

void mzd_xor_uint64_640(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
  mzd_xor_uint64_block(BLOCK(res, 2), CONST_BLOCK(first, 2), CONST_BLOCK(second, 2), 2);
}



/* implementation of mzd_and_* and variants */


static inline void mzd_and_uint64_block(block_t* rblock, const block_t* fblock,
                                        const block_t* sblock, const unsigned int len) {
  for (unsigned int i = 0; i < len; ++i) {
    rblock->w64[i] = fblock->w64[i] & sblock->w64[i];
  }
}

void mzd_and_uint64_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
}

void mzd_and_uint64_192(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
}

void mzd_and_uint64_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}

/* shifts and rotations */

void mzd_shift_left_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = block->w64[1] >> count;
}

void mzd_shift_left_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = block->w64[2] >> count;
}

void mzd_shift_left_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[3] = (block->w64[3] << count) | (block->w64[2] >> right_count);
  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = (block->w64[2] >> count) | (block->w64[3] << left_count);
  rblock->w64[3] = block->w64[3] >> count;
}


void mzd_mul_v_parity_uint64_128_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* vblock = CONST_BLOCK(v, 0);

  word res = 0;
  for (unsigned int i = 15; i; --i) {
    const block_t* Ablock = CONST_BLOCK(At, 15 - i);
    const word parity1 =
        parity64_uint64((vblock->w64[0] & Ablock->w64[0]) ^ (vblock->w64[1] & Ablock->w64[1]));
    const word parity2 =
        parity64_uint64((vblock->w64[0] & Ablock->w64[2]) ^ (vblock->w64[1] & Ablock->w64[3]));
    res |= (parity1 | (parity2 << 1)) << (64 - (2 * i));
  }
  cblock->w64[0] = 0;
  cblock->w64[1] = res;
}




static void clear_uint64_block(block_t* block, const unsigned int idx) {
  for (unsigned int i = 0; i < idx; ++i) {
    block->w64[i] = 0;
  }
}

static void clear_uint64_blocks(block_t* block, unsigned int len) {
  for (; len; --len, ++block) {
    clear_uint64_block(block, 4);
  }
}

static void mzd_xor_mask_uint64_block(block_t* rblock, const block_t* fblock, const word mask,
                                      const unsigned int idx) {
  for (unsigned int i = 0; i < idx; ++i) {
    rblock->w64[i] ^= fblock->w64[i] & mask;
  }
}

void mzd_addmul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 1) {
      const word mask1 = -(idx & 1);
      const word mask2 = -((idx >> 1) & 1);
      cblock->w64[0] ^= (Ablock->w64[0] & mask1) ^ (Ablock->w64[2] & mask2);
      cblock->w64[1] ^= (Ablock->w64[1] & mask1) ^ (Ablock->w64[3] & mask2);
    }
  }
}

void mzd_mul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 2);
  mzd_addmul_v_uint64_128(c, v, A);
}

void mzd_addmul_v_uint64_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  Ablock += 63;
  {
    word idx        = (*vptr) >> 63;
    const word mask = -(idx & 1);
    mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    Ablock++;
    vptr++;
  }

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    }
  }
}

void mzd_mul_v_uint64_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 3);
  mzd_addmul_v_uint64_129(c, v, A);
}

void mzd_addmul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    }
  }
}

void mzd_mul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 3);
  mzd_addmul_v_uint64_192(c, v, A);
}

void mzd_addmul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 4);
    }
  }
}

void mzd_mul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 4);
  mzd_addmul_v_uint64_256(c, v, A);
}

void mzd_mul_v_uint64_128_640(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  clear_uint64_blocks(BLOCK(c, 0), 2);
  clear_uint64_block(BLOCK(c, 2), 2);

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      for (unsigned int j = 0; j < 2; ++j, ++Ablock) {
        mzd_xor_mask_uint64_block(BLOCK(c, j), Ablock, mask, 4);
      }
      mzd_xor_mask_uint64_block(BLOCK(c, 2), Ablock, mask, 2);
    }
  }
}



// specific instances
// bit extract, non-constant time for mask, but mask is public in our calls
static word extract_bits(word in, word mask) {
  word res = 0;
  for (word bb = 1; mask != 0; bb <<= 1, mask &= (mask - 1)) {
    res |= bb & (-((word) !!(in & mask & -mask)));
  }
  return res;
}

static inline void mzd_shuffle_30_idx(mzd_local_t* x, const word mask, unsigned int idx) {
  const word w          = CONST_BLOCK(x, 0)->w64[idx];
  const word a          = extract_bits(w, mask) << 34;
  BLOCK(x, 0)->w64[idx] = a | extract_bits(w, ~mask);
}

void mzd_shuffle_128_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_30_idx(x, mask, 1);
}



// no SIMD
void mzd_addmul_v_uint64_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word idx = CONST_BLOCK(v, 0)->w64[1] >> 34;
  for (unsigned int i = 15; i; --i, idx >>= 2, ++Ablock) {
    const word mask1 = -(idx & 1);
    const word mask2 = -((idx >> 1) & 1);
    for (unsigned int j = 0; j < 2; ++j) {
      cblock->w64[j] ^= (Ablock->w64[j] & mask1) ^ (Ablock->w64[j + 2] & mask2);
    }
  }
}




