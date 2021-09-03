/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */


#include "mpc_lowmc.h"
#include "mzd_additional.h"

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


#include "lowmc_129_129_4.h"

#define MPC_LOOP_CONST(function, result, first, second, sc)                                        \
  do {                                                                                             \
    for (unsigned int e = 0; e < (sc); ++e) {                                                      \
      function((result)[e], (first)[e], (second));                                                 \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_SHARED(function, result, first, second, sc)                                       \
  do {                                                                                             \
    for (unsigned int o = 0; o < (sc); ++o) {                                                      \
      function((result)[o], (first)[o], (second)[o]);                                              \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_SHARED_1(function, result, first, sc)                                             \
  do {                                                                                             \
    for (unsigned int o = 0; o < (sc); ++o) {                                                      \
      function((result)[o], (first)[o]);                                                           \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_CONST_C_0(function, result, first, second, sc)                                    \
  function((result)[0], (first)[0], (second))

#define MPC_LOOP_CONST_C_ch(function, result, first, second, sc, c)                                \
  do {                                                                                             \
    if (!(c)) {                                                                                    \
      MPC_LOOP_CONST_C_0(function, result, first, second, sc);                                     \
    } else if ((c) == (sc)) {                                                                      \
      function((result)[(sc)-1], first[(sc)-1], (second));                                         \
    }                                                                                              \
  } while (0)


/* MPC Sbox implementation for full instances */
#if !defined(NO_UINT64_FALLBACK)
static void mpc_and_uint64_192(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp = {{0}};

  for (unsigned int m = 0; m < SC_PROOF; ++m) {
    const unsigned int j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_uint64_192(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_uint64_192(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_uint64_192(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_uint64_192(&tmp, &r[m], &r[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    if (viewshift) {
      mzd_shift_right_uint64_192(&tmp, &res[m], viewshift);
      mzd_xor_uint64_192(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_uint64_192(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_uint64_192(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                                      const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp = {{0}};

  for (unsigned int m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned int j = m + 1;

    mzd_and_uint64_192(&res[m], &first[m], &second[m]);

    mzd_and_uint64_192(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    mzd_and_uint64_192(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    mzd_xor_uint64_192(&tmp, &r[m], &r[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    if (viewshift) {
      mzd_shift_right_uint64_192(&tmp, &res[m], viewshift);
      mzd_xor_uint64_192(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_uint64_192(&view->s[m], &res[m]);
    }
  }

  if (viewshift) {
    mzd_shift_left_uint64_192(&tmp, &view->s[SC_VERIFY - 1], viewshift);
    mzd_and_uint64_192(&res[SC_VERIFY - 1], &tmp, mask);
  } else {
    mzd_and_uint64_192(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], mask);
  }
}


#define bitsliced_step_1(sc, AND, ROL, MASK_A, MASK_B, MASK_C)                                     \
  mzd_local_t x2m[sc] = {{{0}}};                                                                   \
  mzd_local_t r0m[sc] = {{{0}}}, r1m[sc] = {{{0}}}, r2m[sc] = {{{0}}};                             \
  mzd_local_t x0s[sc] = {{{0}}}, x1s[sc] = {{{0}}}, r0s[sc] = {{{0}}}, r1s[sc] = {{{0}}};          \
                                                                                                   \
  for (unsigned int m = 0; m < (sc); ++m) {                                                        \
    AND(&x0s[m], &in[m], MASK_A);                                                                  \
    AND(&x1s[m], &in[m], MASK_B);                                                                  \
    AND(&x2m[m], &in[m], MASK_C);                                                                  \
                                                                                                   \
    ROL(&x0s[m], &x0s[m], 2);                                                                      \
    ROL(&x1s[m], &x1s[m], 1);                                                                      \
                                                                                                   \
    AND(&r0m[m], &rvec->s[m], MASK_A);                                                             \
    AND(&r1m[m], &rvec->s[m], MASK_B);                                                             \
    AND(&r2m[m], &rvec->s[m], MASK_C);                                                             \
                                                                                                   \
    ROL(&r0s[m], &r0m[m], 2);                                                                      \
    ROL(&r1s[m], &r1m[m], 1);                                                                      \
  }

#define bitsliced_step_2(sc, XOR, ROR)                                                             \
  for (unsigned int m = 0; m < sc; ++m) {                                                          \
    XOR(&r2m[m], &r2m[m], &x0s[m]);                                                                \
    XOR(&x0s[m], &x0s[m], &x1s[m]);                                                                \
    XOR(&r1m[m], &r1m[m], &x0s[m]);                                                                \
    XOR(&r0m[m], &r0m[m], &x0s[m]);                                                                \
    XOR(&r0m[m], &r0m[m], &x2m[m]);                                                                \
                                                                                                   \
    ROR(&x0s[m], &r2m[m], 2);                                                                      \
    ROR(&x1s[m], &r1m[m], 1);                                                                      \
                                                                                                   \
    XOR(&x0s[m], &x0s[m], &r0m[m]);                                                                \
    XOR(&out[m], &x0s[m], &x1s[m]);                                                                \
  }

static void mpc_sbox_prove_uint64_lowmc_129_129_4(mzd_local_t* out, const mzd_local_t* in,
                                                  view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_uint64_192, mzd_shift_left_uint64_192, mask_129_129_43_a,
                   mask_129_129_43_b, mask_129_129_43_c);

  // a & b
  mpc_and_uint64_192(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_uint64_192(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_uint64_192(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_uint64_192, mzd_shift_right_uint64_192);
}

static void mpc_sbox_verify_uint64_lowmc_129_129_4(mzd_local_t* out, const mzd_local_t* in,
                                                   view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_uint64_192, mzd_shift_left_uint64_192, mask_129_129_43_a,
                   mask_129_129_43_b, mask_129_129_43_c);

  // a & b
  mpc_and_verify_uint64_192(r0m, x0s, x1s, r2m, view, mask_129_129_43_c, 0);
  // b & c
  mpc_and_verify_uint64_192(r2m, x1s, x2m, r1s, view, mask_129_129_43_c, 1);
  // c & a
  mpc_and_verify_uint64_192(r1m, x0s, x2m, r0s, view, mask_129_129_43_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_uint64_192, mzd_shift_right_uint64_192);
}


#endif /* NO_UINT_FALLBACK */


/* TODO: get rid of the copies */
#define SBOX(sbox, y, x, views, rvec, n, shares, shares2)                                          \
  {                                                                                                \
    mzd_local_t tmp[shares];                                                                       \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      memcpy(tmp[count].w64, CONST_BLOCK(x[count], 0)->w64, sizeof(mzd_local_t));                  \
    }                                                                                              \
    sbox(tmp, tmp, views, rvec);                                                                   \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      memcpy(BLOCK(y[count], 0)->w64, tmp[count].w64, sizeof(mzd_local_t));                        \
    }                                                                                              \
  }                                                                                                \
  while (0)

#define SBOX_uint64(sbox, y, x, views, rvec, n, shares, shares2)                                   \
  do {                                                                                             \
    uint64_t in[shares];                                                                           \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      in[count] = CONST_BLOCK(x[count], 0)->w64[(n) / (sizeof(word) * 8) - 1];                     \
    }                                                                                              \
    sbox(in, views, rvec->t);                                                                      \
    for (unsigned int count = 0; count < shares2; ++count) {                                       \
      memcpy(BLOCK(y[count], 0)->w64, CONST_BLOCK(x[count], 0)->w64,                               \
             ((n) / (sizeof(word) * 8) - 1) * sizeof(word));                                       \
      BLOCK(y[count], 0)->w64[(n) / (sizeof(word) * 8) - 1] = in[count];                           \
    }                                                                                              \
  } while (0)

#if !defined(NO_UINT64_FALLBACK)
#define IMPL uint64

// uint64 based implementation
#include "lowmc_128_128_20_fns_uint64.h"
#include "mpc_lowmc.c.i"

#include "lowmc_129_129_4_fns_uint64.h"
#include "mpc_lowmc.c.i"

#include "lowmc_192_192_30_fns_uint64.h"
#include "mpc_lowmc.c.i"

#include "lowmc_192_192_4_fns_uint64.h"
#include "mpc_lowmc.c.i"

#include "lowmc_256_256_38_fns_uint64.h"
#include "mpc_lowmc.c.i"

#include "lowmc_255_255_4_fns_uint64.h"
#include "mpc_lowmc.c.i"
#endif


zkbpp_lowmc_implementation_f get_zkbpp_lowmc_implementation(const lowmc_parameters_t* lowmc) {
  assert((lowmc->m == 43 && lowmc->n == 129) || (lowmc->m == 64 && lowmc->n == 192) ||
         (lowmc->m == 85 && lowmc->n == 255) ||
         (lowmc->m == 10 && (lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256)));


#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->m == 10) {
    switch (lowmc->n) {
    }
  }

  if (lowmc->n == 129 && lowmc->m == 43) {
    return mpc_lowmc_prove_uint64_lowmc_129_129_4;
  }
#endif

  return NULL;
}

zkbpp_lowmc_verify_implementation_f
get_zkbpp_lowmc_verify_implementation(const lowmc_parameters_t* lowmc) {
  assert((lowmc->m == 43 && lowmc->n == 129) || (lowmc->m == 64 && lowmc->n == 192) ||
         (lowmc->m == 85 && lowmc->n == 255) ||
         (lowmc->m == 10 && (lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256)));


#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->m == 10) {
    switch (lowmc->n) {
    }
  }

  if (lowmc->n == 129 && lowmc->m == 43) {
    return mpc_lowmc_verify_uint64_lowmc_129_129_4;
  }
#endif

  return NULL;
}

#if !defined(NO_UINT64_FALLBACK)

static void mzd_share_uint64_192(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                                 const mzd_local_t* v3) {
  mzd_xor_uint64_192(r, v1, v2);
  mzd_xor_uint64_192(r, r, v3);
}

#endif


zkbpp_share_implementation_f get_zkbpp_share_implentation(const lowmc_parameters_t* lowmc) {

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n == 129 || lowmc->n == 128) {
    return mzd_share_uint64_192;
  }
#endif

  return NULL;
}
