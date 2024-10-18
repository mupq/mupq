// SPDX-License-Identifier: Apache-2.0
#include "api.h"
#include "hal.h"
#include "randombytes.h"
#include "sendfn.h"

#include <stdint.h>
#include <string.h>

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_BYTES NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_CIPHERTEXTBYTES NAMESPACE(CRYPTO_CIPHERTEXTBYTES)
#define MUPQ_CRYPTO_ALGNAME NAMESPACE(CRYPTO_ALGNAME)

#define MUPQ_crypto_kem_keypair NAMESPACE(crypto_kem_keypair)
#define MUPQ_crypto_kem_enc NAMESPACE(crypto_kem_enc)
#define MUPQ_crypto_kem_dec NAMESPACE(crypto_kem_dec)

#define printcycles(S, U) send_unsignedll((S), (U))

static int cmp_uint64_t(const void *a, const void *b) {
  return (int)((*((const uint64_t *)a)) - (*((const uint64_t *)b)));
}

int main(void) {
  unsigned char key_a[MUPQ_CRYPTO_BYTES], key_b[MUPQ_CRYPTO_BYTES];
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[MUPQ_CRYPTO_CIPHERTEXTBYTES];
  uint64_t cycles_kg[MUPQ_ITERATIONS], cycles_enc[MUPQ_ITERATIONS],
      cycles_dec[MUPQ_ITERATIONS];
  unsigned long long t0, t1;
  int i;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  for (i = 0; i < MUPQ_ITERATIONS; i++) {
    // Key-pair generation
    t0 = hal_get_time();
    MUPQ_crypto_kem_keypair(pk, sk);
    t1 = hal_get_time();
    cycles_kg[i] = t1 - t0;
  }

  for (i = 0; i < MUPQ_ITERATIONS; i++) {
    // Encapsulation
    t0 = hal_get_time();
    MUPQ_crypto_kem_enc(ct, key_a, pk);
    t1 = hal_get_time();
    cycles_enc[i] = t1 - t0;
  }

  for (i = 0; i < MUPQ_ITERATIONS; i++) {
    // Decapsulation
    t0 = hal_get_time();
    MUPQ_crypto_kem_dec(key_b, ct, sk);
    t1 = hal_get_time();
    cycles_dec[i] = t1 - t0;
  }

  qsort(cycles_kg, MUPQ_ITERATIONS, sizeof(uint64_t), cmp_uint64_t);
  qsort(cycles_enc, MUPQ_ITERATIONS, sizeof(uint64_t), cmp_uint64_t);
  qsort(cycles_dec, MUPQ_ITERATIONS, sizeof(uint64_t), cmp_uint64_t);

  printcycles("keypair cycles", cycles_kg[MUPQ_ITERATIONS >> 1]);
  printcycles("encaps cycles", cycles_enc[MUPQ_ITERATIONS >> 1]);
  printcycles("decaps cycles", cycles_dec[MUPQ_ITERATIONS >> 1]);

  hal_send_str("#");

  return 0;
}
