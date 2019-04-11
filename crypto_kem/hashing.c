#include "api.h"
#include "hal.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>

static void printcycles(const char *s, unsigned long long c)
{
  char outs[32];
  hal_send_str(s);
  snprintf(outs,sizeof(outs),"%llu\n",c);
  hal_send_str(outs);
}

unsigned long long hash_cycles;

int main(void)
{
  unsigned char key_a[CRYPTO_BYTES], key_b[CRYPTO_BYTES];
  unsigned char sk[CRYPTO_SECRETKEYBYTES];
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
  unsigned long long t0, t1;

  hal_setup(CLOCK_BENCHMARK);

  hal_send_str("==========================");

  // Key-pair generation
  hash_cycles = 0;
  t0 = hal_get_time();
  crypto_kem_keypair(pk, sk);
  t1 = hal_get_time();
  printcycles("keypair cycles:", t1-t0);
  printcycles("keypair hash cycles:", hash_cycles);

  // Encapsulation
  hash_cycles = 0;
  t0 = hal_get_time();
  crypto_kem_enc(ct, key_a, pk);
  t1 = hal_get_time();
  printcycles("encaps cycles: ", t1-t0);
  printcycles("encaps hash cycles: ", hash_cycles);

  // Decapsulation
  hash_cycles = 0;
  t0 = hal_get_time();
  crypto_kem_dec(key_b, ct, sk);
  t1 = hal_get_time();
  printcycles("decaps cycles: ", t1-t0);
  printcycles("decaps hash cycles: ", hash_cycles);

  if (memcmp(key_a, key_b, CRYPTO_BYTES)) {
    hal_send_str("ERROR KEYS\n");
  }
  else {
    hal_send_str("OK KEYS\n");
  }

  hal_send_str("#");
  while(1);
  return 0;
}
