#include <stdint.h>
#include <string.h>
#include "aes.h"
#include "rijndael.h"

#ifdef PROFILE_HASHING
#include "hal.h"
extern unsigned long long hash_cycles;
#endif

void AES128_load_schedule(const uint8_t *key, uint8_t *schedule) {
#ifdef PROFILE_HASHING
  uint64_t t0 = hal_get_time();
#endif

  memcpy(schedule,key,16);

  AES_128_keyschedule(key,schedule+16);
#ifdef PROFILE_HASHING
  uint64_t t1 = hal_get_time();
  hash_cycles += (t1-t0);
#endif
}

void AES128_ECB_enc_sch(const uint8_t *plaintext, const size_t plaintext_len, const uint8_t *schedule, uint8_t *ciphertext) {
#ifdef PROFILE_HASHING
  uint64_t t0 = hal_get_time();
#endif
  for (size_t block = 0; block < plaintext_len / 16; block++) {
    AES_128_encrypt(schedule, plaintext + (16 * block), ciphertext + (16 * block));
  }
#ifdef PROFILE_HASHING
  uint64_t t1 = hal_get_time();
  hash_cycles += (t1-t0);
#endif
}

