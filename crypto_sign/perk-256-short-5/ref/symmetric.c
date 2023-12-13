
/**
 * @file symmetric.c
 * @brief Implementation of the symmetric functions
 */

#include "symmetric.h"
#include "parameters.h"

#if (SECURITY_BYTES == 16)
#define Keccak_HashInitialize_SHAKE Keccak_HashInitialize_SHAKE128  // SHAKE128
#define Keccak_HashInitialize_SHA3  Keccak_HashInitialize_SHA3_256
#elif (SECURITY_BYTES == 24)
#define Keccak_HashInitialize_SHAKE Keccak_HashInitialize_SHAKE256  // SHAKE256
#define Keccak_HashInitialize_SHA3  Keccak_HashInitialize_SHA3_384
#elif (SECURITY_BYTES == 32)
#define Keccak_HashInitialize_SHAKE Keccak_HashInitialize_SHAKE256  // SHAKE256
#define Keccak_HashInitialize_SHA3  Keccak_HashInitialize_SHA3_512
#endif

void sig_perk_prg_init(sig_perk_prg_state_t *state, const uint8_t domain, const salt_t salt, const seed_t seed) {
    Keccak_HashInitialize_SHAKE(state);
    if (salt != NULL) {
        Keccak_HashUpdate(state, salt, sizeof(salt_t) * 8);
    }
    if (seed != NULL) {
        Keccak_HashUpdate(state, seed, sizeof(seed_t) * 8);
    }
    Keccak_HashUpdate(state, &domain, 1 * 8);
    Keccak_HashFinal(state, NULL);
}

void sig_perk_prg(sig_perk_prg_state_t *state, uint8_t *output, size_t outlen) {
    Keccak_HashSqueeze(state, output, outlen * 8);
}

void sig_perk_hash_init(sig_perk_hash_state_t *state, const salt_t salt, const uint8_t *tau, const uint8_t *n) {
    Keccak_HashInitialize_SHA3(state);
    Keccak_HashUpdate(state, salt, sizeof(salt_t) * 8);

    uint8_t counters[2];
    int j = 0;
    if (tau != NULL) {
        counters[j] = *tau;
        j++;
    }
    if (n != NULL) {
        counters[j] = *n;
        j++;
    }
    if (j != 0) {
        Keccak_HashUpdate(state, counters, j * 8);
    }
}

void sig_perk_hash_update(sig_perk_hash_state_t *state, const uint8_t *message, const size_t message_size) {
    Keccak_HashUpdate(state, message, message_size * 8);
}

void sig_perk_hash_final(sig_perk_hash_state_t *state, digest_t digest, const uint8_t domain) {
    Keccak_HashUpdate(state, &domain, 1 * 8);
    Keccak_HashFinal(state, digest);
}
