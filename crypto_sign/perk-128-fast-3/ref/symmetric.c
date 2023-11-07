
/**
 * @file symmetric.c
 * @brief Implementation of the symmetric functions
 */

#include "symmetric.h"
#include "parameters.h"

#if (SECURITY_BYTES == 16)
#define shake_inc_init shake128_inc_init
#define shake_inc_absorb shake128_inc_absorb
#define shake_inc_finalize shake128_inc_finalize
#define shake_inc_squeeze shake128_inc_squeeze

#define sha3_inc_init sha3_256_inc_init
#define sha3_inc_absorb sha3_256_inc_absorb
#define sha3_inc_finalize sha3_256_inc_finalize
#elif (SECURITY_BYTES == 24) 
#define shake_inc_init shake256_inc_init
#define shake_inc_absorb shake256_inc_absorb
#define shake_inc_finalize shake256_inc_finalize
#define shake_inc_squeeze shake256_inc_squeeze

#define sha3_inc_init sha3_384_inc_init
#define sha3_inc_absorb sha3_384_inc_absorb
#define sha3_inc_finalize sha3_384_inc_finalize
#elif (SECURITY_BYTES == 32)
#define shake_inc_init shake256_inc_init
#define shake_inc_absorb shake256_inc_absorb
#define shake_inc_finalize shake256_inc_finalize
#define shake_inc_squeeze shake256_inc_squeeze

#define sha3_inc_init sha3_512_inc_init
#define sha3_inc_absorb sha3_512_inc_absorb
#define sha3_inc_finalize sha3_512_inc_finalize
#endif

void sig_perk_prg_init(sig_perk_prg_state_t *state, const uint8_t domain, const salt_t salt, const seed_t seed) {
    #if 1 
    shake_inc_init(state);
    if(salt != NULL) {
        shake_inc_absorb(state, salt, sizeof(salt_t));
    }
    if(seed != NULL) {
        shake_inc_absorb(state, seed, sizeof(seed_t));
    }
    shake_inc_absorb(state, &domain, 1);
    shake_inc_finalize(state);
    #else
    Keccak_HashInitialize_SHAKE(state);
    if (salt != NULL) {
        Keccak_HashUpdate(state, salt, sizeof(salt_t) * 8);
    }
    if (seed != NULL) {
        Keccak_HashUpdate(state, seed, sizeof(seed_t) * 8);
    }
    Keccak_HashUpdate(state, &domain, 1 * 8);
    Keccak_HashFinal(state, NULL);
    #endif
}

void sig_perk_prg(sig_perk_prg_state_t *state, uint8_t *output, size_t outlen) {
    #if 1
    shake_inc_squeeze(output, outlen, state);
    #else
    Keccak_HashSqueeze(state, output, outlen * 8);
    #endif
}

void sig_perk_hash_init(sig_perk_hash_state_t *state, const salt_t salt, const uint8_t *tau, const uint8_t *n) {
    #if 1
    sha3_inc_init(state);
    sha3_inc_absorb(state, salt, sizeof(salt_t));
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
        sha3_inc_absorb(state, counters, j);
    }

    #else
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
    #endif
}

void sig_perk_hash_update(sig_perk_hash_state_t *state, const uint8_t *message, const size_t message_size) {
    #if 1 
    sha3_inc_absorb(state, message, message_size);
    #else
    Keccak_HashUpdate(state, message, message_size * 8);
    #endif
}

void sig_perk_hash_final(sig_perk_hash_state_t *state, digest_t digest, const uint8_t domain) {
    #if 1
    sha3_inc_absorb(state, &domain, 1);
    sha3_inc_finalize(digest, state);
    #else
    Keccak_HashUpdate(state, &domain, 1 * 8);
    Keccak_HashFinal(state, digest);
    #endif
}
