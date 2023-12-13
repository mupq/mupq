
/**
 * @file parsing.c
 * @brief Implementation of parsing functions
 */

#include "parsing.h"
#include "parameters.h"

#include <gmp.h>
#include <string.h>
#include "arithmetic.h"
#include "data_structures.h"
#include "rank_unrank_table.h"
#include "symmetric.h"

/**
 * @brief store a 10 bit value in a byte array at bit position i*10
 *        must be called sequentially with increasing index
 *
 * @param sb  byte array pointer
 * @param i   position in the byte array
 * @param val 10 bit value to store
 */
static inline void store_10bit_in_bytearray(uint8_t* sb, int i, uint16_t val) {
    val &= 0x3FF;
    int k = (i / 4) * 5;
    switch (i % 4) {
        case 0:
            sb[k + 0] = val;
            sb[k + 1] = val >> 8;
            break;
        case 1:
            sb[k + 1] |= val << 2;
            sb[k + 2] = val >> 6;
            break;
        case 2:
            sb[k + 2] |= val << 4;
            sb[k + 3] = val >> 4;
            break;
        case 3:
            sb[k + 3] |= val << 6;
            sb[k + 4] = val >> 2;
            break;
        default:
            break;
    }
}

/**
 * @brief load a 10 bit value from a byte array at bit position i*10
 *
 * @param sb byte array pointer
 * @param i  position in the byte array
 * @return   uint16_t loaded value
 */
static inline uint16_t load_10bit_from_bytearray(const uint8_t* sb, int i) {
    int k = (i / 4) * 5;
    uint16_t val = 0;
    switch (i % 4) {
        case 0:
            val = sb[k + 0];
            val |= ((uint16_t)sb[k + 1] & 0x3) << 8;
            break;
        case 1:
            val = sb[k + 1] >> 2;
            val |= ((uint16_t)sb[k + 2] & 0xF) << 6;
            break;
        case 2:
            val = sb[k + 2] >> 4;
            val |= ((uint16_t)sb[k + 3] & 0x3F) << 4;
            break;
        case 3:
            val = sb[k + 3] >> 6;
            val |= ((uint16_t)sb[k + 4] & 0xFF) << 2;
            break;
        default:
            break;
    }

    return val;
}

#define PARAM_RANK_UNRANK_LEN_T \
    ((1 << (PARAM_RANK_UNRANK_K + 1)) - 1) /**< Length of list used in rank/unrank compression */

#if defined(__has_feature)
#if __has_feature(memory_sanitizer)
__attribute__((no_sanitize("memory")))
#endif
#if __has_feature(address_sanitizer)
__attribute__((no_sanitize("address")))
#endif
#endif
void sig_perk_perm_encode(const perm_t in_p, uint8_t* out_buff) {
    mpz_t mul, code;
    mpz_inits(mul, code, NULL);

    uint8_t T[PARAM_RANK_UNRANK_LEN_T] = {0};
    mpz_set_ui(code, 0);
    for (int i = 0; i < PARAM_N1; ++i) {
        uint8_t ctr = in_p[i];
        uint16_t node = (1 << PARAM_RANK_UNRANK_K) + in_p[i];
        for (int j = 0; j < PARAM_RANK_UNRANK_K; ++j) {
            if (node & 0x1) {
                ctr -= T[(node >> 1) << 1];
            }
            T[node] += 1;
            node = node >> 1;
        }
        T[node] += 1;
        mpz_mul_ui(mul, code, PARAM_N1 - i);
        mpz_add_ui(code, mul, ctr);
    }

    if (mpz_sgn(code) == 0) {
        memset(out_buff, 0, PARAM_PERM_COMPRESSION_BYTES);
    } else {
        size_t count;
        mpz_export(out_buff, &count, -1, 1, 1, 0, code);
        memset(out_buff + count, 0, PARAM_PERM_COMPRESSION_BYTES - count);
    }
    mpz_clears(mul, code, NULL);
}

#if defined(__has_feature)
#if __has_feature(memory_sanitizer)
__attribute__((no_sanitize("memory")))
#endif
#if __has_feature(address_sanitizer)
__attribute__((no_sanitize("address")))
#endif
#endif
int sig_perk_perm_decode(const uint8_t* in_buff, perm_t out_p) {
    mpz_t code, fact, tmp;
    mpz_inits(code, fact, tmp, NULL);

    mpz_import(code, PARAM_PERM_COMPRESSION_BYTES, -1, 1, 1, 0, in_buff);

    // validate the permutation encoding to be < n!
    mpz_set_str(fact, factorial[PARAM_N1], 62);
    if (mpz_cmp(fact, code) < 1) {
        mpz_clears(code, fact, tmp, NULL);
        return EXIT_FAILURE;
    }

    mpz_set_str(fact, factorial[PARAM_N1 - 1], 62);
    mpz_divmod(code, tmp, code, fact);

    out_p[0] = mpz_get_ui(code);

    for (int i = 1; i < PARAM_N1 - 1; ++i) {
        mpz_set_str(fact, factorial[PARAM_N1 - i - 1], 62);
        mpz_divmod(code, tmp, tmp, fact);
        out_p[i] = mpz_get_ui(code);
    }

    mpz_set_str(fact, factorial[0], 62);
    mpz_div(code, tmp, fact);
    out_p[PARAM_N1 - 1] = mpz_get_ui(code);

    mpz_clears(code, fact, tmp, NULL);

    uint8_t T[PARAM_RANK_UNRANK_LEN_T] = {0};
    for (int i = 0; i <= PARAM_RANK_UNRANK_K; ++i) {
        for (int j = 0; j < (1 << i); ++j) {
            T[((1 << i)) + j - 1] = 1 << (PARAM_RANK_UNRANK_K - i);
        }
    }

    for (int i = 0; i < PARAM_N1; ++i) {
        int digit = out_p[i];
        uint16_t node = 1;
        for (int j = 0; j < PARAM_RANK_UNRANK_K; ++j) {
            T[node] -= 1;
            node <<= 1;
            if (digit >= T[node]) {
                digit -= T[node];
                node += 1;
            }
        }
        T[node] = 0;
        out_p[i] = node - (1 << PARAM_RANK_UNRANK_K);
    }

    return EXIT_SUCCESS;
}

void sig_perk_private_key_to_bytes(uint8_t sk_bytes[PRIVATE_KEY_BYTES], const perk_private_key_t* sk) {
    memcpy(sk_bytes, sk->seed, SEED_BYTES);
    memcpy(sk_bytes + SEED_BYTES, sk->pk_bytes, PUBLIC_KEY_BYTES);
}

void sig_perk_private_key_from_bytes(perk_private_key_t* sk, const uint8_t sk_bytes[PRIVATE_KEY_BYTES]) {
    memcpy(sk->seed, sk_bytes, SEED_BYTES);
    memcpy(sk->pk_bytes, sk_bytes + SEED_BYTES, PUBLIC_KEY_BYTES);
    sig_perk_perm_set_random(sk->pi, sk->seed);
}

void sig_perk_public_key_to_bytes(uint8_t pk_bytes[PUBLIC_KEY_BYTES], const perk_public_key_t* pk) {
    memcpy(pk_bytes, pk->seed, SEED_BYTES);
    for (int i = 0; i < PARAM_M * PARAM_T; i++) {
        store_10bit_in_bytearray(pk_bytes + SEED_BYTES, i, pk->y[i / PARAM_M][i % PARAM_M]);
    }
}

int sig_perk_public_key_from_bytes(perk_public_key_t* pk, const uint8_t pk_bytes[PUBLIC_KEY_BYTES]) {
    memcpy(pk->seed, pk_bytes, SEED_BYTES);
    for (int i = 0; i < PARAM_M * PARAM_T; i++) {
        uint16_t y = load_10bit_from_bytearray(pk_bytes + SEED_BYTES, i);
        if (y >= PARAM_Q) {
            // y out of range
            return EXIT_FAILURE;
        }
        pk->y[i / PARAM_M][i % PARAM_M] = y;
    }

    sig_perk_prg_state_t prg;
    // initialize prg
    sig_perk_prg_init(&prg, PRG1, NULL, pk->seed);

    // Generate H and x
    sig_perk_mat_set_random(pk->H, &prg);
    sig_perk_vect1_set_random_list(pk->x, &prg);

    return EXIT_SUCCESS;
}

void sig_perk_challenges_from_bytes(challenge_t challenges[PARAM_TAU], const digest_t h1, const digest_t h2) {
    sig_perk_prg_state_t h1_prg_state;
    sig_perk_prg_state_t h2_prg_state;
    uint16_t tmp_kappa;
    uint16_t tmp_alpha;

    // generate first challenge
    sig_perk_prg_init(&h1_prg_state, PRG1, NULL, h1);
    for (int i = 0; i < PARAM_TAU; ++i) {
        uint16_t nonzero_check = 0;
        do {
            for (int j = 0; j < PARAM_T; ++j) {
                do {
                    sig_perk_prg(&h1_prg_state, (uint8_t*)&tmp_kappa, sizeof(tmp_kappa));
                    tmp_kappa = tmp_kappa & PARAM_Q_MASK;
                } while (tmp_kappa >= PARAM_Q);  // 0 <= tmp_kappa < PARAM_Q
                challenges[i].kappa[j] = tmp_kappa;
                nonzero_check |= tmp_kappa;
            }
        } while (!nonzero_check);
    }
    // generate second challenge
    sig_perk_prg_init(&h2_prg_state, PRG1, NULL, h2);
    for (int i = 0; i < PARAM_TAU; ++i) {
        sig_perk_prg(&h2_prg_state, (uint8_t*)&tmp_alpha, sizeof(tmp_alpha));
        tmp_alpha = (tmp_alpha & PARAM_N_MASK) + 1;  // 0 < tmp_alpha <= PARAM_N
        challenges[i].alpha = tmp_alpha;
    }
}

void sig_perk_signature_to_bytes(uint8_t sb[SIGNATURE_BYTES], const perk_signature_t* signature) {
    memcpy(sb, signature->salt, sizeof(salt_t));
    sb += sizeof(salt_t);
    memcpy(sb, signature->h1, sizeof(digest_t));
    sb += sizeof(digest_t);
    memcpy(sb, signature->h2, sizeof(digest_t));
    sb += sizeof(digest_t);

    for (int i = 0; i < PARAM_TAU; i++) {
        memcpy(sb, signature->responses[i].cmt_1_alpha, sizeof(cmt_t));
        sb += sizeof(cmt_t);
        memcpy(sb, signature->responses[i].z2_theta, sizeof(theta_t) * THETA_TREE_LEVELS);
        sb += sizeof(theta_t) * THETA_TREE_LEVELS;
    }

    for (int i = 0; i < PARAM_TAU * PARAM_N1; i++) {
        uint16_t z1 = signature->responses[i / PARAM_N1].z1[i % PARAM_N1];
        store_10bit_in_bytearray(sb, i, z1);
    }
    sb += (PARAM_TAU * PARAM_N1 * PARAM_Q_BITS + 7) / 8;

    for (int i = 0; i < (PARAM_TAU); i++) {
        sig_perk_perm_encode(signature->responses[i].z2_pi, sb);
        sb += PARAM_PERM_COMPRESSION_BYTES;
    }
}

#define Z1_USED_BITS (uint8_t)((PARAM_TAU * PARAM_N1 * PARAM_Q_BITS + 7) % 8 + 1)
static const uint8_t z1_unused_mask = (uint8_t)(((1U << (Z1_USED_BITS)) - 1) ^ 0xFFU);

int sig_perk_signature_from_bytes(perk_signature_t* signature, const uint8_t sb[SIGNATURE_BYTES]) {
    memcpy(signature->salt, sb, sizeof(salt_t));
    sb += sizeof(salt_t);
    memcpy(signature->h1, sb, sizeof(digest_t));
    sb += sizeof(digest_t);
    memcpy(signature->h2, sb, sizeof(digest_t));
    sb += sizeof(digest_t);

    for (int i = 0; i < PARAM_TAU; i++) {
        memcpy(signature->responses[i].cmt_1_alpha, sb, sizeof(cmt_t));
        sb += sizeof(cmt_t);
        memcpy(signature->responses[i].z2_theta, sb, sizeof(theta_t) * THETA_TREE_LEVELS);
        sb += sizeof(theta_t) * THETA_TREE_LEVELS;
    }

    for (int i = 0; i < PARAM_TAU * PARAM_N1; i++) {
        uint16_t z1 = load_10bit_from_bytearray(sb, i);
        if (z1 >= PARAM_Q) {
            // z1 out of range
            return EXIT_FAILURE;
        }
        signature->responses[i / PARAM_N1].z1[i % PARAM_N1] = z1;
    }
    sb += ((PARAM_TAU * PARAM_N1 * PARAM_Q_BITS + 7) / 8) - 1;

    // cppcheck-suppress knownConditionTrueFalse
    // cppcheck-suppress unmatchedSuppression
    if (sb[0] & z1_unused_mask) {
        // unused bits after the z1 != 0
        return EXIT_FAILURE;
    }

    sb += 1;

    for (int i = 0; i < PARAM_TAU; i++) {
        if (sig_perk_perm_decode(sb, signature->responses[i].z2_pi) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
        sb += PARAM_PERM_COMPRESSION_BYTES;
    }

    return EXIT_SUCCESS;
}