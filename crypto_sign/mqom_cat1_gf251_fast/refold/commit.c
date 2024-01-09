#include "commit.h"
#include "hash.h"

// Commitment x1
static void commit_string(uint8_t* digest, const uint8_t* str, size_t len, const uint8_t* salt, uint16_t e, uint16_t i) {
    hash_context ctx;
    hash_init_prefix(&ctx, HASH_PREFIX_COMMITMENT);
    hash_update(&ctx, salt, PARAM_SALT_SIZE);
    hash_update_uint16_le(&ctx, e);
    hash_update_uint16_le(&ctx, i);
    hash_update(&ctx, str, len);
    hash_final(&ctx, digest);
}

void commit_share(uint8_t* digest, const mpc_share_t* share, const uint8_t* salt, uint16_t e, uint16_t i) {
    commit_string(digest, (uint8_t*) share, PARAM_SHARE_SIZE, salt, e, i);
}

void commit_seed(uint8_t* digest, const uint8_t* seed, const uint8_t* salt, uint16_t e, uint16_t i) {
    commit_string(digest, seed, PARAM_SEED_SIZE, salt, e, i);
}

void commit_seed_and_wit(uint8_t* digest, const uint8_t* seed, const mpc_wit_t* wit, const uint8_t* salt, uint16_t e, uint16_t i) {
    hash_context ctx;
    hash_init_prefix(&ctx, HASH_PREFIX_COMMITMENT);
    hash_update(&ctx, salt, PARAM_SALT_SIZE);
    hash_update_uint16_le(&ctx, e);
    hash_update_uint16_le(&ctx, i);
    hash_update(&ctx, seed, PARAM_SEED_SIZE);
    hash_update(&ctx, (uint8_t*) wit, PARAM_WIT_SIZE);
    hash_final(&ctx, digest);
}

void commit_hint(uint8_t* digest, const mpc_hint_t* hint, const uint8_t* salt, uint16_t e, uint16_t i) {
    commit_string(digest, (uint8_t*) hint, PARAM_HINT_SIZE, salt, e, i);
}


// Commitment x4
static void commit_string_x4(uint8_t** digest, uint8_t const* const* str, size_t len, const uint8_t* salt, uint16_t e, const uint16_t* i) {
    hash_context_x4 ctx;
    hash_init_prefix_x4(&ctx, HASH_PREFIX_COMMITMENT);
    hash_update_x4_1(&ctx, salt, PARAM_SALT_SIZE);
    hash_update_x4_uint16_le(&ctx, e);
    hash_update_x4_uint16s_le(&ctx, i);
    hash_update_x4(&ctx, str, len);
    hash_final_x4(&ctx, digest);
}

void commit_share_x4(uint8_t** digest,  mpc_share_t const* const* share, const uint8_t* salt, uint16_t e, const uint16_t* i) {
    commit_string_x4(digest, (uint8_t const* const*) share, sizeof(mpc_share_t), salt, e, i);
}

void commit_seed_x4(uint8_t** digest, uint8_t const* const* seed, const uint8_t* salt, uint16_t e, const uint16_t* i) {
    commit_string_x4(digest, seed, PARAM_SEED_SIZE, salt, e, i);
}

void commit_seed_and_wit_x4(uint8_t** digest, uint8_t const* const* seed, mpc_wit_t const* const* wit, const uint8_t* salt, uint16_t e, const uint16_t* i) {
    hash_context_x4 ctx;
    hash_init_prefix_x4(&ctx, HASH_PREFIX_COMMITMENT);
    hash_update_x4_1(&ctx, salt, PARAM_SALT_SIZE);
    hash_update_x4_uint16_le(&ctx, e);
    hash_update_x4_uint16s_le(&ctx, i);
    hash_update_x4(&ctx, seed, PARAM_SEED_SIZE);
    hash_update_x4(&ctx, (uint8_t const* const*) wit, PARAM_WIT_SIZE);
    hash_final_x4(&ctx, digest);
}

void commit_hint_x4(uint8_t** digest, mpc_hint_t const* const* hint, const uint8_t* salt, uint16_t e, const uint16_t* i) {
    commit_string_x4(digest, (uint8_t const* const*) hint, PARAM_HINT_SIZE, salt, e, i);
}
