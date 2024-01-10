#include "witness.h"
#include "mpc.h"

void hash_update_instance(hash_context* ctx, const instance_t* inst) {
    hash_update(ctx, inst->seed, PARAM_SEED_SIZE);
    hash_update(ctx, inst->y, PARAM_m);
}

void serialize_instance(uint8_t* buf, const instance_t* inst) {
    memcpy(buf, inst->seed, PARAM_SEED_SIZE);
    vec_compress(buf + PARAM_SEED_SIZE, inst->y, PARAM_m);
}

void deserialize_instance(instance_t* inst, const uint8_t* buf) {
    memcpy(inst->seed, buf, PARAM_SEED_SIZE);
    vec_decompress(inst->y, buf + PARAM_SEED_SIZE, PARAM_m);
}

void serialize_instance_solution(uint8_t* buf, const solution_t* sol) {
    vec_compress(buf, sol->x, PARAM_n);
}

void deserialize_instance_solution(solution_t* sol, const uint8_t* buf) {
    vec_decompress(sol->x, buf, PARAM_n);
}

void compress_plain_broadcast(uint8_t* buf, const mpc_broadcast_t* plain_br) {
    vec_compress(buf, plain_br->alpha, sizeof(plain_br->alpha));
}

void uncompress_plain_broadcast(mpc_broadcast_t* plain_br, const uint8_t* buf) {
    vec_decompress(plain_br->alpha, buf, sizeof(plain_br->alpha));
    vec_set_zero(plain_br->v, sizeof(plain_br->v));
}
