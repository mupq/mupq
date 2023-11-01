#ifndef UTILS_H
#define UTILS_H

#include "params.h"
#include "matrixmod.h"

void XOF(uint8_t **buf, size_t *length, const uint8_t *seed, size_t seed_len, int num);

GFq_t rnd_GF(keccak_state *shake);

void rnd_sys_mat(pmod_mat_t *M, int M_r, int M_c, const uint8_t *seed, size_t seed_len);

void rnd_inv_matrix(pmod_mat_t *M, int M_r, int M_c, uint8_t *seed, size_t seed_len);

int parse_hash(uint8_t *digest, int digest_len, uint8_t *h, int len_h);

int solve(pmod_mat_t *A, pmod_mat_t *B_inv, pmod_mat_t *G0prime, GFq_t Amm);

void pi(pmod_mat_t *Gout, pmod_mat_t *A, pmod_mat_t *B, pmod_mat_t *G);

#endif

