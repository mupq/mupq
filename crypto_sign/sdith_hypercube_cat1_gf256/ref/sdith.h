#ifndef SIMPLE_EXAMPLE_SDITH_H
#define SIMPLE_EXAMPLE_SDITH_H

#include <stdint.h>
#include <stdio.h>

#include "rng.h"
#include "param.h"
#include "types.h"

#define L_SUFFIX secl1
// primary parameters
#define PARAM_q 256   // cardinal of F_poly
#if (PARAM_fpoint_size == 3)
#define fpoints_mul_ct gf2p24_mul_ct
#define fpoints_mul gf2p24_mul
#define fpoints_dlog gf2p24_dlog
#define fpoints_dexp gf2p24_dexp
#define fpoints_dlog_pow_log gf2p24_log_pow
#define fpoints_dlog_mul_log_log gf2p24_log_mul_log_log
#elif (PARAM_fpoint_size == 4)
#define fpoints_mul_ct gf2p32_mul_ct
#define fpoints_mul gf2p32_mul
#define fpoints_dlog gf2p32_dlog
#define fpoints_dexp gf2p32_dexp
#define fpoints_dlog_pow_log gf2p32_dlog_pow_l32
#define fpoints_dlog_mul_log_log gf2p32_dlog_mul_l32_l32
#else
#error "Field extension not supported"
#endif

#define LSEC_NAME(NAME) NAME##_secl1

#define sdith_compressed_pubkey_t LSEC_NAME(sdith_compressed_pubkey)
#define sdith_compressed_key_t LSEC_NAME(sdith_compressed_key)
#define sdith_full_pubkey_t LSEC_NAME(sdith_full_pubkey)
#define sdith_full_key_t LSEC_NAME(sdith_full_key)
#define sdith_params_t LSEC_NAME(sdith_params)

/** @brief compressed public key type */
typedef struct {
  seed_t H_a_seed;
  fsd_t y[PAR_y_size];
} sdith_compressed_pubkey_t;

/** @brief compressed secret key type */
typedef struct {
  seed_t m_seed;
} sdith_compressed_key_t;

/** @brief expanded public key type */
typedef struct {
  sdith_compressed_pubkey_t compressed_pubkey;
  ha_slice_t H_a[PAR_ha_nslice][PARAM_k]; // H_a is given by slices of 128 columns
} sdith_full_pubkey_t;

/** @brief expanded secret key type */
typedef struct {
  sdith_compressed_pubkey_t compressed_pubkey;
  fsd_t s_A[PARAM_k];
  fpoly_t q_poly[PARAM_d][PAR_wd];
  fpoly_t p_poly[PARAM_d][PAR_wd];
} sdith_full_key_t;

// TODO(stevenyue): Move types into an internal header.
typedef struct aux_share_struct {
  uint8_t s_A[PARAM_k];
  uint8_t q_poly[PARAM_d][PAR_wd];
  uint8_t p_poly[PARAM_d][PAR_wd];
  uint32_t c[PARAM_d][PARAM_t];
} aux_share_t;

typedef struct signature_struct {
  salt_t salt;
  hash_t h2;
  seed_t tree_prg_seeds[PARAM_tau][PARAM_D];
  commit_t com[PARAM_tau];
  aux_share_t aux[PARAM_tau];
  fpoints_c_t compressed_alpha[PARAM_tau][PARAM_d][PARAM_t];
  fpoints_c_t compressed_beta[PARAM_tau][PARAM_d][PARAM_t];
} signature_t;

void serialize_compressed_pk(FILE* F, sdith_compressed_pubkey_t const* pk);
void deserialize_pk(FILE* F, sdith_compressed_pubkey_t* pk);
void serialize_sk(FILE* F, sdith_compressed_key_t const* sk);
void deserialize_sk(FILE* F, sdith_compressed_key_t* sk);

void field_init();
void keygen(sdith_compressed_pubkey_t* pk, sdith_compressed_key_t* sk);
void uncompress_key(sdith_compressed_pubkey_t const* pk, sdith_compressed_key_t const* sk, sdith_full_pubkey_t* u_pk,
                    sdith_full_key_t* u_sk);
void uncompress_pubkey(sdith_compressed_pubkey_t const* pk, sdith_full_pubkey_t* u_pk);
void sign(sdith_full_pubkey_t const* pk, sdith_full_key_t const* sk, void const* msg, int msgBytes,
          void* sig, int* sigBytes);

typedef struct sdith_ctx_struct sdith_ctx_t;
sdith_ctx_t* new_sdith_ctx();
void free_sdith_ctx(sdith_ctx_t* ctx);
void sign_offline(sdith_ctx_t* ctx, sdith_full_pubkey_t __attribute__((unused)) const* pk, sdith_full_key_t const* sk);
void sign_online(sdith_ctx_t* ctx, sdith_full_pubkey_t const* pk, sdith_full_key_t __attribute__((unused)) const* sk,
                 uint8_t const* msg, int msgBytes, uint8_t* sig, int* sigBytes);
int verify(sdith_full_pubkey_t const* pk, void const* msg, int msgBytes, void const* sig);

int sig_size();

const char* hexmem(const void* addr, uint64_t nbytes);

const char* hashmem(const void* addr, uint64_t nbytes);

#endif  // SIMPLE_EXAMPLE_SDITH_H
