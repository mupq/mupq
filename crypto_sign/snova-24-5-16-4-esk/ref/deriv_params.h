#ifndef CONFIG_H
#define CONFIG_H

#include "params.h"

#define seed_length_public 16
#define seed_length_private 32
#define seed_length (seed_length_public + seed_length_private)

#define n_SNOVA  (v_SNOVA + o_SNOVA)
#define m_SNOVA  (o_SNOVA)
#define lsq_SNOVA  (l_SNOVA * l_SNOVA)

#define GF16s_sha3_public (lsq_SNOVA * (2 * (lsq_SNOVA + l_SNOVA) + m_SNOVA * (n_SNOVA * n_SNOVA - m_SNOVA * m_SNOVA)))
#define GF16s_sha3_private (v_SNOVA*o_SNOVA*l_SNOVA)
#define bytes_sha3_public ((GF16s_sha3_public + 1) >> 1)
#define bytes_sha3_private (seed_length_public + ((GF16s_sha3_private + 1) >> 1))

#define GF16s_hash (o_SNOVA * lsq_SNOVA)
#define GF16s_signature (n_SNOVA * lsq_SNOVA)
#define bytes_hash ((GF16s_hash + 1) >> 1)

#define rank (l_SNOVA)
#define sq_rank (rank * rank)

#define bytes_pk (seed_length_public + ((m_SNOVA * o_SNOVA * o_SNOVA * lsq_SNOVA + 1) >> 1))
#define bytes_sk (((lsq_SNOVA * (4 * lsq_SNOVA + (m_SNOVA * v_SNOVA * v_SNOVA + m_SNOVA * v_SNOVA * o_SNOVA + m_SNOVA * o_SNOVA * v_SNOVA) + v_SNOVA * o_SNOVA) + 1) >> 1) + seed_length_public + seed_length_private)

#define bytes_signature ((GF16s_signature + 1) >> 1)
#define bytes_salt 16

#define GF16s_prng_public (lsq_SNOVA * (2 * (lsq_SNOVA + l_SNOVA) + m_SNOVA * (n_SNOVA * n_SNOVA - m_SNOVA * m_SNOVA)))
#define GF16s_prng_private (v_SNOVA * o_SNOVA * l_SNOVA)
#define bytes_prng_public ((GF16s_prng_public + 1) >> 1)
#define bytes_prng_private ((GF16s_prng_private + 1) >> 1)

#endif
