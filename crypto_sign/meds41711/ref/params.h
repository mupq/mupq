#ifndef PARAMS_H
#define PARAMS_H

#define MEDS_name "MEDS41711"

#define MEDS_digest_bytes 32
#define MEDS_pub_seed_bytes 32
#define MEDS_sec_seed_bytes 32
#define MEDS_st_seed_bytes 24
#define MEDS_st_salt_bytes 32

#define MEDS_p 4093
#define GFq_t uint16_t
#define GFq_bits 12
#define GFq_bytes 2

#define MEDS_m 22
#define MEDS_n 22
#define MEDS_k 22

#define MEDS_s 4
#define MEDS_t 608
#define MEDS_w 26

#define MEDS_seed_tree_height 10
#define SEED_TREE_size 2047
#define MEDS_max_path_len 136

#define MEDS_t_mask 0x000003FF
#define MEDS_t_bytes 2

#define MEDS_s_mask 0x00000003

#define MEDS_PK_BYTES 41711
#define MEDS_SK_BYTES 4420
#define MEDS_SIG_BYTES 41080

#endif

