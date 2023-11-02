#ifndef PARAMS_H
#define PARAMS_H

#define MEDS_name "MEDS134180"

#define MEDS_digest_bytes 32
#define MEDS_pub_seed_bytes 32
#define MEDS_sec_seed_bytes 32
#define MEDS_st_seed_bytes 32
#define MEDS_st_salt_bytes 32

#define MEDS_p 2039
#define GFq_t uint16_t
#define GFq_bits 11
#define GFq_bytes 2

#define MEDS_m 30
#define MEDS_n 30
#define MEDS_k 30

#define MEDS_s 5
#define MEDS_t 192
#define MEDS_w 52

#define MEDS_seed_tree_height 8
#define SEED_TREE_size 511
#define MEDS_max_path_len 116

#define MEDS_t_mask 0x000000FF
#define MEDS_t_bytes 1

#define MEDS_s_mask 0x00000007

#define MEDS_PK_BYTES 134180
#define MEDS_SK_BYTES 9968
#define MEDS_SIG_BYTES 132528

#endif

