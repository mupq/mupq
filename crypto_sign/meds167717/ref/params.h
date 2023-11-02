#ifndef PARAMS_H
#define PARAMS_H

#define MEDS_name "MEDS167717"

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

#define MEDS_s 6
#define MEDS_t 112
#define MEDS_w 66

#define MEDS_seed_tree_height 7
#define SEED_TREE_size 255
#define MEDS_max_path_len 62

#define MEDS_t_mask 0x0000007F
#define MEDS_t_bytes 1

#define MEDS_s_mask 0x00000007

#define MEDS_PK_BYTES 167717
#define MEDS_SK_BYTES 12444
#define MEDS_SIG_BYTES 165464

#endif

