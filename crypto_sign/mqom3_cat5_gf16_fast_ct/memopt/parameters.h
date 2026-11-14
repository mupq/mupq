#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#define MQOM3_PARAM_SECURITY 256
#define MQOM3_PARAM_BASE_FIELD 4
#define MQOM3_PARAM_TRADEOFF 0
#define MQOM3_PARAM_OT_VARIANT 0

/* Fields conf: ref implementation */
#define FIELDS_REF
/* Rijndael conf: bitslice (actually underlying MUPQ implementation for cat1 with the MQOM3_FOR_MUPQ toggle) */
#define RIJNDAEL_BITSLICE
/* Options activated for memory optimization */
#define MEMORY_EFFICIENT_PIOP
#define FIELDS_BITSLICE_COMPOSITE
#define FIELDS_BITSLICE_PUBLIC_JUMP
#define MEMORY_EFFICIENT_BLC
#define MEMORY_EFFICIENT_KEYGEN
#define GGMTREE_NB_ENC_CTX_IN_MEMORY 0
#define SIGN_MEMOPT
#define VERIFY_MEMOPT
#define PRG_ONE_RIJNDAEL_CTX
#define SEED_COMMIT_MEMOPT
#define PIOP_NB_PARALLEL_REPETITIONS_SIGN 6
#define PIOP_NB_PARALLEL_REPETITIONS_VERIFY 4
#define SMALL_GGM_TREE_NB_SIMULTANEOUS_LEAVES_LOG 4
#define BLC_NB_LEAF_SEEDS_IN_PARALLEL 8
/* Specifically target MUPQ */
#define MQOM3_FOR_MUPQ

/* Do not mess with sections as the PQM4 framework uses them */
#define NO_EMBEDDED_SRAM_SECTION

#endif /* __PARAMETERS_H__ */
