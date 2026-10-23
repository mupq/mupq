#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#define MQOM3_PARAM_SECURITY 256
#define MQOM3_PARAM_BASE_FIELD 4
#define MQOM3_PARAM_TRADEOFF 1
#define MQOM3_PARAM_OT_VARIANT 0

/* Fields conf: ref implementation */
#define FIELDS_REF
/* Rijndael conf: bitslice (actually underlying MUPQ implementation for cat1 with the MQOM3_FOR_MUPQ toggle) */
#define RIJNDAEL_BITSLICE
/* Options activated for memory optimization */
#define PIOP_BITSLICE
#define FIELDS_BITSLICE_COMPOSITE
#define FIELDS_BITSLICE_PUBLIC_JUMP
#define MEMORY_EFFICIENT_BLC
#define MEMORY_EFFICIENT_KEYGEN
#define SMALL_GGM_TREE_NB_SIMULTANEOUS_LEAVES_LOG 7
#define BLC_NB_LEAF_SEEDS_IN_PARALLEL 64
#define BLC_SEEDCOMMIT_CACHE
#define BLC_SEEDEXPAND_CACHE
/* Specifically target MUPQ */
#define MQOM3_FOR_MUPQ

/* Do not mess with sections as the PQM4 framework uses them */
#define NO_EMBEDDED_SRAM_SECTION

#endif /* __PARAMETERS_H__ */
