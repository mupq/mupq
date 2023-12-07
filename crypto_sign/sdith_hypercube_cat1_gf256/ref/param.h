#ifndef PARAM_H
#define PARAM_H

#ifndef SUPERCOP
#ifndef CAT_1
#ifndef CAT_3
#ifndef CAT_5
#define CAT_1
#endif
#endif
#endif
#endif

#ifdef CAT_1
#define PARAM_m 230   // code length
#define PARAM_k 126   // code dimension
#define PARAM_w 79    // hammimg weight
#define PARAM_d 1     // twist
#define PARAM_t 3     // number of challenge points per iteration
#define PARAM_tau 17  // number of parallel repetitions
#define PARAM_D   8   // hypercube dimension (2^D leaf parties)
#define PARAM_seed_size 16
#define PARAM_rho_size 16
#define PARAM_salt_size 32
#define PARAM_commit_size 32
#define PARAM_hash_size 32
#define PARAM_fpoint_size 4
#endif

#ifdef CAT_3
#define PARAM_m 352   // code length
#define PARAM_k 193   // code dimension
#define PARAM_w 120   // hammimg weight
#define PARAM_d 2     // twist
#define PARAM_t 3     // number of challenge points per iteration
#define PARAM_tau 26  // number of parallel repetitions
#define PARAM_D   8   // hypercube dimension (2^D leaf parties)
#define PARAM_seed_size 24
#define PARAM_rho_size 24
#define PARAM_salt_size 48
#define PARAM_commit_size 48
#define PARAM_hash_size 48
#define PARAM_fpoint_size 4
#endif

#ifdef CAT_5
#define PARAM_m 480   // code length
#define PARAM_k 278   // code dimension
#define PARAM_w 150   // hammimg weight
#define PARAM_d 2     // twist
#define PARAM_t 4     //
#define PARAM_tau 34  // number of parallel repetitions
#define PARAM_D   8   // hypercube dimension (2^D leaf parties)
#define PARAM_seed_size 32
#define PARAM_rho_size 32
#define PARAM_salt_size 64
#define PARAM_commit_size 64
#define PARAM_hash_size 64
#define PARAM_fpoint_size 4
#endif

// feature defines
// #define BENCHMARK
#define IDS_3_ROUND
// #define PROOF_OF_WORK
#define FULL_TREE

// secondary deduced parameters
#define PAR_wd (PARAM_w/PARAM_d)             // hamming weight per twist
#define PAR_y_size (PARAM_m-PARAM_k)         // code codimension
#define PAR_ha_nslice ((PAR_y_size+127u)/128u) // number of pk slices
#define PAR_md (PARAM_m/PARAM_d)             // dimension per twist
#define PAR_POW_iter (PARAM_D)               // number of proof of work iterations
#define PAR_fpoint_mask ((1UL<<(PARAM_fpoint_size*8u))-1u)

#endif // PARAM_H
