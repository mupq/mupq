/**
 *
 * Reference ISO-C11 Implementation of CROSS.
 *
 * @version 1.1 (March 2023)
 *
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#pragma once
#include <stdint.h>

#define RSDP
#define CATEGORY_5
#define SIG_SIZE
#define SHA3_HASH
#define SHAKE_CSPRNG


/******************************************************************************/
/****************************** RSDP Parameters *******************************/
/******************************************************************************/
#if defined(RSDP)

/* The same base field and restriction are employed for all categories of RSDP */
#define   Q (127)
#define   Z (  7)
/* single-register table representation of E, the value of g^7=1 is also
 * represented to avoid exponent renormalization*/
#define RESTR_G_TABLE ((uint64_t) (0x0140201008040201))
#define RESTR_G_GEN 2
#define FQ_ELEM uint8_t
#define FZ_ELEM uint8_t
#define FQ_DOUBLEPREC uint16_t
#define FQ_TRIPLEPREC uint32_t

/********************************* Category 1 *********************************/
#if defined(CATEGORY_1)
#define SEED_LENGTH_BYTES (16)
#define   N (127)
#define   K ( 76)

#if defined(SPEED)
#define   T (256)
#define   W (216)
#define POSITION_IN_FW_STRING_T uint16_t
/* determined via worst_case_seed_tree_param_computer.py script */
#define TREE_NODES_TO_STORE (108)

#elif defined(SIG_SIZE)
#define   T (871)
#define   W (848)
#define POSITION_IN_FW_STRING_T uint16_t
/* determined via worst_case_seed_tree_param_computer.py script */
#define TREE_NODES_TO_STORE (121)

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define SEED_LENGTH_BYTES (24)
#define   N (187)
#define   K (111)

#if defined(SPEED)
#define   T (256)
#define   W (160)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (136)

#elif defined(SIG_SIZE)
#define   T (1024)
#define   W (987)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (178)

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define SEED_LENGTH_BYTES (32)
#define   N (251)
#define   K (150)

#if defined(SPEED)
#define   T (512)
#define   W (432)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (215)

#elif defined(SIG_SIZE)
#define   T (1024)
#define   W (969)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (233)
#else
#error define optimization corner in Cmakelist
#endif

#else
#error define category for parameters
#endif



/******************************************************************************/
/****************************** RSDP(G) Parameters ****************************/
/******************************************************************************/
#elif defined(RSDPG)

/* The same base field and restriction are employed for all categories of RSDP */
#define   Q (509)
#define   Z (127)
/* Restricted subgroup generator */
#define RESTR_G_GEN 16
#define FZ_ELEM uint8_t
#define FZ_DOUBLEPREC uint16_t
#define FQ_ELEM uint16_t
#define FQ_DOUBLEPREC uint32_t
#define FQ_TRIPLEPREC uint32_t

/********************************* Category 1 *********************************/
#if defined(CATEGORY_1)
#define SEED_LENGTH_BYTES (16)
#define   N ( 42)
#define   K ( 23)
#define   M ( 24)

#if defined(SPEED)
#define   T (243)
#define   W (206)
#define POSITION_IN_FW_STRING_T uint8_t

/* determined via worst_case_seed_tree_param_computer.py script */
#define TREE_NODES_TO_STORE (101)

#elif defined(SIG_SIZE)
#define   T (871)
#define   W (850)
#define POSITION_IN_FW_STRING_T uint16_t
/* determined via worst_case_seed_tree_param_computer.py script */
#define TREE_NODES_TO_STORE (113)

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 3 *********************************/
#elif defined(CATEGORY_3)
#define SEED_LENGTH_BYTES (24)
#define   N ( 63)
#define   K ( 35)
#define   M ( 36)

#if defined(SPEED)
#define   T (255)
#define   W (176)
#define POSITION_IN_FW_STRING_T uint8_t
#define TREE_NODES_TO_STORE (134)

#elif defined(SIG_SIZE)
#define   T (949)
#define   W (914)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (167)

#else
#error define optimization corner in Cmakelist
#endif

/********************************* Category 5 *********************************/
#elif defined(CATEGORY_5)
#define SEED_LENGTH_BYTES (32)
#define   N ( 87)
#define   K ( 47)
#define   M ( 48)

#if defined(SPEED)
#define   T (356)
#define   W (257)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (183)

#elif defined(SIG_SIZE)
#define   T (949)
#define   W (897)
#define POSITION_IN_FW_STRING_T uint16_t
#define TREE_NODES_TO_STORE (218)
#else
#error define optimization corner in Cmakelist
#endif

#else
#error define category for parameters
#endif

#else
#error define either RSDP or RSDPG
#endif


/***************** Derived parameters *****************************************/
/*length of the output of the cryptographic hash, in bytes */
#define HASH_DIGEST_LENGTH (2*SEED_LENGTH_BYTES)
#define SALT_LENGTH_BYTES (SEED_LENGTH_BYTES)
#define NUM_LEAVES_OF_MERKLE_TREE (T)
#define NUM_NODES_OF_MERKLE_TREE (2*NUM_LEAVES_OF_MERKLE_TREE-1)
#define SEED_TREE_MAX_PUBLISHED_BYTES (TREE_NODES_TO_STORE*SEED_LENGTH_BYTES)
#define MERKLE_TREE_MAX_PUBLISHED_BYTES (TREE_NODES_TO_STORE*HASH_DIGEST_LENGTH)

/******************************** Helper macros *******************************/
#define ROUND_UP(amount, round_amt) ( ((amount+round_amt-1)/round_amt)*round_amt )

#define IS_REPRESENTABLE_IN_D_BITS(D, N)                \
  (((unsigned long) N>=(1UL << (D-1)) && (unsigned long) N<(1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (15                                     \
                 + IS_REPRESENTABLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(16, N)    \
                 )                                      \
   )
#define NUM_LEAVES_OF_SEED_TREE ( 1UL << BITS_TO_REPRESENT(T) )
#define NUM_NODES_OF_SEED_TREE ( 2*NUM_LEAVES_OF_SEED_TREE-1 )
#define LOGT BITS_TO_REPRESENT(T)

  /* bitpacking an n-elements vector of  num_bits_for_q-1 bits long values
 * will pack 8 values in num_bits_for_q-1 bytes exactly, leaving the remaining
 * N % 8 as a tail */
#define DENSELY_PACKED_FQ_VEC_SIZE (N/8)*BITS_TO_REPRESENT(Q-1) + \
                                   ROUND_UP( ((N%8)*BITS_TO_REPRESENT(Q-1)),8)/8
#define DENSELY_PACKED_FQ_SYN_SIZE ((N-K)/8)*BITS_TO_REPRESENT(Q-1) + \
                                   ROUND_UP( (((N-K)%8)*BITS_TO_REPRESENT(Q-1)),8)/8
#define DENSELY_PACKED_FZ_VEC_SIZE (N/8)*BITS_TO_REPRESENT(Z-1) + \
                                   ROUND_UP( ((N%8)*BITS_TO_REPRESENT(Z-1)),8)/8
#ifdef RSDPG
#define DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE (M/8)*BITS_TO_REPRESENT(Z-1) + \
                                          ROUND_UP( ((M%8)*BITS_TO_REPRESENT(Z-1)),8)/8
#endif
