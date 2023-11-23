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

#include "parameters.h"
#include "sha2.h"
#include "sha3.h"
#include "aes256_ctr_drbg.h"

/************************* CSPRNG ********************************/
#if defined(AES_CTR_CSPRNG)
#define CSPRNG_STATE_T AES256_CTR_DRBG_struct

static inline
void initialize_csprng(CSPRNG_STATE_T * const csprng_state,
                       const unsigned char * const seed,
                       const uint32_t seed_len_bytes) {
    initialize_aes_csprng(csprng_state, seed, seed_len_bytes);
}

static inline
void csprng_randombytes(unsigned char * const x,
                        unsigned long long xlen,
                        CSPRNG_STATE_T * const csprng_state){
   aes_csprng_randombytes(x,xlen,csprng_state);
}


#elif defined(SHAKE_CSPRNG)

#define CSPRNG_STATE_T SHAKE_STATE_STRUCT
/* initializes a CSPRNG, given the seed and a state pointer */
static inline
void initialize_csprng(CSPRNG_STATE_T * const csprng_state,
                       const unsigned char * const seed,
                       const uint32_t seed_len_bytes) {
   // the second parameter is the security level of the SHAKE instance
   xof_shake_init(csprng_state, SEED_LENGTH_BYTES*8);
   xof_shake_update(csprng_state,seed,seed_len_bytes);
   xof_shake_final(csprng_state);
} /* end initialize_csprng */

/* extracts xlen bytes from the CSPRNG, given the state */
static inline
void csprng_randombytes(unsigned char * const x,
                        unsigned long long xlen,
                        CSPRNG_STATE_T * const csprng_state){
   xof_shake_extract(csprng_state,x,xlen);
}

#else
#error please define either AES_CTR_CSPRNG or SHAKE_CSPRNG
#endif

/******************************************************************************/

/* global csprng state employed to have deterministic randombytes for testing */
extern CSPRNG_STATE_T platform_csprng_state;
/* extracts xlen bytes from the global CSPRNG */
static inline
void randombytes(unsigned char * x,
                 unsigned long long xlen) {
   csprng_randombytes(x,xlen,&platform_csprng_state);
}

/************************* HASH functions ********************************/

#if defined(SHA3_HASH)
void sha3_hash(uint8_t digest[HASH_DIGEST_LENGTH],
               const unsigned char *const m,
               const uint64_t mlen);


/* Opaque algorithm agnostic hash call */
static inline
void hash(uint8_t digest[HASH_DIGEST_LENGTH],
          const unsigned char *const m,
          const uint64_t mlen){
     sha3_hash(digest,m,mlen);
}

#elif defined(SHA2_HASH)
void sha2_hash(uint8_t digest[HASH_DIGEST_LENGTH],
               const unsigned char *const m,
               const uint64_t mlen);

/* Opaque algorithm agnostic hash call */
static inline
void hash(uint8_t digest[HASH_DIGEST_LENGTH],
          const unsigned char *const m,
          const uint64_t mlen){
     sha2_hash(digest,m,mlen);
}

#else
#error please define either SHA3_HASH or SHA2_HASH
#endif

/********************** CSPRNG Sampling functions helpers ********************/

static inline
FQ_ELEM fq_star_rnd_state(CSPRNG_STATE_T * const csprng_state)
{
   const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_TO_REPRESENT(Q-2) ) - 1;
   FQ_ELEM rnd_value;
   do {
      csprng_randombytes((unsigned char *) &rnd_value,
                         sizeof(FQ_ELEM),
                         csprng_state);
      rnd_value = mask & rnd_value;
   } while (rnd_value > Q-2);

   return rnd_value+1;
} /* end fq_star_rnd_state */


static inline
FQ_ELEM fq_rnd_state(CSPRNG_STATE_T * const csprng_state)
{
   const FQ_ELEM mask = ( (FQ_ELEM) 1 << BITS_TO_REPRESENT(Q-1)) - 1;
   FQ_ELEM rnd_value;
   do {
      csprng_randombytes((unsigned char *) &rnd_value,
                         sizeof(FQ_ELEM),
                         csprng_state);
      rnd_value = mask & rnd_value;
   } while (rnd_value >= Q);

   return rnd_value;
} /* end fq_rnd_state */

static inline
FZ_ELEM zz_rnd_state(CSPRNG_STATE_T * const csprng_state)
{
   const FZ_ELEM mask = ( (FZ_ELEM) 1 << BITS_TO_REPRESENT(Z)) - 1;
   FZ_ELEM rnd_value;
   do {
      csprng_randombytes((unsigned char *) &rnd_value,
                         sizeof(FZ_ELEM),
                         csprng_state);
      rnd_value = mask & rnd_value;
   } while (rnd_value >= Z);
   return rnd_value;
} /* end zz_rnd_state */

/***************** Specialized CSPRNGs for non binary domains *****************/

/* CSPRNG sampling fixed weight strings */
void expand_digest_to_fixed_weight(uint8_t fixed_weight_string[T],
                                   const uint8_t digest[HASH_DIGEST_LENGTH]);

static inline
void CSPRNG_zz_vec(FZ_ELEM res[N],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<N; i++){
        res[i] = zz_rnd_state(csprng_state);
    }
}



static inline
void CSPRNG_fq_vec(FQ_ELEM res[N],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<N; i++){
        res[i] = fq_rnd_state(csprng_state);
    }
}

static inline
void CSPRNG_fq_vec_beta(FQ_ELEM res[T],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<T; i++){
        res[i] = fq_star_rnd_state(csprng_state);
    }
}

static inline
void CSPRNG_fq_mat(FQ_ELEM res[N-K][K],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<N-K; i++){
        for(int j = 0; j<K; j++){
            res[i][j] = fq_rnd_state(csprng_state);
        }
    }
}

#if defined(RSDPG)
static inline
void CSPRNG_zz_inf_w(FZ_ELEM res[M],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<M; i++){
        res[i] = zz_rnd_state(csprng_state);
    }
}

static inline
void CSPRNG_fz_mat(FZ_ELEM res[N-M][M],
                   CSPRNG_STATE_T * const csprng_state){
    for(int i = 0; i<N-M; i++){
        for(int j = 0; j<M; j++){
            res[i][j] = zz_rnd_state(csprng_state);
        }
    }
}
#endif

