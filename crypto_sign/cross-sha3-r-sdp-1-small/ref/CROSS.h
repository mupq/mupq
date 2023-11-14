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
#include "pack_unpack.h"
#include <stdint.h>

/* Public key: the parity check matrix is shrunk to a seed, syndrome
 * represented in full */
typedef struct {
   uint8_t seed_pub[SEED_LENGTH_BYTES];
   uint8_t s[DENSELY_PACKED_FQ_SYN_SIZE];
} pubkey_t;

/* Private key: just a single seed*/
typedef struct {
   uint8_t seed[SEED_LENGTH_BYTES];
} prikey_t;

typedef struct {
  uint8_t y[DENSELY_PACKED_FQ_VEC_SIZE];
  /*FZ_elem should be compressed as they are 3b wide */
#if defined(RSDP)
  uint8_t sigma[DENSELY_PACKED_FZ_VEC_SIZE];
#elif defined(RSDPG)
  uint8_t delta[DENSELY_PACKED_FZ_RSDP_G_VEC_SIZE];
#endif
} rsp_0_t;

/* Signature: */
typedef struct {
   uint8_t salt[SALT_LENGTH_BYTES];
   uint8_t digest_01[HASH_DIGEST_LENGTH];
   uint8_t digest_b[HASH_DIGEST_LENGTH];
   /*Merkle tree proof field. Currently all leaves, to be replaced */
   uint8_t mtp[HASH_DIGEST_LENGTH*TREE_NODES_TO_STORE];
   /*seed tree paths storage*/
   uint8_t stp[SEED_TREE_MAX_PUBLISHED_BYTES];
   rsp_0_t rsp_0[T-W];
   uint8_t rsp_1[T-W][HASH_DIGEST_LENGTH];
} sig_t;


/* keygen cannot fail */
void CROSS_keygen(prikey_t *SK,
                 pubkey_t *PK);

/* sign cannot fail */
void CROSS_sign(const prikey_t * const SK,
                const char * const m,
                const uint64_t mlen,
                sig_t * const sig);

/* verify returns 1 if signature is ok, 0 otherwise */
int CROSS_verify(const pubkey_t * const PK,
                 const char * const m,
                 const uint64_t mlen,
                 const sig_t * const sig);
