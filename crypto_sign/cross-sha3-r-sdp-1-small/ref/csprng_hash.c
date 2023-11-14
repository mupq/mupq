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

#include "csprng_hash.h"

/******************************************************************************/
/* Initializes a CSPRNG from either an input seed or the output of
 * clock_gettime. Input seed assumed to be a C convention string */

CSPRNG_STATE_T platform_csprng_state;

void sha3_hash(uint8_t digest[HASH_DIGEST_LENGTH],
               const unsigned char *const m,
               const uint64_t mlen)
{
#if (HASH_DIGEST_LENGTH*8 == 256)
   sha3_256(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 384)
   sha3_384(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 512)
   sha3_512(digest,m, mlen);
#else
#error digest length unsupported by SHA-3
#endif
} /* end compute_digest */

void sha2_hash(uint8_t digest[HASH_DIGEST_LENGTH],
               const unsigned char *const m,
               const uint64_t mlen)
{
#if (HASH_DIGEST_LENGTH*8 == 256)
   sha2_256(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 384)
   sha2_384(digest,m, mlen);
#elif (HASH_DIGEST_LENGTH*8 == 512)
   sha2_512(digest,m, mlen);
#else
#error digest length unsupported by SHA-2
#endif
} /* end compute_digest */

#define  POSITION_MASK (( (uint16_t)1 << BITS_TO_REPRESENT(T-1))-1)
void expand_digest_to_fixed_weight(uint8_t fixed_weight_string[T],
                                   const uint8_t digest[HASH_DIGEST_LENGTH]){
   CSPRNG_STATE_T csprng_state;
   initialize_csprng(&csprng_state,
                     (const unsigned char *) digest,
                     HASH_DIGEST_LENGTH);
   uint16_t rnd_buf;
   int placed_elements = 0;
   while (placed_elements < T-W) {
      POSITION_IN_FW_STRING_T pos;
      do {
         csprng_randombytes((unsigned char *) &rnd_buf,
                             sizeof(uint16_t),
                             &csprng_state);
         pos   = rnd_buf & POSITION_MASK;
      } while ( (  pos >= T) ||             /* rejection sampling */
                (fixed_weight_string[pos] != 0) ); /* skip elements already placed */
      fixed_weight_string[pos] = 1;
      placed_elements += 1;
   }
} /* expand_digest_to_fixed_weight */
