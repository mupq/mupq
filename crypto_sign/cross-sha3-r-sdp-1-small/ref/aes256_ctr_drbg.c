
/**
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

#include <stddef.h>
#include <stdint.h>
#include <string.h> // void *memset(void *s, int c, size_t n);
#include "aes256_ctr_drbg.h"
#include "aes.h"


/*setup-and-encrypt AES ECB wrapper for CTR drbg*/
static inline
void AES256_ECB(unsigned char *key, unsigned char *ptx, unsigned char *ctx) {
   aes256ctx round_key;
   aes256_ecb_keyexp(&round_key, key);
   aes256_ecb(ctx, ptx, 1, &round_key);
   aes256_ctx_release(&round_key);
}

/* aes_256_ctr_drbg_update updates the CTR-DRBG state.
 * If the seed_material, which act as the DRBG seed are a NULL pointer,
 * then the update function simply moves on the DRBG state in a non
 * backtrackable way.
 * The DRBG state is a pair: Vee, Key, where Vee is an AES state, Key is
 * an AES-256 Key.
 * The state is moved on overwriting it with
 *    AES-256-Enc(Key, Vee)|| AES-256-Enc(Key, Vee+1)|| AES-256-Enc(Key, Vee+2)
 * If the seed material is present, it is then xor-added to the updated state
 */
static
void aes_256_ctr_drbg_update(unsigned char *seed_material,
                       unsigned char *Key,
                       unsigned char *Vee) {
   unsigned char temp[48];
   for (int i=0; i<3; i++) {
      //increment V
      for (int j=15; j>=0; j--) {
         if ( Vee[j] == 0xff )
            Vee[j] = 0x00;
         else {
            Vee[j]++;
            break;
         }
      }
      AES256_ECB(Key, Vee, temp+16*i);
   }
   if ( seed_material != NULL ) {
      for (int i=0; i<48; i++) {
         temp[i] ^= seed_material[i];
      }
   }
   memcpy(Key, temp, 32);
   memcpy(Vee, temp+32, 16);
}


void initialize_aes_csprng(AES256_CTR_DRBG_struct * const csprng_state,
                           const unsigned char * const seed,
                           const uint32_t seed_len_bytes){
   memset(csprng_state->Key, 0x00, 32);
   memset(csprng_state->Vee, 0x00, 16);

    /* AES CTR DRBG can be initialized with 48B of seed at a time */
   unsigned char seed_material_buf[48];
   int remaining_seed_len = seed_len_bytes;
   while(remaining_seed_len > 0){
       int seed_len_to_add = remaining_seed_len > 48 ? 48 : remaining_seed_len;
       memset(seed_material_buf, 0x00, 48);
       memcpy(seed_material_buf, seed, seed_len_to_add);
       remaining_seed_len -= remaining_seed_len;
       aes_256_ctr_drbg_update(seed_material_buf,
                               csprng_state->Key,
                               csprng_state->Vee);
   }
   csprng_state->reseed_counter = 1;
   memset(csprng_state->buffer,0,CSPRNG_BUFFER_LENGTH);
   csprng_state->first_cached_byte = CSPRNG_BUFFER_LENGTH-1;
}

static
void aes_csprng_refill_buffer(AES256_CTR_DRBG_struct * const csprng_state){
   for(int i = 0; i < CSPRNG_BUFFER_LENGTH/SEED_LENGTH_BYTES;i++ ){
       //increment Vee as a 16B (128b counter)
       for (int j=15; j>=0; j--) {
          if ( csprng_state->Vee[j] == 0xff )
             csprng_state->Vee[j] = 0x00;
          else {
             csprng_state->Vee[j]++;
             break;
          }
       }
       AES256_ECB(csprng_state->Key,
                  csprng_state->Vee,
                  (csprng_state->buffer)+i*SEED_LENGTH_BYTES);
   }
   csprng_state->first_cached_byte = 0;
   /*after having the actual "call to AES-CTR-DRBG, perform state update */
   aes_256_ctr_drbg_update(NULL, csprng_state->Key, csprng_state->Vee);
}


/* Actual DRBG computation (pseudorandom bytes extraction) */
void aes_csprng_randombytes(unsigned char *const output,
                            const unsigned long long output_len,
                            AES256_CTR_DRBG_struct * const csprng_state){
   int length_left_to_emit;
   length_left_to_emit = output_len;

   int idx_to_write = 0;
   while ( length_left_to_emit > 0 ) {
      /* if there is still computed material in the buffer, use that*/
      int leftover_bytes = CSPRNG_BUFFER_LENGTH-csprng_state->first_cached_byte;
      if(leftover_bytes > 0){
         /*if they're enough, just copy them */
         if(leftover_bytes >= length_left_to_emit){
           memcpy(&output[idx_to_write],
                  &csprng_state->buffer[csprng_state->first_cached_byte],
                  length_left_to_emit);
           csprng_state->first_cached_byte += length_left_to_emit;
           return;
         } else { /*otherwise copy everything that's left*/
           memcpy(&output[idx_to_write],
                  &csprng_state->buffer[csprng_state->first_cached_byte],
                  leftover_bytes);
           idx_to_write += leftover_bytes;
           csprng_state->first_cached_byte = CSPRNG_BUFFER_LENGTH-1;
           length_left_to_emit -= leftover_bytes;
         }
      }
      aes_csprng_refill_buffer(csprng_state);
   }
   csprng_state->reseed_counter++;
} // end aes_csprng_randombytes
