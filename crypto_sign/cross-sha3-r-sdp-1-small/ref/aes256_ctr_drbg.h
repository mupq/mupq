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

#define CSPRNG_BUFFER_LENGTH (3*SEED_LENGTH_BYTES)

typedef struct {
   unsigned char   Key[32];
   unsigned char   Vee[16];
   int             reseed_counter;
   unsigned char   buffer[CSPRNG_BUFFER_LENGTH];
   int             first_cached_byte;
} AES256_CTR_DRBG_struct;


/* initializes a CSPRNG, given the seed and a state pointer */
void initialize_aes_csprng(AES256_CTR_DRBG_struct * const csprng_state,
                           const unsigned char * const seed,
                           const uint32_t seed_len_bytes);

/* extracts xlen bytes from the CSPRNG, given the state */
void aes_csprng_randombytes(unsigned char * const output,
                            unsigned long long output_len,
                            AES256_CTR_DRBG_struct * const csprng_state);
