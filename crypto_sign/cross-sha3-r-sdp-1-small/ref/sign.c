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

#include "api.h"
#include <stddef.h>  //  size_t
#include <string.h>  // memcpy
/*----------------------------------------------------------------------------*/

int crypto_sign_keypair(unsigned char *pk,
                        unsigned char *sk)
{
   /* keygen cannot fail */
   CROSS_keygen((prikey_t *) sk,
               (pubkey_t *) pk);

   return 0; // NIST convention: 0 == zero errors
} // end crypto_sign_keypair

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*... generating a signed message sm[0],sm[1],...,sm[*smlen-1]                */
/*... from original message m[0],m[1],...,m[mlen-1]                           */
/*... under secret key sk[0],sk[1],...                                        */
int crypto_sign(unsigned char *sm,
                size_t *smlen,     // out parameter
                const unsigned char *m, size_t mlen,  // in parameter
                const unsigned char *sk)                          // in parameter
{
   // TODO: this should not be needed, but CROSS_sign seems to not set all bytes of the signature
   memset(sm+mlen, 0, sizeof(sig_t));
   /* sign cannot fail */
   memcpy((unsigned char *) sm, (const unsigned char *) m, (size_t) mlen);
   CROSS_sign((const prikey_t *)
             sk,                                // in parameter
             (const char *const) m, (const uint64_t) mlen,         // in parameter
             (sig_t *) (sm+mlen));                                 // out parameter
   *smlen = mlen + (size_t) sizeof(sig_t);

   return 0;  // NIST convention: 0 == zero errors
} // end crypto_sign

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*.  ... verifying a signed message sm[0],sm[1],...,sm[smlen-1]               */
/*.  ... under public key pk[0],pk[1],...                                     */
/*.  ... and producing original message m[0],m[1],...,m[*mlen-1]              */
int crypto_sign_open(unsigned char *m,
                     size_t *mlen,        // out parameter
                     const unsigned char *sm, size_t smlen, // in parameter
                     const unsigned char *pk)                           // in parameter
{

   /* verify returns 1 if signature is ok, 0 otherwise */
   *mlen = smlen-(size_t) sizeof(sig_t);
   memcpy((unsigned char *) m, (const unsigned char *) sm, (size_t) *mlen);
   int ok = CROSS_verify((const pubkey_t *const)
                        pk,                     // in parameter
                        (const char *const) m, (const uint64_t) *mlen,  // in parameter
                        (const sig_t * const) (sm+*mlen));              // in parameter


   return ok-1; // NIST convention: 0 == zero errors, -1 == error condition
} // end crypto_sign_open

/*----------------------------------------------------------------------------*/
