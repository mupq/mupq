/* Based on the public domain implementation in
 * crypto_hash/keccakc512/simple/ from http://bench.cr.yp.to/supercop.html
 * by Ronny Van Keer
 * and the public domain "TweetFips202" implementation
 * from https://twitter.com/tweetfips202
 * by Gilles Van Assche, Daniel J. Bernstein, and Peter Schwabe */

#include <stddef.h>
#include <stdint.h>

#include "ascon_api.h"

#include "permutations.h"

#define RATE (64 / 8)
#define PA_ROUNDS 12
#define CRYPTO_BYTES 32
#define IV                                            \
  ((u64)(8 * (RATE)) << 48 | (u64)(PA_ROUNDS) << 40 | \
   (u64)(8 * (CRYPTO_BYTES)) << 0)


/*************************************************
 * Name:        ascon_hash
 *
 * Description: ASCON hash function with non-incremental API
 *
 * Arguments:   - uint8_t *output:      pointer to output
                - size_t outlen:        length of output in bytes
 *              - const uint8_t *input: pointer to input
 *              - size_t inlen:         length of input in bytes
 **************************************************/

void ascon_hash(unsigned char *out, size_t outlen, const unsigned char *in, unsigned long long inlen) 
{
  state s;
  size_t outlength;

  // initialization
  s.x0 = IV;
  s.x1 = 0;
  s.x2 = 0;
  s.x3 = 0;
  s.x4 = 0;
  printstate("initial value:", s);
  P12(&s);
  printstate("initialization:", s);

  // absorb plaintext
 
  while (inlen >= RATE) 
  {
    s.x0 ^= BYTES_TO_U64(in, 8);
    P12(&s);
    inlen -= RATE;
    in += RATE;
  }
  s.x0 ^= BYTES_TO_U64(in, inlen);
  s.x0 ^= 0x80ull << (56 - 8 * inlen);
  printstate("absorb plaintext:", s);

  P12(&s);
  printstate("finalization:", s);

  // set hash output

  outlength = outlen;
  while (outlength > RATE) 
  {
    U64_TO_BYTES(out, s.x0, 8);
    P12(&s);
    outlength -= RATE;
    out += RATE;
  }
  U64_TO_BYTES(out, s.x0, 8);
}

