/* Deterministic randombytes by Daniel J. Bernstein */
/* taken from SUPERCOP (https://bench.cr.yp.to)     */

#include "api.h"
#include "randombytes.h"

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define MAXMLEN 2048

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x####y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(MUPQ_NAMESPACE, fun)

// use different names so we can have empty namespaces
#define MUPQ_CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define MUPQ_CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define MUPQ_CRYPTO_BYTES          NAMESPACE(CRYPTO_BYTES)
#define MUPQ_CRYPTO_ALGNAME        NAMESPACE(CRYPTO_ALGNAME)

#define MUPQ_crypto_sign_keypair NAMESPACE(crypto_sign_keypair)
#define MUPQ_crypto_sign NAMESPACE(crypto_sign)
#define MUPQ_crypto_sign_open NAMESPACE(crypto_sign_open)
#define MUPQ_crypto_sign_signature NAMESPACE(crypto_sign_signature)
#define MUPQ_crypto_sign_verify NAMESPACE(crypto_sign_verify)

typedef uint32_t uint32;

static void printbytes(const unsigned char *x, unsigned long long xlen)
{
  char outs[2*xlen+1];
  unsigned long long i;
  for(i=0;i<xlen;i++) {
    sprintf(outs+2*i, "%02x", x[i]);
  }
  outs[2*xlen] = 0;
  puts(outs);
}

static uint32 seed[32] = { 3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5 } ;
static uint32 in[12];
static uint32 out[8];
static int outleft = 0;

#define ROTATE(x,b) (((x) << (b)) | ((x) >> (32 - (b))))
#define MUSH(i,b) x = t[i] += (((x ^ seed[i]) + sum) ^ ROTATE(x,b));

static void surf(void)
{
  uint32 t[12]; uint32 x; uint32 sum = 0;
  int r; int i; int loop;

  for (i = 0;i < 12;++i) t[i] = in[i] ^ seed[12 + i];
  for (i = 0;i < 8;++i) out[i] = seed[24 + i];
  x = t[11];
  for (loop = 0;loop < 2;++loop) {
    for (r = 0;r < 16;++r) {
      sum += 0x9e3779b9;
      MUSH(0,5) MUSH(1,7) MUSH(2,9) MUSH(3,13)
      MUSH(4,5) MUSH(5,7) MUSH(6,9) MUSH(7,13)
      MUSH(8,5) MUSH(9,7) MUSH(10,9) MUSH(11,13)
    }
    for (i = 0;i < 8;++i) out[i] ^= t[i + 4];
  }
}

int randombytes(uint8_t *x, size_t xlen)
{
  while (xlen > 0) {
    if (!outleft) {
      if (!++in[0]) if (!++in[1]) if (!++in[2]) ++in[3];
      surf();
      outleft = 8;
    }
    *x = out[--outleft];
    ++x;
    --xlen;
  }

  return 0;
}

int main(void)
{
  unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
  unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];

  unsigned char mi[MAXMLEN];
  unsigned char sm[MAXMLEN+MUPQ_CRYPTO_BYTES];
  size_t smlen;
  size_t mlen;

  int r;
  unsigned long long i,k;

  for(i=0; i<MAXMLEN; i=(i==0)?i+1:i<<1)
  {
    randombytes(mi,i);

    MUPQ_crypto_sign_keypair(pk, sk);

    printbytes(pk,MUPQ_CRYPTO_PUBLICKEYBYTES);
    printbytes(sk,MUPQ_CRYPTO_SECRETKEYBYTES);

    MUPQ_crypto_sign(sm, &smlen, mi, i, sk);

    printbytes(sm, smlen);

    // By relying on m == sm we prevent having to allocate CRYPTO_BYTES twice
    r = MUPQ_crypto_sign_open(sm, &mlen, sm, smlen, pk);

    if(r)
    {
      puts("ERROR: signature verification failed");
      return -1;
    }
    for(k=0;k<i;k++)
    {
      if(sm[k]!=mi[k])
      {
        puts("ERROR: message recovery failed");
        return -1;
      }
    }
  }
  return 0;
}
