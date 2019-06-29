#include "api.h"
#include "randombytes.h"
#include "hal.h"

#include <string.h>

#define NTESTS 15
#define MLEN 32

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

const uint8_t canary[8] = {
  0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF
};

/* allocate a bit more for all keys and messages and
 * make sure it is not touched by the implementations.
 */
static void write_canary(uint8_t *d) {
  for (size_t i = 0; i < 8; i++) {
    d[i] = canary[i];
  }
}

static int check_canary(const uint8_t *d) {
  for (size_t i = 0; i < 8; i++) {
    if (d[i] != canary[i]) {
      return -1;
    }
  }
  return 0;
}


static int test_sign(void)
{
    unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES+16];
    unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES+16];
    unsigned char sm[MLEN + MUPQ_CRYPTO_BYTES+16];
    unsigned char m[MLEN+16];

    size_t mlen;
    size_t smlen;

    int i;
    write_canary(pk); write_canary(pk+sizeof(pk)-8);
    write_canary(sk); write_canary(sk+sizeof(sk)-8);
    write_canary(sm); write_canary(sm+sizeof(sm)-8);
    write_canary(m); write_canary(m+sizeof(m)-8);

    for (i = 0; i < NTESTS; i++) {
        MUPQ_crypto_sign_keypair(pk+8, sk+8);
        hal_send_str("crypto_sign_keypair DONE.\n");

        randombytes(m+8, MLEN);
        MUPQ_crypto_sign(sm+8, &smlen, m+8, MLEN, sk+8);
        hal_send_str("crypto_sign DONE.\n");

        // By relying on m == sm we prevent having to allocate CRYPTO_BYTES twice
        if (MUPQ_crypto_sign_open(sm+8, &mlen, sm+8, smlen, pk+8))
        {
            hal_send_str("ERROR Signature did not verify correctly!\n");
        }
        else if(check_canary(pk) || check_canary(pk+sizeof(pk)-8) ||
            check_canary(sk) || check_canary(sk+sizeof(sk)-8) ||
            check_canary(sm) || check_canary(sm+sizeof(sm)-8) ||
            check_canary(m) || check_canary(m+sizeof(m)-8))
        {
            hal_send_str("ERROR canary overwritten\n");
        }
        else
        {
            hal_send_str("OK Signature did verify correctly!\n");
        }
        hal_send_str("crypto_sign_open DONE.\n");
    }

    return 0;
}

static int test_wrong_pk(void)
{
    unsigned char pk[MUPQ_CRYPTO_PUBLICKEYBYTES];
    unsigned char pk2[MUPQ_CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[MUPQ_CRYPTO_SECRETKEYBYTES];
    unsigned char sm[MLEN + MUPQ_CRYPTO_BYTES];
    unsigned char m[MLEN];

    size_t mlen;
    size_t smlen;

    int i;

    for (i = 0; i < NTESTS; i++) {
        MUPQ_crypto_sign_keypair(pk2, sk);
        hal_send_str("crypto_sign_keypair DONE.\n");

        MUPQ_crypto_sign_keypair(pk, sk);
        hal_send_str("crypto_sign_keypair DONE.\n");


        randombytes(m, MLEN);
        MUPQ_crypto_sign(sm, &smlen, m, MLEN, sk);
        hal_send_str("crypto_sign DONE.\n");

        // By relying on m == sm we prevent having to allocate CRYPTO_BYTES twice
        if (MUPQ_crypto_sign_open(sm, &mlen, sm, smlen, pk2))
        {
            hal_send_str("OK Signature did not verify correctly under wrong public key!\n");
        }
        else
        {
            hal_send_str("ERROR Signature did verify correctly under wrong public key!\n");
        }
        hal_send_str("crypto_sign_open DONE.\n");
    }

    return 0;
}

int main(void)
{
    hal_setup(CLOCK_FAST);

    // marker for automated testing
    hal_send_str("==========================");
    test_sign();
    test_wrong_pk();
    hal_send_str("#");
    while(1);
    return 0;
}
