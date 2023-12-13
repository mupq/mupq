
/**
 * @file crypto_memset.c
 * @brief cryptographic utility functions
 */

#include "crypto_memset.h"

void *(*volatile memset_volatile)(void *, int, size_t) = memset;

int cmp_const(const void *a, const void *b, const size_t size) {
    const unsigned char *_a = (const unsigned char *)a;
    const unsigned char *_b = (const unsigned char *)b;
    unsigned char result = 0;
    size_t i;

    for (i = 0; i < size; i++) {
        result |= _a[i] ^ _b[i];
    }

    return result; /* returns 0 if equal, nonzero otherwise */
}
