
/**
 * @file crypto_memset.h
 * @brief Header file for cryptographic utility functions
 */

#ifndef CRYPTO_MEMSET_H
#define CRYPTO_MEMSET_H

#include <string.h>

/**
 * safer call to memset https://github.com/veorq/cryptocoding#problem-4
 */
extern void *(*volatile memset_volatile)(void *, int, size_t);

/**
 * constant time call to memcmp
 */
int cmp_const(const void *a, const void *b, const size_t size);

#define memset_zero(ptr, len) memset_volatile((ptr), 0, (len))

#endif  // CRYPTO_MEMSET_H
