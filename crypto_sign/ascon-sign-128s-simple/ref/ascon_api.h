#ifndef ASCON_H
#define ASCON_H

#include <stddef.h>
#include <stdint.h>

void ascon_hash(unsigned char *out, size_t outlen, const unsigned char *in, unsigned long long inlen);

#endif
