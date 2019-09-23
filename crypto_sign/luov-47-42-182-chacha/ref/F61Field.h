/*
	Finite field of order 2^61 implemented as polynomial ring F_2[x]/x^61 + x^5 + x^2 + x + 1
*/

#ifndef F61FIELD_H
#define F61FIELD_H

#include <stdint.h>
#include <stdio.h>
#include "buffer.h"
#include "parameters.h"

typedef uint64_t f61FELT;

/* Field operations */

f61FELT f61multiply(f61FELT a, f61FELT b);
f61FELT f61inverse(f61FELT a);
uint8_t f61log(f61FELT);
f61FELT f61antilog(uint8_t);

uint64_t clmul(uint64_t a, uint64_t b);

/* serialization/deserialization */

void f61serialize_FELT(writer *W, f61FELT a);
f61FELT f61deserialize_FELT(reader *R);

#define f61ZERO ((uint64_t) 0 )
#define f61ONE ((uint64_t) 1)
#define f61add(A,B) (A^B)
#define f61addInPlace(A,B) (*A ^= *B)
#define f61isEqual(A,B) ((A & 0x1fffffffffffffff) == (B & 0x1fffffffffffffff) )

#endif
