/*
	Finite field of order 2^47 implemented as polynomial ring F_2[x]/x^47 + x^5 + 1
*/

#ifndef F47FIELD_H
#define F47FIELD_H

#include <stdint.h>
#include <stdio.h>
#include "buffer.h"
#include "parameters.h"

typedef uint64_t f47FELT;

/* Field operations */

f47FELT f47multiply(f47FELT a, f47FELT b);
f47FELT f47inverse(f47FELT a);
uint8_t f47log(f47FELT);
f47FELT f47antilog(uint8_t);

/* serialization/deserialization */

void f47serialize_FELT(writer *W, f47FELT a);
f47FELT f47deserialize_FELT(reader *R);

#define f47ZERO ((uint64_t) 0 )
#define f47ONE ((uint64_t) 1)
#define f47add(A,B) (A^B)
#define f47addInPlace(A,B) (*A ^= *B)
#define f47isEqual(A,B) ((A & 0x7fffffffffff) == (B & 0x7fffffffffff) )

#endif
