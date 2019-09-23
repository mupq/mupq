/*
	Finite field of order 2^7 implemented as polynomial ring F_2[x] mod x^7 +x + 1
*/

#ifndef F7FIELD_H
#define F7FIELD_H

#include <stdint.h>
#include <stdio.h>
#include "buffer.h"
#include "parameters.h"

enum { twoPow7 = 128, f7units = twoPow7 - 1 };
typedef uint8_t f7FELT;

/* Field operations */

f7FELT f7multiply(f7FELT a, f7FELT b);
f7FELT f7inverse(f7FELT a);
uint8_t f7log(f7FELT);
f7FELT f7antilog(uint8_t);

/* serialization/deserialization */

void f7serialize_FELT(writer *W, f7FELT a);
f7FELT f7deserialize_FELT(reader *R);

#define f7ZERO 0
#define f7ONE 1
#define f7add(A,B) (A^B)
#define f7isEqual(A,B) ((A & 127) == (B & 127) )

#endif
