/*
	Finite field of order 2^47 implemented as polynomial ring F_2[x]/x^79 + x^9 + 1
*/

#ifndef F79FIELD_H
#define F79FIELD_H

#include <stdint.h>
#include <stdio.h>
#include "buffer.h"
#include "parameters.h"
#include "F61Field.h"

typedef struct {
	uint64_t coef[2];
} f79FELT;

static const f79FELT f79ONE = {{ 1,0 }};
static const f79FELT f79ZERO = {{ 0,0}};

/* Field operations */

f79FELT f79add(f79FELT a, f79FELT b);
void f79addInPlace(f79FELT *a, f79FELT *b);
int f79isEqual(f79FELT a, f79FELT b);
f79FELT f79multiply(f79FELT a, f79FELT b);
f79FELT f79inverse(f79FELT a);
uint8_t f79log(f79FELT);
f79FELT f79antilog(uint8_t);

/* serialization/deserialization */

void f79serialize_FELT(writer *W, f79FELT a);
f79FELT f79deserialize_FELT(reader *R);

#endif
