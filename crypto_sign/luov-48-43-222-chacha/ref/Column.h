/*
	Nothing fancy here, just an implementation of a container that stores OIL_VARS bits and supports basic functionalities such as xoring, reading bits and flipping a bit
*/

#ifndef COLUMN_H
#define COLUMN_H

#include <stdint.h>
#include "buffer.h"
#include "prng.h"
#include "parameters.h"

#if OIL_VARS > 64

#define COLUMN_COMPONENTS ((OIL_VARS+63)/64)

typedef struct {
	uint64_t components[COLUMN_COMPONENTS];
} column;

static const column empty = {0};

void serialize_column(writer * Buff, column b);
column deserialize_column(reader *Buff);
column xor(column a, column b);
column random_column(Sponge *sponge);
uint64_t getBit(column container, uint64_t bit);
void flipBit(column *container, uint64_t bit);

#else

#define column uint64_t
#define empty ((uint64_t) 0)
#define xor(a,b) a^b
#define getBit(container,bit) (container & ((uint64_t)1) << bit)
#define flipBit(container,bit) (*container ^= ((uint64_t)1) << bit)
#define random_column(sponge) squeezeuint64_t(sponge,((OIL_VARS+7)/8))
#define serialize_column(W,container) serialize_uint64_t(W, container , OIL_VARS)
#define deserialize_column(R) deserialize_uint64_t(R,OIL_VARS)

#endif

void squeeze_column_array(Sponge *sponge, column *arr, int size);

#endif

