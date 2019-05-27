#include "Column.h"

#if OIL_VARS>64 

/*
	Write a column to a char array
*/
void serialize_column(writer * W, column b) {
	int a = 0;
	int bits = OIL_VARS;
	while (bits >64 ){
		serialize_uint64_t(W, b.components[a++]  , 64);
		bits -= 64;
	}
	serialize_uint64_t(W, b.components[a], bits);
}

/*
	Read a column from a char array
*/
column deserialize_column(reader * R) {
	column out;
	int a = 0;
	int bits = OIL_VARS;
	while (bits >64 ){
		out.components[a++] = deserialize_uint64_t(R, 64);
		bits -= 64;
	}
	out.components[a] = deserialize_uint64_t(R, bits);
	return out;
}

/*
	xor two columns
*/
column xor(column a, column b) {
	int i;
	column BC;
	for(i=0 ; i<COLUMN_COMPONENTS ; i++){
		BC.components[i] = a.components[i] ^ b.components[i];
	}
	return BC;
}

/*
	Randomize column with Keccak Sponge
*/
column random_column(Sponge *sponge) {
	column BC;
	int i;
	for(i=0 ; i<COLUMN_COMPONENTS-1 ; i++){
		BC.components[i] = squeezeuint64_t(sponge,8);
	}
	BC.components[COLUMN_COMPONENTS-1] = squeezeuint64_t(sponge,((OIL_VARS%64)+7)/8); 
	return BC;
}

/*
	Get a bit from the column
*/
uint64_t getBit(column container, uint64_t bit) {
	return (container.components[bit/64] & ((uint64_t)1) << (bit%64) );
}

/*
	Flip a bit from the column
*/
void flipBit(column *container, uint64_t bit) {
	container->components[bit/64] ^= ((uint64_t)1) << (bit%64);
}

#endif

/*
	Generates an array of columns

	sponge : pointer to a Sponge object
	arr    : the array that will receive the generated columns
	size   : the number of columns that is generated
*/
void squeeze_column_array(Sponge *sponge, column *arr, int size) {
	int i;
	for (i = 0; i < size; i++) {
		arr[i] = random_column(sponge);
	}
}