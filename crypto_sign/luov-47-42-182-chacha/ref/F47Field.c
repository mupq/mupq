#include "F47Field.h"

/*
	Write a field element to a char array

	W : writer object
	a : field element to write
*/
void f47serialize_FELT(writer *W, f47FELT a) {
	serialize_uint64_t(W, a , 47);
}

/*
	Read a field element from a char array

	R : reader object

	returns : a field element
*/
f47FELT f47deserialize_FELT(reader *R) {
	return ((f47FELT) deserialize_uint64_t(R,47));
}

static void f47reduce(f47FELT *a){
	f47FELT overflow = *a >> 47;
	*a ^= overflow ^ (overflow << 5);
	*a &= 0x7fffffffffff;
}

/*
	Multiplies two field elements

	a,b : reduced field element to multiply

	return : the reduced product of a and b
*/
f47FELT f47multiply(f47FELT a, f47FELT b) {
	a = a & 0x7fffffffffff;
	b = b & 0x7fffffffffff;
	f47FELT table[16] = {0};
	table[1] = a;
	table[2] = a<<1;
	table[4] = a<<2;
	table[8] = a<<3;

	table[3]  = table[1] ^ table[2];
	table[5]  = table[4] ^ table[1];
	table[6]  = table[4] ^ table[2];
	table[9]  = table[8] ^ table[1];
	table[10] = table[8] ^ table[2];
	table[12] = table[8] ^ table[4];

	table[7]  = table[6] ^ table[1];
	table[11] = table[3] ^ table[8];
	table[13] = table[9] ^ table[4];
	table[14] = table[12] ^ table[2];
	table[15] = table[5] ^ table[10];

	f47FELT out = 0;
	f47FELT out_lower = 0;
	f47FELT out_upper = 0;
	f47FELT b_lower =        b & 0x0f0f0f0f0f0f;
	f47FELT b_upper = (b >> 4) & 0x0f0f0f0f0f0f;

	unsigned char *b_lower_nibbles = (unsigned char *) &b_lower;
	unsigned char *b_upper_nibbles = (unsigned char *) &b_upper;

	out_upper = table[b_upper_nibbles[5]];
	out_lower = table[b_lower_nibbles[5]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[4]];
	out_lower ^= table[b_lower_nibbles[4]];

	out = (out_upper << 4) ^ out_lower;
	f47reduce(&out);

	out_upper = table[b_upper_nibbles[3]];
	out_lower = table[b_lower_nibbles[3]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[2]];
	out_lower ^= table[b_lower_nibbles[2]];

	out = (out <<16) ^ (out_upper << 4) ^ out_lower;
	f47reduce(&out);

	out_upper = table[b_upper_nibbles[1]];
	out_lower = table[b_lower_nibbles[1]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[0]];
	out_lower ^= table[b_lower_nibbles[0]];

	out = (out <<16) ^ (out_upper << 4) ^ out_lower;
	f47reduce(&out);

	return out;
}

/*
	Inverts a field element

	a : field element to invert

	return : the inverse of a, if a is nonzero
*/
f47FELT f47inverse(f47FELT a) {
	uint64_t pow = 0x7ffffffffffe;
	//uint64_t pow = 0x1;
	uint64_t e = a;
	uint64_t ans = 1;
	while(pow != 0){
		if((pow & 1) == 1){
			ans = f47multiply(ans,e);
		}
		e = f47multiply(e,e);
		pow >>= 1;
	}
	return ans;
}

