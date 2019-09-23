#include "F61Field.h"

/*
	Write a field element to a char array

	W : writer object
	a : field element to write
*/
void f61serialize_FELT(writer *W, f61FELT a) {
	serialize_uint64_t(W, a , 61);
}

/*
	Read a field element from a char array

	R : reader object

	returns : a field element
*/
f61FELT f61deserialize_FELT(reader *R) {
	return ((f61FELT) deserialize_uint64_t(R,61));
}

/*
	takes two 32 bit values a, b and computes the carryless multiplication of a and b (a 63 bit value)
*/
uint64_t clmul(uint64_t a, uint64_t b){
	uint64_t table[16] = {0};
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

	uint64_t out = 0;
	uint64_t out_lower = 0;
	uint64_t out_upper = 0;
	uint64_t b_lower =        b & 0x0f0f0f0f0f0f;
	uint64_t b_upper = (b >> 4) & 0x0f0f0f0f0f0f;

	unsigned char *b_lower_nibbles = (unsigned char *) &b_lower;
	unsigned char *b_upper_nibbles = (unsigned char *) &b_upper;

	out_upper = table[b_upper_nibbles[3]];
	out_lower = table[b_lower_nibbles[3]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[2]];
	out_lower ^= table[b_lower_nibbles[2]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[1]];
	out_lower ^= table[b_lower_nibbles[1]];
	out_upper <<= 8;
	out_lower <<= 8;
	out_upper ^= table[b_upper_nibbles[0]];
	out_lower ^= table[b_lower_nibbles[0]];

	out = (out_upper << 4) ^ out_lower;
	return out;
}

/*
	Multiplies two field elements

	a,b : reduced field element to multiply

	return : the reduced product of a and b
*/
f61FELT f61multiply(f61FELT a, f61FELT b) {
	a &= 0x1fffffffffffffff;
	b &= 0x1fffffffffffffff;
	uint64_t a_upper = a >> 32;
	uint64_t a_lower = a & 0xffffffff;
	uint64_t b_upper = b >> 32;
	uint64_t b_lower = b & 0xffffffff;

	// karatsuba
	uint64_t out64  = clmul(a_upper,b_upper);
	uint64_t out0   = clmul(a_lower,b_lower);
	uint64_t out32  = clmul(a_upper ^ a_lower , b_upper ^ b_lower) ^ out64 ^ out0;

//	printf(" out0: %lu \n", out0);
//	printf("out32: %lu \n", out32);
//	printf("out64: %lu \n", out64);

	// reduce
 	out64 ^= (out32 >> 32);
 	out0  ^= (out32 << 32);

//	printf("\n");
//	printf(" out0: %lu \n", out0);
//	printf("out64: %lu \n", out64);

 	out0 ^= (out64 << 8) ^ (out64 << 5) ^ (out64 << 4) ^ (out64 << 3);
 	uint64_t overflow = (out0 >> 61) ^ (( (out64 >> 56) ^ (out64 >> 59) ^ (out64 >> 60) ^ (out64 >> 61) ) << 3 );
 	out0 &= 0x1fffffffffffffff;

// 	printf("\n");
// 	printf("out0: %lu \n", out0);
//	printf("overflow %lu \n", overflow);

 	out0 ^= overflow ^ (overflow << 1) ^ (overflow << 2) ^ (overflow << 5);
 	return out0;
}

/*
	Inverts a field element

	a : field element to invert

	return : the inverse of a, if a is nonzero
*/
f61FELT f61inverse(f61FELT a) {
	uint64_t pow = 0x1ffffffffffffffe;
	uint64_t e = a;
	uint64_t ans = 1;
	while(pow != 0){
		if((pow & 1) == 1){
			ans = f61multiply(ans,e);
		}
		e = f61multiply(e,e);
		pow >>= 1;
	}
	return ans;
}

