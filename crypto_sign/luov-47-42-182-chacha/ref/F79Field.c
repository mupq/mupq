#include "F79Field.h"

f79FELT f79add(f79FELT a, f79FELT b) {
	f79FELT out;
	out.coef[0] = a.coef[0] ^ b.coef[0];
	out.coef[1] = a.coef[1] ^ b.coef[1];
	return out;
}


void f79addInPlace(f79FELT *a, f79FELT *b){
	(*a).coef[0] ^= (*b).coef[0];
	(*a).coef[1] ^= (*b).coef[1];
}

int f79isEqual(f79FELT a, f79FELT b){
	return ((a.coef[0] == b.coef[0]) && ( (a.coef[1] & 0x7fff) == (b.coef[1] & 0x7fff) ) );
}

/*
	Write a field element to a char array

	W : writer object
	a : field element to write
*/
void f79serialize_FELT(writer *W, f79FELT a) {
	serialize_uint64_t(W, a.coef[0] , 64);
	serialize_uint64_t(W, a.coef[1] , 15);
}

/*
	Read a field element from a char array

	R : reader object

	returns : a field element
*/
f79FELT f79deserialize_FELT(reader *R) {
	f79FELT out;
	out.coef[0] = deserialize_uint64_t(R,64);
	out.coef[1] = deserialize_uint64_t(R,15);
	return out;
}

/*
	Multiplies two field elements

	a,b : reduced field element to multiply

	return : the reduced product of a and b
*/
f79FELT f79multiply(f79FELT a, f79FELT b) {
	a.coef[1] &= 0x7fff;
	b.coef[1] &= 0x7fff;
	uint64_t a_upper = a.coef[1];
	uint64_t b_upper = b.coef[1];
	uint64_t a_middle = a.coef[0]>>32;
	uint64_t b_middle = b.coef[0]>>32;
	uint64_t a_lower = a.coef[0] & 0xffffffff;
	uint64_t b_lower = b.coef[0] & 0xffffffff;

	// karatsuba like thing (toom-3 might be better : 5 multiplications insead of 6)
	uint64_t out128 = clmul(a_upper,b_upper);
	uint64_t out0   = clmul(a_lower,b_lower);
	uint64_t out64  = clmul(a_middle,b_middle);
	uint64_t out32  = clmul(a_middle ^ a_lower , b_middle ^ b_lower) ^ out64  ^ out0;
	uint64_t out96  = clmul(a_middle ^ a_upper , b_middle ^ b_upper) ^ out128 ^ out64;
	         out64 ^= clmul(a_lower ^ a_upper , b_lower ^ b_upper) ^ out0 ^ out128;

	// reduce
 	out64  ^= (out32 >> 32);
 	out0   ^= (out32 << 32);
 	out128 ^= (out96 >> 32);
 	out64  ^= (out96 << 32);

 	// x^128 = x^58 + x^49
 	out0  ^= (out128 << 49) ^ (out128 << 58);
 	out64 ^= (out128 >> (64-58)) ^ (out128 >> (64-49));

 	// x^79  =x^9 + 1
 	uint64_t overflow = (out64 >> 15);
 	out64 &= 0x7fff;
 	out0 ^= overflow ^ (overflow << 9);

 	f79FELT out = {{out0, out64}};
 	return out;
}

/*
	Inverts a field element

	a : field element to invert

	return : the inverse of a, if a is nonzero
*/
f79FELT f79inverse(f79FELT a) {
	f79FELT e = f79multiply(a,a);
	f79FELT ans = f79ONE;
	int i;
	for (i = 0; i < 78; ++i)
	{
		ans = f79multiply(ans,e);
		e = f79multiply(e,e);
	}
	return ans;
}

