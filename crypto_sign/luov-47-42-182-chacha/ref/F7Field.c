#include "F7Field.h"

/*
	Write a field element to a char array

	W : writer object
	a : field element to write
*/
void f7serialize_FELT(writer *W, f7FELT a) {
	serialize_uint64_t(W, a , 7);
}

/*
	Read a field element from a char array

	R : reader object

	returns : a field element
*/
f7FELT f7deserialize_FELT(reader *R) {
	return ((f7FELT) deserialize_uint64_t(R,7));
}

/*
	Multiplies two field elements

	a,b : field element to multiply

	return : the product of a and b
*/
f7FELT f7multiply(f7FELT a, f7FELT b) {
	if (((a & 127) == 0) || ((b & 127) == 0))
		return 0;
	return f7antilog((f7log(a) + f7log(b)) % f7units);
}

/*
	Inverts a field element

	a : field element to invert

	return : the inverse of a, if a is nonzero
*/
f7FELT f7inverse(f7FELT a) {
	return f7antilog((f7units - f7log(a)) % f7units);
}

/*
	Calculates the logarithm of a field element with respect to the generator 0x03

	a : a field element

	return : x such that 0x03^x = a or 0 in the case x=0x00
*/
uint8_t f7log(f7FELT a)
{
	static const uint8_t f7LogTable[256] =
	{
		  0,   0,   1,   7,   2,  14,   8,  56,   3,  63,  15,  31,   9,  90,  57,  21,
		  4,  28,  64,  67,  16, 112,  32,  97,  10, 108,  91,  70,  58,  38,  22,  47,
		  5,  54,  29,  19,  65,  95,  68,  45,  17,  43, 113, 115,  33,  77,  98, 117,
		 11,  87, 109,  35,  92,  74,  71,  79,  59, 104,  39, 100,  23,  82,  48, 119,
		  6, 126,  55,  13,  30,  62,  20,  89,  66,  27,  96, 111,  69, 107,  46,  37,
		 18,  53,  44,  94, 114,  42, 116,  76,  34,  86,  78,  73,  99, 103, 118,  81,
		 12, 125,  88,  61, 110,  26,  36, 106,  93,  52,  75,  41,  72,  85,  80, 102,
		 60, 124, 105,  25,  40,  51, 101,  84,  24, 123,  83,  50,  49, 122, 120, 121,
		  0,   0,   1,   7,   2,  14,   8,  56,   3,  63,  15,  31,   9,  90,  57,  21,
		  4,  28,  64,  67,  16, 112,  32,  97,  10, 108,  91,  70,  58,  38,  22,  47,
		  5,  54,  29,  19,  65,  95,  68,  45,  17,  43, 113, 115,  33,  77,  98, 117,
		 11,  87, 109,  35,  92,  74,  71,  79,  59, 104,  39, 100,  23,  82,  48, 119,
		  6, 126,  55,  13,  30,  62,  20,  89,  66,  27,  96, 111,  69, 107,  46,  37,
		 18,  53,  44,  94, 114,  42, 116,  76,  34,  86,  78,  73,  99, 103, 118,  81,
		 12, 125,  88,  61, 110,  26,  36, 106,  93,  52,  75,  41,  72,  85,  80, 102,
		 60, 124, 105,  25,  40,  51, 101,  84,  24, 123,  83,  50,  49, 122, 120, 121
	};
	return f7LogTable[a];
}

/*
	Calculates the inverse of the logarithm with respect to the generator 0x03

	a : an integer

	return : the field element 0x03^a
*/
uint8_t f7antilog(uint8_t a)
{
	static const uint8_t f7AntiLogTable[256] =
	{
		 1,   2,   4,   8,  16,  32,  64,   3,   6,  12,  24,  48,  96,  67,   5,  10,
		 20,  40,  80,  35,  70,  15,  30,  60, 120, 115, 101,  73,  17,  34,  68,  11,
		 22,  44,  88,  51, 102,  79,  29,  58, 116, 107,  85,  41,  82,  39,  78,  31,
		 62, 124, 123, 117, 105,  81,  33,  66,   7,  14,  28,  56, 112,  99,  69,   9,
		 18,  36,  72,  19,  38,  76,  27,  54, 108,  91,  53, 106,  87,  45,  90,  55,
		110,  95,  61, 122, 119, 109,  89,  49,  98,  71,  13,  26,  52, 104,  83,  37,
		 74,  23,  46,  92,  59, 118, 111,  93,  57, 114, 103,  77,  25,  50, 100,  75,
		 21,  42,  84,  43,  86,  47,  94,  63, 126, 127, 125, 121, 113,  97,  65,
		  1,   2,   4,   8,  16,  32,  64,   3,   6,  12,  24,  48,  96,  67,   5,  10,
		 20,  40,  80,  35,  70,  15,  30,  60, 120, 115, 101,  73,  17,  34,  68,  11,
		 22,  44,  88,  51, 102,  79,  29,  58, 116, 107,  85,  41,  82,  39,  78,  31,
		 62, 124, 123, 117, 105,  81,  33,  66,   7,  14,  28,  56, 112,  99,  69,   9,
		 18,  36,  72,  19,  38,  76,  27,  54, 108,  91,  53, 106,  87,  45,  90,  55,
		110,  95,  61, 122, 119, 109,  89,  49,  98,  71,  13,  26,  52, 104,  83,  37,
		 74,  23,  46,  92,  59, 118, 111,  93,  57, 114, 103,  77,  25,  50, 100,  75,
		 21,  42,  84,  43,  86,  47,  94,  63, 126, 127, 125, 121, 113,  97,  65,   1,
 };
	return f7AntiLogTable[a];
}
