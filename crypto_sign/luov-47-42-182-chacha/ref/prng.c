#include "prng.h"

/*
	Initializes a Sponge object, absorbs a seed and finalizes the absorbing phase

	 sponge  : The sponge object
	 seed    : The seed to absorb
	 len     : The length of the seed
*/
void initializeAndAbsorb(Sponge *sponge ,const unsigned char * seed , int len ) {
    shake_inc_init(sponge);
    shake_inc_absorb(sponge, seed, len);
    shake_inc_finalize(sponge);
}

/*
	Squeezes a uint64_t from the sponge object

	sponge : The sponge object
	bytes  : The number of bytes to squeeze from the sponge (should be between 1 and 8)
*/
uint64_t squeezeuint64_t(Sponge *sponge, int bytes){
	unsigned char randombits[8];
    shake_inc_squeeze(randombits, bytes, sponge);
	return  ((uint64_t) randombits[0])         | (((uint64_t) randombits[1]) << 8  ) |
	       (((uint64_t) randombits[2]) << 16 ) | (((uint64_t) randombits[3]) << 24 ) |
	       (((uint64_t) randombits[4]) << 32 ) | (((uint64_t) randombits[5]) << 40 ) |
	       (((uint64_t) randombits[6]) << 48 ) | (((uint64_t) randombits[7]) << 56 );
}
