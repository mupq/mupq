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

void ColumnGenerator_init(ColumnGenerator * col_gen, const unsigned char* key){
	int i;
	unsigned char stream[16] = {0};

	for(i=0; i<STATES; i++){
		stream[0] = i;
		PRNG_INIT(&col_gen->states[i], key, stream);
	}
	col_gen->cols_used = BLOCK_SIZE/2;
}


column Next_Column(ColumnGenerator *col_gen){
	unsigned char i;
	if(col_gen->cols_used == BLOCK_SIZE/2){
		for(i=0; i<STATES; i++){
			PRNG_GET_BLOCK(&col_gen->states[i],col_gen->blocks[i]);
		}
		col_gen->cols_used = 0;
	}

	#if (OIL_VARS <= 64)
		column Out = 0;
		for(i=0; i<STATES; i++){
			Out |= ((column)col_gen->blocks[i][col_gen->cols_used*2    ]) << (16*i);
			Out |= ((column)col_gen->blocks[i][col_gen->cols_used*2 + 1]) << (16*i+8);
		}
	#else
		column Out = {0};
		for(i=0; i<4; i++){
			Out.components[0] |= ((uint64_t)col_gen->blocks[i][col_gen->cols_used*2    ]) << (16*i);
			Out.components[0] |= ((uint64_t)col_gen->blocks[i][col_gen->cols_used*2 + 1]) << (16*i+8);
		}

		for(; i<STATES; i++){
			Out.components[1] |= ((uint64_t)col_gen->blocks[i][col_gen->cols_used*2 + 0]) << (16*(i-4));
			Out.components[1] |= ((uint64_t)col_gen->blocks[i][col_gen->cols_used*2 + 1]) << (16*(i-4)+8);
		}
	#endif

	col_gen->cols_used ++;
	return Out;
}

void PieceGenerator_init(PieceGenerator * col_gen, const unsigned char* key, int IV){
	unsigned char stream[16] = {0};
	stream[0] = IV;
	PRNG_INIT(&col_gen->state, key, stream);
	col_gen->cols_used = BLOCK_SIZE/2;
}

uint16_t Next_Piece(PieceGenerator *col_gen){
	if(col_gen->cols_used == BLOCK_SIZE/2){
		PRNG_GET_BLOCK(&col_gen->state,col_gen->block);
		col_gen->cols_used = 0;
	}

	uint16_t Out = (((uint16_t)col_gen->block[col_gen->cols_used*2    ]) ) | (((uint16_t) col_gen->block[col_gen->cols_used*2 + 1]) << 8);
	col_gen->cols_used ++;
	return Out;
}
