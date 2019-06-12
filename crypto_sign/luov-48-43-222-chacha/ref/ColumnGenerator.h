#include "prng.h"
#include "parameters.h"
#include "Column.h"

#define STATES (((OIL_VARS-1)/16)+1)

typedef struct {
	PRNG_STATE states[STATES];
	unsigned char blocks[STATES][BLOCK_SIZE];
	int cols_used;
} ColumnGenerator;

void ColumnGenerator_init(ColumnGenerator * col_gen, const unsigned char* key);
column Next_Column(ColumnGenerator *col_gen);

typedef struct {
	PRNG_STATE state;
	unsigned char block[BLOCK_SIZE];
	int cols_used;
} PieceGenerator;

void PieceGenerator_init(PieceGenerator * col_gen, const unsigned char* key, int IV);
uint16_t Next_Piece(PieceGenerator *col_gen);