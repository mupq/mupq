#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include "F7Field.h"
#include "F47Field.h"
#include "F61Field.h"
#include "F79Field.h"

#define PRINTMATRIX(M) printf(#M " = \n"); printMatrix(M);

/*Matrix over F_Q*/
typedef struct {
	int rows;
	int cols;
	FELT** array;
} Matrix;

Matrix zeroMatrix(unsigned int rows, unsigned int cols);
void destroy_matrix(Matrix mat);
//void printMatrix(Matrix Mat);
int getUniqueSolution(Matrix A, FELT *solution);

#endif
