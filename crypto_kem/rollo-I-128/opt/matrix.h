#ifndef MATRIX_
#define MATRIX_

#include "typedef.h"

u1 SwapLines(pu4 pu4A, pu4 pu4B, u1 u1Length);
u1 Find_pivot(pu4 pu4Matrix, u1 nu1NbRows, u1 nu1Line, u1 u1Pivot, u1 u1ElementLength);
u1 Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace);
u1 Reduced_Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace);

#endif
