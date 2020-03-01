#ifndef GF2M_H
#define GF2M_H

#include "typedef.h"

u2 Add_GF2m(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1Len);
u2 Modular_Reduction_GF2m(pu4 pu4X, pu4 pu4R);
u2 Mupliply_GF2m_PreCompute(pu4 pu4X, pu4 pu4Storage, u1 u1Len);
u2 Mupliply_GF2m_With_PreCompute(pu4 pu4X, pu4 pu4PreCompute, pu4 pu4R, u1 u1Len, u1 u1Degree);
u2 Modular_inverse_GF2m(pu4 pu4X, pu4 pu4Modulo, pu4 pu4R, pu4 pu4Workspace, u1 u1Len);

#endif
