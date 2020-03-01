#ifndef GF2MN_H
#define GF2MN_H

#include "typedef.h"

u1 Degree_GF2mn(pu4 pu4Polynom, u1 u1CoeffLength, u1 u1DegreeMaxPolynom);
u1 Add_GF2mn(u1 u1NbCoeff, u1 u1CoeffLength, pu4 pu4X, pu4 pu4Y, pu4 pu4R);
u1 Mupliply_GF2mn(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY, u1 u1ModuloDegree, u1 u1CoeffLength);
u1 Euclidean_Division(pu4 pu4Dividend, pu4 pu4Divisor, pu4 pu4Workspace, u1 u1Divisor_degree, u1 u1Dividend_degree, u1 u1NbCoeff, u1 u1CoeffLength, u1 u1ModuloDegree);
pu4 Modular_inverse_GF2mn(pu4 pu4X, pu4 pu4Workspace, u1 u1NbCoeff, u1 u1CoeffLength, u1 u1ModuloDegree);
u2 SchoolMultiply(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY);
u2 Karatsuba(pu4 pu4X, pu4 pu4Y, pu4 pu4R, pu4 pu4Workspace, u1 u1Degree, u1 u1CoeffLength);
u1 Karatsuba_Degree(u1 u1Degree);
u1 S_i_calculus(pu4 pu4F, pu4 pu4S, pu4 pu4R, u1 u1Index, u1 u1SDim);
u1 Intersection(pu4 pu4V1, pu4 pu4V2, pu4 pu4R, u1 pu4V1Dim, u1 pu4V2Dim, u1 u1CoeffLength);
u2 PreCalculusRSR(pu4 pu4F, pu4 pu4S, pu4 pu4Workspace, pu1 u1TabDimInterSi);
u2 Rank_Support_Recover(pu4 pu4E, pu4 pu4F, pu4 pu4Syndrom, pu4 pu4Workspace);
u1 leading_coefficient_div(pu4 pu4R, pu4 pu4V, pu4 pu4U, u1 degree_v, u1 degree_u, pu4 workspace);
u1 mult_by_element(pu4 pu4R, pu4 pu4V, pu4 pu4scalar, u1 nb_coeff, pu4 workspace, u1 nb_shift);
u1 modular_inverse_GF2mn(pu4 pu4R, pu4 pu4X, pu4 pu4modulo);
u1 init_modulo(pu4 pu4v, int n);
u2 ModularReduction(pu4 pu4X, u1 u1NbCoeffX, u1 u1ModuloDegree, u1 u1CoeffLength);

#endif
