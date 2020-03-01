#ifndef UTILS_H
#define UTILS_H

#include "typedef.h"

u2 Shift_32bits_Words_Left(pu4 pu4X, pu4 pu4R, u1 u1Len, u1 u1Shift_Number);
u2 Shift_32bits_5Words_Left_4bits(pu4 pu4X);
u2 Is_Unit(pu4 pu4X, u1 u1Len);
u2 Is_Zero(pu4 pu4X, u1 u1Len);
u2 Degree_GF2m(pu4 pu4X, u1 u1Len);
void ConstToRamMemCpy(pu4 pu4Dest,pu4 pu4Src,u2 u2Length);
u1 word_to_byte(pu1 v_destination, pu4 v_initial, u1 len_v, u1 nb_coeff);
u1 byte_to_word(pu4 v_destination, pu1 v_initial, u1 nb_coeff);

#endif
