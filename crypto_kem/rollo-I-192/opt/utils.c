#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "gf2m.h"
#include "api.h"
#include "utils.h"

u2 Shift_32bits_Words_Left(pu4 pu4X, pu4 pu4R, u1 u1Len, u1 u1Shift_Number)
{
    u4 i;
    u4 u4CarryLSB = 0;
    u4 u4CarryMSB = 0;

    if (u1Shift_Number==0)
    {
    	for(i=0;i<u1Len;i++)
    	{
    		pu4R[i] = *(pu4X + i);
    	}
    	return 0;
    }

    for(i = 0 ; i < u1Len; i++)
    {
        *(pu4R + i)    = *(pu4X + i);
        u4CarryMSB     = *(pu4R + i) >> (32 - u1Shift_Number) ;
        pu4R[i]      <<= u1Shift_Number;
        pu4R[i]       ^= u4CarryLSB;
        u4CarryLSB     = u4CarryMSB;
    }


    *(pu4R) ^= u4CarryLSB;

    return 0;
}


u2 Shift_32bits_5Words_Left_4bits(pu4 pu4X)
{


    u4 u4CarryMSB[6] = { ((*(pu4X    )) & 0xF0000000) >> 28,
                         ((*(pu4X + 1)) & 0xF0000000) >> 28,
                         ((*(pu4X + 2)) & 0xF0000000) >> 28,
                         ((*(pu4X + 3)) & 0xF0000000) >> 28,
                         ((*(pu4X + 4)) & 0xF0000000) >> 28,
                         ((*(pu4X + 5)) & 0xF0000000) >> 28
                       };

    * pu4X      <<= 4;  * pu4X      ^= u4CarryMSB[5];
    *(pu4X + 1) <<= 4;  *(pu4X + 1) ^= u4CarryMSB[0];
    *(pu4X + 2) <<= 4;  *(pu4X + 2) ^= u4CarryMSB[1];
    *(pu4X + 3) <<= 4;  *(pu4X + 3) ^= u4CarryMSB[2];
    *(pu4X + 4) <<= 4;  *(pu4X + 4) ^= u4CarryMSB[3];
    *(pu4X + 5) <<= 4;  *(pu4X + 5) ^= u4CarryMSB[4];

    return 0;

}


/**********************************************************************
*      Is_Unit checks if the polynomial is equals to 1.
*
* pu4X           : Address of the element to compare to 1.
* u1Len          : Length of pu4X in 32-bits words.
*
**********************************************************************/


u2 Is_Unit(pu4 pu4X, u1 u1Len)
{
    u1 u1i;
    if( pu4X[0] == 0x01)
    {
        for(u1i = 1 ; u1i < u1Len; u1i++)
        {
            if( *(pu4X + u1i) != 0)
                return 1;
        }

        return 0;
    }

    return 1;
}



u2 Is_Zero(pu4 pu4X, u1 u1Len)
{
    u1 u1i;
    if( *pu4X == 0)
    {
        for(u1i = 1; u1i<u1Len; u1i++)
        {
            if( *(pu4X + u1i) != 0)
            {
                return 1;
            }

        }
        return 0;
    }

    return 1;

}


/**********************************************************************
*      Degree_GF2m determines the degree of a polynomial in GF(2m).
*
* pu1X           : Address of the element
* Len            : Length in words of 32 bits.
*
**********************************************************************/


u2 Degree_GF2m(pu4 pu4X, u1 u1Len)
{
        u4 u4Last_Word = *(pu4X + u1Len - 1);
        u1 u1Count = 0;
        u1 u1degree = 0;

        while(u4Last_Word != 0)
        {
            u4Last_Word >>= 1;
            u1Count ++;
        }

        u1degree = (u1Len - 1)*32 + u1Count - 1;
        return u1degree;

}


/*********************************************************************
ConstToRamMemCpy() copies a buffer from EEPROM to RAM, while NOT changing the endianness

Input parameters :
  pu1Dest - pu4 - Input
       Pointer to the destination
  pu1Src - pfu1 * - Input
       Pointer to the source
  u2Length - u2 - Input
       Length of the buffer to be copied

Output values :
  none
  *******************************************************************/


void ConstToRamMemCpy(pu4 pu4Dest,pu4 pu4Src,u2 u2Length)
{
  u2  u2Cpt;
  pu4 pu4PtrDest;

  u2Cpt      = 0;
  pu4PtrDest = pu4Dest;
  while (u2Cpt < u2Length)
  {
    *(pu4PtrDest++) = pu4Src[u2Cpt];
    u2Cpt++;
  }

}



u1 word_to_byte(pu1 v_destination, pu4 v_initial, u1 len_v, u1 nb_coeff)
{
  uint8_t byte;
  for(int j=0;j<nb_coeff;j++)
  {

    for(int i=0; i<PARAM_M_LEN_WORD-1; i++)
    {
      byte = *(v_initial+i+PARAM_M_LEN_WORD*j) & 0xff;
      *(v_destination++) = byte;
      byte = (*(v_initial+i+PARAM_M_LEN_WORD*j)>>8 )& 0xff;
      *(v_destination++) = byte;
      byte = (*(v_initial+i+PARAM_M_LEN_WORD*j)>>16 ) & 0xff;
      *(v_destination++) = byte;
      byte = (*(v_initial+i+PARAM_M_LEN_WORD*j)>>24) & 0xff;
      *(v_destination++) = byte;
    }
    byte = *(v_initial+len_v-1+PARAM_M_LEN_WORD*j) & 0xff;
    *(v_destination++) = byte;
    byte = (*(v_initial+len_v-1+PARAM_M_LEN_WORD*j)>>8 )& 0xff;
    *(v_destination++) = byte;

  }
  return 0;

}
