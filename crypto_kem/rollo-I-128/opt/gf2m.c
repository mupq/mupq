#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "api.h"
#include "gf2m.h"
#include "utils.h"


/******************************************
*
*       Operations in GF(2^m)
*
*******************************************/

const u4 MODULO_GF2_M[3] = {0x00000201, 0x00000000, 0x00008000};

/**
 * \fn u2 Add_GF2m(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1Len)
 * \brief This function adds two elements in GF(2^m).
 *
 * \param[in]  pu4X : adress of the first polynomial.
 * \param[in]  pu4Y : adress of the second polynomial.
 * \param[in]  pu4R : adress of the result.
 * \param[in]  u1Len : number of 32-bits words in pu4X and pu4Y.
 */

u2 Add_GF2m(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1Len)
{
    u1 i;

    // XOR word per word.
    for(i = 0 ; i < u1Len ; i++)
        *(pu4R++) = *(pu4X++) ^ *(pu4Y++);
    return 0;
}



/**
 * \fn u2 Mupliply_GF2m_PreCompute(pu4 pu4X, pu4 pu4Storage, u1 u1Len)
 * \brief This function computes the pre-computation for polynomial multiplication in GF(2^m)
 *        as detailed in the next function.
 *
 * \param[in]  pu4X : adress of the polynomial.
 * \param[in]  pu4Storage : adress of the storage.
 *             /!\ Requires at least 16*u1Len 32-bits words.
 * \param[in]  u1Len : number of 32-bits words per element.
 */


u2 Mupliply_GF2m_PreCompute(pu4 pu4X, pu4 pu4Storage, u1 u1Len)
{
    Shift_32bits_Words_Left(pu4X, pu4Storage +     u1Len, u1Len, 0);                      // 0x1
    Shift_32bits_Words_Left(pu4X, pu4Storage +   2*u1Len, u1Len, 1);                      // 0x2
    Add_GF2m(pu4Storage +  u1Len, pu4Storage +   2*u1Len, pu4Storage +  3*u1Len, 3);      // 0x3
    Shift_32bits_Words_Left(pu4X, pu4Storage +   4*u1Len, u1Len, 2);                      // 0x4
    Add_GF2m(pu4Storage +  4*u1Len, pu4Storage +   u1Len, pu4Storage +  5*u1Len, 3);      // 0x5
    Add_GF2m(pu4Storage +  4*u1Len, pu4Storage + 2*u1Len, pu4Storage +  6*u1Len, 3);      // 0x6
    Add_GF2m(pu4Storage +    u1Len, pu4Storage + 2*u1Len, pu4Storage +  7*u1Len, 3);      // 0x7
    Add_GF2m(pu4Storage +  4*u1Len, pu4Storage + 7*u1Len, pu4Storage + 7*u1Len, 3);       // 0x7
    Shift_32bits_Words_Left(pu4X, pu4Storage +   8*u1Len, u1Len, 3);                      // 0x8
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +   u1Len, pu4Storage +  9*u1Len, 3);      // 0x9
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +   2*u1Len, pu4Storage +  10*u1Len, 3);   // 0xA
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +   2*u1Len, pu4Storage +  11*u1Len, 3);   // 0xB
    Add_GF2m(pu4Storage + 11*u1Len, pu4Storage +   u1Len, pu4Storage +  11*u1Len, 3);     // 0xB
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +  4*u1Len, pu4Storage +  12*u1Len, 3);    // 0xC
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +  4*u1Len, pu4Storage +  13*u1Len, 3);    // 0xD
    Add_GF2m(pu4Storage +    u1Len, pu4Storage +  13*u1Len, pu4Storage +  13*u1Len, 3);   // 0xD
    Add_GF2m(pu4Storage +  4*u1Len, pu4Storage +  8*u1Len, pu4Storage +  14*u1Len, 3);    // 0xE
    Add_GF2m(pu4Storage +  2*u1Len, pu4Storage +  14*u1Len, pu4Storage +  14*u1Len, 3);   // 0xE
    Add_GF2m(pu4Storage +  2*u1Len, pu4Storage +   u1Len, pu4Storage +  15*u1Len, 3);     // 0xF
    Add_GF2m(pu4Storage +  4*u1Len, pu4Storage +  15*u1Len, pu4Storage +  15*u1Len, 3);   // 0xF
    Add_GF2m(pu4Storage +  8*u1Len, pu4Storage +  15*u1Len, pu4Storage +  15*u1Len, 3);   // 0xF

    return 0;
}

/**
 * \fn u2 Mupliply_GF2m_With_PreCompute(pu4 pu4X, pu4 pu4PreCompute, pu4 pu4R, u1 u1Len, u1 u1Degree)
 * \brief This function processes the polynomial multiplication in GF(2^m) :
 *        Step 1 : For two polynmials P and Q of degree at most u1Degree, pre-compute C*Q for every value of Q contained in 4 bits.
 *        Step 2 : For every 4 bits value U of P, add U*Q to R with the right shift.
 *
 * \param[in]  pu4X : adress of the first polynomial.
 * \param[in]  pu4PreCompute : adress of the storage.
 *             /!\ Requires at least 16*u1Len 32-bits words.
 * \param[in]  pu4R : adress of the result
 * \param[in]  u1Len : number of 32-bits words in pu4X.
 * \param[in]  u1Degree : maximum degree of the polynomials.
 */

u2 Mupliply_GF2m_With_PreCompute(pu4 pu4X, pu4 pu4PreCompute, pu4 pu4R, u1 u1Len, u1 u1Degree)
{
    int i,j;
    u4 u4Mask = 0xF0000000;
    u1 u1Nb_32bits_Words = (u1Degree + 31) / 32;
    u4 u4Nibble = 0;

    for(i = 7 ; i >= 0 ; i-- )
    {
        for(j = 0 ; j < u1Nb_32bits_Words ; j++)
        {
            u4Nibble = ((*(pu4X + j)) & u4Mask) >> (4 * i);
            Add_GF2m(pu4PreCompute + u4Nibble*u1Len, pu4R + j, pu4R + j, u1Len);
        }

        u4Mask >>= 4;
        if(i != 0)
        {
            Shift_32bits_5Words_Left_4bits(pu4R);
        }
    }

    return 0;
}


/**
 * \fn u2 Modular_Reduction_GF2m(pu4 pu4X, pu4 pu4R)
 * \brief This function efficiently performs the modular reduction modulo P = X^79 + X^9 + 1
 *       This algorithm is designed specifically for ROLLO-I.
 *
 * \param[in]  pu4X : adress of the polynomial to reduce.
 * \param[in]  pu4R : adress of the result.
 * \param[in]  u1Len : number of 32-bits words per element.
 */
u2 Modular_Reduction_GF2m(pu4 pu4X, pu4 pu4R)
{

    u4 u4MSBStorage[8] = { (*(pu4X + 4)) >> 6,
                           (*(pu4X + 4)) >> 15,
                           (*(pu4X + 4)) << 17,
                           (*(pu4X + 4)) << 26,
                           (*(pu4X + 3)) >>  6,
                           (*(pu4X + 3)) >> 15,
                           (*(pu4X + 3)) << 17,
                           (*(pu4X + 3)) << 26
                         };

    *(pu4X + 2) ^=  u4MSBStorage[0] ^ u4MSBStorage[1];
    *(pu4X + 1) ^=  u4MSBStorage[2] ^ u4MSBStorage[3] ^ u4MSBStorage[4] ^ u4MSBStorage[5];
    * pu4X      ^=  u4MSBStorage[6] ^ u4MSBStorage[7];

    // Extract bits 15 to 31 of C[2]

    u4 u4Extract  = *(pu4X + 2) & 0xFFFF8000;  // 0xFFF80000
     * pu4X      ^= u4Extract >> 15;
     *(pu4X ++)  ^= u4Extract >> 6;
     *(pu4X ++)  ^= u4Extract << 22;
     *(pu4X ++)  &= 0x00007FFF;               // Clear MSB


     if((pu4X - 3) != pu4R)
     {
         memcpy(pu4R, pu4X - 3, PARAM_M_LEN_WORD*4);
         memset(pu4X - 3, 0x00, 20);
     }

     *(pu4R + 3) = 0x00;
     *(pu4R + 4) = 0x00;
     *(pu4R + 5) = 0x00;

    return 0;
}


/**
* \fn u2 Modular_inverse_GF2m(pu4 pu4X, pu4 pu4Modulo, pu4 pu4R, pu4 pu4Workspace, u1 u1Len)
* \brief  This function computes the inverse of an element in GF(2^m).
*
* \param[in]  pu4X : adress of the polynomial to inverse.
* \param[in]  pu4Modulo : adress of the modulo.
* \param[in]  pu4R : adress of the result.
* \param[in]  pu4Workspace : adress of workspace.
* \param[in]  u1Len : number of 32-bits words in pu4X.
*/

u2 Modular_inverse_GF2m(pu4 pu4X, pu4 pu4Modulo, pu4 pu4R, pu4 pu4Workspace, u1 u1Len)
{

    ConstToRamMemCpy(pu4Modulo,
                     (pu4) MODULO_GF2_M,
                     MODULO_GF2M_LEN);

    pu4 pu4U;
    pu4 pu4V;
    pu4 pu4G2     = pu4Workspace + u1Len;
    pu4 pu4Temp   = pu4G2 + u1Len;
    pu4 pu4Vtemp  = pu4Temp + u1Len;
    pu4 pu4G2temp = pu4Vtemp + u1Len;

    pu4V = pu4Modulo;
    pu4U = pu4X;
    u1 j;
    u4 u4Count = 0;

    memset(pu4R, 0, u1Len*4);
    memset(pu4G2, 0, u1Len*4);
    *(pu4R) = 0x01;
    *(pu4G2) = 0x00;


    u1 uLen = u1Len;
    u1 vLen = MODULO_GF2M_LEN;
    u1 u1Lentmp;


    u1 uDegree = Degree_GF2m(pu4U, uLen);
    u1 vDegree = Degree_GF2m(pu4V, vLen);


    while(Is_Unit(pu4U, uLen) != 0)
    {

        if(uDegree>=vDegree)
        {
            j = uDegree - vDegree;
        }
        else
        {

            /* Swap between u and v as well as the length in words of u and v */

            pu4Temp = pu4V;
            pu4V = pu4U;
            pu4U = pu4Temp;

            u1Lentmp = uLen;
            uLen = vLen;
            vLen = u1Lentmp;

            /* Swap between g1 and g2  */

            pu4Temp = pu4R;
            pu4R = pu4G2;
            pu4G2 = pu4Temp;
            j = vDegree - uDegree;
            u4Count += 1;
        }

        /* These two shifts amount to multiply v and g2 by x^j */
        Shift_32bits_Words_Left(pu4V,   pu4Vtemp, u1Len, j);
        Shift_32bits_Words_Left(pu4G2, pu4G2temp, u1Len, j);

        /* u = u  x^j * v  and g1 = g1 + x^j * g2*/
        Add_GF2m(pu4U,   pu4Vtemp, pu4U, u1Len);
        Add_GF2m(pu4R, pu4G2temp,  pu4R, u1Len);

        /* Determination of the degree of u and v*/
        if(pu4U[uLen-1] == 0)
        {
            uLen = uLen - 1;
        }
        uDegree = Degree_GF2m(pu4U, uLen);

        if(pu4V[vLen-1] == 0)
        {
            vLen = vLen - 1;
        }
        vDegree = Degree_GF2m(pu4V,vLen);

    }
    if((u4Count %2) == 1)
    {
      memcpy(pu4G2, pu4R, MODULO_GF2M_LEN*4);
      memcpy(pu4X,pu4U,MODULO_GF2M_LEN*4);
    }


    return 0;
}
