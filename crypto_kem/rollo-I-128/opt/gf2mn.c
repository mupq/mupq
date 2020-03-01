#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "api.h"
#include "matrix.h"
#include "gf2mn.h"
#include "gf2m.h"
#include "utils.h"


/*********************************************
*
*       Operations in GF(2^m)^n
*
**********************************************/

const u1 MODULO_GF2_N[2] = {0x00, 0x05};


/**
 * \fn u1 leading_coefficient_div(pu4 pu4R, pu4 pu4V, pu4 pu4U, u1 degree_v, u1 degree_u, pu4 workspace)
 * \brief  This functions is used in the computation of the inverse of a polynomial V modulo a polynomial U,
 *         it computes the multiplication between the inverse of the V leading coefficient and the leading coefficient U.
 * \param[out]  pu4R : adress of the result pointer pu4R = pu4V_(n-1)^-(1)* pu4U_(n-1).
 * \param[in]   pu4V : adress of the divisor polynomial.
 * \param[in]   pu4U : adress of the dividend polynomial.
 * \param[in]   degree_v : degree of the polynomial pu4V in GF(2^m)^n.
 * \param[in]   degree_u : degree of the polynomial pu4U in GF(2^m)^n.
 * \param[in]   workspace : address of the workspace pointer.
 */

u1 leading_coefficient_div(pu4 pu4R, pu4 pu4V, pu4 pu4U, u1 degree_v, u1 degree_u, pu4 workspace)
{
    //definition of tempory buffers
    pu4 element_to_inverse = workspace;
    pu4 inverse = element_to_inverse + PARAM_M_LEN_WORD;
    pu4 base_modulo = inverse + PARAM_M_LEN_WORD;
    pu4 workspace_inverse = base_modulo + MODULO_GF2M_LEN;
    pu4 base_mult = inverse + MODULO_GF2M_LEN;

    // Step 1 : Compute modular inverse of Leading Coefficient (LC) pu4V.
    memcpy(element_to_inverse,
           pu4V + PARAM_M_LEN_WORD*degree_v,
           PARAM_M_LEN_WORD*4);

    Modular_inverse_GF2m(element_to_inverse,
                         base_modulo,
                         inverse,
                         workspace_inverse,
                         MODULO_GF2M_LEN);

    memset(base_modulo, 0x00, PARAM_M_LEN_WORD*4);
    memset(workspace_inverse, 0x00, 5*PARAM_M_LEN_WORD*PARAM_M_LEN_WORD*4);

    // Step 2 :multiplication between the modular inverse of pu4V's LC and the pu4U's LC.
    Mupliply_GF2m_PreCompute(inverse,
                             base_mult,
                             MODULO_GF2M_LEN);



    Mupliply_GF2m_With_PreCompute(pu4U + PARAM_M_LEN_WORD*degree_u,
                                  base_mult,
                                  pu4R,
                                  PARAM_M_LEN_WORD,
                                  PARAM_M);

    // Step 3 : Modular reduction of the result in GF(2^m)^n
    Modular_Reduction_GF2m(pu4R,
                           pu4R);

    memset(workspace,0x00, PARAM_M_LEN_WORD*8*4);
    memset(base_mult, 0x00, PARAM_M_LEN_WORD*16*4);
    return 0;

}

/**
 * \fn u1 mult_by_element(pu4 pu4R, pu4 pu4V, pu4 pu4scalar, u1 nb_coeff, pu4 workspace, u1 nb_shift)
 * \brief  This function computes the mutliplication between a polynomial pu4V in GF(2^m)^n and an element in GF(2^m)
 * \param[out]  pu4R : adress of the resulting multiplication.
 * \param[in]   pu4V : adress of the polynomial to multiply.
 * \param[in]   pu4scalar : adress of the element in GF(2^m).
 * \param[in]   workspace : address of the workspace pointer.
 * \param[in]   nb_shift : number of shift to apply to the polynomial pu4V after the multiplication (= degree of the element pu4scalar)
**/
u1 mult_by_element(pu4 pu4R, pu4 pu4V, pu4 pu4scalar, u1 nb_coeff, pu4 workspace, u1 nb_shift)
{
  pu4 base_mult = workspace;

  //Precomputation required for the multiplication in GF(2^m)
  Mupliply_GF2m_PreCompute(pu4scalar,
                          base_mult,
                          MODULO_GF2M_LEN);

  for (int u1i = 0  ;  u1i < nb_coeff ; u1i++)
    {
        Mupliply_GF2m_With_PreCompute(pu4V + MODULO_GF2M_LEN*u1i,
                                      base_mult,
                                      pu4R+ nb_shift*MODULO_GF2M_LEN+ MODULO_GF2M_LEN*u1i,
                                      MODULO_GF2M_LEN,
                                      PARAM_M);
        //Modular reduction is applied on each coefficient the resulting multiplication
        Modular_Reduction_GF2m(pu4R+ nb_shift*MODULO_GF2M_LEN+ MODULO_GF2M_LEN*u1i,
                               pu4R+ nb_shift*MODULO_GF2M_LEN+ MODULO_GF2M_LEN*u1i);

    }
    memset(workspace,0x00, PARAM_M_LEN_WORD*16*4);
  return 0;
}


/**
 * \fn u1 init_modulo(pu4 pu4v, int n)
 * \brief  This function inits the modulo in memory
 * \param[out]   pu4v : address of the modulo.
 * \param[in]       n : degree of the modulo.
**/
u1 init_modulo(pu4 pu4v, int n)
{
  if(n==47)
  {
    *(pu4v) = 0x1;
    *(pu4v + 15 ) = 0x1;
    *(pu4v + 141) = 0x1;
  }
  return 0;
}


/**
 * \fn u1 modular_inverse_GF2mn(pu4 pu4R, pu4 pu4X, pu4 pu4modulo)
 * \brief  This function computes the modular inverse of an element in GF(2^m)^n
 * \param[out]     pu4R : address of the output inverse.
 * \param[in]      pu4X : address of the element to inverse.
 * \param[in] pu4modulo : address of the modulo in GF(2^m)^n.
**/

u1 modular_inverse_GF2mn(pu4 pu4R, pu4 pu4X, pu4 pu4modulo)
{

  pu4 u  = pu4X;
  memcpy(u, pu4X, PARAM_N*PARAM_M_LEN_WORD*4);
  u1 u_degree   = Degree_GF2mn(pu4X, PARAM_M_LEN_WORD, PARAM_N-1);
  u1 v_degree  = PARAM_N;
  u1 u1NbCoeff_u = u_degree+1;
  u1 u1NbCoeff_v = v_degree+1;
  u1 u1NbCoeff_g1 = 1;
  u1 u1NbCoeff_g2 = u1NbCoeff_g1;
  u1 nb_shift;
  pu4 v  = pu4modulo;
  u1 degree_tmp;
  pu4 g1  = v + u1NbCoeff_v*PARAM_M_LEN_WORD;
  pu4 g2 = g1 + u1NbCoeff_g1*PARAM_M_LEN_WORD;
  pu4 storage_tmp = g2 + u1NbCoeff_g2*PARAM_M_LEN_WORD;
  pu4 pu4vtemp = storage_tmp + 2*PARAM_M_LEN_WORD;
  pu4 pu4g2temp = pu4vtemp + u1NbCoeff_v*PARAM_M_LEN_WORD;
  pu4 workspace_mult = pu4g2temp + u1NbCoeff_g2*PARAM_M_LEN_WORD;
  g1[0] = 0x01;



  while(Is_Zero(u, u1NbCoeff_u*PARAM_M_LEN_WORD)!=0)
  {
    if(u_degree<v_degree)

    {
      memcpy(storage_tmp, v, PARAM_M_LEN_WORD*u1NbCoeff_v*4);
      memset(v , 0x00, u1NbCoeff_v*PARAM_M_LEN_WORD*4);
      v  = u + u1NbCoeff_v*PARAM_M_LEN_WORD;
      memcpy(v, u, PARAM_M_LEN_WORD*u1NbCoeff_u*4);
      memcpy(u, storage_tmp, PARAM_M_LEN_WORD*u1NbCoeff_v*4);
      memset(storage_tmp,0x00, PARAM_M_LEN_WORD*u1NbCoeff_v*4);


      degree_tmp = u_degree;
      u_degree = v_degree;
      v_degree = degree_tmp;
      u1NbCoeff_u = u_degree+1;
      u1NbCoeff_v = v_degree+1;


      memcpy(storage_tmp, g1, PARAM_M_LEN_WORD*u1NbCoeff_g1*4);
      memset(g1,0x00, PARAM_M_LEN_WORD*u1NbCoeff_g1*4);
      g1  = v + u1NbCoeff_v*PARAM_M_LEN_WORD;
      memcpy(g1, g2, PARAM_M_LEN_WORD*u1NbCoeff_g2*4);
      g2 = g1 + u1NbCoeff_g1*PARAM_M_LEN_WORD;
      memcpy(g2, storage_tmp, PARAM_M_LEN_WORD*u1NbCoeff_g1*4);
      memset(storage_tmp,0x00, PARAM_M_LEN_WORD*u1NbCoeff_g1*4);
      u1NbCoeff_g2 = u1NbCoeff_g1;

    }

    nb_shift = u1NbCoeff_u - u1NbCoeff_v;


    u1NbCoeff_g1+= nb_shift;
    storage_tmp = g2 + u1NbCoeff_g1*PARAM_M_LEN_WORD;
    memcpy(storage_tmp, g2, PARAM_M_LEN_WORD*u1NbCoeff_g2*4);
    memset(g2,0x00, PARAM_M_LEN_WORD*u1NbCoeff_g2*4);
    g2 = g1 + u1NbCoeff_g1*PARAM_M_LEN_WORD;
    memcpy(g2, storage_tmp, PARAM_M_LEN_WORD*u1NbCoeff_g2*4);
    memset(storage_tmp,0x00, PARAM_M_LEN_WORD*u1NbCoeff_g2*4);
    pu4g2temp =g2 + u1NbCoeff_g2*PARAM_M_LEN_WORD;
    pu4vtemp = storage_tmp + 2*PARAM_M_LEN_WORD;
    pu4g2temp = pu4vtemp + u1NbCoeff_v*PARAM_M_LEN_WORD;
    workspace_mult = pu4g2temp + (u1NbCoeff_g1+2)*PARAM_M_LEN_WORD;

    leading_coefficient_div(storage_tmp, v, u, v_degree, u_degree, workspace_mult);
    mult_by_element(pu4vtemp, v, storage_tmp,u1NbCoeff_v,workspace_mult,nb_shift);
    Add_GF2m(u, pu4vtemp, u, u1NbCoeff_u*PARAM_M_LEN_WORD);
    memset(pu4vtemp, 0x00, u1NbCoeff_u*PARAM_M_LEN_WORD*4);
    mult_by_element(pu4g2temp,g2, storage_tmp,u1NbCoeff_g2,workspace_mult,nb_shift);
    memset(storage_tmp,0x00,PARAM_M_LEN_WORD*4);
    Add_GF2m(g1, pu4g2temp, g1, u1NbCoeff_g1*PARAM_M_LEN_WORD);

    memset(pu4g2temp, 0x00, u1NbCoeff_g1*PARAM_M_LEN_WORD*4);

    while(!Is_Zero(u + (u1NbCoeff_u*PARAM_M_LEN_WORD-PARAM_M_LEN_WORD),PARAM_M_LEN_WORD))
    {
      u1NbCoeff_u = u1NbCoeff_u - 1;
      if(u1NbCoeff_u==0)
      {
        u_degree = u1NbCoeff_u;
        break;
      }
      else
      {
        u_degree = u1NbCoeff_u - 1;
      }
    }

  }
  if(Is_Unit(v, MODULO_GF2M_LEN))
  {
    pu4 inverse = workspace_mult + 2*MODULO_GF2M_LEN;
    pu4 pu4Rtmp = g2 + u1NbCoeff_g2*PARAM_M_LEN_WORD;
    Modular_inverse_GF2m(v,
                         workspace_mult + MODULO_GF2M_LEN ,
                         inverse,
                         workspace_mult + 7*MODULO_GF2M_LEN,
                         MODULO_GF2M_LEN);
    mult_by_element(pu4Rtmp, g2, inverse, u1NbCoeff_g2,workspace_mult+16*PARAM_M_LEN_WORD*MODULO_GF2M_LEN, 0);
    memcpy(pu4R, pu4Rtmp, u1NbCoeff_g2*MODULO_GF2M_LEN*4);
    memset(pu4Rtmp,0x00, u1NbCoeff_g2*MODULO_GF2M_LEN*4);
    memset(pu4R+u1NbCoeff_g2*PARAM_M_LEN_WORD, 0x00, u1NbCoeff_g2*PARAM_M_LEN_WORD*4);

  }
  else
  {
    memcpy(pu4R, g2, u1NbCoeff_g2*MODULO_GF2M_LEN*4);
  }

 return 0;
}


/**
 * \fn u1 Degree_GF2mn(pu4 pu4Polynom, u1 u1CoeffLength, u1 u1DegreeMaxPolynom)
 * \brief  This function returns the degree of an element in GF(2^m)^n, this amounts to return the highest non zero coefficient.
 * \param[in]         pu4Polynom : address of an element in GF(2^m)^n
 * \param[in]      u1CoeffLength : length of pu4Polynom's coefficient
 * \param[in] u1DegreeMaxPolynom : maximum degree of pu4Polynom
 * \param[out]  degree of  pu4Polynom
**/

u1 Degree_GF2mn(pu4 pu4Polynom, u1 u1CoeffLength, u1 u1DegreeMaxPolynom)
{

  u1 u1i = 0;
  u1 u1Flag;

  // The latest non zero coefficient is finding from the maximum degree of pu4Polynom (=u1DegreeMaxPolynom).
  pu4 pu4Polynom_temp = pu4Polynom + u1CoeffLength*u1DegreeMaxPolynom;

  while(1)
  {
    //If the processing coefficient is the constant term, it returns 0
    if ( (pu4Polynom_temp - u1i* u1CoeffLength) == pu4Polynom)
    {
      return 0;
    }

    // Comparison with the zero element
    u1Flag = Is_Zero(pu4Polynom_temp - u1i* u1CoeffLength, u1CoeffLength);

    // If the processing element is non zero, we thus get the degree
    if (u1Flag == 1)
    {
    	return u1DegreeMaxPolynom - u1i;
    }

    //Case of the zero polynomial (we reach the maximum degree and all coefficients are zero)
    if (u1i == u1DegreeMaxPolynom)
      break;
    u1i += 1;
  }
  return 0;
}


/**
 * \fn u1 Add_GF2mn(u1 u1NbCoeff, u1 u1CoeffLength, pu4 pu4X, pu4 pu4Y, pu4 pu4R)
 * \brief  This function adds two elements in GF(2^m)^n.
 * \param[in]      u1NbCoeff : Number of coefficient in X and Y.
 * \param[in]  u1CoeffLength : Length of coefficient in GF(2^m)^n.
 * \param[in]           pu4X : address of the first element in GF(2^m)^n implied in the addition.
 * \param[in]           pu4Y : address of the second element in GF(2^m)^n.
 * \param[out]          pu4R : address of the addition result R = X + Y in GF(2^m)^n.
**/

u1 Add_GF2mn(u1 u1NbCoeff, u1 u1CoeffLength, pu4 pu4X, pu4 pu4Y, pu4 pu4R)
{
    // Simply add in GF(2^m) every coefficient of the same degree.
    Add_GF2m(pu4X, pu4Y, pu4R, u1NbCoeff*u1CoeffLength);
    return 0;

}


/**
 * \fn u1 Mupliply_GF2mn(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY, u1 u1ModuloDegree, u1 u1CoeffLength)
 * \brief  This function multiplies two elements in GF(2^m)^n. For a memory usage gain, the modular reduction is directly applied during the multiplication.
 * \param[in]         pu4X : address of the first element in GF(2^m)^n implied in the mutliplication.
 * \param[in]         pu4Y : address of the second element in GF(2^m)^n.
 * \param[out]        pu4R : address of the addition result R = X x Y in GF(2^m)^n.
 * \param[in]      u1NbCoeffX : Number of coefficient in X.
 * \param[in]      u1NbCoeffY : Number of coefficient in Y.
 * \param[in]      u1ModuloDegree : Degree of the modulo in GF(2^m)^n.
 * \param[in]  u1CoeffLength : Length of coefficient in GF(2^m)^n.
**/


u1 Mupliply_GF2mn(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY, u1 u1ModuloDegree, u1 u1CoeffLength)
{
    u1 u1i, u1j, u1k, u1l;
    u1 u1IndexAfterReduction;
    u1 u1ReductionDecalage;

    u4 BASE_MULT_PRE[16*PARAM_M_LEN_WORD];
    u4 BASE_MULT[2*PARAM_M_LEN_WORD];
    u4 BASE_ADD[PARAM_M_LEN_WORD];

    for( u1i = 0 ; u1i < u1NbCoeffX ; u1i++)
    {
        memset(BASE_MULT_PRE,
            0x00,
            16*PARAM_M_LEN_WORD*4);
        // Precomputation for pu4X + u1i*u1CoeffLength
        Mupliply_GF2m_PreCompute(pu4X + u1i * MODULO_GF2M_LEN,
                                 BASE_MULT_PRE,
                                 MODULO_GF2M_LEN);


        for( u1j = 0 ; u1j < u1NbCoeffY ; u1j++)
        {
            memset(BASE_MULT,
                   0x00,
                   2*PARAM_M_LEN_WORD*4);
            memset(BASE_ADD,
                   0x00,
                   PARAM_M_LEN_WORD*4);
            // Multiply with result in buffer
            Mupliply_GF2m_With_PreCompute(pu4Y + u1j*MODULO_GF2M_LEN,
                                          BASE_MULT_PRE,
                                          BASE_MULT,
                                          MODULO_GF2M_LEN,
                                          PARAM_M);


            // Apply modular reduction
            Modular_Reduction_GF2m(BASE_MULT,
                                   BASE_MULT);


            // Apply modular reduction in GF(2^u1CoeffLength)^u1NbCoeff
            if( u1i + u1j >= u1ModuloDegree)
            {
                u1IndexAfterReduction = (u1i + u1j) % u1ModuloDegree;

                for (u1k = 0 ; u1k < NB_COEFF_MODULO_GF2_N ; u1k++)
                {
                    u1ReductionDecalage = MODULO_GF2_N[u1k];
                    Add_GF2m(pu4R + (u1IndexAfterReduction + u1ReductionDecalage)*MODULO_GF2M_LEN,
                             BASE_MULT,
                             BASE_ADD,
                             MODULO_GF2M_LEN);


                    // Degree is still higher than modulo's, apply another reduction
                    if ( (((u1i + u1j) % u1ModuloDegree) + MODULO_GF2_N[u1k]) >= u1ModuloDegree )
                    {

                        for( u1l = 0 ; u1l < NB_COEFF_MODULO_GF2_N ; u1l++ )
                        {
                            Add_GF2m(BASE_ADD,
                                     pu4R + (((u1IndexAfterReduction + u1ReductionDecalage) % u1ModuloDegree) + MODULO_GF2_N[u1l])*MODULO_GF2M_LEN,
                                     pu4R + (((u1IndexAfterReduction + u1ReductionDecalage) % u1ModuloDegree) + MODULO_GF2_N[u1l])*MODULO_GF2M_LEN,
                                     MODULO_GF2M_LEN);

                        }

                    }

                    else
                    {
                      memcpy(pu4R + (u1IndexAfterReduction + u1ReductionDecalage)*u1CoeffLength,
                             BASE_ADD,
                             u1CoeffLength*4);

                    }


                }
            }

            // If we not overpass the modulo's degree we only have to copy the value on the right polynomial result's degree.
            else
            {
                Add_GF2m(pu4R + (u1i + u1j)*MODULO_GF2M_LEN,
                         BASE_MULT,
                          pu4R + (u1i + u1j)*MODULO_GF2M_LEN,
                         MODULO_GF2M_LEN);
            }


        }

    }
    return 0;

}



/**
 * \fn u2 SchoolMultiply(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY)
 * \brief  This function multiplies two elements of GF(2^m)^n without modular reduction, the result is in GF(2^m)^(2n).
 * \param[in]         pu4X : address of the first element in GF(2^m)^n implied in the mutliplication.
 * \param[in]         pu4Y : address of the second element in GF(2^m)^n.
 * \param[out]        pu4R : address of the addition result R = X x Y in GF(2^m)^(2n).
 * \param[in]      u1NbCoeffX : Number of coefficient in X.
 * \param[in]      u1NbCoeffY : Number of coefficient in Y.
**/

u2 SchoolMultiply(pu4 pu4X, pu4 pu4Y, pu4 pu4R, u1 u1NbCoeffX, u1 u1NbCoeffY)
{
    u1 u1i, u1j;

    u4 BASE_MULT_PRE[16*PARAM_M_LEN_WORD];
    u4 BASE_MULT[2*PARAM_M_LEN_WORD];

    for( u1i = 0 ; u1i < u1NbCoeffX ; u1i++)
    {
      memset(BASE_MULT_PRE,
            0x00,
            16*PARAM_M_LEN_WORD*4);
        // Precomputation for pu4X + u1i*u1CoeffLength
        Mupliply_GF2m_PreCompute(pu4X + u1i * MODULO_GF2M_LEN,
                                 BASE_MULT_PRE,
                                 MODULO_GF2M_LEN);

        for( u1j = 0 ; u1j < u1NbCoeffY ; u1j++)
        {
            memset(BASE_MULT,
                   0x00,
                   2*PARAM_M_LEN_WORD*4);
            // Multiply with the result in buffer
            Mupliply_GF2m_With_PreCompute(pu4Y + u1j*MODULO_GF2M_LEN,
                                          BASE_MULT_PRE,
                                          BASE_MULT,
                                          MODULO_GF2M_LEN,
                                          PARAM_M);
            // Apply modular reduction
            Modular_Reduction_GF2m(BASE_MULT,
                                   BASE_MULT);


            Add_GF2m(pu4R + (u1i + u1j)*MODULO_GF2M_LEN,
                     BASE_MULT,
                     pu4R + (u1i + u1j)*MODULO_GF2M_LEN,
                     MODULO_GF2M_LEN);

        }

    }

    return 0;
}


/**
 * \fn u2 ModularReduction(pu4 pu4X, pu4 pu4R, u1 u1NbCoeffX, u1 u1ModuloDegree, u1 u1CoeffLength)
 * \brief  This function reduced an element of GF(2^m)^(kn) into an element of GF(2^m)^n.
 * \param[in/out]          pu4X : address of the element X in GF(2^m)^(kn) to reduce.
 * \param[in]      u1NbCoeffX : Number of coefficient in X.
 * \param[in]      u1ModuloDegree : Degree of the modulo in GF(2^m)^n.
* \param[in]       u1CoeffLength : Length of coefficient in GF(2^m)^n.
**/

u2 ModularReduction(pu4 pu4X, u1 u1NbCoeffX, u1 u1ModuloDegree, u1 u1CoeffLength)
{
    u1 u1i;


    for( u1i = u1NbCoeffX-1 ; u1i >= u1ModuloDegree ; u1i--)
    {
       Add_GF2m(pu4X + (u1i + MODULO_GF2_N[0] - u1ModuloDegree)*MODULO_GF2M_LEN,
                pu4X + u1i*MODULO_GF2M_LEN,
                pu4X + (u1i + MODULO_GF2_N[0] - u1ModuloDegree)*MODULO_GF2M_LEN,
                MODULO_GF2M_LEN);

       Add_GF2m(pu4X + (u1i + MODULO_GF2_N[1] -  u1ModuloDegree)*MODULO_GF2M_LEN,
                pu4X + u1i*MODULO_GF2M_LEN,
                pu4X + (u1i+ MODULO_GF2_N[1] - u1ModuloDegree)*MODULO_GF2M_LEN,
                MODULO_GF2M_LEN);
       memset(pu4X + u1i*u1CoeffLength, 0x00, u1CoeffLength*4);
    }



    return 0;

}


/**
 * \fn u2 Karatsuba(pu4 pu4X, pu4 pu4Y, pu4 pu4R, pu4 pu4Workspace, u1 u1Degree, u1 u1CoeffLength)
 * \brief  This function multiplies two element of GF(2^m)^n, it is combined with School book method.
 * \param[in]          pu4X : address of the element X in GF(2^m)^(kn) to reduce.
 * \param[out]         pu4R : address the resulting reduction.
 * \param[out]        pu4R : address of the addition result R = X x Y in GF(2^m)^(2n).
 * \param[in]      u1NbCoeffX : Number of coefficient in X.
 * \param[in]      u1ModuloDegree : Degree of the modulo in GF(2^m)^n.
* \param[in]       u1CoeffLength : Length of coefficient in GF(2^m)^n.
**/


u2 Karatsuba(pu4 pu4X, pu4 pu4Y, pu4 pu4R, pu4 pu4Workspace, u1 u1Degree, u1 u1CoeffLength)
{
  u1 u1i;



  pu4 pu4R2;
  pu4 pu4R3;
  pu4 pu4R4;
  pu4 pu4Workspace2;
  pu4 pu4X2;
  pu4 pu4Y2;

  if( (u1Degree & 1) != 0)
  {
    SchoolMultiply(pu4X, pu4Y, pu4R, u1Degree, u1Degree);
    return 0;
  }

  u1 const Half = u1Degree >> 1;

  pu4R2      = pu4R + Half*u1CoeffLength;
  pu4R3      = pu4R + u1Degree*u1CoeffLength;
  pu4R4      = pu4R + (u1Degree + Half)*u1CoeffLength;
  pu4X2      = pu4X + Half*u1CoeffLength;
  pu4Y2      = pu4Y + Half*u1CoeffLength;
  pu4Workspace2 = pu4Workspace + Half*u1CoeffLength;


   Add_GF2mn(Half, u1CoeffLength, pu4X, pu4X2, pu4R);
   Add_GF2mn(Half, u1CoeffLength, pu4Y, pu4Y2, pu4R2);


  memset(pu4R3, 0x00, Half*u1CoeffLength*4);
  memset(pu4Workspace, 0x00, 2*Half*u1CoeffLength*4);
  Karatsuba(pu4R,  pu4R2, pu4Workspace, pu4R3, Half, u1CoeffLength);

  memset(pu4R, 0x00, Half*u1CoeffLength*4);
  memset(pu4R3, 0x00, 2*Half*u1CoeffLength*4);
  Karatsuba(pu4X2, pu4Y2, pu4R3,        pu4R,  Half, u1CoeffLength);


  for(u1i = 0 ; u1i < Half ; u1i++)
  {
      Add_GF2m(pu4Workspace  + u1i*MODULO_GF2M_LEN,
               pu4R3         + u1i*MODULO_GF2M_LEN,
               pu4Workspace  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);
  }


  memcpy(pu4R2, pu4Workspace, Half*u1CoeffLength*4);


  for(u1i = 0 ; u1i < Half ; u1i++)
  {
      Add_GF2m(pu4Workspace2  + u1i*MODULO_GF2M_LEN,
               pu4R4         + u1i*MODULO_GF2M_LEN,
               pu4Workspace2  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);

      Add_GF2m(pu4R3  + u1i*MODULO_GF2M_LEN,
               pu4Workspace2         + u1i*MODULO_GF2M_LEN,
               pu4R3  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);
  }

  memset(pu4R, 0x00, Half*u1CoeffLength*4);
  memset(pu4Workspace, 0x00, 2*Half*u1CoeffLength*4);
  Karatsuba(pu4X, pu4Y, pu4Workspace, pu4R, Half, u1CoeffLength);


  memcpy(pu4R, pu4Workspace, Half*u1CoeffLength*4);
  for(u1i = 0 ; u1i < Half ; u1i++)
  {
      Add_GF2m(pu4R2  + u1i*MODULO_GF2M_LEN,
               pu4Workspace         + u1i*MODULO_GF2M_LEN,
               pu4R2  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);

      Add_GF2m(pu4R2  + u1i*MODULO_GF2M_LEN,
               pu4Workspace2         + u1i*MODULO_GF2M_LEN,
               pu4R2  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);

      Add_GF2m(pu4R3  + u1i*MODULO_GF2M_LEN,
               pu4Workspace2         + u1i*MODULO_GF2M_LEN,
               pu4R3  + u1i*MODULO_GF2M_LEN,
               MODULO_GF2M_LEN);
  }
  return 0;
}


/**
 * \fn u1 Karatsuba_Degree(u1 u1Degree)
 * \brief  This function returns the padding degree for Karatsuba multiplication.
 * \param[in]          u1Degree : current degree.
 * \param[out]         corresponding padding degree.
**/

u1 Karatsuba_Degree(u1 u1Degree)
{
  switch(u1Degree)
  {
    case 47 :
      return 48;
    case 53 :
      return 56;
  }
  return 0;
}


/**
 * \fn u1 S_i_calculus(pu4 pu4F, pu4 pu4S, pu4 pu4R, u1 u1Index, u1 u1SDim)
 * \brief  This function computes S_i = (f_i)^(-1) S required during RSR algorithm.
 * \param[in]          pu4F : address of the private support F.
 * \param[in]          pu4S : address the syndrome.
 * \param[out]         pu4R : address of the result R = S_i.
 * \param[in]       u1Index : Index i of S_i.
 * \param[in]        u1SDim : syndrome dimension.
**/

u1 S_i_calculus(pu4 pu4F, pu4 pu4S, pu4 pu4R, u1 u1Index, u1 u1SDim)
{
  u1 u1CoeffLength = PARAM_M_LEN_WORD;
  u1 u1i;
  u4 WORKSPACE_INVERSE[7*PARAM_M_LEN_WORD];
  u4 BASE_INVERSE[PARAM_M_LEN_WORD];
  u4 BASE_INVERSE_RESULT[PARAM_M_LEN_WORD];
  u4 BASE_MODULO_M[PARAM_M_LEN_WORD];
  u4 BASE_MULT_PRE[16*PARAM_M_LEN_WORD];
  u4 BASE_MULT[2*PARAM_M_LEN_WORD];

  // Step 1 : Inversion of the coefficient F[u1Index]

  memset(WORKSPACE_INVERSE,0x00,7*PARAM_M_LEN_WORD*4);
  memset(BASE_INVERSE,0x00,PARAM_M_LEN_WORD*4);
  memset(BASE_INVERSE_RESULT,0x00,PARAM_M_LEN_WORD*4);

  memcpy(BASE_INVERSE,
         pu4F + u1CoeffLength*u1Index,
         u1CoeffLength*4);

  Modular_inverse_GF2m(pu4F + MODULO_GF2M_LEN*u1Index,
                       BASE_MODULO_M,
                       BASE_INVERSE_RESULT,
                       WORKSPACE_INVERSE,
                       MODULO_GF2M_LEN);


  // Step 2 : S_i computation

  // Precomputation for modular inverse of divisor's leading coefficient
  memset(BASE_MULT_PRE, 0x00, 16*PARAM_M_LEN_WORD*4);
  memset(BASE_MULT, 0x00, 2*PARAM_M_LEN_WORD*4);
  Mupliply_GF2m_PreCompute(BASE_INVERSE_RESULT,
                           BASE_MULT_PRE,
                           MODULO_GF2M_LEN);

  for(u1i = 0 ; u1i < u1SDim ; u1i++)
  {
      // Multiplication with the result in buffer
      Mupliply_GF2m_With_PreCompute(pu4S + u1i*MODULO_GF2M_LEN,
                                    BASE_MULT_PRE,
                                    BASE_MULT,
                                    MODULO_GF2M_LEN,
                                    PARAM_M);
      // Application of the modular reduction
      Modular_Reduction_GF2m(BASE_MULT,
                              pu4R + u1i*MODULO_GF2M_LEN);
  }

  memcpy(pu4F + u1CoeffLength*u1Index,
         BASE_INVERSE,
         u1CoeffLength*4);

  return 0;
}



/**
 * \fn u1 Intersection(pu4 pu4V1, pu4 pu4V2, pu4 pu4Workspace, u1 u1Dim, u1 u1CoeffLength, u1 u1ModuloDegree)
 * \brief  This function computes the intersection's base between two space vectors - Application of Zassenhaus method.
 * \param[in]          pu4V1 : address of the first space vector.
 * \param[in]          pu4V2 : address of the second space vector.
 * \param[out]          pu4R : address of the result R = pu4V1 inter pu4V2.
 * \param[in]          pu4V1Dim : dimension of pu4V1.
 * \param[in]          pu4V2Dim : dimension of pu4V2.
 * \param[in]     u1CoeffLength : length of a coefficient.
**/

u1 Intersection(pu4 pu4V1, pu4 pu4V2, pu4 pu4R, u1 pu4V1Dim, u1 pu4V2Dim, u1 u1CoeffLength)
{

    u1 u1i, u1j;

    // Step 1 : first part of the Zassenhaus matrix
    for(u1i = 0 ; u1i < pu4V1Dim ; u1i++)
    {
        memcpy(pu4R + (2*u1i + 1)*u1CoeffLength, pu4V1 + u1i*u1CoeffLength, u1CoeffLength*4);
        memcpy(pu4R +  2*u1i     *u1CoeffLength, pu4V1 + u1i*u1CoeffLength, u1CoeffLength*4);
    }

    // Step 2 : Second part of the Zassenhaus matrix
    for(u1i = 0 ; u1i < pu4V2Dim ; u1i++)
        memcpy(pu4R + 2* pu4V1Dim*u1CoeffLength + 2*u1i*u1CoeffLength,
               pu4V2 + u1i*u1CoeffLength,
               u1CoeffLength*4);

    // Step 3 : Row echelon form of the Zassenhaus matrix.
    Echelon_matrix(pu4R, pu4V1Dim+ pu4V2Dim, 2*PARAM_M, 2*u1CoeffLength, pu4R);

    // Intersection's base recovery
    for(u1i = 0 ; u1i <  (pu4V2Dim+pu4V1Dim); u1i++)
    {
        if(Is_Zero(pu4R + 2*u1i* MODULO_GF2M_LEN, MODULO_GF2M_LEN) == 0)
        {
            memset(pu4R, 0x00,((pu4R + 2*u1i* u1CoeffLength) - pu4R)*4);

            // Each element of the intersection is copied at the beginning of pu4R.
            for(u1j = 0 ; u1j < (pu4V1Dim+ pu4V2Dim - u1i)  ; u1j++)
            {
                memcpy(pu4R + u1j*u1CoeffLength, pu4R + (2*(u1i + u1j) + 1)*u1CoeffLength, u1CoeffLength*4);
                memset(pu4R + (2*(u1i + u1j) + 1)*u1CoeffLength, 0x00, u1CoeffLength*4);
            }
            return (pu4V1Dim+ pu4V2Dim - u1i );
        }
    }
    return 0;
}


/**
 * \fn u2 PreCalculusRSR(pu4 pu4F, pu4 pu4S, pu4 pu4Workspace, pu1 u1TabDimInterSi)
 * \brief  This function precomputes the intersections between S_i required during RSR algorithm.
 * \param[in]              pu4F : address of the private support F.
 * \param[in]              pu4S : address of the syndrome S.
 * \param[out]     pu4Workspace : address of the result R = pu4V1 inter pu4V2.
 * \param[out]   u1TabDimInterSi : tab of the intersection's dimension.
**/

u2 PreCalculusRSR(pu4 pu4F, pu4 pu4S, pu4 pu4Workspace, pu1 u1TabDimInterSi)
{
  u1 u1i;
  u1 u1SDimMax = PARAM_D * PARAM_R;
  u1 u1CoeffLength = PARAM_M_LEN_WORD;

  pu4 pu4WorkspaceTop = pu4Workspace;
  pu4 pu4Si           = pu4Workspace;
  pu4 pu4Si1          = pu4Workspace + (u1SDimMax+1)*u1CoeffLength;
  pu4 pu4Si2          = pu4Si1       + (u1SDimMax+1)*u1CoeffLength;
  pu4Workspace        = pu4Si2       + (u1SDimMax+1)*u1CoeffLength;

  //S_0, S_1 ans S_2 computations
  S_i_calculus(pu4F, pu4S, pu4Si  , 0, u1SDimMax);
  S_i_calculus(pu4F, pu4S, pu4Si1 , 1, u1SDimMax);
  S_i_calculus(pu4F, pu4S, pu4Si2 , 2, u1SDimMax);


  // Initialisation
  // S_0 inter S_1
  u1TabDimInterSi[0] = Intersection(pu4Si,
                                    pu4Si1,
                                    pu4Workspace,
                                    u1SDimMax,
                                    u1SDimMax,
                                    u1CoeffLength);

  pu4Workspace += u1TabDimInterSi[0]*u1CoeffLength;
  // S_1 inter S_2
  u1TabDimInterSi[1] =Intersection(pu4Si,
                                   pu4Si2,
                                   pu4Workspace,
                                   u1SDimMax,
                                   u1SDimMax,
                                   u1CoeffLength);

  pu4Workspace += u1TabDimInterSi[1]*u1CoeffLength;
  // S_0 inter S_2
  u1TabDimInterSi[2] = Intersection(pu4Si1,
                                    pu4Si2,
                                    pu4Workspace,
                                    u1SDimMax,
                                    u1SDimMax,
                                    u1CoeffLength);


  pu4Workspace += u1TabDimInterSi[2]*u1CoeffLength;

  for(u1i = 0 ; u1i < PARAM_D - 3 ; u1i++)
  {
    //S_i+1 computation
    S_i_calculus(pu4F, pu4S, pu4Si + (u1i % 3)*(u1SDimMax+1)*u1CoeffLength, u1i + 3, u1SDimMax);

    // S_i+1 inter S_i+2
    u1TabDimInterSi[2*(u1i+1) + 1] =Intersection(pu4Si + ((u1i +1) % 3)*(u1SDimMax+1)*u1CoeffLength,
                                             pu4Si + (u1i % 3)*(u1SDimMax+1)*u1CoeffLength,
                                             pu4Workspace,
                                             u1SDimMax,
                                             u1SDimMax,
                                             u1CoeffLength);

    pu4Workspace += u1TabDimInterSi[2*(u1i+1) + 1]*u1CoeffLength;

    // S_i inter S_i+2
    u1TabDimInterSi[2*(u1i+1) + 2] = Intersection(pu4Si + ((u1i+2) % 3)*(u1SDimMax+1)*u1CoeffLength,
                                             pu4Si + (u1i % 3)*(u1SDimMax+1)*u1CoeffLength,
                                             pu4Workspace,
                                             u1SDimMax,
                                             u1SDimMax,
                                             u1CoeffLength);

    pu4Workspace += u1TabDimInterSi[2*(u1i+1) + 2]*u1CoeffLength;

  }

  // Removal of S_i
  memset(pu4WorkspaceTop, 0x00, 3*(u1SDimMax+1)*u1CoeffLength*4);

  return 0;
}



/**
 * \fn u2 Rank_Support_Recover(pu4 E, pu4 pu4F, pu4 pu4Syndrom, pu4 pu4Workspace, u1 u1EDim)
 * \brief  RSR algorithm recovers the private support E from the support of the private key F and the syndrome.
 * \param[out]             pu4E : address of the output recovered support E.
 * \param[in]              pu4F : address of the private support F.
 * \param[in]        pu4Syndrom : address of the syndrome S.
 * \param[in]     pu4Workspace : address of the workspace.
 * \param[in]           u1EDim : dimension of E.
**/

u2 Rank_Support_Recover(pu4 pu4E, pu4 pu4F, pu4 pu4Syndrom, pu4 pu4Workspace)
{

  pu4 pu4Intersections = pu4Workspace + 3*PARAM_M_LEN_WORD*(PARAM_D*PARAM_R+1);
  pu4 pu4Tmp;

  pu4 pu4Si  = pu4Workspace;
  pu4 pu4Si1 = pu4Si + PARAM_M_LEN_WORD* (PARAM_D* PARAM_R+1);

  u1 inter_dim;
  u1 u1i, u1j;
  u1 u1SDim;
  u1 u1SDimComp;
  u1 u1DirectSumDim;
  u1 u1MultiplyDim;
  u1 pu1TabDimIntersection[PARAM_M_LEN_WORD*(PARAM_D-1)];
  u1 E_dim_expected = PARAM_R;
  u2 nbRows_inter;
  // Step 1 : Computation of the syndrome's support and its rank.
  u1SDim = Echelon_matrix(pu4Syndrom,
                 PARAM_N,
                 PARAM_M,
                 PARAM_M_LEN_WORD,
                 pu4Syndrom);


  // Step 2 : Precomputations of the intersections.
  PreCalculusRSR(pu4F, pu4Syndrom, pu4Workspace, pu1TabDimIntersection);

  u1 count_shift = 0;
  // Step 3 : Expansion of the space vector.
  for( u1i = 0 ; u1i < PARAM_D+2 ; u1i=u1i+2)
  {

    pu4Tmp = pu4Workspace;

    // Copy of the 3 intersections required in the multiplication
    for( u1j = 0 ; u1j < 3 ; u1j++ )
    {
      memcpy(pu4Tmp, pu4Intersections + count_shift, pu1TabDimIntersection[u1i+u1j]*PARAM_M_LEN_WORD*4);
      pu4Tmp += pu1TabDimIntersection[u1i+u1j]*PARAM_M_LEN_WORD;
      count_shift += pu1TabDimIntersection[u1i+u1j]*PARAM_M_LEN_WORD;
     }

    count_shift -= pu1TabDimIntersection[u1i+u1j]*PARAM_M_LEN_WORD;
    nbRows_inter =  pu1TabDimIntersection[u1i] + pu1TabDimIntersection[u1i+1]+pu1TabDimIntersection[u1i+2];
    u1DirectSumDim = Echelon_matrix(pu4Workspace,
                   nbRows_inter,
                   PARAM_M,
                   PARAM_M_LEN_WORD,
                   pu4Workspace);

    pu4Tmp = pu4Workspace + u1DirectSumDim*PARAM_M_LEN_WORD;

    // Computation of F*intersections
    Mupliply_GF2mn(pu4F,
                   pu4Workspace,
                   pu4Tmp,
                   PARAM_D,
                   u1DirectSumDim,
                   PARAM_N,
                   PARAM_M_LEN_WORD);

    u1MultiplyDim   = PARAM_D + u1DirectSumDim;

  // Computation of the basis of the resulting space vector and its dimension.
    u1MultiplyDim   = Echelon_matrix(pu4Tmp,
                   u1MultiplyDim,
                   PARAM_M,
                   PARAM_M_LEN_WORD,
                   pu4Tmp);

    // Computation of S + F*intersections
    memcpy(pu4Tmp + u1MultiplyDim*PARAM_M_LEN_WORD,
           pu4Syndrom,
           (PARAM_M_LEN_WORD * PARAM_R * PARAM_D)*4);

    u1SDimComp = Echelon_matrix(pu4Tmp,
                   u1MultiplyDim + u1SDim,
                   PARAM_M,
                   PARAM_M_LEN_WORD,
                   pu4Tmp);

    // Dimension comparison
    if (u1SDimComp <= PARAM_R * PARAM_D)
    {

      memcpy(pu4Syndrom,
             pu4Tmp,
             u1SDimComp*PARAM_M_LEN_WORD*4);
      u1SDim = u1SDimComp;
    }
    memset(pu4Workspace,
           0x00,
           ((pu4Tmp + u1SDimComp*PARAM_M_LEN_WORD) - pu4Workspace)*4);
  }

  memset(pu4Si, 0x00, (PARAM_M_LEN_WORD* PARAM_D* PARAM_R)*4);
  memset(pu4Si1, 0x00, (PARAM_M_LEN_WORD* PARAM_D* PARAM_R)*4);
  memset(pu4Workspace, 0x00, 4*PARAM_M_LEN_WORD*4);
  // Step  4 : Computation of the space vector E (= intersection of S_i)

  // Initialisation with the first intersection S0 inter S1
   pu4Workspace += 2*PARAM_M_LEN_WORD*(PARAM_D*PARAM_R+1);
   memset(pu4Workspace, 0x00, 3*PARAM_D*PARAM_R*PARAM_M_LEN_WORD*4);
   S_i_calculus(pu4F, pu4Syndrom, pu4Si  , 0, PARAM_D*PARAM_R);
   S_i_calculus(pu4F, pu4Syndrom, pu4Si1 , 1, PARAM_D*PARAM_R);

   inter_dim = Intersection(pu4Si,
                   pu4Si1,
                   pu4Workspace,
                   PARAM_D*PARAM_R,
                   PARAM_D*PARAM_R,
                   PARAM_M_LEN_WORD);
  memset(pu4Si, 0x00, PARAM_R*PARAM_D*PARAM_M_LEN_WORD*4);
  memset(pu4Si1, 0x00, PARAM_R*PARAM_D*PARAM_M_LEN_WORD*4);
  memcpy(pu4Si1, pu4Workspace,  inter_dim*PARAM_M_LEN_WORD*4);
  memset(pu4Workspace, 0x00, inter_dim*PARAM_M_LEN_WORD*4);

  //If the dimension of S0 inter S1 < the expected dimension, then S0 inter S1 is intersects with others S_i
  if(inter_dim > E_dim_expected)
  {
    for(u1i = 2 ; u1i < PARAM_D ; u1i++)
  {

     S_i_calculus(pu4F,
                  pu4Syndrom,
                  pu4Si,
                  u1i,
                  PARAM_D*PARAM_R);

     inter_dim = Intersection(pu4Si1,
                   pu4Si,
                   pu4Workspace,
                   inter_dim,
                   PARAM_D*PARAM_R,
                   PARAM_M_LEN_WORD);


     memset(pu4Si, 0x00, PARAM_R*PARAM_D*PARAM_M_LEN_WORD*4);
     memset(pu4Si1, 0x00, PARAM_R*PARAM_D*PARAM_M_LEN_WORD*4);
     memcpy(pu4Si1, pu4Workspace,  inter_dim*PARAM_M_LEN_WORD*4);
     memset(pu4Workspace, 0x00, inter_dim*PARAM_M_LEN_WORD*4);

     if(inter_dim <= E_dim_expected)
     {
       break;
     }

   }
  }

  //Reduced Gaussian elimination applied to recover the support E

  u1 E_dim_recovered = Reduced_Echelon_matrix(pu4Si1,
                    inter_dim,
                    PARAM_M,
                    PARAM_M_LEN_WORD,
                    pu4Si1);
  memcpy(pu4E,pu4Si1,PARAM_R*PARAM_M_LEN_WORD*4);
  if(E_dim_recovered!=E_dim_expected)
  {
    return 1;
  }

  return 0;
}
