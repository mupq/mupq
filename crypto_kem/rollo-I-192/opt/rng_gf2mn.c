#include <stdio.h>
#include <string.h>

#include "typedef.h"
#include "rng_gf2mn.h"
#include "matrix.h"
#include "api.h"
#include "gf2m.h"
#include "utils.h"
#include "randombytes.h"


/**********************************************************************
*
*       Generation of elements in GF(2^m)^n
*
**********************************************************************/


/**
 * \fn u1 Random_GF2_m(pu1 v, u1 size)
 * \brief Generation of element in GF(2^m)^size
 *
 * \param[out]     v : address of the destination pointer.
 * \param[in]   size : number of coefficients in GF(2^m) to generate.
 **/

u1 Random_GF2_m(pu1 v, u1 size)
{
  u1 mask = (1 << PARAM_M % 8) - 1;
  u1 random[PARAM_M_LEN_BYTES*size];
  randombytes(random, PARAM_M_LEN_BYTES*size);
  for(int i=0;i<size;i++)
  {
    random[(i+1) * PARAM_M_LEN_BYTES - 1] &= mask;
    memcpy(v + PARAM_M_LEN_BYTES_PADDING * i, random + PARAM_M_LEN_BYTES*i, PARAM_M_LEN_BYTES);
  }
  return 0;
}

/**
 * \fn u1 Create_support_using_rng(pu4 pu4Support, u1 u1Rank)
 * \brief Generation of support with a given rank
 *
 * \param[out]     pu4Support : address of the support.
 * \param[in]          u1Rank : Dimmension of the support (= number of independent elements)
 **/

u1 Create_support_using_rng(pu4 pu4Support, u1 u1Rank)
{
  u1 u1RankSupport = -1;
  pu4 Copy_Support = pu4Support + u1Rank*PARAM_M_LEN_WORD;
  pu1 Support = (pu1) pu4Support;


  while ( u1RankSupport != u1Rank)
  {
    Random_GF2_m(Support, u1Rank);
    memcpy(Copy_Support,(pu4)Support,u1Rank*PARAM_M_LEN_WORD*4);
    u1RankSupport = Echelon_matrix(Copy_Support, u1Rank, PARAM_M_LEN_WORD, PARAM_M_LEN_WORD, Copy_Support);
    memset(Copy_Support,0x00, u1Rank*PARAM_M_LEN_WORD*4);
  }
  return 0;
}



/**********************************************************************
*      Cette fonction va créer un élément de rang et de taille voulue.
*
* pu1Random : Adresse où mettre l'élément.
* u1Rank    : Rang de l'élément (dimension de l'espace vectoriel associé).
* u1Length  : Nombre de coefficients.
*
**********************************************************************/
/**
 * \fn u1 Create_gf2mn_using_support(pu4 pu4Element, pu4 pu4Support, u1 u1Rank)
 * \brief Generation of support with a given rank
 *
 * \param[out]     pu4Support : address of the support.
 * \param[in]          u1Rank : Dimmension of the support (= number of independent elements)
 **/

u1 Create_gf2mn_using_support(pu4 pu4Element, pu4 pu4Support, u1 u1Rank)
{
  u1 random_size = 2 * u1Rank;
  u1 random[PARAM_N];
  randombytes(random, random_size);

  u1 i = 0;
  u1 j = 0;
  u1 position;
  while(i != u1Rank) {
    position = random[j];
    if(!Is_Zero(pu4Element + (position % PARAM_N)*PARAM_M_LEN_WORD,PARAM_M_LEN_WORD))
    {
      memcpy(pu4Element + (position % PARAM_N)*PARAM_M_LEN_WORD, pu4Support+i*PARAM_M_LEN_WORD, PARAM_M_LEN_WORD*4);
      i++;
    }
    j++;
    if((j % random_size == 0) && (i != u1Rank))
    {
      randombytes(random, random_size);
      j=0;

    }

  }
  randombytes(random, PARAM_N);

  int k =0;
  int l =0;
  for(i=0;i<PARAM_N;i++)
  {
    if(!Is_Zero(pu4Element+i*PARAM_M_LEN_WORD,PARAM_M_LEN_WORD))
    {
      for(j=0;j<u1Rank;j++)
      {
        if(*(random+k)&0x01)
        {
          Add_GF2m(pu4Element+i*PARAM_M_LEN_WORD, pu4Support+j*PARAM_M_LEN_WORD ,pu4Element+i*PARAM_M_LEN_WORD, PARAM_M_LEN_WORD);
        }
        *(random+k) = *(random+k) >> 1;
        l++;
        if(l == 8) {
          l = 0;
          k++;
        }
      }
    }

  }
  return 0;

}
