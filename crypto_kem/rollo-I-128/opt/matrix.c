#include <string.h>

#include "typedef.h"
#include "api.h"
#include "gf2m.h"
#include "utils.h"
#include "matrix.h"


/***********************************
*
*       Matrix operations
*
************************************/


/**
* \fn     u1 SwapLines(pu4 pu4A, pu4 pu4B, u1 u1Length)
* \brief  Swaps the address of two elements.
*
* \param[in]  pu1A      : Address of the first element to swap.
* \param[in]  pu1B      : Address of the second element.
* \param[in]  u1Length  : Length of elements in 32-bit words.
*
**/

u1 SwapLines(pu4 pu4A, pu4 pu4B, u1 u1Length)
{

 u4 tmp;
 u4   i;

 for(i = 0 ; i < u1Length ; i++)
 {
   tmp = *((pu4)pu4A + i);
   *((pu4) pu4A + i) = *((pu4) pu4B + i);
   *((pu4) pu4B + i) = tmp;
 }

  return 0;
}


/**
* \fn    u1 Find_pivot(pu4 pu4Matrix, u1 nu1NbRows, u1 nu1Line, u1 u1Pivot, u1 u1ElementLength)
* \brief Returns the line containing the pivot.
*
* \param[in]  pu1Matrix       : Address of the matrix.
* \param[in]  nu1NbRows       : Number of rows in the matrix.
* \param[in]  nu1Line         : Line from which the pivot must be searched.
* \param[in]  u1Pivot         : Index of the column of the searched pivot.
* \param[in]  u1ElementLength : Length of an element ( <=> number of columns).
*
**/


u1 Find_pivot(pu4 pu4Matrix, u1 nu1NbRows, u1 nu1Line, u1 u1Pivot, u1 u1ElementLength)
{

  u2 u1i;
  u4 u4Bit;
  u4 u4MaskPivot, u4NbWord;
  u4 decalage;

  // Step 1 : Computation of the mask in order to find the pivot.
  u4NbWord       = (u1Pivot / 32);
  u4MaskPivot    = 1 << (u1Pivot%32);

  // Step 2 : For each row of the matrix, we apply the mask to recover the leading coefficient in each row.
  // The first non zero bit is the pivot.
  for(u1i = nu1Line ; u1i < nu1NbRows ; u1i++)
  {

    decalage    = u1i*u1ElementLength + u4NbWord;
    u4Bit       = *(pu4Matrix + decalage) & u4MaskPivot;
    if ( u4Bit != 0x00)
    {
      return u1i;
    }
  }
  return nu1NbRows;
}



/**
* \fn     u1 Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace)
* \brief  This function puts a matrix under its row echelon form.
*
*
* \param[in] pu1Matrix         : Address of the matrix.
* \param[in] nu1NbRows         : Number of rows in the matrix.
* \param[in] u1NbCols          : Number of columns in the matrix.
* \param[in] u1ElementLength   : Length of an element ( <=> number of columns).
* \param[in] pu1Workspace      : Address of workspace (same size of the matrix).
*                   -> If pu1Matrix == pu1Workspace then pu1Matrix will be remove.
* \param[out] rank_matrix      : rank of pu1Matrix.
*
**/


u1 Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace)
{
  if ( pu4Matrix != pu4Workspace)
  {
	 memcpy(pu4Workspace, pu4Matrix, u1NbRows*u1NbCols*4);
  }

  pu4Matrix = pu4Workspace;
  pu4 pu4BasePivot;
  pu4 pu4ProcessedRow;
  int u1i, u1j, u1Index_pivot;
  u4 u4Mask;
  u4 decalage;
  u2 decalage_bit;
  u4 u4NbWord;
  u1 u1NbLine;
  u1 rank_matrix = 0;


  decalage_bit = 0;


  for( u1i = PARAM_M; u1i >= 0; u1i--)
  {
    u4NbWord = u1i / 32;
    u1NbLine = (PARAM_M - u1i) - decalage_bit;

    // Step 1 : Search the pivot for each column of the matrix.
    u1Index_pivot = Find_pivot(pu4Matrix, u1NbRows, u1NbLine, u1i, u1ElementLength);


    // If no pivot has been found, the next column is processed.
    if( u1Index_pivot == u1NbRows)
    {
      decalage_bit += 1;
      continue;
    }
    else if(u1Index_pivot != u1NbLine)
    {
      SwapLines(pu4Matrix + (u1NbLine)*u1ElementLength,
                pu4Matrix + u1Index_pivot*u1ElementLength,
                u1ElementLength);
    }


    rank_matrix ++;


    // Step 2 : Swap between the row containing the pivot and the row of pivot index.
    pu4BasePivot = (pu4Matrix + (u1NbLine)*u1ElementLength);

    // Compute the mask in order to recover the leading coefficient in each row.
    u4Mask = 1 << (u1i%32);

    for( u1j = u1NbLine+1 ; u1j < u1NbRows ; u1j++)
    {
      decalage = u1j*u1ElementLength + u4NbWord;
       pu4ProcessedRow = pu4Matrix + decalage;

      // If the bit is 1, we add the pivot row with the processed row.
      if( (*(pu4ProcessedRow) & u4Mask))
      {
        Add_GF2m(pu4BasePivot, (pu4ProcessedRow - u4NbWord ), (pu4ProcessedRow - u4NbWord ), u1ElementLength);
      }

    }

  }
  return rank_matrix;
}


/**
* \fn     u1 Reduced_Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace)
* \brief  This function computes the reduced row echelon form of a matrix and returns its rank.
*
*
* \param[in] pu1Matrix         : Address of the matrix.
* \param[in] nu1NbRows         : Number of rows in the matrix.
* \param[in] u1NbCols          : Number of columns in the matrix.
* \param[in] u1ElementLength   : Length of an element ( <=> number of columns).
* \param[in] pu1Workspace      : Address of workspace (same size of the matrix).
*                   -> If pu1Matrix == pu1Workspace then pu1Matrix will be remove.
* \param[out] rank_matrix      : rank of pu1Matrix.
*
**/

u1 Reduced_Echelon_matrix(pu4 pu4Matrix, u1 u1NbRows, u1 u1NbCols,u1 u1ElementLength, pu4 pu4Workspace)
{
  if ( pu4Matrix != pu4Workspace)
  {
	 memcpy(pu4Workspace, pu4Matrix, u1NbRows*u1NbCols*4);
  }

  pu4Matrix = pu4Workspace;
  pu4 pu4BasePivot;
  pu4 pu4ProcessedRow;
  int u1i, u1j, u1Index_pivot;
  u4 u4Mask;
  u4 decalage;
  u2 decalage_bit;
  u4 u4NbWord;
  u1 u1NbLine;
  u1 rank_matrix =0;

  decalage_bit = 0;


  for( u1i = PARAM_M ; u1i >= 0; u1i--)
  {
    u4NbWord = u1i / 32;
    u1NbLine = (PARAM_M - u1i) - decalage_bit;

    // Step 1 : Search the pivot for each column of the matrix.
    u1Index_pivot = Find_pivot(pu4Matrix, u1NbRows, u1NbLine, u1i, u1ElementLength);


    // If no pivot has been found, the next column is processed.
    if( u1Index_pivot == u1NbRows)
    {
      decalage_bit += 1;
      continue;
    }

    // Step 2 : Swap between the row containing the pivot and the row of pivot index.
    else if(u1Index_pivot != u1NbLine)
    {
      SwapLines(pu4Matrix + (u1NbLine)*u1ElementLength,
                pu4Matrix + u1Index_pivot*u1ElementLength,
                u1ElementLength);
    }

    rank_matrix ++;
    pu4BasePivot = (pu4Matrix + (u1NbLine)*u1ElementLength);

    // Compute the mask in order to recover the leading coefficient in each row.
    u4Mask = 1 << (u1i%32);

    for( u1j = 0 ; u1j < u1NbRows ; u1j++)
    {
      decalage = u1j*u1ElementLength + u4NbWord;
      pu4ProcessedRow = pu4Matrix + decalage;

      // If the bit is 1, we add the pivot row with the processed row.
      if( (*(pu4ProcessedRow) & u4Mask)&& (u1j!=u1NbLine))
      {
        Add_GF2m(pu4BasePivot, (pu4ProcessedRow - u4NbWord ), (pu4ProcessedRow - u4NbWord ), u1ElementLength);
      }

    }

  }
  return rank_matrix;
}
