#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixMath.c
* \author       Bill Hill
* \date         June 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for basic arithmatic with matricies.
* \bug          None known.
* \note
* Maintenance log with most recent changes at top of list.
*/

/*!
* \ingroup      Alg
* \defgroup      AlgMatrix
* @{
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \return       <void>
* \brief        Computes the sum of two matricies and returns
*		the result in a third supplied matrix,
*		aM = bM + cM. The dimensions of the all three
*		matricies must be nR, nC. It is safe to supply
*		the same matrix as any combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result.
* \param        bM 			First matrix in the sum.
* \param        cM			Second matrix in the sum.
* \param	nR			Number of rows in matricies.
* \param	nC			Number of columns in matricies.
*/
void		AlgMatrixAdd(double **aM, double **bM, double **cM,
			     int nR, int nC)
{
  int		id0,
  		id1;
  double	*aRowM,
  		*bRowM,
		*cRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    cRowM = cM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ = *bRowM++ + *cRowM++;
    }
  }
}

/*!
* \return       <void>
* \brief        Subtracts on matrix from another and returns the
*		result in a third supplied matrix,
*		aM = bM - cM. The dimensions of the all three
*		matricies must be nR, nC. It is safe to supply
*		the same matrix as any combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result.
* \param        bM 			First matrix in the subtraction.
* \param        cM			Second matrix in the subtraction.
* \param	nR			Number of rows in matricies.
* \param	nC			Number of columns in matricies.
*/
void		AlgMatrixSub(double **aM, double **bM, double **cM,
			     int nR, int nC)
{
  int		id0,
  		id1;
  double	*aRowM,
  		*bRowM,
		*cRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    cRowM = cM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ = *bRowM++ - *cRowM++;
    }
  }
}

/*!
* \return       <void>
* \brief        Computes the product of two matricies and returns
*		the result in a third supplied matrix,
*		aM = bM x cM. The dimensions of the result matrix
*		(aM) must be cR, bC.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result.
* \param        bM 			First matrix in the product.
* \param        cM			Second matrix in the product.
* \param	bR			Number of rows in matrix bM.
* \param	bC			Number of columns in matrix bM.
* \param	cC			Number of columns in matrix cM.
*/
void		AlgMatrixMul(double **aM, double **bM, double **cM,
			     int bR, int bC, int cC)
{
  int		id0,
  		id1,
		id2;
  double	tD0;
  double	*bRowM;

  for(id0 = 0; id0 < bR; ++id0)
  {
    bRowM = bM[id0];
    for(id1 = 0; id1 < cC; ++id1)
    {
      tD0 = 0.0;
      for(id2 = 0; id2 < bC; ++id2)
      {
	tD0 += bRowM[id2] * cM[id2][id1];
      }
      aM[id0][id1] = tD0;
    }
  }
}

/*!
* \return       			The trace of the matrix.
* \brief        Computes the trace of the given matrix.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix.
* \param	nRC			Number of rows and columns in
*					the (square) matrix aM.
*/
double		AlgMatrixTrace(double **aM, int nRC)
{
  int		id0;
  double	trace = 0.0;

  for(id0 = 0; id0 < nRC; ++id0)
  {
    trace += aM[id0][id0];
  }
  return(trace);
}

/*!
* \return       <void>
* \brief        Computes the transpose of the given matrix.
*		aM = transpose(bM). The dimensions of the result
*		matrix (aM) must be aR == bC, aC == bR.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result.
* \param        bM 			Matrix to transpose.
* \param	bR			Number of rows in matrix bM.
* \param	bC			Number of columns in matrix bM.
*/
void		AlgMatrixTranspose(double **aM, double **bM, int bR, int bC)
{
  int		id0,
  		id1;

  for(id0 = 0; id0 < bR; ++id0)
  {
    for(id1 = 0; id1 < bC; ++id1)
    {
      aM[id0][id1] = bM[id1][id0];
    }
  }
}

/*!
* \return       <void>
* \brief        Copies the values of the matrix bM to the result
*		matric aM.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result.
* \param        bM 			Matrix to transpose.
* \param	nR			Number of rows in matricies.
* \param	nC			Number of columns in matricies.
*/
void		AlgMatrixCopy(double **aM, double **bM, int nR, int nC)
{
  int		id0,
  		id1;
  double	*aRowM,
  		*bRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ = *bRowM++;
    }
  }
}


/*!
* @}
*/
