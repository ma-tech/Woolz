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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

/*!
* \ingroup   	Alg
* \defgroup     AlgMatrix
* @{
*/

/*!
* \return       <void>
* \brief        Computes the sum of two matricies and returns
*		the result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} + \mathbf{C}
		\f]
*		The dimensions of the all three matricies must be
*		nR, nC. It is safe to supply the same matrix as any
*		combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			First matrix in the sum,
*					\f$\mathbf{B}\f$.
* \param        cM			Second matrix in the sum,
*					\f$\mathbf{C}\f$.
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
*		result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} - \mathbf{C}
		\f]
*		The dimensions of the all three
*		matricies must be nR, nC. It is safe to supply
*		the same matrix as any combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			First matrix in the subtraction,
					\f$\mathbf{B}\f$.
* \param        cM			Second matrix in the subtraction,
					\f$\mathbf{C}\f$.
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
*		the result in a third supplied matrix:
*               \f[
                  \mathbf{A} = \mathbf{B} \mathbf{C}
		\f]
*		The dimensions of the result matrix (aM) must be cR, bC.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$
* \param        bM 			First matrix in the product,
*					\f$\mathbf{B}\f$
* \param        cM			Second matrix in the product,
*					\f$\mathbf{C}\f$
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
* \brief        Computes the transpose of the given matrix:
*		\f[
		\mathbf{A} = \mathbf{B^T}
		\f]
*		aM = transpose(bM). The dimensions of the result
*		matrix (aM) must be aR == bC, aC == bR.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to transpose,
*					\f$\mathbf{B}\f$.
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
*		matric aM:
*		\f[
		\mathbf{A} = \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to copy, \f$\mathbf{B}\f$.
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
* \return       <void>
* \brief        Multiplies the given matrix by the given scalar:
*		\f[
		\mathbf{A} = s \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Given matrix to scale,
*					\f$\mathbf{B}\f$.
* \param	sv			Scalar value, \f$s\f$.
* \param	nR			Number of rows in matrix aM.
* \param	nC			Number of columns in matrix aM.
*/
void		AlgMatrixScale(double **aM, double **bM, double sv,
			       int nR, int nC)
{
  int		id0,
  		id1;
  double 	*aRowM,
  		*bRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ *= *bRowM++ * sv;
    }
  }
}

/*!
* \return       <void>
* \brief        Multiplies the a matrix by a scalar and then adds
*		another matrix:
*		\f[
		\mathbf{A} = \mathbf{B} + s \mathbf{C}.
		\f]
*		All the matrices must have the same dimensions.
* \note		For efficiency the given parameters are not checked.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Given matrix to scale,
*					\f$\mathbf{B}\f$.
* \param	cM			Matrix too add, \f$\mathbf{C}\f$.
* \param	sv			Scalar value, \f$s\f$.
* \param	nR			Number of rows in each matrix.
* \param	nC			Number of columns in each matrix.
*/
void		AlgMatrixScaleAdd(double **aM, double **bM, double **cM,
				  double sv, int nR, int nC)
{
  int		id0,
  		id1;
  double 	*aRowM,
  		*bRowM,
		*cRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    cRowM = cM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ = *bRowM++ + (*cRowM++ * sv);
    }
  }
}

/*!
* \return       <void>
* \brief        Sets the elements of the given square matrix so that it
*		is a scalar matrix:
*		\f[
                \mathbf{A} = s \mathbf{I}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		This function assumes that the matrix has been allocated
*		by AlcDouble2Malloc().
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param	sv			Scalar value, \f$s\f$.
* \param	nRC			Number of rows and columns in matrix.
*/
void		AlgMatrixScalar(double **aM, double sv, int nRC)
{
  int		id0;

  (void )memset(*aM, 0, nRC * nRC);
  for(id0 = 0; id0 < nRC; ++id0)
  {
    aM[id0][id0] = sv;
  }
}

/*!
* \return       <void>
* \brief        Sets the elements of the given matrix to zero.
* \note		For efficiency the given parameters are not checked.
*		\f[
		\mathbf{A} = \mathbf{0}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by AlcDouble2Malloc().
* \param        aM 			Supplied matrix for result.
* \param	nR			Number of rows in matrix.
* \param	nC			Number of columns in matrix.
*/
void		AlgMatrixZero(double **aM, int nR, int nC)
{
  int		id0,
  		nRC;
  double	*aRowM;

  (void )memset(*aM, 0, nR * nC);
}

/*!
* \return
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B} \mathbf{c}
		\f]
*		This function assumes that the matrix has been allocated
*		by AlcDouble2Malloc().
* \param	aV			Supplied vector for result.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	nR			The number of rows in \f$mathbf{B}\f$.
* \param	nC			The number of columns in
* 					\f$mathbf{c}\f$.
*/
void 		AlgMatrixVectorMul(double *aV, double **bM, double *cV,
				   int nR, int nC)
{
  int		id0,
  		id1;
  double	tD0;
  double	*bRow,
  		*cCol;

  for(id0 = 0; id0 < nR; ++id0)
  {
    tD0 = 0.0;
    cCol = cV;
    bRow = bM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      tD0 += *bRow++ * *cCol++;
    }
    aV[id0] = tD0;
  }
}

/*!
* @}
*/
