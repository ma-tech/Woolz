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
* \ingroup      AlgMatrix
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
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the sum of two matricies and returns
*		the result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} + \mathbf{C}
		\f]
*		The dimensions of the all three matricies must be
*		nR, nC. It is safe to supply the same matrix as any
*		combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
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
			     size_t nR, size_t nC)
{
  size_t        id0,
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
* \return       void
* \ingroup      AlgMatrix
* \brief        Subtracts on matrix from another and returns the
*		result in a third supplied matrix:
*		\f[
		  \mathbf{A} = \mathbf{B} - \mathbf{C}
		\f]
*		The dimensions of the all three
*		matricies must be nR, nC. It is safe to supply
*		the same matrix as any combination of aM, bM and cM.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
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
			     size_t nR, size_t nC)
{
  size_t	id0,
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
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the product of two matricies and returns
*		the result in a third supplied matrix:
*               \f[
                  \mathbf{A} = \mathbf{B} \mathbf{C}
		\f]
*		The dimensions of the result matrix (aM) must be cR, bC.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
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
			     size_t bR, size_t bC, size_t cC)
{
  size_t	id0,
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
* \return       The trace of the matrix.
* \ingroup      AlgMatrix
* \brief        Computes the trace of the given matrix.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix.
* \param	nRC			Number of rows and columns in
*					the (square) matrix aM.
*/
double		AlgMatrixTrace(double **aM, size_t nRC)
{
  size_t	id0;
  double	trace = 0.0;

  for(id0 = 0; id0 < nRC; ++id0)
  {
    trace += aM[id0][id0];
  }
  return(trace);
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Computes the transpose of the given matrix:
*		\f[
		\mathbf{A} = \mathbf{B^T}
		\f]
*		aM = transpose(bM). The dimensions of the result
*		matrix (aM) must be aR == bC, aC == bR.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to transpose,
*					\f$\mathbf{B}\f$.
* \param	bR			Number of rows in matrix bM.
* \param	bC			Number of columns in matrix bM.
*/
void		AlgMatrixTranspose(double **aM, double **bM,
				   size_t bR, size_t bC)
{
  size_t	id0,
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
* \return       void
* \ingroup      AlgMatrix
* \brief        Copies the values of the matrix bM to the result
*		matric aM:
*		\f[
		\mathbf{A} = \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Matrix to copy, \f$\mathbf{B}\f$.
* \param	nR			Number of rows in matricies.
* \param	nC			Number of columns in matricies.
*/
void		AlgMatrixCopy(double **aM, double **bM,
			      size_t nR, size_t nC)
{
  size_t	id0,
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
* \return       void
* \ingroup      AlgMatrix
* \brief        Multiplies the given matrix by the given scalar:
*		\f[
		\mathbf{A} = s \mathbf{B}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result,
*					\f$\mathbf{A}\f$.
* \param        bM 			Given matrix to scale,
*					\f$\mathbf{B}\f$.
* \param	sv			Scalar value, \f$s\f$.
* \param	nR			Number of rows in matrix aM.
* \param	nC			Number of columns in matrix aM.
*/
void		AlgMatrixScale(double **aM, double **bM, double sv,
			       size_t nR, size_t nC)
{
  size_t	id0,
  		id1;
  double 	*aRowM,
  		*bRowM;

  for(id0 = 0; id0 < nR; ++id0)
  {
    aRowM = aM[id0];
    bRowM = bM[id0];
    for(id1 = 0; id1 < nC; ++id1)
    {
      *aRowM++ = *bRowM++ * sv;
    }
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Multiplies the a matrix by a scalar and then adds
*		another matrix:
*		\f[
		\mathbf{A} = \mathbf{B} + s \mathbf{C}.
		\f]
*		All the matrices must have the same dimensions.
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
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
				  double sv, size_t nR, size_t nC)
{
  size_t	id0,
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
* \return       void
* \ingroup      AlgMatrix
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
void		AlgMatrixScalar(double **aM, double sv, size_t nRC)
{
  size_t	id0;

  (void )memset(*aM, 0, nRC * nRC);
  for(id0 = 0; id0 < nRC; ++id0)
  {
    aM[id0][id0] = sv;
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Sets the elements of the given matrix to zero.
* \note		For efficiency the given parameters are not checked.
*		\f[
		\mathbf{A} = \mathbf{0}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by AlcDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param        aM 			Supplied matrix for result.
* \param	nR			Number of rows in matrix.
* \param	nC			Number of columns in matrix.
*/
void		AlgMatrixZero(double **aM, size_t nR, size_t nC)
{
  size_t	id0,
  		nRC;
  double	*aRowM;

  (void )memset(*aM, 0, nR * nC);
}

/*!
* \return
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B} \mathbf{c}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bType			Type of matrix \f$mathbf{B}\f$.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	nR			The number of rows in \f$mathbf{B}\f$.
* \param	nC			The number of columns in
* 					\f$mathbf{c}\f$.
*/
void 		AlgMatrixVectorMul(double *aV,
				   AlgMatrixType bType, double **bM,
				   double *cV, size_t nR, size_t nC)
{
#ifdef _OPENMP
  int		id0,
		id1,
		oNR;
#else
  size_t	id0,
  		id1;
#endif
  double	tD0;
  double	*bRow,
  		*cCol;

#ifdef _OPENMP
  oNR = nR;
#endif
  switch(bType)
  {
    case ALG_MATRIX_RECT:
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(id0,id1,tD0,cCol,bRow)
      for(id0 = 0; id0 < oNR; ++id0)
#else
      for(id0 = 0; id0 < nR; ++id0)
#endif
      {
	tD0 = 0.0;
	cCol = cV;
	bRow = bM[id0];
	for(id1 = 0; id1 < nC; ++id1)
	{
	  tD0 += bRow[id1] * cCol[id1];
	}
	aV[id0] = tD0;
      }
      break;
    case ALG_MATRIX_SYM:
#ifdef _OPENMP
      #pragma omp parallel for default(shared) private(id0,id1,tD0,cCol)
      for(id0 = 0; id0 < oNR; ++id0)
#else
      for(id0 = 0; id0 < nR; ++id0)
#endif
      {
        tD0 = 0.0;
	cCol = cV;
	for(id1 = 0; id1 <= id0; ++id1)
	{
	  tD0 += bM[id0][id1] * cCol[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  tD0 += bM[id1][id0] * cCol[id1];
	}
      }
      break;
  }
}
