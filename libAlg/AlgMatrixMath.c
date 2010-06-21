#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgMatrixMath_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgMatrixMath.c
* \author       Bill Hill
* \date         June 2001
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief        Functions for basic arithmatic with matricies.
* \ingroup      AlgMatrix
* \todo		-
* \bug          None known.
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
void		AlgMatrixCopy(double **aM, double **bM, size_t nR, size_t nC)
{
  size_t	id0;

  for(id0 = 0; id0 < nR; ++id0)
  {
    AlgMatrixVectorCopy(aM[id0], bM[id0], nC);
  }
}

/*!
* \return       void
* \ingroup      AlgMatrix
* \brief        Copies the values of the vector bV to the result
*		vector aV:
*		\f[
		\mathbf{a} = \mathbf{b}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aV 			Supplied vector for result,
*					\f$\mathbf{a}\f$.
* \param        bV 			Vector to copy, \f$\mathbf{b}\f$.
* \param	nV			Number of entries to copy.
*/
void		AlgMatrixVectorCopy(double *aV, double *bV, size_t nV)
{
  size_t	id0;

  for(id0 = 0; id0 < nV; ++id0)
  {
    aV[id0] = bV[id0];
  }
}

/*!
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
  size_t	id0;

  for(id0 = 0; id0 < nR; ++id0)
  {
    AlgMatrixVectorScale(aM[id0], bM[id0], sv, nC);
  }
}

/*!
* \ingroup      AlgMatrix
* \brief        Multiplies the given vector elements by the given scalar:
*		\f[
		\mathbf{a} = s \mathbf{b}
		\f]
* \note		For efficiency the given parameters are not checked.
* \note		Matrix size is limited only by address space.
* \param        aV 			Supplied vector for result,
*					\f$\mathbf{a}\f$.
* \param        bV 			Given vector to scale,
*					\f$\mathbf{b}\f$.
* \param	sv			Scalar value, \f$s\f$.
* \param	nV			Number of vector entries.
*/
void		AlgMatrixVectorScale(double *aV, double *bV, double sv,
			            size_t nV)
{
  size_t	id0;

  for(id0 = 0; id0 < nV; ++id0)
  {
    aV[id0] = sv * bV[id0];
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

  (void )memset(*aM, 0, sizeof(double) * nRC * nRC);
  for(id0 = 0; id0 < nRC; ++id0)
  {
    aM[id0][id0] = sv;
  }
}

/*!
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
  (void )memset(*aM, 0, sizeof(double) * nR * nC);
}

/*!
* \ingroup      AlgMatrix
* \brief        Sets the elements of the given vector to zero.
* \note		For efficiency the given parameters are not checked.
*		\f[
		\mathbf{a} = \mathbf{0}
		\f]
* \note		Vector size is limited only by address space.
* \param        aV 			Supplied Vector for result.
* \param	nV			Number of vector elements.
*/
void		AlgMatrixVectorZero(double *aV, size_t nV)
{
  (void )memset(aV, 0, sizeof(double) * nV);
}

/*!
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
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixVectorMul(double *aV,
				   AlgMatrixType bType, double **bM,
				   double *cV, size_t nR, size_t nC)
{
  size_t	id0,
  		id1;
  double	tD0;
  double	*bRow;

  switch(bType)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nR; ++id0)
      {
	tD0 = 0.0;
	bRow = bM[id0];
	for(id1 = 0; id1 < nC; ++id1)
	{
	  tD0 += bRow[id1] * cV[id1];
	}
	aV[id0] = tD0;
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nR; ++id0)
      {
        tD0 = 0.0;
	for(id1 = 0; id1 <= id0; ++id1)
	{
	  tD0 += bM[id0][id1] * cV[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  tD0 += bM[id1][id0] * cV[id1];
	}
	aV[id0] = tD0;
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$ and adds the vector \f$\mathbf{d}\f$:
*		\f[
		\mathbf{a} = \mathbf{B} \mathbf{c} + \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bType			Type of matrix \f$mathbf{B}\f$.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* \param	nR			The number of rows in \f$mathbf{B}\f$.
* \param	nC			The number of columns in
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixVectorMulAdd(double *aV,
				   AlgMatrixType bType, double **bM,
				   double *cV, double *dV,
				   size_t nR, size_t nC)
{
  size_t	id0,
  		id1;
  double	tD0;
  double	*bRow;

  switch(bType)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nR; ++id0)
      {
	tD0 = 0.0;
	bRow = bM[id0];
	for(id1 = 0; id1 < nC; ++id1)
	{
	  tD0 += bRow[id1] * cV[id1];
	}
	aV[id0] = tD0 + dV[id0];
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nR; ++id0)
      {
        tD0 = 0.0;
	for(id1 = 0; id1 <= id0; ++id1)
	{
	  tD0 += bM[id0][id1] * cV[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  tD0 += bM[id1][id0] * cV[id1];
	}
	aV[id0] = tD0 + dV[id0];
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the matrix \f$\mathbf{B}\f$ by the vector
*		\f$\mathbf{c}\f$ and adds the vector \f$\mathbf{d}\f$
*		using the given weights \f$s\f$ and \f$t\f$:
*		\f[
		\mathbf{a} = s \mathbf{B} \mathbf{c} + t \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bType			Type of matrix \f$mathbf{B}\f$.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* \param	nR			The number of rows in \f$mathbf{B}\f$.
* \param	nC			The number of columns in
* 					\f$mathbf{B}\f$.
* \param	s			First weighting scalar \f$s\f$.
* \param	t			Second weighting scalar \f$t\f$.
*/
void 		AlgMatrixVectorMulWAdd(double *aV,
				   AlgMatrixType bType, double **bM,
				   double *cV, double *dV,
				   size_t nR, size_t nC,
				   double s, double t)
{
  size_t	id0,
  		id1;
  double	tD0;
  double	*bRow;

  switch(bType)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nR; ++id0)
      {
	tD0 = 0.0;
	bRow = bM[id0];
	for(id1 = 0; id1 < nC; ++id1)
	{
	  tD0 += bRow[id1] * cV[id1];
	}
	aV[id0] = (s * tD0) + (t * dV[id0]);
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nR; ++id0)
      {
        tD0 = 0.0;
	for(id1 = 0; id1 <= id0; ++id1)
	{
	  tD0 += bM[id0][id1] * cV[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  tD0 += bM[id1][id0] * cV[id1];
	}
	aV[id0] = (s * tD0) + (t * dV[id0]);
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the transpose of matrix \f$\mathbf{B}\f$ by the
* 		vector \f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B}^T \mathbf{c}
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
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixTVectorMul(double *aV,
				    AlgMatrixType bType, double **bM,
				    double *cV, size_t nR, size_t nC)
{
  size_t	id0,
  		id1;
  double	tD0;

  switch(bType)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nR; ++id0)
      {
	aV[id0] = 0.0;
	for(id1 = 0; id1 < nC; ++id1)
	{
	  aV[id0] += bM[id1][id0] * cV[id1];
	}
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nR; ++id0)
      {
        aV[id0] = bM[id0][0] * cV[0];
	for(id1 = 1; id1 <= id0; ++id1)
	{
	  aV[id0] += bM[id0][id1] * cV[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  aV[id0] += bM[id1][id0] * cV[id1];
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      AlgMatrix
* \brief	Multiplies the transpose of matrix \f$\mathbf{B}\f$ by the
* 		vector \f$\mathbf{c}\f$:
*		\f[
		\mathbf{a} = \mathbf{B}^T \mathbf{c} + \mathbf{d}
		\f]
* \note		This function assumes that the matrix has been allocated
*		by either AlcDouble2Malloc() or AlcSymDouble2Malloc().
* \note		Matrix size is limited only by address space.
* \param	aV			Supplied vector for result.
* \param	bType			Type of matrix \f$mathbf{B}\f$.
* \param	bM			Matrix \f$mathbf{B}\f$.
* \param	cV			Vector \f$\mathbf{c}\f$.
* \param	dV			Vector \f$\mathbf{d}\f$.
* \param	nR			The number of rows in \f$mathbf{B}\f$.
* \param	nC			The number of columns in
* 					\f$mathbf{B}\f$.
*/
void 		AlgMatrixTVectorMulAdd(double *aV,
				       AlgMatrixType bType, double **bM,
				       double *cV, double *dV,
				       size_t nR, size_t nC)
{
  size_t	id0,
  		id1;

  switch(bType)
  {
    case ALG_MATRIX_RECT:
      for(id0 = 0; id0 < nC; ++id0)
      {
	aV[id0] = dV[id0];
        for(id1 = 0; id1 < nR; ++id1)
	{
	  aV[id1] += bM[id0][id1] * cV[id0];
	}
      }
      break;
    case ALG_MATRIX_SYM:
      for(id0 = 0; id0 < nR; ++id0)
      {
        aV[id0] = dV[id0];
	for(id1 = 0; id1 <= id0; ++id1)
	{
	  aV[id0] += bM[id0][id1] * cV[id1];
	}
	for(id1 = id0 + 1; id1 < nR; ++id1)
	{
	  aV[id0] += bM[id1][id0] * cV[id1];
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return	Euclidean or L2 norm of the vector.
* \ingroup	AlgMatrix
* \brief	Computes the Euclidean or L2 norm of the given vector.
* \param	aV			Given vector.
* \param	nV			Number of entries in the vector.
*/
double		AlgMatrixVectorNorm(double *aV, size_t nV)
{
  size_t	id0;
  double	nrm,
  		ssq = 0.0;

  for(id0 = 0; id0 < nV; ++id0)
  {
    ssq += aV[id0] * aV[id0];
  }
  nrm = sqrt(ssq);
  return(nrm);
}
