#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixGauss_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixGauss.c
* \author       John Elder, Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief        Provides a function for solving matrix equations of the
*               form: A.x = b for x using Gaussian elimination with
*               partial pivoting.
* \ingroup      AlgMatrix
*/

#include <Alg.h>
#include <float.h>

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x by Gaussian
*		elimination with partial pivoting. Matrix A is the
*		matrix of coefficients.
* \param	aMat			The augmented matrix of size
*					aSz x (aSz + 1), whose columns
*					0 - (aSz - 1) are the
*					corresponding columns of A, and
*					aSz'th column is b. Overwritten
*					on exit.
* \param	xMat			On exit contains the solution
*					matrix x.
*/
AlgError	AlgMatrixGaussSolve(AlgMatrix aMat, double *xMat)
{
  int		idxI,
		idxK,
		count0,
		pivotRow,
		homogeneous = 0;
  double	tD0,
  		maxPivot;
  double	*tDP0,
  		*tDP1;
  double	**tDPP0;
  size_t	aSz;
  double	**abMat;
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixGaussSolve FE\n"));
  if((aMat.core == NULL) || (aMat.core->type != ALG_MATRIX_RECT) ||
     (aMat.core->nR < 1) || (aMat.core->nC < 1) ||
     (aMat.core->nR >= aMat.core->nC) || (xMat == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    aSz = aMat.rect->nR;
    abMat = aMat.rect->array;
    /* Determine homogeneity */
    tDPP0 = abMat;
    count0 = aSz;
    while((count0-- > 0) && !homogeneous)
    {
      tD0 = *(*tDPP0++ + aSz);
      homogeneous = fabs(tD0) <= DBL_EPSILON;
    }
    /* Perform the elimination */
    for(idxK = 0; (idxK < aSz) && (errCode == ALG_ERR_NONE); ++idxK)
    {
      /* Search for maximum modulus pivot */
      pivotRow = idxK;
      tDPP0 = abMat + idxK;
      maxPivot = fabs(*(*tDPP0 + idxK));
      for(idxI = idxK + 1; idxI < aSz; ++idxI)
      {
	if((tD0 = fabs(*(*++tDPP0 + idxK))) > maxPivot)
	{
	  maxPivot = tD0;
	  pivotRow = idxI;
	}
      }
      if(maxPivot <= DBL_EPSILON)
      {
	errCode = ALG_ERR_MATRIX_SINGULAR;
      }
      else
      {
	/* Interchange rows idxK and pivotRow */
	if(pivotRow != idxK)
	{
	  tDP0 = *(abMat + idxK) + idxK;
	  tDP1 = *(abMat + pivotRow) + idxK;
	  count0 = aSz - idxK;
	  while(count0-- >= 0)	/* for(j = k; j < aSz + 1; ++j) */
	  {
	    tD0 = *tDP0;		/* tmp = ab[k][j]; */
	    *tDP0++ = *tDP1;		/* ab[k][j] = ab[pivotRow][j]; */
	    *tDP1++ = tD0;		/* ab[pivotRow][j] = tmp; */
	  }
	}
	/* Eliminate xMat[idxK] from equations k + 1, ..., aSz */
	tDP0 = *(abMat + idxK);
	tD0 = *(tDP0 + idxK);
	tDP0 += aSz;
	count0 = aSz - idxK;
	while(count0-- >= 0) 		/* for(j = aSz; j >= k; --j) */
	{
	  *tDP0-- /= tD0;		/* ab[k][j] /= ab[k][k]; */
	}
	for(idxI = idxK + 1; idxI < aSz; ++idxI)
	{
	  tDP0 = *(abMat + idxI);
	  tDP1 = *(abMat + idxK) + idxK + 1;
	  tD0 = *(tDP0 + idxK);
	  tDP0 += idxK + 1;
	  count0 = aSz - idxK;
	  while(count0-- > 0)		/* for(j = k + 1; j <= aSz; ++j) */
	  {
	    *tDP0++ -= tD0 * *tDP1++;	/* ab[i][j] -= ab[i][k] * ab[k][j]; */
	  }
	}
      }
    }
    if(errCode == ALG_ERR_NONE)
    {
      if(homogeneous)
      {
	count0 = aSz;
	tDP0 = xMat;
	while(count0-- > 0)
	{
	  *tDP0++ = 0.0;
	}
	errCode = ALG_ERR_MATRIX_HOMOGENEOUS;
      }
      else
      {
	/* Perform backward substitution */
	for(idxI = aSz - 1; idxI >= 0; idxI--)
	{
	  tDPP0 = abMat + idxI - 1;
	  tDP0 = *(abMat + idxI);
	  tD0 = *(tDP0 + aSz) / *(tDP0 + idxI);
	  xMat[idxI] = tD0;		/* x[i] = ab[i][aSz] / ab[i][i]; */
	  count0 = idxI - 1;
	  while(count0-- >= 0)	/* for(k = i - 1; k >= 0; --k) */
	  {
	    tDP1 = *tDPP0--;
	    *(tDP1 + aSz) -= *(tDP1 + idxI) *
			     tD0;	/* ab[k][aSz] -= ab[k][i] * x[i]; */
	  }
	}
      }
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixGaussSolve FX %d\n",
	   (int )errCode));
  return(errCode);
}
