#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixLU.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides functions for solving matrix equations of the
*		form: A.x = b for x, inverting a matrix, calculating
*               the determinant of a matrix and performing LU
*               decomposition.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup	AlgMatrix
* @{
*/

#include <Alg.h>
#include <float.h>

/*!
* \return				Error code.
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		square matrix.
*		On return the matrix A is overwritten with its LU
*		decomposition and the column matrix b is overwitte
*		with the solution column matrix x.
* \param	aMat			Matrix A.
* \param	aSz			Size of matrix A
* \param	bMat			Column matrix b.
* \param	bSz			Size of matrix b (and x).
*/
AlgError	AlgMatrixLUSolve(double **aMat, int aSz,
				 double *bMat, int bSz)
{
  int		count,
  		offset;
  int		*wSpace = NULL;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUSolve FE 0x%lx %d 0x%lx %d\n",
	   (unsigned long )aMat, aSz, (unsigned long )bMat, bSz));
  if(aMat && *aMat && (aSz > 0) && bMat && (bSz > 0))
  {
    errCode = ALG_ERR_NONE;
    if((wSpace = (int *)AlcMalloc(sizeof(int) * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else if((errCode = AlgMatrixLUDecomp(aMat, aSz, wSpace,
					 NULL)) == ALG_ERR_NONE)
    {
      /* Having found the LU decomposition, now perform back substitution */
      count = bSz;
      offset = 0;
      while((errCode == ALG_ERR_NONE) && (count-- > 0))
      {
	errCode = AlgMatrixLUBackSub(aMat, aSz, wSpace, bMat + offset);
	offset += aSz;
      }
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUSolve FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return				Error code.
* \brief	Calculates the inverse of a square matrix.
* \param	aMat			Given matrix A.
* \param	aSz			Size of matrix A.
*/
AlgError	AlgMatrixLUInvert(double **aMat, int aSz)
{
  int		idx0,
  		idx1,
  		offset;
  double	*wSpace = NULL;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUInvert FE 0x%lx %d\n",
	   (unsigned long )aMat, aSz));
  if(aMat && *aMat && (aSz > 0))
  {
    errCode = ALG_ERR_NONE;
    if((wSpace = (double *)AlcCalloc(sizeof(double), aSz * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else
    {
      offset = 0;
      for(idx0=0; idx0 < aSz; ++idx0)
      {
        *(wSpace + offset) = 1.0;	    	    /* Make it a unit matrix */
	offset += aSz + 1;
      }
      errCode = AlgMatrixLUSolve(aMat, aSz, wSpace, aSz);
    }
    if(errCode == ALG_ERR_NONE)
    {
      for(idx0 = 0; idx0 < aSz; ++idx0)
      {
	offset = idx0;
	for(idx1 = 0; idx1 < aSz; ++idx1)
	{
	  aMat[idx0][idx1] = *(wSpace + offset);
	  offset += aSz;
	}
      }
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUInvert FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return				Error code.
* \brief	Calculates the determinant of a matrix. The matrix is
*		overwitten with its LU decomposition on exit.
* \param	aMat			Given matrix A.
* \param	aSz			Size of matrix A.
* \param	determ			Destination ptr for the
*					determinant (may be NULL).
*/
AlgError	AlgMatrixLUDeterm(double **aMat, int aSz, double *determ)
{
  int		idx;
  int		*wSpace = NULL;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUDeterm FE 0x%lx %d 0x%lx\n",
	   (unsigned long )aMat, aSz, (unsigned long )determ));
  if(aMat && *aMat && (aSz > 0))
  {
    errCode = ALG_ERR_NONE;
    if((wSpace = (int *)AlcMalloc(sizeof(int) * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else
    {
        errCode = AlgMatrixLUDecomp(aMat, aSz, wSpace, determ);
    }
    if((errCode == ALG_ERR_NONE) && determ)
    {
        idx = aSz;
	while(idx-- > 0)
	{
	    *determ *= aMat[idx][idx];
	}
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUDeterm FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return				Error code.
* \brief	Replaces the given matrix with the LU decomposition
*		of a row-wise permutation of itself.
*		The given index vector is used to record the
*		permutation effected by the partial pivoting.
* \param	aMat			Given matrix A.
* \param	aSz			Size of matrix A.
* \param	idxVec			Index vector.
* \param	evenOdd			Set to +/-1.0 depending on
*					whether the permutation was
*					even or odd (may be NULL).
*/
AlgError	AlgMatrixLUDecomp(double **aMat, int aSz,
    				  int *idxVec, double *evenOdd)
{
  int		idx0,
  		idx1,
		idx2,
		iMax,
		even;
  double	tD0,
  		aaMax,
  		sum,
  		tiny = 1.0e-20;
  double	*tDP0,
		*tDP1,
  		*wSpace = NULL;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUDecomp FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )aMat, aSz, (unsigned long )idxVec,
	   (unsigned long )evenOdd));
  if(aMat && *aMat && (aSz > 0) && idxVec)
  {
    errCode = ALG_ERR_NONE;
    if((wSpace = (double *)AlcMalloc(sizeof(double) * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else
    {
      even = 1;
      /* Calculate implicit scaling */
      idx0 = 0;
      while((idx0 < aSz) && (errCode == ALG_ERR_NONE))
      {
	aaMax = 0.0;
	tDP0 = *(aMat + idx0);
	for(idx1 = 0; idx1 < aSz; ++idx1)
	{
	  if((tD0 = fabs(*tDP0++)) > aaMax)
	  {
	    aaMax = tD0;
	  }
	}
	if((aaMax >= -(DBL_EPSILON)) && (aaMax <= DBL_EPSILON))
	{
	    errCode = ALG_ERR_SINGULAR;
	}
	else
	{
	  wSpace[idx0++] = 1.0 / aaMax;
	}
      }
    }
    if(errCode == ALG_ERR_NONE)
    {
      for(idx1 = 0; idx1 < aSz; ++idx1)
      {
	if(idx1 > 0)
	{
	  for(idx0 = 0; idx0 < idx1; ++idx0)
	  {
	    tDP0 = *(aMat + idx0);
	    sum = *(tDP0 + idx1);
	    if(idx0 > 0)
	    {
	      for(idx2=0; idx2 < idx0; idx2++)
	      {
		sum -= *(tDP0 + idx2) * aMat[idx2][idx1];
	      }
	      *(tDP0 + idx1) = sum;
	    }
	  }
	}
	aaMax = 0.0;	  	   /* Start search for largest pivot element */
	for(idx0 = idx1; idx0 < aSz; ++idx0)
	{
	  tDP0 = *(aMat + idx0);
	  sum = *(tDP0 + idx1);
	  if(idx1 > 0)
	  {
	    for(idx2=0; idx2 < idx1; ++idx2)
	    {
	      sum -= *(tDP0 + idx2) * aMat[idx2][idx1];
	    }
	    *(tDP0 + idx1) = sum;
	  }
	  tD0 = wSpace[idx0] * fabs(sum);   /* Figure of merit for the pivot */
	  if(tD0 >= aaMax)		  	    /* Check if best so far? */
	  {
	    iMax = idx0;
	    aaMax = tD0;
	  }
	}
	if(idx1 != iMax)				 /* Intechange rows? */
	{
	  tDP0 = *(aMat + iMax);
	  tDP1 = *(aMat + idx1);
	  for(idx2 = 0; idx2 < aSz; ++idx2)
	  {
	    tD0 = *tDP0;
	    *tDP0++ = *tDP1;
	    *tDP1++ = tD0;
	  }
	  even = !even;
	  wSpace[iMax] = wSpace[idx1];
	}
	idxVec[idx1] = iMax;

	if(idx1 != (aSz - 1)) 		      /* Divide by the pivot element */
	{
	  tDP0 =  *(aMat + idx1) + idx1;
	  if(((tD0 = *tDP0) >= -(DBL_EPSILON)) && (tD0 <= DBL_EPSILON))
	  {
	    *tDP0 = tiny;
	  }
	  tD0 = 1.0 / *tDP0;
	  for(idx0=idx1+1; idx0 < aSz; idx0++)
	  {
	    aMat[idx0][idx1] *= tD0;
	  }
	}
      }
      if(((tD0 = aMat[aSz-1][aSz-1]) >= -(DBL_EPSILON)) &&
         (tD0 <= DBL_EPSILON))
      {
	aMat[aSz-1][aSz-1] = tiny;
      }
      if(evenOdd)
      {
        *evenOdd = (even)? 1.0: -1.0;
      }
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUDecomp FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return				Error code.
* \brief	Solves the set of of linear equations A.x = b where
*		A is input as its LU decomposition determined with
*		AlgMatrixLUDecomp()
*		of a row-wise permutation of itself.
*		The given index vector is used to record the
*		permutation effected by the partial pivoting
* \param	aMat			Given matrix A.
* \param	aSz			Size of matrix A.
* \param	idxVec			Index vector.
* \param	bMat			Column matrix b.
*/
AlgError	AlgMatrixLUBackSub(double **aMat, int aSz,
    				  int *idxVec, double *bMat)
{
  int 		idx0,
  		idx1,
		ii = -1,
		tI0;
  double 	sum;
  double	*tDP0,
  		*tDP1;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUBackSub FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )aMat, aSz, (unsigned long )idxVec,
	   (unsigned long )bMat));
  if(aMat && *aMat && (aSz > 0) && bMat)
  {
    for(idx0 = 0; idx0 < aSz; ++idx0)		     /* Forward substitution */
    {
      tI0 = idxVec[idx0];			    /* Uncramble permutation */
      sum = bMat[tI0];
      bMat[tI0] = bMat[idx0];
      if(ii != -1)
      {
	tDP0 = *(aMat + idx0) + ii;
	tDP1 = bMat + ii;
	for(idx1 = ii; idx1 < idx0; ++idx1)
	{
	  sum -= *tDP0++ * *tDP1++;
        }
      }
      else if((sum < -(DBL_EPSILON)) || (sum > DBL_EPSILON))
      {
	ii = idx0;
      }
      bMat[idx0] = sum;
    }
    idx0 = aSz;					      /* Now backsustitution */
    while(idx0-- > 0)
    {
      sum = bMat[idx0];
      if(idx0 <= aSz)
      {
	tDP0 = *(aMat + idx0) + idx0 + 1;
	tDP1 = bMat + idx0 + 1;
	for(idx1= idx0 + 1; idx1 < aSz; ++idx1)
	{
	  sum -= *tDP0++ * *tDP1++;
	}
      }
      bMat[idx0] = sum / aMat[idx0][idx0];
    }
    errCode = ALG_ERR_NONE;
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixLUBackSub FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* @}
*/
