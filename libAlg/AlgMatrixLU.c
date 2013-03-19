#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixLU_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixLU.c
* \author       Richard Baldock, Bill Hill
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
* \brief        Provides functions for solving matrix equations of the
*		form: A.x = b for x, inverting a matrix, calculating
*               the determinant of a matrix and performing LU
*               decomposition.
* \ingroup   	AlgMatrix
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <float.h>


/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		3x3 libAlc double array matrix.
*		On return the matrix A is overwritten with its LU
*		decomposition and the column vector b is overwritten
*		with the solution column matrix x.
* \param	aM			Raw 3x3 array matrix A.
* \param	bV			Column vector b.
* \param	bSz			Size (number of columns) of vector b
* 					(and x).
*/
AlgError	AlgMatrixLUSolveRaw3(double **aM, double *bV, int bSz)
{
  int		k;
  int		wSpace[3];
  AlgError	errCode = ALG_ERR_NONE;

  if((errCode = AlgMatrixLUDecompRaw(aM, 3, wSpace, NULL)) == ALG_ERR_NONE)
  {
    /* Having found the LU decomposition, now perform back substitution */
    for(k = 0; k < bSz; ++k)
    {
      (void )AlgMatrixLUBackSubRaw(aM, 3, wSpace, bV +  (k * 3));
    }
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		4x4 libAlc double array matrix.
*		On return the matrix A is overwritten with its LU
*		decomposition and the column vector b is overwritten
*		with the solution column matrix x.
* \param	aM			Raw 4x4 array matrix A.
* \param	bV			Column vector b.
* \param	bSz			Size (number of columns) of vector b
* 					(and x).
*/
AlgError	AlgMatrixLUSolveRaw4(double **aM, double *bV, int bSz)
{
  int		k;
  int		wSpace[4];
  AlgError	errCode = ALG_ERR_NONE;

  if((errCode = AlgMatrixLUDecompRaw(aM, 4, wSpace, NULL)) == ALG_ERR_NONE)
  {
    /* Having found the LU decomposition, now perform back substitution */
    for(k = 0; k < bSz; ++k)
    {
      (void )AlgMatrixLUBackSubRaw(aM, 4, wSpace, bV +  (k * 4));
    }
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		square matrix.
*		On return the matrix A is overwritten with its LU
*		decomposition and the column vector b is overwritten
*		with the solution column matrix x.
* \param	aM			Square matrix.
* \param	bV			Column vector b.
* \param	bSz			Size (number of columns) of vector b
* 					(and x).
*/
AlgError	AlgMatrixLUSolve(AlgMatrix aM, double *bV, int bSz)
{
  AlgError	errCode = ALG_ERR_NONE;

  if(aM.core == NULL)
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(aM.core->type != ALG_MATRIX_RECT)
  {
    errCode = ALG_ERR_MATRIX_TYPE;
  }
  else
  {
    errCode = AlgMatrixLUSolveRaw(aM.rect->array, aM.rect->nR, bV, bSz);
  }
  return(errCode);
}
/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		square libAlc double array matrix.
*		On return the matrix A is overwritten with its LU
*		decomposition and the column vector b is overwritten
*		with the solution column matrix x.
* \param	aM			Square libalc array matrix.
* \param	aSz			Size of matrix A
* \param	bV			Column vector b.
* \param	bSz			Size (number of columns) of vector b
* 					(and x).
*/
AlgError	AlgMatrixLUSolveRaw(double **aM, int aSz, double *bV, int bSz)
{
  int		count,
  		offset;
  int		*wSpace = NULL;
  AlgError	errCode = ALG_ERR_NONE;

  if((aM == NULL) || (*aM == NULL) || (aSz <= 0) || 
     (bV == NULL) || (bSz <= 0))
  {
    errCode = ALG_ERR_FUNC;
  }
  if(errCode == ALG_ERR_NONE)
  {
    if((wSpace = (int *)AlcMalloc(sizeof(int) * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else if((errCode = AlgMatrixLUDecompRaw(aM, aSz, wSpace,
					    NULL)) == ALG_ERR_NONE)
    {
      /* Having found the LU decomposition, now perform back substitution */
      count = bSz;
      offset = 0;
      while((errCode == ALG_ERR_NONE) && (count-- > 0))
      {
	errCode = AlgMatrixLUBackSubRaw(aM, aSz, wSpace, bV + offset);
	offset += aSz;
      }
    }
    AlcFree(wSpace);
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the inverse of a 3x3 libAlc double array matrix.
* \param	aM			Raw double array of size 3x3.
*/
AlgError	AlgMatrixLUInvertRaw3(double **aM)
{
  double	w[16];
  AlgError	errCode = ALG_ERR_NONE;

  w[0] = 1.0; w[1] = 0.0; w[2] = 0.0;
  w[3] = 0.0; w[4] = 1.0; w[5] = 0.0;
  w[6] = 0.0; w[7] = 0.0; w[8] = 1.0;
  errCode = AlgMatrixLUSolveRaw3(aM, w, 3);
  aM[0][0] = w[0]; aM[0][1] = w[1]; aM[0][2] = w[2];
  aM[1][0] = w[3]; aM[1][1] = w[4]; aM[1][2] = w[5];
  aM[2][0] = w[6]; aM[2][1] = w[7]; aM[2][2] = w[8];
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the inverse of a 4x4 libAlc double array matrix.
* \param	aM			Raw double array of size 4x4.
*/
AlgError	AlgMatrixLUInvertRaw4(double **aM)
{
  double	w[16];
  AlgError	errCode = ALG_ERR_NONE;

  w[ 0] = 1.0; w[ 1] = 0.0; w[ 2] = 0.0; w[ 3] = 0.0;
  w[ 4] = 0.0; w[ 5] = 1.0; w[ 6] = 0.0; w[ 7] = 0.0;
  w[ 8] = 0.0; w[ 9] = 0.0; w[10] = 1.0; w[11] = 0.0;
  w[12] = 0.0; w[13] = 0.0; w[14] = 0.0; w[15] = 1.0;
  errCode = AlgMatrixLUSolveRaw4(aM, w, 4);
  aM[0][0] = w[ 0]; aM[0][1] = w[ 1]; aM[0][2] = w[ 2]; aM[0][3] = w[3];
  aM[1][0] = w[ 4]; aM[1][1] = w[ 5]; aM[1][2] = w[ 6]; aM[1][3] = w[7];
  aM[2][0] = w[ 8]; aM[2][1] = w[ 9]; aM[2][2] = w[10]; aM[2][3] = w[11];
  aM[3][0] = w[12]; aM[3][1] = w[13]; aM[3][2] = w[14]; aM[3][3] = w[15];
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the inverse of a square matrix.
* \param	aM			Square libAlc double array matrix A.
*/
AlgError	AlgMatrixLUInvert(AlgMatrix aM)
{
  AlgError	errCode = ALG_ERR_NONE;

  if(aM.core == NULL)
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(aM.core->type != ALG_MATRIX_RECT)
  {
    errCode = ALG_ERR_MATRIX_TYPE;
  }
  else
  {
    errCode = AlgMatrixLUInvertRaw(aM.rect->array, aM.rect->nR);
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the inverse of a square matrix.
* \param	aM			Square libAlc double array matrix A.
* \param	aSz			Size of matrix A.
*/
AlgError	AlgMatrixLUInvertRaw(double **aM, int aSz)
{
  int		idx0,
  		idx1,
  		offset;
  double	*wSpace = NULL;
  AlgError	errCode = ALG_ERR_NONE;

  if((aM == NULL) || (*aM == NULL) || (aSz <= 0))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
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
      errCode = AlgMatrixLUSolveRaw(aM, aSz, wSpace, aSz);
    }
    if(errCode == ALG_ERR_NONE)
    {
      for(idx0 = 0; idx0 < aSz; ++idx0)
      {
	offset = idx0;
	for(idx1 = 0; idx1 < aSz; ++idx1)
	{
	  aM[idx0][idx1] = *(wSpace + offset);
	  offset += aSz;
	}
      }
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the determinant of a 3x3 double libAlc array matrix.
* 		The matrix is overwitten with its LU decomposition on return.
* \param	aM			Raw 3x3 array matrix.
* \param	det			Destination ptr for the
*					determinant (may be NULL).
*/
AlgError	AlgMatrixLUDetermRaw3(double **aM, double *det)
{
  int		wSpace[3];
  AlgError	errCode = ALG_ERR_NONE;

  errCode = AlgMatrixLUDecompRaw(aM, 3, wSpace, det);
  if(det)
  {
    *det *= aM[0][0] * aM[1][1] * aM[2][2];
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the determinant of a 4x4 double libAlc array matrix.
* 		The matrix is overwitten with its LU decomposition on return.
* \param	aM			Raw 4x4 array matrix.
* \param	det			Destination ptr for the
*					determinant (may be NULL).
*/
AlgError	AlgMatrixLUDetermRaw4(double **aM, double *det)
{
  int		wSpace[4];
  AlgError	errCode = ALG_ERR_NONE;

  errCode = AlgMatrixLUDecompRaw(aM, 4, wSpace, det);
  if(det)
  {
    *det *= aM[0][0] * aM[1][1] * aM[2][2] * aM[3][3];
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the determinant of a square matrix.
* 		The matrix is overwitten with its LU decomposition
* 		on return.
* \param	aM			Matrix A.
* \param	det			Destination ptr for the
*					determinant (may be NULL).
*/
AlgError	AlgMatrixLUDeterm(AlgMatrix aM, double *det)
{
  AlgError	errCode = ALG_ERR_NONE;

  if(aM.core == NULL)
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(aM.core->type != ALG_MATRIX_RECT)
  {
    errCode = ALG_ERR_MATRIX_TYPE;
  }
  else
  {
    errCode = AlgMatrixLUDetermRaw(aM.rect->array, aM.rect->nR, det);
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Calculates the determinant of a square double libAlc array
* 		matrix. The matrix is overwitten with its LU decomposition
* 		on return.
* \param	aM			Raw matrix A.
* \param	aSz			Size of matrix A.
* \param	det			Destination ptr for the
*					determinant (may be NULL).
*/
AlgError	AlgMatrixLUDetermRaw(double **aM, int aSz, double *det)
{
  int		idx;
  int		*wSpace = NULL;
  AlgError	errCode = ALG_ERR_NONE;

  if((aM == NULL) || (*aM == NULL) || (aSz <= 0))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    if((wSpace = (int *)AlcMalloc(sizeof(int) * aSz)) == NULL)
    {
      errCode = ALG_ERR_MALLOC;
    }
    else
    {
        errCode = AlgMatrixLUDecompRaw(aM, aSz, wSpace, det);
    }
    if((errCode == ALG_ERR_NONE) && det)
    {
        idx = aSz;
	while(idx-- > 0)
	{
	    *det *= aM[idx][idx];
	}
    }
    if(wSpace)
    {
      AlcFree(wSpace);
    }
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Replaces the given matrix with the LU decomposition
*		of a row-wise permutation of itself.
*		The given index vector is used to record the
*		permutation effected by the partial pivoting.
* \param	aM			Matrix A.
* \param	iV			Index vector.
* \param	evenOdd			Set to +/-1.0 depending on
*					whether the permutation was
*					even or odd (may be NULL).
*/
AlgError	AlgMatrixLUDecomp(AlgMatrix aM, int *iV, double *evenOdd)
{
  AlgError	errCode = ALG_ERR_NONE;

  if(aM.core == NULL)
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(aM.core->type != ALG_MATRIX_RECT)
  {
    errCode = ALG_ERR_MATRIX_TYPE;
  }
  else
  {
    errCode = AlgMatrixLUDecompRaw(aM.rect->array, aM.rect->nR, iV, evenOdd);
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Replaces the given matrix with the LU decomposition
*		of a row-wise permutation of itself.
*		The given index vector is used to record the
*		permutation effected by the partial pivoting.
* \param	aM			Given raw matrix A.
* \param	aSz			Size of matrix A.
* \param	iV			Index vector.
* \param	evenOdd			Set to +/-1.0 depending on
*					whether the permutation was
*					even or odd (may be NULL).
*/
AlgError	AlgMatrixLUDecompRaw(double **aM, int aSz,
    				     int *iV, double *evenOdd)
{
  int		idx0,
  		idx1,
		idx2,
		even,
		iMax = 0;
  double	tD0,
  		aaMax,
  		sum,
  		tiny = 1.0e-20;
  double	*tDP0,
		*tDP1,
  		*wSpace = NULL;
  double	wSpace4[4];
  AlgError	errCode = ALG_ERR_NONE;

  if((aM == NULL) || (*aM == NULL) || (aSz <= 0) || (iV == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    if(aSz > 4)
    {
      if((wSpace = (double *)AlcMalloc(sizeof(double) * aSz)) == NULL)
      {
	errCode = ALG_ERR_MALLOC;
      }
    }
    else
    {
      wSpace = vSpace4;
    }
    if(errCode == ALG_ERR_NONE)
    {
      even = 1;
      /* Calculate implicit scaling */
      idx0 = 0;
      while((idx0 < aSz) && (errCode == ALG_ERR_NONE))
      {
	aaMax = 0.0;
	tDP0 = *(aM + idx0);
	for(idx1 = 0; idx1 < aSz; ++idx1)
	{
	  if((tD0 = fabs(*tDP0++)) > aaMax)
	  {
	    aaMax = tD0;
	  }
	}
	if((aaMax >= -(DBL_EPSILON)) && (aaMax <= DBL_EPSILON))
	{
	    errCode = ALG_ERR_MATRIX_SINGULAR;
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
	    tDP0 = *(aM + idx0);
	    sum = *(tDP0 + idx1);
	    if(idx0 > 0)
	    {
	      for(idx2=0; idx2 < idx0; idx2++)
	      {
		sum -= *(tDP0 + idx2) * aM[idx2][idx1];
	      }
	      *(tDP0 + idx1) = sum;
	    }
	  }
	}
	aaMax = 0.0;	  	   /* Start search for largest pivot element */
	for(idx0 = idx1; idx0 < aSz; ++idx0)
	{
	  tDP0 = *(aM + idx0);
	  sum = *(tDP0 + idx1);
	  if(idx1 > 0)
	  {
	    for(idx2=0; idx2 < idx1; ++idx2)
	    {
	      sum -= *(tDP0 + idx2) * aM[idx2][idx1];
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
	  tDP0 = *(aM + iMax);
	  tDP1 = *(aM + idx1);
	  for(idx2 = 0; idx2 < aSz; ++idx2)
	  {
	    tD0 = *tDP0;
	    *tDP0++ = *tDP1;
	    *tDP1++ = tD0;
	  }
	  even = !even;
	  wSpace[iMax] = wSpace[idx1];
	}
	iV[idx1] = iMax;
	if(idx1 != (aSz - 1)) 		      /* Divide by the pivot element */
	{
	  tDP0 =  *(aM + idx1) + idx1;
	  if(((tD0 = *tDP0) >= -(DBL_EPSILON)) && (tD0 <= DBL_EPSILON))
	  {
	    *tDP0 = tiny;
	  }
	  tD0 = 1.0 / *tDP0;
	  for(idx0=idx1+1; idx0 < aSz; idx0++)
	  {
	    aM[idx0][idx1] *= tD0;
	  }
	}
      }
      if(((tD0 = aM[aSz-1][aSz-1]) >= -(DBL_EPSILON)) &&
         (tD0 <= DBL_EPSILON))
      {
	aM[aSz-1][aSz-1] = tiny;
      }
      if(evenOdd)
      {
        *evenOdd = (even)? 1.0: -1.0;
      }
    }
    if(wSpace && (vSpace != wSpace4))
    {
      AlcFree(wSpace);
    }
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the set of of linear equations A.x = b where
*		A is input as its LU decomposition determined with
*		AlgMatrixLUDecompRaw() of a row-wise permutation of
*		itself. The given index vector is used to record the
*		permutation effected by the partial pivoting
* \param	aM			Matrix A.
* \param	iV			Index vector.
* \param	bV			Column vector b.
*/
AlgError	AlgMatrixLUBackSub(AlgMatrix aM, int *iV, double *bV)
{
  AlgError	errCode = ALG_ERR_NONE;

  if(aM.core == NULL)
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(aM.core->type != ALG_MATRIX_RECT)
  {
    errCode = ALG_ERR_MATRIX_TYPE;
  }
  else
  {
    errCode = AlgMatrixLUBackSubRaw(aM.rect->array, aM.rect->nR, iV, bV);
  }
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the set of of linear equations A.x = b where
*		A is input as its LU decomposition determined with
*		AlgMatrixLUDecompRaw() of a row-wise permutation of itself.
*		The given index vector is used to record the permutation
*		effected by the partial pivoting
* \param	aM			Given raw libAlc array matrix A.
* \param	aSz			Size of matrix A.
* \param	iV			Index vector.
* \param	bV			Column vector b.
*/
AlgError	AlgMatrixLUBackSubRaw(double **aM, int aSz,
    				      int *iV, double *bV)
{
  int 		idx0,
  		idx1,
		ii = -1,
		tI0;
  double 	sum;
  double	*tDP0,
  		*tDP1;
  AlgError	errCode = ALG_ERR_NONE;

  if((aM == NULL) || (*aM == NULL) || (aSz <= 0) ||
     (iV == NULL) || (bV == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else
  {
    for(idx0 = 0; idx0 < aSz; ++idx0)		     /* Forward substitution */
    {
      tI0 = iV[idx0];			    /* Uncramble permutation */
      sum = bV[tI0];
      bV[tI0] = bV[idx0];
      if(ii != -1)
      {
	tDP0 = *(aM + idx0) + ii;
	tDP1 = bV + ii;
	for(idx1 = ii; idx1 < idx0; ++idx1)
	{
	  sum -= *tDP0++ * *tDP1++;
        }
      }
      else if((sum < -(DBL_EPSILON)) || (sum > DBL_EPSILON))
      {
	ii = idx0;
      }
      bV[idx0] = sum;
    }
    idx0 = aSz;					      /* Now backsustitution */
    while(idx0-- > 0)
    {
      sum = bV[idx0];
      if(idx0 <= aSz)
      {
	tDP0 = *(aM + idx0) + idx0 + 1;
	tDP1 = bV + idx0 + 1;
	for(idx1= idx0 + 1; idx1 < aSz; ++idx1)
	{
	  sum -= *tDP0++ * *tDP1++;
	}
      }
      bV[idx0] = sum / aM[idx0][idx0];
    }
  }
  return(errCode);
}
