#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMatrixSV_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMatrixSV.c
* \author       Bill Hill
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
* \brief        Provides functions for singular value decomposition.
* \ingroup      AlgMatrix
*/

#include <Alg.h>
#include <float.h>

#define ALG_FAST_CODE

static double	AlgMatrixSVPythag(double, double);

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		matrix with at least as many rows as columns.
*		On return the matrix A is overwritten by the matrix U
*		in the singular value decomposition:
*		  A = U.W.V'
* \param	aMat			Matrix A.
* \param	bVec			Column matrix b, overwritten
*					by matrix x on return.
* \param	tol			Tolerance for singular values,
*					1.0e-06 should be suitable as
*					a default value.
* \param	dstIC			Destination pointer for
*					ill-conditioned flag, which may
*					be NULL, but if not NULL is set
*					to the number of singular values
*					which are smaller than the maximum
*					singular value times the given
*					threshold.
*/
AlgError	AlgMatrixSVSolve(AlgMatrix aMat, double *bVec, double tol,
                                 int *dstIC)
{
  int		cnt0 = 0,
  		cntIC = 0;
  size_t	nM,
  		nN;
  double	thresh,
  		wMax = 0.0;
  double	*tDP0,
  		*wVec = NULL;
  AlgMatrix	vMat;
  AlgError	errCode = ALG_ERR_NONE;

  nM = aMat.core->nR;
  nN = aMat.core->nC;
  vMat.core = NULL;
  if((aMat.core == NULL) || (aMat.core->type != ALG_MATRIX_RECT) ||
     (aMat.core->nR <= 0) || (aMat.core->nR < aMat.core->nC) ||
     (bVec == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((wVec = (double *)AlcCalloc(sizeof(double), aMat.rect->nC)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
    vMat.rect = AlgMatrixRectNew(nM, nN, &errCode);
  }
  if(errCode == ALG_ERR_NONE)
  {
    errCode = AlgMatrixSVDecomp(aMat, wVec, vMat);
  }
  if(errCode == ALG_ERR_NONE)
  {
    /* Find maximum singular value. */
    cnt0 = nN;
    cntIC = 0;
    tDP0 = wVec;
    while(cnt0-- > 0)
    {
      if(*tDP0 > wMax)
      {
	++cntIC;
        wMax = *tDP0;
      }
      ++tDP0;
    }
    /* Edit the singular values, replacing any less than tol * max singular
       value with 0.0. */
    cnt0 = nN;
    tDP0 = wVec;
    thresh = tol * wMax;
    while(cnt0-- > 0)
    {
      if(*tDP0 < thresh)
      {
	*tDP0 = 0.0;
      }
      ++tDP0;
    }
    errCode = AlgMatrixSVBackSub(aMat, wVec, vMat, bVec);
  }
  AlcFree(wVec);
  AlgMatrixRectFree(vMat.rect);
  if(dstIC)
  {
    *dstIC = cntIC;
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVSolve FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Performs singular value decomposition of the given
*		matrix (A) and computes two additional matricies
*		(U and V) such that:
*		  A = U.W.V'
*		where V' is the transpose of V.
*		The code for AlgMatrixSVDecomp was derived from:
*		Numerical Recipies function svdcmp(), EISPACK
*		subroutine SVD(), CERN subroutine SVD() and ACM
*		algorithm 358.
*		See AlgMatrixSVSolve() for a usage example.
* \param	aMat			The given matrix A, and U on
*					return.
* \param	wMat			The diagonal matrix of singular
*					values, returned as a vector.
* \param	vMat			The matrix V (not it's
*					transpose).
*/
AlgError	AlgMatrixSVDecomp(AlgMatrix aMat, double *wVec, AlgMatrix vMat)
{
  int		cnt0,
  		flag,
  		idI,
		idJ,
		idK,
		idL = 0,
		its,
		nNL,
		nMI,
		nNM;
  double 	tD0,
  		c,
 		f,
		h,
		s,
		x,
		y,
		z,
		aNorm = 0.0,
		g = 0.0,
		scale = 0.0;
  double	*tDP0,
		*tDP1,
  		*tDVec = NULL;
  double	**tDPP0;
  const int	maxIts = 100;  /* Maximum iterations to find singular value. */
  const double	aScale = 0.01;  /* Used with aNorm to test for small values. */
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVDecomp FE\n"));
  if((aMat.core == NULL) || (aMat.core->type != ALG_MATRIX_RECT) ||
     (vMat.core == NULL) || (aMat.core->type != vMat.core->type) ||
     (aMat.core->nR <= 0) || (aMat.core->nR < aMat.core->nC) ||
     (aMat.core->nR != vMat.core->nR) || (aMat.core->nC != vMat.core->nC) ||
     (wVec == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((tDVec = (double *)AlcCalloc(sizeof(double), aMat.rect->nC)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
    int	   nM,
    	   nN;
    double **aAry,
    	   **vAry;

    nM = aMat.rect->nR;
    nN = aMat.rect->nC;
    aAry = aMat.rect->array;
    vAry = vMat.rect->array;
    /* Householder reduction to bidiagonal form */
    for(idI = 0; idI < nN; ++idI)
    {
      idL = idI + 1;
      nNL = nN - idL;
      nMI = nM - idI;
      tDVec[idI] = scale * g;
      g = s = scale = 0.0;
      if(idI < nM)
      {
        cnt0 = nMI;
	tDPP0 = aAry + idI;
	while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	{
	  tD0 = *(*tDPP0++ + idI);
	  scale += fabs(tD0);	/* scale += fabs(aMat[idK][idI]); */
	}
	if(scale > DBL_EPSILON)			   /* scale must always >= 0 */
	{
	  cnt0 = nMI;
	  tDPP0 = aAry + idI;
	  while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	  {
	    tDP0 = *tDPP0++ + idI;
	    *tDP0 /= scale;	/* aMat[idK][idI] /= scale; */
	    s += *tDP0 * *tDP0;	/* s += aMat[idK][idI] * aMat[idK][idI]; */
	  }
	  f = aAry[idI][idI];
	  g = (f > 0.0)? -(sqrt(s)): sqrt(s);
	  h = (f * g) - s;
	  aAry[idI][idI] = f - g;
	  if(idI != (nN - 1))
	  {
	    for(idJ = idL; idJ < nN; ++idJ)
	    {
	      s = 0.0;
	      cnt0 = nMI;
	      tDPP0 = aAry + idI;
	      while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	      {
	        s += *(*tDPP0 + idI) * *(*tDPP0 + idJ);
		++tDPP0;	/* s += aMat[idK][idI] * aMat[idK][idJ]; */
	      }
	      f = s / h;
	      cnt0 = nMI;
	      tDPP0 = aAry + idI;
	      while(cnt0-- > 0) /* for(idK = idI; idK < nM; ++idK) */
	      {
	        *(*tDPP0 + idJ) += f * *(*tDPP0 + idI);
		++tDPP0;	/* aMat[idK][idJ] += f * aMat[idK][idI]; */
	      }
	    }
	  }
	  cnt0 = nMI;
	  tDPP0 = aAry + idI;
	  while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	  {
	    *(*tDPP0++ + idI) *= scale; /* aMat[idK][idI] *= scale; */
	  }
	}
      }
      wVec[idI] = scale * g;
      g = s = scale = 0.0;
      if((idI < nM) && (idI != (nN - 1)))
      {
	cnt0 = nNL;
	tDP0 = *(aAry + idI) + idL;
	while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	{
	  scale += fabs(*tDP0++); /* scale += fabs(aMat[idI][idK]); */
	}
	if(scale > DBL_EPSILON)				/* scale always >= 0 */
	{
	  cnt0 = nNL;
	  tDP0 = *(aAry + idI) + idL;
	  while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	  {
	    *tDP0 /= scale;	/* aMat[idI][idK] /= scale; */
	    tD0 = *tDP0++;
	    s += tD0 * tD0;	/* s += aMat[idI][idK] * aMat[idI][idK]; */
	  }
	  f = aAry[idI][idL];
	  g = (f > 0.0)? -(sqrt(s)): sqrt(s);
	  h = (f * g) - s;
	  aAry[idI][idL] = f - g;
	  cnt0 = nNL;
	  tDP0 = *(aAry + idI) + idL;
	  tDP1 = tDVec + idL;
	  while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	  {
	    *tDP1++ = *tDP0++ / h; /* tDVec[idK] = aMat[idI][idK] / h; */
	  }
	  if(idI != (nM - 1))
	  {
	    for(idJ = idL; idJ < nM; ++idJ)
	    {
	      s = 0.0;
	      cnt0 = nNL;
	      tDP0 = *(aAry + idI) + idL;
	      tDP1 = *(aAry + idJ) + idL;
	      while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	      {
	        s += *tDP1++ * *tDP0++; /* s += aMat[idJ][idK] *
					        aMat[idI][idK]; */
	      }
	      cnt0 = nNL;
	      tDP0 = tDVec + idL; 
	      tDP1 = *(aAry + idJ) + idL;
	      while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	      {
	        *tDP1++ += s * *tDP0++; /* aMat[idJ][idK] += s * tDVec[idK]; */
	      }
	    }
	  }
	  cnt0 = nNL;
	  tDP0 = *(aAry + idI) + idL;
	  while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	  {
	    *tDP0++ *= scale;	/* aMat[idI][idK] *= scale; */
	  }
	}
      }
      if((tD0 = fabs(wVec[idI]) + fabs(tDVec[idI])) > aNorm)
      {
        aNorm = tD0;
      }
    }
    /* Accumulate right-hand transformations. */
    for(idI = nN - 1; idI >= 0; --idI)
    {
      if(idI < (nN - 1))
      {
	nNL = nN - idL;
	if(fabs(g) > DBL_EPSILON)
	{
	  cnt0 = nNL;
	  tDP0 = *(aAry + idI) + idL;
	  tD0 = *tDP0;
	  tDPP0 = vAry + idL;
	  while(cnt0-- > 0)	/* for(idJ = idL; idJ < nN; ++idJ) */
	  {
	    			/* vMat[idJ][idI] = (aMat[idI][idJ] /
				 		     aMat[idI][idL]) /g */
	    *(*tDPP0++ + idI) = (*tDP0++ / tD0) / g; /* Double division to try
	    						and avoid underflow */
	  }
	  for(idJ = idL; idJ < nN; ++idJ)
	  {
	    s = 0.0;
	    cnt0 = nNL;
	    tDP0 = *(aAry + idI) + idL;
	    tDPP0 = vAry + idL;
	    while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	    {
	      s += *tDP0++ * *(*tDPP0++ + idJ); /* s += aMat[idI][idK] *
	      						vMat[idK][idJ]; */
	    }
	    cnt0 = nNL;
	    tDPP0 = vAry + idL;
	    while(cnt0-- > 0)	
	    {
	      tDP0 = *tDPP0++;
	      *(tDP0 + idJ) += s * *(tDP0 + idI); /* vMat[idK][idJ] += s *
	      						     vMat[idK][idI]; */
	    }
	  }
	}
	cnt0 = nNL;
	tDP0 = *(vAry + idI) + idL;
	tDPP0 = vAry + idL;
	while(cnt0-- > 0)	/* for(idJ = idL; idJ < nN; ++idJ) */
	{
	  *tDP0++ = 0.0;	   /* vMat[idI][idJ] = 0.0 */
	  *(*tDPP0++ + idI) = 0.0; /* vMat[idJ][idI] = 0.0; */
	}
      }
      vAry[idI][idI] = 1.0;
      g = tDVec[idI];
      idL = idI;
    }
    /* Accumulate left-hand transformations. */
    for(idI = nN - 1; idI >= 0; --idI)
    {
      idL = idI + 1;
      nNL = nN - idL;
      nMI = nM - idI;
      g = wVec[idI];
      if(idI < (nN - 1))
      {
	cnt0 = nNL;
	tDP0 = *(aAry + idI) + idL;
	while(cnt0-- > 0)	/* for(idJ = idL; idJ < nN; ++idJ) */
	{
	  *tDP0++ = 0.0;	/* aMat[idI][idJ] = 0.0; */
	}
      }
      if(fabs(g) > DBL_EPSILON)
      {
	g = 1.0 / g;
	if(idI != (nN - 1))
	{
	  for(idJ = idL; idJ < nN; ++idJ)
	  {
	    s = 0.0;
	    cnt0 = nMI - 1;
	    tDPP0 = aAry + idL;
	    while(cnt0-- > 0)	/* for(idK = idL; idK < nM; ++idK) */
	    {
	      tDP0 = *tDPP0++;
	      s +=  *(tDP0 + idI) * *(tDP0 + idJ); /* s += aMat[idK][idI] * 
	      						   aMat[idK][idJ]; */
	    }
	    f = (s / aAry[idI][idI]) * g;
	    cnt0 = nMI;
	    tDPP0 = aAry + idI;
	    while(cnt0-- > 0)   /* for(idK = idI; idK < nM; ++idK) */
	    {
	      tDP0 = *tDPP0++;
	      *(tDP0 + idJ) += f * *(tDP0 + idI);
	    }
	  }
	}
	cnt0 = nMI;
	tDPP0 = aAry + idI;
	while(cnt0-- > 0)   /* for(idJ = idI; idJ < nM; ++idJ) */
	{
	  *(*tDPP0++ + idI) *= g; /* aMat[idJ][idI] *= g; */
	}
      } 
      else
      {
	cnt0 = nMI;
	tDPP0 = aAry + idI;
	while(cnt0-- > 0)	/* for(idJ = idI; idJ < nM; ++idJ) */
	{
	  *(*tDPP0++ + idI) = 0.0; /* aMat[idJ][idI] = 0.0; */
	}
      }
      aAry[idI][idI] += 1.0;
    }
    /* Diagonalize the bidiagonal form. */
    for(idK = nN - 1; (idK >= 0) && (errCode == ALG_ERR_NONE); --idK)
    {
      for(its = 1; its <= maxIts; ++its)
      {
	flag = 1;
	for(idL = idK; idL >= 0; --idL)
	{
	  nNM = idL - 1;
	  if(fabs((tDVec[idL] * aScale) + aNorm) == aNorm)
	  {
	    flag = 0;
	    break;
	  }
	  if((fabs(wVec[nNM] * aScale) + aNorm) == aNorm)
	  {
	    break;
	  }
	}
	if(flag)
	{
	  c = 0.0;
	  s = 1.0;
	  for(idI = idL; idI <= idK; ++idI)
	  {
	    f = s * tDVec[idI];

	    if(fabs(f * aScale) + aNorm != aNorm)
	    {
	      g = wVec[idI];
	      h = AlgMatrixSVPythag(f, g);
	      wVec[idI] = h;
	      h = 1.0 / h;
	      c = g * h;
	      s = (-f * h);
	      cnt0 = nM;
	      tDPP0 = aAry;
	      while(cnt0-- > 0)	/* for(idJ = 0; idJ < nM; ++idJ) */
	      {
	        tDP0 = *tDPP0 + nNM;
		tDP1 = *tDPP0++ + idI;
		y = *tDP0;	/* y = aMat[idJ][nNM]; */
		z = *tDP1;	/* z = aMat[idJ][idI]; */
		*tDP0 = (y * c) + (z * s); /* aMat[idJ][nNM] = (y * c) +
							       (z * s); */
		*tDP1 = (z * c) - (y * s); /* aMat[idJ][idI] = (z * c) -
							       (y * s); */
	      }
	    }
	  }
	}
	/* Test for convergence. */
	z = wVec[idK];
	if(idL == idK)
	{
	  if(z < 0.0)
	  {
	    wVec[idK] = -z;
	    cnt0 = nN;
	    tDPP0 = vAry;
	    while(cnt0-- > 0)	/* for(idJ = 0; idJ < nN; ++idJ) */
	    {
	      tDP0 = *tDPP0++ + idK;
	      *tDP0 = -*tDP0; /* vMat[idJ][idK] = -vMat[idJ][idK]; */
	    }
	  }
	  break;
	}
	if(its >= maxIts)
	{
	  errCode = ALG_ERR_MATRIX_SINGULAR;
	}
	else
	{
	  x = wVec[idL];
	  nNM = idK - 1;
	  y = wVec[nNM];
	  g = tDVec[nNM];
	  h = tDVec[idK];
	  f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
	  g = AlgMatrixSVPythag(f, 1.0);
	  f = (((x - z) * (x + z)) + 
	      (h * ((y / (f + ((f > 0)? g: -g))) - h))) / x;
	  c = s = 1.0;
	  for(idJ = idL; idJ <= nNM; ++idJ)
	  {
	    idI = idJ + 1;
	    g = tDVec[idI];
	    y = wVec[idI];
	    h = s * g;
	    g = c * g;
	    z = AlgMatrixSVPythag(f, h);
	    tDVec[idJ] = z;
	    c = f / z;
	    s = h / z;
	    f = (x * c) + (g * s);
	    g = (g * c) - (x * s);
	    h = y * s;
	    y = y * c;
	    cnt0 = nN;
	    tDPP0 = vAry;
	    while(cnt0-- > 0)	/* for(idM = 0; idM < nN; ++idM) */
	    {
	      tDP0 = *tDPP0 + idJ;
	      tDP1 = *tDPP0++ + idI;
	      x = *tDP0; 	/* x = vMat[idM][idJ]; */
	      z = *tDP1;	/* z = vMat[idM][idI]; */
	      *tDP0 = (x * c) + (z * s); /* vMat[idM][idJ] = (x * c) +
	      						     (z * s); */
	      *tDP1 = (z * c) - (x * s); /* vMat[idM][idI] = (z * c) -
	      						     (x * s); */
	    }
	    z = AlgMatrixSVPythag(f, h);
	    wVec[idJ] = z;
	    if(z > DBL_EPSILON)			       /* Can only be >= 0.0 */
	    {
	      z = 1.0 / z;
	      c = f * z;
	      s = h * z;
	    }
	    f = (c * g) + (s * y);
	    x = (c * y) - (s * g);
	    cnt0 = nM;
	    tDPP0 = aAry;
	    while(cnt0-- > 0)   /* for(idM = 0; idM < nM; ++idM) */
	    {
	      tDP0 = *tDPP0 + idJ;
	      tDP1 = *tDPP0++ + idI;
	      y = *tDP0;	/* y = aMat[idM][idJ]; */
	      z = *tDP1;	/* z = aMat[idM][idI]; */
	      *tDP0 = (y * c) + (z * s); /* aMat[idM][idJ] = (y * c) +
	      						     (z * s); */
	      *tDP1 = (z * c) - (y * s); /* aMat[idM][idI] = (z * c) -
	      						     (y * s); */
	    }
	  }
	  tDVec[idL] = 0.0;
	  tDVec[idK] = f;
	  wVec[idK] = x;
        }
      }
    }
    AlcFree(tDVec);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVDecomp FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the set of of linear equations A.x = b where
*		A is input as its singular value decomposition in
*		the three matricies U, W and V, as returned by
*		AlgMatrixSVDecomp().
*		The code for AlgMatrixSVBackSub was derived from:
*		Numerical Recipies function svbksb().
* \param	uMat			Given matrix U.
* \param	nM 			Number of rows in matrix U and
*       				number of elements in matrix B.
* \param	nN 			Number of columns in matricies
*					U and V, also the number of
*					elements in matricies W and x.
* \param	wVec			The diagonal matrix of singular
*					values, returned as a vector.
* \param	vMat			The matrix V (not it's
*					transpose).
* \param	bVec			Column matrix b, overwritten by
*					column matrix x on return.
*/
AlgError	AlgMatrixSVBackSub(AlgMatrix uMat, double *wVec, AlgMatrix vMat,
				   double *bVec)
{
  int		cnt0,
  		idJ;
  double	s;
  double	*tDP0,
		*tDP1,
  		*tDVec = NULL;
  double	**tDPP0;
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1), ("AlgMatrixSVBackSub FE\n"));
  if((uMat.core == NULL) || (uMat.core->type != ALG_MATRIX_RECT) ||
     (vMat.core == NULL) || (vMat.core->type != vMat.core->type) ||
     (uMat.core->nR <= 0) || (uMat.core->nR < uMat.core->nC) ||
     (uMat.core->nR != vMat.core->nR) || (uMat.core->nC != vMat.core->nC) ||
     (wVec == NULL) || (bVec == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((tDVec = (double *)AlcCalloc(sizeof(double), uMat.rect->nC)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
#ifndef ALG_FAST_CODE
    int	   	idI;
#endif
    int 	nM,
    	   	nN;
    double 	**uAry,
    	   	**vAry;

    nM = uMat.rect->nR;
    nN = uMat.rect->nC;
    uAry = uMat.rect->array;
    vAry = vMat.rect->array;
    for(idJ = 0; idJ < nN; ++idJ) 
    {
      s = 0.0;
      if(fabs(wVec[idJ]) > DBL_EPSILON)
      {
#ifdef ALG_FAST_CODE
	cnt0 = nM;
	tDPP0 = uAry;
	tDP0 = bVec;
	while(cnt0-- > 0)
	{
	  s += *(*tDPP0++ + idJ) * *tDP0++;
	}
#else
	for(idI = 0; idI < nM; ++idI)
	{
	  s += uAry[idI][idJ] * bVec[idI];
        }
#endif
	s /= wVec[idJ];
      }
      tDVec[idJ] = s;
    }
    for(idJ = 0; idJ < nN; ++idJ) 
    {
      s = 0.0;
#ifdef ALG_FAST_CODE
      cnt0 = nN;
      tDP0 = *(vAry + idJ);
      tDP1 = tDVec;
      while(cnt0-- > 0)
      {
        s += *tDP0++ * *tDP1++;
      }
#else
      for(idI = 0; idI < nN; ++idI)
      {
	s += vAry[idJ][idI] * tDVec[idI];
      }
#endif
      bVec[idJ] = s;
    }
    AlcFree(tDVec);
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVBackSub FX %d\n",
	   (int )errCode));
  return(errCode);
}

/*!
* \return	Square root of sum of squares.
* \ingroup      AlgMatrix
* \brief	Computes sqrt(size0^2 + size1^2) without underflow or
*		overflow.
* \param	side0			Length of first side.
* \param	side1			Length of second side1.
*/
static double	AlgMatrixSVPythag(double sd0, double sd1)
{
  double	tD0,
		abs0,
		abs1,
		hyp = 0.0;

  if((abs0 = fabs(sd0)) > (abs1 = fabs(sd1)))
  {
    tD0 = abs1 / abs0;
    hyp = abs0 * sqrt(1.0 + (tD0 * tD0));
  }
  else if(abs0 > DBL_EPSILON)
  {
    tD0 = abs0 / abs1;
    hyp = abs1 * sqrt(1.0 + (tD0 * tD0));
  }
  return(hyp);
}
