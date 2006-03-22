#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgMatrixSV_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgMatrixSV.c
* \author       Bill Hill
* \date         March 1999
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
* \brief        Provides functions for singular value decomposition.
* \ingroup      AlgMatrix
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <float.h>

static double	AlgMatrixSVPythag(double, double);


/*!
* \return	Error code.
* \ingroup      AlgMatrix
* \brief	Solves the matrix equation A.x = b for x, where A is a
*		matrix with at least as many columns as rows.
*		On return the matrix A is overwritten by the matrix U
*		in the singular value decomposition:
*		  A = U.W.V'
* \param	aMat			Matrix A.
* \param	nM			Number of rows in matrix A.
* \param	nN			Number of columns in matrix A.
* \param	bMat			Column matrix b, overwritten
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
AlgError	AlgMatrixSVSolve(double **aMat, int nM, int nN,
				 double *bMat, double tol, int *dstIC)
{
  int		cnt0,
  		cntIC;
  double	thresh,
  		wMax;
  double	*tDP0,
  		*wMat = NULL;
  double	**vMat = NULL;
  AlgError	errCode = ALG_ERR_NONE;

  if((aMat == NULL) || (bMat == NULL) || (nM <= 0) || (nM < nN))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if(((wMat = (double *)AlcCalloc(sizeof(double), nN)) == NULL) ||
          (AlcDouble2Malloc(&vMat, nM, nN) != ALC_ER_NONE))
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
    errCode = AlgMatrixSVDecomp(aMat, nM, nN, wMat, vMat);
  }
  if(errCode == ALG_ERR_NONE)
  {
    /* Find maximum singular value. */
    wMax = 0.0;
    cnt0 = nN;
    cntIC = 0;
    tDP0 = wMat;
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
    tDP0 = wMat;
    thresh = tol * wMax;
    while(cnt0-- > 0)
    {
      if(*tDP0 < thresh)
      {
	*tDP0 = 0.0;
      }
      ++tDP0;
    }
    errCode = AlgMatrixSVBackSub(aMat, nM, nN, wMat, vMat, bMat);
  }
  if(wMat)
  {
    AlcFree(wMat);
  }
  if(vMat)
  {
    AlcDouble2Free(vMat);
  }
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
* \param	nM			Number of rows in matrix A.
* \param	nN			Number of columns in matrix A.
* \param	wMat			The diagonal matrix of singular
*					values, returned as a vector.
* \param	vMat			The matrix V (not it's
*					transpose).
*/
AlgError	AlgMatrixSVDecomp(double **aMat, int nM, int nN,
				  double *wMat, double **vMat)
{
  int		cnt0,
  		flag,
  		idI,
		idJ,
		idK,
		idL,
		idM,
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
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVDecomp FE 0x%lx %d %d 0x%lx 0x%lx\n",
	   (unsigned long )aMat, nM, nN, (unsigned long )wMat,
	   (unsigned long )vMat));
  if((aMat == NULL) || (wMat == NULL) || (vMat == NULL) ||
     (nM <= 0) || (nM < nN))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((tDVec = (double *)AlcCalloc(sizeof(double), nN)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
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
	tDPP0 = aMat + idI;
	while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	{
	  tD0 = *(*tDPP0++ + idI);
	  scale += fabs(tD0);	/* scale += fabs(aMat[idK][idI]); */
	}
	if(scale > DBL_EPSILON)			   /* scale must always >= 0 */
	{
	  cnt0 = nMI;
	  tDPP0 = aMat + idI;
	  while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	  {
	    tDP0 = *tDPP0++ + idI;
	    *tDP0 /= scale;	/* aMat[idK][idI] /= scale; */
	    s += *tDP0 * *tDP0;	/* s += aMat[idK][idI] * aMat[idK][idI]; */
	  }
	  f = aMat[idI][idI];
	  g = (f > 0.0)? -(sqrt(s)): sqrt(s);
	  h = (f * g) - s;
	  aMat[idI][idI] = f - g;
	  if(idI != (nN - 1))
	  {
	    for(idJ = idL; idJ < nN; ++idJ)
	    {
	      s = 0.0;
	      cnt0 = nMI;
	      tDPP0 = aMat + idI;
	      while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	      {
	        s += *(*tDPP0 + idI) * *(*tDPP0 + idJ);
		++tDPP0;	/* s += aMat[idK][idI] * aMat[idK][idJ]; */
	      }
	      f = s / h;
	      cnt0 = nMI;
	      tDPP0 = aMat + idI;
	      while(cnt0-- > 0) /* for(idK = idI; idK < nM; ++idK) */
	      {
	        *(*tDPP0 + idJ) += f * *(*tDPP0 + idI);
		++tDPP0;	/* aMat[idK][idJ] += f * aMat[idK][idI]; */
	      }
	    }
	  }
	  cnt0 = nMI;
	  tDPP0 = aMat + idI;
	  while(cnt0-- > 0)	/* for(idK = idI; idK < nM; ++idK) */
	  {
	    *(*tDPP0++ + idI) *= scale; /* aMat[idK][idI] *= scale; */
	  }
	}
      }
      wMat[idI] = scale * g;
      g = s = scale = 0.0;
      if((idI < nM) && (idI != (nN - 1)))
      {
	cnt0 = nNL;
	tDP0 = *(aMat + idI) + idL;
	while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	{
	  scale += fabs(*tDP0++); /* scale += fabs(aMat[idI][idK]); */
	}
	if(scale > DBL_EPSILON)				/* scale always >= 0 */
	{
	  cnt0 = nNL;
	  tDP0 = *(aMat + idI) + idL;
	  while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	  {
	    *tDP0 /= scale;	/* aMat[idI][idK] /= scale; */
	    tD0 = *tDP0++;
	    s += tD0 * tD0;	/* s += aMat[idI][idK] * aMat[idI][idK]; */
	  }
	  f = aMat[idI][idL];
	  g = (f > 0.0)? -(sqrt(s)): sqrt(s);
	  h = (f * g) - s;
	  aMat[idI][idL] = f - g;
	  cnt0 = nNL;
	  tDP0 = *(aMat + idI) + idL;
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
	      tDP0 = *(aMat + idI) + idL;
	      tDP1 = *(aMat + idJ) + idL;
	      while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	      {
	        s += *tDP1++ * *tDP0++; /* s += aMat[idJ][idK] *
					        aMat[idI][idK]; */
	      }
	      cnt0 = nNL;
	      tDP0 = tDVec + idL; 
	      tDP1 = *(aMat + idJ) + idL;
	      while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	      {
	        *tDP1++ += s * *tDP0++; /* aMat[idJ][idK] += s * tDVec[idK]; */
	      }
	    }
	  }
	  cnt0 = nNL;
	  tDP0 = *(aMat + idI) + idL;
	  while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	  {
	    *tDP0++ *= scale;	/* aMat[idI][idK] *= scale; */
	  }
	}
      }
      if((tD0 = fabs(wMat[idI]) + fabs(tDVec[idI])) > aNorm)
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
	  tDP0 = *(aMat + idI) + idL;
	  tD0 = *tDP0;
	  tDPP0 = vMat + idL;
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
	    tDP0 = *(aMat + idI) + idL;
	    tDPP0 = vMat + idL;
	    while(cnt0-- > 0)	/* for(idK = idL; idK < nN; ++idK) */
	    {
	      s += *tDP0++ * *(*tDPP0++ + idJ); /* s += aMat[idI][idK] *
	      						vMat[idK][idJ]; */
	    }
	    cnt0 = nNL;
	    tDPP0 = vMat + idL;
	    while(cnt0-- > 0)	
	    {
	      tDP0 = *tDPP0++;
	      *(tDP0 + idJ) += s * *(tDP0 + idI); /* vMat[idK][idJ] += s *
	      						     vMat[idK][idI]; */
	    }
	  }
	}
	cnt0 = nNL;
	tDP0 = *(vMat + idI) + idL;
	tDPP0 = vMat + idL;
	while(cnt0-- > 0)	/* for(idJ = idL; idJ < nN; ++idJ) */
	{
	  *tDP0++ = 0.0;	   /* vMat[idI][idJ] = 0.0 */
	  *(*tDPP0++ + idI) = 0.0; /* vMat[idJ][idI] = 0.0; */
	}
      }
      vMat[idI][idI] = 1.0;
      g = tDVec[idI];
      idL = idI;
    }
    /* Accumulate left-hand transformations. */
    for(idI = nN - 1; idI >= 0; --idI)
    {
      idL = idI + 1;
      nNL = nN - idL;
      nMI = nM - idI;
      g = wMat[idI];
      if(idI < (nN - 1))
      {
	cnt0 = nNL;
	tDP0 = *(aMat + idI) + idL;
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
	    tDPP0 = aMat + idL;
	    while(cnt0-- > 0)	/* for(idK = idL; idK < nM; ++idK) */
	    {
	      tDP0 = *tDPP0++;
	      s +=  *(tDP0 + idI) * *(tDP0 + idJ); /* s += aMat[idK][idI] * 
	      						   aMat[idK][idJ]; */
	    }
	    f = (s / aMat[idI][idI]) * g;
	    cnt0 = nMI;
	    tDPP0 = aMat + idI;
	    while(cnt0-- > 0)   /* for(idK = idI; idK < nM; ++idK) */
	    {
	      tDP0 = *tDPP0++;
	      *(tDP0 + idJ) += f * *(tDP0 + idI);
	    }
	  }
	}
	cnt0 = nMI;
	tDPP0 = aMat + idI;
	while(cnt0-- > 0)   /* for(idJ = idI; idJ < nM; ++idJ) */
	{
	  *(*tDPP0++ + idI) *= g; /* aMat[idJ][idI] *= g; */
	}
      } 
      else
      {
	cnt0 = nMI;
	tDPP0 = aMat + idI;
	while(cnt0-- > 0)	/* for(idJ = idI; idJ < nM; ++idJ) */
	{
	  *(*tDPP0++ + idI) = 0.0; /* aMat[idJ][idI] = 0.0; */
	}
      }
      aMat[idI][idI] += 1.0;
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
	  if(fabs((tDVec[idL]) + aNorm) == aNorm)
	  {
	    flag = 0;
	    break;
	  }
	  if((fabs(wMat[nNM]) + aNorm) == aNorm)
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
	    if(fabs(f) + aNorm != aNorm)
	    {
	      g = wMat[idI];
	      h = AlgMatrixSVPythag(f, g);
	      wMat[idI] = h;
	      h = 1.0 / h;
	      c = g * h;
	      s = (-f * h);
	      cnt0 = nM;
	      tDPP0 = aMat;
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
	z = wMat[idK];
	if(idL == idK)
	{
	  if(z < 0.0)
	  {
	    wMat[idK] = -z;
	    cnt0 = nN;
	    tDPP0 = vMat;
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
	  errCode = ALG_ERR_SINGULAR;
	}
	else
	{
	  x = wMat[idL];
	  nNM = idK - 1;
	  y = wMat[nNM];
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
	    y = wMat[idI];
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
	    tDPP0 = vMat;
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
	    wMat[idJ] = z;
	    if(z > DBL_EPSILON)			       /* Can only be >= 0.0 */
	    {
	      z = 1.0 / z;
	      c = f * z;
	      s = h * z;
	    }
	    f = (c * g) + (s * y);
	    x = (c * y) - (s * g);
	    cnt0 = nM;
	    tDPP0 = aMat;
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
	  wMat[idK] = x;
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
* \param	wMat			The diagonal matrix of singular
*					values, returned as a vector.
* \param	vMat			The matrix V (not it's
*					transpose).
* \param	bMat			Column matrix b, overwritten by
*					column matrix x on return.
*/
AlgError	AlgMatrixSVBackSub(double **uMat, int nM, int nN,
				   double *wMat, double **vMat,
				   double *bMat)
{
  int		cnt0,
  		idJ;
  double	s;
  double	*tDP0,
		*tDP1,
  		*tDVec = NULL;
  double	**tDPP0;
  AlgError	errCode = ALG_ERR_NONE;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgMatrixSVBackSub FE 0x%lx %d %d 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )uMat, nM, nN,
	   (unsigned long )wMat, (unsigned long)vMat,
	   (unsigned long )bMat));

  if((uMat == NULL) || (wMat == NULL) || (vMat == NULL) ||
     (bMat == NULL) || (vMat == NULL) || (nM <= 0) || (nM < nN))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((tDVec = (double *)AlcCalloc(sizeof(double), nN)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  else
  {
    for(idJ = 0; idJ < nN; ++idJ) 
    {
      s = 0.0;
      if(fabs(wMat[idJ]) > DBL_EPSILON)
      {
	cnt0 = nM;
	tDPP0 = uMat;
	tDP0 = bMat;
	while(cnt0-- > 0)	/* for(idI = 0; idI < nM; ++idI) */
	{
	  s += *(*tDPP0++ + idJ) * *tDP0++; /* s += uMat[idI][idJ] *
	  					    bMat[idI]; */
	}
	s /= wMat[idJ];
      }
      tDVec[idJ] = s;
    }
    for(idJ = 0; idJ < nN; ++idJ) 
    {
      s = 0.0;
      cnt0 = nN;
      tDP0 = *(vMat + idJ);
      tDP1 = tDVec;
      while(cnt0-- > 0)		/* for(idI = 0; idI < nN; ++idI) */
      {
        s += *tDP0++ * *tDP1++; /* s += vMat[idJ][idI] * tDVec[idI]; */
      }
      bMat[idJ] = s;
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
