#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgPolyLSQ.c
* \author       John Elder, Richard Baldock, Bill Hill
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
* \brief        Provides functions for fitting a polynomial using
*		least squares.
* \ingroup      AlgFit
* \todo         -
* \bug          None known.
*/

#include <Alg.h>

/*!
* \return	Error code.
* \ingroup      AlgFit
* \brief	Attempts to fit a polynomial to the given data using
*		a least squares approach.
* \param	xVec			Data vector x of size vecSz.
* \param	yVec			Data vector y of size vecSz.
* \param	vecSz			Size of data vectors.
* \param	polyDeg			Degree of ploynomial.
* \param	cVec			Destination vector for the
*					polynomial coefficients, which
*					must have at least polyDeg + 1
*					elements.
*/
AlgError	AlgPolynomialLSq(double *xVec, double *yVec,
				 int vecSz, int polyDeg, double *cVec)
{
  int		tI0,
  		tI1,
		idxI,
		idxJ,
		idxK,
		count0,
		polyDeg2;
  double	tD0,
  		tD1,
		tD2;
  double	*tDP0,
  		*tDP1,
		*sig = NULL,
  		*xPow = NULL,
		*sigXKy = NULL;
  double	**tDPP0,
  		**aMat = NULL;
  AlgError	errCode = ALG_ERR_FUNC;

  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgPolynomialLSq FE 0x%lx 0x%lx %d %d 0x%lx\n",
	   (unsigned long )xVec, (unsigned long )yVec, vecSz, polyDeg,
	   (unsigned long )cVec));
  if(xVec && yVec && cVec && (vecSz > 0) && (polyDeg > 0) && (vecSz > polyDeg))
  {
    polyDeg2 = polyDeg * 2;
    tI0 = sizeof(double) * (polyDeg + 1);
    tI1 = sizeof(double) * (polyDeg2 + 1);
    if(((sig = (double *)AlcMalloc(tI1)) == NULL) ||
       ((xPow = (double *)AlcMalloc(tI1)) == NULL) ||
       ((sigXKy  = (double *)AlcMalloc(tI0)) == NULL) ||
       (AlcDouble2Malloc(&aMat, polyDeg2 + 1, polyDeg2 + 2) != ALC_ER_NONE))
    {
      errCode = ALG_ERR_MALLOC;
    }
    else
    {
      *xPow = 1.0;
      for(idxI = 0; idxI < vecSz; ++idxI)
      {
	tD0 = *(xVec + idxI);
	tD1 = *(yVec + idxI);
	tDP0 = xPow;
	tDP1 = sig;
	count0 = polyDeg2;
	while(count0-- > 0)		/* for(k = 1; k <= polyDeg2; ++k) */
	{
	  tD2 = *tDP0 * tD0;
	  *++tDP0 = tD2; 		/* xPow[k] = xPow[k - 1] * tD0; */
	  *++tDP1 += tD2; 		/* sig[k] += xPow[k]; */
	}
	tDP0 = xPow + 2; 
	tDP1 = sigXKy;
	*tDP1++ += tD1;			/* sigxky[0] += yVec[i] */
	*tDP1++ += tD0 * tD1;		/* sigxky[1] += xVec[i] * yVec[i] */ 
	count0 = polyDeg - 1;
	while(--count0 > 0)		/* for(k = 2; k <= polyDeg; ++k) */
	{
	  *tDP1++ += *tDP0++ * tD1;	/* sigXKy[k] += xPow[k] * yVec[i]; */
	}
      }
      *sig = vecSz;
      for(idxI = 0; idxI <= polyDeg; ++idxI)
      {
	tDP0 = *(aMat + idxI);
	tDP1 = sig + idxI;
	count0 = polyDeg;
	while(count0-- >= 0)		/*  for(j = 0; j <= polyDeg; ++j) */
	{
	  *tDP0++ = *tDP1++;		/* aMat[i][j] = sig[i + j]; */
	}
      }
      tI0 = polyDeg + 1;
      tDP0 = sigXKy;
      tDPP0 = aMat;
      count0 = polyDeg;
      while(count0-- >= 0)		/* for(k = 0; k <= polyDeg; ++k) */
      {
	*(*(tDPP0++) + tI0) = *tDP0++;	/* aMat[k][polyDeg + 1] = sigXKy[k]; */
      }
      errCode = AlgMatrixGaussSolve(aMat, polyDeg + 1, cVec);
    }
    if(sig)
    {
      AlcFree(sig);
    }
    if(xPow)
    {
      AlcFree(xPow);
    }
    if(sigXKy)
    {
      AlcFree(sigXKy);
    }
    if(aMat)
    {
      AlcDouble2Free(aMat);
    }
  }
  ALG_DBG((ALG_DBG_LVL_FN|ALG_DBG_LVL_1),
	  ("AlgPolynomialLSq FX %d\n",
	   (int )errCode));
  return(errCode);
}
