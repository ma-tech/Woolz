#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgLinearFit_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgLinearFit.c
* \author       Bill Hill
* \date         May 2000
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
* \brief	Provides functions for fitting linear models to data,
*		ie linear regression.
* \ingroup      AlgFit
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

/*!
* \return	Error code.
* \ingroup	AlgFit
* \brief	Computes the least squares best fit straight line
*		(y = a + bx) through the given data, ie linear
*		regression.
*		This function is based on the function fit():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
*  \param	datSz			Number of elements in given
*					data arrays array.
* \param	datXA			Data array with 'x' values.
* \param	datYA			Data array with 'y' values.
* \param	dstA			Destination ptr for intercept
*					'a', may be NULL.
* \param	dstB			Destination ptr for gradient
*					'b', may be NULL.
* \param	dstSigA			Destination ptr for std dev
*					of 'a', may be NULL.
* \param	dstSigB			Destination ptr for std dev
*					of 'b', may be NULL.
* \param	dstQ			Destination ptr for goodness
*					of fit, may be NULL.
*/
AlgError	AlgLinearFit1D(int datSz, double *datXA, double *datYA,
			       double *dstA, double *dstB,
			       double *dstSigA, double *dstSigB,
			       double *dstQ)
{
  int		cnt;
  double	tD0,
  		a = 0.0,
  		b = 0.0,
		t,
		sTSq = 0.0,
		chiSq = 0.0,
  		sX = 0.0,
		sY = 0.0,
		sigA,
		sigB,
		q;
  double	*datXP,
  		*datYP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((datSz < 3) || (datXA == NULL) || (datYA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      sX += *datXP++;
      sY += *datYP++;
    }
    tD0 = sX / datSz;
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      t = *datXP++ - tD0;
      sTSq += t * t;
      b += t * *datYP++;
    }
    b /= sTSq;
    a = (sY - (sX * b)) / datSz;
    sigA = sqrt((1.0 + ((sX * sX) / (datSz * sTSq))) / datSz);
    sigB = sqrt(1.0 / sTSq);
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      tD0 = *datYP++ - a - (b * *datXP++);
      chiSq += tD0 * tD0;
    }
    tD0 = sqrt(chiSq / (datSz - 2));
    sigA *= tD0;
    sigB *= tD0;
    q = AlgGammaP(0.5 * (datSz - 2), 0.5 * (chiSq), &algErr);
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(dstA)
    {
      *dstA = a;
    }
    if(dstB)
    {
      *dstB = b;
    }
    if(dstSigA)
    {
      *dstSigA = sigA;
    }
    if(dstSigB)
    {
      *dstSigB = sigB;
    }
    if(dstQ)
    {
      *dstQ = q;
    }
  }
  return(algErr);
}

/*!
* \return	Error code.
* \ingroup	AlgFit
* \brief	Computes the least squares best fit straight line
*		(y = a + bx) through the given data, ie linear
*		regression.
* \param	datXA			Data array with 'x' values.
* \param	datYA			Data array with 'y' values.
* \param	idxXA			Index array with indicies into
*					the 'x' data buffer.
*					for the values use examine.
* \param	idxYA			Index array with indicies into
*					the 'y' data buffer.
*					for the values use examine.
* \param	idxASz			Number of elements in each of
*					the given index arrays.
* \param	dstA			Destination ptr for intercept
*					'a', may be NULL.
* \param	dstB			Destination ptr for gradient
*					'b', may be NULL.
* \param	dstSigA			Destination ptr for std dev
*					of 'a', may be NULL.
* \param	dstSigB			Destination ptr for std dev
*					of 'b', may be NULL.
* \param	dstQ			Destination ptr for goodness
*					of fit, may be NULL.
*/
AlgError	AlgLinearFitIdx1D(double *datXA, double *datYA,
				  int *idxXA, int *idxYA, int idxASz,
				  double *dstA, double *dstB,
				  double *dstSigA, double *dstSigB,
				  double *dstQ)
{
  int		cnt;
  double	tD0,
  		a = 0.0,
  		b = 0.0,
		t,
		sTSq = 0.0,
		chiSq = 0.0,
  		sX = 0.0,
		sY = 0.0,
		sigA,
		sigB,
		q;
  int		*idxXP,
  		*idxYP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((idxASz < 3) || (datXA == NULL) || (datYA == NULL) ||
     (idxXA == NULL) || (idxXA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      sX += *(datXA + *idxXP++);
      sY += *(datYA + *idxYP++);
    }
    tD0 = sX / idxASz;
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      t = *(datXA + *idxXP++) - tD0;
      sTSq += t * t;
      b += t * *(datYA + *idxYP++);
    }
    b /= sTSq;
    a = (sY - (sX * b)) / idxASz;
    sigA = sqrt((1.0 + ((sX * sX) / (idxASz * sTSq))) / idxASz);
    sigB = sqrt(1.0 / sTSq);
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      tD0 = *(datYA + *idxYP++) - a - (b * *(datXA + *idxXP++));
      chiSq += tD0 * tD0;
    }
    tD0 = sqrt(chiSq / (idxASz - 2));
    sigA *= tD0;
    sigB *= tD0;
    q = AlgGammaP(0.5 * (idxASz - 2), 0.5 * (chiSq), &algErr);
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(dstA)
    {
      *dstA = a;
    }
    if(dstB)
    {
      *dstB = b;
    }
    if(dstSigA)
    {
      *dstSigA = sigA;
    }
    if(dstSigB)
    {
      *dstSigB = sigB;
    }
    if(dstQ)
    {
      *dstQ = q;
    }
  }
  return(algErr);
}
