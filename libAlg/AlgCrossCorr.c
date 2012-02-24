#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgCrossCorr_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgCrossCorr.c
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
* \brief	Frequency domain cross correlation functions.
* \ingroup	AlgCorr
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <math.h>
#include <float.h>


/*!
* \return	Error code.
* \ingroup	AlgCorr
* \brief	Cross correlates the given 2D double arrays leaving
*		the result in the first of the two arrays.
*		The cross correlation data are un-normalized.
* \param	data0			Data for/with obj0's FFT 
*					(source: AlcDouble2Malloc)
*					which holds the cross	
*					correlation data on return.
* \param	data1			Data for/with obj1's FFT 
*					(source: AlcDouble2Malloc).
* \param	nX			Number of columns in each of the
*					data arrays.
* \param	nY			Number of lines in each of the
*					data arrays.
*/
AlgError	AlgCrossCorrelate2D(double **data0, double **data1,
			            int nX, int nY)
{
  int		tI0,
		tI1,
		idX,
		idY,
  		nX2,
		nY2;
  double	tD1,
		tD2,
		tD3,
		tD4;
  double	*tDP1,
		*tDP2;
  double	*reBuf = NULL,
		*imBuf = NULL;
  AlgError	errNum = ALG_ERR_NONE;
  const int	cThr = 1,
  		minN = 8,
  		maxN = 1048576;

  if((data0 == NULL) || (data1 == NULL) ||
     (nX < minN) || (nX > maxN) || (nY < minN) || (nY > maxN))
  {
     errNum = ALG_ERR_FUNC;
  }
  else
  {
    (void )AlgBitNextPowerOfTwo((unsigned int *)&tI0, nX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&tI1, nY);
    if((tI0 != nX) || (tI1 != nY))
    {
      errNum = ALG_ERR_FUNC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    if(((reBuf = (double *)AlcMalloc(sizeof(double) * nY)) == NULL) ||
       ((imBuf = (double *)AlcMalloc(sizeof(double) * nY)) == NULL))
    {
      errNum = ALG_ERR_MALLOC;
    }
  }
  if(errNum == ALG_ERR_NONE)
  {
    AlgFourReal2D(data0, reBuf, imBuf, nX, nY, cThr);
    AlgFourReal2D(data1, reBuf, imBuf, nX, nY, cThr);
    nX2 = nX / 2;
    nY2 = nY / 2;
    for(idY = 0; idY < nY; ++idY)
    {
      tDP1 = *(data0 + idY) + 1;
      tDP2 = *(data1 + idY) + 1;
      for(idX = 1; idX < nX2; ++idX)
      {
	tD1 = *tDP1;
	tD2 = *(tDP1 + nX2);
	tD3 = *tDP2;
	tD4 = -*(tDP2 + nX2);
	*tDP1 = tD1 * tD3 - tD2 * tD4;
	*(tDP1 + nX2) = tD1 * tD4 + tD2 * tD3;
	++tDP1;
	++tDP2;
      }
    }
    for(idX = 0; idX < nX; idX += nX2)
    {
      for(idY = 1; idY < nY2; ++idY)
      {
	tDP1 = *(data0 + idY) + idX;
	tDP2 = *(data0 + nY2 + idY) + idX;
	tD1 = *tDP1;
	tD2 = *tDP2;
	tD3 = *(*(data1 + idY) + idX);
	tD4 = -*(*(data1 + nY2 + idY) + idX);
	*tDP1 = tD1 * tD3 - tD2 * tD4;
	*tDP2 = tD1 * tD4 + tD2 * tD3;
      }
    }
    **data0 *= **data1;
    **(data0 + nY2) *= **(data1 + nY2);
    *(*data0 + nX2) *= *(*data1 + nX2);
    *(*(data0 + nY2) + nX2) *= *(*(data1 + nY2) + nX2);
    AlgFourRealInv2D(data0, reBuf, imBuf, nX, nY, cThr);
  }
  if(reBuf)
  {
    AlcFree(reBuf);
  }
  if(imBuf)
  {
    AlcFree(imBuf);
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup	AlgCorr
* \brief	Find the maximum correlation value in the given two
*		dimensional array. Because the correlation data are 
*		stored in wrap-around order and the maximum is known
*		to lie within a limited search range, only the four
*		corners of the array are searched.		
* \param	dstMaxX			Destination ptr for column
*					coordinate with maximum value.
* \param	dstMaxY			Destination ptr for line
*					coordinate with maximum value.		
* \param	dstMaxVal		Destination ptr for maximum value.	
* \param	data			Data to search for maximum.
* \param	nX			Number of colimns in data.
* \param	nY			Number of lines in data.
* \param	searchX			Maximum number of columns to
*					search.
* \param	searchY			Maximum number of lines to
*					search.
*/
void		AlgCrossCorrPeakXY(int *dstMaxX, int *dstMaxY,
				   double *dstMaxVal,
				   double **data, int nX, int nY,
				   int searchX, int searchY)
{
  int		yIdx,
		xIdx,
		xMax,
		yMax;
  double	maxVal;
  double	*left,
		*right;

  xMax = 0;
  yMax = 0;
  maxVal = **data;
  for(yIdx = 0; yIdx < searchY; ++yIdx) 	/* Search two bottom corners */
  {
    left = *(data + yIdx);
    right = left + nX - 1;
    for(xIdx = 0; xIdx < searchX; ++xIdx)
    {
      if(*left > maxVal)
      {
	maxVal = *left;
	xMax = xIdx;
	yMax = yIdx;
      }
      if(*right > maxVal)
      {
	maxVal = *right;
	xMax = -(xIdx + 1);
	yMax = yIdx;
      }
      ++left;
      --right;
    }
  }
  for(yIdx = 1; yIdx <= searchY; ++yIdx)	   /* Search two top corners */
  {
    left = *(data + nY - yIdx);
    right = left + nX - 1;
    for(xIdx= 0; xIdx < searchX; ++xIdx)
    {
      if(*left > maxVal)
      {
	maxVal = *left;
	xMax = xIdx;
	yMax = -yIdx;
      }
      if(*right > maxVal)
      {
	maxVal = *right;
	xMax = -(xIdx + 1);
	yMax = -yIdx;
      }
      ++left;
      --right;
    }
  }
  if(dstMaxVal)
  {
    *dstMaxVal = maxVal;
  }
  if(dstMaxX)
  {
    *dstMaxX = xMax;
  }
  if(dstMaxY)
  {
    *dstMaxY = yMax;
  }
}

/*!
* \return	void
* \brief	Finds peak value in cross correlation data, only
*		searching the first column.
*		The search is particularly simple because only the
*		first column of the cross correlation data needs to be
*		searched for the angle of rotation.		
*		Data are in wrap around order.			
* \param	dstMaxY			Destination ptr for line
*					coordinate with maximum value.
*					correlation maximum.	
* \param	dstMaxVal		Cross correlation maximum.
* \param	data			Cross correlation data.	
* \param	nY			Number of lines in data array.
*/
void		AlgCrossCorrPeakY(int *dstMaxY, double *dstMaxVal,
				   double **data, int nY)
{
  int		idx,
		maxIdx;
  double	maxVal;

  maxIdx = 0;
  maxVal = (**data);
  for(idx = 1; idx < nY; ++idx)
  {
    if(**(data + idx) > maxVal)
    {
      maxVal = **(data + idx);
      maxIdx = idx;
    }
  }
  if(maxIdx > (nY / 2))
  {
    maxIdx -= nY;
  }
  if(dstMaxY)
  {
    *dstMaxY = maxIdx;
  }
  if(dstMaxVal)
  {
    *dstMaxVal = maxVal;
  }
}

