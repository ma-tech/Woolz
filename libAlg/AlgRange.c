#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgRange_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgRange.c
* \author       Bill Hill
* \date         May 2000
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
* \brief        Provides functions for computing the range of values
*		within a given array.
* \ingroup      AlgRange
* \todo         -
* \bug          None known.
*/


#include <Alg.h>
#include <math.h>
#include <float.h>

/*!
* \return	Error code.
* \ingroup	AlgRange
* \brief	Computes the range of the given data, ie it's minimum
*		and maximum values.
* \param	datASz			Number of elements in given data
*					array.
* \param	datA			Data array to examine for minimum
*					and maximum values.
* \param	dstMin			Destination ptr for minimum
*					value, may be NULL.
* \param	dstMax			Destination ptr for maximum
*					value, may be NULL.
*/
AlgError	AlgRange1D(int datASz, double *datA,
			   double *dstMin, double *dstMax)
{
  int		cnt;
  double	datV,
  		minV,
  		maxV;
  double	*datP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((datASz < 1) || (datA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = datASz;
    datP = datA;
    maxV = minV = *datP;
    while(--cnt > 0)
    {
      if((datV = *++datP) > maxV)
      {
        maxV = datV;
      }
      else if (datV < minV)
      {
        minV = datV;
      }
    }
    if(dstMin)
    {
      *dstMin = minV;
    }
    if(dstMax)
    {
      *dstMax = maxV;
    }
  }
  return(algErr);
}

/*!
* \return	Error code.
* \ingroup	AlgRange
* \brief	Computes the range of the given indexed data, ie it's
*		minimum and maximum values.
* \param	datA			Data array to examine for minimum
*					and maximum values.
* \param	idxASz			Number of elements in given
*					index array.
* \param	idxA			Index array with indicies into
*					the data buffer for the values
*					to examine.
* \param	dstMin			Destination ptr for minimum
*					value, may be NULL.
* \param	dstMax			Destination ptr for maximum
*					value, may be NULL.
*/
AlgError	AlgRangeIdx1D(double *datA, int idxASz, int *idxA,
			      double *dstMin, double *dstMax)
{
  int		cnt;
  double	datV,
  		minV,
  		maxV;
  int		*idxP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((idxASz < 1) || (datA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = idxASz;
    idxP = idxA;
    maxV = minV = *(datA + *idxP);
    while(--cnt > 0)
    {
      if((datV = *(datA + *++idxP)) > maxV)
      {
        maxV = datV;
      }
      else if (datV < minV)
      {
        minV = datV;
      }
    }
    if(dstMin)
    {
      *dstMin = minV;
    }
    if(dstMax)
    {
      *dstMax = maxV;
    }
  }
  return(algErr);
}

