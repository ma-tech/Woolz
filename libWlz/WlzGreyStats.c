#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyStats_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGreyStats.c
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
* \brief	Calculates simple statistics about an object's 
* 		grey values.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Number of bytes required to store a single grey value or zero
* 		on error (invalid grey type).
* \ingroup	WlzFeatures
* \brief	Given a grey type, computes the number of bytes required to
* 		to store a single grey value of that type. Note the grey
* 		type WLZ_GREY_BIT is able to store multiple values in a
* 		single byte but one byte is required to store a single bit.
* \param	gType			Given grey type.
*/
size_t		WlzGreySize(WlzGreyType gType)
{
  size_t	sz = 0;

  switch(gType)
  {
    case WLZ_GREY_LONG:
      sz = sizeof(long);
      break;
    case WLZ_GREY_INT:
      sz = sizeof(int);
      break;
    case WLZ_GREY_SHORT:
      sz = sizeof(short);
      break;
    case WLZ_GREY_UBYTE:
      sz = sizeof(WlzUByte);
      break;
    case WLZ_GREY_FLOAT:
      sz = sizeof(float);
      break;
    case WLZ_GREY_DOUBLE:
      sz = sizeof(double);
      break;
    case WLZ_GREY_BIT:
      sz = sizeof(WlzUByte);
      break;
    case WLZ_GREY_RGBA:
      sz = sizeof(WlzUInt);
      break;
    default:
      break;
  }
  return(sz);
}

/*!
* \return	Object area or -1 on error.
* \ingroup	WlzFeatures
* \brief	Calculates simple quick statistics for the given
*               2D domain object with grey values.
* \param	srcObj			2D domain object from which to
*                                       calculate the statistics.
* \param	dstGType		Pointer for grey type.
* \param	dstMin			Pointer for minimum value.
* \param	dstMax			Pointer for maximum value.
* \param	dstSum			Pointer for sum of values.
* \param	dstSumSq		Pointer for sum of squares of
*                                       values.
* \param	dstErr			Destination pointer for error, may
*					be NULL.
*/
static int	WlzGreyStats2D(WlzObject *srcObj,
			       WlzGreyType *dstGType,
			       double *dstMin, double *dstMax,
			       double *dstSum, double *dstSumSq,
			       WlzErrorNum *dstErr)
{
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace gWSp;
  WlzGreyP	gPix;
  int		count,
  		area = 0;
  double	gVal,
  		min,
		max,
		sum = 0.0,
		sumSq = 0.0;
  WlzUInt	rgbVal;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  if((errNum = WlzInitGreyScan(srcObj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
  {
    while((WlzNextGreyInterval(&iWSp) == 0) && (errNum == WLZ_ERR_NONE))
    {
      gPix = gWSp.u_grintptr;
      gType = gWSp.pixeltype;
      count = iWSp.rgtpos - iWSp.lftpos + 1;
      while(count-- > 0)
      {
	switch(gType)
	{
	  case WLZ_GREY_INT:
	    gVal = *(gPix.inp)++;
	    break;
	  case WLZ_GREY_SHORT:
	    gVal = *(gPix.shp)++;
	    break;
	  case WLZ_GREY_UBYTE:
	    gVal = *(gPix.ubp)++;
	    break;
	  case WLZ_GREY_FLOAT:
	    gVal = *(gPix.flp)++;
	    break;
	  case WLZ_GREY_DOUBLE:
	    gVal = *(gPix.dbp)++;
	    break;
	  case WLZ_GREY_RGBA: /* RGBA - make this the modulus
			         For RGB stats call WlzRGBAGreyStats() */
	    rgbVal = *(gPix.rgbp)++;
	    gVal = WLZ_RGBA_MODULUS(rgbVal);
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	}
	if(area == 0)
	{
	  min = gVal;
	  max = gVal;
	}
	else
	{
	  if(gVal < min)
	  {
	    min = gVal;
	  }
	  else if(gVal > max)
	  {
	    max = gVal;
	  }
	}
	sum += gVal;
	sumSq += gVal * gVal;
        ++area;
      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
      *dstGType = gType;
      *dstMin = min;
      *dstMax = max;
      *dstSum = sum;
      *dstSumSq = sumSq;
  }
  else
  {
    area = -1;
  }
  *dstErr = errNum;
  return(area);
}

/*!
* \return	Object area or -1 on error.
* \ingroup	WlzFeatures
* \brief	Calculates simple quick statistics for given 2D or
*               3D domain object with grey values.
*               Pointers provided for results may be NULL without
*               causing an error.
* \param	srcObj			Object from which to calculate
*                                       the statistics.
* \param	dstGType		Pointer for grey type.
* \param	dstMin			Pointer for minimum value.
* \param	dstMax			Pointer for maximum value.
* \param	dstSum			Pointer for sum of values.
* \param	dstSumSq		Pointer for sum of squares of
*                                       values.
* \param	dstMean			Mean value.
* \param	dstStdDev		Standard deviation of values.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL if not
*                                       required.
*/
int		WlzGreyStats(WlzObject *srcObj,
			     WlzGreyType *dstGType,
			     double *dstMin, double *dstMax,
			     double *dstSum, double *dstSumSq,
			     double *dstMean, double *dstStdDev,
			     WlzErrorNum *dstErr)
{
  int		area = 0,
		pIdx,
		pCnt;
  double	min,
		max,
  		mean = -1.0,
		stdDev = -1.0,
		sum = 0.0,
		sumSq = 0.0;
  WlzGreyType	gType;
  WlzDomain	*dom;
  WlzValues	*val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzGreyStats FE 0x%lx "
	  "0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	  (unsigned long )srcObj,
	  (unsigned long )dstGType,
	  (unsigned long )dstMin, (unsigned long )dstMax,
	  (unsigned long )dstSum, (unsigned long )dstSumSq,
	  (unsigned long )dstMean, (unsigned long )dstStdDev,
	  (unsigned long )dstErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        area = WlzGreyStats2D(srcObj, &gType, &min, &max,
			      &sum, &sumSq, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
	if(srcObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if((dom = srcObj->domain.p->domains) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else if(((val = srcObj->values.vox->values) == NULL))
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  pCnt = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	  for(pIdx =  0; pIdx < pCnt; ++pIdx)
	  {
	    if((errNum == WLZ_ERR_NONE) &&
	       (dom[pIdx].core != NULL) &&
	       (val[pIdx].core != NULL))
	    {
	      int	area2D;
	      double	min2D,
			max2D,
			sum2D,
			sumSq2D;
	      WlzObject	*srcObj2D = NULL;
	      WlzErrorNum errNum2D = WLZ_ERR_NONE;

	      srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom[pIdx], val[pIdx],
				     NULL, NULL, &errNum2D);
	      if(errNum == WLZ_ERR_NONE)
	      {
		area2D = WlzGreyStats2D(srcObj2D, &gType,
					&min2D, &max2D, &sum2D, &sumSq2D,
					&errNum2D);
	      }
	      WlzFreeObj(srcObj2D);
#ifdef _OPENMP
#pragma omp critical
	      {
#endif
		if(errNum2D == WLZ_ERR_NONE)
		{
		  if(area == 0)
		  {
		    min = min2D;
		    max = max2D;
		    area = area2D;
		    sum = sum2D;
		    sumSq = sumSq2D;
		  }
		  else
		  {
		    min = (min < min2D)? min: min2D;
		    max = (max > max2D)? max: max2D;
		    area += area2D;
		    sum += sum2D;
		    sumSq += sumSq2D;
		  }
		}
		else
		{
		  if(errNum == WLZ_ERR_NONE)
		  {
		    errNum = errNum2D;
		  }
		}
#ifdef _OPENMP
	      }
#endif
	    }
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
          ("WlzGreyStats 01 %d %d %g %g %g %g\n",
	  area, (int )gType, min, max, sum, sumSq));
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstGType)
    {
      *dstGType = gType;
    }
    if(dstMin)
    {
      *dstMin = min;
    }
    if(dstMax)
    {
      *dstMax = max;
    }
    if(dstSum)
    {
      *dstSum = sum;
    }
    if(dstSumSq)
    {
      *dstSumSq = sumSq;
    }
    if(dstMean)
    {
      if(area > 0)
      {
        mean = sum / area;
      }
      *dstMean = mean;
    }
    if(dstStdDev)
    {
      if(area > 1)
      {
	stdDev = sqrt((sumSq - (sum * sum / area)) / (area - 1));
      }
      else
      {
	stdDev = 0.0;
      }
      *dstStdDev = stdDev;
    }
  }
  else
  {
    area = -1;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzGreyStats FX %d\n",
	  area));
  return(area);
}
