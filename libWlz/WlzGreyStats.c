#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyStats.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates simple statistics about a Woolz object's
*		grey values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzGreyStats2D						*
* Returns:	int:			Object area or -1 on error.	*
* Purpose:	Calculates simple quick statistics for the given	*
*		2D domain object with grey values.			*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	2D domain object from which to	*
*					calculate the statistics.	*
*		WlzGreyType *dstGType:	Pointer for grey type.		*
*		double *dstMin:		Pointer for minimum value.	*
*		double *dstMax:		Pointer for maximum value.	*
*		double *dstSum:		Pointer for sum of values.	*
*		double *dstSumSq:	Pointer for sum of squares of	*
*					values.				*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number.				*
************************************************************************/
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

/************************************************************************
* Function:	WlzGreyStats						*
* Returns:	int:			Object area or -1 on error.	*
* Purpose:	Calculates simple quick statistics for given 2D or	*
*		3D domain object with grey values.			*
*		Pointers provided for results may be NULL without 	*
*		causing an error.					*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Object from which to calculate	*
*					the statistics.			*
*		WlzGreyType *dstGType:	Pointer for grey type.		*
*		double *dstMin:		Pointer for minimum value.	*
*		double *dstMax:		Pointer for maximum value.	*
*		double *dstSum:		Pointer for sum of values.	*
*		double *dstSumSq:	Pointer for sum of squares of	*
*					values.				*
*		double *dstMean:	Mean value.			*
*		double *dstStdDev:	Standard deviation of values.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number, may be NULL if not	*
*					required.			*
************************************************************************/
int		WlzGreyStats(WlzObject *srcObj,
			     WlzGreyType *dstGType,
			     double *dstMin, double *dstMax,
			     double *dstSum, double *dstSumSq,
			     double *dstMean, double *dstStdDev,
			     WlzErrorNum *dstErr)
{
  int		area = 0,
  		area2D,
		planeIdx,
		planeCount;
  double	min,
		min2D,
		max,
		max2D,
  		mean = -1.0,
		stdDev = -1.0,
		sum = 0.0,
		sum2D,
		sumSq = 0.0,
		sumSq2D;
  WlzGreyType	gType;
  WlzDomain	dummyDom;
  WlzValues	dummyValues;
  WlzDomain	*srcDomains;
  WlzValues	*srcValues;
  WlzObject	*srcObj2D;
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
	dummyDom.core = NULL;
	dummyValues.core = NULL;
	if(srcObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if((srcDomains = srcObj->domain.p->domains) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else if(((srcValues = srcObj->values.vox->values) == NULL))
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
	else
	{
	  srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
	      		         dummyValues, NULL, NULL, &errNum);
	  if((srcObj2D == NULL) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WLZ_ERR_UNSPECIFIED;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  planeIdx =  0;
	  planeCount = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
	  while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	  {
	    srcObj2D->domain = *(srcDomains + planeIdx);
	    srcObj2D->values = *(srcValues + planeIdx);
	    if(srcObj2D->domain.core && srcObj2D->values.core)
	    {
	      area2D = WlzGreyStats2D(srcObj2D, &gType,
				      &min2D, &max2D, &sum2D, &sumSq2D,
				      &errNum);
	      if(errNum == WLZ_ERR_NONE)
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
	    }
	  }
	  srcObj2D->domain.core = NULL;
	  srcObj2D->values.core = NULL;
	  WlzFreeObj(srcObj2D);
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
