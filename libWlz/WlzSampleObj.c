#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzSampleObj.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Subsamples an object through a convolution kernel.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

#define WLZ_SAMPLE_KERNEL_INORM	(0x000100)

static WlzObject *WlzSampleObj2D(WlzObject *, WlzIVertex2,
				WlzSampleFn, WlzIVertex2, WlzErrorNum *),
		*WlzSampleObjIDom(WlzObject *, WlzIVertex2, WlzErrorNum *),
		*WlzSampleObjPoint(WlzObject *, WlzIVertex2, WlzErrorNum *),
		*WlzSampleObjConvI(WlzObject *, int **, WlzIVertex2, int,
				WlzIVertex2, WlzErrorNum *),
		*WlzSampleObjConvD(WlzObject *, double **,
				WlzIVertex2, WlzIVertex2, WlzErrorNum *),
	        *WlzSampleObjRankI(WlzObject *, WlzIVertex2, WlzSampleFn,
				WlzIVertex2, WlzErrorNum *),
		*WlzSampleObjRankD(WlzObject *, WlzIVertex2, WlzSampleFn,
				WlzIVertex2, WlzErrorNum *);
static WlzValues WlzSampleObjConstructRectValues(void **, WlzGreyType,
				WlzIBox2, WlzPixelV, WlzErrorNum *);
static int	WlzSampleObjEstMaxIntervals(WlzDomain, int, int,
			    	WlzIVertex2),
		WlzSampleObjGaussKernelD(double **, WlzIVertex2, WlzIVertex2),
		WlzSampleObjGaussKernelI(int **, WlzIVertex2, int *,
					 WlzIVertex2),
		WlzSampleObjMeanKernelD(double **, WlzIVertex2),
		WlzSampleObjMeanKernelI(int **, WlzIVertex2, int *);

/************************************************************************
* Function:	WlzSampleObj						*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using the given sampling	*
*		factor and sampling method.				*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzSampleFn samFn:	Sampling method.		*
*		WlzErrorNUm *dstErrNum: Destination pointer for error	*
*					number, may be NULL if not	*
*					required.			*
************************************************************************/
WlzObject	*WlzSampleObj(WlzObject *srcObj, WlzIVertex2 samFac,
			      WlzSampleFn samFn,
			      WlzErrorNum *dstErrNum)
{
  WlzObject	*dstObj = NULL;
  WlzIVertex2 	kernelSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzSampleObj FE 0x%lx {%d %d} %d 0x%lx\n",
	   (unsigned long )srcObj, samFac.vtX, samFac.vtY,
	   (int )samFn, (unsigned long )dstErrNum));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((samFac.vtX < 1) || (samFac.vtY < 1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	if(samFac.vtX == 1)
	{
	  kernelSz.vtX = 1;
	}
	else
	{
	  kernelSz.vtX = (samFac.vtX % 2)? samFac.vtX + 2: samFac.vtX + 1;
	}
	if(samFac.vtY == 1)
	{
	  kernelSz.vtY = 1;
	}
	else
	{
	  kernelSz.vtY = (samFac.vtY % 2)? samFac.vtY + 2: samFac.vtY + 1;
	}
	dstObj = WlzSampleObj2D(srcObj, samFac, samFn, kernelSz,
				   &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzSampleObj FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObj2D						*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given 2D domain object using the given	*
*		sampling factor kernel size and sampling method.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzSampleFn samFn:	Sampling method.		*
*		WlzIVertex2 kernelSz:	Size of the convolution kernel.	*
*		WlzErrorNUm *dstErrNum: Destination pointer for error	*
*					number, may be NULL if not	*
*					required.			*
************************************************************************/
static WlzObject *WlzSampleObj2D(WlzObject *srcObj, WlzIVertex2 samFac,
			         WlzSampleFn samFn, WlzIVertex2 kernelSz,
				 WlzErrorNum *dstErrNum)
{
  int		kernelSum,
  		integralGrey;
  WlzGreyType	greyType;
  WlzObject	*dstObj = NULL;
  int		**kernelI;
  double	**kernelD;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj->values.core == NULL)
  {
    dstObj = WlzSampleObjIDom(srcObj, samFac, &errNum);
  }
  else
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type,
					  &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if((samFac.vtX <= 0) || (samFac.vtY <= 0) ||
	 ((samFn != WLZ_SAMPLEFN_POINT) &&
	  ((kernelSz.vtX <= 0) || (kernelSz.vtY <= 0) ||
	   ((kernelSz.vtX % 2) == 0) || ((kernelSz.vtY % 2) == 0))))
      {
	errNum = WLZ_ERR_PARAM_DATA;
      }
      else if((samFn ==  WLZ_SAMPLEFN_POINT) ||
	      ((samFac.vtX == 1) && (samFac.vtY == 1)))
      {
	dstObj = WlzSampleObjPoint(srcObj, samFac, &errNum);
      }
      else
      {
	switch(greyType)
	{
	  case WLZ_GREY_INT:
	  case WLZ_GREY_SHORT:
	  case WLZ_GREY_UBYTE:
	    integralGrey = 1;
	    break;
	  case WLZ_GREY_FLOAT:
	  case WLZ_GREY_DOUBLE:
	    integralGrey = 0;
	    break;
	}
	switch(samFn)
	{
	  case WLZ_SAMPLEFN_GAUSS:
	  case WLZ_SAMPLEFN_MEAN:
	    if(integralGrey)
	    {
	      if(AlcInt2Malloc(&kernelI,
			       kernelSz.vtY, kernelSz.vtX) != ALC_ER_NONE)
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		if(((samFn ==  WLZ_SAMPLEFN_GAUSS) &&
		    WlzSampleObjGaussKernelI(kernelI, kernelSz, 
		    			     &kernelSum, samFac)) ||
		   ((samFn ==  WLZ_SAMPLEFN_MEAN) &&
		    WlzSampleObjMeanKernelI(kernelI, kernelSz, &kernelSum)))
		{
		  dstObj = WlzSampleObjConvI(srcObj, kernelI, kernelSz,
					     kernelSum, samFac, &errNum);
		}
		AlcInt2Free(kernelI);
	      }
	    }
	    else
	    {
	      if(AlcDouble2Malloc(&kernelD,
				  kernelSz.vtY, kernelSz.vtX) != ALC_ER_NONE)
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		if(((samFn ==  WLZ_SAMPLEFN_GAUSS) && 
		    WlzSampleObjGaussKernelD(kernelD, kernelSz, samFac)) ||
		    ((samFn ==  WLZ_SAMPLEFN_MEAN) && 
		     WlzSampleObjMeanKernelD(kernelD, kernelSz)))
		{
		  dstObj = WlzSampleObjConvD(srcObj, kernelD, kernelSz,
					     samFac, &errNum);
		}
		AlcDouble2Free(kernelD);
	      }
	    }
	    break;
	  case WLZ_SAMPLEFN_MIN:
	  case WLZ_SAMPLEFN_MAX:
	  case WLZ_SAMPLEFN_MEDIAN:
	    if(integralGrey)
	    {
	      dstObj = WlzSampleObjRankI(srcObj, samFac, samFn, kernelSz,
	      				 &errNum);
	    }
	    else
	    {
	      dstObj = WlzSampleObjRankD(srcObj, samFac, samFn, kernelSz,
	      				 &errNum);
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
    }
    dstObj = NULL;
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjIDom					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object's interval domain only using	*
*		the given sampling factor.				*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjIDom(WlzObject *srcObj, WlzIVertex2 samFac,
				   WlzErrorNum *dstErrNum)
{
  int		itvCount,
		totItvCount,
  		maxItvCount,
		dstLinePos,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth,
		srcInvWidth;
  WlzObject	*dstObj = NULL;
  WlzIBox2	srcBox,
  		dstBox;
  WlzInterval	*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzIntervalWSpace srcIWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  srcDom = srcObj->domain;
  srcBox.xMin = srcDom.i->kol1;
  srcBox.yMin = srcDom.i->line1;
  srcBox.xMax = srcDom.i->lastkl;
  srcBox.yMax = srcDom.i->lastln;
  dstBox.xMin = (srcBox.xMin < 0) ?
		(srcBox.xMin  - samFac.vtX + 1) / samFac.vtX :
		(srcBox.xMin  + samFac.vtX - 1) / samFac.vtX;
  dstBox.yMin = (srcBox.yMin < 0) ?
		(srcBox.yMin - samFac.vtY + 1) / samFac.vtY :
		(srcBox.yMin + samFac.vtY - 1) / samFac.vtY;
  dstBox.xMax = srcBox.xMax / samFac.vtX;
  dstBox.yMax = srcBox.yMax / samFac.vtY;
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItv0 = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					      sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
	if(dstItv0)
	{
	  AlcFree(dstItv0);
	}
	if(errNum == WLZ_ERR_NONE)
	{
          errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      else
      {
        dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr, (void *)dstItv0,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstDom.core)
  {
    totItvCount = itvCount = 0;
    dstItv1 = dstItv0;
    errNum = WlzInitRasterScan(srcObj, &srcIWsp, WLZ_RASTERDIR_ILIC);
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextInterval(&srcIWsp)) == WLZ_ERR_NONE))
    {
      if((srcIWsp.linpos % samFac.vtY) == 0)
      {
	dstInvLeftPos = (srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX;
	dstInvRgtPos = srcIWsp.rgtpos / samFac.vtX;
	srcInvWidth = srcIWsp.rgtpos - srcIWsp.lftpos + 1;
	dstInvWidth = (srcInvWidth >= (srcIWsp.rgtpos % samFac.vtX)) ?
		      dstInvRgtPos - dstInvLeftPos + 1 : 0;
	if(dstInvWidth > 0)
	{
	  dstLinePos = srcIWsp.linpos / samFac.vtY;
	  dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	  dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	  ++dstItv1;
	  totItvCount += ++itvCount;
	  if(totItvCount >= maxItvCount)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (itvCount > 0) && (srcIWsp.intrmn == 0))
	{
	  WlzMakeInterval(dstLinePos, dstDom.i, itvCount, dstItv0);
	  dstItv0 = dstItv1;
	  itvCount = 0;
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzStandardIntervalDomain(dstDom.i);
      dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, srcObj->values,
			NULL, NULL, &errNum);
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjPoint					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using a simple point sampling	*
*		method and the given sampling factor.			*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjPoint(WlzObject *srcObj, WlzIVertex2 samFac,
				    WlzErrorNum *dstErrNum)
{
  int		tI0,
		itvCount,
		totItvCount,
  		maxItvCount,
		dstWidth,
		dstOffset,
		srcOffset,
		dstLinePos,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth,
		srcInvWidth;
  WlzGreyType	greyType;
  WlzObject	*dstObj = NULL;
  WlzIBox2	srcBox,
  		dstBox;
  void		*dstGreyValues = NULL;
  WlzInterval	*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzValues	dstValues;
  WlzGreyP	srcPix,
		dstPix;
  WlzPixelV	backgroundPix;
  WlzIntervalWSpace srcIWsp;
  WlzGreyWSpace	srcGWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstValues.core = NULL;
  srcDom = srcObj->domain;
  srcBox.xMin = srcDom.i->kol1;
  srcBox.yMin = srcDom.i->line1;
  srcBox.xMax = srcDom.i->lastkl;
  srcBox.yMax = srcDom.i->lastln;
  dstBox.xMin = (srcBox.xMin < 0) ?
		(srcBox.xMin  - samFac.vtX + 1) / samFac.vtX :
		(srcBox.xMin  + samFac.vtX - 1) / samFac.vtX;
  dstBox.yMin = (srcBox.yMin < 0) ?
		(srcBox.yMin - samFac.vtY + 1) / samFac.vtY :
		(srcBox.yMin + samFac.vtY - 1) / samFac.vtY;
  dstBox.xMax = srcBox.xMax / samFac.vtX;
  dstBox.yMax = srcBox.yMax / samFac.vtY;
  dstWidth = dstBox.xMax - dstBox.xMin + 1;
  backgroundPix = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues = WlzSampleObjConstructRectValues(&dstGreyValues, greyType,
						dstBox,
						backgroundPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItv0 = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					      sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
	if(dstItv0)
	{
	  AlcFree(dstItv0);
	}
	if(errNum == WLZ_ERR_NONE)
	{
          errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      else
      {
        dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr, (void *)dstItv0,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core)
  {
    totItvCount = itvCount = 0;
    dstItv1 = dstItv0;
    errNum = WlzInitGreyScan(srcObj, &srcIWsp, &srcGWsp);
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&srcIWsp)) == WLZ_ERR_NONE))
    {
      if((srcIWsp.linpos % samFac.vtY) == 0)
      {
	dstInvLeftPos = (srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX;
	dstInvRgtPos = srcIWsp.rgtpos / samFac.vtX;
	srcInvWidth = srcIWsp.rgtpos - srcIWsp.lftpos + 1;
	dstInvWidth = (srcInvWidth >= (srcIWsp.rgtpos % samFac.vtX)) ?
		      dstInvRgtPos - dstInvLeftPos + 1 : 0;
	if(dstInvWidth > 0)
	{
	  dstLinePos = srcIWsp.linpos / samFac.vtY;
	  dstOffset = ((dstLinePos - dstBox.yMin) * dstWidth) +
	  	      dstInvLeftPos - dstBox.xMin;
	  srcOffset = srcIWsp.lftpos % samFac.vtX;
	  dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	  dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	  ++dstItv1;
	  totItvCount += ++itvCount;
	  if(totItvCount >= maxItvCount)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    switch(greyType)
	    {
	      case WLZ_GREY_INT:
		srcPix.inp = srcGWsp.u_grintptr.inp + srcOffset;
		dstPix.inp = (int *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.inp)++ = *(srcPix.inp);
		  srcPix.inp += samFac.vtX;
		}
		break;
	      case WLZ_GREY_SHORT:
		srcPix.shp = srcGWsp.u_grintptr.shp + srcOffset;
		dstPix.shp = (short *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.shp)++ = *(srcPix.shp);
		  srcPix.shp += samFac.vtX;
		}
		break;
	      case WLZ_GREY_UBYTE:
		srcPix.ubp = srcGWsp.u_grintptr.ubp + srcOffset;
		dstPix.ubp = (UBYTE *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.ubp)++ = *(srcPix.ubp);
		  srcPix.ubp += samFac.vtX;
		}
		break;
	      case WLZ_GREY_FLOAT:
		srcPix.flp = srcGWsp.u_grintptr.flp + srcOffset;
		dstPix.flp = (float *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.flp)++ = *(srcPix.flp);
		  srcPix.flp += samFac.vtX;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		srcPix.dbp = srcGWsp.u_grintptr.dbp + srcOffset;
		dstPix.dbp = (double *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.dbp)++ = *(srcPix.dbp);
		  srcPix.dbp += samFac.vtX;
		}
		break;
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (itvCount > 0) && (srcIWsp.intrmn == 0))
	{
	  WlzMakeInterval(dstLinePos, dstDom.i, itvCount, dstItv0);
	  dstItv0 = dstItv1;
	  itvCount = 0;
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzStandardIntervalDomain(dstDom.i);
      dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstValues,
			NULL, NULL, &errNum);
      if(dstObj)
      {
	(void )WlzSetBackground(dstObj, WlzGetBackground(srcObj, NULL));
      }
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjConvI					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using the given integer	*
*		convolution kernel and the given sampling factor.	*
*		The domain of the the new object is smaller than the	*
*		source because of the convolution.			*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		int **kernel:		Integer convolution kernel.	*
*		WlzIVertex2 kernelSz:	Size of the convolution kernel.	*
*		int kernelSum:		Kernel sum for normalization.	*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjConvI(WlzObject *srcObj, int **kernel,
				    WlzIVertex2 kernelSz, int kernelSum,
				    WlzIVertex2 samFac, WlzErrorNum *dstErrNum)
{
  int		tI0,
  		tI1,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		itvCount,
		totItvCount,
  		maxItvCount,
  		srcWidth,
		dstWidth,
		dstOffset,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth,
		backgroundVal;
  WlzGreyType	greyType;
  int		**bufData = NULL;
  int		*tIP0,
		*tIP1,
  		*bufMarks = NULL;
  WlzObject	*dstObj = NULL;
  WlzIVertex2	bufPos,
		dstPos;
  WlzIBox2	dstBox;
  void		*dstGreyValues = NULL;
  WlzInterval	*dstItvBase = NULL,
  		*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzValues	dstValues;
  WlzGreyP	tGP0;
  WlzPixelV	backgroundPix;
  WlzIntervalWSpace bufIWsp,
		srcIWsp;
  WlzGreyWSpace	bufGWsp,
		srcGWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstValues.core = NULL;
  srcDom = srcObj->domain;
  dstBox.xMin = (srcDom.i->kol1 + (kernelSz.vtX / 2) + samFac.vtX - 1) / 
		samFac.vtX;
  dstBox.yMin = (srcDom.i->line1 + (kernelSz.vtY / 2) + samFac.vtY - 1) /
		samFac.vtY;
  dstBox.xMax = (srcDom.i->lastkl - (kernelSz.vtX / 2)) / samFac.vtX;
  dstBox.yMax = (srcDom.i->lastln - (kernelSz.vtY / 2)) / samFac.vtY;
  dstWidth = dstBox.xMax - dstBox.xMin + 1;
  srcWidth = srcDom.i->lastkl - srcDom.i->kol1 + 1;
  backgroundPix = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues = WlzSampleObjConstructRectValues(&dstGreyValues, greyType,
						dstBox,
						backgroundPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					       sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
        if(dstItvBase)
	{
	  AlcFree(dstItvBase);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
        }
      }
      else
      {
	dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr,
					   (void *)dstItvBase,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core)
  {
    if(AlcInt2Malloc(&bufData, kernelSz.vtY, srcWidth) == ALC_ER_NONE)
    {
      if((bufMarks = (int *)AlcMalloc((unsigned long )(kernelSz.vtY) *
      				      sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	tIP0 = bufMarks;
	tI0 = kernelSz.vtY;
	tI1 = srcDom.i->line1 - (2 * kernelSz.vtY); 	 /* Use invalid line */
	while(tI0-- > 0)        /* Mark buffer lines stale with invalid line */
	{
	  *tIP0++ = tI1;
	}
	switch(backgroundPix.type)
	{
	  case WLZ_GREY_INT:
	    backgroundVal = backgroundPix.v.inv;
	    break;
	  case WLZ_GREY_SHORT:
	    backgroundVal = (int )(backgroundPix.v.shv);
	    break;
	  case WLZ_GREY_UBYTE:
	    backgroundVal = (int )(backgroundPix.v.ubv);
	    break;
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core &&
     bufData && bufMarks)
  {
    totItvCount = itvCount = 0;
    dstItv0 = dstItvBase;
    dstItv1 = dstItvBase;
    bufBase = WLZ_ABS(srcDom.i->line1) + kernelSz.vtY;  /* Make sure not -ve */
    if(((errNum = WlzInitGreyScan(srcObj, &srcIWsp,
    				  &srcGWsp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(srcObj, &bufIWsp,
    				  &bufGWsp)) == WLZ_ERR_NONE))
    {
      bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
    }
    while((errNum == WLZ_ERR_NONE) && (WlzNextGreyInterval(&srcIWsp) == 0))
    {
      if(((srcIWsp.linpos % samFac.vtY) == 0) &&
      	 ((tI0 = (srcIWsp.linpos / samFac.vtY)) >= dstBox.yMin) &&
	 (tI0 <= dstBox.yMax))
      {
	dstInvLeftPos = WLZ_MAX((srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX,
				dstBox.xMin);
	dstInvRgtPos = WLZ_MIN(srcIWsp.rgtpos / samFac.vtX, dstBox.xMax);
	dstInvWidth = dstInvRgtPos - dstInvLeftPos + 1;
	dstPos.vtY = srcIWsp.linpos / samFac.vtY;
	dstOffset = ((dstPos.vtY - dstBox.yMin) * dstWidth) +
		    dstInvLeftPos - dstBox.xMin;
	dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	++dstItv1;
	totItvCount += ++itvCount;
	if(totItvCount >= maxItvCount)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  for(idY = 0; idY < kernelSz.vtY; ++idY)    /* Check and update buf */
	  {
	    bufPos.vtX = srcDom.i->kol1;
	    idB = (bufBase + bufPos.vtY) % kernelSz.vtY;
	    if(*(bufMarks + idB) != bufPos.vtY)  /* Check if update required */
	    {
	      *(bufMarks + idB) = bufPos.vtY;
	      if((bufPos.vtY < srcDom.i->line1) ||
		 (bufPos.vtY > srcDom.i->lastln))
	      {
		WlzValueSetInt(*(bufData + idB), backgroundVal, srcWidth);
	      }
	      else
	      {
		while((bufIWsp.linpos <= bufPos.vtY) && (bufIwspFlag == 0))
		{
		  if(bufIWsp.linpos == bufPos.vtY)
		  {
		    if((tI1 = bufIWsp.lftpos - bufPos.vtX) > 0)
		    {
		      tIP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		      WlzValueSetInt(tIP0, backgroundVal, tI1);
		    }
		    tI1 = bufIWsp.rgtpos - bufIWsp.lftpos + 1;
		    tIP0 = *(bufData + idB) + bufIWsp.lftpos - srcDom.i->kol1;
		    switch(greyType)
		    {
		      case WLZ_GREY_INT:
			WlzValueCopyIntToInt(tIP0, bufGWsp.u_grintptr.inp,
					     tI1);
			break;
		      case WLZ_GREY_SHORT:
			WlzValueCopyShortToInt(tIP0, bufGWsp.u_grintptr.shp,
					       tI1);
			break;
		      case WLZ_GREY_UBYTE:
			WlzValueCopyUByteToInt(tIP0, bufGWsp.u_grintptr.ubp,
					       tI1);
			break;
		    }
		    bufPos.vtX = bufIWsp.rgtpos + 1;
		  }
		  bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
		}
		if((tI1 = srcDom.i->lastkl - bufPos.vtX) > 0)
		{
		  tIP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  WlzValueSetInt(tIP0, backgroundVal, tI1);
		}
	      }
	    }
	    ++(bufPos.vtY);
	  }
	  bufPos.vtX = dstInvLeftPos *
		       samFac.vtX;   /* Sample and convolve through interval */
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  switch(greyType)
	  {
	    case WLZ_GREY_INT:
	      tGP0.inp = (int *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_SHORT:
	      tGP0.shp = (short *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_UBYTE:
	      tGP0.ubp = (UBYTE *)dstGreyValues + dstOffset;
	      break;
	  }
	  tI0 = dstInvWidth;
	  while(tI0-- > 0)
	  {
	    tI1 = 0;
	    for(idY = 0; idY < kernelSz.vtY; ++idY)
	    {
	      tIP0 = *(kernel + idY);
	      idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
	      tIP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
	      idX = kernelSz.vtX / 2;
	      tI1 += *tIP0++ * *tIP1++;
	      while(idX-- > 0)
	      {
		tI1 += (*tIP0++ * *tIP1++);
		tI1 += (*tIP0++ * *tIP1++);
	      }
	    }
	    tI1 /= kernelSum;
	    switch(greyType)
	    {
	      case WLZ_GREY_INT:
		*(tGP0.inp)++ = tI1;
		break;
	      case WLZ_GREY_SHORT:
		*(tGP0.shp)++ = tI1;
		break;
	      case WLZ_GREY_UBYTE:
		if(tI1 < 0)
		{
		  tI1 = 0;
		}
		else if(tI1 > 255)
		{
		  tI1 = 255;
		}
		*(tGP0.ubp)++ = (UBYTE )tI1;
		break;
	    }
	    bufPos.vtX += samFac.vtX;
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     (itvCount > 0) && (srcIWsp.intrmn == 0))
	  {
	    WlzMakeInterval(dstPos.vtY, dstDom.i, itvCount, dstItv0);
	    dstItv0 = dstItv1;
	    itvCount = 0;
	  }
	}
      }
    }
    (void )WlzStandardIntervalDomain(dstDom.i);
    dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstValues, NULL, NULL,
    			 &errNum);
    if(dstObj)
    {
      (void )WlzSetBackground(dstObj, WlzGetBackground(srcObj, NULL));
    }
  }
  if(bufData) 					      /* Free up buffer data */
  {
    AlcInt2Free(bufData);
  }
  if(bufMarks)
  {
    AlcFree(bufMarks);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(dstObj == NULL)				        /* Clear up on error */
  {
    (void )WlzFreeValues(dstValues);
    if(dstDom.core)
    {
      (void )WlzFreeDomain(dstDom);
    }
    else if(dstItvBase)
    {
      AlcFree(dstItvBase);
    }
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjConvD					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using the given double		*
*		convolution kernel and the given sampling factor.	*
*		The domain of the the new object is smaller than the	*
*		source because of the convolution.			*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		double **kernel:	Double convolution kernel.	*
*		WlzIVertex2 kernelSz:	Size of the convolution kernel.	*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjConvD(WlzObject *srcObj, double **kernel,
				    WlzIVertex2 kernelSz, WlzIVertex2 samFac,
				    WlzErrorNum *dstErrNum)
{
  int		tI0,
  		tI1,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		itvCount,
		totItvCount,
  		maxItvCount,
  		srcWidth,
		dstWidth,
		dstOffset,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth;
  WlzGreyType 	greyType;
  double	tD0,
		backgroundVal;
  double	**bufData = NULL;
  double	*tDP0,
  		*tDP1;
  int		*tIP0,
  		*bufMarks = NULL;
  WlzObject	*dstObj = NULL;
  WlzIVertex2	bufPos,
		dstPos;
  WlzIBox2	dstBox;
  void		*dstGreyValues = NULL;
  WlzInterval	*dstItvBase = NULL,
  		*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzValues	dstValues;
  WlzGreyP	tGP0;
  WlzPixelV	backgroundPix;
  WlzIntervalWSpace bufIWsp,
		srcIWsp;
  WlzGreyWSpace	bufGWsp,
		srcGWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstValues.core = NULL;
  srcDom = srcObj->domain;
  dstBox.xMin = (srcDom.i->kol1 + (kernelSz.vtX / 2) + samFac.vtX - 1) / 
		samFac.vtX;
  dstBox.yMin = (srcDom.i->line1 + (kernelSz.vtY / 2) + samFac.vtY - 1) /
		samFac.vtY;
  dstBox.xMax = (srcDom.i->lastkl - (kernelSz.vtX / 2)) / samFac.vtX;
  dstBox.yMax = (srcDom.i->lastln - (kernelSz.vtY / 2)) / samFac.vtY;
  dstWidth = dstBox.xMax - dstBox.xMin + 1;
  srcWidth = srcDom.i->lastkl - srcDom.i->kol1 + 1;
  backgroundPix = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues = WlzSampleObjConstructRectValues(&dstGreyValues, greyType,
						dstBox,
						backgroundPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					       sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
        if(dstItvBase)
	{
	  AlcFree(dstItvBase);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      else
      {
	dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr,
					   (void *)dstItvBase,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core)
  {
    if(AlcDouble2Malloc(&bufData, kernelSz.vtY, srcWidth) == ALC_ER_NONE)
    {
      if((bufMarks = (int *)AlcMalloc((unsigned long )kernelSz.vtY *
      				      sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	tIP0 = bufMarks;
	tI0 = kernelSz.vtY;
	tI1 = srcDom.i->line1 - (2 * kernelSz.vtY); 	 /* Use invalid line */
	while(tI0-- > 0)        /* Mark buffer lines stale with invalid line */
	{
	  *tIP0++ = tI1;
	}
	switch(backgroundPix.type)
	{
	  case WLZ_GREY_FLOAT:
	    backgroundVal = backgroundPix.v.flv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    backgroundVal = backgroundPix.v.dbv;
	    break;
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core &&
     bufData && bufMarks)
  {
    totItvCount = itvCount = 0;
    dstItv0 = dstItvBase;
    dstItv1 = dstItvBase;
    bufBase = WLZ_ABS(srcDom.i->line1) + kernelSz.vtY;  /* Make sure never -ve */
    if(((errNum = WlzInitGreyScan(srcObj, &srcIWsp,
    				  &srcGWsp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(srcObj, &bufIWsp,
    				  &bufGWsp)) == WLZ_ERR_NONE))
    {
      bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
    }
    while((errNum == WLZ_ERR_NONE) && (WlzNextGreyInterval(&srcIWsp) == 0))
    {
      if(((srcIWsp.linpos % samFac.vtY) == 0) &&
      	 ((tI0 = (srcIWsp.linpos / samFac.vtY)) >= dstBox.yMin) &&
	 (tI0 <= dstBox.yMax))
      {
	dstInvLeftPos = WLZ_MAX((srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX,
				dstBox.xMin);
	dstInvRgtPos = WLZ_MIN(srcIWsp.rgtpos / samFac.vtX, dstBox.xMax);
	dstInvWidth = dstInvRgtPos - dstInvLeftPos + 1;
	dstPos.vtY = srcIWsp.linpos / samFac.vtY;
	dstOffset = ((dstPos.vtY - dstBox.yMin) * dstWidth) +
		    dstInvLeftPos - dstBox.xMin;
	dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	++dstItv1;
	totItvCount += ++itvCount;
	if(totItvCount >= maxItvCount)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	{
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  for(idY = 0; idY < kernelSz.vtY; ++idY) /* Check and update buffer */
	  {
	    bufPos.vtX = srcDom.i->kol1;
	    idB = (bufBase + bufPos.vtY) % kernelSz.vtY;
	    if(*(bufMarks + idB) != bufPos.vtY)  /* Check if update required */
	    {
	      *(bufMarks + idB) = bufPos.vtY;
	      if((bufPos.vtY < srcDom.i->line1) ||
		 (bufPos.vtY > srcDom.i->lastln))
	      {
		WlzValueSetDouble(*(bufData + idB), backgroundVal, srcWidth);
	      }
	      else
	      {
		while((bufIWsp.linpos <= bufPos.vtY) && (bufIwspFlag == 0))
		{
		  if(bufIWsp.linpos == bufPos.vtY)
		  {
		    if((tI1 = bufIWsp.lftpos - bufPos.vtX) > 0)
		    {
		      tDP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		      WlzValueSetDouble(tDP0, backgroundVal, tI1);
		    }
		    tI1 = bufIWsp.rgtpos - bufIWsp.lftpos + 1;
		    tDP0 = *(bufData + idB) + bufIWsp.lftpos - srcDom.i->kol1;
		    switch(greyType)
		    {
		      case WLZ_GREY_FLOAT:
			WlzValueCopyFloatToDouble(tDP0, bufGWsp.u_grintptr.flp,
						  tI1);
			break;
		      case WLZ_GREY_DOUBLE:
			WlzValueCopyDoubleToDouble(tDP0, bufGWsp.u_grintptr.dbp,
						   tI1);
			break;
		    }
		    bufPos.vtX = bufIWsp.rgtpos + 1;
		  }
		  bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
		}
		if((tI1 = srcDom.i->lastkl - bufPos.vtX) > 0)
		{
		  tDP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  WlzValueSetDouble(tDP0, backgroundVal, tI1);
		}
	      }
	    }
	    ++(bufPos.vtY);
	  }
	  bufPos.vtX = dstInvLeftPos *
		       samFac.vtX;   /* Sample and convolve through interval */
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  switch(greyType)
	  {
	    case WLZ_GREY_FLOAT:
	      tGP0.flp = (float *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_DOUBLE:
	      tGP0.dbp = (double *)dstGreyValues + dstOffset;
	      break;
	  }
	  tI0 = dstInvWidth;
	  while(tI0-- > 0)
	  {
	    tD0 = 0;
	    for(idY = 0; idY < kernelSz.vtY; ++idY)
	    {
	      tDP0 = *(kernel + idY);
	      idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
	      tDP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
	      idX = kernelSz.vtX / 2;
	      tD0 += *tDP0++ * *tDP1++;
	      while(idX-- > 0)
	      {
		tD0 += (*tDP0++ * *tDP1++);
		tD0 += (*tDP0++ * *tDP1++);
	      }
	    }
	    switch(greyType)
	    {
	      case WLZ_GREY_FLOAT:
		*(tGP0.flp)++ = tD0;
		break;
	      case WLZ_GREY_DOUBLE:
		*(tGP0.dbp)++ = tD0;
		break;
	    }
	    bufPos.vtX += samFac.vtX;
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     (itvCount > 0) && (srcIWsp.intrmn == 0))
	  {
	    WlzMakeInterval(dstPos.vtY, dstDom.i, itvCount, dstItv0);
	    dstItv0 = dstItv1;
	    itvCount = 0;
	  }
	}
      }
    }
    (void )WlzStandardIntervalDomain(dstDom.i);
    dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstValues, NULL, NULL,
    			 &errNum);
    if(dstObj)
    {
      (void )WlzSetBackground(dstObj, WlzGetBackground(srcObj, NULL));
    }
  }
  if(bufData) 					      /* Free up buffer data */
  {
    AlcDouble2Free(bufData);
  }
  if(bufMarks)
  {
    AlcFree(bufMarks);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(dstObj == NULL)				        /* Clear up on error */
  {
    (void )WlzFreeValues(dstValues);
    if(dstDom.core)
    {
      (void )WlzFreeDomain(dstDom);
    }
    else if(dstItvBase)
    {
      AlcFree(dstItvBase);
    }
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjRankI					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using the given (rank) 	*
*		sampling method and sampling factor.			*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzSampleFn samFn:	Rank sampling method, eg min,	*
*					max, median,....		*
*		WlzIVertex2 kernelSz:	Size of the convolution kernel.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjRankI(WlzObject *srcObj, WlzIVertex2 samFac,
				    WlzSampleFn samFn, WlzIVertex2 kernelSz,
				    WlzErrorNum *dstErrNum)
{
  int		tI0,
  		tI1,
		tI2,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		bufMedianSz,
		itvCount,
		totItvCount,
  		maxItvCount,
  		srcWidth,
		dstWidth,
		dstOffset,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth,
		backgroundVal;
  WlzGreyType	greyType;
  int		**bufData = NULL;
  int		*tIP0,
		*tIP1,
  		*bufMarks = NULL,
		*bufMedian = NULL;
  WlzObject	*dstObj = NULL;
  WlzIVertex2	bufPos,
		dstPos;
  WlzIBox2	dstBox;
  void		*dstGreyValues = NULL;
  WlzInterval	*dstItvBase = NULL,
  		*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzValues	dstValues;
  WlzGreyP	tGP0;
  WlzPixelV	backgroundPix;
  WlzIntervalWSpace bufIWsp,
		srcIWsp;
  WlzGreyWSpace	bufGWsp,
		srcGWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstValues.core = NULL;
  srcDom = srcObj->domain;

  dstBox.xMin = (srcDom.i->kol1 + (kernelSz.vtX / 2) + samFac.vtX - 1) / 
		samFac.vtX;
  dstBox.yMin = (srcDom.i->line1 + (kernelSz.vtY / 2) + samFac.vtY - 1) /
		samFac.vtY;
  dstBox.xMax = (srcDom.i->lastkl - (kernelSz.vtX / 2)) / samFac.vtX;
  dstBox.yMax = (srcDom.i->lastln - (kernelSz.vtY / 2)) / samFac.vtY;
  dstWidth = dstBox.xMax - dstBox.xMin + 1;
  srcWidth = srcDom.i->lastkl - srcDom.i->kol1 + 1;
  backgroundPix = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues = WlzSampleObjConstructRectValues(&dstGreyValues, greyType,
						dstBox,
						backgroundPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					       sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
        if(dstItvBase)
	{
	  AlcFree(dstItvBase);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
        }
      }
      else
      {
	dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr,
					   (void *)dstItvBase,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core)
  {
    if(AlcInt2Malloc(&bufData, kernelSz.vtY, srcWidth) == ALC_ER_NONE)
    {
      if((bufMarks = (int *)AlcMalloc((unsigned long )(kernelSz.vtY) *
      				      sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	tIP0 = bufMarks;
	tI0 = kernelSz.vtY;
	tI1 = srcDom.i->line1 - (2 * kernelSz.vtY); 	 /* Use invalid line */
	while(tI0-- > 0)        /* Mark buffer lines stale with invalid line */
	{
	  *tIP0++ = tI1;
	}
	switch(backgroundPix.type)
	{
	  case WLZ_GREY_INT:
	    backgroundVal = backgroundPix.v.inv;
	    break;
	  case WLZ_GREY_SHORT:
	    backgroundVal = (int )(backgroundPix.v.shv);
	    break;
	  case WLZ_GREY_UBYTE:
	    backgroundVal = (int )(backgroundPix.v.ubv);
	    break;
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (samFn == WLZ_SAMPLEFN_MEDIAN))
  {
    bufMedianSz = kernelSz.vtX * kernelSz.vtY;
    if((bufMedian = (int *)AlcMalloc((unsigned long )bufMedianSz *
    				     sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core &&
     bufData && bufMarks)
  {
    totItvCount = itvCount = 0;
    dstItv0 = dstItvBase;
    dstItv1 = dstItvBase;
    bufBase = WLZ_ABS(srcDom.i->line1) + kernelSz.vtY;  /* Make sure not -ve */
    if(((errNum = WlzInitGreyScan(srcObj, &srcIWsp,
    				  &srcGWsp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(srcObj, &bufIWsp,
    				  &bufGWsp)) == WLZ_ERR_NONE))
    {
      bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
    }
    while((errNum == WLZ_ERR_NONE) && (WlzNextGreyInterval(&srcIWsp) == 0))
    {
      if(((srcIWsp.linpos % samFac.vtY) == 0) &&
      	 ((tI0 = (srcIWsp.linpos / samFac.vtY)) >= dstBox.yMin) &&
	 (tI0 <= dstBox.yMax))
      {
	dstInvLeftPos = WLZ_MAX((srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX,
				dstBox.xMin);
	dstInvRgtPos = WLZ_MIN(srcIWsp.rgtpos / samFac.vtX, dstBox.xMax);
	dstInvWidth = dstInvRgtPos - dstInvLeftPos + 1;
	dstPos.vtY = srcIWsp.linpos / samFac.vtY;
	dstOffset = ((dstPos.vtY - dstBox.yMin) * dstWidth) +
		    dstInvLeftPos - dstBox.xMin;
	dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	++dstItv1;
	totItvCount += ++itvCount;
	if(totItvCount >= maxItvCount)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	{
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  for(idY = 0; idY < kernelSz.vtY; ++idY) /* Check and update buffer */
	  {
	    bufPos.vtX = srcDom.i->kol1;
	    idB = (bufBase + bufPos.vtY) % kernelSz.vtY;
	    if(*(bufMarks + idB) != bufPos.vtY)  /* Check if update required */
	    {
	      *(bufMarks + idB) = bufPos.vtY;
	      if((bufPos.vtY < srcDom.i->line1) ||
		 (bufPos.vtY > srcDom.i->lastln))
	      {
		WlzValueSetInt(*(bufData + idB), backgroundVal, srcWidth);
	      }
	      else
	      {
		while((bufIWsp.linpos <= bufPos.vtY) && (bufIwspFlag == 0))
		{
		  if(bufIWsp.linpos == bufPos.vtY)
		  {
		    if((tI1 = bufIWsp.lftpos - bufPos.vtX) > 0)
		    {
		      tIP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		      WlzValueSetInt(tIP0, backgroundVal, tI1);
		    }
		    tI1 = bufIWsp.rgtpos - bufIWsp.lftpos + 1;
		    tIP0 = *(bufData + idB) + bufIWsp.lftpos - srcDom.i->kol1;
		    switch(greyType)
		    {
		      case WLZ_GREY_INT:
			WlzValueCopyIntToInt(tIP0, bufGWsp.u_grintptr.inp,
					     tI1);
			break;
		      case WLZ_GREY_SHORT:
			WlzValueCopyShortToInt(tIP0, bufGWsp.u_grintptr.shp,
					       tI1);
			break;
		      case WLZ_GREY_UBYTE:
			WlzValueCopyUByteToInt(tIP0, bufGWsp.u_grintptr.ubp,
					       tI1);
			break;
		    }
		    bufPos.vtX = bufIWsp.rgtpos + 1;
		  }
		  bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
		}
		if((tI1 = srcDom.i->lastkl - bufPos.vtX) > 0)
		{
		  tIP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  WlzValueSetInt(tIP0, backgroundVal, tI1);
		}
	      }
	    }
	    ++(bufPos.vtY);
	  }
	  /* Sample by rank through the interval */
	  bufPos.vtX = dstInvLeftPos * samFac.vtX;
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  switch(greyType)
	  {
	    case WLZ_GREY_INT:
	      tGP0.inp = (int *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_SHORT:
	      tGP0.shp = (short *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_UBYTE:
	      tGP0.ubp = (UBYTE *)dstGreyValues + dstOffset;
	      break;
	  }
	  for(tI0 = 0; tI0 < dstInvWidth; ++tI0)
	  {
	    tI2 = 0;
	    switch(samFn)
	    {
	      case WLZ_SAMPLEFN_MIN:
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tIP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  if(tI2)
		  {
		    if(*tIP1 < tI1)
		    {
		      tI1 = *tIP1;
		    }
		  }
		  else
		  {
		    tI1 = *tIP1;
		    tI2 = 1;
		  }
		  while(idX-- > 0)
		  {
		    if(*++tIP1 < tI1)
		    {
		      tI1 = *tIP1;
		    }
		    if(*++tIP1 < tI1)
		    {
		      tI1 = *tIP1;
		    }
		  }
		}
		break;
	      case WLZ_SAMPLEFN_MAX:
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tIP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  if(tI2)
		  {
		    if(*tIP1 > tI1)
		    {
		      tI1 = *tIP1;
		    }
		  }
		  else
		  {
		    tI1 = *tIP1;
		    tI2 = 1;
		  }
		  while(idX-- > 0)
		  {
		    if(*++tIP1 > tI1)
		    {
		      tI1 = *tIP1;
		    }
		    if(*++tIP1 > tI1)
		    {
		      tI1 = *tIP1;
		    }
		  }
		}
		break;
	      case WLZ_SAMPLEFN_MEDIAN:
		tIP0 = bufMedian;
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tIP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  /* First copy to the median buffer. */
		  *tIP0++ = *tIP1++;
		  while(idX-- > 0)
		  {
		    *tIP0++ = *tIP1++;
		    *tIP0++ = *tIP1++;
		  }
		  /* Then find median. */
		  tI1 = WlzValueMedianInt(bufMedian, bufMedianSz);
		}
		break;
	    }
	    switch(greyType)
	    {
	      case WLZ_GREY_INT:
		*(tGP0.inp)++ = tI1;
		break;
	      case WLZ_GREY_SHORT:
		*(tGP0.shp)++ = tI1;
		break;
	      case WLZ_GREY_UBYTE:
		*(tGP0.ubp)++ = (UBYTE )tI1; /* Fits because rank operation. */
		break;
	    }
	    bufPos.vtX += samFac.vtX;
	  }

	  if((errNum == WLZ_ERR_NONE) &&
	     (itvCount > 0) && (srcIWsp.intrmn == 0))
	  {
	    WlzMakeInterval(dstPos.vtY, dstDom.i, itvCount, dstItv0);
	    dstItv0 = dstItv1;
	    itvCount = 0;
	  }
	}
      }
    }
    (void )WlzStandardIntervalDomain(dstDom.i);
    dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstValues, NULL, NULL,
    			 &errNum);
    if(dstObj)
    {
      (void )WlzSetBackground(dstObj, WlzGetBackground(srcObj, NULL));
    }
  }
  if(bufData) 					      /* Free up buffer data */
  {
    AlcInt2Free(bufData);
  }
  if(bufMarks)
  {
    AlcFree(bufMarks);
  }
  if(bufMedian)
  {
    AlcFree(bufMedian);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(dstObj == NULL)				        /* Clear up on error */
  {
    (void )WlzFreeValues(dstValues);
    if(dstDom.core)
    {
      (void )WlzFreeDomain(dstDom);
    }
    else if(dstItvBase)
    {
      AlcFree(dstItvBase);
    }
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjRankD					*
* Returns:	WlzObject *:		New sampled object.		*
* Purpose:	Samples the given object using the given (rank) 	*
*		sampling method and sampling factor.			*
*		This function assumes it's parameters to be valid.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given source object.		*
*		WlzIVertex2 samFac:	Sampling factor for both rows	*
*					and columns. Every pixel == 1,	*
*					every other pixel == 2, ....	*
*		WlzSampleFn samFn:	Rank sampling method, eg min,	*
*					max, median,....		*
*		WlzIVertex2 kernelSz:	Size of the convolution kernel.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number (may be NULL).		*
************************************************************************/
static WlzObject *WlzSampleObjRankD(WlzObject *srcObj, WlzIVertex2 samFac,
				    WlzSampleFn samFn, WlzIVertex2 kernelSz,
				    WlzErrorNum *dstErrNum)
{
  int		tI0,
  		tI1,
		tI2,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		bufMedianSz,
		itvCount,
		totItvCount,
  		maxItvCount,
  		srcWidth,
		dstWidth,
		dstOffset,
		dstInvLeftPos,
		dstInvRgtPos,
		dstInvWidth;
  double	tD0,
		backgroundVal;
  WlzGreyType	greyType;
  double	**bufData = NULL;
  double	*tDP0,
  		*tDP1,
		*bufMedian = NULL;
  int		*tIP0,
  		*bufMarks = NULL;
  WlzObject	*dstObj = NULL;
  WlzIVertex2	bufPos,
		dstPos;
  WlzIBox2	dstBox;
  void		*dstGreyValues = NULL;
  WlzInterval	*dstItvBase = NULL,
  		*dstItv0,
  		*dstItv1;
  WlzDomain	dstDom,
  		srcDom;
  WlzValues	dstValues;
  WlzGreyP	tGP0;
  WlzPixelV	backgroundPix;
  WlzIntervalWSpace bufIWsp,
		srcIWsp;
  WlzGreyWSpace	bufGWsp,
		srcGWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstValues.core = NULL;
  srcDom = srcObj->domain;

  dstBox.xMin = (srcDom.i->kol1 + (kernelSz.vtX / 2) + samFac.vtX - 1) / 
		samFac.vtX;
  dstBox.yMin = (srcDom.i->line1 + (kernelSz.vtY / 2) + samFac.vtY - 1) /
		samFac.vtY;
  dstBox.xMax = (srcDom.i->lastkl - (kernelSz.vtX / 2)) / samFac.vtX;
  dstBox.yMax = (srcDom.i->lastln - (kernelSz.vtY / 2)) / samFac.vtY;
  dstWidth = dstBox.xMax - dstBox.xMin + 1;
  srcWidth = srcDom.i->lastkl - srcDom.i->kol1 + 1;
  backgroundPix = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues = WlzSampleObjConstructRectValues(&dstGreyValues, greyType,
						dstBox,
						backgroundPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((maxItvCount = WlzSampleObjEstMaxIntervals(srcDom, 
						  dstBox.yMin * samFac.vtY,
						  dstBox.yMax * samFac.vtY,
						  samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc((unsigned long )maxItvCount *
      					       sizeof(WlzInterval))) == NULL) ||
         ((dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
	 				    dstBox.yMin, dstBox.yMax,
					    dstBox.xMin,
					    dstBox.xMax,
					    &errNum)) == NULL))
      {
        if(dstItvBase)
	{
	  AlcFree(dstItvBase);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
        }
      }
      else
      {
	dstDom.i->freeptr = WlzPushFreePtr(dstDom.i->freeptr,
					   (void *)dstItvBase,
					   &errNum);
      }
    }
    else
    {
      dstObj = WlzMakeEmpty(&errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core)
  {
    if(AlcDouble2Malloc(&bufData, kernelSz.vtY, srcWidth) == ALC_ER_NONE)
    {
      if((bufMarks = (int *)AlcMalloc((unsigned long )(kernelSz.vtY) *
      				      sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	tIP0 = bufMarks;
	tI0 = kernelSz.vtY;
	tI1 = srcDom.i->line1 - (2 * kernelSz.vtY); 	 /* Use invalid line */
	while(tI0-- > 0)        /* Mark buffer lines stale with invalid line */
	{
	  *tIP0++ = tI1;
	}
	switch(backgroundPix.type)
	{
	  case WLZ_GREY_DOUBLE:
	    backgroundVal = backgroundPix.v.dbv;
	    break;
	  case WLZ_GREY_FLOAT:
	    backgroundVal = (double )(backgroundPix.v.flv);
	    break;
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (samFn == WLZ_SAMPLEFN_MEDIAN))
  {
    bufMedianSz = kernelSz.vtX * kernelSz.vtY;
    if((bufMedian = (double *)AlcMalloc((unsigned long )bufMedianSz *
    				        sizeof(double))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstValues.core && dstDom.core &&
     bufData && bufMarks)
  {
    totItvCount = itvCount = 0;
    dstItv0 = dstItvBase;
    dstItv1 = dstItvBase;
    bufBase = WLZ_ABS(srcDom.i->line1) + kernelSz.vtY;  /* Make sure not -ve */
    if(((errNum = WlzInitGreyScan(srcObj, &srcIWsp,
    				  &srcGWsp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(srcObj, &bufIWsp,
    				  &bufGWsp)) == WLZ_ERR_NONE))
    {
      bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
    }
    while((errNum == WLZ_ERR_NONE) && (WlzNextGreyInterval(&srcIWsp) == 0))
    {
      if(((srcIWsp.linpos % samFac.vtY) == 0) &&
      	 ((tI0 = (srcIWsp.linpos / samFac.vtY)) >= dstBox.yMin) &&
	 (tI0 <= dstBox.yMax))
      {
	dstInvLeftPos = WLZ_MAX((srcIWsp.lftpos + samFac.vtX - 1) / samFac.vtX,
				dstBox.xMin);
	dstInvRgtPos = WLZ_MIN(srcIWsp.rgtpos / samFac.vtX, dstBox.xMax);
	dstInvWidth = dstInvRgtPos - dstInvLeftPos + 1;
	dstPos.vtY = srcIWsp.linpos / samFac.vtY;
	dstOffset = ((dstPos.vtY - dstBox.yMin) * dstWidth) +
		    dstInvLeftPos - dstBox.xMin;
	dstItv1->ileft = dstInvLeftPos - dstBox.xMin;
	dstItv1->iright = dstInvRgtPos - dstBox.xMin;
	++dstItv1;
	totItvCount += ++itvCount;
	if(totItvCount >= maxItvCount)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	{
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  for(idY = 0; idY < kernelSz.vtY; ++idY) /* Check and update buffer */
	  {
	    bufPos.vtX = srcDom.i->kol1;
	    idB = (bufBase + bufPos.vtY) % kernelSz.vtY;
	    if(*(bufMarks + idB) != bufPos.vtY)  /* Check if update required */
	    {
	      *(bufMarks + idB) = bufPos.vtY;
	      if((bufPos.vtY < srcDom.i->line1) ||
		 (bufPos.vtY > srcDom.i->lastln))
	      {
		WlzValueSetDouble(*(bufData + idB), backgroundVal, srcWidth);
	      }
	      else
	      {
		while((bufIWsp.linpos <= bufPos.vtY) && (bufIwspFlag == 0))
		{
		  if(bufIWsp.linpos == bufPos.vtY)
		  {
		    if((tI1 = bufIWsp.lftpos - bufPos.vtX) > 0)
		    {
		      tDP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		      WlzValueSetDouble(tDP0, backgroundVal, tI1);
		    }
		    tI1 = bufIWsp.rgtpos - bufIWsp.lftpos + 1;
		    tDP0 = *(bufData + idB) + bufIWsp.lftpos - srcDom.i->kol1;
		    switch(greyType)
		    {
		      case WLZ_GREY_FLOAT:
			WlzValueCopyFloatToDouble(tDP0, bufGWsp.u_grintptr.flp,
						  tI1);
			break;
		      case WLZ_GREY_DOUBLE:
			WlzValueCopyDoubleToDouble(tDP0, bufGWsp.u_grintptr.dbp,
						   tI1);
			break;
		    }
		    bufPos.vtX = bufIWsp.rgtpos + 1;
		  }
		  bufIwspFlag = WlzNextGreyInterval(&bufIWsp);
		}
		if((tI1 = srcDom.i->lastkl - bufPos.vtX) > 0)
		{
		  tDP0 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  WlzValueSetDouble(tDP0, backgroundVal, tI1);
		}
	      }
	    }
	    ++(bufPos.vtY);
	  }
	  /* Sample by rank through the interval */
	  bufPos.vtX = dstInvLeftPos * samFac.vtX;
	  bufPos.vtY = srcIWsp.linpos - (kernelSz.vtY / 2);
	  switch(greyType)
	  {
	    case WLZ_GREY_FLOAT:
	      tGP0.flp = (float *)dstGreyValues + dstOffset;
	      break;
	    case WLZ_GREY_DOUBLE:
	      tGP0.dbp = (double *)dstGreyValues + dstOffset;
	      break;
	  }
	  for(tI0 = 0; tI0 < dstInvWidth; ++tI0)
	  {
	    tI2 = 0;
	    switch(samFn)
	    {
	      case WLZ_SAMPLEFN_MIN:
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tDP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  if(tI2)
		  {
		    if(*tDP1 < tD0)
		    {
		      tD0 = *tDP1;
		    }
		  }
		  else
		  {
		    tD0 = *tDP1;
		    tI2 = 1;
		  }
		  while(idX-- > 0)
		  {
		    if(*++tDP1 < tD0)
		    {
		      tD0 = *tDP1;
		    }
		    if(*++tDP1 < tD0)
		    {
		      tD0 = *tDP1;
		    }
		  }
		}
		break;
	      case WLZ_SAMPLEFN_MAX:
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tDP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  if(tI2)
		  {
		    if(*tDP1 > tD0)
		    {
		      tD0 = *tDP1;
		    }
		  }
		  else
		  {
		    tD0 = *tDP1;
		    tI2 = 1;
		  }
		  while(idX-- > 0)
		  {
		    if(*++tDP1 > tD0)
		    {
		      tD0 = *tDP1;
		    }
		    if(*++tDP1 > tD0)
		    {
		      tD0 = *tDP1;
		    }
		  }
		}
		break;
	      case WLZ_SAMPLEFN_MEDIAN:
		tDP0 = bufMedian;
		for(idY = 0; idY < kernelSz.vtY; ++idY)
		{
		  idB = (bufBase + (bufPos.vtY)++) % kernelSz.vtY;
		  tDP1 = *(bufData + idB) + bufPos.vtX - srcDom.i->kol1;
		  idX = kernelSz.vtX / 2;
		  /* First copy to the median buffer. */
		  *tDP0++ = *tDP1++;
		  while(idX-- > 0)
		  {
		    *tDP0++ = *tDP1++;
		    *tDP0++ = *tDP1++;
		  }
		  /* Then find median. */
		  tD0 = WlzValueMedianDouble(bufMedian, bufMedianSz);
		}
		break;
	    }
	    switch(greyType)
	    {
	      case WLZ_GREY_FLOAT:
		*(tGP0.flp)++ = (float )tD0; /* Fits because rank operation. */
		break;
	      case WLZ_GREY_DOUBLE:
		*(tGP0.dbp)++ = tD0;
		break;
	    }
	    bufPos.vtX += samFac.vtX;
	  }

	  if((errNum == WLZ_ERR_NONE) &&
	     (itvCount > 0) && (srcIWsp.intrmn == 0))
	  {
	    WlzMakeInterval(dstPos.vtY, dstDom.i, itvCount, dstItv0);
	    dstItv0 = dstItv1;
	    itvCount = 0;
	  }
	}
      }
    }
    (void )WlzStandardIntervalDomain(dstDom.i);
    dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstValues, NULL, NULL,
    			 &errNum);
    if(dstObj)
    {
      (void )WlzSetBackground(dstObj, WlzGetBackground(srcObj, NULL));
    }
  }
  if(bufData) 					      /* Free up buffer data */
  {
    AlcDouble2Free(bufData);
  }
  if(bufMarks)
  {
    AlcFree(bufMarks);
  }
  if(bufMedian)
  {
    AlcFree(bufMedian);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(dstObj == NULL)				        /* Clear up on error */
  {
    (void )WlzFreeValues(dstValues);
    if(dstDom.core)
    {
      (void )WlzFreeDomain(dstDom);
    }
    else if(dstItvBase)
    {
      AlcFree(dstItvBase);
    }
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzSampleObjGaussKernelD				*
* Returns:	int:			Non zero if kernel computed ok.	*
* Purpose:	Computes a (double) floating point gaussian convolution	*
*		kernel for WlzSampleObj().				*
* Global refs:	-							*
* Parameters:	int **kernel:		Space for kernel allocated with	*
*					AlcDouble2Malloc().		*
*		WlzIVertex2 kernelSz:	Kernel size (3, 5, 7, ...).	*
*		WlzIVertex2 samFac:	Sampling factor (2, 3, 4, ...).	*
************************************************************************/
static int	WlzSampleObjGaussKernelD(double **kernel, WlzIVertex2 kernelSz,
				         WlzIVertex2 samFac)
{
  double	tD0,
		kX,
		kY,
		kR,
		min,
  		sum = 0;
  int		ok = 0,
  		idX,
		idY;
  const double	widFac = 1.0 / (16.0 * WLZ_M_LN2);

  if(kernel)
  {
    kR = -widFac * ((samFac.vtX * samFac.vtX) + (samFac.vtY * samFac.vtY));
    for(idY = 0; idY < kernelSz.vtY; ++idY)
    {
      kY = idY - (kernelSz.vtY / 2);
      kY *= kY;
      for(idX = 0; idX < kernelSz.vtX; ++idX)
      {
	kX = idX - (kernelSz.vtX / 2);
	kX *= kX;
        tD0 = exp((kX + kY) / kR);
	if(((idX == 0) && (idY == 0)) || (tD0 < min))
	{
	  min = tD0;
	}
        *(*(kernel + idY) + idX) = tD0;
	sum += tD0;
      }
    }
    sum -= min * kernelSz.vtY * kernelSz.vtX;
    for(idY = 0; idY < kernelSz.vtY; ++idY)
    {
      for(idX = 0; idX < kernelSz.vtX; ++idX)
      {
	*(*(kernel + idY) + idX) = (*(*(kernel + idY) + idX) - min) / sum;
      }
    }
    ok = 1;
  }
  return(ok);
}

/************************************************************************
* Function:	WlzSampleObjGaussKernelI				*
* Returns:	int:			Non zero if kernel computed ok.	*
* Purpose:	Computes an integral gaussian convolution kernel for	*
*		WlzSampleObj().						*
* Global refs:	-							*
* Parameters:	int **kernel:		Space for kernel allocated with	*
*					AlcInt2Malloc().		*
*		WlzIVertex2 kernelSz:	Kernel size (3, 5, 7, ...).	*
*		int *kernelSum:		Destination pointer for kernel 	*
*					sum.				*
*		WlzIVertex2 samFac:	Sampling factor (2, 3, 4, ...).	*
************************************************************************/
static int	WlzSampleObjGaussKernelI(int **kernel, WlzIVertex2 kernelSz,
				         int *kernelSum, WlzIVertex2 samFac)
{
  int		sum = 0,
  		ok = 0,
  		idX,
		idY;
  double	tD0;
  double	**kernelD = NULL;

  if((AlcDouble2Malloc(&kernelD,
  		       kernelSz.vtY, kernelSz.vtX) == ALC_ER_NONE) &&
     WlzSampleObjGaussKernelD(kernelD, kernelSz, samFac))
  {
    for(idY = 0; idY < kernelSz.vtY; ++idY)
    {
      for(idX = 0; idX < kernelSz.vtX; ++idX)
      {
	tD0 = WLZ_SAMPLE_KERNEL_INORM * *(*(kernelD + idY) + idX);
        sum += *(*(kernel + idY) + idX) = WLZ_NINT(tD0);
      }
    }
    ok = 1;
    *kernelSum = sum;
  }
  if(kernelD)
  {
    AlcDouble2Free(kernelD);
  }
  return(ok);
}

/************************************************************************
* Function:	WlzSampleObjMeanKernelD					*
* Returns:	int:			Non zero if kernel computed ok.	*
* Purpose:	Computes a (double) floating point mean convolution	*
*		kernel for WlzSampleObj().				*
* Global refs:	-							*
* Parameters:	double **kernel:	Space for kernel allocated with	*
*					AlcDouble2Malloc().		*
*		WlzIVertex2 kernelSz:	Kernel size (3, 5, 7, ...).	*
************************************************************************/
static int	WlzSampleObjMeanKernelD(double **kernel, WlzIVertex2 kernelSz)
{
  double 	tD0;
  int		ok = 0,
  		idX,
		idY;

  tD0 = 1.0 / (kernelSz.vtX * kernelSz.vtY);
  for(idY = 0; idY < kernelSz.vtY; ++idY)
  {
    for(idX = 0; idX < kernelSz.vtX; ++idX)
    {
      *(*(kernel + idY) + idX) = tD0;
    }
  }
  ok = 1;
  return(ok);
}

/************************************************************************
* Function:	WlzSampleObjMeanKernelI					*
* Returns:	int:			Non zero if kernel computed ok.	*
* Purpose:	Computes an integral mean convolution kernel for	*
*		WlzSampleObj().						*
* Global refs:	-							*
* Parameters:	int **kernel:		Space for kernel allocated with	*
*					AlcInt2Malloc().		*
*		WlzIVertex2 kernelSz:	Kernel size (3, 5, 7, ...).	*
*		int *kernelSum:		Destination pointer for kernel 	*
*					sum.				*
************************************************************************/
static int	WlzSampleObjMeanKernelI(int **kernel, WlzIVertex2 kernelSz,
					int *kernelSum)
{
  int		ok = 0,
  		tI0,
  		idX,
		idY;

  tI0 = WLZ_SAMPLE_KERNEL_INORM / (kernelSz.vtX * kernelSz.vtY);
  for(idY = 0; idY < kernelSz.vtY; ++idY)
  {
    for(idX = 0; idX < kernelSz.vtX; ++idX)
    {
      *(*(kernel + idY) + idX) = tI0;
    }
  }
  *kernelSum = tI0 * kernelSz.vtX * kernelSz.vtY;
  ok = 1;
  return(ok);
}

/************************************************************************
* Function:     WlzSampleObjConstructRectValues				*
* Returns:      WlzValues :         	New value table (rectangular). 	*
* Purpose:      Constructs a new (rectangular) value table with grey	*
*		values allocated but not set.			        *
* Global refs:  -							*
* Parameters:   void **dstValues:       Destination pointer for grey	*
*                                       values.				*
*               WlzGreyType greyType:   Grey type.			*
*               WlzIBox2 rBox:           Box encoding the rectangle	*
*                                       verticies.			*
*		WlzPixelV bgdPix:	Background value.		*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number, may be NULL if not	*
*					required.			*
************************************************************************/
static WlzValues WlzSampleObjConstructRectValues(void **dstValues,
						  WlzGreyType greyType,
						  WlzIBox2 rBox,
						  WlzPixelV bgdPix,
						  WlzErrorNum *dstErrNum)
{
  int           bCount,
                rSzX;
  WlzObjectType	gTabType;
  void		*greyValues = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzValues  	values;
 
  values.core = NULL;
  if((rBox.xMin > rBox.xMax) || (rBox.yMin > rBox.yMax))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    rSzX = rBox.xMax - rBox.xMin + 1;
    bCount = rSzX * (rBox.yMax - rBox.yMin + 1);
    switch(greyType)
    {
      case WLZ_GREY_INT:
        bCount *= sizeof(int);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_INT,
				    &errNum);
        break;
      case WLZ_GREY_SHORT:
        bCount *= sizeof(short);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_SHORT,
				    &errNum);
        break;
      case WLZ_GREY_UBYTE:
        bCount *= sizeof(char);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_UBYTE,
				    &errNum);
        break;
      case WLZ_GREY_FLOAT:
	bCount *= sizeof(float);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_FLOAT,
				    &errNum);
        break;
      case WLZ_GREY_DOUBLE:
	bCount *= sizeof(double);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_DOUBLE,
			            &errNum);
        break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((greyValues = AlcMalloc(bCount)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      values.r = WlzMakeRectValueTb(gTabType, rBox.yMin, rBox.yMax,
				    rBox.xMin, rSzX, bgdPix,
				    (int *)greyValues, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    values.r->freeptr = WlzPushFreePtr(values.r->freeptr,
    				       greyValues, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstValues)
    {
      *dstValues = greyValues;
    }
  }
  else
  {
    if(values.core)
    {
      (void )WlzFreeValues(values);
    }
    else if(greyValues)
    {
      AlcFree(greyValues);
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  return(values);
}

/************************************************************************
* Function:	WlzSampleObjEstMaxIntervals				*
* Returns:	int:			Number of intervals.		*
* Purpose:	Estimates the maximum number of intervals required	*
*		for the interval domain of a subsampled object.  	*
*		The number returned may be much bigger than required.	*
* Global refs:	-							*
* Parameters:	WlzDomain srcDom:	Source objects domain.		*
*		int line1:		First line.			*
*		int lastLn:		Last line.			*
*		WlzIVertex2 samFac:	Sampling factor.		*
************************************************************************/
static int	WlzSampleObjEstMaxIntervals(WlzDomain srcDom,
					    int line1, int lastLn,
					    WlzIVertex2 samFac)
{
    int		srcLine,
    		srcLastLine,
		itvCount = 0;

  if(srcDom.core &&
     ((srcDom.core->type == WLZ_INTERVALDOMAIN_INTVL) ||
      (srcDom.core->type == WLZ_INTERVALDOMAIN_RECT)) &&
     (samFac.vtX > 0) && (samFac.vtY > 0))
  {
    line1 -= samFac.vtY;
    lastLn += samFac.vtY;
    srcLine = WLZ_MAX(srcDom.i->line1, line1);
    srcLastLine = WLZ_MIN(srcDom.i->lastln, lastLn);
    switch(srcDom.i->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
	if(srcLine <= srcLastLine)
	{
	  while(srcLine <= srcLastLine)
	  {
	    itvCount += (srcDom.i->intvlines + srcLine -
	    		 srcDom.i->line1)->nintvs;
	    ++srcLine;
	  }
	  itvCount += ((srcLastLine - srcLine) /  samFac.vtY) + samFac.vtY + 1;
	}
        break;
      case WLZ_INTERVALDOMAIN_RECT:
	itvCount = ((srcLastLine - srcLine) /  samFac.vtY) + samFac.vtY + 1;
        break;
    }
  }
  return(itvCount);
}

#undef WLZ_SAMPLE_KERNEL_INORM
