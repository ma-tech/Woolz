#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzSampleObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzSampleObj.c
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
* \brief	Subsamples objects using an integer sampling factor
* 		and a convolution kernel.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

#define WLZ_SAMPLE_KERNEL_INORM	(0x000100)

static WlzObject 		*WlzSampleObj2D(
				  WlzObject *,
				  WlzIVertex2,
				  WlzSampleFn,
				  WlzIVertex2,
				  WlzErrorNum *);
static WlzObject 		*WlzSampleObj3D(
				  WlzObject *srcObj,
				  WlzIVertex3 samFac,
			          WlzSampleFn samFn,
				  WlzIVertex3 kernelSz,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzSampleObjIDom(
				  WlzObject *srcObj,
				  WlzIVertex2 samFac,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzSampleObjConvI(
				  WlzObject *,
				  int **,
				  WlzIVertex2,
				  int,
				  WlzIVertex2,
				  WlzErrorNum *);
static WlzObject		*WlzSampleObjConvD(
				  WlzObject *,
				  double **,
				  WlzIVertex2,
				  WlzIVertex2,
				  WlzErrorNum *);
static WlzObject	        *WlzSampleObjRankI(
				  WlzObject *,
				  WlzIVertex2,
				  WlzSampleFn,
				  WlzIVertex2,
				  WlzErrorNum *);
static WlzObject		*WlzSampleObjRankD(
				  WlzObject *,
				  WlzIVertex2,
				  WlzSampleFn,
				  WlzIVertex2,
				  WlzErrorNum *);
static WlzValues 		WlzSampleObjConstructRectValues(
				  void **,
				  WlzGreyType,
				  WlzIBox2,
				  WlzPixelV,
				  WlzErrorNum *);
static int			WlzSampleObjEstIntervals(
				  WlzDomain,
				  int,
				  int,
				  WlzIVertex2);
static int			WlzSampleObjGaussKernelD(
				  double **,
				  WlzIVertex2,
				  WlzIVertex2);
static int			WlzSampleObjGaussKernelI(
				  int **,
				  WlzIVertex2,
				  int *,
				  WlzIVertex2);
static int			WlzSampleObjMeanKernelD(
				  double **,
				  WlzIVertex2);
static int			WlzSampleObjMeanKernelI(
				  int **,
				  WlzIVertex2,
				  int *);
static WlzErrorNum 		WlzSampleObjMoreIntervals(
				  WlzDomain dstDom,
				  int *delItvCount,
				  int itvCount,
				  WlzInterval **dstItv0,
				  WlzInterval **dstItv1);

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using the given sampling
*               factor and sampling method.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	samFn			Sampling method.
* \param	dstErr			Destination pointer for error,
					may be NULL.
*/
WlzObject	*WlzSampleObj(WlzObject *srcObj, WlzIVertex3 samFac,
			      WlzSampleFn samFn,
			      WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzIVertex2 	samFac2,
  		kernelS2;
  WlzIVertex3 	kernelS3;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzSampleObj FE 0x%lx {%d %d} %d 0x%lx\n",
	   (unsigned long )srcObj, samFac.vtX, samFac.vtY,
	   (int )samFn, (unsigned long )dstErr));
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
	samFac2.vtX = samFac.vtX;
	samFac2.vtY = samFac.vtY;
	if(samFac2.vtX == 1)
	{
	  kernelS2.vtX = 1;
	}
	else
	{
	  kernelS2.vtX = (samFac2.vtX % 2)? samFac2.vtX + 2: samFac2.vtX + 1;
	}
	if(samFac2.vtY == 1)
	{
	  kernelS2.vtY = 1;
	}
	else
	{
	  kernelS2.vtY = (samFac2.vtY % 2)? samFac2.vtY + 2: samFac2.vtY + 1;
	}
	dstObj = WlzSampleObj2D(srcObj, samFac2, samFn, kernelS2,
				   &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	if(samFac.vtX == 1)
	{
	  kernelS3.vtX = 1;
	}
	else
	{
	  kernelS3.vtX = (samFac.vtX % 2)? samFac.vtX + 2: samFac.vtX + 1;
	}
	if(samFac.vtY == 1)
	{
	  kernelS3.vtY = 1;
	}
	else
	{
	  kernelS3.vtY = (samFac.vtY % 2)? samFac.vtY + 2: samFac.vtY + 1;
	}
	if(samFac.vtZ == 1)
	{
	  kernelS3.vtZ = 1;
	}
	else
	{
	  kernelS3.vtZ = (samFac.vtZ % 2)? samFac.vtZ + 2: samFac.vtZ + 1;
	}
	dstObj = WlzSampleObj3D(srcObj, samFac, samFn, kernelS3,
				   &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzSampleObj FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given 2D domain object using the given
*               sampling factor kernel size and sampling method.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	samFn			Sampling method.
* \param	kernelSz		Size of the convolution kernel.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzSampleObj2D(WlzObject *srcObj, WlzIVertex2 samFac,
			         WlzSampleFn samFn, WlzIVertex2 kernelSz,
				 WlzErrorNum *dstErr)
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
	dstObj = WlzSampleObjPoint2D(srcObj, samFac, &errNum);
      }
      else
      {
	switch(greyType)
	{
	  case WLZ_GREY_INT:
	  case WLZ_GREY_SHORT:
	  case WLZ_GREY_UBYTE:
	  case WLZ_GREY_RGBA:
	    integralGrey = 1;
	    break;
	  case WLZ_GREY_FLOAT:
	  case WLZ_GREY_DOUBLE:
	    integralGrey = 0;
	    break;
	  default:
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given 3D domain object using the given
*               sampling factor kernel size and sampling method.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	samFn			Sampling method.
* \param	kernelSz		Size of the convolution kernel.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzSampleObj3D(WlzObject *srcObj, WlzIVertex3 samFac,
			         WlzSampleFn samFn, WlzIVertex3 kernelSz,
				 WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj->values.core == NULL)
  {
    dstObj = WlzSampleObjPoint3D(srcObj, samFac, &errNum);
  }
  else
  {
    switch(samFn)
    {
      case WLZ_SAMPLEFN_POINT:
	dstObj = WlzSampleObjPoint3D(srcObj, samFac, &errNum);
        break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
        break;
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New sampled object.
* \ingroup      WlzTransform
* \brief	Samples the given object's interval domain only using
*               the given sampling factor.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzObject *WlzSampleObjIDom(WlzObject *srcObj, WlzIVertex2 samFac,
				   WlzErrorNum *dstErr)
{
  int		itvCount,
		totItvCount,
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;

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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItv0 = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
        dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr, (void *)dstItv0,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	    errNum = WlzSampleObjMoreIntervals(dstDom,
	    				       &delItvCount, itvCount,
					       &dstItv0, &dstItv1);
	    totItvCount = itvCount;
	    maxItvCount = delItvCount;
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New sampled object.
* \ingroup      WlzTransform
* \brief	Samples the given object's plane domain only using
*               the given sampling factor.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every voxel == 1,
*                                       every other voxel == 2, ....
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
WlzObject 	*WlzSampleObjPoint3D(WlzObject *srcObj, WlzIVertex3 samFac,
				     WlzErrorNum *dstErr)
{
  int		dPlCnt,
		dPlIdx,
  		sPlIdx;
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
  		*dstObj = NULL;
  WlzIVertex2	samFac2;
  WlzIBox3	srcBox,
  		dstBox;
  WlzValues	dumVal,
  		dstVal,
  		srcVal;
  WlzDomain	dstDom,
  		srcDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dumVal.core = NULL;
  dstVal.core = NULL;
  srcDom = srcObj->domain;
  srcVal = srcObj->values;
  if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    samFac2.vtX = samFac.vtX;
    samFac2.vtY = samFac.vtY;
    srcBox.xMin = srcDom.p->kol1;
    srcBox.xMax = srcDom.p->lastkl;
    srcBox.yMin = srcDom.p->line1;
    srcBox.yMax = srcDom.p->lastln;
    srcBox.zMin = srcDom.p->plane1;
    srcBox.zMax = srcDom.p->lastpl;
    dstBox.xMin = (srcBox.xMin < 0) ?
		  (srcBox.xMin  - samFac.vtX + 1) / samFac.vtX :
		  (srcBox.xMin  + samFac.vtX - 1) / samFac.vtX;
    dstBox.xMax = srcBox.xMax / samFac.vtX;
    dstBox.yMin = (srcBox.yMin < 0) ?
		  (srcBox.yMin - samFac.vtY + 1) / samFac.vtY :
		  (srcBox.yMin + samFac.vtY - 1) / samFac.vtY;
    dstBox.yMax = srcBox.yMax / samFac.vtY;
    dstBox.zMin = (srcBox.zMin < 0) ?
		  (srcBox.zMin - samFac.vtZ + 1) / samFac.vtZ :
		  (srcBox.zMin + samFac.vtZ - 1) / samFac.vtZ;
    dstBox.zMax = srcBox.zMax / samFac.vtZ;
    dstDom.p = WlzMakePlaneDomain(srcDom.p->type, dstBox.zMin, dstBox.zMax,
    				  dstBox.yMin, dstBox.yMax,
				  dstBox.xMin, dstBox.xMax,
				  &errNum);
    if(srcVal.core && (srcVal.core->type != WLZ_EMPTY_OBJ))
    {
      dstVal.vox = WlzMakeVoxelValueTb(srcObj->values.vox->type,
      				       dstBox.zMin, dstBox.zMax,
				       WlzGetBackground(srcObj, NULL),
				       NULL, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sPlIdx = 0;
    dPlIdx = sPlIdx * samFac.vtZ;
    dPlCnt = dstBox.zMax - dstBox.zMin + 1;
    while((errNum == WLZ_ERR_NONE) && (dPlCnt-- > 0))
    {
      if(dstVal.vox)
      {
	tObj0 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			    *(srcDom.p->domains + sPlIdx),
			    *(srcVal.vox->values + sPlIdx),
			    NULL, NULL, &errNum);
      }
      else
      {
        tObj0 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			    *(srcDom.p->domains + sPlIdx),
			    dumVal,
			    NULL, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(dstVal.vox)
	{
          tObj1 = WlzSampleObjPoint2D(tObj0, samFac2, &errNum);
	}
	else
	{
          tObj1 = WlzSampleObjIDom(tObj0, samFac2, &errNum);
	}
      }
      if(tObj0)
      {
	(void )WlzFreeObj(tObj0);
	tObj0 = NULL;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        *(dstDom.p->domains + dPlIdx) = tObj1->domain;
	if(dstVal.vox)
	{
	  *(dstVal.vox->values + dPlIdx) = tObj1->values;
	}
      }
      if(tObj1)
      {
	tObj1->domain.core = NULL;
	tObj1->values.core = NULL;
	(void )WlzFreeObj(tObj1);
	tObj1 = NULL;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        ++dPlIdx;
	sPlIdx += samFac.vtZ;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dstDom.p, dstVal.vox);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstDom.p->voxel_size[0] = srcDom.p->voxel_size[0];
    dstDom.p->voxel_size[1] = srcDom.p->voxel_size[1];
    dstDom.p->voxel_size[2] = srcDom.p->voxel_size[2];
    dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstVal,
    		         NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstDom.core = NULL;
    dstVal.core = NULL;
  }
  if(dstDom.core)
  {
    (void )WlzFreeDomain(dstDom);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}


/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using a simple point sampling
*               method and the given sampling factor.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
WlzObject 	*WlzSampleObjPoint2D(WlzObject *srcObj, WlzIVertex2 samFac,
				     WlzErrorNum *dstErr)
{
  int		tI0,
		itvCount,
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;
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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItv0 = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
        dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr, (void *)dstItv0,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	    errNum = WlzSampleObjMoreIntervals(dstDom,
	    				       &delItvCount, itvCount,
					       &dstItv0, &dstItv1);
	    totItvCount = itvCount;
	    maxItvCount = delItvCount;
	  }
	  if(errNum == WLZ_ERR_NONE)
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
		dstPix.ubp = (WlzUByte *)dstGreyValues + dstOffset;
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
	      case WLZ_GREY_RGBA:
		srcPix.rgbp = srcGWsp.u_grintptr.rgbp + srcOffset;
		dstPix.rgbp = (WlzUInt *)dstGreyValues + dstOffset;
		tI0 = dstInvWidth;
		while(tI0-- > 0)
		{
		  *(dstPix.rgbp)++ = *(srcPix.rgbp);
		  srcPix.rgbp += samFac.vtX;
		}
		break;
	      default:
	        errNum = WLZ_ERR_GREY_TYPE;
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using the given integer
*               convolution kernel and the given sampling factor.
*               The domain of the the new object is smaller than the
*               source because of the convolution.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	kernel			Integer convolution kernel.
* \param	kernelSz		Size of the convolution kernel.
* \param	kernelSum		Kernel sum for normalization.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzObject *WlzSampleObjConvI(WlzObject *srcObj, int **kernel,
				    WlzIVertex2 kernelSz, int kernelSum,
				    WlzIVertex2 samFac, WlzErrorNum *dstErr)
{
  int		tI0,
  		tI1,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		itvCount,
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;
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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
	dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
					     (void *)dstItvBase,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	  default:
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
	  errNum = WlzSampleObjMoreIntervals(dstDom,
	  				     &delItvCount, itvCount,
					     &dstItv0, &dstItv1);
	  totItvCount = itvCount;
	  maxItvCount = delItvCount;
	}
	if(errNum == WLZ_ERR_NONE)
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
		      default:
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
	      tGP0.ubp = (WlzUByte *)dstGreyValues + dstOffset;
	      break;
	    default:
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
		*(tGP0.ubp)++ = (WlzUByte )tI1;
		break;
	      default:
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
  if(dstErr)
  {
    *dstErr = errNum;
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

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using the given double
*               convolution kernel and the given sampling factor.
*               The domain of the the new object is smaller than the
*               source because of the convolution.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	kernel			Double convolution kernel.
* \param	kernelSz		Size of the convolution kernel.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzObject *WlzSampleObjConvD(WlzObject *srcObj, double **kernel,
				    WlzIVertex2 kernelSz, WlzIVertex2 samFac,
				    WlzErrorNum *dstErr)
{
  int		tI0,
  		tI1,
		idB,
		idX,
		idY,
		bufBase,
		bufIwspFlag,
		itvCount,
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;

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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
	dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
					     (void *)dstItvBase,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	  default:
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
	  errNum = WlzSampleObjMoreIntervals(dstDom,
	  				     &delItvCount, itvCount,
					     &dstItv0, &dstItv1);
	  totItvCount = itvCount;
	  maxItvCount = delItvCount;
	}
	if(errNum == WLZ_ERR_NONE)
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
		      default:
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
	    default:
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
	      default:
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
  if(dstErr)
  {
    *dstErr = errNum;
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

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using the given (rank)
*               sampling method and sampling factor.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	samFn			Rank sampling method, eg min,
*                                       max, median,....
* \param	kernelSz		Size of the convolution kernel.
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzObject *WlzSampleObjRankI(WlzObject *srcObj, WlzIVertex2 samFac,
				    WlzSampleFn samFn, WlzIVertex2 kernelSz,
				    WlzErrorNum *dstErr)
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
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;

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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
	dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
					     (void *)dstItvBase,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	  default:
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
	  errNum = WlzSampleObjMoreIntervals(dstDom,
	  				     &delItvCount, itvCount,
					     &dstItv0, &dstItv1);
	  totItvCount = itvCount;
	  maxItvCount = delItvCount;
	}
	if(errNum == WLZ_ERR_NONE)
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
		      default:
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
	      tGP0.ubp = (WlzUByte *)dstGreyValues + dstOffset;
	      break;
	    default:
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
	      default:
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
		*(tGP0.ubp)++ = (WlzUByte )tI1; /* Fits because rank op. */
		break;
	      default:
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
  if(dstErr)
  {
    *dstErr = errNum;
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

/*!
* \return	New sampled object.
* \ingroup	WlzTransform
* \brief	Samples the given object using the given (rank)
*               sampling method and sampling factor.
*               This function assumes it's parameters to be valid.
* \param	srcObj			Given source object.
* \param	samFac			Sampling factor for both rows
*                                       and columns. Every pixel == 1,
*                                       every other pixel == 2, ....
* \param	samFn			Rank sampling method, eg min,
*                                       max, median,....
* \param	kernelSz		Size of the convolution kernel.
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzObject *WlzSampleObjRankD(WlzObject *srcObj, WlzIVertex2 samFac,
				    WlzSampleFn samFn, WlzIVertex2 kernelSz,
				    WlzErrorNum *dstErr)
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
		delItvCount,
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
  AlcErrno	alcErr = ALC_ER_NONE;

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
    if((delItvCount = WlzSampleObjEstIntervals(srcDom, 
					       dstBox.yMin * samFac.vtY,
					       dstBox.yMax * samFac.vtY,
					       samFac)) > 0)
    {
      if(((dstItvBase = (WlzInterval *)AlcMalloc(delItvCount *
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
	maxItvCount = delItvCount;
	dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
					     (void *)dstItvBase,
					     &alcErr);
        if(alcErr != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
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
	  default:
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
	  errNum = WlzSampleObjMoreIntervals(dstDom,
	  				     &delItvCount, itvCount,
					     &dstItv0, &dstItv1);
	  totItvCount = itvCount;
	  maxItvCount = delItvCount;
	}
	if(errNum == WLZ_ERR_NONE)
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
		      default:
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
	    default:
	      break;
	  }
	  for(tI0 = 0; tI0 < dstInvWidth; ++tI0)
	  {
	    tI2 = 0;
	    switch(samFn)
	    {
	      case WLZ_SAMPLEFN_MIN:
		tD0 = DBL_MAX;
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
	      default:
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
	      default:
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
  if(dstErr)
  {
    *dstErr = errNum;
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

/*!
* \return	Non-zero if kernel computed without error.
* \ingroup	WlzTransform
* \brief	Computes a (double) floating point gaussian convolution
*               kernel for WlzSampleObj().
* \param	kernel			Space for kernel as allocated with
*                                       AlcDouble2Malloc().
* \param	kernelSz		Kernel size (3, 5, 7, ...).
* \param	samFac			Sampling factor (2, 3, 4, ...).
*/
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
    min = 0.0;			      /* Just to keep lint/compillers happy. */
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

/*!
* \return	Non-zero if the kernel is computed without error.
* \ingroup	WlzTransform
* \brief	Computes an integral gaussian convolution kernel for
*               WlzSampleObj().
* \param	kernel			Space for kernel as allocated with
*                                       AlcInt2Malloc().
* \param	kernelSz		Kernel size (3, 5, 7, ...).
* \param	kernelSum		Destination pointer for kernel sum.
* \param	samFac			Sampling factor (2, 3, 4, ...).
*/
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

/*!
* \return	Non-zero if computed without error.
* \ingroup	WlzTransform
* \brief	Computes a (double) floating point mean convolution
*               kernel for WlzSampleObj().
* \param	kernel			Space for kernel as allocated with
*                                       AlcDouble2Malloc().
* \param	kernelSz		Kernel size (3, 5, 7, ...).
*/
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

/*!
* \return	Non-zero if computed without error.
* \ingroup	WlzTransform
* \brief	Computes an integral mean convolution kernel for
*               WlzSampleObj().
* \param	kernel			Space for kernel as allocated with
*                                       AlcInt2Malloc().
* \param	kernelSz		Kernel size (3, 5, 7, ...).
* \param	kernelSum		Destination pointer for kernel sum.
*/
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

/*!
* \return	New (rectangular) value table.
* \ingroup	WlzTransform
* \brief	Constructs a new (rectangular) value table with grey
*               values allocated but not set.
* \param	dstValues		Destination pointer for grey
*                                       values.
* \param	greyType		Grey type.
* \param	rBox			Box encoding the rectangle
*                                       verticies.
* \param	bgdPix			Background value.
* \param	dstErr			Destination pointer for error,
                                        may be NULL.
*/
static WlzValues WlzSampleObjConstructRectValues(void **dstValues,
						  WlzGreyType greyType,
						  WlzIBox2 rBox,
						  WlzPixelV bgdPix,
						  WlzErrorNum *dstErr)
{
  int           bCount,
                rSzX;
  WlzObjectType	gTabType;
  void		*greyValues = NULL;
  WlzValues  	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  AlcErrno	alcErr = ALC_ER_NONE;
 
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
      case WLZ_GREY_RGBA:
	bCount *= sizeof(WlzUInt);
	gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, WLZ_GREY_RGBA,
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
    values.r->freeptr = AlcFreeStackPush(values.r->freeptr,
					 greyValues, &alcErr);
    if(alcErr != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(values);
}

/*!
* \return	Number of intervals.
* \ingroup	WlzTransform
* \brief	Estimates the maximum number of intervals required
*               for the interval domain of a subsampled object.
*               The number returned may be <, = or > the number of
*               intervals actualy required.
* \param	srcDom			Source objects domain.
* \param	line1			First line.
* \param	lastLn			Last line.
* \param	samFac			Sampling factor.
*/
static int	WlzSampleObjEstIntervals(WlzDomain srcDom,
					    int line1, int lastLn,
					    WlzIVertex2 samFac)
{
    int		line,
		itvCount = 0;

  if(srcDom.core &&
     ((srcDom.core->type == WLZ_INTERVALDOMAIN_INTVL) ||
      (srcDom.core->type == WLZ_INTERVALDOMAIN_RECT)) &&
     (samFac.vtX > 0) && (samFac.vtY > 0) && (line1 <= lastLn))
  {
    switch(srcDom.i->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
        line = line1;
	while(line <= lastLn)
	{
	  itvCount += (srcDom.i->intvlines + line - srcDom.i->line1)->nintvs;
	  line += samFac.vtY;
	}
	break;
      case WLZ_INTERVALDOMAIN_RECT:
	itvCount = ((lastLn - line1) /  samFac.vtY) + samFac.vtY + 1;
        break;
      default:
        break;
    }
  }
  return(itvCount);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Gets more intervals for the sampling functions.
* \param	dstDom			Destination domain.
* \param	delItvCount		Estimate of the number of
*                                       intervals to allocate.
* \param	itvCount		Number of intervals not yet
*                                       committed to the domain.
* \param	dstItv0			First uncommitted interval.
* \param	dstItv1			Next uncommitted interval.
*/
static WlzErrorNum WlzSampleObjMoreIntervals(WlzDomain dstDom,
					     int *delItvCount, int itvCount,
  					     WlzInterval **dstItv0,
  					     WlzInterval **dstItv1)
{
  int		idx;
  WlzInterval	*newItv;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(*delItvCount < itvCount)
  {
    *delItvCount += itvCount;
  }
  if((newItv = (WlzInterval *)AlcMalloc(*delItvCount *
  				        sizeof(WlzInterval))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    dstDom.i->freeptr = AlcFreeStackPush(dstDom.i->freeptr,
					 (void *)newItv, &alcErr);
    if(alcErr != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Copy existing uncommitted intervals. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idx = 0; idx < itvCount; ++idx)
    {
      *(newItv + idx) = *(*dstItv0 + idx);
    }
    *dstItv0 = newItv;
    *dstItv1 = *dstItv0 + itvCount;
  }
  return(errNum);
}

#undef WLZ_SAMPLE_KERNEL_INORM
