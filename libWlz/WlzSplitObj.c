#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSplitObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzSplitObj.c
* \author       Bill Hill
* \date         November 2004
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
* \brief	Functions to split a single object into component
* 		objects.
* \ingroup	WlzBinaryOps
*/

#include <float.h>
#include <Wlz.h>

typedef struct _WlzSplitObjData
{
  int		nLComp;
  int		*compI;
  WlzLong	*lCompSz;
  WlzObject	**lComp;
} WlzSplitObjData;

static int			WlzSplitObjSortSzFn(
				  const void *cData,
				  const void *a,
				  const void *b);

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Splits the reference object into component objects cliped
*		from the reference object, with the bounding box of each
*		of the component objects determined using the pre-processed
*		object. The component objects are returned in size order.
* \param	refObj			Reference object.
* \param	ppObj			Pre-processed object which is
*					normalised to values in the range
*					0 - 255 as WlzUByte greys.
* \param	bWidth			Border width.
* \param	bgdFrac			Minimum fraction of values which are
* 					background values, with range
*					[0.0+ - 1.0-].
* \param	sigma			Histogram smoothing parameter used
*					by WlzHistogramCnvGauss().
* \param	compThrMethod		Method for computing threshold, used
*					in call to WlzCompThresholdVT().
* \param	nReqComp		Number of required components.
* \param	dstNComp		Destination pointer for the number of
*					components extracted, must not be NULL.
* \param	dstComp			Destination pointer for the extracted
*					components, must not be NULL.
*/
WlzErrorNum	WlzSplitObj(WlzObject *refObj, WlzObject *ppObj,
			      int bWidth, double bgdFrac, double sigma,
			      WlzCompThreshType compThrMethod,
			      int nReqComp, int *dstNComp,
			      WlzObject ***dstComp)
{
  int		idC,
		dim = 0;
  WlzObject	*hObj = NULL,
  		*tObj = NULL;
  WlzObject 	**comp = NULL;
  WlzBox	box;
  WlzPixelV	tV;
  WlzSplitObjData split;
  WlzThresholdType tType;
  WlzConnectType lCon = WLZ_0_CONNECTED;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxComp = 1024;

  split.nLComp = 0;
  split.compI = NULL;
  split.lCompSz = NULL;
  split.lComp = NULL;
  if((refObj == NULL) || (ppObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((refObj->domain.core == NULL) || (ppObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((refObj->values.core == NULL) || (ppObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(refObj->type != ppObj->type)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((dstNComp == NULL) || (dstComp == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((bgdFrac < DBL_EPSILON) || (bgdFrac > (1.0 - DBL_EPSILON)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(refObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dim = 2;
        lCon = WLZ_8_CONNECTED;
	break;
      case WLZ_3D_DOMAINOBJ:
	dim = 3;
        lCon = WLZ_26_CONNECTED;
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Compute threshold value and type from histogram. */
  if(errNum == WLZ_ERR_NONE)
  {
    hObj = WlzAssignObject(
    	   WlzHistogramObj(ppObj, 256, 0.0, 1.0, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzHistogramCnvGauss(hObj, sigma, 0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCompThresholdVT(hObj, compThrMethod, bgdFrac, 0.0,
    			        0.0, &tV, &tType);
  }
  (void )WlzFreeObj(hObj); hObj = NULL;
  /* Threshold object. */
  if(errNum == WLZ_ERR_NONE)
  {
    tObj = WlzAssignObject(
     	   WlzThreshold(ppObj, tV, tType, &errNum), NULL);
  }
  /* Label to get connected components. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzLabel(tObj, &(split.nLComp), &(split.lComp), maxComp, 0, lCon);
  }
  /* Sort connected components by size. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(split.nLComp < nReqComp)
    {
      nReqComp = split.nLComp;
    }
    if(((split.compI = (int *)AlcMalloc(sizeof(int) *
                                        split.nLComp)) == NULL) ||
       ((split.lCompSz = (WlzLong *)AlcMalloc(sizeof(WlzLong) *
                                              split.nLComp)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idC = 0; 
    while((errNum == WLZ_ERR_NONE) && (idC < split.nLComp))
    {
      split.compI[idC] = idC;
      split.lCompSz[idC] = (dim == 2)? WlzArea(split.lComp[idC], &errNum):
      				       WlzVolume(split.lComp[idC], &errNum);
      ++idC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Sort component indices by component size. */
    AlgQSort(split.compI, split.nLComp, sizeof(int), &split,
    	     WlzSplitObjSortSzFn);
    /* Allocate array for cliped component objects. */
    if((comp = (WlzObject **)AlcCalloc(sizeof(WlzObject *),
    				       split.nLComp)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute bounding box and clip objects from the reference object. */
  if(errNum == WLZ_ERR_NONE)
  {
    idC = 0;
    while((errNum == WLZ_ERR_NONE) && (idC < nReqComp))
    {
      if(dim == 2)
      {
        box.i2 = WlzBoundingBox2I(split.lComp[split.compI[idC]], &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  box.i2.xMin -= bWidth;
	  box.i2.yMin -= bWidth;
	  box.i2.xMax += bWidth;
	  box.i2.yMax += bWidth;
	  comp[idC] = WlzClipObjToBox2D(refObj, box.i2, &errNum);
	}
      }
      else /* dim == 3 */
      {
        box.i3 = WlzBoundingBox3I(split.lComp[split.compI[idC]], &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  box.i3.xMin -= bWidth;
	  box.i3.yMin -= bWidth;
	  box.i3.zMin -= bWidth;
	  box.i3.xMax += bWidth;
	  box.i3.yMax += bWidth;
	  box.i3.zMax += bWidth;
	  comp[idC] = WlzClipObjToBox3D(refObj, box.i3, &errNum);
	}
      }
      ++idC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNComp = nReqComp;
    *dstComp = comp;
  }
  /* Free temporary storage. */
  if(split.lComp)
  {
    for(idC = 0; idC < split.nLComp; ++idC)
    {
      (void )WlzFreeObj(split.lComp[idC]);
    }
    AlcFree(split.lComp);
  }
  AlcFree(split.compI);
  AlcFree(split.lCompSz);
  (void )WlzFreeObj(tObj);
  return(errNum);
}

/*!
* \return	Result of comparison.
* \ingroup	WlzBinaryOps
* \brief	Sort function for AlgQSort() to sort components
*		by size (area or volume).
* \param	cData			Used to pass the split image data.
* \param	a			Used to pass first component.
* \param	b			Used to pass second component.
*/
static int	WlzSplitObjSortSzFn(const void *cData,
				    const void *a, const void *b)
{
  int		cmp;
  WlzLong	dif;
  WlzSplitObjData *split;

  split = (WlzSplitObjData *)cData;
  dif = *(split->lCompSz + *(int *)b) - *(split->lCompSz + *(int *)a);
  cmp = (dif > 0) - (dif < 0);
  return(cmp);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinaryOps
* \brief	Splits the given montage object into component objects
*		clipped from the montage object.  The montage object
*		must be composed of component images embedded in a
*		background, with little variation in the background
*		values.
* \param	mObj			Montage object, which must be either
*					a WLZ_2D_DOMAINOBJ or a
*					WLZ_3D_DOMAINOBJ with values.
* \param	gapV			Value for the uniform background.
*					Must be either WLZ_GREY_INT or
*					WLZ_GREY_RGBA.
* \param	tol			Tolerance (fraction) for the
*					variation in background values.
* \param	bWidth			Additional boundary width added
*					to detected images before they are
*					clipped.
* \param	minArea			Minimum area for a valid component
*					image, must be greater than zero.
* \param	maxComp			Maximum number of components.
* \param	dstNComp		Destination pointer for the number of
*					components extracted, must not be NULL.
* \param	dstComp			Destination pointer for the extracted
*					components, must not be NULL.
*/
WlzErrorNum 	WlzSplitMontageObj(WlzObject *mObj, WlzPixelV gapV,
				      double tol, int bWidth,
				      WlzLong minArea, int maxComp,
				      int *dstNComp, WlzObject ***dstComp)
{
  int		id0,
  		id1,
		nLComp = 0;
  WlzLong	area = 0;
  WlzObject	*gObj = NULL,
  		*tObj = NULL;
  WlzObject	**lComp;
  WlzGreyType	objG;
  WlzBox	box;
  WlzPixelV	gapLV,
  		gapHV;
  WlzConnectType lCon;
  int		tI[8];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tol = WLZ_CLAMP(tol, 0.0, 1.0);
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(minArea < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(mObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	lCon = WLZ_4_CONNECTED;
        break;
      case WLZ_3D_DOMAINOBJ:
	lCon = WLZ_6_CONNECTED;
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    objG = WlzGreyTypeFromObj(mObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gapV.type)
    {
      case WLZ_GREY_INT: /* FALLTHROUGH */
      case WLZ_GREY_RGBA:
        break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(objG == WLZ_GREY_RGBA)
    {
      if(gapV.type != WLZ_GREY_RGBA)
      {
        (void )WlzValueConvertPixel(&gapV, gapV, WLZ_GREY_RGBA);
      }
    }
    else
    {
      if(gapV.type != WLZ_GREY_INT)
      {
        (void )WlzValueConvertPixel(&gapV, gapV, WLZ_GREY_INT);
      }
    }
    gapLV.type = gapHV.type = gapV.type;
    if(gapV.type == WLZ_GREY_INT)
    {
      tI[0] = gapV.v.inv * tol;
      gapLV.v.inv = gapV.v.inv - tI[0];
      gapHV.v.inv = gapV.v.inv + tI[0];
      tObj = WlzThreshold(mObj, gapLV, WLZ_THRESH_HIGH, &errNum);
      if((errNum == WLZ_ERR_NONE) && (tObj != NULL))
      {
	gObj = WlzThreshold(tObj, gapHV, WLZ_THRESH_LOW, &errNum);
      }
      (void )WlzFreeObj(tObj);
      tObj = NULL;
    }
    else /* gapV.type == WLZ_GREY_RGBA */
    {
	tI[0] = WLZ_RGBA_RED_GET(gapV.v.rgbv);
	tI[1] = (int )floor((double )(tI[0]) * tol);
	tI[2] = tI[0] - tI[1];
	tI[5] = tI[0] + tI[1];
	tI[0] = WLZ_RGBA_GREEN_GET(gapV.v.rgbv);
	tI[1] = (int )floor((double )(tI[0]) * tol);
	tI[3] = tI[0] - tI[1];
	tI[6] = tI[0] + tI[1];
	tI[0] = WLZ_RGBA_BLUE_GET(gapV.v.rgbv);
	tI[1] = (int )floor((double )(tI[0]) * tol);
	tI[4] = tI[0] - tI[1];
	tI[7] = tI[0] + tI[1];
	tI[2] = WLZ_CLAMP(tI[2], 0, 255);
	tI[3] = WLZ_CLAMP(tI[3], 0, 255);
	tI[4] = WLZ_CLAMP(tI[4], 0, 255);
	WLZ_RGBA_RGBA_SET(gapLV.v.rgbv, tI[2], tI[3], tI[4], 255);
	tI[5] = WLZ_CLAMP(tI[5], 0, 255);
	tI[6] = WLZ_CLAMP(tI[6], 0, 255);
	tI[7] = WLZ_CLAMP(tI[7], 0, 255);
	WLZ_RGBA_RGBA_SET(gapHV.v.rgbv, tI[5], tI[6], tI[7], 255);
        gObj = WlzRGBABoxThreshold(mObj, gapLV, gapHV, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj = WlzDiffDomain(mObj, gObj, &errNum);
  }
  (void )WlzFreeObj(gObj);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzLabel(tObj, &nLComp, &lComp, maxComp, 0, lCon);
  }
  (void )WlzFreeObj(tObj);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Get rid of small objects using minArea as the threshold. */
    id0 = 0;
    id1 = 0;
    while(id0 < nLComp)
    {
      switch((*(lComp + id0))->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  area = WlzArea(*(lComp + id0), NULL);
	  break;
        case WLZ_3D_DOMAINOBJ:
	  area = WlzVolume(*(lComp + id0), NULL);
	  break;
        default:
          area = 0;
	  break;
      }
      if(area >= minArea)
      {
        *(lComp + id1) = *(lComp + id0);
        ++id1;
      }
      else
      {
        (void )WlzFreeObj(*(lComp + id0));
	*(lComp + id0) = NULL;
      }
      ++id0;
    }
    nLComp = id1;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Clip rectangular objects from the montage object. */
    id0 = 0;
    while((errNum == WLZ_ERR_NONE) && (id0 < nLComp))
    {
      if(tObj->type == WLZ_2D_DOMAINOBJ)
      {
        box.i2 = WlzBoundingBox2I(*(lComp + id0), &errNum);
	box.i2.xMin -= bWidth;
	box.i2.yMin -= bWidth;
	box.i2.xMax += bWidth;
	box.i2.yMax += bWidth;
	(void )WlzFreeObj(*(lComp + id0));
	*(lComp + id0) = WlzClipObjToBox2D(mObj, box.i2, &errNum);
      }
      else /* tObj->type == WLZ_3D_DOMAINOBJ */
      {
        box.i3 = WlzBoundingBox3I(*(lComp + id0), &errNum);
	box.i3.xMin -= bWidth;
	box.i3.yMin -= bWidth;
	box.i3.zMin -= bWidth;
	box.i3.xMax += bWidth;
	box.i3.yMax += bWidth;
	box.i3.zMax += bWidth;
	(void )WlzFreeObj(*(lComp + id0));
	*(lComp + id0) = WlzClipObjToBox3D(mObj, box.i3, &errNum);
      }
      ++id0;
    }
  }
  *dstNComp = nLComp;
  *dstComp = lComp;
  return(errNum);
}
