#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzSplitObj.c
* \author       Bill Hill
* \date         November 2004
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions to split a single object into component objects.
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/
#include <float.h>
#include <Wlz.h>

typedef struct _WlzSplitObjData
{
  int		nLComp;
  int		*compI;
  int		*lCompSz;
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
*					0 - 255 as UBYTE greys.
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
  int		tI,
  		dim,
  		idC;
  WlzObject	*hObj = NULL,
  		*tObj = NULL;
  WlzObject 	**comp = NULL;
  WlzBox	box;
  WlzPixelV	tV;
  WlzSplitObjData split;
  WlzThresholdType tType;
  WlzConnectType lCon;
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
       ((split.lCompSz = (int *)AlcMalloc(sizeof(int) *
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
  int		idx,
  		aSz,
  		bSz,
		cmp;
  WlzSplitObjData *split;

  split = (WlzSplitObjData *)cData;
  idx = *(int *)a;
  aSz = *(split->lCompSz + idx);
  idx = *(int *)b;
  bSz = *(split->lCompSz + idx);
  cmp = bSz - aSz;
  return(cmp);
}
