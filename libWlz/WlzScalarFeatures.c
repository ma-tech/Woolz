#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzScalarFeatures.c
* \author       Bill Hill
* \date         November 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for extracting scalar features from objects.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Finds scalar features within the given object. Where the
*		feature values are all either above or below the given
*		threshold value, depending on the value of thrHL, and
*		the features have the given minimum seperation distance.
*
*		The features are found using the following algorithm:
*		<ul>
*		  <li>
*		  Compute feature values from given object.
*		  </li>
*		  <li>
*		  Threshold and extract a list of features.
*		  </li>
*		  <li>
*		  Sort the list of features by value.
*		  </li>
*		  <li>
*		  Create a kD-tree.
*		  </li>
*		  <li>
*		  While the list of feature values is not empty.
*		  <ul>
*		    <li>
*		    Remove item from list of features.
*		    </li>
*		    <li>
*		    Find distance to nearest neighbour in the kD-tree.
*		    </li>
*		    <li>
*		    If the distance is greater than the minimum.
*		    <ul>
*		      <li>
*		      Add the feature to the kD-tree.
*		      </li>
*		    </ul></li>
*		  </ul></li>
*		  <li>
*		  Output features of the kD-tree.
*		  </li>
*		</ul>
*		The cost of this function depends on the threshold value
*		and to some extent on the minimum distance.
*
*		If feature values are computed then they are normalised
*		using WlzGreyNormalise(), to allow threshold values to be
*		set.
* \param	obj			Given object.
* \param	fType			Type of feature to find.
* \param	thrHL			High or low feature values.
* \param	thrV			Threshold feature value.
* \param	minDist			Minimum distance between features.
* \param	dstNFeat		Destination pointer for the number
*					of scalar features found.
* \param	dstFeat			Destination pointer for coordinates
*					of the features found.
* \param	filterV			Filter value for computing features.
* \param	dstFeat			Destination pointer for the
*					coordinates of the scalar features.
*/
WlzErrorNum	WlzScalarFeatures2D(WlzObject *obj,
				    int *dstNFeat,
				    WlzIVertex2 **dstFeat,
				    WlzScalarFeatureType fType,
				    WlzThresholdType thrHL, WlzPixelV thrV,
				    double filterV, double minDist)
{
  int		tI0,
  		cIdx,
  		fCnt,
  		rIdx,
  		nVal = 0;
  double	dist;
  WlzGreyType	vType;
  WlzGreyP	valP;
  WlzPixelV	minV,
  		maxV;
  int		*vRank = NULL,
  		*dRank = NULL;
  WlzObject	*fObj = NULL,
  		*tObj = NULL;
  WlzVertexP	vCoords,
  		fCoords;
  AlcKDTTree	*tree = NULL;
  AlcKDTNode	*node;
  int 		(*sortFn)(void *, int *, int, int);
  double	pos[2];
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  valP.inp = NULL;
  vCoords.v = fCoords.v = NULL;
  if(minDist < (1.0 - DBL_EPSILON))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* Compute object with feature values from given object. */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(fType)
    {
      case WLZ_SCALARFEATURE_VALUE:
	fObj = WlzAssignObject(obj, NULL);
        break;
      case WLZ_SCALARFEATURE_GRADIENT:
        fObj = WlzAssignObject(
	       WlzGreyModGradient(obj, filterV, &errNum) ,NULL);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzGreyNormalise(fObj);
	}
        break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
        break;
    }
  }
  /* Threshold the feature values object and extract a list of features. */
  if(errNum == WLZ_ERR_NONE)
  {
    tObj = WlzAssignObject(WlzThreshold(fObj, thrV, thrHL, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((tI0 = (int )floor((minDist / 3.0) + DBL_EPSILON)) < 1)
    {
      tI0 = 1;
    }
    errNum = WlzSampleValuesAndCoords(tObj, &vType, &nVal, &valP, &vCoords,
    				(thrHL == WLZ_THRESH_LOW)?
				WLZ_SAMPLEFN_MIN: WLZ_SAMPLEFN_MAX, tI0);
  }
  if((errNum == WLZ_ERR_NONE) && (nVal > 0))
  {
    if(((vRank = (int *)AlcMalloc(nVal * sizeof(int))) == NULL) ||
       ((dRank = (int *)AlcMalloc(nVal * sizeof(int))) == NULL))
    {
      errNum == WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      for(rIdx = 0; rIdx < nVal; ++rIdx)
      {
        *(vRank + rIdx) = rIdx;
      }
    }
  }
  /* Clear up the temporary objects. */
  (void )WlzFreeObj(fObj);
  (void )WlzFreeObj(tObj);
  if((errNum == WLZ_ERR_NONE) && (nVal > 0))
  {
    /* Compute rank indices for sorted grey values. */
    switch(vType)
    {
      case WLZ_GREY_LONG:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxLFn: AlgHeapSortCmpIdxLFn;
	break;
      case WLZ_GREY_INT:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxIFn: AlgHeapSortCmpIdxIFn;
	break;
      case WLZ_GREY_SHORT:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxSFn: AlgHeapSortCmpIdxSFn;
	break;
      case WLZ_GREY_UBYTE:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxUFn: AlgHeapSortCmpIdxUFn;
	break;
      case WLZ_GREY_FLOAT:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxFFn: AlgHeapSortCmpIdxFFn;
	break;
      case WLZ_GREY_DOUBLE:
        sortFn = (thrHL == WLZ_THRESH_HIGH)?
		 AlgHeapSortInvCmpIdxDFn: AlgHeapSortCmpIdxDFn;
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )AlgHeapSortIdx(valP.inp, vRank, nVal, sortFn);
    }
  }
  /* Clear up the values since only their rank is needed. */
  AlcFree(valP.inp);
  /* Create a kD-tree and then try adding the coordinates maintaining
   * the required minimum distance. */
  if((errNum == WLZ_ERR_NONE) && (nVal > 0))
  {
    if((tree = AlcKDTTreeNew(ALC_POINTTYPE_DBL, 2, 0.0, nVal, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      /* Insert highest ranked value into the kD-tree. */
      fCnt = 0;
      rIdx = 0;
      cIdx = *(vRank + 0);
      *(dRank + fCnt++) = cIdx;
      pos[0] = (vCoords.i2 + cIdx)->vtX;
      pos[1] = (vCoords.i2 + cIdx)->vtY;
      (void )AlcKDTInsert(tree, pos, NULL, &alcErr);
      rIdx = 1;
      while((alcErr == ALC_ER_NONE) && (rIdx < nVal))
      {
	cIdx = *(vRank + rIdx);
	pos[0] = (vCoords.i2 + cIdx)->vtX;
	pos[1] = (vCoords.i2 + cIdx)->vtY;
	node = AlcKDTGetNN(tree, pos, INT_MAX, &dist, NULL);
	if(dist > minDist)
	{
	  *(dRank + fCnt++) = cIdx;
	  (void )AlcKDTInsert(tree, pos, NULL, &alcErr);
	}
        ++rIdx;
      }
      (void )AlcKDTTreeFree(tree);
      if(alcErr != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  /* Value rank table no longer needed. */
  AlcFree(vRank);
  if((errNum == WLZ_ERR_NONE) && (fCnt > 0))
  {
    if((fCoords.v = AlcMalloc(fCnt * sizeof(WlzIVertex2))) == NULL)
    {
      errNum = WLZ_ERR_NONE;
    }
    else
    {
      for(rIdx = 0; rIdx < fCnt; ++rIdx)
      {
        *(fCoords.i2 + rIdx) = *(vCoords.i2 + *(dRank + rIdx));
      }
    }
  }
  if(errNum == ALC_ER_NONE)
  {
    *dstNFeat = fCnt;
    *dstFeat = fCoords.i2;
  }
  AlcFree(dRank);
  AlcFree(vCoords.v);
  return(errNum);
}
