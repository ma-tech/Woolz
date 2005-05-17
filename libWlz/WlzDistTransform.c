#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzDistTransform.c
* \author       Konstantinos Liakos, Bill Hill
* \date         October 2004
* \version      $Id$
* \note
*               Copyright
*               2004 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Distance transform functions which calculate the distance of
* 		every pixel in a foreground object from a reference object.
* \ingroup	WlzMorphologyOps
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Distance object which shares the given foreground object's
*		domain and has integer distance values, null on error.
* \ingroup	WlzMorphologyOps
* \brief	Computes the distance of every pixel/voxel in the foreground
* 		object from the reference object.
*		For octagonal distance: 4 and 8 connectivities are alternated
*		in 2D and 6 and 26 connectivities are alternated in 3D. See:
*		G. Borgefors. "Distance Transformations in Arbitrary
*		Dimensions" CVGIP 27:321-345, 1984.
* \param	forObj			Foreground object.
* \param	refObj			Reference object.
* \param	dFn			Distance function which must be
*					appropriate to the dimension of
*					the foreground and reference objects.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*DistanceTransform(WlzObject *forObj, WlzObject *refObj,
				   WlzDistanceType dFn, WlzErrorNum *dstErr)
{
  int 		idP,
		lastP,
		dim,
  		notDone = 1;
  double	dInc = 1.0;
  WlzObject	*dilObj,
  		*dstObj,
		*difObj,
		*prvItrObj,
  		*curItrObj = NULL;
  WlzObject 	*bothObj[2];
  WlzDomain	*difDoms;
  WlzPixelV 	dstV,
  		bgdV;
  WlzValues 	*difVals;
  WlzConnectType con;
  WlzObjectType dstGType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzValues 	difVal,
  		dstVal;

  if((forObj == NULL) || (refObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((forObj->domain.core == NULL) || (refObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  /* Create new values for the computed distances. */
  if(errNum == WLZ_ERR_NONE)
  {
    bgdV.type = WLZ_GREY_INT;
    bgdV.v.ubv = 0;
    dstV.type = WLZ_GREY_DOUBLE;
    dstV.v.dbv = 0.0;
    dstGType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_SHORT, NULL);
    switch(forObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	switch(dFn)
	{
	  case WLZ_4_DISTANCE: /* FALLTHROUGH */
	  case WLZ_8_DISTANCE: /* FALLTHROUGH */
	  case WLZ_OCTAGONAL_DISTANCE:
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dim = 2;
          dstVal.v = WlzNewValueTb(forObj, dstGType, bgdV, &errNum);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	switch(dFn)
	{
	  case WLZ_6_DISTANCE:  /* FALLTHROUGH */
	  case WLZ_18_DISTANCE: /* FALLTHROUGH */
	  case WLZ_26_DISTANCE: /* FALLTHROUGH */
	  case WLZ_OCTAGONAL_DISTANCE:
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dim = 3;
	  dstVal.vox = WlzNewValuesVox(forObj, dstGType, bgdV, &errNum);
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(dFn)
    {
      case WLZ_4_DISTANCE:
        con = WLZ_4_CONNECTED;
	break;
      case WLZ_6_DISTANCE:
        con = WLZ_6_CONNECTED;
	break;
      case WLZ_8_DISTANCE:
        con = WLZ_8_CONNECTED;
	break;
      case WLZ_18_DISTANCE:
        con = WLZ_18_CONNECTED;
	break;
      case WLZ_26_DISTANCE:
        con = WLZ_26_CONNECTED;
	break;
      case WLZ_OCTAGONAL_DISTANCE:
        con = (dim == 2)? WLZ_8_DISTANCE: WLZ_26_CONNECTED;
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  /* Create a distance object using the foreground object's domain and
   * the new distance values. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzMakeMain(forObj->type, forObj->domain, dstVal,
			 NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bothObj[0] = forObj;
    errNum = WlzGreySetValue(dstObj, dstV);
  }
  /* Dilate the reference object while setting the distances in each
   * dilated shell. */
  while((errNum == WLZ_ERR_NONE) && notDone)
  {
    dstV.v.dbv += dInc;
    prvItrObj = curItrObj;
    dilObj = WlzDilation(refObj, con, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(forObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  curItrObj = WlzIntersect2(dilObj, forObj, &errNum);
	  break;
        case WLZ_3D_DOMAINOBJ:
	  bothObj[1] = dilObj;
	  curItrObj = WlzIntersectN(2, bothObj, 1, &errNum);
	  break;
      }
      (void)WlzFreeObj(dilObj);
    }
    /* Create difference object for the expanding shell. */
    if(errNum == WLZ_ERR_NONE)
    {
      difObj = WlzDiffDomain(curItrObj, refObj, &errNum);
    }
    if((difObj == NULL) || WlzIsEmpty(difObj, &errNum))
    {
      notDone = 0;
    }
    else
    {
      /* Assign the distance object's values to the difference object
       * and set all it's values to the current distance. */
      if(errNum == WLZ_ERR_NONE)
      {
	switch(forObj->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    difObj->values = WlzAssignValues(dstObj->values, NULL);
	    errNum = WlzGreySetValue(difObj, dstV);
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    /* 3D is more complex than 2D: Need to create a temporary
	     * voxel valuetable and assign the individual 2D values. */
	    difVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					     difObj->domain.p->plane1,
					     difObj->domain.p->lastpl,
					     bgdV, NULL, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      difObj->values = WlzAssignValues(difVal, NULL);
	      difDoms = difObj->domain.p->domains;
	      difVals = difObj->values.vox->values;
	      idP = difObj->domain.p->plane1;
	      lastP = difObj->domain.p->lastpl;
	      while(idP <= lastP)
	      {
		if((*difDoms).core)
		{
		  dstVal = dstObj->values.vox->values[idP - 
						      dstObj->domain.p->plane1];
		  *difVals = WlzAssignValues(dstVal, NULL);
		}
		++idP;
		++difDoms;
		++difVals;
	      }
	      if(difObj->domain.p->lastpl > difObj->domain.p->plane1)
	      {
		errNum = WlzGreySetValue(difObj, dstV);
	      }
	    }
	    break;
	}
      }
      refObj = curItrObj;
    }
    (void )WlzFreeObj(difObj); difObj = NULL;
    (void)WlzFreeObj(prvItrObj); prvItrObj = NULL;
    /* Alternate connectivities for octagonal distance. */
    if(dFn == WLZ_OCTAGONAL_DISTANCE)
    {
      if(dim == 2)
      {
        con = (con == WLZ_4_DISTANCE)? WLZ_8_DISTANCE: WLZ_4_DISTANCE;
      }
      else /* dim == 3 */
      {
        con = (con == WLZ_6_DISTANCE)? WLZ_26_DISTANCE: WLZ_6_DISTANCE;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(dstObj); dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}
