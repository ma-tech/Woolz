#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzHyThreshold.c
* Date:         May 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A hysteresis threshold filter for Woolz.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzHyThreshold
* Returns:	WlzObject:		Thresholded object or
*					NULL on error.
* Purpose:	Hysteresis thresholds the given Woolz object.
*		Values are in the domain of the hysteresis threshold'd
*		object if they are above/below the primary threshold
*		or above/below the secondary threshold and connected
*		to values above/below the primary threshold.
*		hilo se
*		
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Object to be thresholded.
*		WlzPixelV *pThrV:	Primary hysteresis threshold
*					value
*		WlzPixelV *sThrV:	Secondary hysteresis threshold
*					value
*		WlzThresholdType hilo:	Threshold for above or below
*					values.
*		WlzConnectType con:	Connectivity to examine for
*					hysteresis.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzObject	*WlzHyThreshold(WlzObject *srcObj,
			       WlzPixelV pThrV, WlzPixelV sThrV,
			       WlzThresholdType hilo, WlzConnectType con,
			       WlzErrorNum *dstErr)
{
  int		simpleThr = 0;
  WlzPixelV	tmpV;
  WlzObject	*dstObj = NULL,
  		*pThrObj = NULL,
  		*sThrObj = NULL,
		*dThrObj = NULL,
		*iThrObj = NULL,
		*uThrObj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    if(con == WLZ_0_CONNECTED)
    {
      simpleThr = 1;
    }
    else
    {
      if((errNum = WlzValueConvertPixel(&tmpV, sThrV,
      					pThrV.type)) == WLZ_ERR_NONE)
      {
        switch(pThrV.type)
	{
	  case WLZ_GREY_INT:
	    if(tmpV.v.inv == pThrV.v.inv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(tmpV.v.shv == pThrV.v.shv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(tmpV.v.ubv == pThrV.v.ubv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(fabs(tmpV.v.flv - pThrV.v.flv) <= FLT_EPSILON)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(fabs(tmpV.v.dbv - pThrV.v.dbv) <= DBL_EPSILON)
	    {
	      simpleThr = 1;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(simpleThr)
    {
      dstObj = WlzThreshold(srcObj, pThrV, hilo, &errNum);
    }
    else
    {
      pThrObj = WlzThreshold(srcObj, pThrV, hilo, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	if(pThrObj->type == WLZ_EMPTY_OBJ)
	{
	  dstObj = pThrObj;
	  pThrObj = NULL;
	}
	else
	{
	  sThrObj = WlzThreshold(srcObj, sThrV, hilo, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(sThrObj->type == WLZ_EMPTY_OBJ)
	    {
	      dstObj = pThrObj;
	      pThrObj = NULL;
	    }
	    else
	    {
	      dThrObj = WlzDilation(pThrObj, con, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		iThrObj = WlzIntersect2(dThrObj, sThrObj, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		uThrObj = WlzUnion2(pThrObj, iThrObj, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		dstObj = WlzMakeMain(uThrObj->type,
				     uThrObj->domain, srcObj->values,
				     srcObj->plist, srcObj, &errNum);
	      }
	      if(dThrObj)
	      {
		WlzFreeObj(dThrObj);
	      }
	      if(iThrObj)
	      {
		WlzFreeObj(iThrObj);
	      }
	      if(uThrObj)
	      {
		WlzFreeObj(uThrObj);
	      }
	    }
	    if(sThrObj)
	    {
	      WlzFreeObj(sThrObj);
	    }
	    if(pThrObj)
	    {
	      WlzFreeObj(pThrObj);
	    }
	  }
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}
