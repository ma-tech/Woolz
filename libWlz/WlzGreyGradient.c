#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyGradient.c
* Date:         May 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes a new grey valued object where the grey values
*		are the gradient of the gray values in the original
*		image.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzGreyGradient2D
* Returns:	WlzObject:		New Woolz domain object with
*					gradient grey values or NULL
*					on error.
* Purpose:	Computes the gradient of the gray values in the
*		given Woolz 2D domain object.
* Global refs:	-
* Parameters:	WlzObject **dstGrdY:	Destination pointer for the
*					gradient (partial derivative)
*					through lines, may be NULL.
*		WlzObject **dstGrdX:	Destination pointer for the
*					gradient (partial derivative)
*					through columns, may be NULL.
*		WlzObject *srcObj:	Given source object.
*		WlzRsvFilter *flt:	Recursive filter to use.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzObject *WlzGreyGradient2D(WlzObject **dstGrdY, WlzObject **dstGrdX,
				    WlzObject *srcObj,
				    WlzRsvFilter *flt, WlzErrorNum *dstErr)
{
  WlzObject	*grdX,
  		*grdY,
		*dstObj = NULL;
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
  if(errNum == WLZ_ERR_NONE)
  {
    grdX = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_X, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grdY = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_Y, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzImageArithmetic(grdY, grdX, WLZ_MAGNITUDE, 0, &errNum);
  }
  if(dstGrdX && (errNum == WLZ_ERR_NONE))
  {
    *dstGrdX = grdX;
  }
  else if(grdX)
  {
    WlzFreeObj(grdX);
  }
  if(dstGrdY && (errNum == WLZ_ERR_NONE))
  {
    *dstGrdY = grdY;
  }
  else if(grdY)
  {
    WlzFreeObj(grdY);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzGreyGradient3D
* Returns:	WlzObject:		New Woolz domain object with
*					gradient grey values or NULL
*					on error.
* Purpose:	Computes the gradient of the gray values in the
*		given Woolz 3D domain object.
* Global refs:	-
* Parameters:	WlzObject **dstGrdZ:	Destination pointer for the
*					gradient (partial derivative)
*					through planes, may be NULL.
*		WlzObject **dstGrdY:	Destination pointer for the
*					gradient (partial derivative)
*					through lines, may be NULL.
*		WlzObject **dstGrdX:	Destination pointer for the
*					gradient (partial derivative)
*					through columns, may be NULL.
*		WlzObject *srcObj:	Given source object.
*		WlzRsvFilter *flt:	Recursive filter to use.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzObject *WlzGreyGradient3D(WlzObject **dstGrdZ, WlzObject **dstGrdY,
				    WlzObject **dstGrdX, WlzObject *srcObj,
				    WlzRsvFilter *flt, WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WLZ_ERR_OBJECT_TYPE; /* TODO implement this function */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/************************************************************************
* Function:	WlzGreyGradient
* Returns:	WlzObject:		New Woolz domain object with
*					gradient grey values or NULL
*					on error.
* Purpose:	Computes the gradient of the gray values in the
*		given Woolz domain object.
* Global refs:	-
* Parameters:	WlzObject **dstGrdZ:	Destination pointer for the
*					gradient (partial derivative)
*					through planes, may be NULL.
*		WlzObject **dstGrdY:	Destination pointer for the
*					gradient (partial derivative)
*					through lines, may be NULL.
*		WlzObject **dstGrdX:	Destination pointer for the
*					gradient (partial derivative)
*					through columns, may be NULL.
*		WlzObject *srcObj:	Given source object.
*		WlzRsvFilter *flt:	Recursive filter to use.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzObject	*WlzGreyGradient(WlzObject **dstGrdZ, WlzObject **dstGrdY,
				 WlzObject **dstGrdX, WlzObject *srcObj,
				 WlzRsvFilter *flt, WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(flt == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dstObj = WlzGreyGradient2D(dstGrdY, dstGrdX, srcObj,
				   flt, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
        dstObj = WlzGreyGradient3D(dstGrdZ, dstGrdY, dstGrdX, srcObj,
				   flt, &errNum);
        break;
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
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
  return(dstObj);
}
