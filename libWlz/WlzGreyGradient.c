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
* 05-06-2000 bill Removed unused variables.
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>

static WlzObject *WlzGreyGradient2D(WlzObject **dstGrdY, WlzObject **dstGrdX,
			WlzObject *srcObj, WlzRsvFilter *flt,
			WlzErrorNum *dstErr);
static WlzObject *WlzGreyGradient3D(WlzObject **dstGrdZ, WlzObject **dstGrdY,
			WlzObject **dstGrdX, WlzObject *srcObj,
			WlzRsvFilter *flt, WlzErrorNum *dstErr);
static 	WlzObject *WlzGreyMagnitude3D(WlzObject *srcObj0, WlzObject *srcObj1,
			WlzObject *srcObj2, WlzErrorNum *dstErr);
static WlzObject *WlzGreyMagnitude2D3(WlzObject *srcObj0, WlzObject *srcObj1,
			WlzObject *srcObj2, WlzErrorNum *dstErr);
static void	WlzBufMagD3(double *dP0, double *dP1, double *dP2, int cnt);

/************************************************************************
* Function:	WlzBufMagD3
* Returns:	void
* Purpose:	Computes the magnitudes of the 3 sets of components
*		given in the double vectors and returns the magnitudes
*		in the first vector.
* Global refs:	-
* Parameters:	double *dP0:		First vector of components, also
*					used to return magnitudes.
*		double *dP1:		Second vector of components.
*		double *dP2:		Third vector of components.
*		int cnt:		Number of components.
************************************************************************/
static void	WlzBufMagD3(double *dP0, double *dP1, double *dP2, int cnt)
{
  double	tD0;

  while(cnt-- > 0)
  {
    if((tD0 = (*dP0 * *dP0) + (*dP1 * *dP1) + (*dP2 *  *dP2)) > DBL_EPSILON)
    {
      *dP0 = sqrt(tD0);
    }
    else
    {
      *dP0 = 0.0;
    }
    ++dP0;
    ++dP1;
    ++dP2;
  }
}

/************************************************************************
* Function:	WlzGreyGradient2D
* Returns:	WlzObject:		New Woolz domain object with
*					gradient grey values or NULL
*					on error.
* Purpose:	Computes the magnitude of the gray values in the
*		3 given Woolz 2D domain objects.
* Global refs:	-
* Parameters:	WlzObject *srcObj0:	First object.
*		WlzObject *srcObj1:	Second object.
*		WlzObject *srcObj2:	Third object.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzObject *WlzGreyMagnitude2D3(WlzObject *srcObj0, WlzObject *srcObj1,
				      WlzObject *srcObj2,
				      WlzErrorNum *dstErr)
{
  int		idN,
		itvLen,
  		bufSz;
  double	**iBufA = NULL;
  WlzObject     *tObj0,
		*istObj = NULL,
                *dstObj = NULL;
  WlzObject	*iObjA[3],
  		*tObjA[3];
  WlzGreyP	tGP0;
  WlzGreyType	dstGType;
  WlzGreyType	gTypeA[3];
  WlzPixelV	dstBgd;
  WlzPixelV	bgdA[3];
  WlzValues	dstVal;
  WlzIntervalWSpace dstIWSp;
  WlzGreyWSpace	dstGWSp;
  WlzIntervalWSpace iIWSpA[3];
  WlzGreyWSpace	iGWSpA[3];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  *(iObjA + 0) = *(iObjA + 1) = *(iObjA + 2) = NULL;
  /* Check source objects. */
  if((srcObj0 == NULL) || (srcObj1 == NULL) ||(srcObj2 == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;;
  }
  else if((srcObj0->type != WLZ_2D_DOMAINOBJ) ||
          (srcObj1->type != WLZ_2D_DOMAINOBJ) ||
          (srcObj2->type != WLZ_2D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((srcObj0->domain.core == NULL) ||
          (srcObj1->domain.core == NULL) ||
          (srcObj2->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((srcObj0->values.core == NULL) ||
          (srcObj1->values.core == NULL) ||
          (srcObj2->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  /* Compute the intersection of the source objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    *(tObjA + 0) = srcObj0;
    *(tObjA + 1) = srcObj1;
    *(tObjA + 2) = srcObj2;
    istObj = WlzIntersectN(3, tObjA, 0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(iObjA + 0) = WlzMakeMain(WLZ_2D_DOMAINOBJ,
    			       istObj->domain, srcObj0->values,
    			       NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(iObjA + 1) = WlzMakeMain(WLZ_2D_DOMAINOBJ,
    			       istObj->domain, srcObj1->values,
    			       NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(iObjA + 2) = WlzMakeMain(WLZ_2D_DOMAINOBJ,
    			       istObj->domain, srcObj2->values,
    			       NULL, NULL, &errNum);
  }
  /* Get background value and grey types */
  idN = 0;
  while((errNum == WLZ_ERR_NONE) && (idN < 3))
  {
    tObj0 = *(iObjA + idN);
    *(gTypeA + idN) = WlzGreyTableTypeToGreyType(tObj0->values.core->type,
				                 &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      *(bgdA + idN) = WlzGetBackground(tObj0, &errNum);
    }
    ++idN;
  }
  /* Promote grey types. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((*(gTypeA + 0) == WLZ_GREY_DOUBLE) ||
       (*(gTypeA + 1) == WLZ_GREY_DOUBLE) ||
       (*(gTypeA + 2) == WLZ_GREY_DOUBLE))
    {
      dstGType = WLZ_GREY_DOUBLE;
    }
    else if((*(gTypeA + 0) == WLZ_GREY_FLOAT) ||
	    (*(gTypeA + 1) == WLZ_GREY_FLOAT) ||
	    (*(gTypeA + 2) == WLZ_GREY_FLOAT))
    {
      dstGType = WLZ_GREY_FLOAT;
    }
    else if((*(gTypeA + 0) == WLZ_GREY_INT) ||
	    (*(gTypeA + 1) == WLZ_GREY_INT) ||
	    (*(gTypeA + 2) == WLZ_GREY_INT))
    {
      dstGType = WLZ_GREY_INT;
    }
    else if((*(gTypeA + 0) == WLZ_GREY_SHORT) ||
	    (*(gTypeA + 1) == WLZ_GREY_SHORT) ||
	    (*(gTypeA + 2) == WLZ_GREY_SHORT))
    {
      dstGType = WLZ_GREY_SHORT;
    }
    else if((*(gTypeA + 0) == WLZ_GREY_UBYTE) ||
	    (*(gTypeA + 1) == WLZ_GREY_UBYTE) ||
	    (*(gTypeA + 2) == WLZ_GREY_UBYTE))
    {
      dstGType = WLZ_GREY_SHORT;
    }
    else
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
  }
  /* Make destination object with intersection domain and new values. */
  if(errNum == WLZ_ERR_NONE)
  {
    (void )WlzValueConvertPixel(&dstBgd, *(bgdA + 0), dstGType);
    dstVal.v = WlzNewValueTb(*(iObjA + 0),
    			     WlzGreyTableType(WLZ_GREY_TAB_RAGR,
			     dstGType, NULL),
			     dstBgd, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, istObj->domain, dstVal,
      			   NULL, NULL, &errNum);
    }
  }
  if(istObj)
  {
    WlzFreeObj(istObj);
  }
  /* Make buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz = dstObj->domain.i->lastkl - dstObj->domain.i->kol1 + 1;
    if(AlcDouble2Malloc(&iBufA, 3, bufSz) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Scan through the objects computing the magnitude. */
  if(errNum == WLZ_ERR_NONE)
  {
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < 3))
    {
      errNum = WlzInitGreyScan(*(iObjA + idN), iIWSpA + idN, iGWSpA + idN);
      ++idN;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(dstObj, &dstIWSp, &dstGWSp);
    }
    while((errNum == WLZ_ERR_NONE) &&
    	  ((errNum = WlzNextGreyInterval(iIWSpA + 0)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(iIWSpA + 1)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(iIWSpA + 2)) == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&dstIWSp)) == WLZ_ERR_NONE))
    {
      itvLen = dstIWSp.rgtpos - dstIWSp.lftpos + 1;
      /* Copy intervals to double buffers. */
      idN = 0;
      while(idN < 3)
      {
	tGP0.dbp = *(iBufA + idN);
        WlzValueCopyGreyToGrey(tGP0, 0,
			       WLZ_GREY_DOUBLE,
			       (iGWSpA + idN)->u_grintptr, 0,
			       (iGWSpA + idN)->pixeltype,
			       itvLen);
        ++idN;
      }
      /* Compute magnitude. */
      WlzBufMagD3(*(iBufA + 0), *(iBufA + 1), *(iBufA + 2), itvLen);
      /* Clamp into destination interval. */
      tGP0.dbp = *(iBufA + 0);
      WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0, dstGWSp.pixeltype,
				tGP0, 0, WLZ_GREY_DOUBLE,
				itvLen);
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  /* Free intersection objects. */
  idN = 0;
  while(idN < 3)
  {
    if(iObjA[idN])
    {
      WlzFreeObj(iObjA[idN]);
    }
    ++idN;
  }
  /* Free buffers. */
  if(iBufA)
  {
    Alc2Free((void **)iBufA);
  }
  /* Tidy up on error. */
  if(dstObj && (errNum != WLZ_ERR_NONE))
  {
    WlzFreeObj(dstObj);
    dstObj = NULL;
  }
  /* Pass back error status. */
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
* Purpose:	Computes the magnitude of the gray values in the
*		3 given Woolz 3D domain objects.
* Global refs:	-
* Parameters:	WlzObject *srcObj0:	First object.
*		WlzObject *srcObj1:	Second object.
*		WlzObject *srcObj2:	Third object.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static 	WlzObject *WlzGreyMagnitude3D(WlzObject *srcObj0, WlzObject *srcObj1,
				      WlzObject *srcObj2, WlzErrorNum *dstErr)
{
  int		idN,
		idP,
  		nPlanes;
  WlzDomain	dstDom;
  WlzValues	dstVal;
  WlzObject	*dstObj2D,
  		*istObj = NULL,
		*dstObj = NULL;
  WlzObject	*srcObjA[3],
  		*iObj3DA[3],
  		*iObj2DA[3];
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  dstVal.core = NULL;
  iObj3DA[0] = iObj3DA[1] = iObj3DA[2] = NULL;
  if((srcObj0->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN) ||
     (srcObj1->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN) ||
     (srcObj2->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if ((srcObj0->values.core->type != WLZ_VOXELVALUETABLE_GREY) ||
           (srcObj1->values.core->type != WLZ_VOXELVALUETABLE_GREY) ||
	   (srcObj2->values.core->type != WLZ_VOXELVALUETABLE_GREY))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    srcObjA[0] = srcObj0;
    srcObjA[1] = srcObj1;
    srcObjA[2] = srcObj2;
    istObj = WlzIntersectN(3, srcObjA, 0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(istObj->type == WLZ_EMPTY_OBJ)
    {
      dstObj = istObj;
    }
    else
    {
      dstDom = WlzAssignDomain(istObj->domain, NULL);
      WlzFreeObj(istObj);
      idN = 0;
      while((errNum == WLZ_ERR_NONE) && (idN < 3))
      {
        iObj3DA[idN] = WlzMakeMain(WLZ_3D_DOMAINOBJ,
				   dstDom, srcObjA[idN]->values,
				   NULL, NULL, &errNum);
        ++idN;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	idP = 0;
        nPlanes = dstDom.p->lastpl - dstDom.p->plane1 + 1;
	while((errNum == WLZ_ERR_NONE) &&
	      (idP < nPlanes))
        {
	  if((dstDom.p->domains + idN)->core == NULL)
	  {
	    (dstVal.vox->values + idP)->core = NULL;
	  }
	  else
	  {
	    dstObj2D = NULL;
	    iObj2DA[0] = iObj2DA[1] = iObj2DA[2] = NULL;
	    idN = 0;
	    /* Make 2D objects. */
	    while((errNum == WLZ_ERR_NONE) && (idN < 3))
	    {
	      iObj2DA[idN] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				    *(iObj3DA[idN]->domain.p->domains + idP),
				    *(iObj3DA[idN]->values.vox->values + idP),
				    NULL, NULL, &errNum);
	      ++idN;
	    }
	    /* Compute magnitude. */
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstObj2D = WlzGreyMagnitude2D3(iObj2DA[0], iObj2DA[1],
					     iObj2DA[2], &errNum);
	    }
	    /* If first plane get grey type and background, then create new
	     * 3D object with new voxel values. */
	    if((idP == 0) && (errNum == WLZ_ERR_NONE))
	    {
      	      bgdV = WlzGetBackground(iObj2DA[0], &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		dstVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
						 dstDom.p->plane1,
						 dstDom.p->lastpl,
						 bgdV, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstVal,
				     NULL, NULL, &errNum);
	      }
	    }
	    /* Add 2D values to the 3D objects values. */
	    if(errNum == WLZ_ERR_NONE)
	    {
	      *(dstVal.vox->values + idP) = WlzAssignValues(dstObj2D->values,
	      						    NULL);
	    }
	    /* Free temporary 2D object. */
	    if(dstObj2D)
	    {
	      WlzFreeObj(dstObj2D);
	    }
	    idN = 0;
	    while(idN < 3)
	    {
	      if(iObj2DA[idN])
	      {
		WlzFreeObj(iObj2DA[idN]);
	      }
	      ++idN;
	    }
	  }
	  ++idP;
	}
      }
      /* Free 3D intersection gradient objects. */
      idN = 0;
      while(idN < 3)
      {
	if(iObj3DA[idN])
	{
	  WlzFreeObj(iObj3DA[idN]);
	}
	++idN;
      }
      /* Decrement link count to the destination 3D domain. */
      (void )WlzFreeDomain(dstDom);
    }
  }
  if((errNum != WLZ_ERR_NONE) && (dstObj != NULL))
  {
    WlzFreeObj(dstObj);
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

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
    grdX = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_X, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grdY = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_Y, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzImageArithmetic(grdY, grdX, WLZ_BO_MAGNITUDE, 0, &errNum);
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
  WlzObject	*grdX = NULL,
  		*grdY = NULL,
		*grdZ = NULL,
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
    grdX = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_X, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grdY = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_Y, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grdZ = WlzRsvFilterObj(srcObj, flt, WLZ_RSVFILTER_ACTION_Z, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzGreyMagnitude3D(grdZ, grdY, grdX, &errNum);
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
  if(dstGrdZ && (errNum == WLZ_ERR_NONE))
  {
    *dstGrdZ = grdZ;
  }
  else if(grdZ)
  {
    WlzFreeObj(grdZ);
  }
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
