#pragma ident "MRC HGU $Id$"
/*!
* \file         Wlz3DSection.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for cutting 2D sections from 3D objects.
* \ingroup	WlzSectionTransform
* \todo         -
* \bug          None known.
*/
#include <limits.h>
#include <float.h>

#include <Wlz.h>

static WlzObject 		*WlzGetSectionFrom3DDomObj(
  				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzInterpolationType	interp,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzGetMaskedSectionFrom3DDomObj(
  				  WlzObject *obj,
				  WlzThreeDViewStruct *viewStr,
				  WlzInterpolationType	interp,
				  WlzErrorNum *dstErr);
static WlzContour 		*WlzGetSectionFromCtr(
				  WlzContour *ctr,
				  WlzThreeDViewStruct *view,
				  WlzInterpolationType	interp,
				  WlzErrorNum *dstErr);
static WlzPixelP 		WlzGetSectionConvertGreyType(
				  WlzPixelP pixptr,
				  WlzGreyType grey_type);


/*!
* \return				A new 2D object cut from the given
*					3D object.
* \ingroup	WlzSectionTransform
* \brief	Cuts the 2D object which lies on the plane specified by
*		the given view structure from the given 3D object.
*		If the given object is a 3D domain object with grey values
*		then a new 2D object is created with the same grey-type
*		as the given object. Only grey values within the area
*		defined by the view structure reference object are
*		extracted. Returns a rectangular object and value table
*		with size determined by the bounding box which encloses
* 		the bounding box of the original but in the viewing direction.
* \param	obj			Given 3D object.
* \param	view			The given view structure.
* \param	interp			Interpolation type - nearest neighbour
 or linear
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
WlzObject 	*WlzGetSectionFromObject(
  WlzObject 		*obj,
  WlzThreeDViewStruct 	*view,
  WlzInterpolationType	interp,
  WlzErrorNum 		*dstErr)
{
  WlzObject	*newObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_3D_DOMAINOBJ:
        newObj = WlzGetSectionFrom3DDomObj(obj, view, interp, &errNum);
        break;
      case WLZ_CONTOUR:
	dom.core = NULL;
	val.core = NULL;
        dom.ctr = WlzGetSectionFromCtr(obj->domain.ctr, view, interp, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  newObj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    WlzFreeContour(dom.ctr);
	  }
	}
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
  return(newObj);
}

/*!
* \return				A new 2D object cut from the given
*					3D object.
* \ingroup	WlzSectionTransform
* \brief	Cuts the 2D object which lies on the plane specified by
*		the given view structure from the given 3D object.
*		If the given object is a 3D domain object with grey values
*		then a new 2D object is created with the same grey-type
*		as the given object. Only grey values within the area
*		defined by the view structure reference object are
*		extracted. Returns an object with domain defined by the
*		section cut through the reference object. The value
*		table is rectangular and the same size as from
*		WlzGetSectionFromObject.
* \param	obj			Given 3D object.
* \param	view			The given view structure.
* \param	interp			Interpolation type - nearest neighbour
 or linear
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
WlzObject 	*WlzGetMaskedSectionFromObject(WlzObject *obj,
					       WlzThreeDViewStruct *view,
					       WlzInterpolationType	interp,
					       WlzErrorNum *dstErr)
{
  WlzObject	*newObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_3D_DOMAINOBJ:
        newObj = WlzGetMaskedSectionFrom3DDomObj(obj, view, interp, &errNum);
        break;
      case WLZ_CONTOUR:
	dom.core = NULL;
	val.core = NULL;
        dom.ctr = WlzGetSectionFromCtr(obj->domain.ctr, view, interp, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  newObj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    WlzFreeContour(dom.ctr);
	  }
	}
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
  return(newObj);
}

/*!
* \return				A new 2D contour cut from the given
*					3D contour.
* \ingroup	WlzSectionTransform
* \brief	Cuts the 2D contour which lies on the plane specified by
*		the given view structure from the given 3D contour.
* \param	ctr			Given contour.
* \param	view			The given view structure.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
static WlzContour *WlzGetSectionFromCtr(WlzContour *ctr,
				        WlzThreeDViewStruct *view,
					WlzInterpolationType	interp,
					WlzErrorNum *dstErr)
{
  WlzGMModel	*newModel = NULL;
  WlzContour	*newCtr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr->model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(ctr->model->type)
    {
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
      case WLZ_GMMOD_3N:
	newModel = WlzGetSectionFromGMModel(ctr->model, view, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  newCtr = WlzMakeContour(&errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    (void )WlzGMModelFree(newModel);
	  }
	  else
	  {
	    newCtr->model = newModel;
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newCtr);
}

/*!
* \return				The section as a new 2D woolz object.
* \ingroup	WlzSectionTransform
* \brief	Cuts a 2D object of the same grey-type as the given
*		object corresponding to the planar section defined by the given
*		view structure. Only an area defined by the view-structure
*		reference object is filled. The  object has a domain which
*		matches the domain of the 3D reference object. This is different
*		to WlzGetSectionFromObj which returns a rectangular object is
*		the reference object has a grey-table.
*		Currently binary images are handled
*		by filling an array and thresholding.
* \param	obj			Given 3D object.
* \param	viewStr			The given view structure.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
static WlzObject *WlzGetMaskedSectionFrom3DDomObj(
  WlzObject		*obj,
  WlzThreeDViewStruct	*viewStr,
  WlzInterpolationType	interp,
  WlzErrorNum	*dstErr)
{
  WlzObject		*newobj, *tmp_obj, *mask;
  WlzDomain		domain;
  WlzValues		values;
  WlzVoxelValues	*voxvals;
  WlzGreyType		grey_type;
  WlzPixelV		pixval;
  WlzPixelP		pixptr;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzIntervalWSpace	iwsp, iwsp1;
  WlzGreyWSpace		gwsp, gwsp1;
  WlzFVertex3		vtx;
  int			k, xp, yp, p;
  WlzDVertex3		tDV0, tDV1;
  double		tD0;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  newobj = NULL;
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) && (obj->type != WLZ_3D_DOMAINOBJ) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  if( (errNum == WLZ_ERR_NONE) && (obj->domain.core == NULL) ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) &&
     (obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN) ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  /* set local pointers and get the grey-type */
  if( errNum == WLZ_ERR_NONE ){
    grey_type = WLZ_GREY_UBYTE;
    voxvals = NULL;
    pixval.type = WLZ_GREY_UBYTE;
    pixval.v.ubv = (UBYTE) 0;
    if( obj->values.core ){
      voxvals = obj->values.vox;
      pixval = voxvals->bckgrnd;
      if( voxvals->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VOXELVALUES_TYPE;
      }
      else {
	for(p=0; p < (voxvals->lastpl - voxvals->plane1 + 1); p++){
	  if( voxvals->values[p].core ){
	    grey_type = 
	      WlzGreyTableTypeToGreyType(voxvals->values[p].core->type,
					 &errNum);
	    break;
	  }
	}
      }
    }
  }
      

  /* check the view structure */
  if( (errNum == WLZ_ERR_NONE) && (viewStr == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) && (viewStr->type != WLZ_3D_VIEW_STRUCT) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /* check the interpolation parameter */
  if( errNum == WLZ_ERR_NONE ){
    switch( interp ){
    case WLZ_INTERPOLATION_NEAREST:
    case WLZ_INTERPOLATION_LINEAR:
      break;

    default:
      errNum = WLZ_ERR_INTERPOLATION_TYPE;
      break;
    }
  }

  /* create a new rectangular object */
  if( errNum == WLZ_ERR_NONE ){
    domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     WLZ_NINT(viewStr->minvals.vtY),
				     WLZ_NINT(viewStr->maxvals.vtY),
				     WLZ_NINT(viewStr->minvals.vtX),
				     WLZ_NINT(viewStr->maxvals.vtX),
				     &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    newobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
			 &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    values.v = WlzNewValueTb(newobj,
			     WlzGreyTableType(WLZ_GREY_TAB_RECT,
					      grey_type, NULL),
			     pixval, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    newobj->values = WlzAssignValues(values, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && (voxvals)){
    pixval.type = WLZ_GREY_UBYTE;
    pixval.v.ubv = (UBYTE) 0;
    values.v = WlzNewValueTb(newobj,
			     WlzGreyTableType(WLZ_GREY_TAB_RECT,
					      grey_type, NULL),
			     pixval, &errNum);
    mask = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
		       &errNum);
  }


  /* scan object setting values */
  if( errNum == WLZ_ERR_NONE ){
    if( voxvals ){
      gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
      errNum = WlzInitGreyScan(mask, &iwsp1, &gwsp1);
    }
    if( errNum == WLZ_ERR_NONE ){
      errNum = WlzInitGreyScan(newobj, &iwsp, &gwsp);
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	if( voxvals ){
	  (void) WlzNextGreyInterval(&iwsp1);
	}
	yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
	for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	  xp = k - WLZ_NINT(viewStr->minvals.vtX);
	  vtx.vtX = viewStr->xp_to_x[xp] + viewStr->yp_to_x[yp];
	  vtx.vtY = viewStr->xp_to_y[xp] + viewStr->yp_to_y[yp];
	  vtx.vtZ = viewStr->xp_to_z[xp] + viewStr->yp_to_z[yp];

	  /* apply interpolation */
	  if( voxvals ){
	    switch( interp ){
	    case WLZ_INTERPOLATION_NEAREST:
	      
	      WlzGreyValueGet(gVWSp, WLZ_NINT(vtx.vtZ),
			      WLZ_NINT(vtx.vtY), WLZ_NINT(vtx.vtX));
	      pixptr.p = *(gVWSp->gPtr);
	      pixptr.type = gVWSp->gType;
	      break;

	    case WLZ_INTERPOLATION_LINEAR:
	      /* set the background to current value */
	      WlzGreyValueGetCon(gVWSp, floor(vtx.vtZ),
				 floor(vtx.vtY), floor(vtx.vtX));
	      pixptr.p = *(gVWSp->gPtr);
	      pixptr.type = gVWSp->gType;

	      tDV0.vtX = vtx.vtX - floor(vtx.vtX);
	      tDV0.vtY = vtx.vtY - floor(vtx.vtY);
	      tDV0.vtZ = vtx.vtZ - floor(vtx.vtZ);
	      tDV1.vtX = 1.0 - tDV0.vtX;
	      tDV1.vtY = 1.0 - tDV0.vtY;
	      tDV1.vtZ = 1.0 - tDV0.vtZ;
	      switch(gVWSp->gType)
	      {
	      case WLZ_GREY_INT:
		tD0 = ((gVWSp->gVal[0]).inv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).inv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).inv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).inv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).inv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).inv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).inv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).inv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, INT_MIN, INT_MAX);
		*(pixptr.p.inp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_SHORT:
		tD0 = ((gVWSp->gVal[0]).shv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).shv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).shv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).shv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).shv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).shv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).shv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).shv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, SHRT_MIN, SHRT_MAX);
		*(pixptr.p.shp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_UBYTE:
		tD0 = ((gVWSp->gVal[0]).ubv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).ubv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).ubv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).ubv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).ubv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).ubv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).ubv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).ubv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, 0, 255);
		*(pixptr.p.ubp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_FLOAT:
		tD0 = ((gVWSp->gVal[0]).flv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).flv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).flv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).flv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).flv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).flv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).flv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).flv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, FLT_MIN, FLT_MAX);
		*(pixptr.p.flp) = tD0;
		break;
	      case WLZ_GREY_DOUBLE:
		tD0 = ((gVWSp->gVal[0]).dbv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).dbv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).dbv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).dbv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).dbv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).dbv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).dbv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).dbv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		*(pixptr.p.dbp) = tD0;
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	      }
	      break;
	    }

	    /* now copy the value */
	    if( pixptr.type != grey_type ){
	      pixptr = WlzGetSectionConvertGreyType(pixptr, grey_type);
	    }
	    switch( grey_type ){
	    case WLZ_GREY_INT:
	      *(gwsp.u_grintptr.inp) = *(pixptr.p.inp);
	      gwsp.u_grintptr.inp++;
	      break;
	    case WLZ_GREY_SHORT:
	      *(gwsp.u_grintptr.shp) = *(pixptr.p.shp);
	      gwsp.u_grintptr.shp++;
	      break;
	    case WLZ_GREY_UBYTE:
	      *(gwsp.u_grintptr.ubp) = *(pixptr.p.ubp);
	      gwsp.u_grintptr.ubp++;
	      break;
	    case WLZ_GREY_FLOAT:
	      *(gwsp.u_grintptr.flp) = *(pixptr.p.flp);
	      gwsp.u_grintptr.flp++;
	      break;
	    case WLZ_GREY_DOUBLE:
	      *(gwsp.u_grintptr.dbp) = *(pixptr.p.dbp);
	      gwsp.u_grintptr.dbp++;
	      break;
	    case WLZ_GREY_RGBA:
	      *(gwsp.u_grintptr.rgbp) = *(pixptr.p.rgbp);
	      gwsp.u_grintptr.rgbp++;
	      break;
	    default:
	      break;
	    }

	    /* set the mask value */
	    if( WlzInsideDomain(obj, WLZ_NINT(vtx.vtZ), WLZ_NINT(vtx.vtY),
				WLZ_NINT(vtx.vtX), NULL) ){
	      *(gwsp1.u_grintptr.ubp) = 0;
	    }
	    else {
	      *(gwsp1.u_grintptr.ubp) = 128;
	    }
	    gwsp1.u_grintptr.ubp++;
	  }
	  else {
	    if( WlzInsideDomain(obj, WLZ_NINT(vtx.vtZ), WLZ_NINT(vtx.vtY),
				WLZ_NINT(vtx.vtX), NULL) ){
	      *(gwsp.u_grintptr.ubp) = 128;
	    }
	    else {
	      *(gwsp.u_grintptr.ubp) = 0;
	    }
	    gwsp.u_grintptr.ubp++;
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
      if( gVWSp ){
	WlzGreyValueFreeWSp(gVWSp);
      }
    }
  }

  /* if binary then threshold and free valuetable */
  if( (errNum == WLZ_ERR_NONE) ){
    pixval.type = WLZ_GREY_INT;
    pixval.v.inv = 128;
    if( voxvals ){
      mask = WlzAssignObject(mask, NULL);
      newobj = WlzAssignObject(newobj, NULL);
      tmp_obj = WlzThreshold(mask, pixval, WLZ_THRESH_HIGH, &errNum);
      WlzFreeObj(mask);
      mask = WlzMakeMain(newobj->type, tmp_obj->domain, newobj->values,
			 NULL, NULL, &errNum);
      WlzFreeObj(newobj);
      WlzFreeObj(tmp_obj);
      newobj = mask;
    }
    else {
      newobj = WlzAssignObject(newobj, NULL);
      tmp_obj = WlzThreshold(newobj, pixval, WLZ_THRESH_HIGH, &errNum);
      WlzFreeObj(newobj);
      newobj = tmp_obj;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return newobj;
}

/*!
* \return				The section as a new 2D woolz object.
* \ingroup	WlzSectionTransform
* \brief	Cuts a 2D object of the same grey-type as the given
*		object corresponding to the planar section defined by the given
*		view structure. Only an area defined by the view-structure
*		reference object is filled. Currently binary images are handled
*		by filling an array and thresholding.
* \param	obj			Given 3D object.
* \param	viewStr			The given view structure.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
static WlzObject *WlzGetSectionFrom3DDomObj(
  WlzObject		*obj,
  WlzThreeDViewStruct	*viewStr,
  WlzInterpolationType	interp,
  WlzErrorNum	*dstErr)
{
  WlzObject		*newobj, *tmp_obj;
  WlzDomain		domain;
  WlzValues		values;
  WlzVoxelValues	*voxvals;
  WlzGreyType		grey_type;
  WlzPixelV		pixval;
  WlzPixelP		pixptr;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzFVertex3		vtx;
  int			k, xp, yp, p;
  WlzDVertex3		tDV0, tDV1;
  double		tD0;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  newobj = NULL;
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) && (obj->type != WLZ_3D_DOMAINOBJ) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  if( (errNum == WLZ_ERR_NONE) && (obj->domain.core == NULL) ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) &&
     (obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN) ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  /* set local pointers and get the grey-type */
  if( errNum == WLZ_ERR_NONE ){
    grey_type = WLZ_GREY_UBYTE;
    voxvals = NULL;
    pixval.type = WLZ_GREY_UBYTE;
    pixval.v.ubv = (UBYTE) 0;
    if( obj->values.core ){
      voxvals = obj->values.vox;
      pixval = voxvals->bckgrnd;
      if( voxvals->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VOXELVALUES_TYPE;
      }
      else {
	for(p=0; p < (voxvals->lastpl - voxvals->plane1 + 1); p++){
	  if( voxvals->values[p].core ){
	    grey_type = 
	      WlzGreyTableTypeToGreyType(voxvals->values[p].core->type,
					 &errNum);
	    break;
	  }
	}
      }
    }
  }
      

  /* check the view structure */
  if( (errNum == WLZ_ERR_NONE) && (viewStr == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (errNum == WLZ_ERR_NONE) && (viewStr->type != WLZ_3D_VIEW_STRUCT) ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /* check the interpolation parameter */
  if( errNum == WLZ_ERR_NONE ){
    switch( interp ){
    case WLZ_INTERPOLATION_NEAREST:
    case WLZ_INTERPOLATION_LINEAR:
      break;

    default:
      errNum = WLZ_ERR_INTERPOLATION_TYPE;
      break;
    }
  }

  /* create a new rectangular object */
  if( errNum == WLZ_ERR_NONE ){
    domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     WLZ_NINT(viewStr->minvals.vtY),
				     WLZ_NINT(viewStr->maxvals.vtY),
				     WLZ_NINT(viewStr->minvals.vtX),
				     WLZ_NINT(viewStr->maxvals.vtX),
				     &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    newobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
			 &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    values.v = WlzNewValueTb(newobj,
			     WlzGreyTableType(WLZ_GREY_TAB_RECT,
					      grey_type, NULL),
			     pixval, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    newobj->values = WlzAssignValues(values, &errNum);
  }

  /* scan object setting values */
  if( errNum == WLZ_ERR_NONE ){
    if( voxvals ){
      gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
    }
    if( errNum == WLZ_ERR_NONE ){
      errNum = WlzInitGreyScan(newobj, &iwsp, &gwsp);
      while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
	yp = iwsp.linpos - newobj->domain.i->line1;
	for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	  xp = k - newobj->domain.i->kol1;
	  vtx.vtX = viewStr->xp_to_x[xp] + viewStr->yp_to_x[yp];
	  vtx.vtY = viewStr->xp_to_y[xp] + viewStr->yp_to_y[yp];
	  vtx.vtZ = viewStr->xp_to_z[xp] + viewStr->yp_to_z[yp];

	  /* apply interpolation */
	  if( voxvals ){
	    switch( interp ){
	    case WLZ_INTERPOLATION_NEAREST:
	      
	      WlzGreyValueGet(gVWSp, WLZ_NINT(vtx.vtZ),
			      WLZ_NINT(vtx.vtY), WLZ_NINT(vtx.vtX));
	      pixptr.p = *(gVWSp->gPtr);
	      pixptr.type = gVWSp->gType;
	      break;

	    case WLZ_INTERPOLATION_LINEAR:
	      WlzGreyValueGetCon(gVWSp, floor(vtx.vtZ),
				 floor(vtx.vtY), floor(vtx.vtX));
	      pixptr.p = *(gVWSp->gPtr);
	      pixptr.type = gVWSp->gType;

	      tDV0.vtX = vtx.vtX - floor(vtx.vtX);
	      tDV0.vtY = vtx.vtY - floor(vtx.vtY);
	      tDV0.vtZ = vtx.vtZ - floor(vtx.vtZ);
	      tDV1.vtX = 1.0 - tDV0.vtX;
	      tDV1.vtY = 1.0 - tDV0.vtY;
	      tDV1.vtZ = 1.0 - tDV0.vtZ;
	      switch(gVWSp->gType)
	      {
	      case WLZ_GREY_INT:
		tD0 = ((gVWSp->gVal[0]).inv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).inv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).inv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).inv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).inv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).inv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).inv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).inv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, INT_MIN, INT_MAX);
		*(pixptr.p.inp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_SHORT:
		tD0 = ((gVWSp->gVal[0]).shv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).shv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).shv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).shv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).shv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).shv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).shv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).shv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, SHRT_MIN, SHRT_MAX);
		*(pixptr.p.shp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_UBYTE:
		tD0 = ((gVWSp->gVal[0]).ubv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).ubv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).ubv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).ubv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).ubv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).ubv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).ubv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).ubv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, 0, 255);
		*(pixptr.p.ubp) = WLZ_NINT(tD0);
		break;
	      case WLZ_GREY_FLOAT:
		tD0 = ((gVWSp->gVal[0]).flv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).flv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).flv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).flv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).flv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).flv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).flv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).flv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		tD0 = WLZ_CLAMP(tD0, FLT_MIN, FLT_MAX);
		*(pixptr.p.flp) = tD0;
		break;
	      case WLZ_GREY_DOUBLE:
		tD0 = ((gVWSp->gVal[0]).dbv *
		       tDV1.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[1]).dbv *
		   tDV0.vtX * tDV1.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[2]).dbv *
		   tDV1.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[3]).dbv *
		   tDV0.vtX * tDV0.vtY * tDV1.vtZ) +
		  ((gVWSp->gVal[4]).dbv *
		   tDV1.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[5]).dbv *
		   tDV0.vtX * tDV1.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[6]).dbv *
		   tDV1.vtX * tDV0.vtY * tDV0.vtZ) +
		  ((gVWSp->gVal[7]).dbv *
		   tDV0.vtX * tDV0.vtY * tDV0.vtZ);
		*(pixptr.p.dbp) = tD0;
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	      }
	      break;
	    }

	    /* now copy the value */
	    if( pixptr.type != grey_type ){
	      pixptr = WlzGetSectionConvertGreyType(pixptr, grey_type);
	    }
	    switch( grey_type ){
	    case WLZ_GREY_INT:
	      *(gwsp.u_grintptr.inp) = *(pixptr.p.inp);
	      gwsp.u_grintptr.inp++;
	      break;
	    case WLZ_GREY_SHORT:
	      *(gwsp.u_grintptr.shp) = *(pixptr.p.shp);
	      gwsp.u_grintptr.shp++;
	      break;
	    case WLZ_GREY_UBYTE:
	      *(gwsp.u_grintptr.ubp) = *(pixptr.p.ubp);
	      gwsp.u_grintptr.ubp++;
	      break;
	    case WLZ_GREY_FLOAT:
	      *(gwsp.u_grintptr.flp) = *(pixptr.p.flp);
	      gwsp.u_grintptr.flp++;
	      break;
	    case WLZ_GREY_DOUBLE:
	      *(gwsp.u_grintptr.dbp) = *(pixptr.p.dbp);
	      gwsp.u_grintptr.dbp++;
	      break;
	    case WLZ_GREY_RGBA:
	      *(gwsp.u_grintptr.rgbp) = *(pixptr.p.rgbp);
	      gwsp.u_grintptr.rgbp++;
	      break;
	    default:
	      break;
	    }
	  }
	  else {
	    if( WlzInsideDomain(obj, WLZ_NINT(vtx.vtZ), WLZ_NINT(vtx.vtY),
				WLZ_NINT(vtx.vtX), NULL) ){
	      *(gwsp.u_grintptr.ubp) = 128;
	    }
	    else {
	      *(gwsp.u_grintptr.ubp) = 0;
	    }
	    gwsp.u_grintptr.ubp++;
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
      if( gVWSp ){
	WlzGreyValueFreeWSp(gVWSp);
      }
    }
  }

  /* if binary then threshold and free valuetable */
  if( (errNum == WLZ_ERR_NONE) && (voxvals == NULL) ){
    pixval.type = WLZ_GREY_INT;
    pixval.v.inv = 128;
    newobj = WlzAssignObject(newobj, NULL);
    tmp_obj = WlzThreshold(newobj, pixval, WLZ_THRESH_HIGH, &errNum);
    WlzFreeObj(newobj);
    newobj = tmp_obj;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return newobj;
}

/*!
* \return				Converted pixel.
* \brief	Convert the type of the given pixel.
* \param	pixptr			Given pixel.
* \param	grey_type		Required grey type.
*/
static WlzGreyV GetSectionStaticGreyVal;
static WlzPixelP WlzGetSectionConvertGreyType(
  WlzPixelP	pixptr,
  WlzGreyType	grey_type)
{
  /*WlzGreyV	val;*/
  WlzPixelP	pix;
  UINT		uval;

  pix.type = grey_type;
  pix.p.inp = &(GetSectionStaticGreyVal.inv);

  switch( pixptr.type ){
  case WLZ_GREY_INT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (UBYTE) *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_RGBA:
      uval = WLZ_CLAMP(*(pixptr.p.inp), 0, 255);
      GetSectionStaticGreyVal.rgbv = uval + (uval<<8) + (uval<<16) + 0xff000000;
      return pix;
    }
  case WLZ_GREY_SHORT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (UBYTE) *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_RGBA:
      uval = WLZ_CLAMP(*(pixptr.p.shp), 0, 255);
      GetSectionStaticGreyVal.rgbv = uval + (uval<<8) + (uval<<16) + 0xff000000;
      return pix;
    }
  case WLZ_GREY_UBYTE:
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (UBYTE) *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_RGBA:
      uval = *(pixptr.p.ubp);
      GetSectionStaticGreyVal.rgbv = uval + (uval<<8) + (uval<<16) + 0xff000000;
      return pix;
    }
  case WLZ_GREY_FLOAT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (UBYTE) *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_RGBA:
      uval = WLZ_CLAMP(*(pixptr.p.flp), 0, 255);
      GetSectionStaticGreyVal.rgbv = uval + (uval<<8) + (uval<<16) + 0xff000000;
      return pix;
    }
  case WLZ_GREY_DOUBLE:
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (UBYTE) *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_RGBA:
      uval = WLZ_CLAMP(*(pixptr.p.dbp), 0, 255);
      GetSectionStaticGreyVal.rgbv = uval + (uval<<8) + (uval<<16) + 0xff000000;
      return pix;
    }
  case WLZ_GREY_RGBA:
    uval = WLZ_RGBA_MODULUS(*(pixptr.p.rgbp));
    switch( grey_type ){
    case WLZ_GREY_INT:
      GetSectionStaticGreyVal.inv = uval;
      return pix;
    case WLZ_GREY_SHORT:
      GetSectionStaticGreyVal.shv = uval;
      return pix;
    case WLZ_GREY_UBYTE:
      GetSectionStaticGreyVal.ubv = (uval/sqrt(3.0));
      return pix;
    case WLZ_GREY_FLOAT:
      GetSectionStaticGreyVal.flv = uval;
      return pix;
    case WLZ_GREY_DOUBLE:
      GetSectionStaticGreyVal.dbv = uval;
      return pix;
    case WLZ_GREY_RGBA:
      GetSectionStaticGreyVal.rgbv = *(pixptr.p.rgbp);
      return pix;
    }
  }
  return(pix);
}
