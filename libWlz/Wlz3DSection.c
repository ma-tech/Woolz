#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        Wlz3DSection.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for cutting 2D sections from 3D objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

static WlzGreyV		val;
static WlzPixelP	pix;

static WlzPixelP convert_grey_type(
  WlzPixelP	pixptr,
  WlzGreyType	grey_type)
{
  pix.type = grey_type;
  pix.p.inp = &(val.inv);

  switch( pixptr.type ){
  case WLZ_GREY_INT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      val.inv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_SHORT:
      val.shv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_UBYTE:
      val.ubv = (UBYTE) *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_FLOAT:
      val.flv = *(pixptr.p.inp);
      return pix;
    case WLZ_GREY_DOUBLE:
      val.dbv = *(pixptr.p.inp);
      return pix;
    }
  case WLZ_GREY_SHORT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      val.inv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_SHORT:
      val.shv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_UBYTE:
      val.ubv = (UBYTE) *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_FLOAT:
      val.flv = *(pixptr.p.shp);
      return pix;
    case WLZ_GREY_DOUBLE:
      val.dbv = *(pixptr.p.shp);
      return pix;
    }
  case WLZ_GREY_UBYTE:
    switch( grey_type ){
    case WLZ_GREY_INT:
      val.inv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_SHORT:
      val.shv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_UBYTE:
      val.ubv = (UBYTE) *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_FLOAT:
      val.flv = *(pixptr.p.ubp);
      return pix;
    case WLZ_GREY_DOUBLE:
      val.dbv = *(pixptr.p.ubp);
      return pix;
    }
  case WLZ_GREY_FLOAT:
    switch( grey_type ){
    case WLZ_GREY_INT:
      val.inv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_SHORT:
      val.shv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_UBYTE:
      val.ubv = (UBYTE) *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_FLOAT:
      val.flv = *(pixptr.p.flp);
      return pix;
    case WLZ_GREY_DOUBLE:
      val.dbv = *(pixptr.p.flp);
      return pix;
    }
  case WLZ_GREY_DOUBLE:
    switch( grey_type ){
    case WLZ_GREY_INT:
      val.inv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_SHORT:
      val.shv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_UBYTE:
      val.ubv = (UBYTE) *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_FLOAT:
      val.flv = *(pixptr.p.dbp);
      return pix;
    case WLZ_GREY_DOUBLE:
      val.dbv = *(pixptr.p.dbp);
      return pix;
    }
  }
  return(pix);
}

/************************************************************************
*   Function   : WlzGetSectionFromObject				*
*   Date       : Thu Feb  6 15:20:23 1997				*
*************************************************************************
*   Synopsis   :Return a 2D object of the same grey-type as the given	*
*		object corresponding to the planar section defined by	*
*		the given view structure. Only an area defined by the	*
*		view-structure reference object is filled. Currently	*
*		binary images are handled by filling an array and      	*
*		thresholding.						*
*   Returns    :WlzObject *: the section as a woolz 2D object		*
*   Parameters :WlzObject           *obj: source object			*
*		WlzThreeDViewStruct *viewStr: structure defining the	*
*			section.					*
*   Global refs:None.							*
************************************************************************/
WlzObject *WlzGetSectionFromObject(
  WlzObject		*obj,
  WlzThreeDViewStruct	*viewStr,
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
	yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
	for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	  xp = k - WLZ_NINT(viewStr->minvals.vtX);
	  vtx.vtX = viewStr->xp_to_x[xp] + viewStr->yp_to_x[yp];
	  vtx.vtY = viewStr->xp_to_y[xp] + viewStr->yp_to_y[yp];
	  vtx.vtZ = viewStr->xp_to_z[xp] + viewStr->yp_to_z[yp];
	  vtx.vtX = WLZ_NINT(vtx.vtX);
	  vtx.vtY = WLZ_NINT(vtx.vtY);
	  vtx.vtZ = WLZ_NINT(vtx.vtZ);
	  if( voxvals ){
	    WlzGreyValueGet(gVWSp, vtx.vtZ, vtx.vtY, vtx.vtX);
	    pixptr.p = *(gVWSp->gPtr);
	    pixptr.type = gVWSp->gType;
	    if( pixptr.type != grey_type ){
	      pixptr = convert_grey_type(pixptr, grey_type);
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
	    default:
	      break;
	    }
	  }
	  else {
	    if( WlzInsideDomain(obj, vtx.vtZ, vtx.vtY, vtx.vtX, NULL) ){
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
    tmp_obj = WlzThreshold(newobj, pixval, WLZ_THRESH_HIGH, &errNum);
    WlzFreeObj(newobj);
    newobj = tmp_obj;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return newobj;
}
