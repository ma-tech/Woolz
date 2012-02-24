#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DSubSection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/Wlz3DSubSection.c
* \author       Richard Baldock
* \date         May 2008
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
* \brief        Return a sub-region of a 3D section via a 3D section transform.
* \ingroup      WlzSectionTransform
*/

#include <limits.h>
#include <float.h>

#include <Wlz.h>

static WlzObject *WlzGetSubSectionFrom3DDomObj(
  WlzObject 		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct 	*viewStr,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum 		*dstErr);


WlzObject 	*WlzGetSubSectionFromObject(
  WlzObject 		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct 	*view,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum 		*dstErr)
{
  WlzObject	*newObj = NULL;
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
	newObj = WlzGetSubSectionFrom3DDomObj(obj, subDomain, view,
					      interp, maskRtn, &errNum);
        break;

      default:
	newObj = WlzGetSectionFromObject(obj, view, interp, &errNum);
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newObj);
}


static WlzObject *WlzGetSubSectionFrom3DDomObj(
  WlzObject		*obj,
  WlzObject		*subDomain,
  WlzThreeDViewStruct	*viewStr,
  WlzInterpolationType	interp,
  WlzObject		**maskRtn,
  WlzErrorNum		*dstErr)
{
  WlzObject		*newObj, *tmpObj, *mask;
  WlzDomain		domain;
  WlzValues		values;
  WlzPixelV		pixval;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzFVertex3		vtx;
  int			k, xp, yp;
  WlzDVertex3		tDV0, tDV1;
  double		tD0;
  int			maskFlg, greyFlg;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check the object */
  newObj = NULL;
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

  /* sort out return object requirements */
  maskFlg = 0;
  greyFlg = 0;
  if( errNum == WLZ_ERR_NONE ){
    if((obj->values.core == NULL)){
      /* object without values therefore mask object only */
      maskFlg = 1;
      greyFlg = 0;
    }
    else {
      greyFlg = 1;
      if( maskRtn ){
	maskFlg = 1;
      }
    }
  }

  /* check the view structure */
  if( (errNum == WLZ_ERR_NONE) && (viewStr == NULL) ){
    if( viewStr == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if( viewStr->type != WLZ_3D_VIEW_STRUCT ){
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if( !(viewStr->initialised) ){
      errNum = WLZ_ERR_OBJECT_DATA;
    }
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

  /* create a new return object - domain only */
  if( errNum == WLZ_ERR_NONE ){
    /* create maximum sized rectangular object */
    if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     WLZ_NINT(viewStr->minvals.vtY),
				     WLZ_NINT(viewStr->maxvals.vtY),
				     WLZ_NINT(viewStr->minvals.vtX),
				     WLZ_NINT(viewStr->maxvals.vtX),
					 &errNum))){
      values.core = NULL;
      tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain,
			   values, NULL, NULL, &errNum);
    }
  }

  /* check for given subdomain - note require intersection so view structure
     LUTs are valid */
  if( errNum == WLZ_ERR_NONE ){
    if( subDomain ){
      switch( subDomain->type ){
      case WLZ_2D_DOMAINOBJ:
	newObj = WlzIntersect2(tmpObj, subDomain, &errNum);
	WlzFreeObj(tmpObj);
	break;

      case WLZ_2D_POLYGON:
      case WLZ_BOUNDLIST:
      case WLZ_RECTANGLE:
	errNum = WLZ_ERR_UNIMPLEMENTED;
	break;

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
    else {
      newObj = tmpObj;
    }
  }

  /* now add a value table */
  if((errNum == WLZ_ERR_NONE) && greyFlg ){
    pixval = WlzGetBackground(obj, &errNum);
    if((values.v = WlzNewValueTb
	(newObj, WlzGreyTableType(WLZ_GREY_TAB_RECT,
				  WlzGreyTypeFromObj(obj, NULL), NULL),
	 pixval, &errNum))){
      newObj->values = WlzAssignValues(values, &errNum);
    }
  }

  /* check if mask required */
  if( errNum == WLZ_ERR_NONE ){
    if( maskFlg ){
      pixval.type = WLZ_GREY_UBYTE;
      pixval.v.ubv = (WlzUByte) 0;
      if((values.v = WlzNewValueTb(newObj,
				   WlzGreyTableType(WLZ_GREY_TAB_RECT,
						    WLZ_GREY_UBYTE, NULL),
				   pixval, &errNum))){
        mask = WlzMakeMain(WLZ_2D_DOMAINOBJ, newObj->domain, values, NULL, NULL,
                              &errNum);  //Changed 10/02/10 by Zsolt Husz
      }

      /* newObj is redundant if mask only required */
      if( !greyFlg ){
	WlzFreeObj(newObj);
	newObj = NULL;
      }
    }
  }

  /* scan object setting values */
  if((errNum == WLZ_ERR_NONE) && greyFlg ){
    if((gVWSp = WlzGreyValueMakeWSp(obj, &errNum))){
      errNum = WlzInitGreyScan(newObj, &iwsp, &gwsp);
      while((errNum == WLZ_ERR_NONE) && 
	    ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE)){
	yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
	for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	  xp = k - WLZ_NINT(viewStr->minvals.vtX);
	  vtx.vtX = (float )(viewStr->xp_to_x[xp] + viewStr->yp_to_x[yp]);
	  vtx.vtY = (float )(viewStr->xp_to_y[xp] + viewStr->yp_to_y[yp]);
	  vtx.vtZ = (float )(viewStr->xp_to_z[xp] + viewStr->yp_to_z[yp]);

	  /* apply interpolation */
	  switch( interp ){
	  case WLZ_INTERPOLATION_NEAREST:
	      
	    WlzGreyValueGet(gVWSp, WLZ_NINT(vtx.vtZ),
			    WLZ_NINT(vtx.vtY), WLZ_NINT(vtx.vtX));
	    WlzValueSetGrey(gwsp.u_grintptr, 0, gVWSp->gVal[0],
			    gVWSp->gType, 1);
	    break;

	  case WLZ_INTERPOLATION_LINEAR:
	    /* set the background to current value */
	    WlzGreyValueGetCon(gVWSp, floor(vtx.vtZ),
			       floor(vtx.vtY), floor(vtx.vtX));

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
	      *(gwsp.u_grintptr.inp) = WLZ_NINT(tD0);
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
	      *(gwsp.u_grintptr.shp) = (short )WLZ_NINT(tD0);
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
	      *(gwsp.u_grintptr.ubp) = (WlzUByte )WLZ_NINT(tD0);
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
	      *(gwsp.u_grintptr.flp) = (float )tD0;
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
	      *(gwsp.u_grintptr.dbp) = tD0;
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    break;
	  }

	  /* increment destination pointer */
	  switch( gVWSp->gType ){
	  case WLZ_GREY_INT:
	    gwsp.u_grintptr.inp++;
	    break;
	  case WLZ_GREY_SHORT:
	    gwsp.u_grintptr.shp++;
	    break;
	  case WLZ_GREY_UBYTE:
	    gwsp.u_grintptr.ubp++;
	    break;
	  case WLZ_GREY_FLOAT:
	    gwsp.u_grintptr.flp++;
	    break;
	  case WLZ_GREY_DOUBLE:
	    gwsp.u_grintptr.dbp++;
	    break;
	  case WLZ_GREY_RGBA:
	    gwsp.u_grintptr.rgbp++;
	    break;
	  default:
	    break;
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
      WlzGreyValueFreeWSp(gVWSp);
    }
  }

  /* check if mask required */
  if((errNum == WLZ_ERR_NONE) && maskFlg ){
    errNum = WlzInitGreyScan(mask, &iwsp, &gwsp);
    while((errNum == WLZ_ERR_NONE) && 
	  ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE) ){
	yp = iwsp.linpos - WLZ_NINT(viewStr->minvals.vtY);
	for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	  xp = k - WLZ_NINT(viewStr->minvals.vtX);
	  vtx.vtX = (float )(viewStr->xp_to_x[xp] + viewStr->yp_to_x[yp]);
	  vtx.vtY = (float )(viewStr->xp_to_y[xp] + viewStr->yp_to_y[yp]);
	  vtx.vtZ = (float )(viewStr->xp_to_z[xp] + viewStr->yp_to_z[yp]);

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
    if(errNum == WLZ_ERR_EOO)	   /* Reset EOO error from end of intervals */ 
    {
      errNum = WLZ_ERR_NONE;
    }

    /* threshold to determine the mask */
    if( errNum == WLZ_ERR_NONE ){
      pixval.type = WLZ_GREY_INT;
      pixval.v.inv = 128;
      mask = WlzAssignObject(mask, NULL);
      tmpObj = WlzThreshold(mask, pixval, WLZ_THRESH_HIGH, &errNum);
      WlzFreeObj(mask);
      values.core = NULL;
      mask = WlzMakeMain(tmpObj->type, tmpObj->domain, values,
			 NULL, NULL, &errNum);
      WlzFreeObj(tmpObj);

      /* set returns */
      if( !greyFlg ){
	/* new object if mask is returned twice */
	if( maskRtn ){
	  *maskRtn = WlzMakeMain(mask->type, mask->domain, mask->values,
			 NULL, NULL, &errNum);
	}
	newObj = mask;
      }
      else if ( maskRtn ){
	*maskRtn = mask;
      }
      else {
	/* should not get here */
	WlzFreeObj(mask);
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return newObj;
}
