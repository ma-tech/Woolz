#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRGBAThreshold_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzRGBAThreshold.c
* \author       Richard Baldock
* \date         July 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Threshold functions for RGBA objects.
* \ingroup	WlzThreshold
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzRGBAMultiThreshold    */
/*! 
* \ingroup      WlzThreshold
* \brief        Apply independent thresholds to each colour channel
 independently and combine according to the settings encoded in
 combineMode. Each channel can have one of two modes: WLZ_BO_AND
 and WLZ_BO_OR. These are encoded into a single mode variable using
 the RGBA macro, e.g.:
 WLZ_RGBA_RGBA_SET(combineMode, redMode, greenMode, blueMode, 255);

 The macro WLZ_RGBA_RED_GET(combineMode) will return redMode and
 similarly for green and blue.
*
* \return       Thresholded object
* \param    obj	Object to be thresholded
* \param    lowVal	RGB low values
* \param    highVal	RGB high values
* \param    combineMode	Combination rules as an RGBA encoded unsigned integer
* \param    dstErr	Error return
* \par      Source:
*                WlzRGBAThreshold.c
*/
WlzObject *WlzRGBAMultiThreshold(
  WlzObject	*obj,
  WlzPixelV	lowVal,
  WlzPixelV	highVal,
  WlzUInt	combineMode,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzObject	*obj1, *obj2;
  WlzValues	values;
  WlzCompoundArray	*cobj=NULL;
  int		low[3], high[3];
  WlzUInt	mode[3]; 

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    /* must be grey-type RGBA or compound with at least three channels */
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else {
	/* create compound object */
	if((cobj = WlzRGBAToCompound(obj, WLZ_RGBA_SPACE_RGB,
	                             &errNum)) != NULL){
	  cobj = (WlzCompoundArray *) WlzAssignObject((WlzObject *) cobj,
						      &errNum);
	}
      }
      break;

    case WLZ_TRANS_OBJ:
      if((obj1 = WlzRGBAMultiThreshold(obj->values.obj, lowVal, highVal,
				       combineMode, &errNum)) != NULL){
	values.obj = WlzAssignObject(obj1, NULL);
	rtnObj = WlzMakeMain(obj->type, obj->domain, values,
			     NULL, obj, &errNum);
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      cobj = (WlzCompoundArray *) WlzAssignObject(obj, &errNum);
      if( cobj->n < 3 ){
	errNum = WLZ_ERR_OBJECT_DATA;
	WlzFreeObj((WlzObject *) cobj);
      }
      else if((cobj->o[0]->values.core == NULL) ||
	      (cobj->o[1]->values.core == NULL) ||
	      (cobj->o[2]->values.core == NULL)){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      break;

    case WLZ_EMPTY_OBJ:
      rtnObj = WlzMakeEmpty(&errNum);
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if((errNum == WLZ_ERR_NONE) && (rtnObj == NULL)){
    if((lowVal.type != WLZ_GREY_RGBA) ||
       (highVal.type != WLZ_GREY_RGBA)){
      errNum = WLZ_ERR_PARAM_TYPE;
    }
    else {
      low[0] = WLZ_RGBA_RED_GET(lowVal.v.rgbv);
      low[1] = WLZ_RGBA_GREEN_GET(lowVal.v.rgbv);
      low[2] = WLZ_RGBA_BLUE_GET(lowVal.v.rgbv);
      high[0] = WLZ_RGBA_RED_GET(highVal.v.rgbv);
      high[1] = WLZ_RGBA_GREEN_GET(highVal.v.rgbv);
      high[2] = WLZ_RGBA_BLUE_GET(highVal.v.rgbv);
    }
  }

  if((errNum == WLZ_ERR_NONE) && (rtnObj == NULL)){
    mode[0] = WLZ_RGBA_RED_GET(combineMode);
    mode[1] = WLZ_RGBA_GREEN_GET(combineMode);
    mode[2] = WLZ_RGBA_BLUE_GET(combineMode);
  }

  /* get thresholded channels */
  if((errNum == WLZ_ERR_NONE) &&
     (cobj != NULL) && (rtnObj == NULL)){
    WlzObject	*objs[3];
    int		i;
    WlzPixelV	threshV;

    for(i=0; i < 3; i++){
      threshV.type = WLZ_GREY_INT;
      threshV.v.inv = low[i];
      if((obj1 = WlzThreshold(cobj->o[i], threshV, WLZ_THRESH_HIGH,
      			      &errNum)) != NULL){
	obj1 = WlzAssignObject(obj1, &errNum);
	threshV.v.inv = high[i] + 1;
	if((obj2 = WlzThreshold(obj1, threshV, WLZ_THRESH_LOW,
	                        &errNum)) != NULL){
	  objs[i] = WlzAssignObject(obj2, &errNum);
	}
	else {
	  objs[i] = NULL;
	}
	WlzFreeObj(obj1);
      }
      else {
	objs[i] = NULL;
      }
    }

    /* combine according to mode
       what to do here? AND against a channel implies that the threshold
       constraint must be satisfied, OR implies it may be satisfied so
       all AND implies intersection, all OR implies union. Otherwise union
       of ORs and intersect with the ANDs. What about XOR? */
    /* find union of all then intersect with ANDs */
    obj1 = WlzAssignObject(WlzMakeEmpty(&errNum), NULL);
    for(i=0; (i < 3) && (errNum == WLZ_ERR_NONE); i++){
      if( objs[i] ){
	obj2 = WlzUnion2(obj1, objs[i], &errNum);
	WlzFreeObj(obj1);
	obj1 = WlzAssignObject(obj2, &errNum);
      }
    }
    for(i=0; i < 3; i++){
      if( objs[i] ){
	if( (mode[i] == WLZ_BO_AND) && (errNum == WLZ_ERR_NONE) ){
	  obj2 = WlzIntersect2(obj1, objs[i], &errNum);
	  WlzFreeObj(obj1);
	  obj1 = WlzAssignObject(obj2, &errNum);
	}
	WlzFreeObj(objs[i]);
      }
    }

    /* create the return object and add grey-table if possible */
    if( obj1 ){
      if( WlzIsEmpty(obj1, &errNum) ){
	rtnObj = WlzMakeMain(obj1->type, obj1->domain, obj1->values,
			     NULL, NULL, &errNum);
      }
      else {
	if( obj1->type == obj->type ){
	  rtnObj = WlzMakeMain(obj1->type, obj1->domain, obj->values,
			       NULL, NULL, &errNum);
	}
	else {
	  rtnObj = WlzMakeMain(obj1->type, obj1->domain, obj1->values,
			       NULL, NULL, &errNum);
	}
      }
      WlzFreeObj(obj1);
    }
  }

  if( cobj ){
    WlzFreeObj((WlzObject *) cobj);
  }

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

typedef struct _WlzVectorThresholdStruct
{
  WlzRGBAThresholdType	type;
  double	       	origin[3];
  double		direction[3];
  double		dist;
  double		eccen;
} WlzVectorThresholdStruct;

int WlzVectorThreshCb(
  WlzObject	*obj,
  void		*clientData,
  WlzThreshCbStr	*cbs)
{
  WlzVectorThresholdStruct *vts=(WlzVectorThresholdStruct *) clientData;
  int		i, rtnVal=-1;
  double	d, d1, vect[3];
  WlzUInt	val;

  switch( cbs->pix.type ){
  case WLZ_GREY_INT:
    vect[0] = vect[1] = vect[2] = *(cbs->pix.p.inp);
    break;
  case WLZ_GREY_SHORT:
    vect[0] = vect[1] = vect[2] = *(cbs->pix.p.shp);
    break;
  case WLZ_GREY_UBYTE:
    vect[0] = vect[1] = vect[2] = *(cbs->pix.p.ubp);
    break;
  case WLZ_GREY_FLOAT:
    vect[0] = vect[1] = vect[2] = *(cbs->pix.p.flp);
    break;
  case WLZ_GREY_DOUBLE:
    vect[0] = vect[1] = vect[2] = *(cbs->pix.p.dbp);
    break;
  case WLZ_GREY_RGBA:
    val = *(cbs->pix.p.rgbp);
    vect[0] = WLZ_RGBA_RED_GET(val);
    vect[1] = WLZ_RGBA_GREEN_GET(val);
    vect[2] = WLZ_RGBA_BLUE_GET(val);
    break;
  default:
    break;
  }

  /* calculate threshold rule */
  switch( vts->type ){
  case WLZ_RGBA_THRESH_SPHERE:
    for(i=0, d=0.0; i < 3; i++){
      d1 = vect[i] - vts->origin[i];
      d += d1 * d1;
    }
    break;

  case WLZ_RGBA_THRESH_SLICE:
    for(i=0, d=0.0; i < 3; i++){
      d += (vect[i] - vts->origin[i]) * vts->direction[i];
    }
    break;

  default:
    rtnVal = -1;
    break;
  }
  if((d >= 0.0) && (d <= vts->dist)){
    rtnVal = 1;
  }
  else {
    rtnVal = 0;
  }

  return rtnVal;
}

WlzObject *WlzRGBASliceThreshold(
  WlzObject	*obj,
  WlzPixelV	lowVal,
  WlzPixelV	highVal,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzVectorThresholdStruct clientData;
  double	x, y, z, s;

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  
  if( errNum == WLZ_ERR_NONE ){

    /* set up the callback data structure */
    clientData.type = WLZ_RGBA_THRESH_SLICE;
    clientData.origin[0] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_RED,
					    &errNum);
    clientData.origin[1] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_GREEN,
					    &errNum);
    clientData.origin[2] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_BLUE,
					    &errNum);
    x = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_RED, &errNum) -
      clientData.origin[0];
    y = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_GREEN, &errNum) -
      clientData.origin[1];
    z = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_BLUE, &errNum) -
      clientData.origin[2];
    s = sqrt(x*x + y*y + z*z + 0.0000001);
    if( s < 0.001 ){
      x = 1.0; y = z = 0.0;
    }
    else {
      x /= s; y /= s; z /= s;
    }
    clientData.direction[0] = x;
    clientData.direction[1] = y;
    clientData.direction[2] = z;
    clientData.dist = s;
    
    /* now call the threshold function */
    rtnObj = WlzCbThreshold(obj, WlzVectorThreshCb, &clientData, &errNum);
  }

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzObject *WlzRGBABoxThreshold(
  WlzObject	*obj,
  WlzPixelV	lowVal,
  WlzPixelV	highVal,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  int		h, l;
  WlzUInt	combineMode;

  /* check and reset low annd high values */
  l = WLZ_RGBA_RED_GET(lowVal.v.rgbv);
  h = WLZ_RGBA_RED_GET(highVal.v.rgbv);
  if( l > h ){
    WLZ_RGBA_RED_SET(highVal.v.rgbv, l);
    WLZ_RGBA_RED_SET(lowVal.v.rgbv, h);
  }
  l = WLZ_RGBA_GREEN_GET(lowVal.v.rgbv);
  h = WLZ_RGBA_GREEN_GET(highVal.v.rgbv);
  if( l > h ){
    WLZ_RGBA_GREEN_SET(highVal.v.rgbv, l);
    WLZ_RGBA_GREEN_SET(lowVal.v.rgbv, h);
  }
  l = WLZ_RGBA_BLUE_GET(lowVal.v.rgbv);
  h = WLZ_RGBA_BLUE_GET(highVal.v.rgbv);
  if( l > h ){
    WLZ_RGBA_BLUE_SET(highVal.v.rgbv, l);
    WLZ_RGBA_BLUE_SET(lowVal.v.rgbv, h);
  }

  /* set all AND for combine value and call multi-threshold */
  WLZ_RGBA_RGBA_SET(combineMode,
		    WLZ_BO_AND, WLZ_BO_AND, WLZ_BO_AND, 255);
  rtnObj = WlzRGBAMultiThreshold(obj, lowVal, highVal,
				 combineMode, &errNum);

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

WlzObject *WlzRGBAEllipsoidThreshold(
  WlzObject	*obj,
  WlzPixelV	lowVal,
  WlzPixelV	highVal,
  double	eccentricity,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzObject	*rtnObj=NULL;
  WlzVectorThresholdStruct clientData;
  double	x, y, z, s;

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  
  if( errNum == WLZ_ERR_NONE ){

    /* set up the callback data structure */
    clientData.type = WLZ_RGBA_THRESH_SPHERE;
    clientData.origin[0] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_RED,
					    &errNum);
    clientData.origin[1] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_GREEN,
					    &errNum);
    clientData.origin[2] = WlzRGBAPixelValue(lowVal, WLZ_RGBA_CHANNEL_BLUE,
					    &errNum);
    x = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_RED, &errNum) -
      clientData.origin[0];
    y = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_GREEN, &errNum) -
      clientData.origin[1];
    z = WlzRGBAPixelValue(highVal, WLZ_RGBA_CHANNEL_BLUE, &errNum) -
      clientData.origin[2];
    s = sqrt(x*x + y*y + z*z + 0.0000001);
    if( s < 0.001 ){
      x = 1.0; y = z = 0.0;
    }
    else {
      x /= s; y /= s; z /= s;
    }
    clientData.direction[0] = x;
    clientData.direction[1] = y;
    clientData.direction[2] = z;
    clientData.dist = s*s;
    clientData.eccen = eccentricity;
    
    /* now call the threshold function */
    rtnObj = WlzCbThreshold(obj, WlzVectorThreshCb, &clientData, &errNum);
  }

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

