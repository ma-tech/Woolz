#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzRGBAConvert.c
* \author       Richard Baldock
* \date         May 2003
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
* \brief	Conversion routines for RGBA data including conversion
* 		to modulus.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/* static functions defined later */
static WlzCompoundArray *WlzRGBAToCompound3D(
  WlzObject	*obj,
  WlzRGBAColorSpace	colSpc,
  WlzErrorNum	*dstErr);

static WlzObject *WlzCompoundToRGBA3D(
  WlzCompoundArray	*cmpnd,
  WlzRGBAColorSpace	colSpc,
  int			clipFlg,
  WlzErrorNum		*dstErr);

static WlzObject *WlzIndexToRGBA3D(
  WlzObject	*obj,
  unsigned char	colormap[3][256],
  WlzErrorNum	*dstErr);

static WlzObject *WlzRGBAToModulus3D(
  WlzObject	*obj,
  WlzErrorNum	*dstErr);

static WlzObject *WlzRGBAToChannel3D(
  WlzObject	*obj,
  WlzRGBAColorChannel	chan,
  WlzErrorNum	*dstErr);

/* function:     WlzRGBAToCompound    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Convert a RGBA image to a compound object. The RGBA
 channels are at array indices 0,1,2,3 respectively and the sub-object
 grey types will be WLZ_GREY_UBYTE.
*
* \return       Compound array of rgba values
* \param    obj	Input domain object with value type WLZ_GREY_RGBA
* \param    colSpc	The colour space.
* \param    dstErr	error return
* \par      Source:
*                WlzRGBAConvert.c
*/
WlzCompoundArray *WlzRGBAToCompound(
  WlzObject	*obj,
  WlzRGBAColorSpace	colSpc,
  WlzErrorNum	*dstErr)
{
  WlzCompoundArray	*cobj=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object type, and value type */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      return WlzRGBAToCompound3D(obj, colSpc, dstErr);

    case WLZ_TRANS_OBJ:
      /* not difficult, do it later */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      /* bit recursive this ! */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_EMPTY_OBJ:
      return (WlzCompoundArray *) WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* check colour space */
  if( errNum == WLZ_ERR_NONE ){
    switch( colSpc ){
    case WLZ_RGBA_SPACE_RGB:
    case WLZ_RGBA_SPACE_HSB:
    case WLZ_RGBA_SPACE_CMY:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* create compound object return */
  if( errNum == WLZ_ERR_NONE ){
    WlzValues	values;
    WlzObject	*objs[4];
    WlzObjectType	type;
    WlzPixelV	oldBck, newBck;

    type = WlzGreyTableType(
      WlzGreyTableTypeToTableType(obj->values.core->type, &errNum),
      WLZ_GREY_UBYTE, &errNum);
    oldBck = WlzGetBackground(obj, &errNum);

    /* red */
    newBck.type = WLZ_GREY_UBYTE;
    newBck.v.ubv = WLZ_RGBA_RED_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[0] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* green */
    newBck.v.ubv = WLZ_RGBA_GREEN_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[1] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* blue */
    newBck.v.ubv = WLZ_RGBA_BLUE_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[2] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* alpha */
    newBck.v.ubv = WLZ_RGBA_ALPHA_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[3] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);

    /* create compound object, object pointers are assigned for mode=3
       so no need to free objects */
    cobj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 4, &(objs[0]),
				 obj->type, &errNum);
  }

  /* iterate through objects setting values */
  if( errNum == WLZ_ERR_NONE ){
    WlzIntervalWSpace	iwsp0, iwsp[4];
    WlzGreyWSpace	gwsp0, gwsp[4];
    int			i, j, k;
    int			a, col[3];

    errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0);
    for(i=0; i < 4; i++){
      errNum = WlzInitGreyScan(cobj->o[i], &(iwsp[i]), &(gwsp[i]));
    }
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
      for(i=0; i < 4; i++){
	errNum = WlzNextGreyInterval(&(iwsp[i]));
      }
      switch( colSpc ){
      case WLZ_RGBA_SPACE_RGB:
	for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	      gwsp0.u_grintptr.rgbp++){
	  *(gwsp[0].u_grintptr.ubp++) =
	    WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[1].u_grintptr.ubp++) =
	    WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[2].u_grintptr.ubp++) =
	    WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[3].u_grintptr.ubp++) =
	    WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
	}
	break;
      
      case WLZ_RGBA_SPACE_HSB: /* each normalised to [0,255] */
	for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	      gwsp0.u_grintptr.rgbp++){
	  col[0] = WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	  col[1] = WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	  col[2] = WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	  a = WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
	  WlzRGBAConvertRGBToHSV_UBYTENormalised(col);
	  *(gwsp[0].u_grintptr.ubp++) = col[0];
	  *(gwsp[1].u_grintptr.ubp++) = col[1];
	  *(gwsp[2].u_grintptr.ubp++) = col[2];
	  *(gwsp[3].u_grintptr.ubp++) = a;
	    }
	break;
      
      case WLZ_RGBA_SPACE_CMY:
	for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	      gwsp0.u_grintptr.rgbp++){
	  col[0] = WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	  col[1] = WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	  col[2] = WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	  a = WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[0].u_grintptr.ubp++) = (col[1] + col[2]) / 2;
	  *(gwsp[1].u_grintptr.ubp++) = (col[2] + col[0]) / 2;
	  *(gwsp[2].u_grintptr.ubp++) = (col[0] + col[1]) / 2;
	  *(gwsp[3].u_grintptr.ubp++) = a;
	    }
	break;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }
      

  if( dstErr ){
    *dstErr = errNum;
  }
  return cobj;
}

static WlzCompoundArray *WlzRGBAToCompound3D(
  WlzObject	*obj,
  WlzRGBAColorSpace	colSpc,
  WlzErrorNum	*dstErr)
{
  WlzCompoundArray	*cobj=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  
  if( dstErr ){
    *dstErr = errNum;
  }
  return NULL;
}

WlzObject *WlzCompoundToRGBA(
  WlzCompoundArray	*cmpnd,
  WlzRGBAColorSpace	colSpc,
  int			clipFlg,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzPixelV	bckgrnd;
  int		i, j;
  UINT		b[4];
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object - must have at least 3 grey-level objects
     all of the same type. A fourth is used as the alpha-channel
     which is otherwise set to 255.
     The union of domains is returned with values external to
     any input domain set using the respective background */
  if( cmpnd == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( cmpnd->n < 3 ){
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else {
    switch( cmpnd->o[0]->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if((cmpnd->o[1]->type != cmpnd->o[0]->type) ||
	 (cmpnd->o[2]->type != cmpnd->o[0]->type)){
	errNum = WLZ_ERR_OBJECT_DATA;
      }
      break;

    default:
      errNum = WLZ_ERR_OBJECT_DATA;
      break;
    }
  }

  /* check the colour space parameter */
  if( errNum == WLZ_ERR_NONE ){
    switch( colSpc ){
    case WLZ_RGBA_SPACE_RGB:
    case WLZ_RGBA_SPACE_HSB:
    case WLZ_RGBA_SPACE_CMY:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* check for 3D */
  if( errNum == WLZ_ERR_NONE ){
    if( cmpnd->o[0]->type == WLZ_3D_DOMAINOBJ ){
      return WlzCompoundToRGBA3D(cmpnd, colSpc, clipFlg, dstErr);
    }
  }

  /* 2D case */
  if( errNum == WLZ_ERR_NONE ){
    /* build the 2D object - union of the rgb channels */
    if( rtnObj = WlzUnionN(3, cmpnd->o, 0, &errNum) ){
      WlzValues	values;
      WlzObjectType	vType;

      /* add an RGBA valuetable, extract background for each channel */
      vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_RGBA, NULL);
      for(i=0; i < 3; i++){
	bckgrnd = WlzGetBackground(cmpnd->o[i], NULL);
	WlzValueConvertPixel(&bckgrnd, bckgrnd, WLZ_GREY_UBYTE);
	b[i] = bckgrnd.v.ubv;
      }
      bckgrnd.type = WLZ_GREY_RGBA;
      WLZ_RGBA_RGBA_SET(bckgrnd.v.rgbv, b[0], b[1], b[2], 255);
      if( values.v = WlzNewValueTb(rtnObj, vType, bckgrnd, &errNum) ){
	rtnObj->values = WlzAssignValues(values, &errNum);
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }

  /* transfer values */
  if( errNum == WLZ_ERR_NONE ){
    WlzGreyValueWSpace	*gValWSpc[4];
    WlzIntervalWSpace	iwsp;
    WlzGreyWSpace	gwsp;
    WlzGreyV		gval;

    /* do it dumb fashion for now, rgb only */
    for(i=0; i < 3; i++){
      gValWSpc[i] = WlzGreyValueMakeWSp(cmpnd->o[i], &errNum);
    }
    
    errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
    while((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE){
      WlzPixelV	pix;

      for(j = iwsp.lftpos; j <= iwsp.rgtpos; j++){
	for(i=0; i < 3; i++){
	  WlzGreyValueGet(gValWSpc[i], 0, iwsp.linpos, j);
	  pix.type = gValWSpc[i]->gType;
	  pix.v = gValWSpc[i]->gVal[0];
	  WlzValueConvertPixel(&pix, pix, WLZ_GREY_UBYTE);
	  b[i] = pix.v.ubv;
	}
	WLZ_RGBA_RGBA_SET(gval.rgbv,
			  b[0], b[1], b[2], 255);
	*gwsp.u_grintptr.rgbp = gval.rgbv;
	gwsp.u_grintptr.rgbp++;
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }

    for(i=0; i < 3; i++){
      if( gValWSpc[i] ){
	WlzGreyValueFreeWSp(gValWSpc[i]);
      }
    }
  }

  /* set error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

static WlzObject *WlzCompoundToRGBA3D(
  WlzCompoundArray	*cmpnd,
  WlzRGBAColorSpace	colSpc,
  int			clipFlg,
  WlzErrorNum		*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/* function:     WlzRGBAToModulus    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Calculate the modulus of the rgb values and return
 in an image of grey type WLZ_GREY_SHORT
*
* \return       Grey-level object of modulus values
* \param    obj	Input rgba object
* \param    dstErr	error return
* \par      Source:
*                WlzRGBAConvert.c
*/
WlzObject *WlzRGBAToModulus(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object type, and value type */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      return WlzRGBAToModulus3D(obj, dstErr);

    case WLZ_TRANS_OBJ:
      /* not difficult, do it later */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      /* bit recursive this ! */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* create object return */
  if( errNum == WLZ_ERR_NONE ){
    WlzValues	values;
    WlzObjectType	type;
    WlzPixelV	oldBck, newBck;

    type = WlzGreyTableType(
      WlzGreyTableTypeToTableType(obj->values.core->type, &errNum),
      WLZ_GREY_SHORT, &errNum);
    oldBck = WlzGetBackground(obj, &errNum);
    newBck.type = WLZ_GREY_SHORT;
    newBck.v.shv = WLZ_RGBA_MODULUS(oldBck.v.rgbv);

    /* make values table and return object */
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    rtnObj = WlzMakeMain(obj->type, obj->domain, values,
			 NULL, NULL, &errNum);
  }

  /* iterate through objects setting values */
  if( errNum == WLZ_ERR_NONE ){
    WlzIntervalWSpace	iwsp0, iwsp1;
    WlzGreyWSpace	gwsp0, gwsp1;
    int			i, j, k;

    errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0);
    errNum = WlzInitGreyScan(rtnObj, &iwsp1, &gwsp1);
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
      errNum = WlzNextGreyInterval(&iwsp1);
      for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	    gwsp0.u_grintptr.rgbp++){
	*(gwsp1.u_grintptr.shp++) = WLZ_RGBA_MODULUS(*(gwsp0.u_grintptr.rgbp));
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

static WlzObject *WlzRGBAToModulus3D(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}


/* function:     WlzIndexToRGBA    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Convert a grey-level woolz object to RGBA via
 a colourmap look-up table. Values are clamped to [0,255] and the
 LUT is assumed to be unsigned byte 3x256 
*
* \return       Woolz object
* \param    obj	Input object to be converted
* \param    colormap	Colourmap array
* \param    dstErr	Error return
* \par      Source:
*                WlzRGBAConvert.c
*/
WlzObject *WlzIndexToRGBA(
  WlzObject	*obj,
  unsigned char	colormap[3][256],
  WlzErrorNum	*dstErr)
{
  WlzObject		*rtnObj=NULL;
  WlzGreyType		oldpixtype;
  WlzGreyP		go, gn;
  WlzIntervalWSpace	oldiwsp, newiwsp;
  WlzGreyWSpace		oldgwsp, newgwsp;
  WlzObjectType		newvtbltype;
  WlzPixelV		bg;
  WlzValues		newvalues,
  			values;
  int 			k, greyVal, redVal, greenVal, blueVal;
  unsigned int		rgbaVal;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object - must be domain object with a values table */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core ){
	if( obj->values.core == NULL ){
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else {
	  oldpixtype =
	    WlzGreyTableTypeToGreyType(obj->values.core->type, NULL);
	  if( oldpixtype == WLZ_GREY_RGBA ){
	    return WlzMakeMain(obj->type, obj->domain, obj->values,
			       NULL, NULL, dstErr);
	  }
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzIndexToRGBA3D(obj, colormap, dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /*
   * Set type of new value table so as to preserve
   * rectangular/single interval/multiple interval
   * type.
   */
  if( errNum == WLZ_ERR_NONE ){
    newvtbltype = WlzGreyTableTypeToTableType(obj->values.core->type,
					      &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    newvtbltype = WlzGreyTableType(newvtbltype, WLZ_GREY_RGBA, &errNum);
  }

  /* get the background  - note background now carries its own type */
  if( errNum == WLZ_ERR_NONE ){
    bg = WlzGetBackground(obj, &errNum);
    switch( bg.type ){
    case WLZ_GREY_INT:
      greyVal = WLZ_CLAMP(bg.v.inv, 0, 255);
      break;

    case WLZ_GREY_SHORT:
      greyVal = WLZ_CLAMP(bg.v.shv, 0, 255);
      break;

    case WLZ_GREY_UBYTE:
      greyVal = bg.v.ubv;
      break;

    case WLZ_GREY_FLOAT:
      greyVal = WLZ_CLAMP(bg.v.flv, 0, 255);
      break;

    case WLZ_GREY_DOUBLE:
      greyVal = WLZ_CLAMP(bg.v.dbv, 0, 255);
      break;

    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    }
    bg.type = WLZ_GREY_RGBA;
    WLZ_RGBA_RGBA_SET(bg.v.rgbv, colormap[0][greyVal],
		      colormap[1][greyVal], colormap[2][greyVal],
		      255);
  }

  /*
   * Make the new object with new value table type and value table
   * allocated (but blank). Share original idom.
   */
  if( errNum == WLZ_ERR_NONE ){
    values.v = WlzNewValueTb(obj, newvtbltype, bg, &errNum);
  }
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj->domain, values,
			 obj->plist, obj->assoc, &errNum);
  }

  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(obj, &oldiwsp, &oldgwsp);
  }
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(rtnObj, &newiwsp, &newgwsp);
  }

  while( ((errNum = WlzNextGreyInterval(&oldiwsp)) == WLZ_ERR_NONE)
	 && ((errNum = WlzNextGreyInterval(&newiwsp)) == WLZ_ERR_NONE) ){
    go = oldgwsp.u_grintptr;	
    gn = newgwsp.u_grintptr;

    for(k=0; k <= oldiwsp.colrmn; k++){
      switch( oldgwsp.pixeltype ){
      case WLZ_GREY_INT:
	greyVal = WLZ_CLAMP(*(go.inp), 0, 255);
	go.inp++;
	break;

      case WLZ_GREY_SHORT:
	greyVal = WLZ_CLAMP(*(go.shp), 0, 255);
	go.shp++;
	break;

      case WLZ_GREY_UBYTE:
	greyVal = *(go.ubp);
	go.ubp++;
	break;

      case WLZ_GREY_FLOAT:
	greyVal = WLZ_CLAMP(*(go.flp), 0, 255);
	go.flp++;
	break;

      case WLZ_GREY_DOUBLE:
	greyVal = WLZ_CLAMP(*(go.dbp), 0, 255);
	go.dbp++;
	break;

      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
      }
      redVal = colormap[0][greyVal];
      greenVal = colormap[1][greyVal];
      blueVal = colormap[2][greyVal];
      WLZ_RGBA_RGBA_SET(rgbaVal, redVal, greenVal, blueVal, 0xff);
      *(gn.rgbp) = rgbaVal;
      gn.rgbp++;
    }
  } /* while */
  if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
  {
    errNum = WLZ_ERR_NONE;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

static WlzObject *WlzIndexToRGBA3D(
  WlzObject	*obj,
  unsigned char	colormap[3][256],
  WlzErrorNum	*dstErr)
{
  WlzObject		*obj1, *temp;
  WlzPlaneDomain	*pdom, *npdom;
  WlzVoxelValues	*voxtab, *nvoxtab;
  WlzDomain		*domains, *ndomains, domain;
  WlzValues		*values, *nvalues, vals;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* no need to check the object pointer or type because this procedure
     can only be accessed via WlzIndexToRGBA. The domain and valuetable
     must be checked however */
  obj1 = NULL;
  if( obj->domain.p == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->values.vox == NULL ){
    errNum = WLZ_ERR_VALUES_NULL;
  }

  /* check types */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->domain.p->type ){

    case WLZ_PLANEDOMAIN_DOMAIN:
      break;

    default:
      errNum = WLZ_ERR_PLANEDOMAIN_TYPE;
      break;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    switch( obj->values.vox->type ){

    case WLZ_VOXELVALUETABLE_GREY:
      break;

    default:
      errNum = WLZ_ERR_VOXELVALUES_TYPE;
      break;
    }
  }

  /* make new planedomain and voxelvaluetable */
  if( errNum == WLZ_ERR_NONE ){
    pdom = obj->domain.p;
    voxtab = obj->values.vox;
    npdom = WlzMakePlaneDomain(pdom->type,
			       pdom->plane1, pdom->lastpl,
			       pdom->line1, pdom->lastln,
			       pdom->kol1, pdom->lastkl, &errNum);
  }
    
  if((errNum == WLZ_ERR_NONE) &&
     ((nvoxtab = WlzMakeVoxelValueTb(voxtab->type, voxtab->plane1,
				     voxtab->lastpl, voxtab->bckgrnd,
				     NULL, &errNum)) == NULL) ){
    WlzFreePlaneDomain(npdom);
  }

  if( errNum == WLZ_ERR_NONE ){
    /* set up variables */
    domains = pdom->domains;
    ndomains = npdom->domains;
    values = voxtab->values;
    nvalues = nvoxtab->values;
    nplanes = pdom->lastpl - pdom->plane1 + 1;

    /* copy voxel_sizes */
    for(i=0; i < 3; i++){
      npdom->voxel_size[i] = pdom->voxel_size[i];
    }
  }

  /* convert each plane */
#ifdef _OPENMP
  {
    int		pnIdx;
    WlzErrorNum pvErrNum;
    /* (*(ndomains + pnIdx)).core and (*(nvalues + pnIdx)).core are
     * initialized to NULL by WlzMakePlaneDomain() and WlzMakeVoxelValueTb()
     * when they are created. */
    #pragma omp parallel for default(shared) private(pvErrNum,temp,obj1)
    for(pnIdx = 0; pnIdx < nplanes; ++pnIdx)
    {
      if((errNum == WLZ_ERR_NONE) &&
         (*(domains + pnIdx)).core && (*(values + pnIdx)).core)
      {
	temp = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			   *(domains + pnIdx), *(values + pnIdx),
			   NULL, NULL, &pvErrNum);
	obj1 = NULL;
	if(temp && temp->domain.core)
	{
	  obj1 = WlzIndexToRGBA(temp, colormap, &pvErrNum);
	}
	if(obj1)
	{
	  if( obj1->type == WLZ_2D_DOMAINOBJ ){
	    *(ndomains + pnIdx) = WlzAssignDomain(obj1->domain, NULL);
	    *(nvalues + pnIdx) = WlzAssignValues(obj1->values, NULL);
	  }
	  else {
	    (*(ndomains + pnIdx)).core = NULL;
	    (*(nvalues + pnIdx)).core = NULL;
	  }
	  WlzFreeObj(obj1);
	}
	#pragma omp critical
	{
	  if((pvErrNum != WLZ_ERR_NONE) && (errNum == WLZ_ERR_NONE))
	  {
	    errNum = pvErrNum;
	  }
	}
      }
    }
  }
#else
  while( (errNum == WLZ_ERR_NONE) && nplanes-- ){
    if(((*domains).core == NULL) || ((*values).core == NULL)){
      (*ndomains).core = NULL;
      (*nvalues).core = NULL;
    }
    else if( temp = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				NULL, NULL, &errNum) ){

      if( temp->domain.i != NULL ){
	if( obj1 = WlzIndexToRGBA(temp, colormap, &errNum) ){
	  if( obj1->type == WLZ_2D_DOMAINOBJ ){
	    *ndomains = WlzAssignDomain(obj1->domain, NULL);
	    *nvalues = WlzAssignValues(obj1->values, NULL);
	  }
	  else {
	    (*ndomains).core = NULL;
	    (*nvalues).core = NULL;
	  }
	  WlzFreeObj(obj1);
	}
      } else {
	(*ndomains).core = NULL;
	(*nvalues).core = NULL;
      }
    }
    else {
      WlzFreePlaneDomain(npdom);
      WlzFreeVoxelValueTb( nvoxtab );
      break;
    }

    domains++;
    ndomains++;
    values++;
    nvalues++;
  }
#endif
  /* return a new object */
  if( errNum == WLZ_ERR_NONE ){
    domain.p = npdom;
    vals.vox = nvoxtab;
    if( obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			   NULL, obj, &errNum) ){
      /*	nvoxtab->original = obj1; */
      nvoxtab->original_table = WlzAssignValues(obj->values, NULL);
    }
    else {
      WlzFreePlaneDomain( npdom );
      WlzFreeVoxelValueTb( nvoxtab );
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return obj1;
}

WlzObject *WlzRGBAToChannel(
  WlzObject	*obj,
  WlzRGBAColorChannel	chan,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzValues		values;
  WlzObjectType	type;
  WlzPixelV		pixVal, oldBck, newBck;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object and channel */
  if( obj ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else if ( obj->values.core == NULL ){
	errNum = WLZ_ERR_VALUES_NULL;
      }
      else if( obj->values.vox->type != WLZ_VOXELVALUETABLE_GREY ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if( WlzGreyTypeFromObj(obj, &errNum) != WLZ_GREY_RGBA ){
	errNum = WLZ_ERR_VALUES_TYPE;
      }
      return WlzRGBAToChannel3D(obj, chan, dstErr);

    case WLZ_TRANS_OBJ:
      /* not difficult, do it later */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      /* bit recursive this ! */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( chan ){
    case WLZ_RGBA_CHANNEL_RED:
    case WLZ_RGBA_CHANNEL_GREEN:
    case WLZ_RGBA_CHANNEL_BLUE:
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
    case WLZ_RGBA_CHANNEL_BRIGHTNESS:
    case WLZ_RGBA_CHANNEL_CYAN:
    case WLZ_RGBA_CHANNEL_MAGENTA:
    case WLZ_RGBA_CHANNEL_YELLOW:
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* now extract data */
  if( errNum == WLZ_ERR_NONE ){

    type = WlzGreyTableType(
      WlzGreyTableTypeToTableType(obj->values.core->type, &errNum),
      WLZ_GREY_UBYTE, &errNum);
    oldBck = WlzGetBackground(obj, &errNum);
    newBck.type = WLZ_GREY_UBYTE;
    newBck.v.ubv = (UBYTE) WlzRGBAPixelValue(oldBck, chan, &errNum);
  }

  /* make values table and return object */
  if( errNum == WLZ_ERR_NONE ){
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    rtnObj = WlzMakeMain(obj->type, obj->domain, values,
			 NULL, NULL, &errNum);
  }

  /* iterate through objects setting values */
  if( errNum == WLZ_ERR_NONE ){
    WlzIntervalWSpace	iwsp0, iwsp1;
    WlzGreyWSpace	gwsp0, gwsp1;
    int			i, j, k;

    errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0);
    errNum = WlzInitGreyScan(rtnObj, &iwsp1, &gwsp1);
    pixVal.type = WLZ_GREY_RGBA;
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
      errNum = WlzNextGreyInterval(&iwsp1);
      for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	    gwsp0.u_grintptr.rgbp++){
	pixVal.v.rgbv = (*(gwsp0.u_grintptr.rgbp));
	*(gwsp1.u_grintptr.ubp++) = (UBYTE)
	  WlzRGBAPixelValue(pixVal, chan, &errNum);
      }
    }
    if( errNum == WLZ_ERR_EOO ){
      errNum = WLZ_ERR_NONE;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}


static WlzObject *WlzRGBAToChannel3D(
  WlzObject	*obj,
  WlzRGBAColorChannel	chan,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
