#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzRGBAConvert.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Wed May 14 07:47:16 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzValuesUtils
* \brief        Conversion routines for RGBA data including conversion to modulus.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
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

static WlzObject *WlzRGBAToModulus3D(
  WlzObject	*obj,
  WlzErrorNum	*dstErr);

static void WlzRGBAConvertRGBToHSV_UBYTENormalised(
  int		*col)
{
  int	h, s, b;
  int	max, min;

  /* algorithm from Foley, van Dam, Feiner, Hughes,
     Computer Graphics. Modified so each value is in the range
     [0,255]. If saturation is zero thgen the hue is undefined
     and in this case set to zero */
  max = WLZ_MAX(col[0],col[1]);
  max = WLZ_MAX(max, col[2]);
  min = WLZ_MIN(col[0],col[1]);
  min = WLZ_MIN(min, col[2]);

  b = max;
  if( max > 0 ){
    s = (max - min) * 255 / max;
  }
  else {
    s = 0;
  }
  if( s == 0 ){
    h = 0;
  }
  else {
    if( col[0] == max ){
      h = (col[1] - col[2]) * 42.5 / (max - min);
    }
    else if( col[1] == max ){
      h = 85 + (col[2] - col[0]) * 42.5 / (max - min);
    }
    else if( col[2] == max ){
      h = 170 + (col[0] - col[1]) * 42.5 / (max - min);
    }
    if( h < 0 ){
      h += 255;
    }
  }
  col[0] = h;
  col[1] = s;
  col[2] = b;

  return;
}



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
