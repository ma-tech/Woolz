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
  WlzErrorNum	*dstErr);

static WlzObject *WlzRGBAToModulus3D(
  WlzObject	*obj,
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
* \param    dstErr	error return
* \par      Source:
*                WlzRGBAConvert.c
*/
WlzCompoundArray *WlzRGBAToCompound(
  WlzObject	*obj,
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
      return WlzRGBAToCompound3D(obj, dstErr);

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

    errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0);
    for(i=0; i < 4; i++){
      errNum = WlzInitGreyScan(cobj->o[i], &(iwsp[i]), &(gwsp[i]));
    }
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
      for(i=0; i < 4; i++){
	errNum = WlzNextGreyInterval(&(iwsp[i]));
      }
      for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	    gwsp0.u_grintptr.rgbp++){
	*(gwsp[0].u_grintptr.ubp++) = WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	*(gwsp[1].u_grintptr.ubp++) = WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	*(gwsp[2].u_grintptr.ubp++) = WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	*(gwsp[3].u_grintptr.ubp++) = WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
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
  WlzErrorNum	*dstErr)
{
  WlzCompoundArray	*cobj=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  
  if( dstErr ){
    *dstErr = errNum;
  }
  return NULL;
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
