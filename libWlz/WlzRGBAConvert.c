#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBAConvert_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRGBAConvert.c
* \author       Richard Baldock
* \date         May 2003
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
* \brief	Conversion routines for RGBA data including conversion
* 		to modulus.
* \ingroup	WlzValuesUtils
*/

#include <Wlz.h>

/* static functions defined later */
static WlzCompoundArray 	*WlzRGBAToCompound3D(
  				  WlzObject *obj,
  				  WlzRGBAColorSpace colSpc,
  				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCompoundToRGBA2D(
  				  WlzCompoundArray *cmpnd,
  				  WlzRGBAColorSpace colSpc,
  				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCompoundToRGBA3D(
  				  WlzCompoundArray *cmpnd,
  				  WlzRGBAColorSpace colSpc,
  				  WlzErrorNum *dstErr);
static WlzObject 		*WlzIndexToRGBA3D(
  				  WlzObject *obj,
				  unsigned char colormap[3][256],
				  WlzErrorNum *dstErr);

static WlzObject		*WlzRGBAToModulus3D(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);

static WlzObject		*WlzRGBAToChannel3D(
				  WlzObject *obj,
				  WlzRGBAColorChannel chan,
				  WlzErrorNum *dstErr);

/*! 
* \ingroup      WlzValuesUtils
* \brief        Convert a RGBA image to a compound object. The RGBA
*		 channels are at array indices 0,1,2,3 respectively
*		 and the sub-object grey types will be WLZ_GREY_UBYTE.
*
* \return       Compound array of rgba values
* \param    obj				Input domain object with value type
* 					WLZ_GREY_RGBA
* \param    colSpc			The colour space.
* \param    dstErr			Destination error ponyer, may be NULL.
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
      else if (WlzGreyTableIsTiled(obj->values.core->type)) {
	errNum = WLZ_ERR_VALUES_TYPE;
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
    newBck.v.ubv = (WlzUByte )WLZ_RGBA_RED_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[0] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* green */
    newBck.v.ubv = (WlzUByte )WLZ_RGBA_GREEN_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[1] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* blue */
    newBck.v.ubv = (WlzUByte )WLZ_RGBA_BLUE_GET(oldBck.v.rgbv);
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    objs[2] = WlzMakeMain(obj->type, obj->domain, values,
			  NULL, NULL, &errNum);
    /* alpha */
    newBck.v.ubv = (WlzUByte )WLZ_RGBA_ALPHA_GET(oldBck.v.rgbv);
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
	    (WlzUByte )WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[1].u_grintptr.ubp++) =
	    (WlzUByte )WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[2].u_grintptr.ubp++) =
	    (WlzUByte )WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[3].u_grintptr.ubp++) =
	    (WlzUByte )WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
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
	  *(gwsp[0].u_grintptr.ubp++) = (WlzUByte )(col[0]);
	  *(gwsp[1].u_grintptr.ubp++) = (WlzUByte )(col[1]);
	  *(gwsp[2].u_grintptr.ubp++) = (WlzUByte )(col[2]);
	  *(gwsp[3].u_grintptr.ubp++) = (WlzUByte )a;
	    }
	break;
      
      case WLZ_RGBA_SPACE_CMY:
	for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	      gwsp0.u_grintptr.rgbp++){
	  col[0] = WLZ_RGBA_RED_GET(*(gwsp0.u_grintptr.rgbp));
	  col[1] = WLZ_RGBA_GREEN_GET(*(gwsp0.u_grintptr.rgbp));
	  col[2] = WLZ_RGBA_BLUE_GET(*(gwsp0.u_grintptr.rgbp));
	  a = WLZ_RGBA_ALPHA_GET(*(gwsp0.u_grintptr.rgbp));
	  *(gwsp[0].u_grintptr.ubp++) = (WlzUByte )((col[1] + col[2]) / 2);
	  *(gwsp[1].u_grintptr.ubp++) = (WlzUByte )((col[2] + col[0]) / 2);
	  *(gwsp[2].u_grintptr.ubp++) = (WlzUByte )((col[0] + col[1]) / 2);
	  *(gwsp[3].u_grintptr.ubp++) = (WlzUByte )a;
	    }
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
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
  WlzErrorNum		errNum = WLZ_ERR_UNIMPLEMENTED;
  
  if( dstErr ){
    *dstErr = errNum;
  }
  return NULL;
}

/*!
* \return	New object with WLZ_GREY_RGBA values or possible an empty
* 		object.
* \ingroup      WlzValuesUtils
* \brief	Creates a WLZ_GREY_RGBA valued object from the given compound
* 		array. This is a static function which will always be called
* 		with valid parameters so they aren't checked. If all members of
* 		the compound array are empty then the returned object will be
* 		empty too.
* \param	cObj			Compound array object.
* \param	cSpc 			The colour space.
* \param	dstErr			Destination error pointer may be NULL.
*/
WlzObject	*WlzCompoundToRGBA(WlzCompoundArray *cmpnd,
				WlzRGBAColorSpace colSpc,
  				WlzErrorNum *dstErr)
{
  int		i;
  WlzObject	*rObj = NULL;
  WlzObjectType	oType = WLZ_EMPTY_OBJ;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* Check the compound array type and it's object types which must all
   * be either WLZ_2D_DOMAINOBJ and WLZ_EMPTY_OBJ or WLZ_3D_DOMAINOBJ
   * and WLZ_EMPTY_OBJ. If the number of channels is less than four
   * then the components will be assigned in order of R, G, B and then
   * A. Components not supplied will be treated as if they were WLZ_EMPTY_OBJ. */
  if(cmpnd == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cmpnd->n < 3)
  {
    errNum = WLZ_ERR_OBJECT_DATA;
  }
  else
  {
    for(i = 0; i < cmpnd->n; ++i)
    {
      WlzObject *o;
      if((o = cmpnd->o[i]) != NULL)
      {
	switch(o->type)
	{
	  case WLZ_EMPTY_OBJ:
	    break;
	  case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
	  case WLZ_3D_DOMAINOBJ:
	    if(oType == WLZ_EMPTY_OBJ)
	    {
	      oType = o->type;
	    }
	    else if(oType != o->type)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
    }
  }
  /* Check the colour space parameter */
  if(errNum == WLZ_ERR_NONE)
  {
    switch(colSpc)
    {
      case WLZ_RGBA_SPACE_RGB: /* FALLTHROUGH */
      case WLZ_RGBA_SPACE_HSB: /* FALLTHROUGH */
      case WLZ_RGBA_SPACE_CMY:
	break;
      default:
	errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(oType)
    {
      case WLZ_EMPTY_OBJ:
        rObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
        rObj =  WlzCompoundToRGBA2D(cmpnd, colSpc, &errNum);
        break;
      case WLZ_3D_DOMAINOBJ:
        rObj =  WlzCompoundToRGBA3D(cmpnd, colSpc, &errNum);
        break;
      default:
        break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New 3D domain object with corresponding WLZ_GREY_RGBA values.
* \ingroup      WlzValuesUtils
* \brief	Creates a WLZ_GREY_RGBA valued object from the given compound
* 		array. This is a static function which will always be called
* 		with valid parameters so they aren't checked.
* \param	cObj			Compound array object.
* \param	cSpc 			The colour space.
* \param	dstErr			Destination error pointer may be NULL.
*/
static WlzObject *WlzCompoundToRGBA2D(WlzCompoundArray *cObj,
  				WlzRGBAColorSpace cSpc, WlzErrorNum *dstErr)
{
  int		i,
  		j,
		numObjs = 3;
  WlzObject	*rtnObj=NULL;
  WlzPixelV	bckgrnd;
  WlzObject	*objs[4];
  WlzObjectType vType;
  WlzUInt	b[4];
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* Make a copy of the object pointers because WlzUnionN() modifies the
   * array if it contains empty objects. */
  if(cObj->n >= 4)
  {
    numObjs = 4;
  }
  for(i = 0; i < numObjs; ++i)
  {
    objs[i] = cObj->o[i];
  }
  rtnObj = WlzUnionN(numObjs, objs, 0, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Add an RGBA valuetable, extract background for each channel */
    vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_RGBA, &errNum);
    b[0] = b[1] = b[2] = b[3] = 255;
    for(i=0; (errNum == WLZ_ERR_NONE) && (i < numObjs); i++)
    {
      bckgrnd = WlzGetBackground(cObj->o[i], &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzValueConvertPixel(&bckgrnd, bckgrnd, WLZ_GREY_UBYTE);
        b[i] = bckgrnd.v.ubv;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	values;

    bckgrnd.type = WLZ_GREY_RGBA;
    WLZ_RGBA_RGBA_SET(bckgrnd.v.rgbv, b[0], b[1], b[2], b[3]);
    values.v = WlzNewValueTb(rtnObj, vType, bckgrnd, &errNum);
    if(values.v != NULL)
    {
      rtnObj->values = WlzAssignValues(values, &errNum);
    }
    else
    {
      (void )WlzFreeObj(rtnObj);
      rtnObj = NULL;
    }
  }
  /* Transfer values */
  if( errNum == WLZ_ERR_NONE)
  {
    WlzGreyValueWSpace	*gValWSpc[4];
    WlzIntervalWSpace	iwsp = {0};
    WlzGreyWSpace	gwsp;
    WlzGreyV		gval;

    /* do it dumb fashion for now, rgb only */
    gValWSpc[0] = gValWSpc[1] = gValWSpc[2] = gValWSpc[3] = NULL;
    for(i=0; i < numObjs; i++)
    {
      if((cObj->o[i] != NULL) && (cObj->o[i]->type != WLZ_EMPTY_OBJ))
      {
        gValWSpc[i] = WlzGreyValueMakeWSp(cObj->o[i], &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
    }
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE))
    {
      WlzPixelV	pix;

      for(j = iwsp.lftpos; j <= iwsp.rgtpos; j++)
      {
	for(i = 0; i < numObjs; i++)
	{
	  if(gValWSpc[i] == NULL)
	  {
	    pix.v.ubv = (i < 3)? 0: 255;
	  }
	  else
	  {
	    WlzGreyValueGet(gValWSpc[i], 0, iwsp.linpos, j);
	    pix.type = gValWSpc[i]->gType;
	    pix.v = gValWSpc[i]->gVal[0];
	    WlzValueConvertPixel(&pix, pix, WLZ_GREY_UBYTE);
	  }
	  b[i] = pix.v.ubv;
	}
	WLZ_RGBA_RGBA_SET(gval.rgbv, b[0], b[1], b[2], b[3]);
	*gwsp.u_grintptr.rgbp = gval.rgbv;
	gwsp.u_grintptr.rgbp++;
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
    if(iwsp.gryptr == &gwsp)
    {
      (void )WlzEndGreyScan(&gwsp);
    }
    for(i=0; i < numObjs; i++)
    {
      WlzGreyValueFreeWSp(gValWSpc[i]);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	New 3D domain object with corresponding WLZ_GREY_RGBA values.
* \ingroup      WlzValuesUtils
* \brief	Creates a WLZ_GREY_RGBA valued object from the given compound
* 		array. This is a static function which will always be called
* 		with valid parameters so they aren't checked.
* \param	cObj			Compound array object.
* \param	cSpc 			The colour space.
* \param	dstErr			Destination error pointer may be NULL.
*/
static WlzObject *WlzCompoundToRGBA3D(WlzCompoundArray *cObj,
  				WlzRGBAColorSpace cSpc, WlzErrorNum *dstErr)
{
  int		i;
  WlzIBox3 	bBox;
  WlzDomain	dom;
  WlzValues	val;
  WlzPixelV	bgd;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  
  dom.core = NULL;
  val.core = NULL;
  bgd.v.rgbv = 0;
  bgd.type = WLZ_GREY_RGBA;
  /* Check there are no tiled value tables as these aren't supported for
   * 3D yet. */
  for(i = 0; i < cObj->n; ++i)
  {
    WlzObject *obj;

    obj = cObj->o[i];
    if(obj && obj->values.core)
    {
      if(WlzGreyTableIsTiled(obj->values.core->type))
      {
	errNum = WLZ_ERR_VALUES_TYPE;
        break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox = WlzBoundingBox3I((WlzObject *)cObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN, bBox.zMin, bBox.zMax,
		               bBox.yMin, bBox.yMax, bBox.xMin, bBox.xMax,
			       &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
                                  bBox.zMin, bBox.zMax, bgd, NULL,
			          &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP,
    		nPl;

    nPl = bBox.zMax - bBox.zMin + 1;
    for(idP = 0; idP < nPl; ++idP)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	int	idC;
	WlzObject *rObj2 = NULL;
	WlzCompoundArray *cObj2 = NULL;
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

	cObj2 = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_2, 1, cObj->n,
	                             NULL, WLZ_2D_DOMAINOBJ, &errNum2);

	idC = 0;
	while((errNum2 == WLZ_ERR_NONE) && (idC < cObj->n))
	{
	  WlzDomain	dom2;
	  WlzValues	val2;

	  dom2.core = NULL;
	  val2.core = NULL;
	  if((cObj->o[idC] != NULL) &&
	     (cObj->o[idC]->type == WLZ_3D_DOMAINOBJ))
	  {
	    if((bBox.zMin + idP >= cObj->o[idC]->domain.p->plane1) &&
	       (bBox.zMin + idP <= cObj->o[idC]->domain.p->lastpl))
	    {
	      int	idP2;

	      idP2 = bBox.zMin + idP - cObj->o[idC]->domain.p->plane1;
	      dom2 = *(cObj->o[idC]->domain.p->domains + idP2);
	      val2 = *(cObj->o[idC]->values.vox->values + idP2);
	    }
	  }
	  cObj2->o[idC] = (dom2.core == NULL)?
	                  WlzMakeEmpty(&errNum2):
	                  WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2, val2,
	                              NULL, NULL, &errNum2);
	  ++idC;
	}
	if(errNum2 == WLZ_ERR_NONE)
	{
	  rObj2 = WlzCompoundToRGBA2D(cObj2, cSpc, &errNum2);
	}
	if(errNum2 == WLZ_ERR_NONE)
	{
	  dom.p->domains[idP] = WlzAssignDomain(rObj2->domain, NULL);
	  val.vox->values[idP] = WlzAssignValues(rObj2->values, NULL);
	}
	(void )WlzFreeObj(rObj2);
	(void )WlzFreeObj((WlzObject *)cObj2);
	if(errNum2 != WLZ_ERR_NONE)
	{
	  errNum = errNum2;
	}
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj != NULL)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      (void )WlzFreePlaneDomain(dom.p);
      (void )WlzFreeVoxelValueTb(val.vox);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
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
    newBck.v.shv = (short )WLZ_RGBA_MODULUS(oldBck.v.rgbv);

    /* make values table and return object */
    values.v = WlzNewValueTb(obj, type, newBck, &errNum);
    rtnObj = WlzMakeMain(obj->type, obj->domain, values,
			 NULL, NULL, &errNum);
  }

  /* iterate through objects setting values */
  if( errNum == WLZ_ERR_NONE ){
    WlzIntervalWSpace	iwsp0, iwsp1;
    WlzGreyWSpace	gwsp0, gwsp1;
    int			j, k;

    errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0);
    errNum = WlzInitGreyScan(rtnObj, &iwsp1, &gwsp1);
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
      errNum = WlzNextGreyInterval(&iwsp1);
      for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	    gwsp0.u_grintptr.rgbp++){
	*(gwsp1.u_grintptr.shp++) = (short )
	                            WLZ_RGBA_MODULUS(*(gwsp0.u_grintptr.rgbp));
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

static WlzObject *WlzRGBAToModulus3D(WlzObject *obj, WlzErrorNum *dstErr)
{
  int           pln;
  WlzPlaneDomain *gDom,
  		 *rDom = NULL;
  WlzVoxelValues *gVal,
                 *rVal = NULL;
  WlzObject	*rObj=NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  gDom = obj->domain.p;
  gVal = obj->values.vox;
  if(errNum == WLZ_ERR_NONE)
  {
    rDom = WlzMakePlaneDomain(gDom->type,
			      gDom->plane1, gDom->lastpl,
			      gDom->line1, gDom->lastln,
			      gDom->kol1, gDom->lastkl,
			      &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rVal = WlzMakeVoxelValueTb(gVal->type, gVal->plane1, gVal->lastpl,
    			       gVal->bckgrnd, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain dom;
    WlzValues val;

    dom.p = rDom;
    val.vox = rVal;
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, val, NULL, obj, &errNum);
  }
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(pln = gDom->plane1; pln <= gDom->lastpl; ++pln)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain *gDom2,
      		*rDom2;
      WlzValues *gVal2,
                *rVal2;
      WlzErrorNum errNum2 = WLZ_ERR_NONE;

      if(((gDom2 = gDom->domains + pln - gDom->plane1) != NULL) &&
	 ((*gDom2).core != NULL) &&
         ((rDom2 = rDom->domains + pln - rDom->plane1) != NULL) &&
         ((gVal2 = gVal->values  + pln - gDom->plane1) != NULL) &&
         ((rVal2 = rVal->values  + pln - rDom->plane1) != NULL))
      {
        WlzObject *gObj2 = NULL,
                  *rObj2 = NULL;

        if((gObj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *gDom2, *gVal2, NULL, NULL,
                                &errNum2)) != NULL)
        {
          rObj2 = WlzRGBAToModulus(gObj2, &errNum2);
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    *rDom2 = WlzAssignDomain(rObj2->domain, NULL);
	    *rVal2 = WlzAssignValues(rObj2->values, NULL);
	  }
        }
        (void )WlzFreeObj(gObj2);
        (void )WlzFreeObj(rObj2);
      }
#ifdef _OPENMP
#pragma omp critical
      {
#endif
        if(errNum2 != WLZ_ERR_NONE)
        {
          errNum = errNum2;
        }
#ifdef _OPENMP
      }
#endif
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj)
    {
      WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      (void )WlzFreePlaneDomain(rDom);
      (void )WlzFreeVoxelValueTb(rVal);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
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
  WlzIntervalWSpace	oldiwsp = {0},
  			newiwsp = {0};
  WlzGreyWSpace		oldgwsp,
  			newgwsp;
  WlzObjectType		newvtbltype;
  WlzPixelV		bg;
  WlzValues		values;
  int 			k, redVal, greenVal, blueVal, greyVal = 0;
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
      greyVal = (int )WLZ_CLAMP(bg.v.dbv, 0, 255);
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
	greyVal = (int )WLZ_CLAMP(*(go.dbp), 0, 255);
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
  if(oldiwsp.gryptr == &oldgwsp)
  {
    (void )WlzEndGreyScan(&oldgwsp);
  }
  if(newiwsp.gryptr == &newgwsp)
  {
    (void )WlzEndGreyScan(&newgwsp);
  }
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
  while( (errNum == WLZ_ERR_NONE) && nplanes-- ){
    if(((*domains).core == NULL) || ((*values).core == NULL)){
      (*ndomains).core = NULL;
      (*nvalues).core = NULL;
    }
    else if((temp = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, *values,
				NULL, NULL, &errNum)) != NULL){

      if( temp->domain.i != NULL ){
	if((obj1 = WlzIndexToRGBA(temp, colormap, &errNum)) != NULL){
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
  /* return a new object */
  if( errNum == WLZ_ERR_NONE ){
    domain.p = npdom;
    vals.vox = nvoxtab;
    if((obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			   NULL, obj, &errNum)) != NULL){
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
    newBck.v.ubv = (WlzUByte )WlzRGBAPixelValue(oldBck, chan, &errNum);
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
    int			j, k;

    if((errNum = WlzInitGreyScan(obj, &iwsp0, &gwsp0)) == WLZ_ERR_NONE)
    {
      if((errNum = WlzInitGreyScan(rtnObj, &iwsp1, &gwsp1)) == WLZ_ERR_NONE)
      {
	pixVal.type = WLZ_GREY_RGBA;
	while((errNum == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iwsp0)) == WLZ_ERR_NONE)){
	  errNum = WlzNextGreyInterval(&iwsp1);
	  for(j=0, k=iwsp0.lftpos; k <= iwsp0.rgtpos; j++, k++,
	      gwsp0.u_grintptr.rgbp++){
	    pixVal.v.rgbv = (*(gwsp0.u_grintptr.rgbp));
	    *(gwsp1.u_grintptr.ubp++) = (WlzUByte )
	      WlzRGBAPixelValue(pixVal, chan, &errNum);
	  }
	}
	(void )WlzEndGreyScan(&gwsp1);
      }
      (void )WlzEndGreyScan(&gwsp0);
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
  WlzErrorNum	errNum=WLZ_ERR_UNIMPLEMENTED;

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
