#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzShadeCorrect.c
* \author       Bill Hill
* \date         January 2001
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	
* \ingroup	WlzValueFilters
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

static WlzObject 		*WlzShadeCorrect2DG(
				  WlzObject *srcObj,
				  WlzObject *shdObj,
				  double nrmVal,
				  int inPlace,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzShadeCorrect2D(
				  WlzObject *srcObj,
				  WlzObject *shdObj,
				  double nrmVal,
				  int inPlace,
				  WlzErrorNum *dstErr);

/*!
* \return	Shade corrected object or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Shade corrects the given domain object with grey
*               values.
*		\f[
		  p_i = \frac{n  o_i}{s_i}
		\f]
*               The shade corrected image P with values \f$p_i\f$ is
*		created by applying a correction factor to each image
*		value of the given image O with values \f$o_i\f$.
* \param	srcObj			Given object to be shade corrected.
* \param	shdObj			Given bright field object.
* \param	nrmVal			Normalization value.
* \param	inPlace			Modify the grey values of the given
*					object in place if non-zero.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzObject	*WlzShadeCorrect(WlzObject *srcObj, WlzObject *shdObj,
			         double nrmVal, int inPlace,
			         WlzErrorNum *dstErr)
{
  return WlzShadeCorrectBFDF(srcObj, shdObj, NULL, nrmVal, inPlace, dstErr);
}

/*!
* \return	Shade corrected object or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Shade corrects the given domain object with grey
*               values.
*		\f[
		  p_i = n \frac{(o_i - d_i)}{(s_i - d_i}
		\f]
*               The shade corrected image P with values \f$p_i\f$ is
*		created by applying a correction factor to each image
*		value of the given image O with values \f$o_i\f$.
\f$s_i\f$ and \f$d_i\f$ are the bright and dark-field shade image values
respectively.
* \param	srcObj			Given object to be shade corrected.
* \param	shdObj			Given bright field object.
* \param	shdDFObj		Given dark field object (may be NULL).
* \param	nrmVal			Normalization value.
* \param	inPlace			Modify the grey values of the given
*					object in place if non-zero.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzObject	*WlzShadeCorrectBFDF(
  WlzObject 	*srcObj,
  WlzObject 	*shdObj,
  WlzObject 	*shdDFObj,
  double 	nrmVal,
  int 		inPlace,
  WlzErrorNum 	*dstErr)
{
  WlzObject	*rtnObj = NULL;
  WlzObject	*obj1, *obj2;
  WlzCompoundArray	*inCobj, *outCobj;
  WlzObject	**outArray;
  WlzValues	values;
  int		i;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcObj == NULL) || (shdObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((srcObj->domain.core == NULL) || (shdObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((srcObj->values.core == NULL) || (shdObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
    case WLZ_2D_DOMAINOBJ:
      if(shdObj->type != srcObj->type)
      {
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if( shdDFObj ){
	if( WlzGreyTypeFromObj(srcObj, &errNum) == WLZ_GREY_RGBA ){
	  obj1 = WlzRGBAImageArithmetic(srcObj, shdDFObj, WLZ_BO_SUBTRACT,
					0, &errNum);
	  obj2 = WlzRGBAImageArithmetic(shdObj, shdDFObj, WLZ_BO_SUBTRACT,
					0, &errNum);
	}
	else {
	  obj1 = WlzImageArithmetic(srcObj, shdDFObj, WLZ_BO_SUBTRACT,
				    0, &errNum);
	  obj2 = WlzImageArithmetic(shdObj, shdDFObj, WLZ_BO_SUBTRACT,
				    0, &errNum);
	}
	rtnObj = WlzShadeCorrect2D(obj1, obj2, nrmVal, inPlace, &errNum);
	WlzFreeObj(obj1);
	WlzFreeObj(obj2);;
      }
      else {
	rtnObj = WlzShadeCorrect2D(srcObj, shdObj, nrmVal, inPlace, &errNum);
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      inCobj = (WlzCompoundArray *) srcObj;
      if( inCobj->n > 0 ){
	if( (outArray = (WlzObject **)
	     AlcCalloc(inCobj->n, sizeof(WlzObject *))) == NULL){
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else {
	  for(i=0; (i < inCobj->n) && (errNum == WLZ_ERR_NONE); i++){
	    /* allow special behaviour of a non-compound shade image
	       to be used for each object in a compound but must align */
	    switch( shdObj->type ){
	    case WLZ_2D_DOMAINOBJ:
	    case WLZ_3D_DOMAINOBJ:
	      outArray[i] = WlzShadeCorrectBFDF(inCobj->o[i],
						shdObj, shdDFObj,
						nrmVal, inPlace, &errNum);
	      break;

	    case WLZ_COMPOUND_ARR_1:
	    case WLZ_COMPOUND_ARR_2:
	      outArray[i] =
		WlzShadeCorrectBFDF(inCobj->o[i],
				    ((WlzCompoundArray *) shdObj)->o[i],
				    shdDFObj?
				    ((WlzCompoundArray *) shdDFObj)->o[i]:NULL,
				    nrmVal, inPlace, &errNum);
	      break;
	    }
	      
	  }
	  if( errNum == WLZ_ERR_NONE ){
	    outCobj = WlzMakeCompoundArray(srcObj->type, 3, inCobj->n,
					   outArray, outArray[0]->type,
					   &errNum);
	    rtnObj = (WlzObject *) outCobj;
	  }
	  else {
	    while( i >= 0 ){
	      WlzFreeObj(outArray[i--]);
	    }
	  }
	  AlcFree(outArray);
	}
      }
      else {
	errNum = WLZ_ERR_OBJECT_DATA;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;

    case WLZ_TRANS_OBJ:
      if(rtnObj =
	 WlzShadeCorrectBFDF(srcObj->values.obj, shdObj, shdDFObj,
			     nrmVal, inPlace, &errNum) ){
	values.obj = rtnObj;
	return WlzMakeMain(WLZ_TRANS_OBJ, srcObj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	Shade corrected object or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Shade corrects the given 2D domain object with grey
*               values. This function just checks that the given and
*               shade objects have the same grey types, it then calls
*               WlzShadeCorrect2DG() to shade correct the given object.
* \param	srcObj			Given object to be shade corrected.
* \param	shdObj			Given bright field object.
* \param	nrmVal			Normalization value.
* \param	inPlace			Modify the grey values of the
*                                       given object if non-zero.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzShadeCorrect2D(WlzObject *srcObj, WlzObject *shdObj,
				    double nrmVal, int inPlace,
				    WlzErrorNum *dstErr)
{
  WlzGreyType	srcG,
  		shdG;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  srcG = WlzGreyTableTypeToGreyType(srcObj->values.v->type, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    shdG = WlzGreyTableTypeToGreyType(shdObj->values.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(srcG != shdG)
    {
      errNum = WLZ_ERR_GREY_TYPE;
    }
    else
    {
      rtnObj = WlzShadeCorrect2DG(srcObj, shdObj, nrmVal, inPlace, &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	Shade corrected object or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Shade corrects the given 2D domain object with grey
*               values. Grey value types known to be the same.
* \param	srcObj			Given object to be shade
*                                       corrected.
* \param	shdObj			Given bright field object.
* \param	nrmVal			Normalization value.
* \param	inPlace			Modify the grey values of the
*                                       given object if non-zero.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzShadeCorrect2DG(WlzObject *srcObj, WlzObject *shdObj,
				     double nrmVal, int inPlace,
				     WlzErrorNum *dstErr)
{
  int		tI0,
  		iCnt, red, green, blue;
  double	tD0;
  UINT		tUI0, tUI1;
  WlzObject	*uObj = NULL,
		*uSrcObj = NULL,
		*uShdObj = NULL,
  		*rtnObj = NULL;
  WlzGreyP	srcPix,
  		shdPix,
		rtnPix;
  WlzValues	newVal;
  WlzIntervalWSpace srcIWSp,
  		shdIWSp,
		rtnIWSp;
  WlzGreyWSpace	srcGWSp,
  		shdGWSp,
		rtnGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find intersection of the given and shade objects. */
  uObj = WlzIntersect2(srcObj, shdObj, &errNum);
  /* Make new objects with the values of the given and shade objects
   * but the domain of their intersection. */
  if(errNum == WLZ_ERR_NONE)
  {
    uSrcObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, uObj->domain, srcObj->values,
    			  NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    uShdObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, uObj->domain, shdObj->values,
    			  NULL, NULL, &errNum);
  }
  /* Make a new object, again using the union for the domain, but this time
   * either sharing the given objects values or creating a new value table. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(inPlace)
    {
      rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, uObj->domain, srcObj->values,
      			   NULL, NULL, &errNum);
    }
    else
    {
      newVal.v = WlzNewValueTb(uObj, srcObj->values.core->type,
      			       WlzGetBackground(srcObj, NULL), &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        if((rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, uObj->domain, newVal,
			         NULL, NULL, &errNum)) == NULL)
        {
	  (void )WlzFreeValueTb(newVal.v);
	}
      }
    }
  }

  /* Work through the intervals setting the grey values. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((errNum = WlzInitGreyScan(uSrcObj, &srcIWSp,
    				  &srcGWSp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(uShdObj, &shdIWSp,
    				  &shdGWSp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(rtnObj, &rtnIWSp,
    				  &rtnGWSp)) == WLZ_ERR_NONE))
    {
      while(((errNum = WlzNextGreyInterval(&srcIWSp)) == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&shdIWSp)) == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&rtnIWSp)) == WLZ_ERR_NONE))
      {
	srcPix = srcGWSp.u_grintptr;
	shdPix = shdGWSp.u_grintptr;
	rtnPix = rtnGWSp.u_grintptr;
	iCnt = rtnIWSp.rgtpos - rtnIWSp.lftpos + 1;
        switch(rtnGWSp.pixeltype)
	{
	  case WLZ_GREY_INT:
	    while(iCnt-- > 0)
	    {
	      tD0 = (*(srcPix.inp)++ * nrmVal) / (*(shdPix.inp)++ + 1.0);
	      *(rtnPix.inp)++ = WLZ_NINT(tD0);
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    while(iCnt-- > 0)
	    {
	      tD0 = (*(srcPix.shp)++ * nrmVal) / (*(shdPix.shp)++ + 1.0);
	      tI0 = WLZ_NINT(tD0);
	      *(rtnPix.shp)++ = WLZ_CLAMP(tI0, SHRT_MIN, SHRT_MAX);
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    while(iCnt-- > 0)
	    {
	      tD0 = (*(srcPix.ubp)++ * nrmVal) / (*(shdPix.ubp)++ + 1.0);
	      tI0 = WLZ_NINT(tD0);
	      *(rtnPix.ubp)++ = WLZ_CLAMP(tI0, 0, 255);
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    while(iCnt-- > 0)
	    {
	      tD0 = (*(srcPix.flp)++ * nrmVal) / (*(shdPix.flp)++ + 1.0);
	      *(rtnPix.flp)++ = tD0;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    while(iCnt-- > 0)
	    {
	      tD0 = (*(srcPix.dbp)++ * nrmVal) / (*(shdPix.dbp)++ + 1.0);
	      *(rtnPix.dbp)++ = tD0;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    while(iCnt-- > 0)
	    {
	      /* slightly different logic here. Avoid divide by zero
		 by explicit check */
	      tUI0 = *(srcPix.rgbp)++;
	      tUI1 = *(shdPix.rgbp)++;
	      red = WLZ_RGBA_RED_GET(tUI1);
	      red = (red)?((WLZ_RGBA_RED_GET(tUI0) * nrmVal)) / red : nrmVal;
	      red = WLZ_CLAMP(red, 0, 255);
	      green = WLZ_RGBA_GREEN_GET(tUI1);
	      green = (green)?((WLZ_RGBA_GREEN_GET(tUI0) * nrmVal)) / green : nrmVal;
	      green = WLZ_CLAMP(green, 0, 255);
	      blue = WLZ_RGBA_BLUE_GET(tUI1);
	      blue = (blue)?((WLZ_RGBA_BLUE_GET(tUI0) * nrmVal)) / blue : nrmVal;
	      blue = WLZ_CLAMP(blue, 0, 255);
	      WLZ_RGBA_RGBA_SET(tUI0, red, green, blue, 255);
	      *(rtnPix.rgbp)++ = tUI0;
	    }
	    break;
	}
      }
      if(errNum == WLZ_ERR_EOO)         /* Reset error from end of intervals */
      {
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  (void )WlzFreeObj(uObj);
  (void )WlzFreeObj(uSrcObj);
  (void )WlzFreeObj(uShdObj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rtnObj);
    rtnObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}
