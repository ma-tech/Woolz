#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzStringTypes.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for converting between Woolz data types
*		and a string representation.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 14-03-01 bill Add WLZ_ERR_UNIMPLEMENTED.
* 15-08-00 bill	Remove obsolete types WLZ_VECTOR_(INT)|(FLOAT) and
*		WLZ_POINT_(INT)|(FLOAT). Add WLZ_CONTOUR and the geometric
*		model types.
* 28-01-00 bill Add ALG_ERR_CONVERGENCE.
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzStringFromObjType				
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given object's type.	
* Global refs:	-						
* Parameters:	WlzObject *obj:		Given object.		
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromObjType(WlzObject *obj,
				      WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*oTypeStr = NULL;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        oTypeStr = "WLZ_2D_DOMAINOBJ";
	break;
      case WLZ_3D_DOMAINOBJ:
        oTypeStr = "WLZ_3D_DOMAINOBJ";
	break;
      case WLZ_TRANS_OBJ:
        oTypeStr = "WLZ_TRANS_OBJ";
	break;
      case WLZ_3D_WARP_TRANS:
        oTypeStr = "WLZ_3D_WARP_TRANS";
	break;
      case WLZ_2D_POLYGON:
        oTypeStr = "WLZ_2D_POLYGON";
	break;
      case WLZ_BOUNDLIST:
        oTypeStr = "WLZ_BOUNDLIST";
	break;
      case WLZ_CONV_HULL:
        oTypeStr = "WLZ_CONV_HULL";
	break;
      case WLZ_HISTOGRAM:
        oTypeStr = "WLZ_HISTOGRAM";
	break;
      case WLZ_3D_POLYGON:
        oTypeStr = "WLZ_3D_POLYGON";
	break;
      case WLZ_CONTOUR:
        oTypeStr = "WLZ_CONTOUR";
	break;
      case WLZ_RECTANGLE:
        oTypeStr = "WLZ_RECTANGLE";
	break;
      case WLZ_CONVOLVE_INT:
        oTypeStr = "WLZ_CONVOLVE_INT";
	break;
      case WLZ_CONVOLVE_FLOAT:
        oTypeStr = "WLZ_CONVOLVE_FLOAT";
	break;
      case WLZ_AFFINE_TRANS:
        oTypeStr = "WLZ_AFFINE_TRANS";
	break;
      case WLZ_WARP_TRANS:
        oTypeStr = "WLZ_WARP_TRANS";
	break;
      case WLZ_FMATCHOBJ:
        oTypeStr = "WLZ_FMATCHOBJ";
	break;
      case WLZ_TEXT:
        oTypeStr = "WLZ_TEXT";
	break;
      case WLZ_COMPOUND_ARR_1:
        oTypeStr = "WLZ_COMPOUND_ARR_1";
	break;
      case WLZ_COMPOUND_ARR_2:
        oTypeStr = "WLZ_COMPOUND_ARR_2";
	break;
      case WLZ_COMPOUND_LIST_1:
        oTypeStr = "WLZ_COMPOUND_LIST_1";
	break;
      case WLZ_COMPOUND_LIST_2:
        oTypeStr = "WLZ_COMPOUND_LIST_2";
	break;
      case WLZ_PROPERTY_OBJ:
        oTypeStr = "WLZ_PROPERTY_OBJ";
	break;
      case WLZ_EMPTY_OBJ:
        oTypeStr = "WLZ_EMPTY_OBJ";
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
  return(oTypeStr);
}

/************************************************************************
* Function:	WlzStringToObjType				
* Returns:	WlzObjectType:		Woolz object type or WLZ_NULL
*					on error.		
* Purpose:	Finds an enumerated type for the given object type
*		string.						
* Global refs:	-						
* Parameters:	const char *oTypeStr:	Given object type string.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzObjectType	WlzStringToObjType(const char *oTypeStr,
				   WlzErrorNum *dstErr)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_OBJECT_TYPE;
  WlzObjectType	oType = WLZ_NULL;

  if(WlzStringMatchValue(&tI0, oTypeStr,
		"WLZ_2D_DOMAINOBJ", WLZ_2D_DOMAINOBJ,
		"WLZ_3D_DOMAINOBJ", WLZ_3D_DOMAINOBJ,
		"WLZ_TRANS_OBJ", WLZ_TRANS_OBJ,
		"WLZ_3D_WARP_TRANS", WLZ_3D_WARP_TRANS,
		"WLZ_2D_POLYGON", WLZ_2D_POLYGON,
		"WLZ_BOUNDLIST", WLZ_BOUNDLIST,
		"WLZ_CONV_HULL", WLZ_CONV_HULL,
		"WLZ_CONTOUR", WLZ_CONTOUR,
		"WLZ_HISTOGRAM", WLZ_HISTOGRAM,
		"WLZ_3D_POLYGON", WLZ_3D_POLYGON,
		"WLZ_RECTANGLE", WLZ_RECTANGLE,
		"WLZ_CONVOLVE_INT", WLZ_CONVOLVE_INT,
		"WLZ_CONVOLVE_FLOAT", WLZ_CONVOLVE_FLOAT,
		"WLZ_AFFINE_TRANS", WLZ_AFFINE_TRANS,
		"WLZ_WARP_TRANS", WLZ_WARP_TRANS,
		"WLZ_FMATCHOBJ", WLZ_FMATCHOBJ,
		"WLZ_TEXT", WLZ_TEXT,
		"WLZ_COMPOUND_ARR_1", WLZ_COMPOUND_ARR_1,
		"WLZ_COMPOUND_ARR_2", WLZ_COMPOUND_ARR_2,
		"WLZ_COMPOUND_LIST_1", WLZ_COMPOUND_LIST_1,
		"WLZ_COMPOUND_LIST_2", WLZ_COMPOUND_LIST_2,
		"WLZ_PROPERTY_OBJ", WLZ_PROPERTY_OBJ,
		"WLZ_EMPTY_OBJ", WLZ_EMPTY_OBJ,
		NULL))
  {
    oType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oType);
}

/************************************************************************
* Function:	WlzStringFromObjDomainType			
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given object's domain type.
* Global refs:	-						
* Parameters:	WlzObject *obj:		Given object.		
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromObjDomainType(WlzObject *obj,
					    WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*oDomTypeStr = NULL;

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
      case WLZ_2D_DOMAINOBJ:
	switch(obj->domain.core->type)
	{
          case WLZ_INTERVALDOMAIN_INTVL:
	    oDomTypeStr = "WLZ_INTERVALDOMAIN_INTVL";
	    break;
          case WLZ_INTERVALDOMAIN_RECT:
	    oDomTypeStr = "WLZ_INTERVALDOMAIN_RECT";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
        }
	break;
      case WLZ_3D_DOMAINOBJ:
	switch(obj->domain.core->type)
	{
          case WLZ_PLANEDOMAIN_DOMAIN:
            oDomTypeStr = "WLZ_PLANEDOMAIN_DOMAIN";
	    break;
          case WLZ_PLANEDOMAIN_POLYGON:
            oDomTypeStr = "WLZ_PLANEDOMAIN_POLYGON";
	    break;
          case WLZ_PLANEDOMAIN_BOUNDLIST:
            oDomTypeStr = "WLZ_PLANEDOMAIN_BOUNDLIST";
	    break;
          case WLZ_PLANEDOMAIN_HISTOGRAM:
            oDomTypeStr = "WLZ_PLANEDOMAIN_HISTOGRAM";
	    break;
          case WLZ_PLANEDOMAIN_AFFINE:
            oDomTypeStr = "WLZ_PLANEDOMAIN_AFFINE";
	    break;
          case WLZ_PLANEDOMAIN_WARP:
            oDomTypeStr = "WLZ_PLANEDOMAIN_WARP";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
        }
	break;
      case WLZ_TRANS_OBJ:
	switch(obj->domain.core->type)
	{
	  case WLZ_TRANSFORM_2D_AFFINE:
	    oDomTypeStr = "WLZ_TRANSFORM_2D_AFFINE";
	    break;
	  case WLZ_TRANSFORM_3D_AFFINE:
	    oDomTypeStr = "WLZ_TRANSFORM_3D_AFFINE";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_3D_WARP_TRANS:
	switch(obj->domain.core->type)
	{
	  case WLZ_WARP_TRANS:
	    oDomTypeStr = "WLZ_WARP_TRANS";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_2D_POLYGON:
	switch(obj->domain.core->type)
	{
	  case WLZ_POLYGON_INT:
	    oDomTypeStr = "WLZ_POLYGON_INT";
	    break;
	  case WLZ_POLYGON_FLOAT:
	    oDomTypeStr = "WLZ_POLYGON_FLOAT";
	    break;
	  case WLZ_POLYGON_DOUBLE:
	    oDomTypeStr = "WLZ_POLYGON_DOUBLE";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_BOUNDLIST:
	switch(obj->domain.core->type)
	{
	  case WLZ_BOUNDLIST_PIECE:
	    oDomTypeStr = "WLZ_BOUNDLIST_PIECE";
	    break;
	  case WLZ_BOUNDLIST_HOLE:
	    oDomTypeStr = "WLZ_BOUNDLIST_HOLE";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_CONTOUR:
        oDomTypeStr = "WLZ_CONTOUR";
	break;
      case WLZ_HISTOGRAM:
	switch(obj->domain.core->type)
	{
	  case WLZ_HISTOGRAMDOMAIN_INT:
	    oDomTypeStr = "WLZ_HISTOGRAMDOMAIN_INT";
	    break;
	  case WLZ_HISTOGRAMDOMAIN_FLOAT:
	    oDomTypeStr = "WLZ_HISTOGRAMDOMAIN_FLOAT";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
        break;
      case WLZ_RECTANGLE:
	switch(obj->domain.core->type)
	{
	  case WLZ_RECTANGLE_DOMAIN_INT:
	    oDomTypeStr = "WLZ_RECTANGLE_DOMAIN_INT";
	    break;
	  case WLZ_RECTANGLE_DOMAIN_FLOAT:
	    oDomTypeStr = "WLZ_RECTANGLE_DOMAIN_FLOAT";
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
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
  return(oDomTypeStr);
}

/************************************************************************
* Function:	WlzStringToObjDomainType			
* Returns:	WlzObjectType:		Woolz object type or WLZ_NULL
*					on error.		
* Purpose:	Finds an enumerated type for the given object domain
*		type string.					
* Global refs:	-						
* Parameters:	const char *oDomTypeStr:Given object domain type
*					string.			
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzObjectType	WlzStringToObjDomainType(const char *oDomTypeStr,
				         WlzErrorNum *dstErr)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_OBJECT_TYPE;
  WlzObjectType	oDomType = WLZ_NULL;

  if(WlzStringMatchValue(&tI0, oDomTypeStr,
		"WLZ_INTERVALDOMAIN_INTVL", WLZ_INTERVALDOMAIN_INTVL,
		"WLZ_INTERVALDOMAIN_RECT", WLZ_INTERVALDOMAIN_RECT,
		"WLZ_PLANEDOMAIN_DOMAIN", WLZ_PLANEDOMAIN_DOMAIN,
		"WLZ_PLANEDOMAIN_POLYGON", WLZ_PLANEDOMAIN_POLYGON,
		"WLZ_PLANEDOMAIN_BOUNDLIST", WLZ_PLANEDOMAIN_BOUNDLIST,
		"WLZ_PLANEDOMAIN_HISTOGRAM", WLZ_PLANEDOMAIN_HISTOGRAM,
		"WLZ_PLANEDOMAIN_AFFINE", WLZ_PLANEDOMAIN_AFFINE,
		"WLZ_PLANEDOMAIN_WARP", WLZ_PLANEDOMAIN_WARP,
		"WLZ_TRANSFORM_2D_AFFINE", WLZ_TRANSFORM_2D_AFFINE,
		"WLZ_TRANSFORM_3D_AFFINE", WLZ_TRANSFORM_3D_AFFINE,
		"WLZ_WARP_TRANS", WLZ_WARP_TRANS,
		"WLZ_POLYGON_INT", WLZ_POLYGON_INT,
		"WLZ_POLYGON_FLOAT", WLZ_POLYGON_FLOAT,
		"WLZ_POLYGON_DOUBLE", WLZ_POLYGON_DOUBLE,
		"WLZ_BOUNDLIST_PIECE", WLZ_BOUNDLIST_PIECE,
		"WLZ_BOUNDLIST_HOLE", WLZ_BOUNDLIST_HOLE,
		"WLZ_CONTOUR", WLZ_CONTOUR,
		"WLZ_HISTOGRAMDOMAIN_INT", WLZ_HISTOGRAMDOMAIN_INT,
		"WLZ_HISTOGRAMDOMAIN_FLOAT", WLZ_HISTOGRAMDOMAIN_FLOAT,
		"WLZ_RECTANGLE_DOMAIN_INT", WLZ_RECTANGLE_DOMAIN_INT,
		"WLZ_RECTANGLE_DOMAIN_FLOAT", WLZ_RECTANGLE_DOMAIN_FLOAT,
		NULL))
  {
    oDomType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oDomType);
}

/************************************************************************
* Function:	WlzStringFromObjValuesType			
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given object's values type.
* Global refs:	-						
* Parameters:	WlzObject *obj:		Given object.		
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromObjValuesType(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzObjectType	valType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*oValTypeStr = NULL;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    valType = WlzGreyTableTypeToTableType(obj->values.core->type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(obj->type)
      {
	case WLZ_2D_DOMAINOBJ:
	  switch(valType)
	  {
	    case WLZ_GREY_TAB_RAGR:
	      oValTypeStr = "WLZ_GREY_TAB_RAGR";
	      break;
	    case WLZ_GREY_TAB_RECT:
	      oValTypeStr = "WLZ_GREY_TAB_RECT";
	      break;
	    case WLZ_GREY_TAB_INTL:
	      oValTypeStr = "WLZ_GREY_TAB_INTL";
	      break;
	    default:
	      errNum = WLZ_ERR_VALUES_TYPE;
	      break;
	  }
	  break;
	case WLZ_3D_DOMAINOBJ:
	  switch(obj->values.core->type)
	  {
	    case WLZ_VOXELVALUETABLE_GREY:
	      oValTypeStr = "WLZ_VOXELVALUETABLE_GREY";
	      break;
	    default:
	      errNum = WLZ_ERR_VOXELVALUES_TYPE;
	      break;
	  }
	  break;
	case WLZ_TRANS_OBJ:
	  oValTypeStr = WlzStringFromObjType(obj->values.obj, &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oValTypeStr);
}

/************************************************************************
* Function:	WlzStringToObjValuesType			
* Returns:	WlzObjectType:		Woolz object type or WLZ_NULL
*					on error.		
* Purpose:	Finds an enumerated type for the given object values
*		type string.					
* Global refs:	-						
* Parameters:	const char *oValTypeStr:Given object values type
*					string.			
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzObjectType	WlzStringToObjValuesType(const char *oValTypeStr,
				         WlzErrorNum *dstErr)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_OBJECT_TYPE;
  WlzObjectType	oValType = WLZ_NULL;

  if(WlzStringMatchValue(&tI0, oValTypeStr,
		"WLZ_GREY_TAB_RAGR", WLZ_GREY_TAB_RAGR,
		"WLZ_GREY_TAB_RECT", WLZ_GREY_TAB_RECT,
		"WLZ_GREY_TAB_INTL", WLZ_GREY_TAB_INTL,
		"WLZ_VOXELVALUETABLE_GREY", WLZ_VOXELVALUETABLE_GREY,
		NULL))
  {
    oValType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(oValType);
}

/************************************************************************
* Function:	WlzStringFromTransformType			
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given transform type.	
* Global refs:	-						
* Parameters:	WlzTransformType tType:	Given transform type.	
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromTransformType(WlzTransformType tType,
				            WlzErrorNum *dstErr)
{
  const char	*tStr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(tType)
  {
    case WLZ_TRANSFORM_2D_AFFINE:
      tStr = "WLZ_TRANSFORM_2D_AFFINE";
      break;
    case WLZ_TRANSFORM_3D_AFFINE:
      tStr = "WLZ_TRANSFORM_3D_AFFINE";
      break;
    default:
      errNum = WLZ_ERR_TRANSFORM_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tStr);
}

/************************************************************************
* Function:	WlzStringToTransformType			
* Returns:	WlzTransformType: 	Woolz transform type.	
* Purpose:	Finds an enumerated type for the given transform
*		type string.					
* Global refs:	-						
* Parameters:	const char *tStr:	Given transform type string.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzTransformType WlzStringToTransformType(const char *tStr,
				          WlzErrorNum *dstErr)
{
  int		tI0;
  WlzTransformType	tType = WLZ_TRANSFORM_2D_AFFINE;
  WlzErrorNum	errNum = WLZ_ERR_TRANSFORM_TYPE;

  if(WlzStringMatchValue(&tI0, tStr,
  		         "WLZ_TRANSFORM_2D_AFFINE", WLZ_TRANSFORM_2D_AFFINE,
			 "WLZ_TRANSFORM_3D_AFFINE", WLZ_TRANSFORM_3D_AFFINE,
			 NULL))
  {
    tType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tType);
}

/************************************************************************
* Function:	WlzStringFromGMModelType			
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given transform type.	
* Global refs:	-						
* Parameters:	WlzGMModelType mType:	Given model type.	
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromGMModelType(WlzGMModelType mType,
				          WlzErrorNum *dstErr)
{
  const char	*tStr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(mType)
  {
    case WLZ_GMMOD_2I:
      tStr = "WLZ_GMMOD_2I";
      break;
    case WLZ_GMMOD_2D:
      tStr = "WLZ_GMMOD_2D";
      break;
    case WLZ_GMMOD_3I:
      tStr = "WLZ_GMMOD_3I";
      break;
    case WLZ_GMMOD_3D:
      tStr = "WLZ_GMMOD_3D";
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tStr);
}

/************************************************************************
* Function:	WlzStringToGMModelType			
* Returns:	WlzGMModelType: 	Woolz GM model type.	
* Purpose:	Finds an enumerated type for the given GM model
*		type string.					
* Global refs:	-						
* Parameters:	const char *tStr:	Given GM Model type string.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzGMModelType	WlzStringToGMModelType(const char *tStr,
					WlzErrorNum *dstErr)
{
  int		tI0;
  WlzGMModelType mType = WLZ_GMMOD_2I;
  WlzErrorNum	errNum = WLZ_ERR_DOMAIN_TYPE;

  if(WlzStringMatchValue(&tI0, tStr,
			 "WLZ_GMMOD_2I", WLZ_GMMOD_2I,
			 "WLZ_GMMOD_2D", WLZ_GMMOD_2D,
			 "WLZ_GMMOD_3I", WLZ_GMMOD_3I,
			 "WLZ_GMMOD_3D", WLZ_GMMOD_3D,
			 NULL))
  {
    mType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mType);
}

/************************************************************************
* Function:	WlzStringFromGreyType				
* Returns:	const char*:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given grey type.		
* Global refs:	-						
* Parameters:	WlzGreyType gType:	Given grey type.	
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
const char	*WlzStringFromGreyType(WlzGreyType gType,
				       WlzErrorNum *dstErr)
{
  const char	*gStr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_INT:
      gStr = "WLZ_GREY_INT";
      break;
    case WLZ_GREY_SHORT:
      gStr = "WLZ_GREY_SHORT";
      break;
    case WLZ_GREY_UBYTE:
      gStr = "WLZ_GREY_UBYTE";
      break;
    case WLZ_GREY_FLOAT:
      gStr = "WLZ_GREY_FLOAT";
      break;
    case WLZ_GREY_DOUBLE:
      gStr = "WLZ_GREY_DOUBLE";
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gStr);
}

/************************************************************************
* Function:	WlzStringToGreyType				
* Returns:	WlzGreyType:		Woolz grey type or	
*					WLZ_GREY_ERROR on error.
* Purpose:	Finds an enumerated type for the given grey type
*		string.						
* Global refs:	-						
* Parameters:	const char *gStr:	Given grey type string.	
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.		
************************************************************************/
WlzGreyType	WlzStringToGreyType(const char *gStr,
				    WlzErrorNum *dstErr)
{
  int		tI0;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_GREY_TYPE;

  if(WlzStringMatchValue(&tI0, gStr,
  		         "WLZ_GREY_INT", WLZ_GREY_INT,
			 "WLZ_GREY_SHORT", WLZ_GREY_SHORT,
			 "WLZ_GREY_UBYTE", WLZ_GREY_UBYTE,
			 "WLZ_GREY_FLOAT", WLZ_GREY_FLOAT,
			 "WLZ_GREY_DOUBLE", WLZ_GREY_DOUBLE,
			 NULL))
  {
    gType = tI0;
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gType);
}

/************************************************************************
* Function:	WlzStringFromErrorNum				
* Returns:	const char *:		Pointer to read only string or
*					NULL on error.		
* Purpose:	Finds a string for the given error.		
* Global refs:	-						
* Parameters:	WlzErrorNum wlzErr:	Given error.		
*		const char **dstMsgStr:	Destination pointer for a
*					'meaningful' message string,
*					may be NULL.		
************************************************************************/
const char	*WlzStringFromErrorNum(WlzErrorNum wlzErr,
				       const char **dstMsgStr)
{
  const char	*errStr = NULL,
    *msgStr = NULL;

  switch(wlzErr)
  {
  case WLZ_ERR_NONE:
    errStr = "WLZ_ERR_NONE",
      msgStr = "No error";
    break;
  case WLZ_ERR_EOO:
    errStr = "WLZ_ERR_EOO",
      msgStr = "End of object";
    break;
  case WLZ_ERR_OBJECT_NULL:
    errStr = "WLZ_ERR_OBJECT_NULL",
      msgStr = "NULL object";
    break;
  case WLZ_ERR_OBJECT_TYPE:
    errStr = "WLZ_ERR_OBJECT_TYPE",
      msgStr = "Invalid object type";
    break;
  case WLZ_ERR_OBJECT_DATA:
    errStr = "WLZ_ERR_OBJECT_DATA",
      msgStr = "Object data";
    break;
  case WLZ_ERR_DOMAIN_NULL :
    errStr = "WLZ_ERR_DOMAIN_NULL ",
      msgStr = "NULL domain";
    break;
  case WLZ_ERR_DOMAIN_TYPE:
    errStr = "WLZ_ERR_DOMAIN_TYPE",
      msgStr = "Invalid domain type";
    break;
  case WLZ_ERR_DOMAIN_DATA:
    errStr = "WLZ_ERR_DOMAIN_DATA",
      msgStr = "Domain data";
    break;
  case WLZ_ERR_VALUES_NULL:
    errStr = "WLZ_ERR_VALUES_NULL",
      msgStr = "NULL values";
    break;
  case WLZ_ERR_VALUES_TYPE:
    errStr = "WLZ_ERR_VALUES_TYPE",
      msgStr = "Invalid values type";
    break;
  case WLZ_ERR_VALUES_DATA:
    errStr = "WLZ_ERR_VALUES_DATA",
      msgStr = "Values data";
    break;
  case WLZ_ERR_PROPERTY_NULL:
    errStr = "WLZ_ERR_PROPERTY_NULL",
      msgStr = "NULL property";
    break;
  case WLZ_ERR_GMELM_NULL:
    errStr = "WLZ_ERR_GMELM_NULL",
      msgStr = "Geometric model element NULL";
    break;
  case WLZ_ERR_GMELM_TYPE:
    errStr = "WLZ_ERR_GMELM_TYPE",
      msgStr = "Geometric model element type";
    break;
  case WLZ_ERR_GMELM_DATA:
    errStr = "WLZ_ERR_GMELM_DATA",
      msgStr = "Geometric model element data";
    break;
  case WLZ_ERR_PARAM_NULL:
    errStr = "WLZ_ERR_PARAM_NULL",
      msgStr = "NULL parameter";
    break;
  case WLZ_ERR_PARAM_TYPE:
    errStr = "WLZ_ERR_PARAM_TYPE",
      msgStr = "Invalid parameter type";
    break;
  case WLZ_ERR_PARAM_DATA:
    errStr = "WLZ_ERR_PARAM_DATA",
      msgStr = "Parameter data";
    break;
  case WLZ_ERR_INT_DATA:
    errStr = "WLZ_ERR_INT_DATA",
      msgStr = "Int data";
    break;
  case WLZ_ERR_SHORT_DATA:
    errStr = "WLZ_ERR_SHORT_DATA",
      msgStr = "Short int data";
    break;
  case WLZ_ERR_UBYTE_DATA:
    errStr = "WLZ_ERR_UBYTE_DATA",
      msgStr = "Unsigned byte data";
    break;
  case WLZ_ERR_FLOAT_DATA:
    errStr = "WLZ_ERR_FLOAT_DATA",
      msgStr = "Float data";
    break;
  case WLZ_ERR_DOUBLE_DATA:
    errStr = "WLZ_ERR_DOUBLE_DATA",
      msgStr = "Double data";
    break;
  case WLZ_ERR_GREY_TYPE:
    errStr = "WLZ_ERR_GREY_TYPE",
      msgStr = "Invalid grey type";
    break;
  case WLZ_ERR_GREY_DATA:
    errStr = "WLZ_ERR_GREY_DATA",
      msgStr = "Grey data";
    break;
  case WLZ_ERR_PLANEDOMAIN_TYPE:
    errStr = "WLZ_ERR_PLANEDOMAIN_TYPE",
      msgStr = "Invalid planedomain type";
    break;
  case WLZ_ERR_PLANEDOMAIN_DATA:
    errStr = "WLZ_ERR_PLANEDOMAIN_DATA",
      msgStr = "";
    break;
  case WLZ_ERR_INTERVALDOMAIN_NULL:
    errStr = "WLZ_ERR_INTERVALDOMAIN_NULL",
      msgStr = "NULL intervaldomain";
    break;
  case WLZ_ERR_INTERVALDOMAIN_TYPE:
    errStr = "WLZ_ERR_INTERVALDOMAIN_TYPE",
      msgStr = "Invalid intervaldomain type";
    break;
  case WLZ_ERR_INTERVALLINE_NULL:
    errStr = "WLZ_ERR_INTERVALLINE_NULL",
      msgStr = "NULL intervalline";
    break;
  case WLZ_ERR_INTERVAL_NULL:
    errStr = "WLZ_ERR_INTERVAL_NULL",
      msgStr = "NULL interval";
    break;
  case WLZ_ERR_INTERVAL_DATA:
    errStr = "WLZ_ERR_INTERVAL_DATA",
      msgStr = "Interval data";
    break;
  case WLZ_ERR_INTERVAL_ADJACENT:
    errStr = "WLZ_ERR_INTERVAL_ADJACENT",
      msgStr = "Adjacent interval";
    break;
  case WLZ_ERR_INTERVAL_BOUND:
    errStr = "WLZ_ERR_INTERVAL_BOUND",
      msgStr = "Interval bounds";
    break;
  case WLZ_ERR_INTERVAL_NUMBER:
    errStr = "WLZ_ERR_INTERVAL_NUMBER",
      msgStr = "Interval number";
    break;
  case WLZ_ERR_TRANSFORM_DATA:
    errStr = "WLZ_ERR_TRANSFORM_DATA",
      msgStr = "Transform data";
    break;
  case WLZ_ERR_TRANSFORM_TYPE:
    errStr = "WLZ_ERR_TRANSFORM_TYPE",
      msgStr = "Invalid transform type";
    break;
  case WLZ_ERR_VOXELVALUES_TYPE:
    errStr = "WLZ_ERR_VOXELVALUES_TYPE",
      msgStr = "Invalid voxelvalues";
    break;
  case WLZ_ERR_COLUMN_DATA:
    errStr = "WLZ_ERR_COLUMN_DATA",
      msgStr = "Column data";
    break;
  case WLZ_ERR_LINE_DATA:
    errStr = "WLZ_ERR_LINE_DATA",
      msgStr = "Line data";
    break;
  case WLZ_ERR_PLANE_DATA:
    errStr = "WLZ_ERR_PLANE_DATA",
      msgStr = "Plane data";
    break;
  case WLZ_ERR_BINARY_OPERATOR_TYPE:
    errStr = "WLZ_ERR_BINARY_OPERATOR_TYPE",
      msgStr = "Invalid binary operator type";
    break;
  case WLZ_ERR_COMPTHRESH_TYPE:
    errStr = "WLZ_ERR_COMPTHRESH_TYPE",
      msgStr = "Invalid computed threshold type";
    break;
  case WLZ_ERR_CONNECTIVITY_TYPE:
    errStr = "WLZ_ERR_CONNECTIVITY_TYPE",
      msgStr = "Invalid connectivity type";
    break;
  case WLZ_ERR_INTERPOLATION_TYPE:
    errStr = "WLZ_ERR_INTERPOLATION_TYPE",
      msgStr = "Invalid interpolation type";
    break;
  case WLZ_ERR_POLYGON_TYPE:
    errStr = "WLZ_ERR_POLYGON_TYPE",
      msgStr = "Invalid polygon type";
    break;
  case WLZ_ERR_RASTERDIR_TYPE:
    errStr = "WLZ_ERR_RASTERDIR_TYPE",
      msgStr = "Invalid raster direction type";
    break;
  case WLZ_ERR_LINKCOUNT_DATA:
    errStr = "WLZ_ERR_LINKCOUNT_DATA",
      msgStr = "Linkcount data";
    break;
  case WLZ_ERR_MEM_ALLOC:
    errStr = "WLZ_ERR_MEM_ALLOC",
      msgStr = "Memory allocation";
    break;
  case WLZ_ERR_MEM_FREE:
    errStr = "WLZ_ERR_MEM_FREE",
      msgStr = "Memory de-allocation";
    break;
  case WLZ_ERR_READ_EOF:
    errStr = "WLZ_ERR_READ_EOF",
      msgStr = "End of file when reading an object";
    break;
  case WLZ_ERR_READ_INCOMPLETE:
    errStr = "WLZ_ERR_READ_INCOMPLETE",
      msgStr = "Incomplete object read";
    break;
  case WLZ_ERR_WRITE_EOF:
    errStr = "WLZ_ERR_WRITE_EOF",
      msgStr = "End of file when writing an object";
    break;
  case WLZ_ERR_WRITE_INCOMPLETE:
    errStr = "WLZ_ERR_WRITE_INCOMPLETE",
      msgStr = "Incomplete object write";
    break;
  case WLZ_ERR_ALG:
    errStr = "WLZ_ERR_ALG",
      msgStr = "Numerical algorithm error";
    break;
  case WLZ_ERR_ALG_SINGULAR:
    errStr = "WLZ_ERR_ALG_SINGULAR",
      msgStr = "Numerical algorithm singular matrix";
    break;
  case WLZ_ERR_ALG_HOMOGENEOUS:
    errStr = "WLZ_ERR_ALG_HOMOGENEOUS",
      msgStr = "Numerical algorithm homogeneous matrix";
    break;
  case WLZ_ERR_ALG_CONVERGENCE:
    errStr = "WLZ_ERR_ALG_CONVERGENCE",
      msgStr = "Numerical algorithm failed to converge";
    break;
  case WLZ_ERR_UNIMPLEMENTED:
    errStr = "WLZ_ERR_UNIMPLEMENTED",
      msgStr = "Unimplemented feature";
    break;
  case WLZ_ERR_FILE_OPEN:
    errStr = "WLZ_ERR_FILE_OPEN",
      msgStr = "File open failure";
    break;
  case WLZ_ERR_FILE_FORMAT:
    errStr = "WLZ_ERR_FILE_FORMAT",
      msgStr = "Invalid file format";
    break;
  case WLZ_ERR_IMAGE_TYPE:
    errStr = "WLZ_ERR_IMAGE_TYPE",
      msgStr = "Attempt to read an image type not-supported in Woolz";
    break;
  default:	/* FALLTHROUGH */
  case WLZ_ERR_UNSPECIFIED:
    errStr = "WLZ_ERR_UNSPECIFIED",
      msgStr = "Unspecified error";
    break;
  }
  if(dstMsgStr)
  {
    *dstMsgStr = msgStr;
  }
  return(errStr);
}

/************************************************************************
* Function:	WlzStringToErrorNum				
* Returns:	WlzErrorNum:		Matched error number.	
* Purpose:	Finds an error number for the given error number
*		string.						
* Global refs:	-						
* Parameters:	const char *errStr:	Given error number string.
************************************************************************/
WlzErrorNum	WlzStringToErrorNum(const char *errStr)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_UNSPECIFIED;

  if(WlzStringMatchValue(&tI0, errStr,
	      "WLZ_ERR_NONE", WLZ_ERR_NONE,
	      "WLZ_ERR_EOO", WLZ_ERR_EOO,
	      "WLZ_ERR_OBJECT_NULL", WLZ_ERR_OBJECT_NULL,
	      "WLZ_ERR_OBJECT_TYPE", WLZ_ERR_OBJECT_TYPE,
	      "WLZ_ERR_OBJECT_DATA", WLZ_ERR_OBJECT_DATA,
	      "WLZ_ERR_DOMAIN_NULL", WLZ_ERR_DOMAIN_NULL ,
	      "WLZ_ERR_DOMAIN_TYPE", WLZ_ERR_DOMAIN_TYPE,
	      "WLZ_ERR_DOMAIN_DATA", WLZ_ERR_DOMAIN_DATA,
	      "WLZ_ERR_VALUES_NULL", WLZ_ERR_VALUES_NULL,
	      "WLZ_ERR_VALUES_TYPE", WLZ_ERR_VALUES_TYPE,
	      "WLZ_ERR_VALUES_DATA", WLZ_ERR_VALUES_DATA,
	      "WLZ_ERR_PROPERTY_NULL", WLZ_ERR_PROPERTY_NULL,
	      "WLZ_ERR_GMELM_NULL", WLZ_ERR_GMELM_NULL,
	      "WLZ_ERR_GMELM_TYPE", WLZ_ERR_GMELM_TYPE,
	      "WLZ_ERR_GMELM_DATA", WLZ_ERR_GMELM_DATA,
	      "WLZ_ERR_PARAM_NULL", WLZ_ERR_PARAM_NULL,
	      "WLZ_ERR_PARAM_TYPE", WLZ_ERR_PARAM_TYPE,
	      "WLZ_ERR_PARAM_DATA", WLZ_ERR_PARAM_DATA,
	      "WLZ_ERR_INT_DATA", WLZ_ERR_INT_DATA,
	      "WLZ_ERR_SHORT_DATA", WLZ_ERR_SHORT_DATA,
	      "WLZ_ERR_UBYTE_DATA", WLZ_ERR_UBYTE_DATA,
	      "WLZ_ERR_FLOAT_DATA", WLZ_ERR_FLOAT_DATA,
	      "WLZ_ERR_DOUBLE_DATA", WLZ_ERR_DOUBLE_DATA,
	      "WLZ_ERR_GREY_TYPE", WLZ_ERR_GREY_TYPE,
	      "WLZ_ERR_GREY_DATA", WLZ_ERR_GREY_DATA,
	      "WLZ_ERR_PLANEDOMAIN_TYPE", WLZ_ERR_PLANEDOMAIN_TYPE,
	      "WLZ_ERR_PLANEDOMAIN_DATA", WLZ_ERR_PLANEDOMAIN_DATA,
	      "WLZ_ERR_INTERVALDOMAIN_NULL", WLZ_ERR_INTERVALDOMAIN_NULL,
	      "WLZ_ERR_INTERVALDOMAIN_TYPE", WLZ_ERR_INTERVALDOMAIN_TYPE,
	      "WLZ_ERR_INTERVALLINE_NULL", WLZ_ERR_INTERVALLINE_NULL,
	      "WLZ_ERR_INTERVAL_NULL", WLZ_ERR_INTERVAL_NULL,
	      "WLZ_ERR_INTERVAL_DATA", WLZ_ERR_INTERVAL_DATA,
	      "WLZ_ERR_INTERVAL_ADJACENT", WLZ_ERR_INTERVAL_ADJACENT,
	      "WLZ_ERR_INTERVAL_BOUND", WLZ_ERR_INTERVAL_BOUND,
	      "WLZ_ERR_INTERVAL_NUMBER", WLZ_ERR_INTERVAL_NUMBER,
	      "WLZ_ERR_TRANSFORM_DATA", WLZ_ERR_TRANSFORM_DATA,
	      "WLZ_ERR_TRANSFORM_TYPE", WLZ_ERR_TRANSFORM_TYPE,
	      "WLZ_ERR_VOXELVALUES_TYPE", WLZ_ERR_VOXELVALUES_TYPE,
	      "WLZ_ERR_COLUMN_DATA", WLZ_ERR_COLUMN_DATA,
	      "WLZ_ERR_LINE_DATA", WLZ_ERR_LINE_DATA,
	      "WLZ_ERR_PLANE_DATA", WLZ_ERR_PLANE_DATA,
	      "WLZ_ERR_BINARY_OPERATOR_TYPE", WLZ_ERR_BINARY_OPERATOR_TYPE,
	      "WLZ_ERR_COMPTHRESH_TYPE", WLZ_ERR_COMPTHRESH_TYPE,
	      "WLZ_ERR_CONNECTIVITY_TYPE", WLZ_ERR_CONNECTIVITY_TYPE,
	      "WLZ_ERR_INTERPOLATION_TYPE", WLZ_ERR_INTERPOLATION_TYPE,
	      "WLZ_ERR_POLYGON_TYPE", WLZ_ERR_POLYGON_TYPE,
	      "WLZ_ERR_RASTERDIR_TYPE", WLZ_ERR_RASTERDIR_TYPE,
	      "WLZ_ERR_LINKCOUNT_DATA", WLZ_ERR_LINKCOUNT_DATA,
	      "WLZ_ERR_MEM_ALLOC", WLZ_ERR_MEM_ALLOC,
	      "WLZ_ERR_MEM_FREE", WLZ_ERR_MEM_FREE,
	      "WLZ_ERR_READ_EOF", WLZ_ERR_READ_EOF,
	      "WLZ_ERR_READ_INCOMPLETE", WLZ_ERR_READ_INCOMPLETE,
	      "WLZ_ERR_WRITE_EOF", WLZ_ERR_WRITE_EOF,
	      "WLZ_ERR_WRITE_INCOMPLETE", WLZ_ERR_WRITE_INCOMPLETE,
	      "WLZ_ERR_ALG", WLZ_ERR_ALG,
	      "WLZ_ERR_ALG_SINGULAR", WLZ_ERR_ALG_SINGULAR,
	      "WLZ_ERR_ALG_HOMOGENEOUS", WLZ_ERR_ALG_HOMOGENEOUS,
	      "WLZ_ERR_ALG_CONVERGENCE", WLZ_ERR_ALG_CONVERGENCE,
	      "WLZ_ERR_UNIMPLEMENTED", WLZ_ERR_UNIMPLEMENTED,
	      "WLZ_ERR_FILE_OPEN", WLZ_ERR_FILE_OPEN,
	      "WLZ_ERR_FILE_FORMAT", WLZ_ERR_FILE_FORMAT,
	      "WLZ_ERR_IMAGE_TYPE", WLZ_ERR_IMAGE_TYPE,
	      "WLZ_ERR_UNSPECIFIED", WLZ_ERR_UNSPECIFIED,
	      NULL))
  {
    errNum = tI0;
  }
  return(errNum);
}
