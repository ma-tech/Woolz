#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzValueTableUtils.c
* \author       Bill Hill, Richard Baldock
* \date         March 1999
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
* \ingroup	WlzValuesUtils
* \brief	Functions for computing value amd value table types.
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Type of grey table.
* \ingroup	WlzValuesUtils
* \brief	Computes a grey table type from table and grey types.
* \param	tableType		The basic table type.
* \param	greyType		The grey type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzGreyTableType(WlzObjectType tableType,
			         WlzGreyType greyType,
				 WlzErrorNum *dstErr)
{
  WlzObjectType gTabType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(tableType)
  {
    case WLZ_GREY_TAB_RAGR:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RAGR_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RAGR_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RAGR_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RAGR_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RAGR_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_RAGR_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_RECT:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_RECT_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_RECT_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_RECT_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_RECT_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_RECT_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_RECT_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_GREY_TAB_INTL:
      switch(greyType)
      {
        case WLZ_GREY_INT:
	  gTabType = WLZ_VALUETABLE_INTL_INT;
	  break;
        case WLZ_GREY_SHORT:
	  gTabType = WLZ_VALUETABLE_INTL_SHORT;
	  break;
        case WLZ_GREY_UBYTE:
	  gTabType = WLZ_VALUETABLE_INTL_UBYTE;
	  break;
        case WLZ_GREY_FLOAT:
	  gTabType = WLZ_VALUETABLE_INTL_FLOAT;
	  break;
        case WLZ_GREY_DOUBLE:
	  gTabType = WLZ_VALUETABLE_INTL_DOUBLE;
	  break;
        case WLZ_GREY_RGBA:
	  gTabType = WLZ_VALUETABLE_INTL_RGBA;
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    default:
      errNum = WLZ_ERR_VALUES_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(gTabType);
}

/*!
* \return	Type of grey value.
* \ingroup	WlzValuesUtils
* \brief	Computes the type of grey from a grey table type.
* \param	gTabType		The basic table type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGreyType	WlzGreyTableTypeToGreyType(WlzObjectType gTabType,
				           WlzErrorNum *dstErr)
{
  WlzGreyType	greyType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_INT:
      greyType = WLZ_GREY_INT;
      break;
    case WLZ_VALUETABLE_RAGR_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT:
      greyType = WLZ_GREY_SHORT;
      break;
    case WLZ_VALUETABLE_RAGR_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE:
      greyType = WLZ_GREY_UBYTE;
      break;
    case WLZ_VALUETABLE_RAGR_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT:
      greyType = WLZ_GREY_FLOAT;
      break;
    case WLZ_VALUETABLE_RAGR_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE:
      greyType = WLZ_GREY_DOUBLE;
      break;
    case WLZ_VALUETABLE_RAGR_RGBA: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_RGBA: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_RGBA:
      greyType = WLZ_GREY_RGBA;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(greyType);
}

/*!
* \return	Type of grey value.
* \ingroup	WlzValuesUtils
* \brief	Computes the type of table from a  grey table type.
* \param	gTabType		The basic table type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType WlzGreyTableTypeToTableType(WlzObjectType gTabType,
				 	  WlzErrorNum *dstErr)
{
  WlzObjectType	tableType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(gTabType)
  {
    case WLZ_VALUETABLE_RAGR_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RAGR_RGBA:
      tableType = WLZ_GREY_TAB_RAGR;
      break;
    case WLZ_VALUETABLE_RECT_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_RECT_DOUBLE /* FALLTHROUGH */:
    case WLZ_VALUETABLE_RECT_RGBA:
      tableType = WLZ_GREY_TAB_RECT;
      break;
    case WLZ_VALUETABLE_INTL_INT:   /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_SHORT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_UBYTE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_FLOAT: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_DOUBLE: /* FALLTHROUGH */
    case WLZ_VALUETABLE_INTL_RGBA:
      tableType = WLZ_GREY_TAB_INTL;
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tableType);
}

/*!
* \return	Grey type, WLZ_GREY_ERROR on error.
* \ingroup	WlzValuesUtils
* \brief	Gets the grey type of the values in a Woolz object.
*               If the object is not a domain object with values
*               an error is returned. If the object is a 3D domain
*               object with values, all 2D value tables are checked
*               and an error is returned if they don't all have the
*               same grey type.
* \param	obj			Given object.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
WlzGreyType	WlzGreyTypeFromObj(WlzObject *obj, WlzErrorNum *dstErr)
{
  int		pCnt,
    pIdx,
    firstDomFlg;
  WlzDomain	dom2D;
  WlzValues	val2D;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(obj->type)
    {
    case WLZ_2D_DOMAINOBJ:
      gType = WlzGreyTableTypeToGreyType(obj->values.core->type, &errNum);
      break;

    case WLZ_3D_DOMAINOBJ:
      pIdx = 1;
      firstDomFlg = 1;
      pCnt = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
      while((errNum == WLZ_ERR_NONE) && (pIdx < pCnt))
      {
	if(((dom2D = *(obj->domain.p->domains + pIdx)).core != NULL) &&
	   (dom2D.core->type != WLZ_EMPTY_DOMAIN))
	{
	  if(((val2D = *(obj->values.vox->values + pIdx)).core == NULL) ||
	     (val2D.core->type == WLZ_EMPTY_VALUES))
	  {
	    errNum = WLZ_ERR_VALUES_NULL;
	  }
	  else
	  {
	    if(firstDomFlg)
	    {
	      gType = WlzGreyTableTypeToGreyType(val2D.core->type, &errNum);
	    }
	    else
	    {
	      if(gType != WlzGreyTableTypeToGreyType(val2D.core->type,
						     &errNum))
	      {
		errNum = WLZ_ERR_VALUES_DATA;
	      }
	    }
	  }
	}
	++pIdx;
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
  return(gType);
}

/*!
* \return	Voxel size.
* \ingroup	WlzValuesUtils
* \brief	Gets the given 3D domain objects voxel size.
* \param	obj			Given 3D domain object.
* \param	dstErr			Destination pointer for error code,
* 					may be NULL.
*/
WlzDVertex3	WlzVozelSz(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzDVertex3	sz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  sz.vtX = sz.vtY = sz.vtZ = 0;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    sz.vtX = obj->domain.p->voxel_size[0];
    sz.vtY = obj->domain.p->voxel_size[1];
    sz.vtZ = obj->domain.p->voxel_size[2];
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(sz);
}
