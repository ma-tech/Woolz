#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzArea.c
* \author       Richard Baldock
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
* \brief	Computes the area of an object.
* \ingroup	WlzFeatures
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>

/*!
* \return	Area of object or a negative value on error.
* \ingroup	WlzFeatures
* \brief	Computes the area of an object.
* \todo		Implement for objects other than 2D domain objects.
* \param	obj			Input object.
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
int WlzArea(
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  int 			size;
  WlzIntervalWSpace 	iwsp;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == (WlzObject *) NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
    size = -1;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	size = -1;
      }
      break;
    
    case WLZ_EMPTY_OBJ:
      size = 0;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      size = -1;
      break;

    }
  }

  /* calculate size */
  if( errNum == WLZ_ERR_NONE ){
    size = 0;
    errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC);
    if( errNum == WLZ_ERR_NONE ){
      while( (errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	size += iwsp.colrmn;
      }
      if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return size;
}
