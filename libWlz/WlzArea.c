#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzArea.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes the area of an object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzArea						*
*   Returns    : int - area (number of pixels) of an object		*
*		       -1 on error					*
*   Parameters : WlzObject *obj: input object				*
*   Date       : Fri Oct 11 12:10:26 1996				*
*   Synopsis   : return the area of an object currently only implemented*
*	     	for 2d domain objects - to be extended			*
************************************************************************/

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
