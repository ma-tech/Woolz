#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGreyNormalise.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Normalise the grey-values of a woolz object to the
*		range 0-255.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>


/************************************************************************
*   Function   : WlzGreyNormalise					*
*   Date       : Wed Sep 17 11:21:15 1997				*
*************************************************************************
*   Synopsis   :normalise grey-table of object to lie in range 1-255	*
*   Returns    :							*
*   Parameters :							*
*   Global refs:							*
************************************************************************/

WlzErrorNum WlzGreyNormalise(
  WlzObject	*obj)
{
  WlzPixelV	min, max, Min, Max;
  WlzErrorNum	wlzErrno=WLZ_ERR_NONE;

  /* get then set the grey-range of the object */
  wlzErrno = WlzGreyRange(obj, &min, &max);
  Min = min;
  Max = max;
  if( wlzErrno == WLZ_ERR_NONE ){
    switch( min.type ){
    case WLZ_GREY_INT:
      Min.v.inv = 0;
      Max.v.inv = 255;
      break;
    case WLZ_GREY_SHORT:
      Min.v.shv = 0;
      Max.v.shv = 255;
      break;
    case WLZ_GREY_UBYTE:
      Min.v.ubv = 0;
      Max.v.ubv = 255;
      break;
    case WLZ_GREY_FLOAT:
      Min.v.flv = 0.0;
      Max.v.flv = 255.0;
      break;
    case WLZ_GREY_DOUBLE:
      Min.v.dbv = 0.0;
      Max.v.dbv = 255.0;
      break;
    }
    wlzErrno = WlzGreySetRange(obj, min, max, Min, Max);
  }

  return wlzErrno;
}
