#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntersect2.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the set intersection between two Woolz
*		domain objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzIntersect2						*
*   Date       : Mon Oct 28 18:23:36 1996				*
*************************************************************************
*   Synopsis   :Calculate the set intersection between two domain 	*
*		objects.						*
*   Returns    :WlzObject *: return from WlzIntersectN. NULL on error.	*
*   Parameters :WlzObject *obj1: first object pointer			*
*		WlzObject *obj2: second object pointer			*
*		Note WLZ_EMPTY_OBJ is legal but a NULL pointer is not.	*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzIntersect2(WlzObject *obj1,
			 WlzObject *obj2,
			 WlzErrorNum	*dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;
  return(WlzIntersectN(2, objs, 0, dstErr));
}
