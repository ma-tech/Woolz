#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzUnion2.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Procedure to calculate the union of two Woolz domain
*		objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzUnion2						*
*   Date       : Mon Oct 28 18:26:14 1996				*
*************************************************************************
*   Synopsis   :Convenience procedure to calculate the union of two	*
*		woolz domain objects. Note uvt set to zero		*
*   Returns    :WlzObject *: object pointer, NULL on error.		*
*   Parameters :WlzObject *obj1, *obj2: two domain objects.		*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzUnion2(WlzObject *obj1,
		     WlzObject *obj2,
		     WlzErrorNum *dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;

  return WlzUnionN(2, objs, 0, dstErr);
}
