#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntervalCount.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Counts the number of intervals (or equivalent) in a
*		Woolz object's domain.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>

/************************************************************************
*   Function   : WlzIntervalCount					*
*   Returns    :int: number of intervals if WLZ_INTERVALDOMAIN_INTVL	*
*		else number of lines					*
*   Parameters :WlzIntervalDomain *idom: woolz interval domain		*
*		WlzErrorNum *wlzErr:     Woolz error return.		*
*   Date       : Mon Oct 14 13:18:00 1996				*
*   Synopsis   :Count the number of intervals or equivalent if 		*
*		rectangular						*
************************************************************************/

int WlzIntervalCount(WlzIntervalDomain *idom, WlzErrorNum *wlzErr)
{
  int		l,
		ll,
		intcount = -1;
  WlzIntervalLine *intl;
  WlzErrorNum	errNum = WLZ_ERR_UNSPECIFIED;

  if(idom == NULL)
  {
    errNum = WLZ_ERR_INTERVALDOMAIN_NULL;
  }
  else
  {
    switch(idom->type)
    {
      case WLZ_INTERVALDOMAIN_INTVL:
	intl = idom->intvlines;
	intcount = 0;
	ll = idom->lastln;
	for (l = idom->line1; l <= ll; ++l)
	{
	  intcount += (intl++)->nintvs;
	}
	errNum = WLZ_ERR_NONE;
	break;
      case WLZ_INTERVALDOMAIN_RECT:
	errNum = WLZ_ERR_NONE;
	intcount = idom->lastln - idom->line1 + 1;
	break;
      default:
	errNum = WLZ_ERR_INTERVALDOMAIN_TYPE;
	break;
    }
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  return(intcount);
}
