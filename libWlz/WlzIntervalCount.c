#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzIntervalCount.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Sep 26 11:15:09 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzDomainOps
* \brief        Counts the number of intervals (or equivalent)
 in a Woolz object domain.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

/* function:     WlzIntervalCount    */
/*! 
* \ingroup      WlzDomainOps
* \brief        Count the number of intervals or equivalent if rectangular.
*
* \return       Number of intervals if the domain type is
 <tt>WLZ_INTERVALDOMAIN_INTVL</tt> else the number of lines.
* \param    idom	Input domain.
* \param    wlzErr	Error return.
* \par      Source:
*                WlzIntervalCount.c
*/
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
