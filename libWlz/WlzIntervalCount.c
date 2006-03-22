#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzIntervalCount_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzIntervalCount.c
* \author       Bill Hill
* \date         September 2003
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Counts the number of intervals (or equivalent)
* 		in an object's domain.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
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
