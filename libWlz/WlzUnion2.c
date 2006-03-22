#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzUnion2_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzUnion2.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Convenience function to calculate the union
* 		of two domain objects.
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <Wlz.h>

/* function:     WlzUnion2    */
/*! 
* \return       Union object pointer.
* \ingroup      WlzBinaryOps
* \brief        Convenience procedure to calculate the union of two
* 		woolz domain objects. This calls WlzUnnionN() with
* 		uvt=0. Objects must be of the same type.
* \par      Source:
*                WlzUnion2.c
* \param    obj1	First inout object.
* \param    obj2	Second input object.
* \param    dstErr	errro return.
*/
WlzObject *WlzUnion2(WlzObject *obj1,
		     WlzObject *obj2,
		     WlzErrorNum *dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;

  return WlzUnionN(2, objs, 0, dstErr);
}
