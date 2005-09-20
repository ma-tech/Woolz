#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzIntersect2.c
* \author       Richard Baldock
* \date         August 2003
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
* \brief	Calculates the intersection between two domain objects.
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/* function:     WlzIntersect2    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Calculate the set intersection between two domain
 objects. This is a convenience routine calling WlzIntersectN with
 uvt=0. Input objects must be domain objects of the same type (2D
 or 3D) and non-NULL. Type WLZ_EMPTY_OBJ is legal, clearly an empty
 domain will be returned.
*
* \return       Intersection object with NULL value table
* \param    obj1	first input object
* \param    obj2	second input object
* \param    dstErr	error return.
* \par      Source:
*                WlzIntersect2.c
*/
WlzObject *WlzIntersect2(
  WlzObject *obj1,
  WlzObject *obj2,
  WlzErrorNum	*dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;
  return(WlzIntersectN(2, objs, 0, dstErr));
}
