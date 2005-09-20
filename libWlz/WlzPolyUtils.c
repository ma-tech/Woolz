#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzPolyUtils.c
* \author       Bill Hill
* \date         September 2002
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
* \brief	Functions for manipulating polygon domains.
* \ingroup	WlzPolyline
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzPolyline
* \brief	Access to an polygon domain's integer vertices.
* \param	poly			Given polygon domain.
* \param	dstArySz		Destination pointer for the number of
* 					vertices.
* \param	dstPolyAry		Destination pointer for the vertices.
*/
WlzErrorNum	WlzPolyVertices2I(WlzPolygonDomain *poly,
				  int *dstArySz, WlzIVertex2 **dstPolyAry)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(poly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstArySz = poly->nvertices;
    *dstPolyAry = poly->vtx;
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzPolyline
* \brief	Access to an polygon domain's double precision vertices.
* \param	poly			Given polygon domain.
* \param	dstArySz		Destination pointer for the number of
* 					vertices.
* \param	dstPolyAry		Destination pointer for the vertices.
*/
WlzErrorNum	WlzPolyVertices2D(WlzPolygonDomain *poly,
				  int *dstArySz, WlzDVertex2 **dstPolyAry)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(poly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstArySz = poly->nvertices;
    *dstPolyAry = (WlzDVertex2 *)(poly->vtx);
  }
}


