#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzPolyReverse.c
* \author       Jim Piper, Richard Baldock
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
* \brief	Functions to reverse the vertex ordering in a polygon
* 		domain.
* \ingroup	WlzPolyline
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <Wlz.h>

/* function:     WlzPolyReverse   */
/*! 
* \ingroup	WlzPolyline
* \brief        Reverse the vertex ordering in a polygon domain.
*
* \return       New polygon domain with vertices reversed
* \param    poly	input polygon domain
* \param    dstErr	error return
* \par      Source:
*               WlzPolyReverse.c
*/
WlzPolygonDomain *WlzPolyReverse(
  WlzPolygonDomain 	*poly,
  WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*rtnPoly=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzIVertex2		*iVtxs, *jVtxs;
  WlzFVertex2		*fVtxs, *gVtxs;
  WlzDVertex2		*dVtxs, *eVtxs;
  int			i, n;
  double		x, y;

  /* check object and parameters */
  if( poly == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    if( rtnPoly = WlzMakePolygonDomain(poly->type, poly->nvertices, NULL,
				 poly->nvertices, 1, &errNum) ){
      switch( rtnPoly->type ){
      case WLZ_POLYGON_INT:
	n = rtnPoly->nvertices;
	iVtxs = rtnPoly->vtx;
	jVtxs = iVtxs + n - 1;
	n = (n+1)/2;
	for(i=0; i < n; i++, iVtxs++, jVtxs--) {
	  x = iVtxs->vtX;
	  y = iVtxs->vtY;
	  iVtxs->vtX = jVtxs->vtX;
	  iVtxs->vtY = jVtxs->vtY;
	  jVtxs->vtX = x;
	  jVtxs->vtY = y;
	}
	break;

      case WLZ_POLYGON_FLOAT:
	n = rtnPoly->nvertices;
	fVtxs = (WlzFVertex2 *) rtnPoly->vtx;
	gVtxs = fVtxs + n - 1;
	n = (n+1)/2;
	for(i=0; i < n; i++, fVtxs++, gVtxs--) {
	  x = fVtxs->vtX;
	  y = fVtxs->vtY;
	  fVtxs->vtX = gVtxs->vtX;
	  fVtxs->vtY = gVtxs->vtY;
	  gVtxs->vtX = x;
	  gVtxs->vtY = y;
	}
	break;

      case WLZ_POLYGON_DOUBLE:
	n = rtnPoly->nvertices;
	dVtxs = (WlzDVertex2 *) rtnPoly->vtx;
	eVtxs = dVtxs + n - 1;
	n = (n+1)/2;
	for(i=0; i < n; i++, dVtxs++, eVtxs--) {
	  x = dVtxs->vtX;
	  y = dVtxs->vtY;
	  dVtxs->vtX = eVtxs->vtX;
	  dVtxs->vtY = eVtxs->vtY;
	  eVtxs->vtX = x;
	  eVtxs->vtY = y;
	}
	break;

      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnPoly;
}
