#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzPolyReverse.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Mon Jul 30 12:07:51 2001
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions to reverse the vertex ordering in a polygon domain
*               
* \todo         -
* \bug          None known
* \ingroup      WlzPolyline
*
* Maintenance log with most recent changes at top of list.
*/

/*
 * polyinvert.c		Jim Piper	November 30 1983
 * polygon domain operations
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
