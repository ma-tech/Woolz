#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzPolyEquispace.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Thu Jul 12 11:23:23 2001
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
* \brief        Build a new polygon domain with equi-spaced vertices from the
old.
*               
* \todo         put in keepVertices option in WlzPolyEquispace
* \bug          None known
* \ingroup      WlzPolyline
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <Wlz.h>


/* function:     WlzPolyLength    */
/*! 
* \ingroup      WlzPolyline
* \brief        Calculate the length of the input polyline.
*
* \return       length of the input polygon domain
* \param    poly	input polygon domain
* \param    wrap	wrap value of the p[olygon
* \param    dstErr	error return
* \par      Source:
*               WlzPolyEquispace.c
*/
double WlzPolyLength(
  WlzPolygonDomain	*poly,
  int			wrap,
  WlzErrorNum		*dstErr)
{
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzIVertex2		*iVtxs;
  WlzFVertex2		*fVtxs;
  WlzDVertex2		*dVtxs; 
  double		length=-1.0;
  double		x0, y0, dx, dy;
  int			i, nVtxs;

  /* check inputs */
  if( poly == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if( wrap < 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  if( errNum == WLZ_ERR_NONE ){
    length = 0.0;
    switch( poly->type ){
    case WLZ_POLYGON_INT:
      iVtxs = poly->vtx;
      nVtxs = wrap ? poly->nvertices - wrap + 1 : poly->nvertices;
      x0 = iVtxs->vtX;
      y0 = iVtxs->vtY;
      iVtxs++;
      for(i=1; i < nVtxs; i++, iVtxs++){
	dx = iVtxs->vtX - x0;
	x0 = iVtxs->vtX;
	dy = iVtxs->vtY - y0;
	y0 = iVtxs->vtY;
	length += sqrt(dx*dx + dy*dy);
      }
      break;

    case WLZ_POLYGON_FLOAT:
      fVtxs = (WlzFVertex2 *) poly->vtx;
      nVtxs = wrap ? poly->nvertices - wrap + 1 : poly->nvertices;
      x0 = fVtxs->vtX;
      y0 = fVtxs->vtY;
      fVtxs++;
      for(i=1; i < nVtxs; i++, fVtxs++){
	dx = fVtxs->vtX - x0;
	x0 = fVtxs->vtX;
	dy = fVtxs->vtY - y0;
	y0 = fVtxs->vtY;
	length += sqrt(dx*dx + dy*dy);
      }
      break;

    case WLZ_POLYGON_DOUBLE:
      dVtxs = (WlzDVertex2 *) poly->vtx;
      nVtxs = wrap ? poly->nvertices - wrap + 1 : poly->nvertices;
      x0 = dVtxs->vtX;
      y0 = dVtxs->vtY;
      dVtxs++;
      for(i=1; i < nVtxs; i++, dVtxs++){
	dx = dVtxs->vtX - x0;
	x0 = dVtxs->vtX;
	dy = dVtxs->vtY - y0;
	y0 = dVtxs->vtY;
	length += sqrt(dx*dx + dy*dy);
      }
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return length;
}


/* function:     WlzPolyEquispace    */
/*! 
* \ingroup      WlzPolyline
* \brief        Create a new polygon domain with vertices qually spaced
along the original polyline. The wrap value is preserved and if keepOrigVtxs
is non-zero then the original vertices will be kept. This results in non-equal
spacing but the new line will be faithfull to the old and not "cut-corners".
*
* \return       new polygon domain with equally spaced vertices and 
the given wrap value
* \param    poly	input polygon domain
* \param    wrap	wrap value of the input polygon
* \param    spacing	required spacing
* \param    keepOrigVtxs	flag to retain the original vertices, 
currently ignored
* \param    dstErr	error return
* \par      Source:
*               WlzPolyEquispace.c
*/
WlzPolygonDomain  *WlzPolyEquispace(
  WlzPolygonDomain	*poly,
  int			wrap,
  double		spacing,
  int			keepOrigVtxs,
  WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*rtnPoly=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			i, n, nVtxs;
  WlzIVertex2		*iVtxs, *jVtxs;
  WlzFVertex2		*fVtxs, *gVtxs;
  WlzDVertex2		*dVtxs, *eVtxs; 
  double		length, dx, dy, del, dist;
  double		x, y;

  /* check inputs */
  if( poly == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if( wrap < 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if( spacing < 0.5 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if( (poly->type == WLZ_POLYGON_INT) && (spacing < 1.0) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* calculate length and generate new polydomain
     leave space for kept vertices if required - not yet implemented */
  if( errNum == WLZ_ERR_NONE ){
    length = WlzPolyLength(poly, wrap, &errNum);
    n = length / spacing + wrap + poly->nvertices;
    rtnPoly = WlzMakePolygonDomain(poly->type, 0, NULL, n, 1, &errNum);
  }

  /* set values, interpolate as required */
  if( errNum == WLZ_ERR_NONE ){
    switch( poly->type ){
    case WLZ_POLYGON_INT:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      iVtxs = (WlzIVertex2 *) poly->vtx;
      jVtxs = (WlzIVertex2 *) rtnPoly->vtx;
      jVtxs[0] = iVtxs[0];
      n = 1;
      i = 1;
      del = spacing;
      dx = iVtxs[i].vtX - iVtxs[i-1].vtX;
      dy = iVtxs[i].vtY - iVtxs[i-1].vtY;
      dist = sqrt(dx*dx + dy*dy);
      while( i < nVtxs ){
	if( del <= dist ){
	  /* new point in this span */
	  x = iVtxs[i-1].vtX + dx*del/dist;
	  y = iVtxs[i-1].vtY + dy*del/dist;
	  jVtxs[n].vtX = WLZ_NINT(x);
	  jVtxs[n].vtY = WLZ_NINT(y);
	  del += spacing;
	  n++;
	}
	else {
	  /* next span */
	  i++;
	  del -= dist;
	  if( i < nVtxs ){
	    /* should put in last point here */
	    dx = iVtxs[i].vtX - iVtxs[i-1].vtX;
	    dy = iVtxs[i].vtY - iVtxs[i-1].vtY;
	    dist = sqrt(dx*dx + dy*dy);
	  }
	}
      }
      /* add last point if needed, reset wrap */
      if( del < (spacing/2.0) ){
	jVtxs[n] = iVtxs[nVtxs-1];
	n++;
      }
      for(i=1; i < wrap; i++, n++){
	jVtxs[n] = jVtxs[i];
      }
      rtnPoly->nvertices = n;
      break;

    case WLZ_POLYGON_FLOAT:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      fVtxs = (WlzFVertex2 *) poly->vtx;
      gVtxs = (WlzFVertex2 *) rtnPoly->vtx;
      gVtxs[0] = fVtxs[0];
      n = 1;
      i = 1;
      del = spacing;
      dx = fVtxs[i].vtX - fVtxs[i-1].vtX;
      dy = fVtxs[i].vtY - fVtxs[i-1].vtY;
      dist = sqrt(dx*dx + dy*dy);
      while( i < nVtxs ){
	if( del <= dist ){
	  /* new point in this span */
	  gVtxs[n].vtX = fVtxs[i-1].vtX + dx*del/dist;
	  gVtxs[n].vtY = fVtxs[i-1].vtY + dy*del/dist;
	  del += spacing;
	  n++;
	}
	else {
	  /* next span */
	  i++;
	  del -= dist;
	  if( i < nVtxs ){
	    /* should put in last point here */
	    dx = fVtxs[i].vtX - fVtxs[i-1].vtX;
	    dy = fVtxs[i].vtY - fVtxs[i-1].vtY;
	    dist = sqrt(dx*dx + dy*dy);
	  }
	}
      }
      /* add last point if needed, reset wrap */
      if( del < (spacing/2.0) ){
	gVtxs[n] = fVtxs[nVtxs-1];
	n++;
      }
      for(i=1; i < wrap; i++, n++){
	gVtxs[n] = gVtxs[i];
      }
      rtnPoly->nvertices = n;
      break;

    case WLZ_POLYGON_DOUBLE:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      dVtxs = (WlzDVertex2 *) poly->vtx;
      eVtxs = (WlzDVertex2 *) rtnPoly->vtx;
      eVtxs[0] = dVtxs[0];
      n = 1;
      i = 1;
      del = spacing;
      dx = dVtxs[i].vtX - dVtxs[i-1].vtX;
      dy = dVtxs[i].vtY - dVtxs[i-1].vtY;
      dist = sqrt(dx*dx + dy*dy);
      while( i < nVtxs ){
	if( del <= dist ){
	  /* new point in this span */
	  eVtxs[n].vtX = dVtxs[i-1].vtX + dx*del/dist;
	  eVtxs[n].vtY = dVtxs[i-1].vtY + dy*del/dist;
	  del += spacing;
	  n++;
	}
	else {
	  /* next span */
	  i++;
	  del -= dist;
	  if( i < nVtxs ){
	    /* should put in last point here */
	    dx = dVtxs[i].vtX - dVtxs[i-1].vtX;
	    dy = dVtxs[i].vtY - dVtxs[i-1].vtY;
	    dist = sqrt(dx*dx + dy*dy);
	  }
	}
      }
      /* add last point if needed, reset wrap */
      if( del < (spacing/2.0) ){
	eVtxs[n] = dVtxs[nVtxs-1];
	n++;
      }
      for(i=1; i < wrap; i++, n++){
	eVtxs[n] = dVtxs[i];
      }
      rtnPoly->nvertices = n;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnPoly;
}

