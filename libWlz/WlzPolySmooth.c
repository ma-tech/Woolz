#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzPolySmooth_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzPolySmooth.c
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
* \brief	Smoothing operations for polylines.
* \ingroup	WlzPolyline
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <Wlz.h>


/*!
* \ingroup WlzPolyline
* \def MAX_POLYSMOOTH_ITERATIONS
* 	Maximum number of iterations of 1-2-1 smoothing in WlzPolySmooth
*/

#define MAX_POLYSMOOTH_ITERATIONS 10

/* function:     WlzPolySmooth   */
/*! 
* \ingroup WlzPolyline
* \return       smoothed polygon domain, NULL on error
* \brief        performs iterative 1-2-1 smoothing on the polyline
*		treating the x & y values independently. The smoothing
*		of an integer polyline is done by converting first to
*		double vertices then cnverting back at the end.
*
* \param    poly	input polygon domain, same type as the input
*			polyline but after passing through
*			WlzPolyEquispace
* \param    wrap	wrap parameter for the input polyline,
*			the new will have the same wrap
* \param    iterations	number of ierations of 1-2-1 smoothing
* \param    dstErr	error return
* \par      Source:
*               WlzPolySmooth.c
*/
WlzPolygonDomain *WlzPolySmooth(
  WlzPolygonDomain	*poly,
  int			wrap,
  int 			iterations,
  WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*rtnPoly=NULL, *tmp1Poly, *tmp2Poly;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzFVertex2		*fVtxs;
  WlzDVertex2		*dVtxs;
  int			n, i;
  double		x0, x1, y0, y1, x2, y2, factor;

  /* check object and parameters */
  if( poly == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( (iterations < 0) || (iterations > MAX_POLYSMOOTH_ITERATIONS) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  
  if( errNum == WLZ_ERR_NONE ){
    switch( poly->type ){
    case WLZ_POLYGON_INT:
      /* smooth using double vertices */
      if((tmp1Poly = WlzConvertPolyType(poly, WLZ_POLYGON_DOUBLE,
                                        &errNum)) != NULL){
	if((tmp2Poly = WlzPolySmooth(tmp1Poly, wrap, iterations,
	                             &errNum)) != NULL){
	  rtnPoly = WlzConvertPolyType(tmp2Poly, WLZ_POLYGON_INT, &errNum);
	  WlzFreePolyDmn(tmp2Poly);
	}
	WlzFreePolyDmn(tmp1Poly);
      }
      break;

    case WLZ_POLYGON_FLOAT:
      if((rtnPoly = WlzPolyEquispace(poly, wrap, 1.0, 0, &errNum)) != NULL){

	/* 1-2-1 smoothing, iterated, leaving end-points alone unless wrapped */
	for(n=0, factor=1.0; n < iterations; n++){
	  fVtxs = (WlzFVertex2 *) rtnPoly->vtx;

	  /* start values for smoothing */
	  if( wrap ){
	    x1 = (fVtxs+rtnPoly->nvertices-wrap-1)->vtX;
	    y1 = (fVtxs+rtnPoly->nvertices-wrap-1)->vtY;
	    x2 = fVtxs->vtX;
	    y2 = fVtxs->vtY;
	  }
	  else {
	    x1 = 2*fVtxs->vtX - (fVtxs+1)->vtX;
	    y1 = 2*fVtxs->vtY - (fVtxs+1)->vtY;
	  }

	  /* loop though polyline */
	  for (i=0; i < (rtnPoly->nvertices - wrap - 1); i++) {
	    x0 = fVtxs->vtX;
	    fVtxs->vtX = (float )(x1 + 2 * x0 + (fVtxs+1)->vtX);
	    x1 = x0;
	    y0 = fVtxs->vtY;
	    fVtxs->vtY = (float )(y1 + 2 * y0 + (fVtxs+1)->vtY);
	    y1 = y0;
	    fVtxs++;
	  }

	  /* last vertex, take care of wrap */
	  if( wrap ){
	    x0 = fVtxs->vtX;
	    fVtxs->vtX = (float )(x1 + 2 * x0 + x2);
	    y0 = fVtxs->vtY;
	    fVtxs->vtY = (float )(y1 + 2 * y0 + y2);
	    fVtxs++;
	    i++;
	    while( i < rtnPoly->nvertices ){
	      *fVtxs = *(fVtxs - rtnPoly->nvertices + wrap);
	      fVtxs++;
	      i++;
	    }
	  }
	  else {
	    fVtxs->vtX *= 4.0;
	    fVtxs->vtY *= 4.0;
	  }
	  factor *= 4.0;
	}

	/* re-normalise */
	fVtxs = (WlzFVertex2 *) rtnPoly->vtx;
	for(i=0; i < rtnPoly->nvertices; i++, fVtxs++) {
	  fVtxs->vtX /= factor;
	  fVtxs->vtY /= factor;
	}
      }
      break;

    case WLZ_POLYGON_DOUBLE:
      if((rtnPoly = WlzPolyEquispace(poly, wrap, 1.0, 0, &errNum)) != NULL){

	/* 1-2-1 smoothing, iterated, leaving end-points alone unless wrapped */
	for(n=0, factor=1.0; n < iterations; n++){
	  dVtxs = (WlzDVertex2 *) rtnPoly->vtx;

	  /* start values for smoothing */
	  if( wrap ){
	    x1 = (dVtxs+rtnPoly->nvertices-wrap-1)->vtX;
	    y1 = (dVtxs+rtnPoly->nvertices-wrap-1)->vtY;
	    x2 = dVtxs->vtX;
	    y2 = dVtxs->vtY;
	  }
	  else {
	    x1 = 2*dVtxs->vtX - (dVtxs+1)->vtX;
	    y1 = 2*dVtxs->vtY - (dVtxs+1)->vtY;
	  }

	  /* loop though polyline */
	  for (i=0; i < (rtnPoly->nvertices - wrap - 1); i++) {
	    x0 = dVtxs->vtX;
	    dVtxs->vtX = x1 + 2 * x0 + (dVtxs+1)->vtX;
	    x1 = x0;
	    y0 = dVtxs->vtY;
	    dVtxs->vtY = y1 + 2 * y0 + (dVtxs+1)->vtY;
	    y1 = y0;
	    dVtxs++;
	  }

	  /* last vertex, take care of wrap */
	  if( wrap ){
	    x0 = dVtxs->vtX;
	    dVtxs->vtX = x1 + 2 * x0 + x2;
	    y0 = dVtxs->vtY;
	    dVtxs->vtY = y1 + 2 * y0 + y2;
	    dVtxs++;
	    i++;
	    while( i < rtnPoly->nvertices ){
	      *dVtxs = *(dVtxs - rtnPoly->nvertices + wrap);
	      dVtxs++;
	      i++;
	    }
	  }
	  else {
	    dVtxs->vtX *= 4.0;
	    dVtxs->vtY *= 4.0;
	  }
	  factor *= 4.0;
	}

	/* re-normalise */
	dVtxs = (WlzDVertex2 *) rtnPoly->vtx;
	for(i=0; i < rtnPoly->nvertices; i++, dVtxs++) {
	  dVtxs->vtX /= factor;
	  dVtxs->vtY /= factor;
	}
      }
      break;

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnPoly;
}

/* function:     WlzBoundSmooth   */
/*! 
* \ingroup WlzBoundary
* \return       new boundlist, NULL on error
* \brief        Smooth a boundary list using WlzPolySmooth().
*
* \param    bound	input boundary list
* \param    iterations	number of iterations for smoothing the polylines
* \param    dstErr	error return
* \par      Source:
*               WlzPolySmooth.c
*/
WlzBoundList *WlzBoundSmooth(
  WlzBoundList	*bound,
  int 		iterations,
  WlzErrorNum	*dstErr)
{
  WlzBoundList	*rtnBound=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object */
  if( bound == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((rtnBound = WlzMakeBoundList(bound->type, bound->wrap, NULL,
				       &errNum)) != NULL){
    if( (errNum == WLZ_ERR_NONE) && bound->next ){
      rtnBound->next = 
	WlzAssignBoundList(WlzBoundSmooth(bound->next, iterations, &errNum),
			   NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->down ){
      rtnBound->down = 
	WlzAssignBoundList(WlzBoundSmooth(bound->down, iterations, &errNum),
			   NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->poly ){
      rtnBound->poly =
	WlzAssignPolygonDomain(WlzPolySmooth(bound->poly, bound->wrap,
					     iterations, &errNum), NULL);
    }
    if( errNum != WLZ_ERR_NONE ){
      WlzFreeBoundList(rtnBound);
      rtnBound = NULL;
    }
  }

  if(dstErr){
    *dstErr = errNum;
  }
  return rtnBound;
}
