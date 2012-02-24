#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPolyDecimate_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzPolyDecimate.c
* \author       Richard Baldock
* \date         July 2001
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
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
* \brief	Functions to decimate polyline and boundary domains
* 		The functions remove vertices that are parts of straight
* 		lines as defined by a maximum distance.
* 		The algorithm starts at vertex 1, walks along the line
* 		until at least one vertex between the start and current
* 		position is more than max-dist from the straight line
* 		between vertex 1 and current. All vertices between
* 		position 1 and current-1 are removed and position 1
* 		is incremented (to what was current-1). The process
* 		is then repeated.
* \ingroup	WlzPolyline
*/

#include <stdio.h>
#include <Wlz.h>


/* function:    WlzIVtx2TriangleHeight */
/*! 
* \ingroup      WlzPolyline
* \brief        Calculate the height of a triangle from the last vertex
*		to the line defined by the first two. Uses formula that
*		can be derived from Faux and Pratt p57-65 or see my lab
*		notebook #1, p17. For three vertices \f$\mathbf{v_1}, 
*		\mathbf{v_2}, \mathbf{v_3},\f$ we define the vectors
*		\f$\mathbf{r} = \mathbf{v_2} - \mathbf{v_1}\f$,
*		\f$\mathbf{s} = \mathbf{v_3} - \mathbf{v_1}\f$ then the
*		the required height \f$h\f$ is given by:
*		\f[
*		h = \frac{|\mathbf{r}|^2|\mathbf{s}|^2 -
*			(\mathbf{r}\cdot\mathbf{s})^2}
*			{|\mathbf{r}|^2}
*		\f]
*
* \return       height of the triangle to vertex k from baseline
*		defined by vertices i, j
* \param    vtxs	vertex array
* \param    i	index of first vertex
* \param    j	index of second vertex
* \param    k	index of third vertex
* \par      Source:
*               WlzPolyDecimate.c
*/
double WlzIVtx2TriangleHeight(
  WlzIVertex2	*vtxs,
  int		i,
  int		j,
  int		k)
{
  double val=-1.0;
  double rx, ry;
  double sx, sy;
  double m1, m2, m3;

  rx = vtxs[j].vtX - vtxs[i].vtX;
  ry = vtxs[j].vtY - vtxs[i].vtY;
  sx = vtxs[k].vtX - vtxs[i].vtX;
  sy = vtxs[k].vtY - vtxs[i].vtY;

  m1 = rx*rx + ry*ry;
  m2 = sx*sx + sy*sy;
  m3 = rx*sx + ry*sy;
  if( m1 > 0.0 ){
    val = sqrt((m1*m2 - m3*m3) / m1);
  }

  return val;
}

static double WlzFVtx2TriangleHeight(
  WlzFVertex2	*vtxs,
  int		i,
  int		j,
  int		k)
{
  double val=-1.0;
  double rx, ry;
  double sx, sy;
  double m1, m2, m3;

  rx = vtxs[j].vtX - vtxs[i].vtX;
  ry = vtxs[j].vtY - vtxs[i].vtY;
  sx = vtxs[k].vtX - vtxs[i].vtX;
  sy = vtxs[k].vtY - vtxs[i].vtY;

  m1 = rx*rx + ry*ry;
  m2 = sx*sx + sy*sy;
  m3 = rx*sx + ry*sy;
  if( m1 > 0.0 ){
    val = sqrt((m1*m2 - m3*m3) / m1);
  }

  return val;
}

static double WlzDVtx2TriangleHeight(
  WlzDVertex2	*vtxs,
  int		i,
  int		j,
  int		k)
{
  double val=-1.0;
  double rx, ry;
  double sx, sy;
  double m1, m2, m3;

  rx = vtxs[j].vtX - vtxs[i].vtX;
  ry = vtxs[j].vtY - vtxs[i].vtY;
  sx = vtxs[k].vtX - vtxs[i].vtX;
  sy = vtxs[k].vtY - vtxs[i].vtY;

  m1 = rx*rx + ry*ry;
  m2 = sx*sx + sy*sy;
  m3 = rx*sx + ry*sy;
  if( m1 > 0.0 ){
    val = sqrt((m1*m2 - m3*m3) / m1);
  }

  return val;
}

/* function:     WlzPolyDecimate */
/*! 
* \ingroup    WlzPolyline
* \return       decimated polygon domain, NULL on error
* \brief        Decimate a polyline by removing vertices that are
*		within straight line segments as defined by a maximum
*		distance.
*		The algorithm starts at vertex 1, walks along the line	
*		until at least one vertex between the start and current	
*		position is more than maxDist from the straight line	
*		between vertex 1 and current. All vertices between	
*		position 1 and current-1 are removed and position 1 	
*		is incremented (to what was current-1). The process	
*		is then repeated.
* \param    poly	input polygon domain
* \param    wrap	wrap value of the input polyline, the returned
*			line will have wrap=1 if the input wrap >= 1
* \param    maxDist	distance parameter to test vertex removal
* \param    dstErr	error return
* \par      Source:
*               WlzPolyDecimate.c
*/
WlzPolygonDomain  *WlzPolyDecimate(
  WlzPolygonDomain	*poly,
  int			wrap,
  double		maxDist,
  WlzErrorNum		*dstErr)
{
  WlzPolygonDomain	*rtnPoly=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			i, j, k, n, nVtxs;
  WlzIVertex2		*iVtxs, *jVtxs;
  WlzFVertex2		*fVtxs, *gVtxs;
  WlzDVertex2		*dVtxs, *eVtxs;
  

  /* check object and input parameters */
  if( poly == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }

  if( maxDist < 0.0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  if( wrap < 0 ){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* create return polyline and start testing
     Note with this version the first point can never be removed
     Wrap is used to detect that the line is closed but is always
     reset to 1 (one) i.e. the last point set equal to the first.
  */
  if((errNum == WLZ_ERR_NONE) &&
     (rtnPoly = WlzMakePolygonDomain(poly->type, 0, NULL,
			       poly->nvertices, 1, &errNum))){
    switch( poly->type ){
    case WLZ_POLYGON_INT:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      n = 1;
      iVtxs = poly->vtx;
      jVtxs = rtnPoly->vtx;
      jVtxs[0] = iVtxs[0];
      i=0;
      j=1;
      while( j < nVtxs ){
	for(k=(i+1); k < j; k++){
	  if( WlzIVtx2TriangleHeight(iVtxs, i, j, k) > maxDist ){
	    /* set next vertex, reset i, j */
	    jVtxs[n] = iVtxs[j-1];
	    n++;
	    i=j-1;
	    break;
	  }
	}
	j++;
      }
      if( i < (nVtxs-1) ){
	jVtxs[n] = iVtxs[nVtxs-1];
	n++;
      }
      rtnPoly->nvertices = n;
      break;

    case WLZ_POLYGON_FLOAT:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      n = 1;
      fVtxs = (WlzFVertex2 *) poly->vtx;
      gVtxs = (WlzFVertex2 *) rtnPoly->vtx;
      gVtxs[0] = fVtxs[0];
      i=0;
      j=1;
      while( j < nVtxs ){
	for(k=(i+1); k < j; k++){
	  if( WlzFVtx2TriangleHeight(fVtxs, i, j, k) > maxDist ){
	    /* set next vertex, reset i, j */
	    gVtxs[n] = fVtxs[j-1];
	    n++;
	    i=j-1;
	    break;
	  }
	}
	j++;
      }
      if( i < (nVtxs-1) ){
	gVtxs[n] = fVtxs[nVtxs-1];
	n++;
      }
      rtnPoly->nvertices = n;
      break;

    case WLZ_POLYGON_DOUBLE:
      nVtxs = (wrap > 0) ? poly->nvertices - wrap + 1 : poly->nvertices;
      n = 1;
      dVtxs = (WlzDVertex2 *) poly->vtx;
      eVtxs = (WlzDVertex2 *) rtnPoly->vtx;
      eVtxs[0] = dVtxs[0];
      i=0;
      j=1;
      while( j < nVtxs ){
	for(k=(i+1); k < j; k++){
	  if( WlzDVtx2TriangleHeight(dVtxs, i, j, k) > maxDist ){
	    /* set next vertex, reset i, j */
	    eVtxs[n] = dVtxs[j-1];
	    n++;
	    i=j-1;
	    break;
	  }
	}
	j++;
      }
      if( i < (nVtxs-1) ){
	eVtxs[n] = dVtxs[nVtxs-1];
	n++;
      }
      rtnPoly->nvertices = n;
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


/* function:     WlzBoundDecimate */
/*! 
* \ingroup    WlzBoundary
* \return       decimated boundary list, NULL on error
* \brief        Decimate a boundary list using WlzPolyDecimate() on the
*		boundary polylines.
*
* \param    bound	input boundary list
* \param    maxDist	distance parameter to test for vertex removal
* \param    dstErr	error return
* \par      Source:
*               WlzPolyDecimate.c
*/
WlzBoundList *WlzBoundDecimate(
  WlzBoundList	*bound,
  double 	maxDist,
  WlzErrorNum	*dstErr)
{
  WlzBoundList	*rtnBound=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check object */
  if( bound == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((rtnBound = WlzMakeBoundList(bound->type, 1, NULL,
				       &errNum)) != NULL){
    if( (errNum == WLZ_ERR_NONE) && bound->next ){
      rtnBound->next = 
	WlzAssignBoundList(WlzBoundDecimate(bound->next, maxDist, &errNum),
			   NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->down ){
      rtnBound->down = 
	WlzAssignBoundList(WlzBoundDecimate(bound->down, maxDist, &errNum),
			   NULL);
    }

    if( (errNum == WLZ_ERR_NONE) && bound->poly ){
      rtnBound->poly =
	WlzAssignPolygonDomain(WlzPolyDecimate(bound->poly, bound->wrap,
					       maxDist, &errNum), NULL);
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

