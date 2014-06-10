#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMwrAngle_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzMwrAngle.c
* \author       Jim Piper, Bill Hill
* \date         February 2002
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
* \brief	Computes the minimum width rectangle from a
* 		convex hull. This code has been extensively rewritten
* 		for the new convex hull domains, so that only the original
* 		interface and algorithm remain.
* \ingroup	WlzConvexHull
*/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

static double			WlzMwrAngleConvHull2I(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr);
static double			WlzMwrAngleConvHull2D(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr);

/*! 
* \return       Angle of the minimum width rectangles longest side with
* 	`	respect to the x-axis.
* \ingroup      WlzConvexHull
* \brief        Given a 2D convex hull object, computes the angle of the
* 		minimum width rectangle which encloses the convex hull.
* 		Where the angle is in radians and with respect to the x-axis.
* 		The minimum width rectangle long side must be parallel to an
* 		edge of convex hull and all sides must have at least one vertex
* 		of convex hull lying within them.
* \param    	cObj			Given object with a 2D convex hull
* 					domain.
* \param    	dstErr			Destination error pointer, may be NULL.
*/
double 				WlzMwrAngle(
				  WlzObject *cObj,
				  WlzErrorNum *dstErr)
{
  double	a = 0.0;
  WlzConvHullDomain2 *cvh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->type != WLZ_CONV_HULL)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((cvh = cObj->domain.cvh2) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cvh->type != WLZ_CONVHULL_DOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    switch(cvh->vtxType)
    {
      case WLZ_VERTEX_I2:
	a = WlzMwrAngleConvHull2I(cvh, &errNum);
	break;
      case WLZ_VERTEX_D2:
	a = WlzMwrAngleConvHull2D(cvh, &errNum);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_DATA;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(a);
}

/*!
* \return       Angle of the minimum width rectangles longest side with
* 	`	respect to the x-axis.
* \ingroup      WlzConvexHull
* \brief        Given a 2D convex hull object with integer vertices, computes
* 		the angle of the minimum width rectangle which encloses the
* 		convex hull (see WlzMwrAngle()).
* \param	cvh			Given 2D convex hull domain with
* 					integer vertices.
* \param    	dstErr			Destination error pointer, may be NULL.
*/
static double			WlzMwrAngleConvHull2I(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr)
{
  int		i,
  		iMin = -1;
  double	a = 0.0,
  		dMin = 0.0;
  WlzIVertex2 	p0,
    		p1;
  WlzDVertex2	eMin;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  p1 = cvh->vertices.i2[0];
  for(i = 0; i < cvh->nVertices; ++i)
  {
    int		j;

    /* Compute edge vector e, from current vertex of the convex hull to
     * the next. */
    j = (i + 1) % cvh->nVertices;
    p0 = p1;
    p1 = cvh->vertices.i2[j];
    if((p0.vtX - p1.vtX != 0) || (p0.vtY - p1.vtY != 0))
    {
      double	l,
      		dMax = 0.0;
      WlzDVertex2 e,
      		  p;

      WLZ_VTX_2_SUB(e, p1, p0);
      /* Compute normalised vector p perpendicular to edge e, knowing that
       * the length of edge e is non-zero. */
      l = 1.0 / (WLZ_VTX_2_LENGTH(e));
      WLZ_VTX_2_SET(p, -(e.vtY) * l, e.vtX * l);
      /* Find the vertex of the convex hull with the greatest perpendicular
       * distance from the edge e. */
      for(j = 0; j < cvh->nVertices; ++j)
      {
	double	d;
        WlzIVertex2 q;
        WlzDVertex2 f;

        q = cvh->vertices.i2[j];
	WLZ_VTX_2_SUB(f, q, p0);
	d = fabs(WLZ_VTX_2_DOT(p, f));
	if(d > dMax)
	{
	  dMax = d;
	}
      }
      if((iMin < 0) || (dMax < dMin))
      {
        dMin = dMax;
	eMin = e;
	iMin = i;
      }
    }
  }
  if(iMin < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    a = atan2(eMin.vtY, eMin.vtX);
  }
  *dstErr = errNum;
  return(a);
}

/*!
* \return       Angle of the minimum width rectangles longest side with
* 	`	respect to the x-axis.
* \ingroup      WlzConvexHull
* \brief        Given a 2D convex hull object with double vertices, computes
* 		the angle of the minimum width rectangle which encloses the
* 		convex hull (see WlzMwrAngle()).
* \param	cvh			Given 2D convex hull domain with
* 					double vertices.
* \param    	dstErr			Destination error pointer, may be NULL.
*/
static double			WlzMwrAngleConvHull2D(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr)
{
  int		i,
  		iMin = -1;
  double	a = 0.0,
  		dMin = 0.0;
  WlzDVertex2 	p0,
    		p1,
             	eMin;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-06;

  p1 = cvh->vertices.d2[0];
  for(i = 0; i < cvh->nVertices; ++i)
  {
    int		j;

    /* Compute edge vector e, from current vertex of the convex hull to
     * the next. */
    j = (i + 1) % cvh->nVertices;
    p0 = p1;
    p1 = cvh->vertices.d2[j];
    if((fabs(p0.vtX - p1.vtX) > eps) || (fabs(p0.vtY - p1.vtY) > eps))
    {
      double	l,
      		dMax = 0.0;
      WlzDVertex2 e,
      		  p;

      WLZ_VTX_2_SUB(e, p1, p0);
      /* Compute normalised vector p perpendicular to edge e, knowing that
       * the length of edge e is non-zero. */
      l = 1.0 / (WLZ_VTX_2_LENGTH(e));
      WLZ_VTX_2_SET(p, -(e.vtY) * l, e.vtX * l);
      /* Find the vertex of the convex hull with the greatest perpendicular
       * distance from the edge e. */
      for(j = 0; j < cvh->nVertices; ++j)
      {
	double	d;
        WlzDVertex2 q,
		    f;

        q = cvh->vertices.d2[j];
	WLZ_VTX_2_SUB(f, q, p0);
	d = fabs(WLZ_VTX_2_DOT(p, f));
	if(d > dMax)
	{
	  dMax = d;
	}
      }
      if((iMin < 0) || (dMin < dMax))
      {
        dMin = dMax;
	eMin = e;
	iMin = i;
      }
    }
  }
  if(iMin < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    a = atan2(eMin.vtY, eMin.vtX);
  }
  *dstErr = errNum;
  return(a);
}

