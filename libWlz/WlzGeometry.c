#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzGeometry.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Provides geometry utility functions.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static int		WlzGeomVtxSortRadialFn(
			  void *p0,
			  void *p1);

/************************************************************************
* Function:	WlzGeomTriangleCircumcentre			
* Returns:	int:			Zero if the circumcentre of the
*					triangle lies at infinity,
*					else non-zero.		
* Purpose:	Computes the circumcentre of the given triangle.
*		Given a triangle (a0, a1), (b0, b1), (c0, c1) then the
*		circumcentre (p0, 1) is given by:		
*		  p0 = (a0^2*b1 - a0^2*c1 - b1^2*a1 + c1^2*a1 +	
*		        b0^2*c1	+ a1^2*b1 + c0^2*a1 - c1^2*b1 -	
*		        c0^2*b1 - b0^2*a1 + b1^2*c1 - a1^2*c1) / D
*		  p1 = (a0^2*c0 + a1^2*c0 + b0^2*a0 - b0^2*c0 +	
*		        b1^2*a0 - b1^2*c0 - a0^2*b0 - a1^2*b0 -	
*		        c0^2*a0 + c0^2b0 - c1^2*a0 + c1^2*b0) / D
*		Where:						
*		  D = 2 * (a1*c0 + b1*a0 - b1*c0 - a1*b0 -	
*			   c1*a0 + c1*b0)			
*		This is taken from J. O'Rourke: Computational Geometry
*		in C, p201.					
* Global refs:	-						
* Parameters:	WlzDVertex2 *ccVx:	Destination ptr for the	
*					circumcentre.		
*		WlzDVertex2 vx0:	First vertex of triangle.
*		WlzDVertex2 vx1:	Second vertex of triangle.
*		WlzDVertex2 vx2:	Third vertex of triangle.
************************************************************************/
int		WlzGeomTriangleCircumcentre(WlzDVertex2 *ccVx,
					    WlzDVertex2 vx0,
					    WlzDVertex2 vx1,
					    WlzDVertex2 vx2)
{
  int		finite = 0;
  double	tD[13];

  tD[0]  = vx0.vtX * vx1.vtY;
  tD[1]  = vx0.vtX * vx2.vtY;
  tD[2]  = vx0.vtY * vx1.vtX;
  tD[3]  = vx0.vtY * vx2.vtX;
  tD[4]  = vx1.vtX * vx2.vtY;
  tD[5]  = vx1.vtY * vx2.vtX;
  tD[6]  = vx0.vtX * vx1.vtX;
  tD[7]  = vx0.vtX * vx2.vtX;
  tD[8]  = vx1.vtX * vx2.vtX;
  tD[9]  = vx0.vtY * vx1.vtY;
  tD[10] = vx0.vtY * vx2.vtY;
  tD[11] = vx1.vtY * vx2.vtY;
  tD[12] = tD[0] + tD[3] - tD[5] - tD[2] - tD[1] + tD[4];
  if((tD[12] * tD[12]) > DBL_EPSILON)
  {
    finite = 1;
    tD[12] = 0.5 / tD[12];
    ccVx->vtX = ((vx0.vtX * (tD[0] - tD[1])) +
		 (vx0.vtY * (tD[9] - tD[10])) +
		 (vx1.vtX * (tD[4] - tD[2])) +
		 (vx1.vtY * (tD[11] - tD[9])) +
		 (vx2.vtX * (tD[3] - tD[5])) +
		 (vx2.vtY * (tD[10] - tD[11]))) * tD[12];
    ccVx->vtY = ((vx0.vtX * (tD[7] - tD[6])) +
		 (vx0.vtY * (tD[3] - tD[2])) +
		 (vx1.vtX * (tD[6] - tD[8])) +
		 (vx1.vtY * (tD[0] - tD[5])) +
		 (vx2.vtX * (tD[8] - tD[7])) +
		 (vx2.vtY * (tD[4] - tD[1]))) * tD[12];
  }
  return(finite);
}

/************************************************************************
* Function:	WlzGeomVxInTriangle				
* Returns:	int:			Non-zero if vertex is inside
*					the triangle.		
* Purpose:	Test's to set if the given vertex lies within the given
*		triangle using a barycentric coordinates test.	
*		If a triangle has vertices p0, p1, p2, then any point
*		in the plane containing the triangle can be represented
*		by:						
* 		  p = alpha*p0 + beta*p2 + gamma*p3		
*		subject to the constraint:			
*		  alpha + beta + gamma = 1			
*		If p is inside the triangle at least one of alpha, beta
*		and gamma is -ve.				
* Global refs:	-						
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.
*		WlzDVertex2 vx1:	Second vertex of triangle.
*		WlzDVertex2 vx2:	Third vertex of triangle.
*		WlzDVertex2 vxP:	Given vertex.		
************************************************************************/
int		 WlzGeomVxInTriangle(WlzDVertex2 vx0, WlzDVertex2 vx1,
				     WlzDVertex2 vx2, WlzDVertex2 vxP)
{
  int		isInside = 0;
  double	tD0,
  		tD1,
		x0,
  		y0,
		x1,
		y1,
		x2,
		y2,
		alpha,
		beta,
		gamma;

  tD0 = vx2.vtX;
  tD1 = vx2.vtY;
  x0 = vxP.vtX - tD0;
  y0 = vxP.vtY - tD1;
  x1 = vx0.vtX - tD0;
  y1 = vx0.vtY - tD1;
  x2 = vx1.vtX - tD0;
  y2 = vx1.vtY - tD1;
  if((x2 * x2) > DBL_EPSILON)
  {
    isInside = (fabs(tD0 = (x1 * y2) - (x2 * y1)) > DBL_EPSILON) &&
	       ((alpha = ((x0 * y2) - (x2 * y0)) / tD0) >= 0.0) &&
	       (alpha <= 1.0) &&
	       ((beta = (x0 - (alpha * x1)) / x2) >= 0.0) &&
	       (beta <= 1.0) &&
	       ((gamma = (1 - (alpha + beta))) >= 0.0) &&
	       (gamma <= 1.0);
  }
  else if((y2 * y2)  > DBL_EPSILON)
  {
    isInside = (fabs(tD0 = (x1 * y2) - (x2 * y1)) > DBL_EPSILON) &&
	       ((alpha = ((x0 * y2) - (x2 * y0)) / tD0) >= 0.0) &&
	       (alpha <= 1.0) &&
	       ((beta = (y0 - (alpha * y1)) / y2) >= 0.0) &&
	       (beta <= 1.0) &&
	       ((gamma = (1 - (alpha + beta))) >= 0.0) &&
	       (gamma <= 1.0);
  }
  return(isInside);
}

/************************************************************************
* Function:	WlzGeomTriangleSnArea2				
* Returns:	double:			Twice the signed area of the
*					given triangle.		
* Purpose:	Computes twice the signed area of the given triangle.
*		The determinant is NOT computed with:
*		  (x0 - x1)(y1 - y2) - (y0 - y1)(x1 - x2)
*		instead the factorized form is used because it is
*		more robust numericaly.
* Global refs:	-						
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.
*		WlzDVertex2 vx1:	Second vertex of triangle.
*		WlzDVertex2 vx2:	Third vertex of triangle.
************************************************************************/
double		WlzGeomTriangleSnArea2(WlzDVertex2 vx0, WlzDVertex2 vx1,
				       WlzDVertex2 vx2)
{
  double	area2;

  area2 = (vx0.vtX * vx1.vtY) - (vx0.vtY * vx1.vtX) +
          (vx0.vtY * vx2.vtX) - (vx0.vtX * vx2.vtY) +
          (vx1.vtX * vx2.vtY) - (vx2.vtX * vx1.vtY);
  return(area2);
}

/************************************************************************
* Function:	WlzGeomInTriangleCircumcircle			
* Returns:	int:			Non zero if there's a conflict.
* Purpose:	Tests to see if the given vertex is inside the	
*		circumcircle of the given triangle.		
* Global refs:	-						
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.
*		WlzDVertex2 vx1:	Second vertex of triangle.
*		WlzDVertex2 vx2:	Third vertex of triangle.
*		WlzDVertex2 gVx:	Given vertex to test.	
************************************************************************/
int		WlzGeomInTriangleCircumcircle(WlzDVertex2 vx0, WlzDVertex2 vx1,
					      WlzDVertex2 vx2, WlzDVertex2 gVx)
{
  int			conflict;
  double		xp,
			yp,
			x0,
			y0,
			x1,
			y1,
			x2,
			y2,
			z1,
			z2,
			alpha,
			beta,
			gamma;

  x0 = vx0.vtX;
  y0 = vx0.vtY;
  x1 = vx1.vtX - x0;
  y1 = vx1.vtY - y0;
  x2 = vx2.vtX - x0;
  y2 = vx2.vtY - y0;
  xp = gVx.vtX - x0;
  yp = gVx.vtY - y0;
  z1 = (x1 * x1) + (y1 * y1);
  z2 = (x2 * x2) + (y2 * y2);
  alpha = (y1 * z2) - (z1 * y2);
  beta = (x2 * z1) - (x1 * z2);
  gamma = (x1 * y2) - (y1 * x2);
  conflict = ((alpha * xp) + (beta * yp) +
	      (gamma * ((xp * xp) + (yp * yp))) < 0.0);
  return(conflict);
}

/************************************************************************
* Function:	WlzGeomLineSegmentsIntersect			
* Returns:	int:			Non zero if line segments
*					intersect.
* Purpose:	Tests to see if the two given line segments intersect.
*		This is taken from J. O'Rourke: Computational Geometry
*		in C, p250.					
* Global refs:	-						
* Parameters:	WlzDVertex2 p0:		1st vertex of 1st line segment.
*		WlzDVertex2 p1:		2nd vertex 1st line segment.
*		WlzDVertex2 q0:		1st vertex of 2nd line segment.
*		WlzDVertex2 q1:		2nd vertex of 2nd line segment.
*		WlzDVertex2 *dstN	Destination ptr for intersection
*					vertex, may be NULL.
************************************************************************/
int		WlzGeomLineSegmentsIntersect(WlzDVertex2 p0, WlzDVertex2 p1,
					     WlzDVertex2 q0, WlzDVertex2 q1,
					     WlzDVertex2 *dstN)
{
  int		intersect = 0;
  double	sp,
  		tp,
		dnm;

  dnm = (p0.vtX * (q1.vtY - q0.vtY)) + (p1.vtX * (q0.vtY - q1.vtY)) +
        (q0.vtX * (p0.vtY - p1.vtY)) + (q1.vtX * (p1.vtY - p0.vtY));
  if(fabs(dnm) > DBL_EPSILON)
  {
    sp = ((p0.vtX * (q1.vtY - q0.vtY)) +
          (q0.vtX * (p0.vtY - q1.vtY)) +
	  (q1.vtX * (q0.vtY - p0.vtY))) / dnm;
    tp = -(p0.vtX * (q0.vtY - p1.vtY) +
           p1.vtX * (p0.vtY - q0.vtY) +
	   q0.vtX * (p1.vtY - p0.vtY)) / dnm;
  }
  if((sp > 0.0) && (sp <= 1.0) && (tp > 0.0) && (tp <= 1.0))
  {
    if(dstN)
    {
      dstN->vtX = p0.vtX + (sp * (p1.vtX - p0.vtX));
      dstN->vtY = p0.vtY + (sp * (p1.vtY - p0.vtY));
    }
  }
  return(intersect);
}

/************************************************************************
* Function:	WlzGeomCmpAngle
* Returns:	int:			Result of comparison: -ve, 0
*					or +ve. Only the sign is
*					meaningful.
* Purpose:	Given two end connected 2D line segments: (p0, O) and
*		(p1, O), compares the CCW angle of the segments,
*		where O is the origin (0,0).
* Global refs:	-
* Parameters:	WlzDVertex2 p0:		1st segment endpoint vertex.
*		WlzDVertex2 p1:		2nd segment endpoint vertex.
************************************************************************/
int		WlzGeomCmpAngle(WlzDVertex2 p0, WlzDVertex2 p1)
{
  int		q0,
  		q1,
		o0,
		o1,
		cmp = 0;
  double	tst = 0.0;
  WlzDVertex2	s0,
  		s1;
  const int	quadTbl[4] = {2, 3, 1, 0},
  		octTbl[8] = {4, 0, 1, 5, 6, 2, 3, 7};

  /* Compute relative endpoints. */
  /* Find quadrants:
   *
   *        ^ Y
   *        |
   *   q=1  | q=0
   * -------O--------> X
   *   q=2  | q=3
   *        |
   */
  q0 = quadTbl[((p0.vtY > 0.0) << 1) | (p0.vtX > 0.0)]; 
  q1 = quadTbl[((p1.vtY > 0.0) << 1) | (p1.vtX > 0.0)]; 
  /* Compare quadrants. */
  if((cmp = q0 - q1) == 0)
  {
    /* Find octants:
     *       ^ Y
     *       |
     *  \o=2 |o=1 /
     *    \  |  /
     *  o=3 \|/ o=0
     * ------O--------> X
     *  o=4 /|\ o=7
     *    /  |  \
     *  / 0=5|o=6 \
     *       |
     */
    s0.vtX = p0.vtX * p0.vtX; s0.vtY = p0.vtY * p0.vtY;
    s1.vtX = p1.vtX * p1.vtX; s1.vtY = p1.vtY * p1.vtY;
    o0 = octTbl[((s0.vtX > s0.vtY) << 2) | q0];
    o1 = octTbl[((s1.vtX > s1.vtY) << 2) | q1];
    /* Compare octants. */
    if((cmp = o0 - o1) == 0)
    {
      /* The octants are the same need separate case for each octant. */
      switch(o0)
      {
        case 0:
	  if((s0.vtX > DBL_EPSILON) && (s1.vtX > DBL_EPSILON))
	  {
	    tst = (s0.vtY / s0.vtX) - (s1.vtY / s1.vtX);
	  }
	  break;
        case 1:
	  if((s0.vtY > DBL_EPSILON) && (s1.vtY > DBL_EPSILON))
	  {
	    tst = (s1.vtX / s1.vtY) - (s0.vtX / s0.vtY);
	  }
	  break;
        case 2:
	  if((s0.vtY > DBL_EPSILON) && (s1.vtY > DBL_EPSILON))
	  {
	    tst = (s0.vtX / s0.vtY) - (s1.vtX / s1.vtY);
	  }
	  break;
        case 3:
	  if((s0.vtX > DBL_EPSILON) && (s1.vtX > DBL_EPSILON))
	  {
	    tst = (s1.vtY / s1.vtX) - (s0.vtY / s0.vtX);
	  }
	  break;
        case 4:
	  if((s0.vtX > DBL_EPSILON) && (s1.vtX > DBL_EPSILON))
	  {
	    tst = (s0.vtY / s0.vtX) - (s1.vtY / s1.vtX);
	  }
	  break;
        case 5:
	  if((s0.vtY > DBL_EPSILON) && (s1.vtY > DBL_EPSILON))
	  {
	    tst = (s1.vtX / s1.vtY) - (s0.vtX / s0.vtY);
	  }
	  break;
        case 6:
	  if((s0.vtY > DBL_EPSILON) && (s1.vtY > DBL_EPSILON))
	  {
	    tst = (s0.vtX / s0.vtY) - (s1.vtX / s1.vtY);
	  }
	  break;
        case 7:
	  if((s0.vtX > DBL_EPSILON) && (s1.vtX > DBL_EPSILON))
	  {
	    tst = (s1.vtY / s1.vtX) - (s0.vtY / s0.vtX);
	  }
	  break;
      }
      if(tst > DBL_EPSILON)
      {
        cmp = 1.0;
      }
      else if(tst < -(DBL_EPSILON))
      {
        cmp = -1.0;
      }
    }
  }
  return(cmp);
}

/************************************************************************
* Function:	WlzGeomVtxEqual2D
* Returns:	int:			1 if node positions are equal,
*					else 0.
* Purpose:	Checks to see if two verticies are the same
*		within some tollerance.
* Global refs:	-
* Parameters:	WlzDVertex2 pos0:	First node position.
*		WlzDVertex2 pos1:	Second node position.
*		double tolSq:		Square of tollerance value.
************************************************************************/
int		WlzGeomVtxEqual2D(WlzDVertex2 pos0, WlzDVertex2 pos1,
				  double tolSq)
{
  int		equal;

  pos0.vtX -= pos1.vtX;
  pos0.vtY -= pos1.vtY;
  equal = ((pos0.vtX * pos0.vtX) + (pos0.vtY * pos0.vtY)) < tolSq;
  return(equal);
}

/************************************************************************
* Function:	WlzGeomVtxSortRadialFn
* Returns:	int:			Result of comparison.
* Purpose:	Simple wrapper for WlzGeomCmpAngle().
* Global refs:	-
* Parameters:	void *p0:		Ptr to first vertex.
*		void *p1:		Ptr to second vertex.
************************************************************************/
static int	WlzGeomVtxSortRadialFn(void *p0, void *p1)
{
  return(WlzGeomCmpAngle(*(WlzDVertex2 *)p0, *(WlzDVertex2 *)p1));
}

/************************************************************************
* Function:	WlzGeomVtxSortRadial
* Returns:	void
* Purpose:	Sorts the given 3D verticies, which lie in a plane
*		perpendicular to the radial vector, in order of their
*		angle the radial vector.
*		No checks are made of the given parameters validity,
*		it's assumed that:
*		  (nV > 0) &&
*		  (vP != NULL) && (wP != NULL) && (iP != NULL)
*		  (|rV| > 0) && (rV.(uV = *vP) == 0)
* Global refs:	-
* Parameters:	int nV:			Number of 3D verticies.
*		WlzDVertex3 *vP:	The 3D verticies.
*		WlzDVertex2 *wP:	Workspace with nV 2D verticies.
*		WlzDVertex3 rV:		The radial vector.
************************************************************************/
void		WlzGeomVtxSortRadial(int nV, WlzDVertex3 *vP,
				     WlzDVertex2 *wP, WlzDVertex3 rV)
{
  int		idI;
  double	tD0;
  WlzDVertex3	uV,
  		vV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Compute basis vectors orthogonal to rV and in the plane of the
   * verticies. */
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(rV);
  WLZ_VTX_3_SCALE(rV, rV, tD0);
  uV = *vP;
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(uV);
  WLZ_VTX_3_SCALE(uV, uV, tD0);
  WLZ_VTX_3_CROSS(vV, rV, uV);
  /* TODO check I need to scale vV. */
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(vV);
  WLZ_VTX_3_SCALE(vV, vV, tD0);
  /* Compute the projections of the verticies onto the basis vectors. */
  for(idI = 0; idI < nV; ++idI)
  {
    (wP + idI)->vtX = WLZ_VTX_3_DOT(uV, *(vP + idI));
    (wP + idI)->vtY = WLZ_VTX_3_DOT(vV, *(vP + idI));
  }
  /* Sort the vericies. */
  (void )AlgHeapSort(wP, nV, sizeof(WlzDVertex2), WlzGeomVtxSortRadialFn);

}
