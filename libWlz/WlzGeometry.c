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


/************************************************************************
* Function:	WlzGeomTriangleCircumcentre				*
* Returns:	int:			Zero if the circumcentre of the	*
*					triangle lies at infinity,	*
*					else non-zero.			*
* Purpose:	Computes the circumcentre of the given triangle.	*
*		Given a triangle (a0, a1), (b0, b1), (c0, c1) then the	*
*		circumcentre (p0, 1) is given by:			*
*		  p0 = (a0^2*b1 - a0^2*c1 - b1^2*a1 + c1^2*a1 +		*
*		        b0^2*c1	+ a1^2*b1 + c0^2*a1 - c1^2*b1 -		*
*		        c0^2*b1 - b0^2*a1 + b1^2*c1 - a1^2*c1) / D	*
*		  p1 = (a0^2*c0 + a1^2*c0 + b0^2*a0 - b0^2*c0 +		*
*		        b1^2*a0 - b1^2*c0 - a0^2*b0 - a1^2*b0 -		*
*		        c0^2*a0 + c0^2b0 - c1^2*a0 + c1^2*b0) / D	*
*		Where:							*
*		  D = 2 * (a1*c0 + b1*a0 - b1*c0 - a1*b0 -		*
*			   c1*a0 + c1*b0)				*
*		This is taken from J. O'Rourke: Computational Geometry	*
*		in C, p201.						*
* Global refs:	-							*
* Parameters:	WlzDVertex2 *ccVx:	Destination ptr for the		*
*					circumcentre.			*
*		WlzDVertex2 vx0:	First vertex of triangle.	*
*		WlzDVertex2 vx1:	Second vertex of triangle.	*
*		WlzDVertex2 vx2:	Third vertex of triangle.	*
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
* Function:	WlzGeomVxInTriangle					*
* Returns:	int:			Non-zero if vertex is inside	*
*					the triangle.			*
* Purpose:	Test's to set if the given vertex lies within the given	*
*		triangle using a barycentric coordinates test.		*
*		If a triangle has vertices p0, p1, p2, then any point	*
*		in the plane containing the triangle can be represented	*
*		by:							*
* 		  p = alpha*p0 + beta*p2 + gamma*p3			*
*		subject to the constraint:				*
*		  alpha + beta + gamma = 1				*
*		If p is inside the triangle at least one of alpha, beta	*
*		and gamma is -ve.					*
* Global refs:	-							*
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.	*
*		WlzDVertex2 vx1:	Second vertex of triangle.	*
*		WlzDVertex2 vx2:	Third vertex of triangle.	*
*		WlzDVertex2 vxP:	Given vertex.			*
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
* Function:	WlzGeomTriangleSnArea2					*
* Returns:	double:			Twice the signed area of the	*
*					given triangle.			*
* Purpose:	Computes twice the signed area of the given triangle.	*
*		The determinant is NOT computed the factorized form:	*
*		(x0 - x1)(y1 - y2) - (y0 - y1)(x1 - x2) because it is	*
*		less numericaly robust, although faster.		*
* Global refs:	-							*
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.	*
*		WlzDVertex2 vx1:	Second vertex of triangle.	*
*		WlzDVertex2 vx2:	Third vertex of triangle.	*
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
* Function:	WlzGeomInTriangleCircumcircle				*
* Returns:	int:			Non zero if there's a conflict.	*
* Purpose:	Tests to see if the given vertex is inside the		*
*		circumcircle of the given triangle.			*
* Global refs:	-							*
* Parameters:	WlzDVertex2 vx0:	First vertex of triangle.	*
*		WlzDVertex2 vx1:	Second vertex of triangle.	*
*		WlzDVertex2 vx2:	Third vertex of triangle.	*
*		WlzDVertex2 gVx:	Given vertex to test.		*
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
