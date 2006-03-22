#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGeometry_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGeometry.c
* \author       Bill Hill
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
* \brief	Geometric utility functions.
* \ingroup	WlzGeometry
* \todo         -
* \bug          None known.
*/

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>

static int		WlzGeomVtxSortRadialFn(
			  void *p0,
			  int *idxP,
			  int idx0,
			  int idx1);

/*!
* \return	Zero if the circumcentre of the triangle lies at infinity,
*		else non-zero.
* \ingroup	WlzGeometry
* \brief	Computes the circumcentre of the given triangle.
*
*		Given a triangle \f$(a_0, a_1), (b_0, b_1), (c_0, c_1)\f$
*		then the circumcentre \f$(p_0, p_1)\f$ is given by:		
*		  \f[
		  p0 = (a_0^2 b_1 - a_0^2 c_1 - b_1^2 a_1 + c_1^2 a_1 +	
 		        b_0^2 c_1	+ a_1^2 b_1 + c_0^2 a_1 - c_1^2 b_1 -	
 		        c_0^2 b_1 - b_0^2 a_1 + b_1^2 c_1 - a_1^2 c_1) / D
*		  \f]
*		  \f[
 		  p1 = (a_0^2 c_0 + a_1^2 c_0 + b_0^2 a_0 - b_0^2 c_0 +	
 		        b_1^2 a_0 - b_1^2 c_0 - a_0^2 b_0 - a_1^2 b_0 -	
 		        c_0^2 a_0 + c_0^2 b_0 - c_1^2 a_0 + c_1^2 b_0) / D
*		  \f]
*		Where:						
*		  \f[
 		  D = 2 (a_1 c_0 + b_1 a_0 - b_1 c_0 - a_1 b_0 -	
 			 c_1 a_0 + c_1 b_0)			
*		  \f]
*		This is taken from J. O'Rourke: Computational Geometry
*		in C, p201.					
* \param	ccVx			Destination pointer for the
*					circumcentre.
* \param	vx0			First vertex of triangle.
* \param	vx1			Second vertex of triangle.
* \param	vx2			Third vertex of triangle.
*/
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

/*!
* \return	Value indicating the position of the vertex with respect
*               to the triangle:
*		  +ve if the vertex is inside the triangle,
*		  0   if the vertex is on an edge of the triangle and
*		  -ve if the vertex is outside the triangle.
* \ingroup	WlzGeometry
* \brief	Test's to set if the given vertex lies within the given
*		triangle using a barycentric coordinates test.
*
*		If a triangle has vertices \f$p_0, p_1, p_2\f$, then any point
*		in the plane containing the triangle can be represented
*		by: \f[p = \alpha*p_0 + \beta*p_2 + \gamma*p_3\f]
*		subject to the constraint: \f[\alpha + \beta + \gamma = 1\f]
*		If \f$p\f$ is outside the triangle at one or more of 
*		\f$alpha\f$, \f$beta\f$ and \f$gamma\f$ is -ve. It
*               is inside if all are -ve and on an edge of the triangle
*		if any are close to zero (ie < DBL_EPSILON).
* \param	vx0			First vertex of triangle.
* \param	vx1			Second vertex of triangle.
* \param	vx2			Third vertex of triangle.
* \param	vxP			Given vertex.
*/
int		 WlzGeomVxInTriangle(WlzDVertex2 vx0, WlzDVertex2 vx1,
				     WlzDVertex2 vx2, WlzDVertex2 vxP)
{
  int		inside = 0;
  double	tA,
  		tB,
		tC,
		tD,
		tE,
		tF,
		alpha,
		beta,
		gamma,
		delta;

  tA = vx0.vtX - vx2.vtX;
  tB = vx1.vtX - vx2.vtX;
  tD = vx0.vtY - vx2.vtY;
  tE = vx1.vtY - vx2.vtY;
  delta = (tA * tE) - (tB * tD);
  if(fabs(delta) > DBL_EPSILON)
  {
    tC = vx2.vtX - vxP.vtX;
    tF = vx2.vtY - vxP.vtY;
    alpha = ((tB * tF) - (tC * tE)) / delta;
    beta  = ((tC * tD) - (tA * tF)) / delta;
    gamma = 1.0 - (alpha + beta);
    if((alpha < -DBL_EPSILON) || (beta < -DBL_EPSILON) ||
       (gamma < -DBL_EPSILON))
    {
      inside = -1;
    }
    else if((alpha > DBL_EPSILON) && (beta > DBL_EPSILON) &&
            (gamma > DBL_EPSILON))
    {
      inside = 1;
    }
  }
  return(inside);

}

/*!
* \return	Twice the signed area of the given triangle.
* \ingroup	WlzGeometry
* \brief	Computes twice the signed area of the given triangle.
*
*		Computes twice the signed area of the given triangle.
*		The determinant is NOT computed with:
*		\f[(x_0 - x_1)(y_1 - y_2) - (y_0 - y_1)(x_1 - x_2)\f]
*		instead the factorized form is used because it is more robust
*		numericaly.
* \param	vx0			First vertex of triangle.
* \param	vx1			Second vertex of triangle.
* \param	vx2			Third vertex of triangle.
*/
double		WlzGeomTriangleSnArea2(WlzDVertex2 vx0, WlzDVertex2 vx1,
				       WlzDVertex2 vx2)
{
  double	area2;

  area2 = (vx0.vtX * vx1.vtY) - (vx0.vtY * vx1.vtX) +
          (vx0.vtY * vx2.vtX) - (vx0.vtX * vx2.vtY) +
          (vx1.vtX * vx2.vtY) - (vx2.vtX * vx1.vtY);
  return(area2);
}

/*!
* \return	Twice the square of the area of the given triangle.
* \ingroup	WlzGeometry
* \brief	Computes twice the square of the area of the given
*		3D triangle.
*
*		A nieve approach is used in which the area \f$A\f$ is
*		computed using:
*		\f[
		2 A^2 = \left|\left| \: \left|
			\begin{array}{ccc}
			\mathbf{i} & \mathbf{j} & \mathbf{k} \\
			a_x & a_y & a_z \\
			b_x & b_y & b_z
			\end{array}
		        \right| \: \right|\right|^2
		\f]
*		Where \f$\mathbf{a} = \mathbf{v_0} - \mathbf{v_1}\f$ and
*		\f$\mathbf{b} = \mathbf{v_2} - \mathbf{v_1}\f$.
* \param	vx0			First vertex of triangle
* 					\f$\mathbf{v_0}\f$.
* \param	vx1			Second vertex of triangle
*					\f$\mathbf{v_1}\f$.
* \param	vx2			Third vertex of triangle
*					\f$\mathbf{v_2}\f$.
*/
double		WlzGeomTriangleArea2Sq3(WlzDVertex3 vx0, WlzDVertex3 vx1,
				        WlzDVertex3 vx2)
{
  WlzDVertex3	a,
  		b,
		t;
  double	area3;

  WLZ_VTX_3_SUB(a, vx1, vx0);
  WLZ_VTX_3_SUB(b, vx2, vx0);
  WLZ_VTX_3_CROSS(t, a, b);
  area3 = WLZ_VTX_3_SQRLEN(t);
  return(area3);
}

/*!
* \return	Non zero if the given vertex is inside the circumcircle of
*		the given triangle.
* \ingroup      WlzGeometry
* \brief	Tests to see if the given vertex is inside the circumcircle of
* 		the given triangle.
* \param	vx0			First vertex of triangle.
* \param	vx1			Second vertex of triangle.
* \param	vx2			Third vertex of triangle.
* \param	gVx			Given vertex to test.
*/
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

/*!
* \return	Integer value which classifies the intersection of
*		the line segments, with values:
*		<ul>
*		  <li>0 no intersection.</li>
*		  <li>1 intersection along line segments or all points
*			are coincident.</li>
*                 <li>2 intersection at end points.</li>
*                 <li>3 intersection at a single point (not end points).</li>
*		</ul>
* \ingroup	WlzGeometry
* \brief	Tests to see if the two given line segments intersect.
*
*		Tests to see if the two given line segments intersect
*		using the DBL_EPSILON tollerance value.
*               This is taken from J. O'Rourke: Computational Geometry
*               in C, p250, but has ben modified to include the use of
*		DBL_EPSILON.
* \param	p0			1st vertex of 1st line segment.
* \param	p1			2nd vertex 1st line segment.
* \param	q0			1st vertex of 2nd line segment.
* \param	q1			2nd vertex of 2nd line segment.
* \param	dstN			Destination ptr for intersection
*					vertex, may be NULL. The intersection
*					value is not set if there is no
*					intersection or the intersection is
*					along the line segmants.
*/
int		WlzGeomLineSegmentsIntersect(WlzDVertex2 p0, WlzDVertex2 p1,
					     WlzDVertex2 q0, WlzDVertex2 q1,
					     WlzDVertex2 *dstN)
{
  int		ic,
  		intersect = 0;
  double	sgn,
		spd,
  		dna,
		dn,
  		sp,
  		tp;
  WlzDVertex2	ict;

  /* To minimize numerical problems don't factorize the expresion for dn. */
  dn = (p0.vtX * (q1.vtY - q0.vtY)) + (p1.vtX * (q0.vtY - q1.vtY)) +
       (q0.vtX * (p0.vtY - p1.vtY)) + (q1.vtX * (p1.vtY - p0.vtY));
  sgn = (dn < 0)? -1.0: 1.0;
  sp = ((p0.vtX * (q1.vtY - q0.vtY)) +
	(q0.vtX * (p0.vtY - q1.vtY)) +
	(q1.vtX * (q0.vtY - p0.vtY))) * sgn;
  tp = ((p0.vtX * (q0.vtY - p1.vtY)) +
	(p1.vtX * (p0.vtY - q0.vtY)) +
	(q0.vtX * (p1.vtY - p0.vtY))) * sgn * -1.0;
  if(fabs(dn) < DBL_EPSILON)
  {
    /* Line segments are parallel. */
    if((fabs(sp) < DBL_EPSILON) && (fabs(tp) < DBL_EPSILON))
    {
      /* Line segments are coincident. */
      ic = 0;
      if((fabs(p0.vtX - q0.vtX) < DBL_EPSILON) &&
	 (fabs(p0.vtY - q0.vtY) < DBL_EPSILON))
      {
        ict = p0;
	++ic;
      }
      if((fabs(p0.vtX - q1.vtX) < DBL_EPSILON) &&
	 (fabs(p0.vtY - q1.vtY) < DBL_EPSILON))
      {
        ict = p0;
	++ic;
      }
      if((fabs(p1.vtX - q0.vtX) < DBL_EPSILON) &&
	 (fabs(p1.vtY - q0.vtY) < DBL_EPSILON))
      {
        ict = p1;
	++ic;
      }
      if((fabs(p1.vtX - q1.vtX) < DBL_EPSILON) &&
	 (fabs(p1.vtY - q1.vtY) < DBL_EPSILON))
      {
        ict = p1;
	++ic;
      }
      if(ic != 1)
      {
        intersect = 1;
      }
      else
      {
        intersect = 2;
	if(dstN != NULL)
	{
	  *dstN = ict;
	}
      }
    }
  }
  else
  {
    /* Line segments are not parallel. */
    dna = dn * sgn;
    if((sp >= -(DBL_EPSILON)) && (sp < dna + DBL_EPSILON) &&
       (tp >= -(DBL_EPSILON)) && (tp < dna + DBL_EPSILON))
    {
      /* Line segments intersect. */
      intersect = ((sp > DBL_EPSILON) && (sp < (dna - DBL_EPSILON)) &&
		   (tp > DBL_EPSILON) && (tp < (dna - DBL_EPSILON)))? 3: 2;
      if(dstN != NULL)
      {
	spd = sp / dna;
	dstN->vtX = p0.vtX + (spd * (p1.vtX - p0.vtX));
	dstN->vtY = p0.vtY + (spd * (p1.vtY - p0.vtY));
      }
    }
  }
  return(intersect);
}

#ifdef WLZ_GEOM_LINESEGMENTSINTERSECT_MAIN
int		main(int argc, char *argv[])
{
  int		ic;
  WlzDVertex2	is,
  		p0,
  		p1,
		q0,
		q1;
  char		buf[1000];
  const int	bufSz = 1000;

  while((fgets(buf, bufSz, stdin) != NULL) &&
        (sscanf(buf, "%lg %lg %lg %lg %lg %lg %lg %lg",
	        &(p0.vtX), &(p0.vtY),
	        &(p1.vtX), &(p1.vtY),
	        &(q0.vtX), &(q0.vtY),
	        &(q1.vtX), &(q1.vtY)) == 8))
  {
    (void )printf("((%g, %g), (%g, %g)), ((%g, %g), (%g, %g)) -> ",
                  p0.vtX, p0.vtY, p1.vtX, p1.vtY,
                  q0.vtX, q0.vtY, q1.vtX, q1.vtY);
    ic = WlzGeomLineSegmentsIntersect(p0, p1, q0, q1, &is);
    (void )printf("%d", ic);
    if(ic > 1)
    {
      (void )printf(" (%g, %g)", is.vtX, is.vtY);
    }
    (void )printf("\n");
  }
  exit(0);
}

#endif /* WLZ_GEOM_LINESEGMENTSINTERSECT_MAIN */

/*!
* \return	Result of comparison: -ve, 0 or +ve. Only the sign is
*               meaningful.
* \ingroup	WlzGeometry
* \brief	Given two end connected 2D line segments this function
*		compares the CCW angle of the segments.
*
*		Given two end connected 2D line segments: \f$(p_0, O)\f$ and
*               \f$(p_1, O)\f$, compares the CCW angle of the segments,
*               where \f$O\f$ is the origin \f$(0, 0)\f$.
* \param	p0			1st segment endpoint vertex.
* \param	p1			2nd segment endpoint vertex.
*/
int		WlzGeomCmpAngle(WlzDVertex2 p0, WlzDVertex2 p1)
{
  int		i0,
  		i1,
		q0,
  		q1,
		o0,
		o1,
		cmp = 0;
  double	tst = 0.0;
  WlzDVertex2	s0,
  		s1;
  const int	quadTbl[4] = {2, 3, 1, 0},
  		octTbl[8] = {5, 6, 2, 1, 4, 7, 3, 0};

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
  i0 = ((p0.vtY > 0.0) << 1) | (p0.vtX > 0.0);
  i1 = ((p1.vtY > 0.0) << 1) | (p1.vtX > 0.0);
  q0 = quadTbl[i0]; 
  q1 = quadTbl[i1]; 
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
    i0 |= ((s0.vtX > s0.vtY) << 2);
    i1 |= ((s1.vtX > s1.vtY) << 2);
    o0 = octTbl[i0];
    o1 = octTbl[i1];
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

/*!
* \return	1 if node positions are equal, else 0.
* \ingroup	WlzGeometry
* \brief	Checks to see if two verticies are the same
*               within some tollerance.
* \param	pos0			First node position.
* \param	pos1			Second node position.
* \param	tolSq			Square of tollerance value.
*/
int		WlzGeomVtxEqual2D(WlzDVertex2 pos0, WlzDVertex2 pos1,
				  double tolSq)
{
  int		equal;

  pos0.vtX -= pos1.vtX;
  pos0.vtY -= pos1.vtY;
  equal = ((pos0.vtX * pos0.vtX) + (pos0.vtY * pos0.vtY)) < tolSq;
  return(equal);
}

/*!
* \return	Result of comparison.
* \ingroup	WlzGeometry
* \brief	Simple wrapper for WlzGeomCmpAngle().
* \param	p0			Ptr to verticies.
* \param	idxP			Ptr to indicies of verticies.
* \param	idx0			Index to index of first vertex.
* \param	idx1			Index to index of second vertex.
*/
static int	WlzGeomVtxSortRadialFn(void *p0, int *idxP, int idx0, int idx1)
{
  return(WlzGeomCmpAngle(*((WlzDVertex2 *)p0 + *(idxP + idx0)),
  			 *((WlzDVertex2 *)p0 + *(idxP + idx1))));
}

/*!
* \return	void
* \ingroup	WlzGeometry
* \brief	Sorts the given 3D verticies, which lie in a plane
*               perpendicular to the radial vector, in order of their
*               angle the radial vector.
*
*		Sorts the given 3D verticies, which lie in a plane
*               perpendicular to the radial vector, in order of their
*               angle the radial vector.
*               No checks are made of the given parameters validity,
*               it's assumed that:
*                 (nV > 0) &&
*                 (vP != NULL) && (wP != NULL) && (iP != NULL)
*                 (|rV| > 0) && (rV.(uV = *vP) == 0)
*               Note that it is the indicies that are sorted NOT the
*               verticies themselves.
* \param	nV			Number of 3D verticies.
* \param	vP			The 3D verticies.
* \param	idxBuf			Buffer of nV indicies used
*                                       for sorting the verticies.
* \param	wP			Workspace with nV 2D verticies.
* \param	rV			The radial vector.
*/
void		WlzGeomVtxSortRadial(int nV, WlzDVertex3 *vP,
				     int *idxBuf, WlzDVertex2 *wP,
				     WlzDVertex3 rV)
{
  int		idI;
  double	tD0;
  WlzDVertex3	uV,
  		vV;

  /* Compute basis vectors orthogonal to rV and in the plane of the
   * verticies. */
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(rV);
  WLZ_VTX_3_SCALE(rV, rV, tD0);
  uV = *vP;
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(uV);
  WLZ_VTX_3_SCALE(uV, uV, tD0);
  WLZ_VTX_3_CROSS(vV, rV, uV);
  tD0 = 1.0 / WLZ_VTX_3_LENGTH(vV);
  WLZ_VTX_3_SCALE(vV, vV, tD0);
  /* Compute the projections of the verticies onto the basis vectors. */
  for(idI = 0; idI < nV; ++idI)
  {
    *(idxBuf + idI) = idI;
    (wP + idI)->vtX = WLZ_VTX_3_DOT(uV, *(vP + idI));
    (wP + idI)->vtY = WLZ_VTX_3_DOT(vV, *(vP + idI));
  }
  /* Sort the vericies. */
  (void )AlgHeapSortIdx(wP, idxBuf, nV, WlzGeomVtxSortRadialFn);
}

/*!
* \return	Normal vector.
* \ingroup	WlzGeometry
* \brief	Computes the unit normal vector perpendicular to the
*               triangle \f$v_0, v_1, v_2\f$.
* \param	v0			First vertex of triangle.
* \param	v1			Second vertex of triangle.
* \param	v2			Third vertex of triangle.
*/
WlzDVertex3	WlzGeomTriangleNormal(WlzDVertex3 v0, WlzDVertex3 v1,
				      WlzDVertex3 v2)
{
  double	len;
  WlzDVertex3	nrm;

  WLZ_VTX_3_SUB(v1, v1, v0); 
  WLZ_VTX_3_SUB(v2, v2, v0); 
  WLZ_VTX_3_CROSS(nrm, v1, v2);
  len = WLZ_VTX_3_LENGTH(nrm);
  if(len > DBL_EPSILON)
  {
    WLZ_VTX_3_SCALE(nrm, nrm, 1.0 / len);
  }
  else
  {
    nrm.vtX = 0.0;
    nrm.vtY = 0.0;
    nrm.vtZ = 0.0;
  }
  return(nrm);
}

/*!
* \return	Non-zero if there is an intersection.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the plane defined by the
*		equation: \f$ax + by + cz + d = 0\f$ and the given axis
*		aligned bounding box.
* \param	a			Plane X parameter.
* \param	b			Plane Y parameter.
* \param	c			Plane Z parameter.
* \param	d			Other plane parameter.
* \param	box			Axis aligned bounding box.
*/
int		WlzGeomPlaneAABBIntersect(double a, double b,
					  double c, double d,
					  WlzDBox3 box)
{

  int           idI,
                maxP,
                intersect = 0;
  double        iVal;
  double        aP[3];
  WlzDVertex3   bV[4];

  /* Check for approximate direction of the plane. */
  aP[0] = fabs(a); aP[1] = fabs(b); aP[2] = fabs(c);
  if(aP[0] > aP[1])
  {
    maxP = (aP[0] > aP[2])? 0: 2;
  }
  else
  {
    maxP = (aP[1] > aP[2])? 1: 2;
  }
  if(aP[maxP] > DBL_EPSILON)
  {
    /* Get the verticies of the bounding box and check for an intersection
     * between four edges (not parallel to the plane) and the plane. */
    switch(maxP)
    {
      case 0:
	idI = 0;
        bV[0].vtX = -(box.xMin * a);     /* -ve implies inverted comparison. */
        bV[1].vtX = -(box.xMax * a);
        bV[0].vtY = bV[3].vtY = box.yMin;
        bV[0].vtZ = bV[1].vtZ = box.zMin;
        bV[1].vtY = bV[2].vtY = box.yMax;
        bV[2].vtZ = bV[3].vtZ = box.zMax;
        do
        {
          iVal = (b * bV[idI].vtY) + (c * bV[idI].vtZ) + d;
          intersect = (iVal <= bV[0].vtX) && (iVal >= bV[1].vtX);
        }
        while((intersect == 0) && (++idI < 4));
        break;
      case 1:
	idI = 0;
        bV[0].vtY = -(box.yMin * b);
        bV[1].vtY = -(box.yMax * b);
        bV[0].vtZ = bV[3].vtZ = box.zMin;
        bV[0].vtX = bV[1].vtX = box.xMin;
        bV[1].vtZ = bV[2].vtZ = box.zMax;
        bV[2].vtX = bV[3].vtX = box.xMax;
        do
        {
          iVal = (c * bV[idI].vtZ) + (a * bV[idI].vtX) + d;
          intersect = (iVal <= bV[0].vtY) && (iVal >= bV[1].vtY);
        }
        while((intersect == 0) && (++idI < 4));
        break;
      case 2:
	idI = 0;
        bV[0].vtZ = -(box.zMin * c);
        bV[1].vtZ = -(box.zMax * c);
        bV[0].vtX = bV[3].vtX = box.xMin;
        bV[0].vtY = bV[1].vtY = box.yMin;
        bV[1].vtX = bV[2].vtX = box.xMax;
        bV[2].vtY = bV[3].vtY = box.yMax;
        do
        {
          iVal = (a * bV[idI].vtX) + (b * bV[idI].vtY) + d;
          intersect = (iVal <= bV[0].vtZ) && (iVal >= bV[1].vtZ);
        }
        while((intersect == 0) && (++idI < 4));
        break;
    }
  }
  return(intersect);
}

/*!
* \return	0 if the plane and line do not intersect, 1 if the line
*		segment intersects the plane at a single point or 2 if the
*		line segment is wholly on the plane.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the plane defined by the
*		equation: \f$ax + by + cz + d = 0\f$ and the line segment
*		with end points \f$p_0\f$ and \f$p_1\f$.
* \param	a			Plane X parameter.
* \param	b			Plane Y parameter.
* \param	c			Plane Z parameter.
* \param	d			Other plane parameter.
* \param	p0			First end point of the line segment.
* \param	p1			Second end point of the line segment.
* \param	dstIsn			Destination pointer for point of
*					intersection, may be NULL.
*/
int		WlzGeomPlaneLineIntersect(double a, double b,
					  double c, double d,
					  WlzDVertex3 p0,
					  WlzDVertex3 p1,
					  WlzDVertex3 *dstIsn)
{
  int		intersect = 0;
  double	tD0,
  		tD1,
		tD2;
  
  tD0 = (a * p0.vtX) + (b * p0.vtY) + (c * p0.vtZ) + d;
  if(fabs(tD0) < DBL_EPSILON)
  {
    /* The first line segment end point lies on the plane if the second end
     * point also lies on the plane then there is no unique intersection
     * point. */
    tD1 = (a * p1.vtX) + (b * p1.vtY) + (c * p1.vtZ) + d;
    if(fabs(tD1) < DBL_EPSILON)
    {
      intersect = 2;
    }
    else
    {
      intersect = 1;
      if(dstIsn)
      {
        *dstIsn = p0;
      }
    }
  }
  else
  {
    tD1 = (a * (p0.vtX - p1.vtX)) + (b * (p0.vtY - p1.vtY)) +
	  (c * (p0.vtZ - p1.vtZ));
    if(fabs(tD1) > DBL_EPSILON)
    {
      /* The line segment does not have zero length and the plane equation is
       * valid. */
      tD2 = tD0 / tD1;
      if((tD2 >= 0.0) && (tD2 <= 1.0))
      {
        intersect = 1;
	if(dstIsn)
	{
	  dstIsn->vtX = p0.vtX + (tD2 * (p1.vtX - p0.vtX));
	  dstIsn->vtY = p0.vtY + (tD2 * (p1.vtY - p0.vtY));
	  dstIsn->vtZ = p0.vtZ + (tD2 * (p1.vtZ - p0.vtZ));
	}
      }
    }
  }
  return(intersect);
}

/*!
* \return	An intersection code:
*			- 0 if the plane and triangle do not intersect;
*			- 1 if a single vertex of the triangle is on the
*		            plane;
*		        - 2 if a single edge or a line through the triangle
*			    is on the plane;
*		  	- 3 if the triangle is wholly on the plane.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between a plane and a triangle.
*
*		Tests for an intersection between the plane defined by the
*		equation: \f$ax + by + cz + d = 0\f$ and the triangle
*		with end verticies \f$p_0\f$, \f$p_1\f$ and \f$p_2\f$.
*		If the destination pointers for the intersection points
*		are not NULL and the intersection code is 1 then the single
*		point of intersection is returned in dstIsn0. If the
*		destination pointers are not NULL and the intersection code
*		is either 1 or 2 then the twther the single intersection
*		point is returned in dstIsn0 or the two intersection points
*		are returned in dstIsn0 and dstIsn1.
* \param	a			Plane X parameter.
* \param	b			Plane Y parameter.
* \param	c			Plane Z parameter.
* \param	d			Other plane parameter.
* \param	p0			First triangle vertex.
* \param	p1			Second triangle vertex.
* \param	p2			Third triangle vertex.
* \param	dstIsn0			Destination pointer for first point
*					of intersection, may be NULL.
* \param	dstIsn1			Destination pointer for second point
*					of intersection, may be NULL.
*/
int		WlzGeomPlaneTriangleIntersect(double a, double b,
					      double c, double d,
					      WlzDVertex3 p0,
					      WlzDVertex3 p1,
					      WlzDVertex3 p2,
					      WlzDVertex3 *dstIsn0,
					      WlzDVertex3 *dstIsn1)
{
  int		iCode,
  		valFlg,
  		intersect = 0;
  double	tD0,
  		tD1;
  WlzDVertex3	tDV0,
  		tDV1;
  WlzDVertex3	isn0,
  		isn1;
  int		isnFlg[3];
  WlzDVertex3	isnVal[3];

  valFlg = dstIsn0 && dstIsn1;
  isnFlg[0] = WlzGeomPlaneLineIntersect(a, b, c, d, p0, p1, isnVal + 0);
  isnFlg[1] = WlzGeomPlaneLineIntersect(a, b, c, d, p1, p2, isnVal + 1);
  isnFlg[2] = WlzGeomPlaneLineIntersect(a, b, c, d, p2, p0, isnVal + 2);
  iCode = (isnFlg[2] * 100) + (isnFlg[1] * 10) + isnFlg[0];
  switch(iCode)
  {
    case 000:
      intersect = 0;
      break;
    case 001:
      intersect = 1;
      isn0 = isnVal[0];
      break;
    case 010:
      intersect = 1;
      isn0 = isnVal[1];
      break;
    case 100:
      intersect = 1;
      isn0 = isnVal[2];
      break;
    case 002:
    case 012:
    case 102:
    case 112:
      intersect = 2;
      isn0 = p0;
      isn1 = p1;
      break;
    case 020:
    case 021:
    case 120:
    case 121:
      intersect = 2;
      isn0 = p1;
      isn1 = p2;
      break;
    case 200:
    case 201:
    case 210:
    case 211:
      intersect = 2;
      isn0 = p2;
      isn1 = p0;
      break;
    case 022:
    case 122:
    case 202:
    case 212:
    case 220:
    case 221:
    case 222:
      intersect = 3;
      break;
    case 011:
      WLZ_VTX_3_SUB(tDV0, isnVal[1], isnVal[0]);
      tD0 = WLZ_VTX_3_SQRLEN(tDV0);
      if(tD0 < DBL_EPSILON)
      {
	intersect = 1;
	isn0 = p1;
      }
      else
      {
	intersect = 2;
	isn0 = isnVal[0];
	isn1 = isnVal[1];
      }
      break;
    case 101:
      WLZ_VTX_3_SUB(tDV0, isnVal[2], isnVal[0]);
      tD0 = WLZ_VTX_3_SQRLEN(tDV0);
      if(tD0 < DBL_EPSILON)
      {
	intersect = 1;
	isn0 = p0;
      }
      else
      {
	intersect = 2;
	isn0 = isnVal[2];
	isn1 = isnVal[0];
      }
      break;
    case 110:
      WLZ_VTX_3_SUB(tDV0, isnVal[1], isnVal[2]);
      tD0 = WLZ_VTX_3_SQRLEN(tDV0);
      if(tD0 < DBL_EPSILON)
      {
	intersect = 1;
	isn0 = p2;
      }
      else
      {
	intersect = 2;
	isn0 = isnVal[1];
	isn1 = isnVal[2];
      }
      break;
    case 111:
      WLZ_VTX_3_SUB(tDV0, isnVal[1], isnVal[0]);
      tD0 = WLZ_VTX_3_SQRLEN(tDV0);
      WLZ_VTX_3_SUB(tDV1, isnVal[2], isnVal[1]);
      tD1 = WLZ_VTX_3_SQRLEN(tDV1);
      if((tD0 > DBL_EPSILON) && (tD1 > DBL_EPSILON))
      {
	intersect = 3;
      }
      else
      {
	intersect = 2;
	if(tD0 < DBL_EPSILON)
	{
	  isn0 = p1;
	  isn1 = p2;
	}
	else
	{
	  isn0 = p0;
	  isn1 = p1;
	}
      }
      break;
  }
  if(valFlg)
  {
    *dstIsn0 = isn0;
    *dstIsn1 = isn1;
  }
  return(intersect);
}

/*!
* \return	Distance ratio.
* \ingroup	WlzGeometry
* \brief	Given an ellipse defined by it's centre \f$\mathbf{c}\f$
*		and it's semi axes \f$\mathbf{a}\f$. This function computes
*		the square of the ratio of the distances from the centre
*		of the ellipse to the given point and from the centre of
*		the ellipse in the direction of the given point to the
*		ellipse.
*
*		Equation of ellipse is:
* 		\f[
		{(\frac{x}{a})}^2 + {(\frac{y}{b})}^2 = 1
		\f]
*		and a straight line through the origin:
*		\f[
		y = m x
		\f]
*		Solving for \f$x\f$ and \f$y\f$ at the ellipse gives the
*		square of the distance ratio for given point \f$\mathbf{p}\f$
*		is
*		\f[
		d = \frac{({(p_x - c_x)}^2 + {(p_y - c_y)}^2)
		          (a_x^2 m^2 + a_y^2)}
		         {a_x^2 a_y^2 ( 1 + m^2)}
		\f]
*		with \f$m = \frac{p_y - c_y}{p_x - c_x}\f$.
* 
* \param	centre			Centre of ellipse.
* \param	sAx			Ellipse semi axes, both \f$> 0.0\f$.
* \param	gPnt			Given point.
*/
double		WlzGeomEllipseVxDistSq(WlzDVertex2 centre, WlzDVertex2 sAx,
				       WlzDVertex2 gPnt)
{
  double 	grdSq,
  		dRatSq;
  WlzDVertex2	sAxSq,
  		rPntSq;

  WLZ_VTX_2_SUB(rPntSq, gPnt, centre);
  rPntSq.vtX = rPntSq.vtX * rPntSq.vtX;
  rPntSq.vtY = rPntSq.vtY * rPntSq.vtY;
  sAxSq.vtX = sAx.vtX * sAx.vtX;
  sAxSq.vtY = sAx.vtY * sAx.vtY;
  if(rPntSq.vtX < DBL_EPSILON)
  {
    dRatSq = rPntSq.vtY / sAxSq.vtY;
  }
  else if(rPntSq.vtY < DBL_EPSILON)
  {
    dRatSq = rPntSq.vtX / sAxSq.vtX;
  }
  else
  {
    grdSq = rPntSq.vtY / rPntSq.vtX;
    dRatSq = ((rPntSq.vtX + rPntSq.vtY) * ((sAxSq.vtX * grdSq) + sAxSq.vtY)) /
             (sAxSq.vtX * sAxSq.vtY * (1.0 + grdSq));
  }
  return(dRatSq);
}

/*!
* \return	Hash value.
* \ingroup      WlzGeometry
* \brief	Computes a hash value from a given 3D double precision
*               position.
* \param	pos			Given position.
* \param	tol			Tolerance, \f$x = x \pm tol\f$.
*/
unsigned int	WlzGeomHashVtx3D(WlzDVertex3 pos, double tol)
{
  unsigned int	hVal;
  double	fF,
  		fI;
  const unsigned int pX = 399989, /* These are just different 6 digit primes */
  		pY = 599999,
		pZ = 999983;
  
  fF = modf(pos.vtX, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pX;
  fF *= pY * pZ;
  hVal = ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pY;
  fF *= pZ * pX;
  hVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtZ, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pZ;
  fF *= pX * pY;
  hVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  return(hVal);
}

/*!
* \return	Hash value.
* \ingroup      WlzGeometry
* \brief	Computes a hash value from a given 2D double precision
*               position.
* \param	pos			Given position.
* \param	tol			Tolerance, \f$x = x \pm tol\f$.
*/
unsigned int	WlzGeomHashVtx2D(WlzDVertex2 pos, double tol)
{
  unsigned int	hVal;
  double	fF,
  		fI;
  const unsigned int pX = 399989, /* These are just different 6 digit primes */
  		pY = 599999,
		pZ = 999983;
  
  fF = modf(pos.vtX, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pX;
  fF *= pY * pZ;
  hVal = ((long long )fI + (long long )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pY;
  fF *= pZ * pX;
  hVal ^= ((long long )fI + (long long )fF) & UINT_MAX;
  return(hVal);
}

/*!
* \return	The result of the vertex comparison: -1, 0 or +1.
* \ingroup      WlzGeometry
* \brief	Compares the coordinates of the given 3D double precision
*		vertices to find a signed value for sorting.
* \param	pos0			First vertex.
* \param	pos1			Second vertex.
* \param	tol			Tolerance, \f$x = x \pm tol\f$.
*/
int		WlzGeomCmpVtx3D(WlzDVertex3 pos0, WlzDVertex3 pos1, double tol)
{
  int		cmp;
  WlzDVertex3	cmp3D;

  WLZ_VTX_3_SUB(cmp3D, pos0, pos1);
  if(cmp3D.vtZ < -tol)
  {
    cmp = -1;
  }
  else if(cmp3D.vtZ > tol)
  {
    cmp = 1;
  }
  else
  {
    if(cmp3D.vtY < -tol)
    {
      cmp = -1;
    }
    else if(cmp3D.vtY > tol)
    {
      cmp = 1;
    }
    else
    {
      if(cmp3D.vtX < -tol)
      {
	cmp = -1;
      }
      else if(cmp3D.vtX > tol)
      {
	cmp = 1;
      }
      else
      {
	cmp = 0;
      }
    }
  }
  return(cmp);
}

/*!
* \return	The sign of the vertex comparison: -1, 0 or +1.
* \ingroup      WlzGeometry
* \brief	Compares the coordinates of the given 2D double precision
*		vertices to find a signed value for sorting.
* \param	pos0			First vertex.
* \param	pos1			Second vertex.
* \param	tol			Tolerance, \f$x = x \pm tol\f$.
*/
int		WlzGeomCmpVtx2D(WlzDVertex2 pos0, WlzDVertex2 pos1, double tol)
{
  int		cmp;
  WlzDVertex2	cmp2D;

  WLZ_VTX_2_SUB(cmp2D, pos0, pos1);
  if(cmp2D.vtY < -tol)
  {
    cmp = -1;
  }
  else if(cmp2D.vtY > tol)
  {
    cmp = 1;
  }
  else
  {
    if(cmp2D.vtX < -tol)
    {
      cmp = -1;
    }
    else if(cmp2D.vtX > tol)
    {
      cmp = 1;
    }
    else
    {
      cmp = 0;
    }
  }
  return(cmp);
}

/*!
* \return	Unit vector or zero vector if vertices are coincident.
* \ingroup	WlzGeometry
* \brief	Computes the unit vector
*		\f$\frac{1}{|\mathbf{v}|} \mathbf{v}\f$.
* \param	vec			Given vector, \f$\mathbf{v}\f$.
*/
WlzDVertex2	WlzGeomUnitVector2D(WlzDVertex2 vec)
{
  double	len;

  if((len = WLZ_VTX_2_LENGTH(vec)) > DBL_EPSILON)
  {
    len = 1.0 / len;
    WLZ_VTX_2_SCALE(vec, vec, len);
  }
  else
  {
    vec.vtX = vec.vtY = 0.0;
  }
  return(vec);
}

/*!
* \return	Unit vector or zero vector if vertices are coincident.
* \ingroup	WlzGeometry
* \brief	Computes the unit vector with the direction given by
*		\f$\mathbf{p}_1 - \mathbf{p}_0\f$.
*		If the two given vertices are coincident then a
*		zero vector is returned instead of a unit vector.
* \param	pos1			Position of vertex, \f$\mathbf{p}_1\f$.
* \param	pos0			Position of vertex, \f$\mathbf{p}_0\f$.
*/
WlzDVertex2	WlzGeomUnitVector2D2(WlzDVertex2 pos1, WlzDVertex2 pos0)
{
  double	len;
  WlzDVertex2	vec;

  WLZ_VTX_2_SUB(vec, pos1, pos0);
  if((len = WLZ_VTX_2_LENGTH(vec)) > DBL_EPSILON)
  {
    len = 1.0 / len;
    WLZ_VTX_2_SCALE(vec, vec, len);
  }
  else
  {
    vec.vtX = vec.vtY = 0.0;
  }
  return(vec);
}

/*!
* \return	Non-zero if point is inside diametral circle.
* \ingroup	WlzGeometry
* \brief	Determines whether a point is inside the diametral circle
*		of a line segment.
*		If two vectors \f$\mathbf{v_0}\f$ and \f$\mathbf{v_1}\f$ are
*		directed from the line segment end points to the point,
*		then the angle between the vectors is \f$<\f$ \f$90^\circ\f$
*		if the point is inside the diametral circle, \f$0\f$ if
*		it lies on the circle and \f$>\f$ \f$90^{\circ}\f$ if it
*		lies outside the circle. This is easily tested by
*		\f$\mathbf{v_0} \cdot \mathbf{v_1} < 0\f$.
* \param	lPos0			Vertex at one end of line segment.
* \param	lPos1			Vertex at other end of line segment.
* \param	pos			Position of point.
*/
int		WlzGeomVertexInDiamCircle(WlzDVertex2 lPos0, WlzDVertex2 lPos1,
					  WlzDVertex2 pos)
{
  int		inside;
  double	prod;
  WlzDVertex2	v0,
  		v1;

  WLZ_VTX_2_SUB(v0, pos, lPos0);
  WLZ_VTX_2_SUB(v1, pos, lPos1);
  prod = WLZ_VTX_2_DOT(v0, v1);
  inside = prod < 0.0;
  return(inside);
}

/*!
* \return	Incremented spiral step count.
* \ingroup	WlzGeometry
* \brief	Iterates the given positions coordinates through an
*		expanding integer spiral.
* \param	step			Spiral step count, must be zero
*					when this function is called for the
*					for first step.
* \param	pX			Destination pointer for column
*					coordinate.
* \param	pY			Destination pointer for line
*					coordinate.
*/
int             WlzGeomItrSpiral2I(int step, int *pX, int *pY)
{
  int           ring,
                ring2,
                square;
  const int	lutX[9] = { 1,  0, -1, -1,  0,  0,  1,  1,  1},
  		lutY[9] = { 0,  1,  0,  0, -1, -1,  0,  0,  0};
  if(step <= 0)
  {
    step = 1;
    ++*pX;
  }
  else if(step < 9)
  {
    *pX += lutX[step];
    *pY += lutY[step];
    ++step;
  }
  else
  {
    ++step;
    ring = (int )floor(sqrt(step) + 1) / 2;
    ring2 = 2 * ring;
    square = (ring2 - 1) * (ring2 - 1);
    if(step == square)
    {
      ++*pX;
    }
    else if(step < (square + ring2))
    {
      ++*pY;
    }
    else if(step < (square + (2 * ring2)))
    {
      --*pX;
    }
    else if(step < (square + (3 * ring2)))
    {
      --*pY;
    }
    else
    {
      ++*pX;
    }
  }
  return(step);
}

/*!
* \return	Euclidean distance between the given vertices.
* \ingroup	WlzGeometry
* \brief	Computes square of the Euclidean distance between the given
* 		two vertices.
* \param	v0			First of the given vertices.
* \param	v1			Second of the given vertices.
*/
double		WlzGeomDistSq2D(WlzDVertex2 v0, WlzDVertex2 v1)
{
  double	dst;

  WLZ_VTX_2_SUB(v0, v0, v1);
  dst = WLZ_VTX_2_SQRLEN(v0);
  return(dst);
}

/*!
* \return       Non zero if the area of the triangle is very small.
* \ingroup      WlzGeometry
* \brief        If the unsigned area of the triangle is very small
*               then the only the transform translation coefficients
*               are computed with the other coefficients being set to
*               zero.
*               If the unsigned area of the triangle is not very small
*               then a system of linear equations is solved for the
*               coefficients of the 2D affine transform from the source
*               triangle to the destination triangle.
* \param        xTr                     Transform coordinates for x.
* \param        yTr                     Transform coordinates for y.
* \param        dd                      Twice the area of the source triangle.
* \param        sVx                     Source triangle vertices.
* \param        dVx                     Destination triangle vertices.
* \param        thresh                  Threshold value for twice the area.
*/
int             WlzGeomTriangleAffineSolve(double *xTr, double *yTr, double dd,
                                        WlzDVertex2 *sVx, WlzDVertex2 *dVx,
                                        double thresh)
{
  int           squashed = 0;
  double        tD0,
                tD1,
                tD2;

  if(fabs(dd) < thresh)
  {
    squashed = 1;
    xTr[0] = 0.0;
    xTr[1] = 0.0;
    xTr[2] = (dVx[0].vtX + dVx[1].vtX + dVx[2].vtX -
              sVx[0].vtX - sVx[1].vtX - sVx[2].vtX) / 3.0;
    yTr[0] = 0.0;
    yTr[1] = 0.0;
    yTr[2] = (dVx[0].vtY + dVx[1].vtY + dVx[2].vtY -
              sVx[0].vtY - sVx[1].vtY - sVx[2].vtY) / 3.0;
  }
  else
  {
    squashed = 0;
    dd = 1.0 / dd;
    tD0 = sVx[1].vtY - sVx[2].vtY;
    tD1 = sVx[2].vtY - sVx[0].vtY;
    tD2 = sVx[0].vtY - sVx[1].vtY;
    xTr[0] = ((dVx[0].vtX * tD0) + (dVx[1].vtX * tD1) +
              (dVx[2].vtX * tD2)) * dd;
    yTr[0] = ((dVx[0].vtY * tD0) + (dVx[1].vtY * tD1) +
              (dVx[2].vtY * tD2)) * dd;
    tD0 = sVx[2].vtX - sVx[1].vtX;
    tD1 = sVx[0].vtX - sVx[2].vtX;
    tD2 = sVx[1].vtX - sVx[0].vtX;
    xTr[1] = ((dVx[0].vtX * tD0) + (dVx[1].vtX * tD1) +
              (dVx[2].vtX * tD2)) * dd;
    yTr[1] = ((dVx[0].vtY * tD0) + (dVx[1].vtY * tD1) +
              (dVx[2].vtY * tD2)) * dd;
    tD0 = (sVx[1].vtX * sVx[2].vtY) - (sVx[2].vtX * sVx[1].vtY);
    tD1 = (sVx[2].vtX * sVx[0].vtY) - (sVx[0].vtX * sVx[2].vtY);
    tD2 = (sVx[0].vtX * sVx[1].vtY) - (sVx[1].vtX * sVx[0].vtY);
    xTr[2] = ((dVx[0].vtX * tD0) + (dVx[1].vtX * tD1) +
              (dVx[2].vtX * tD2)) * dd;
    yTr[2] = ((dVx[0].vtY * tD0) + (dVx[1].vtY * tD1) +
              (dVx[2].vtY * tD2)) * dd;
  }
  return(squashed);
}
