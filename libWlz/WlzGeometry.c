#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGeometry_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzGeometry.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Geometric utility functions.
* \ingroup	WlzGeometry
*/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>

extern double			cbrt(double c);
#ifdef _BORLAND_FOR_JATLASVIEWER
/*!
 *  Dummy function to allow compilation with borland bcc32 under Windows
 */
double cbrt(double c) {
   return (0.0);
}
#endif

static int			WlzGeomVtxSortRadialFn(
				  void *p0,
				  int *idxP,
				  int idx0,
				  int idx1);
static int			WlzGeomLineTriangleIntersectEdge3D(
				  WlzDVertex3 org,
				  WlzDVertex3 dir,
				  WlzDVertex3 v0,
				  WlzDVertex3 v1,
				  WlzDVertex3 v2);
static int			WlzGeomTriAABBIsnDir(
				  WlzDVertex3 d,
				  WlzDVertex3 *t,
				  WlzDVertex3 *b);
static int			WlzGeomTetAABBIsnDir(
				  WlzDVertex3 d,
                                  WlzDVertex3 *t,
				  WlzDVertex3 *b);
static int 			WlzGeomTriTri3DCoplanar(
				  WlzDVertex3 n,
				  WlzDVertex3 s[],
				  WlzDVertex3 t[]);
static int			WlzGeomTriTri3DIsn(
				  WlzDVertex3 s[],
				  double sp[],
				  double d[],
				  double d0d1,
				  double d0d2,
				  double is[]);
static int			WlzGeomTriTriPlaneTest(
				  double ds[],
				  WlzDVertex3 *n,
				  double *d,
				  double *d0d1,
				  double *d0d2,
				  WlzDVertex3 s[],
				  WlzDVertex3 t[]);
static double			WlzGeomCot2D3(
				  WlzDVertex2 a,
				  WlzDVertex2 b,
				  WlzDVertex2 c);

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
  const double  tol = ALG_DBL_TOLLERANCE;

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
  if((tD[12] * tD[12]) > tol)
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
* \brief	Tests to set if the given vertex lies within the given
*		triangle using a barycentric coordinates test.
*
*		If a triangle has vertices \f$p_0, p_1, p_2\f$, then any point
*		in the plane containing the triangle can be represented
*		by: \f$p = \lambda_0 p_0 + \lambda_1 p_2 + \lambda_2 p_3\f$
*		subject to the constraint:
*		\f$\lambda_0 + \lambda_1 + \lambda_2 = 1\f$
*		\f$p\f$ is outside the triangle at one or more of 
*		\f$\lambda_0\f$, \f$\lambda_1\f$ and \f$\lambda_2\f$ is -ve.
*		It is inside if all are +ve and on an edge of the
*		triangle if any are close to zero (ie < ALG_DBL_TOLLERANCE).
* \param	p0			First vertex of triangle.
* \param	p1			Second vertex of triangle.
* \param	p2			Third vertex of triangle.
* \param	pP			Given vertex.
*/
int		 WlzGeomVxInTriangle2D(WlzDVertex2 p0, WlzDVertex2 p1,
				       WlzDVertex2 p2, WlzDVertex2 pP)
{
  int		inside = 0;
  double	l0,
		l1,
		l2,
		delta;
  WlzDVertex2	q0,
  		q1,
		qP;
  const double	eps = 1.0e-10;

  WLZ_VTX_2_SUB(q0, p0, p2);
  WLZ_VTX_2_SUB(q1, p1, p2);
  delta = (q0.vtX * q1.vtY) - (q1.vtX * q0.vtY);
  if(fabs(delta) > eps)
  {
    delta = 1.0 / delta;
    WLZ_VTX_2_SUB(qP, pP, p2);
    l0 = delta * ((qP.vtX * q1.vtY) - (q1.vtX * qP.vtY));
    l1 = delta * ((q0.vtX * qP.vtY) - (qP.vtX * q0.vtY));
    l2 = 1.0 - (l0 + l1);
    if((l0 < -eps) || (l1 < -eps) || (l2 < -eps))
    {
      inside = -1;
    }
    else if((l0 > eps) && (l1 > eps) && (l2 > eps))
    {
      inside = 1;
    }
  }
  return(inside);
}

/*!
* \return	Value indicating the position of the vertex with respect
*               to the triangle:
*		  +ve if the vertex is inside the triangle,
*		  0   if the vertex is on an edge of the triangle and
*		  -ve if the vertex is outside the triangle.
* \ingroup	WlzGeometry
* \brief	First finds the closest point on the plane of the triangle
* 		to the given point. Then if the distance from the point
* 		to the plane is less than the given tolerance vvalue tests
* 		to set if the given vertex lies within the given triangle
* 		using a barycentric coordinates test (see
* 		WlzGeomVxInTriangle2D()).
* \param	v0			First vertex of triangle.
* \param	v1			Second vertex of triangle.
* \param	v2			Third vertex of triangle.
* \param	vQ			Given query vertex.
* \param	vvMax			Maximum plane vertex distance.
*/
int		 WlzGeomVxInTriangle3D(WlzDVertex3 v0, WlzDVertex3 v1,
				       WlzDVertex3 v2, WlzDVertex3 vQ,
				       double vPMax)
{
  int		inside = -1;
  double	lnn;
  WlzDVertex3	n,
  		u0,
  		u1,
		uQ;
  const double	eps = 1.0e-10;

  WLZ_VTX_3_SUB(u0, v0, v2);
  WLZ_VTX_3_SUB(u1, v1, v2);
  WLZ_VTX_3_CROSS(n, u0, u1);
  lnn = WLZ_VTX_3_SQRLEN(n);
  if(lnn > eps)
  {
    double	d,
    		ln;
    WlzDVertex3 t0,
		t1,
		t2;

    lnn = 1.0 / lnn;
    WLZ_VTX_3_SUB(uQ, vQ, v2);
    /* Make uQ closest point in plane of triangle and compute distance d
     * from uQ to vQ. */
    ln = sqrt(lnn);
    WLZ_VTX_3_SCALE(t0, n, ln);
    d = WLZ_VTX_3_DOT(uQ, t0);
    if(fabs(d) < vPMax + eps)
    {
      double	l0,
		l1,
		l2;

      inside = 0;
      WLZ_VTX_3_SCALE(t1, t0, d);
      WLZ_VTX_3_SUB(uQ, uQ, t1);
      /* Now have triangle (O, u0, u1) and query vertex uQ all on a plane.
       * Compute the barycentric ccordinates \f$\lambda_0\f$ \f£lambda_1\f$
       * and \f$\lambda_2\f$.
       */
      WLZ_VTX_3_SUB(t0, u1, u0);
      WLZ_VTX_3_SUB(t1, uQ, u0);
      WLZ_VTX_3_CROSS(t2, t0, t1);
      l0 = WLZ_VTX_3_DOT(n, t2) * lnn;
      WLZ_VTX_3_CROSS(t2, u0, uQ);
      l1 = WLZ_VTX_3_DOT(n, t2) * lnn;
      l2 = 1.0 - (l0 + l1);
      if((l0 < -eps) || (l1 < -eps) || (l2 < -eps))
      {
	inside = -1;
      }
      else if((l0 > eps) && (l1 > eps) && (l2 > eps))
      {
	inside = 1;
      }
    }
  }
  return(inside);
}

/*!
* \return	Value indicating the position of the vertex with respect
*               to the tetrahedron:
*		  +ve if the vertex is inside the tetrahedron,
*		  0   if the vertex is on an edge of the tetrahedron and
*		  -ve if the vertex is outside the tetrahedron.
* \ingroup	WlzGeometry
* \brief	Tests to set if the given vertex lies within the given
*		tetrahedron using a barycentric coordinates test.
*
*		If a tetrahedron has vertices \f$p_0, p_1, p_2, p_3\f$,
*		then any point in the 3D space containing the 
*		tetrahedron can be represented by:
*		\f[p = \lambda_0 p_0 + \lambda_1 p_1 + \lambda_2 p_2 +
                       \lambda_3 p_3\f]
*		subject to the constraint:
*		\f[\lambda_0 + \lambda_1 + \lambda_2 + \lambda_3 = 1\f]
*		\f$p\f$ is outside the tetrahedron at one or more of 
*		\f$\lambda_0\f$, \f$\lambda_1\f$, \f$\lambda_2\f$ and
*		\f$\lambda_3\f$ is -ve. It is inside if all are +ve and
*		on an edge of the tetrahedron if any are close to
*		zero (ie fabs(x) < ALG_DBL_TOLLERANCE).
*
* 		The barycentric coordinates are computed by inverting the
*		The tetrahedron vertices
*		\f[
		V = \left[
		    \begin{array}{cccc}
		    vx_0 & vx_1 & vx_2 & vx_3 \\
		    vy_0 & vy_1 & vy_2 & vy_3 \\
		    vz_0 & vz_1 & vz_2 & vz_3 \\
		       1 &    1 &    1 &    1
		    \end{array}
		    \right]
		\f]
*		the barycentric coordinates
*		\f[
                L = \left[
		    \begin{array}{c}
		    \lambda_0 \\
		    \lambda_1 \\
		    \lambda_2 \\
		    \lambda_3
		    \end{array}
		    \right]
                \f]
*              and the point to be queried
*		\f[
                P = \left[
		    \begin{array}{c}
		    px \\
		    py \\
		    pz \\
		    1
		    \end{array}
		    \right]
                \f]
*		can be written \f[V L = P \f] and solved for the
*		the barycentric coordinates using \f[L = V^{-1} P\f].
* \param	v0			First vertex of tetrahedron.
* \param	v1			Second vertex of tetrahedron.
* \param	v2			Third vertex of tetrahedron.
* \param	v3			Fourth vertex of tetrahedron.
* \param	vP			Given vertex.
*/
int		 WlzGeomVxInTetrahedron(WlzDVertex3 v0, WlzDVertex3 v1,
				        WlzDVertex3 v2, WlzDVertex3 v3,
				        WlzDVertex3 vP)
{
  int		inside = 0;
  double	delta,
		l0,
		l1,
		l2,
		l3,
		u0,
		u1,
		u2,
		u3,
		vy01,
		vy02,
		vy03,
		vy12,
		vy13,
		vy23,
		vyz01,
		vyz02,
		vyz03,
		vyz12,
		vyz13,
		vyz23,
                vz01,
  		vz02,
		vz03,
		vz12,
		vz13,
		vz23;
  const double	eps = 1.0e-10;

  vz01 =  v0.vtZ - v1.vtZ;
  vz02 =  v0.vtZ - v2.vtZ;
  vz03 =  v0.vtZ - v3.vtZ;
  vz12 =  v1.vtZ - v2.vtZ;
  vz13 =  v1.vtZ - v3.vtZ;
  vz23 =  v2.vtZ - v3.vtZ;
  delta =   v0.vtX * ( v1.vtY * vz23 - v2.vtY * vz13 + v3.vtY * vz12)
	  + v1.vtX * (-v0.vtY * vz23 + v2.vtY * vz03 - v3.vtY * vz02)
	  + v2.vtX * ( v0.vtY * vz13 - v1.vtY * vz03 + v3.vtY * vz01)
	  + v3.vtX * (-v0.vtY * vz12 + v1.vtY * vz02 - v2.vtY * vz01);
  if(fabs(delta) > eps)
  {
    delta = 1.0 / delta;
    vy01 =  v0.vtY - v1.vtY;
    vy02 =  v0.vtY - v2.vtY;
    vy03 =  v0.vtY - v3.vtY;
    vy12 =  v1.vtY - v2.vtY;
    vy13 =  v1.vtY - v3.vtY;
    vy23 =  v2.vtY - v3.vtY;
    vyz01 = v0.vtY * v1.vtZ - v1.vtY * v0.vtZ;
    vyz02 = v0.vtY * v2.vtZ - v2.vtY * v0.vtZ;
    vyz03 = v0.vtY * v3.vtZ - v3.vtY * v0.vtZ;
    vyz12 = v1.vtY * v2.vtZ - v2.vtY * v1.vtZ;
    vyz13 = v1.vtY * v3.vtZ - v3.vtY * v1.vtZ;
    vyz23 = v2.vtY * v3.vtZ - v3.vtY * v2.vtZ;
    u0 =   v1.vtY * vz23  - v2.vtY * vz13  + v3.vtY * vz12;
    u1 = - v1.vtX * vz23  + v2.vtX * vz13  - v3.vtX * vz12;
    u2 =   v1.vtX * vy23  - v2.vtX * vy13  + v3.vtX * vy12;
    u3 = - v1.vtX * vyz23 + v2.vtX * vyz13 - v3.vtX * vyz12;
    l0 =  delta * (u0 * vP.vtX + u1 * vP.vtY + u2 * vP.vtZ + u3);
    if(l0 < -eps)
    {
      inside = -1;
    }
    else
    {
      u0 = - v0.vtY * vz23  + v2.vtY * vz03  - v3.vtY * vz02;
      u1 =   v0.vtX * vz23  - v2.vtX * vz03  + v3.vtX * vz02;
      u2 = - v0.vtX * vy23  + v2.vtX * vy03  - v3.vtX * vy02;
      u3 =   v0.vtX * vyz23 - v2.vtX * vyz03 + v3.vtX * vyz02;
      l1 =  delta * (u0 * vP.vtX + u1 * vP.vtY + u2 * vP.vtZ + u3);
      if(l1 < -eps)
      {
        inside = -1;
      }
      else
      {
	u0 =   v0.vtY * vz13  - v1.vtY * vz03  + v3.vtY * vz01;
	u1 = - v0.vtX * vz13  + v1.vtX * vz03  - v3.vtX * vz01;
	u2 =   v0.vtX * vy13  - v1.vtX * vy03  + v3.vtX * vy01;
	u3 = - v0.vtX * vyz13 + v1.vtX * vyz03 - v3.vtX * vyz01;
	l2 =  delta * (u0 * vP.vtX + u1 * vP.vtY + u2 * vP.vtZ + u3);
	if(l2 < -eps)
	{
	  inside = -1;
	}
	else
	{
          l3 = 1.0 - (l0 + l1 + l2);
	  if(l3 < -eps)
	  {
	    inside = -1;
	  }
	  else if((l0 > eps) && (l1 > eps) && (l2 > eps) && (l3 > eps))
	  {
	    inside = 1;
	  }
        }
      }
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
* \return	Six times the signed volume of the given tetrahedron.
* \ingroup	WlzGeometry
* \brief	Computes six times the signed volume of the given tetrahedron.
*
* 		The signed volume is computed using simple determinant
* 		evaluation:
*		\f[
		area \times 6 = \left|
		     \begin{array}{cccc}
		     x_0 & y_0 & z_0 & 1 \\
		     x_1 & y_1 & z_1 & 1 \\
		     x_2 & y_2 & z_2 & 1 \\
		     x_3 & y_3 & z_3 & 1
		     \end{array}
		     \right|
		 \f]
*		which can be written
*		\f[
                area \times 6 =
  x_0 ((y_2 z_3 + y_1 (z_2 - z_3)) - (y_3 z_2 + z_1 (y_2 - y_3))) 
- y_0 ((x_2 z_3 + x_1 (z_2 - z_3)) - (x_3 z_2 + z_1 (x_2 - x_3)))
+ z_0 ((x_2 y_3 + x_1 (y_2 - y_3)) - (x_3 y_2 + y_1 (x_2 - x_3)))
- x_1 (y_2 z_3 - y_3 z_2) + y_1 (x_2 z_3 - x_3 z_2) - z_1 (x_2 y_3 - x_3 y_2)
		\f]
*		Simple evaluation of this determinant is not robust.
* \param	vx0			First vertex of tetrahedron.
* \param	vx1			Second vertex of tetrahedron.
* \param	vx2			Third vertex of tetrahedron.
* \param	vx3			Forth vertex of tetrahedron.
*/
double		WlzGeomTetraSnVolume6(WlzDVertex3 vx0, WlzDVertex3 vx1,
				       WlzDVertex3 vx2, WlzDVertex3 vx3)
{
  double	vol6;

  vol6 = vx0.vtX * ((vx2.vtY * vx3.vtZ + vx1.vtY * (vx2.vtZ - vx3.vtZ)) -
                    (vx3.vtY * vx2.vtZ + vx1.vtZ * (vx2.vtY - vx3.vtY))) -
	 vx0.vtY * ((vx2.vtX * vx3.vtZ + vx1.vtX * (vx2.vtZ - vx3.vtZ)) -
	            (vx3.vtX * vx2.vtZ + vx1.vtZ * (vx2.vtX - vx3.vtX))) +
	 vx0.vtZ * ((vx2.vtX * vx3.vtY + vx1.vtX * (vx2.vtY - vx3.vtY)) -
	            (vx3.vtX * vx2.vtY + vx1.vtY * (vx2.vtX - vx3.vtX))) -
	 vx1.vtX * (vx2.vtY * vx3.vtZ - vx3.vtY * vx2.vtZ) +
	 vx1.vtY * (vx2.vtX * vx3.vtZ - vx3.vtX * vx2.vtZ) -
	 vx1.vtZ * (vx2.vtX * vx3.vtY - vx3.vtX * vx2.vtY);
  return(vol6);
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
*		using the ALG_DBL_TOLLERANCE tollerance value.
*               This is taken from J. O'Rourke: Computational Geometry
*               in C, p250, but has ben modified to include the use of
*		ALG_DBL_TOLLERANCE.
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
  const double  tol = ALG_DBL_TOLLERANCE;

  ict = p0;
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
  if(fabs(dn) < tol)
  {
    /* Line segments are parallel. */
    if((fabs(sp) < tol) && (fabs(tp) < tol))
    {
      /* Line segments are coincident. */
      ic = 0;
      if((fabs(p0.vtX - q0.vtX) < tol) && (fabs(p0.vtY - q0.vtY) < tol))
      {
        ict = p0;
	++ic;
      }
      if((fabs(p0.vtX - q1.vtX) < tol) && (fabs(p0.vtY - q1.vtY) < tol))
      {
        ict = p0;
	++ic;
      }
      if((fabs(p1.vtX - q0.vtX) < tol) && (fabs(p1.vtY - q0.vtY) < tol))
      {
        ict = p1;
	++ic;
      }
      if((fabs(p1.vtX - q1.vtX) < tol) && (fabs(p1.vtY - q1.vtY) < tol))
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
    if((sp >= -tol) && (sp - dna < tol) && (tp >= -tol) && (tp - dna < tol))
    {
      /* Line segments intersect. */
      intersect = ((sp > tol) && (sp - dna < tol) &&
		   (tp > tol) && (tp - dna < tol))? 3: 2;
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
  const double  tol = ALG_DBL_TOLLERANCE;

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
	  if((s0.vtX > tol) && (s1.vtX > tol))
	  {
	    tst = (s0.vtY / s0.vtX) - (s1.vtY / s1.vtX);
	  }
	  break;
        case 1:
	  if((s0.vtY > tol) && (s1.vtY > tol))
	  {
	    tst = (s1.vtX / s1.vtY) - (s0.vtX / s0.vtY);
	  }
	  break;
        case 2:
	  if((s0.vtY > tol) && (s1.vtY > tol))
	  {
	    tst = (s0.vtX / s0.vtY) - (s1.vtX / s1.vtY);
	  }
	  break;
        case 3:
	  if((s0.vtX > tol) && (s1.vtX > tol))
	  {
	    tst = (s1.vtY / s1.vtX) - (s0.vtY / s0.vtX);
	  }
	  break;
        case 4:
	  if((s0.vtX > tol) && (s1.vtX > tol))
	  {
	    tst = (s0.vtY / s0.vtX) - (s1.vtY / s1.vtX);
	  }
	  break;
        case 5:
	  if((s0.vtY > tol) && (s1.vtY > tol))
	  {
	    tst = (s1.vtX / s1.vtY) - (s0.vtX / s0.vtY);
	  }
	  break;
        case 6:
	  if((s0.vtY > tol) && (s1.vtY > tol))
	  {
	    tst = (s0.vtX / s0.vtY) - (s1.vtX / s1.vtY);
	  }
	  break;
        case 7:
	  if((s0.vtX > tol) && (s1.vtX > tol))
	  {
	    tst = (s1.vtY / s1.vtX) - (s0.vtY / s0.vtX);
	  }
	  break;
      }
      if(tst > tol)
      {
        cmp = 1.0;
      }
      else if(tst < -tol)
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
* \return	Non-zero if node positions are equal, else 0.
* \ingroup	WlzGeometry
* \brief	Checks to see if two verticies are the same
*               within some tollerance.
* \param	pos0			First node position.
* \param	pos1			Second node position.
* \param	tol			Tollerance value.
*/
int		WlzGeomVtxEqual3D(WlzDVertex3 pos0, WlzDVertex3 pos1,
				  double tol)
{
  int		equal;

  equal = (fabs(pos0.vtX - pos1.vtX) < tol) &&
          (fabs(pos0.vtY - pos1.vtY) < tol) &&
	  (fabs(pos0.vtZ - pos1.vtZ) < tol);
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
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_3_SUB(v1, v1, v0); 
  WLZ_VTX_3_SUB(v2, v2, v0); 
  WLZ_VTX_3_CROSS(nrm, v1, v2);
  len = WLZ_VTX_3_LENGTH(nrm);
  if(len > tol)
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
* \return	The result of the intersection test: 0 - no intersection,
* 		1 - triangle and box are touch or intersect.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the given triangle and
* 		the axis aligned bounding box using the Separating Axis
* 		Theorem (SAT).
*
* 		Given an axis aligned bounding box and a triangle this
* 		function tests for an intersection using the Separating
* 		Axis Theorem (SAT) which states : Two 2D convex domains do
* 		not intersect iff there exists a line, called a separating
* 		axis, on which projection intervals of the domains do not
* 		intersect.  The minimal set of axes that need to be considered
* 		is formed by the normals to all edges of the (polygonal)
* 		domains. For an axis aligned bounding box and a triangle
* 		in 2D, this is equivalent to testing for the intersection
* 		of the given axis aligned bounding box with the axis
* 		aligned bounding box of the triangle and the axes
* 		normal to the faces of the triangle. The mathematics
* 		are simplified by the box being axis aligned.
*
* 		The algorithm may return false positives when the domains
* 		are very close to touching.
* \param	t0			First vertex of triangle.
* \param	t1			Second vertex of triangle.
* \param	t2			Third vertex of triangle.
* \param	b0			Minimum coordinates of axis aligned
* 					bounding box.
* \param	b1			Maximum coordinates of axis aligned
* 					bounding box.
* \param	tst			Determines the actual intersection
* 					tests used:
* 					0 - AABB / triangle.
* 					1 - AABB / AABB(triangle) only.
* 					2 - AABB / triangle omitting the
* 					    AABB / AABB(triangle) test
* 					    this is probably only useful if
* 					    the AABB / AABB(triangle) are
* 					    known to intersect.
*/
int		WlzGeomTriangleAABBIntersect2D(WlzDVertex2 t0,
				WlzDVertex2 t1, WlzDVertex2 t2,
				WlzDVertex2 b0, WlzDVertex2 b1,
				int tst)
{
  int		idx,
  		isn = 1;
  WlzDVertex2	b,
  		c;
  WlzDVertex2	t[3];
  const double	tol = 10.0 * ALG_DBL_TOLLERANCE;

  /* Make origin centroid of the AABB. */
  c.vtX = (b0.vtX + b1.vtX) * 0.5;
  c.vtY = (b0.vtY + b1.vtY) * 0.5;
  WLZ_VTX_2_SUB(b, b1, c);
  WLZ_VTX_2_SUB(t[0], t0, c);
  WLZ_VTX_2_SUB(t[1], t1, c);
  WLZ_VTX_2_SUB(t[2], t2, c);
  /* Check AABB and the AABB of the triangle intersect. */
  if((tst == 0) || (tst == 1))
  {
    /* Compute the AABB of the triangle. */
    WlzDVertex2	bT[2];

    bT[0] = bT[1] = t[0];
    for(idx = 1; idx <= 2; ++idx)
    {
      if(t[idx].vtX < bT[0].vtX)
      {
	bT[0].vtX = t[idx].vtX;
      }
      else if(t[idx].vtX > bT[1].vtX)
      {
	bT[1].vtX = t[idx].vtX;
      }
      if(t[idx].vtY < bT[0].vtY)
      {
	bT[0].vtY = t[idx].vtY;
      }
      else if(t[idx].vtY > bT[1].vtY)
      {
	bT[1].vtY = t[idx].vtY;
      }
    }
    /* Compare AABB of triangle with given AABB. If is an intersection
     * when using tolerance then there may be an intersection. Set
     * intersection (will keep looking below). */
    if((-b.vtX - bT[1].vtX > tol) || (bT[0].vtX - b.vtX > tol) ||
       (-b.vtY - bT[1].vtY > tol) || (bT[0].vtY - b.vtY > tol))
    {
      isn = 0;
    }
  }
  if((tst == 0) || (tst == 2))
  {
    if(isn != 0)
    {
      WlzDVertex2	e,
		  f;
      double	p[2],
		  q[3];

      for(idx = 0; idx < 3; ++idx)
      {
	/* Compute an edge vector for the triangle. */
	WLZ_VTX_2_SUB(e, t[(idx + 1) % 3], t[idx]);
	/* Project two vertices onto perpendicular to first edge. */
	p[0] = t[idx].vtX * e.vtY - t[idx].vtY * e.vtX;
	p[1] = t[(idx + 2) % 3].vtX * e.vtY - t[(idx + 2) % 3].vtY * e.vtX;
	if(p[0] > p[1])
	{
	  q[2] = p[0]; p[0] = p[1]; p[1] = q[2];
	}
	/* Project AABB vertices onto perpendicular and find limits. */
	f.vtX = b.vtX * e.vtY;
	f.vtY = b.vtY * e.vtX;
	q[0] = -f.vtX + f.vtY;
	q[1] =  f.vtX + f.vtY;
	if(q[1] < q[0])
	{
	  q[2] = q[0]; q[0] = q[1]; q[1] = q[2];
	}
	q[2] =  f.vtX - f.vtY;
	if(q[2] < q[0])
	{
	  q[0] = q[2];
	}
	else if(q[2] > q[1])
	{
	  q[1] = q[2];
	}
	q[2] = -f.vtX - f.vtY;
	if(q[2] < q[0])
	{
	  q[0] = q[2];
	}
	else if(q[2] > q[1])
	{
	  q[1] = q[2];
	}
	/* Look for intersection of projections. */
	if((p[0] - q[1] > -tol) || (q[0] + p[1] > -tol))
	{
	  isn = 0;
	  break;
	}
      }
    }
  }
  return(isn);
}

/*!
* \return	The result of the intersection test: 0 - no intersection,
* 		1 - tetrahedron and box touch or intersect.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the given tetrahedron and
* 		the axis aligned bounding box using the Separating Axis
* 		Theorem (SAT).
*
* 		Given an axis aligned bounding box and a tetrahedron this
* 		function tests for an intersection using the Separating
* 		Axis Theorem (SAT) which states : Two 3D convex domains do
* 		not intersect iff there exists a line, called a separating
* 		axis, on which projection intervals of the domains do not
* 		intersect.  The minimal set of axes that need to be considered
* 		is formed by the normals to all faces of the (polyhedral)
* 		domains and the cross product of all edge combinations in
* 		which one edge is from each polyhedron. For an axis aligned
* 		bounding box and a tetrahedron in 3D, this is equivalent to
* 		testing for the intersection of the given axis aligned
* 		bounding box with the axis aligned bounding box of the
* 		tetrahedron, testing for intersections on the axes
* 		normal to the faces of the tetrahedron and testing
* 		for intersection along the cross product of the axis
* 		aligned bounding box - tetrahedron edges. The mathematics
* 		are simplified by the box being axis aligned.
*
* 		The algorithm may return false positives when the domains
* 		are very close to touching.
* \param	t0			First vertex of tetrahedron.
* \param	t1			Second vertex of tetrahedron.
* \param	t2			Third vertex of tetrahedron.
* \param	t3			Fourth vertex of tetrahedron.
* \param	b0			Minimum coordinates of axis aligned
* 					bounding box.
* \param	b1			Maximum coordinates of axis aligned
* 					bounding box.
* \param	tst			Determines the actual intersection
* 					tests used:
* 					0 - AABB / tetrahedron.
* 					1 - AABB / AABB(tetrahedron) only.
* 					2 - AABB / tetrahedron omitting the
* 					    AABB / AABB(tetrahedron) test
* 					    this is probably only useful if
* 					    the AABB / AABB(tetrahedron) are
* 					    known to intersect.
*/
int		WlzGeomTetrahedronAABBIntersect3D(WlzDVertex3 t0,
				WlzDVertex3 t1, WlzDVertex3 t2,
				WlzDVertex3 t3, WlzDVertex3 b0,
				WlzDVertex3 b1, int tst)
{
  int		idx,
  		isn = 1;
  double	l;
  WlzDVertex3	b,
  		c,
		d;
  WlzDVertex3	e[6],
  		t[4],
		x[8];
  const double	tol = 10.0 * ALG_DBL_TOLLERANCE;

  /* Make origin centroid of the AABB. */
  c.vtX = (b0.vtX + b1.vtX) * 0.5;
  c.vtY = (b0.vtY + b1.vtY) * 0.5;
  c.vtZ = (b0.vtZ + b1.vtZ) * 0.5;
  WLZ_VTX_3_SUB(b, b1, c);
  WLZ_VTX_3_SUB(t[0], t0, c);
  WLZ_VTX_3_SUB(t[1], t1, c);
  WLZ_VTX_3_SUB(t[2], t2, c);
  WLZ_VTX_3_SUB(t[3], t3, c);
  /* Check AABB and the AABB of the tetrahedron intersect. This is equivalent
   * to checking for an intersection using projections of the tetrahedron
   * onto vectors perpendicular to the faces of the AABB. */
  if((tst == 0) || (tst == 1))
  {
    /* Compute the AABB of the tetrahedron. */
    WlzDVertex3	bT[2];

    bT[0] = bT[1] = t[0];
    for(idx = 1; idx <= 3; ++idx)
    {
      if(t[idx].vtX < bT[0].vtX)
      {
	bT[0].vtX = t[idx].vtX;
      }
      else if(t[idx].vtX > bT[1].vtX)
      {
	bT[1].vtX = t[idx].vtX;
      }
      if(t[idx].vtY < bT[0].vtY)
      {
	bT[0].vtY = t[idx].vtY;
      }
      else if(t[idx].vtY > bT[1].vtY)
      {
	bT[1].vtY = t[idx].vtY;
      }
      if(t[idx].vtZ < bT[0].vtZ)
      {
	bT[0].vtZ = t[idx].vtZ;
      }
      else if(t[idx].vtZ > bT[1].vtZ)
      {
	bT[1].vtZ = t[idx].vtZ;
      }
    }
    /* Compare AABB of triangle with given AABB. */
    if((-b.vtX - bT[1].vtX > tol) || (bT[0].vtX - b.vtX > tol) ||
       (-b.vtY - bT[1].vtY > tol) || (bT[0].vtY - b.vtY > tol) ||
       (-b.vtZ - bT[1].vtZ > tol) || (bT[0].vtZ - b.vtZ > tol))
    {
      isn = 0;        /* No intersection of the AABB with AABB(tetrahedron). */
    }
  }
  /* Check for intersection using projections of the AABB onto vectors
   * perpendicular to the faces of the tetrahedron. */
  if((tst == 0) || (tst == 2))
  {
    if(isn != 0)
    {
      /* Compute the 6 edge vectors for the tetrahedron with the first 3 being
       * along the edges of one face. */
      WLZ_VTX_3_SUB(e[0], t[1], t[0]);
      WLZ_VTX_3_SUB(e[1], t[2], t[1]);
      WLZ_VTX_3_SUB(e[2], t[0], t[2]);
      WLZ_VTX_3_SUB(e[3], t[3], t[0]);
      WLZ_VTX_3_SUB(e[4], t[3], t[1]);
      WLZ_VTX_3_SUB(e[5], t[3], t[2]);
      x[0].vtX = -b.vtX; x[0].vtY = -b.vtY; x[0].vtZ = -b.vtZ;
      x[1].vtX = -b.vtX; x[1].vtY = -b.vtY; x[1].vtZ =  b.vtZ;
      x[2].vtX = -b.vtX; x[2].vtY =  b.vtY; x[2].vtZ = -b.vtZ;
      x[3].vtX = -b.vtX; x[3].vtY =  b.vtY; x[3].vtZ =  b.vtZ;
      x[4].vtX =  b.vtX; x[4].vtY = -b.vtY; x[4].vtZ = -b.vtZ;
      x[5].vtX =  b.vtX; x[5].vtY = -b.vtY; x[5].vtZ =  b.vtZ;
      x[6].vtX =  b.vtX; x[6].vtY =  b.vtY; x[6].vtZ = -b.vtZ;
      x[7].vtX =  b.vtX; x[7].vtY =  b.vtY; x[7].vtZ =  b.vtZ;
      /* For each direction vector normal to a tetrahedron face project
       * the AABB and check for intersection. */
      idx = 0;
      WLZ_VTX_3_CROSS(d, e[0], e[1]); 
      l = WLZ_VTX_3_SQRLEN(d);
      if(l > tol)
      {
	isn = WlzGeomTetAABBIsnDir(d, t, x);
      }
      while((isn != 0) && (idx < 3))
      {
	WLZ_VTX_3_CROSS(d, e[idx], e[idx + 3]); 
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTetAABBIsnDir(d, t, x);
	}
	++idx;
      }
    }
    /* Check for intersection along vectors which are the cross product of
     * all possible edge combinations where one edge is from the AABB and
     * the other is from the tetrahedron. */
    if(isn != 0)
    {
      /* Cross product of tetrahedron edges with the AABB edges can be
       * computed avoiding a full cross product as (0, -z, y), (z, 0, -x),
       * and (-y, x, 0) where (x, y, z) is the tetrahedron edge vector.
       * Edge - edge intersections with the edges parallel and consequent
       * zero length cross product are ignored because the AABB - AABB
       * intersection test above will find these intersections. */
      idx = 0;
      do
      {
	d.vtX =  0.0;
	d.vtY = -e[idx].vtZ;
	d.vtZ =  e[idx].vtY;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTetAABBIsnDir(d, t, x);
	  if(isn == 0)
	  {
	    break;
	  }
	}
	d.vtX =  e[idx].vtZ;
	d.vtY =  0.0;
	d.vtZ = -e[idx].vtX;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTetAABBIsnDir(d, t, x);
	  if(isn == 0)
	  {
	    break;
	  }
	}
	d.vtX = -e[idx].vtY;
	d.vtY =  e[idx].vtX;
	d.vtZ =  0.0;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTetAABBIsnDir(d, t, x);
	}
      } while((isn >= 0) && (++idx < 6));
    }
  }
  return(isn);
}

/*!
* \return	Intersection code : 0 - no intersection,
* 		1 - tetrahedron and box touch or intersect.
* \ingroup	WlzGeometry
* \brief	Intersection interval test code for
* 		WlzGeomTetrahedronAABBIntersect3D().
* \param	d			Direction vector.
* \param	t			Array of four vertices of the
* 					tetrahedron.
* \param	b			Array of 8 vertices of the box.
*/
static int	WlzGeomTetAABBIsnDir(WlzDVertex3 d,
                                     WlzDVertex3 t[], WlzDVertex3 b[])
{
  int		idx;
  double	f;
  double	p[2],
  		q[2];
  int		isn = 1;
  const double  tol = ALG_DBL_TOLLERANCE;

  p[0] = p[1] = WLZ_VTX_3_DOT(d, t[0]);
  for(idx = 1; idx < 4; ++idx)
  {
    f = WLZ_VTX_3_DOT(d, t[idx]);
    if(f < p[0])
    {
      p[0] = f;
    }
    else if(f > p[1])
    {
      p[1] = f;
    }
  }
  q[0] = q[1] = WLZ_VTX_3_DOT(d, b[0]);
  for(idx = 1; idx < 8; ++idx)
  {
    f = WLZ_VTX_3_DOT(d, b[idx]);
    if(f < q[0])
    {
      q[0] = f;
    }
    else if(f > q[1])
    {
      q[1] = f;
    }
  }
  if((p[0] - q[1] > -tol) || (q[0] - p[1] > -tol))
  {
    isn = 0;              /* No intersection of the AABB with tetrahedron. */
  }
  return(isn);
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
  const double  tol = ALG_DBL_TOLLERANCE;

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
  if(aP[maxP] > tol)
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
  const double  tol = ALG_DBL_TOLLERANCE;
  
  tD0 = (a * p0.vtX) + (b * p0.vtY) + (c * p0.vtZ) + d;
  if(fabs(tD0) < tol)
  {
    /* The first line segment end point lies on the plane if the second end
     * point also lies on the plane then there is no unique intersection
     * point. */
    tD1 = (a * p1.vtX) + (b * p1.vtY) + (c * p1.vtZ) + d;
    if(fabs(tD1) < tol)
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
    if(fabs(tD1) > tol)
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
  const double  tol = ALG_DBL_TOLLERANCE;

  isn0 = isn1 = p0;
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
      if(tD0 < tol)
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
      if(tD0 < tol)
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
      if(tD0 < tol)
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
      if((tD0 > tol) && (tD1 > tol))
      {
	intersect = 3;
      }
      else
      {
	intersect = 2;
	if(tD0 < tol)
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
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_2_SUB(rPntSq, gPnt, centre);
  rPntSq.vtX = rPntSq.vtX * rPntSq.vtX;
  rPntSq.vtY = rPntSq.vtY * rPntSq.vtY;
  sAxSq.vtX = sAx.vtX * sAx.vtX;
  sAxSq.vtY = sAx.vtY * sAx.vtY;
  if(rPntSq.vtX < tol)
  {
    dRatSq = rPntSq.vtY / sAxSq.vtY;
  }
  else if(rPntSq.vtY < tol)
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
  hVal = ((WlzLong )fI + (WlzLong )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pY;
  fF *= pZ * pX;
  hVal ^= ((WlzLong )fI + (WlzLong )fF) & UINT_MAX;
  fF = modf(pos.vtZ, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pZ;
  fF *= pX * pY;
  hVal ^= ((WlzLong )fI + (WlzLong )fF) & UINT_MAX;
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
  hVal = ((WlzLong )fI + (WlzLong )fF) & UINT_MAX;
  fF = modf(pos.vtY, &fI);
  fF = floor(fF / tol) * tol;
  fI *= pY;
  fF *= pZ * pX;
  hVal ^= ((WlzLong )fI + (WlzLong )fF) & UINT_MAX;
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
* \return	Unit vector or zero vector if components are both zero.
* \ingroup	WlzGeometry
* \brief	Computes the 2D unit vector
*		\f$\frac{1}{|\mathbf{v}|} \mathbf{v}\f$.
* \param	vec			Given vector, \f$\mathbf{v}\f$.
*/
WlzDVertex2	WlzGeomUnitVector2D(WlzDVertex2 vec)
{
  double	len;
  const double  tol = ALG_DBL_TOLLERANCE;

  if((len = WLZ_VTX_2_LENGTH(vec)) > tol)
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
* \return	Unit vector or zero vector if components are both zero.
* \ingroup	WlzGeometry
* \brief	Computes the 3D unit vector
*		\f$\frac{1}{|\mathbf{v}|} \mathbf{v}\f$.
* \param	vec			Given vector, \f$\mathbf{v}\f$.
*/
WlzDVertex3	WlzGeomUnitVector3D(WlzDVertex3 vec)
{
  double	len;
  const double  tol = ALG_DBL_TOLLERANCE;

  if((len = WLZ_VTX_3_LENGTH(vec)) > tol)
  {
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(vec, vec, len);
  }
  else
  {
    vec.vtX = vec.vtY = vec.vtZ = 0.0;
  }
  return(vec);
}

/*!
* \return	Unit vector or zero vector if vertices are coincident.
* \ingroup	WlzGeometry
* \brief	Computes the unit 2D vector with the direction given by
*		\f$\mathbf{p}_1 - \mathbf{p}_0\f$.
*		If the two given vertices are coincident then a
*		zero vector is returned instead of a unit vector.
* \param	v1			Position of vertex, \f$\mathbf{p}_1\f$.
* \param	v0			Position of vertex, \f$\mathbf{p}_0\f$.
*/
WlzDVertex2	WlzGeomUnitVector2D2(WlzDVertex2 v1, WlzDVertex2 v0)
{
  double	len;
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_2_SUB(v1, v1, v0);
  if((len = WLZ_VTX_2_LENGTH(v1)) > tol)
  {
    len = 1.0 / len;
    WLZ_VTX_2_SCALE(v1, v1, len);
  }
  else
  {
    v1.vtX = v1.vtY = 0.0;
  }
  return(v1);
}

/*!
* \return	Unit vector or zero vector if vertices are coincident.
* \ingroup	WlzGeometry
* \brief	Computes the unit 3D vector with the direction given by
*		\f$\mathbf{p}_1 - \mathbf{p}_0\f$.
*		If the two given vertices are coincident then a
*		zero vector is returned instead of a unit vector.
* \param	v1			Position of vertex, \f$\mathbf{p}_1\f$.
* \param	v0			Position of vertex, \f$\mathbf{p}_0\f$.
*/
WlzDVertex3	WlzGeomUnitVector3D2(WlzDVertex3 v1, WlzDVertex3 v0)
{
  double	len;
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_3_SUB(v1, v1, v0);
  if((len = WLZ_VTX_3_LENGTH(v1)) > tol)
  {
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(v1, v1, len);
  }
  else
  {
    v1.vtX = v1.vtY = v1.vtZ = 0.0;
  }
  return(v1);
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
* \return	Ring of spiral.
* \ingroup	WlzGeometry
* \brief	Computes the ring of a spiral. If two rings differ by more
* 		than one then at least one itteration outwards on the spiral
* 		has been performed between the rings.
* \param	step			Spiral step count.
*/
int		WlzGeomItrSpiralRing(int step)
{
  int ring;

  ring = (int )floor(sqrt(step) + 1) / 2;
  return(ring);
}

/*!
* \return	Incremented spiral step count.
* \ingroup	WlzGeometry
* \brief	Iterates the given positions coordinates through a 2D
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
* \return	Shell of spiral.
* \ingroup	WlzGeometry
* \brief	Computes the shell of a spiral. If two shells differ by more
* 		than one then at least one itteration outwards on the spiral
* 		has been performed between the shells.
* \param	step			Spiral step count.
*/
int		WlzGeomItrSpiralShell(int step)
{
  int 		shell;

  shell = ((int )(floor(cbrt(step))) + 1) / 2;
  return(shell);
}

/*!
* \return	Incremented spiral step count.
* \ingroup	WlzGeometry
* \brief	Iterates the given positions coordinates through a 3D
*		expanding integer spiral.
* \param	step			Spiral step count, must be zero
*					when this function is called for the
*					for first step.
* \param	pX			Destination pointer for column
*					coordinate.
* \param	pY			Destination pointer for line
*					coordinate.
* \param	pZ			Destination pointer for plane
*					coordinate.
*/
int             WlzGeomItrSpiral3I(int step, int *pX, int *pY, int *pZ)
{
  int           tI0,
  		intermediate,
		ring,
  		shell,
		stepInPn,
		stepInPn1,
		stepInRing,
		stepInShell;
  const int	lutX[27] = { 0,
                             0,  1,  1,  0, -1, -1, -1,  0,  1,
			         1,  1,  0, -1, -1, -1,  0,  1,
			     0,  1,  1,  0, -1, -1, -1,  0,  1},
  		lutY[27] = { 0,
		             0,  0,  1,  1,  1,  0, -1, -1, -1,
		                 0,  1,  1,  1,  0, -1, -1, -1,
		             0,  0,  1,  1,  1,  0, -1, -1, -1},
  		lutZ[27] = { 0,
		            -1, -1, -1, -1, -1, -1, -1, -1, -1,
		                 0,  0,  0,  0,  0,  0,  0,  0,
		             1,  1,  1,  1,  1,  1,  1,  1,  1};
  if(step < 27)
  {
    if(step < 0)
    {
      step = 0;
    }
    *pX = lutX[step];
    *pY = lutY[step];
    *pZ = lutZ[step];
  }
  else
  {
    intermediate = 0;
    /* Each shell + all sub shells has (2n + 1)^3 steps, where n is the shell
     * number. */
    shell = ((int )(floor(cbrt(step))) + 1) / 2;
    tI0 = 2 * shell - 1;
    stepInShell = step - tI0 * tI0 * tI0; 
    /* First and last planes of a shell have (2n +1)^2 steps.
     * The (2n - 1) intermediate planes each have 8n steps. */
    tI0 = 2 * shell + 1;
    stepInPn1 = tI0 * tI0;
    if(stepInShell <= stepInPn1)
    {
      /* In first plane of shell. */
      stepInPn = stepInShell;
      *pZ = -shell;
    }
    else
    {
      tI0 = stepInPn1 + 8 * shell * (2 * shell - 1);
      if(stepInShell > tI0)
      {
        /* In last plane of shell. */
	stepInPn = stepInShell - tI0; 
	*pZ = shell;
      }
      else
      {
        /* In intermediate plane of shell. */
        intermediate = 1;
	tI0 = stepInShell - stepInPn1;
	stepInPn = tI0 % (8 * shell);
	*pZ = (tI0 / (8 * shell)) + 1 - shell;
      }
    }
    if(intermediate)
    {
      ring = shell;
      stepInRing = stepInPn;
    }
    else
    {
      if(stepInPn > 0)
      {
        ring =  ((int )floor(sqrt(stepInPn)) + 1) / 2;
	tI0 = (2 * ring - 1);
	stepInRing = stepInPn - tI0 * tI0;
      }
      else
      {
        ring = 0;
	stepInRing = 0;
      }
    }
    if(ring == 0)
    {
      *pX = 0;
      *pY = 0;
    }
    else
    {
      tI0 = stepInRing % ring;
      switch(tI0)
      {
	case 0:
	  *pX = ring;
	  *pY = stepInRing;
	  break;
	case 1: /* FALLTHROUGH */
	case 2:
	  *pX = (2 * ring) - stepInRing;
	  *pY = ring;
	  break;
	case 3: /* FALLTHROUGH */
	case 4:
	  *pX =  -ring;
	  *pY = (4 * ring) - stepInRing;
	  break;
	case 5: /* FALLTHROUGH */
	case 6:
	  *pX = (6 * ring) - stepInRing;
	  *pY = -ring;
	  break;
	case 7:
	  *pX =  ring;
	  *pY = (8 * ring) - stepInRing;
	  break;
      }
    }
  }
  ++step;
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
* \return	Euclidean distance between the given vertices.
* \ingroup	WlzGeometry
* \brief	Computes square of the Euclidean distance between the given
* 		two vertices.
* \param	v0			First of the given vertices.
* \param	v1			Second of the given vertices.
*/
double		WlzGeomDistSq3D(WlzDVertex3 v0, WlzDVertex3 v1)
{
  double	dst;

  WLZ_VTX_3_SUB(v0, v0, v1);
  dst = WLZ_VTX_3_SQRLEN(v0);
  return(dst);
}

/*!
* \return	Euclidean distance between the given vertices.
* \ingroup	WlzGeometry
* \brief	Computes the Euclidean distance between the given
* 		two vertices.
* \param	v0			First of the given vertices.
* \param	v1			Second of the given vertices.
*/
double		WlzGeomDist2D(WlzDVertex2 v0, WlzDVertex2 v1)
{
  double	dst;

  WLZ_VTX_2_SUB(v0, v0, v1);
  dst = WLZ_VTX_2_LENGTH(v0);
  return(dst);
}

/*!
* \return	Euclidean distance between the given vertices.
* \ingroup	WlzGeometry
* \brief	Computes the Euclidean distance between the given
* 		two vertices.
* \param	v0			First of the given vertices.
* \param	v1			Second of the given vertices.
*/
double		WlzGeomDist3D(WlzDVertex3 v0, WlzDVertex3 v1)
{
  double	dst;

  WLZ_VTX_3_SUB(v0, v0, v1);
  dst = WLZ_VTX_3_LENGTH(v0);
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

/*!
* \return       Non zero if the volume of the tetrahedron is very small.
* \ingroup      WlzGeometry
* \brief        Computes the affine transform coefficients from the
*		source to target tetrahedron.
*
*		This is done by solving for general 3D affine transform
*		matrix which transforms the source tetrahedrons vertices
*		to those of the destination tetrahedron
*		If the destination tetrahedron vertices are given by
*		\f[
                D = \left(
		    \begin{array}{cccc}
                    d_{x0} & d_{x0} & d_{x1} & d_{x2} \\
                    d_{y0} & d_{y0} & d_{y1} & d_{y2} \\
                    d_{z0} & d_{z0} & d_{z1} & d_{z2} \\
                    1      & 1      & 1      & 1
		    \end{array}
		    \right)
                \f]
*		and the source vertices by
*		\f[
                S = \left(
		    \begin{array}{cccc}
                    s_{x0} & s_{x0} & s_{x1} & s_{x2} \\
                    s_{y0} & s_{y0} & s_{y1} & s_{y2} \\
                    s_{z0} & s_{z0} & s_{z1} & s_{z2} \\
                    1      & 1      & 1      & 1
		    \end{array}
		    \right)
                \f]
*		then the transform being sought satisfies
*		\f[
                D = T S
                \f]
*		Solving for \f$T\f$
*		\f[
                T = D S^{-1}
                \f]
*		Setting
*		\f[
		U = S^{-1}
		\f]
		with
*		\f[
                U = \left(
		    \begin{array}{cccc}
                    u_{11} & u_{12} & u_{13} & u_{14} \\
                    u_{21} & u_{22} & u_{23} & u_{24} \\
                    u_{31} & u_{32} & u_{33} & u_{34} \\
                    u_{41} & u_{42} & u_{43} & u_{44}
		    \end{array}
		    \right)
                \f]
*		and
		\f[
		d = \det(S)
		\f]
*		For efficiency the followingcollect common sub-expressions
*		are used
		\f{eqnarray*}
		cz2z3 = s_{z2} - s_{z3} \\
		cz1z3 = s_{z1} - s_{z3} \\
		cz1z2 = s_{z1} - s_{z2} \\
		cy2y3 = s_{y2} - s_{y3} \\
		cy1y3 = s_{y1} - s_{y3} \\
		cy1y2 = s_{y1} - s_{y2} \\
		cy2z3y3z2 = s_{y2} s_{z3} - s_{y3} s_{z2} \\
		cy1z3y3z1 = s_{y1} s_{z3} - s_{y3} s_{z1} \\
		cy1z2y2z1 = s_{y1} s_{z2} - s_{y2} s_{z1} \\
		cz0z3 = s_{z0} - s_{z3} \\
		cz0z2 = s_{z0} - s_{z2} \\
		cy0y3 = s_{y0} - s_{y3} \\
		cy0y2 = s_{y0} - s_{y2} \\
		cy0z3y3z0 = s_{y0} s_{z3} - s_{y3} s_{z0} \\
		cy0z2y2z0 = s_{y0} s_{z2} - s_{y2} s_{z0} \\
		cz0z1 = s_{z0} - s_{z1} \\
		cy0y1 = s_{y0} - s_{y1} \\
		cy0z1y1z0 = s_{y0} s_{z1} - s_{y1} s_{z0}
		\f}
*		giving
*		\f[
		d =   s_{x0}( s_{y1} cz2z3 - s_{y2} cz1z3 + s_{y3} (cz1z2)) 
                    + s_{x1}(-s_{y0} cz2z3 + s_{y2} cz0z3 - s_{y3} (cz0z2))
                    + s_{x2}( s_{y0} cz1z3 - s_{y1} cz0z3 + s_{y3} (cz0z1))
                    + s_{x3}(-s_{y0} cz1z2 + s_{y1} cz0z2 - s_{y2} (cz0z1))
		\f]
*		and for the non-degenerate case, when \f$d \not= 0\f$
*		\f{eqnarray*}
  u11 = \frac{1}{d}( s_{y1} cz2z3 - s_{y2} cz1z3 + s_{y3} cz1z2) \\
  u12 = \frac{1}{d}(-s_{x1} cz2z3 + s_{x2} cz1z3 - s_{x3} cz1z2) \\
  u13 = \frac{1}{d}( s_{x1} cy2y3 - s_{x2} cy1y3 + s_{x3} cy1y2) \\
  u14 = \frac{1}{d}(-s_{x1} cy2z3y3z2 + s_{x2} cy1z3y3z1 - s_{x3} cy1z2y2z1) \\
  u21 = \frac{1}{d}(-s_{y0} cz2z3 + s_{y2} cz0z3 - s_{y3} cz0z2)  \\
  u22 = \frac{1}{d}( s_{x0} cz2z3 - s_{x2} cz0z3 + s_{x3} cz0z2) \\
  u23 = \frac{1}{d}(-s_{x0} cy2y3 + s_{x2} cy0y3 - s_{x3} cy0y2)  \\
  u24 = \frac{1}{d}( s_{x0} cy2z3y3z2 - s_{x2} cy0z3y3z0 + s_{x3} cy0z2y2z0) \\
  u31 = \frac{1}{d}( s_{y0} cz1z3 - s_{y1} cz0z3 + s_{y3} cz0z1) \\
  u32 = \frac{1}{d}(-s_{x0} cz1z3 + s_{x1} cz0z3 - s_{x3} cz0z1) \\
  u33 = \frac{1}{d}( s_{x0} cy1y3 - s_{x1} cy0y3 + s_{x3} cy0y1)  \\
  u34 = \frac{1}{d}(-s_{x0} cy1z3y3z1 + s_{x1} cy0z3y3z0 - s_{x3} cy0z1y1z0) \\
  u41 = \frac{1}{d}(-s_{y0} cz1z2 + s_{y1} cz0z2 - s_{y2} cz0z1) \\
  u42 = \frac{1}{d}( s_{x0} cz1z2 - s_{x1} cz0z2 + s_{x2} cz0z1)  \\
  u43 = \frac{1}{d}(-s_{x0} cy1y2 + s_{x1} cy0y2 - s_{x2} cy0y1) \\
  u44 = \frac{1}{d}( s_{x0} cy1z2y2z1 - s_{x1} cy0z2y2z0 + s_{x2} cy0z1y1z0)
		\f}
*		and so (\f$T = D U\f$)
                \f{eqnarray*}
		t11 = d_{x3} u41 + d_{x2} u31 + d_{x1} u21 + d_{x0} u11 \\
		t12 = d_{x3} u42 + d_{x2} u32 + d_{x1} u22 + d_{x0} u12 \\
		t13 = d_{x3} u43 + d_{x2} u33 + d_{x1} u23 + d_{x0} u13 \\
		t14 = d_{x3} u44 + d_{x2} u34 + d_{x1} u24 + d_{x0} u14 \\
		t21 = d_{y3} u41 + d_{y2} u31 + d_{y1} u21 + d_{y0} u11 \\
		t22 = d_{y3} u42 + d_{y2} u32 + d_{y1} u22 + d_{y0} u12 \\
		t23 = d_{y3} u43 + d_{y2} u33 + d_{y1} u23 + d_{y0} u13 \\
		t24 = d_{y3} u44 + d_{y2} u34 + d_{y1} u24 + d_{y0} u14 \\
		t31 = d_{z3} u41 + d_{z2} u31 + d_{z1} u21 + d_{z0} u11 \\
		t32 = d_{z3} u42 + d_{z2} u32 + d_{z1} u22 + d_{z0} u12  \\
		t33 = d_{z3} u43 + d_{z2} u33 + d_{z1} u23 + d_{z0} u13 \\
		t34 = d_{z3} u44 + d_{z2} u34 + d_{z1} u24 + d_{z0} u14 \\
		t41 = u41 + u31 + u21 + u11 (0) \\
		t42 = u42 + u32 + u22 + u12 (0) \\
		t43 = u43 + u33 + u23 + u13 (0) \\
		t44 = u44 + u34 + u24 + u14 (1)
		\f}
*		For the degenerate cases in which \f$d \approx 0\f$
*		then the transformation is given as a simple translation.
* \param        tr                      Transform matrix with 4x4 contiguous
* 					coefficients which are equivalent
* 					to the base storage of the matrix
*					in a WlzAffineTransform.
* \param        sVx                     Source tetrahedron vertices.
* \param        dVx                     Destination tetrahedron vertices.
* \param	thresh			Threshold value which is the lower
* 					limit of 6 x the source tetrahedron's
*					volume.
*/
int             WlzGeomTetraAffineSolve(double *tr,
					WlzDVertex3 *sVx,
					WlzDVertex3 *dVx,
					double thresh)
{
  int           squashed = 0;
  double	d,
  		cz2z3,
		cz1z3,
		cz1z2,
		cy2y3,
		cy1y3,
		cy1y2,
		cy2z3y3z2,
		cy1z3y3z1,
		cy1z2y2z1,
		cz0z3,
		cz0z2,
		cy0y3,
		cy0y2,
		cy0z3y3z0,
		cy0z2y2z0,
		cz0z1,
		cy0y1,
		cy0z1y1z0;
  double	u[4];

  cz2z3 = sVx[2].vtZ - sVx[3].vtZ;
  cz1z3 = sVx[1].vtZ - sVx[3].vtZ;
  cz1z2 = sVx[1].vtZ - sVx[2].vtZ;
  cy2y3 = sVx[2].vtY - sVx[3].vtY;
  cy1y3 = sVx[1].vtY - sVx[3].vtY;
  cy1y2 = sVx[1].vtY - sVx[2].vtY;
  cy2z3y3z2 = sVx[2].vtY * sVx[3].vtZ - sVx[3].vtY * sVx[2].vtZ;
  cy1z3y3z1 = sVx[1].vtY * sVx[3].vtZ - sVx[3].vtY * sVx[1].vtZ;
  cy1z2y2z1 = sVx[1].vtY * sVx[2].vtZ - sVx[2].vtY * sVx[1].vtZ;
  cz0z3 = sVx[0].vtZ - sVx[3].vtZ;
  cz0z2 = sVx[0].vtZ - sVx[2].vtZ;
  cy0y3 = sVx[0].vtY - sVx[3].vtY;
  cy0y2 = sVx[0].vtY - sVx[2].vtY;
  cy0z3y3z0 = sVx[0].vtY * sVx[3].vtZ - sVx[3].vtY * sVx[0].vtZ;
  cy0z2y2z0 = sVx[0].vtY * sVx[2].vtZ - sVx[2].vtY * sVx[0].vtZ;
  cz0z1 = sVx[0].vtZ - sVx[1].vtZ;
  cy0y1 = sVx[0].vtY - sVx[1].vtY;
  cy0z1y1z0 = sVx[0].vtY * sVx[1].vtZ - sVx[1].vtY * sVx[0].vtZ;
  d =  sVx[0].vtX * ( sVx[1].vtY * cz2z3 - sVx[2].vtY * cz1z3 +
                      sVx[3].vtY * (cz1z2))  +
       sVx[1].vtX * (-sVx[0].vtY * cz2z3 + sVx[2].vtY * cz0z3 -
                      sVx[3].vtY * (cz0z2)) +
       sVx[2].vtX * ( sVx[0].vtY * cz1z3 - sVx[1].vtY * cz0z3 +
                      sVx[3].vtY * (cz0z1)) +
       sVx[3].vtX * (-sVx[0].vtY * cz1z2 + sVx[1].vtY * cz0z2 -
                      sVx[2].vtY * (cz0z1));
  if(fabs(d) < thresh)
  {
    squashed = 1;
    tr[ 0] = tr[ 1] = tr[ 2] = 0.0;
    tr[ 4] = tr[ 5] = tr[ 6] = 0.0;
    tr[ 8] = tr[ 9] = tr[10] = 0.0;
    tr[11] = tr[12] = tr[13] = 0.0;
    tr[ 3] = 0.25 * (dVx[0].vtX + dVx[1].vtX + dVx[2].vtX + dVx[3].vtX -
                     sVx[0].vtX - sVx[1].vtX - sVx[2].vtX - sVx[3].vtX);
    tr[ 7] = 0.25 * (dVx[0].vtY + dVx[1].vtY + dVx[2].vtY + dVx[3].vtY -
                     sVx[0].vtY - sVx[1].vtY - sVx[2].vtY - sVx[3].vtY);
    tr[11] = 0.25 * (dVx[0].vtZ + dVx[1].vtZ + dVx[2].vtZ + dVx[3].vtZ -
                     sVx[0].vtZ - sVx[1].vtZ - sVx[2].vtZ - sVx[3].vtZ);
    tr[15] = 1.0;
  }
  else
  {
    squashed = 0;
    d = 1.0 / d;
    u[0] =  sVx[1].vtY * cz2z3 - sVx[2].vtY * cz1z3 + sVx[3].vtY * cz1z2;
    u[1] = -sVx[0].vtY * cz2z3 + sVx[2].vtY * cz0z3 - sVx[3].vtY * cz0z2;
    u[2] =  sVx[0].vtY * cz1z3 - sVx[1].vtY * cz0z3 + sVx[3].vtY * cz0z1;
    u[3] = -sVx[0].vtY * cz1z2 + sVx[1].vtY * cz0z2 - sVx[2].vtY * cz0z1;
    tr[ 0] = d * (dVx[3].vtX * u[3] + dVx[2].vtX * u[2] +
                  dVx[1].vtX * u[1] + dVx[0].vtX * u[0]);
    tr[ 4] = d * (dVx[3].vtY * u[3] + dVx[2].vtY * u[2] +
                  dVx[1].vtY * u[1] + dVx[0].vtY * u[0]);
    tr[ 8] = d * (dVx[3].vtZ * u[3] + dVx[2].vtZ * u[2] +
                  dVx[1].vtZ * u[1] + dVx[0].vtZ * u[0]);
    tr[12] = 0.0;
    u[0] = -sVx[1].vtX * cz2z3 + sVx[2].vtX * cz1z3 - sVx[3].vtX * cz1z2;
    u[1] =  sVx[0].vtX * cz2z3 - sVx[2].vtX * cz0z3 + sVx[3].vtX * cz0z2;
    u[2] = -sVx[0].vtX * cz1z3 + sVx[1].vtX * cz0z3 - sVx[3].vtX * cz0z1;
    u[3] =  sVx[0].vtX * cz1z2 - sVx[1].vtX * cz0z2 + sVx[2].vtX * cz0z1;
    tr[ 1] = d * (dVx[3].vtX * u[3] + dVx[2].vtX * u[2] +
                  dVx[1].vtX * u[1] + dVx[0].vtX * u[0]);
    tr[ 5] = d * (dVx[3].vtY * u[3] + dVx[2].vtY * u[2] +
                  dVx[1].vtY * u[1] + dVx[0].vtY * u[0]);
    tr[ 9] = d * (dVx[3].vtZ * u[3] + dVx[2].vtZ * u[2] +
                  dVx[1].vtZ * u[1] + dVx[0].vtZ * u[0]);
    tr[13] = 0.0;
    u[0] =  sVx[1].vtX * cy2y3 - sVx[2].vtX * cy1y3 + sVx[3].vtX * cy1y2;
    u[1] = -sVx[0].vtX * cy2y3 + sVx[2].vtX * cy0y3 - sVx[3].vtX * cy0y2;
    u[2] =  sVx[0].vtX * cy1y3 - sVx[1].vtX * cy0y3 + sVx[3].vtX * cy0y1;
    u[3] = -sVx[0].vtX * cy1y2 + sVx[1].vtX * cy0y2 - sVx[2].vtX * cy0y1;
    tr[ 2] = d * (dVx[3].vtX * u[3] + dVx[2].vtX * u[2] +
                  dVx[1].vtX * u[1] + dVx[0].vtX * u[0]);
    tr[ 6] = d * (dVx[3].vtY * u[3] + dVx[2].vtY * u[2] +
                  dVx[1].vtY * u[1] + dVx[0].vtY * u[0]);
    tr[10] = d * (dVx[3].vtZ * u[3] + dVx[2].vtZ * u[2] +
                  dVx[1].vtZ * u[1] + dVx[0].vtZ * u[0]);
    tr[14] = 0.0;
    u[0] = -sVx[1].vtX * cy2z3y3z2 + sVx[2].vtX * cy1z3y3z1 -
            sVx[3].vtX * cy1z2y2z1;
    u[1] =  sVx[0].vtX * cy2z3y3z2 - sVx[2].vtX * cy0z3y3z0 +
            sVx[3].vtX * cy0z2y2z0;
    u[2] = -sVx[0].vtX * cy1z3y3z1 + sVx[1].vtX * cy0z3y3z0 -
            sVx[3].vtX * cy0z1y1z0;
    u[3] =  sVx[0].vtX * cy1z2y2z1 - sVx[1].vtX * cy0z2y2z0 +
            sVx[2].vtX * cy0z1y1z0;
    tr[ 3] = d * (dVx[3].vtX * u[3] + dVx[2].vtX * u[2] +
                  dVx[1].vtX * u[1] + dVx[0].vtX * u[0]);
    tr[ 7] = d * (dVx[3].vtY * u[3] + dVx[2].vtY * u[2] +
                  dVx[1].vtY * u[1] + dVx[0].vtY * u[0]);
    tr[11] = d * (dVx[3].vtZ * u[3] + dVx[2].vtZ * u[2] +
                  dVx[1].vtZ * u[1] + dVx[0].vtZ * u[0]);
    tr[15] = 1.0;
  }
  return(squashed);
}

/*!
* \return      	Position of intersection.
* \ingroup      WlzGeometry
* \brief        Given a Woolz object and two vertices, finds the
*               position along a line segment between the two
*               vertices which is just inside/outside the boundary
*		of the object. The destination pointer is used to
*		return the status of the vertices, using the
*		following
*               code:
*                 0 - One of the given verticies was inside and the
*                     other outside,
*                 1 - Both the given verticies were inside,
*                 2 - Both the given verticies were outside.
*		This function assumes that the line segment only crosses
*		the object's boundary once.
*
*               Given the line segment \f$p_0\f$, \f$p_1\f$ any position
*               along the segment can be given by a parameter \f$\alpha\f$
*               (range [0-1]), where \f$p_x = p_0 + \alpha(p_1 - p_0)\f$.
* \param        obj                     Given object. Object must have a
*                                       valid domain.
* \param        p0                      First vertex.
* \param        p1                      Second vertex.
* \param        tol                     Acceptable placement error.
* \param	inside			Non-zero if the returned position
*					should be inside or on the
*					boundary, if zero it will be
*					outside or on the boundary.
* \param        dstStat                 Destination pointer for status,
*                                       may be NULL.
*/
WlzDVertex2	WlzGeomObjLineSegIntersect2D(WlzObject *obj,
					WlzDVertex2 p0, WlzDVertex2 p1,
					double tol, int inside, int *dstStat)
{
  int           s0,
                s1,
                s2,
                stat;
  double        dErr;
  WlzDVertex2   p2;

  tol *= tol;
  s0 = WlzInsideDomain(obj, 0.0, p0.vtY, p0.vtX, NULL);
  s1 = WlzInsideDomain(obj, 0.0, p1.vtY, p1.vtX, NULL);
  if(s0 != s1)
  {
    stat = 0;
    if(s0 != 0)
    {
      /* Ensure that p0 is outside the domain. */
      p2 = p0; p0 = p1; p1 = p2;
      s2 = s0; s0 = s1; s1 = s2;
    }
    do
    {
      /* Find midpoint of p0 and p1. */
      p2.vtX = 0.5 * (p0.vtX + p1.vtX);
      p2.vtY = 0.5 * (p0.vtY + p1.vtY);
      /* Check if the midpoint is within the object and update end points. */
      s2 = WlzInsideDomain(obj, 0.0, p2.vtY, p2.vtX, NULL);
      if(s2 != 0)
      {
        p1 = p2;
      }
      else
      {
        p0 = p2;
      }
      /* Check distance error. */
      WLZ_VTX_2_SUB(p2, p0, p1);
      dErr = WLZ_VTX_2_SQRLEN(p2);
    }
    while(dErr > tol);
  }
  else if(s0 != 0)
  {
    stat = 1;
  }
  else
  {
    stat = 2;
  }
  p2 = (inside)? p1: p0;
  if(dstStat)
  {
    *dstStat = stat;
  }
  return(p2);
}

/*!
* \return      	Position of intersection.
* \ingroup      WlzGeometry
* \brief        Given a Woolz object and two vertices, finds the
*               position along a line segment between the two
*               vertices which is just inside/outside the boundary
*		of the object. The destination pointer is used to
*		return the status of the vertices, using the
*		following
*               code:
*                 0 - One of the given verticies was inside and the
*                     other outside,
*                 1 - Both the given verticies were inside,
*                 2 - Both the given verticies were outside.
*		This function assumes that the line segment only crosses
*		the object's boundary once.
*
*               Given the line segment \f$p_0\f$, \f$p_1\f$ any position
*               along the segment can be given by a parameter \f$\alpha\f$
*               (range [0-1]), where \f$p_x = p_0 + \alpha(p_1 - p_0)\f$.
* \param        obj                     Given object. Object must have a
*                                       valid domain.
* \param        p0                      First vertex.
* \param        p1                      Second vertex.
* \param        tol                     Acceptable placement error.
* \param	inside			Non-zero if the returned position
*					should be inside or on the
*					boundary, if zero it will be
*					outside or on the boundary.
* \param        dstStat                 Destination pointer for status,
*                                       may be NULL.
*/
WlzDVertex3	WlzGeomObjLineSegIntersect3D(WlzObject *obj,
					WlzDVertex3 p0, WlzDVertex3 p1,
					double tol, int inside, int *dstStat)
{
  int           s0,
                s1,
                s2,
                stat;
  double        dErr;
  WlzDVertex3   p2;

  tol *= tol;
  s0 = WlzInsideDomain(obj, p0.vtZ, p0.vtY, p0.vtX, NULL);
  s1 = WlzInsideDomain(obj, p1.vtZ, p1.vtY, p1.vtX, NULL);
  if(s0 != s1)
  {
    stat = 0;
    if(s0 != 0)
    {
      /* Ensure that p0 is outside the domain. */
      p2 = p0; p0 = p1; p1 = p2;
      s2 = s0; s0 = s1; s1 = s2;
    }
    do
    {
      /* Find midpoint of p0 and p1. */
      p2.vtX = 0.5 * (p0.vtX + p1.vtX);
      p2.vtY = 0.5 * (p0.vtY + p1.vtY);
      p2.vtZ = 0.5 * (p0.vtZ + p1.vtZ);
      /* Check if the midpoint is within the object and update end points. */
      s2 = WlzInsideDomain(obj, p2.vtZ, p2.vtY, p2.vtX, NULL);
      if(s2 != 0)
      {
        p1 = p2;
      }
      else
      {
        p0 = p2;
      }
      /* Check distance error. */
      WLZ_VTX_3_SUB(p2, p0, p1);
      dErr = WLZ_VTX_3_SQRLEN(p2);
    }
    while(dErr > tol);
  }
  else if(s0 != 0)
  {
    stat = 1;
  }
  else
  {
    stat = 2;
  }
  p2 = (inside)? p1: p0;
  if(dstStat)
  {
    *dstStat = stat;
  }
  return(p2);
}

/*!
* \return	Maximum diameter sphere inscribed within the tetrahedron.
* 		
* \ingroup	WlzGeometry
* \brief	Given the coordinates of the four vertices of a tetrahedron
* 		the function computes the maximum diameter of an inscribed
* 		sphere.
* 		
* 		Diameter of the (maximum) inscribed sphere \f$d\f$ is given
* 		by:
* 		\f[
                    d = \frac{6V}{A}
 	        \f]
*               where \f$V\f$ is the volume of the tetrahedron and \f$A\f$
*               is it's surface area. "Encyclopedia of Mathematics"
*               ISBN 1402006098 (http://eom.springer.de).
* \param	vx0			First vertex of tetrahedron.
* \param	vx1			Second vertex of tetrahedron.
* \param	vx2			Third vertex of tetrahedron.
* \param	vx3			Forth vertex of tetrahedron.
*/
double		WlzGeomTetraInSphereDiam(WlzDVertex3 vx0, WlzDVertex3 vx1,
				         WlzDVertex3 vx2, WlzDVertex3 vx3)
{
  double	diam = 0.0,
		a0,
		a1,
		a2,
		a3,
  		vol6;
  const double  tol = ALG_DBL_TOLLERANCE;

  vol6 = fabs(WlzGeomTetraSnVolume6(vx0, vx1, vx2, vx3));
  if(vol6 > tol)
  {
    a0 = WlzGeomTriangleArea2Sq3(vx0, vx1, vx2);
    a1 = WlzGeomTriangleArea2Sq3(vx1, vx2, vx3);
    a2 = WlzGeomTriangleArea2Sq3(vx2, vx3, vx0);
    a3 = WlzGeomTriangleArea2Sq3(vx3, vx0, vx1);
    diam = (2.0 * vol6) /
           (sqrt(fabs(a0)) + sqrt(fabs(a1)) +
	    sqrt(fabs(a2)) + sqrt(fabs(a3)));
  }
  return(diam);
}

/*!
* \return	Maximum diameter sphere inscribed within the regular
* 		tetrahedron.
* 		
* \ingroup	WlzGeometry
* \brief	Given the side length of a regular tetrahedron this function
* 		computes the maximum diameter of an inscribed sphere.
* 		
* 		Diameter of the (maximum) inscribed sphere \f$d\f$ is given
* 		by:
* 		\f[
                    d = \frac{6V}{A}
 	        \f]
*               where \f$V\f$ is the volume of the tetrahedron and \f$A\f$
*               is it's surface area. "Encyclopedia of Mathematics"
*               ISBN 1402006098 (http://eom.springer.de).
*               Standard formulae are used for the area and volume.
* \param	side			Length of an edge.
*/
double		WlzGeomTetraInSphereRegDiam(double side)
{
  double	diam = 0.0,
		area,
		sideSq,
  		vol12;
  const double  tol = ALG_DBL_TOLLERANCE;

  sideSq = side * side;
  if(sideSq > tol)
  {
    area = sideSq * ALG_M_SQRT3;
    vol12 = side * sideSq * ALG_M_SQRT2;
    diam = 0.5 * vol12 / area;
  }
  return(diam);
}

/*!
* \return       Angle, range \f$[0 - 2\pi]\f$.
* \ingroup      WlzGeometry
* \brief        Computes the angle of a ray from the origin to the
*               destination vertex along with the length of the ray.
*               Angles are:
*               \verbatim
                      (y)
                       ^  pi/2
                       |
                       |
                 pi    |
                <------+-------> (x)
                       |       0
                       |
                  3pi/2|
                       V
                \endverbatim
* \param        org                     Position of the origin.
* \param        dst                     Position of the destination.
* \param        dstRad                  Destination pointer for the length
*                                       of the ray, may be NULL.
*/
double          WlzGeomPolar2D(WlzDVertex2 org, WlzDVertex2 dst,
                                double *dstRad)
{
  double        ang;
  WlzDVertex2   del;

  WLZ_VTX_2_SUB(del, dst, org);
  ang = atan2(del.vtY, del.vtX);
  if(ang < 0.0)
  {
    ang += 2.0 * ALG_M_PI;
  }
  if(dstRad)
  {
    *dstRad = WLZ_VTX_2_LENGTH(del);
  }
  return(ang);
}


/*!
* \return	Cosine of angle between line segments (v0, v1) and (v1, v2)
* 		with value in the range [0-1].
* \ingroup	WlzGeometry
* \brief	Computes the cosine of angle between line segments (v0, v1)
* 		and (v1, v2). If any of these vertices are coincident then
* 		zero is returned.
* \param	v0			First vertex.
* \param	v1			Second vertex (the common one).
* \param	v2			Third vertex.
*/
double		WlzGeomCos3V(WlzDVertex2 v0, WlzDVertex2 v1, WlzDVertex2 v2)
{
  double	c = 0.0,
  		l0,
  		l1,
		l2;
  WlzDVertex2	s0,
  		s1,
		s2;
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_2_SUB(s0, v1, v2);
  WLZ_VTX_2_SUB(s1, v2, v0);
  WLZ_VTX_2_SUB(s2, v0, v1);
  l0 = WLZ_VTX_2_SQRLEN(s0);
  l1 = WLZ_VTX_2_SQRLEN(s1);
  l2 = WLZ_VTX_2_SQRLEN(s2);
  if((l0 > tol) && (l2 > tol))
  {
    c = (0.5 * (l0 + l2 - l1)) / sqrt(l0 * l2);
  }
  return(c);
}

/*!
* \return	Non-zero value only if test vertex is on the line segment.
* \ingroup	WlzGeometry
* \brief	Tests whether the given test vertex is on the given line
* 		segment.
* 		If all three vertices are coincident then the test vertex
* 		is considered to line on the line segment.
* \param	tst				Test vertex.
* \param	seg0				First vertex of line segment.
* \param	seg1				Second vertex of line segment.
* \param	tol				Tollerance.
*/
int             WlzGeomVtxOnLineSegment2D(WlzDVertex2 tst,
                                          WlzDVertex2 seg0, WlzDVertex2 seg1,
                                          double tol)
{
  int           onSeg = 0;
  double        mS,
                mT,
                tolSq;
  WlzDVertex2   delS,
                delT;
  WlzDBox2	box;

  /* 1. Simple in box test. */
  box.xMin = box.xMax = seg0.vtX;
  box.yMin = box.yMax = seg0.vtY;
  if(seg1.vtX < box.xMin)
  {
    box.xMin = seg1.vtX;
  }
  else if(seg1.vtX > box.xMax)
  {
    box.xMax = seg1.vtX;
  }
  if(seg1.vtY < box.yMin)
  {
    box.yMin = seg1.vtY;
  }
  else if(seg1.vtY > box.yMax)
  {
    box.yMax = seg1.vtY;
  }
  if(((tst.vtX - box.xMin) > -tol) && ((tst.vtY - box.yMin) > -tol) &&
     ((box.xMax - tst.vtX) > -tol) && ((box.yMax - tst.vtY) > -tol))
  {
    /* 2. Test vertex coincident with either segment vertex. */
    tolSq = tol * tol;
    if(WlzGeomVtxEqual2D(seg0, tst, tolSq) ||
       WlzGeomVtxEqual2D(seg1, tst, tolSq))
    {
      onSeg = 1;
    }
    else
    {
      /* 3. Test gradients of the lines (seg0, seg1) and (seg0, tst) are
       * equal. */
      WLZ_VTX_2_SUB(delS, seg1, seg0);
      WLZ_VTX_2_SUB(delT, tst, seg0);
      if(delS.vtX - delS.vtY > 0)
      {
	if(delT.vtX - delT.vtY > -tol)
	{
	  mS = delS.vtY / delS.vtX;
	  mT = delT.vtY / delT.vtX;
	  if(fabs(mS - mT) < tol)
	  {
	    onSeg = 1;
	  }
	}
      }
      else
      {
	if(delT.vtY - delT.vtX > -tol)
	{
	  mS = delS.vtX / delS.vtY;
	  mT = delT.vtX / delT.vtY;
	  if(fabs(mS - mT) < tol)
	  {
	    onSeg = 1;
	  }
	}
      }
    }
  }
  return(onSeg);
}

/*!
* \return	Integer code corresponding to position of the test
* 		vertex with respect to the line segment:
*		<ul>
*		  <li>0 no intersection.</li>
*		  <li>1 intersection at an end point.</li>
*                 <li>2 intersection at a single point (not end points).</li>
*		</ul>
* \ingroup	WlzGeometry
* \brief	Tests whether the given test vertex is on the given line
* 		segment. If all three vertices are coincident then the
* 		test vertex is considered to be coincident with an end
* 		point on the line segment.
*
* 		Consider a line segment from a vertex at \f$\mathbf{p_0}\f$
* 		to another vertex at \f$\mathbf{p_1}\f$, with a third test
* 		vertex at \f$\mathbf{p_x}\f$, the shortest path from
* 		\f$\mathbf{p_x}\f$ to the line segment will be perpendicular
* 		to the line segment. Let the position of the intersection
* 		of this perpendicular with the line segment be at
* 		\f$\mathbf{p}\f$, with distance \f$d\f$ from \mathbf{p_x},
* 		then:
* 		\f[
		\mathbf{p} = \mathbf{p_0} + s (\mathbf{p_1} - \mathbf{p_0},
		\f]
* 		\f[
		d^2 = \| (\mathbf{p_0} - \mathbf{p_x}) +
		         (\mathbf{p} - \mathbf{p_0}) \|^2,
		\f]
* 		\f[
                s = \frac{ (\mathbf{p_0} - \mathbf{p_x}) \cdot
		           (\mathbf{p_1} - \mathbf{p_0})}
		         {\| \mathbf{p_1} - \mathbf{p_0} \|^2},
		\f]
* 		and after substitution:
* 		\f[
		d^2 = \frac{\| \mathbf{p_0} - \mathbf{p_x} \|^2
		            \| \mathbf{p_1} - \mathbf{p_0} \|^2 -
			    \left[(\mathbf{p_1} - \mathbf{p_0})\right]^2 }
		           {\| \mathbf{p_1} - \mathbf{p_0} \|^2},
		\f]
* 		Observing that the numerator is a vector quad product
* 		\f[
		\mathbf{A} \times \mathbf{B} = 
		\|\mathbf{A}\|^2\|\mathbf{A}\|^2 -
		(\mathbf{A} \cdot \mathbf{Ba})^2
		\f]
		gives
* 		\f[
		d^2 = \frac{\| \mathbf{p_1} - \mathbf{p_0}  \times
		               \mathbf{p_0} - \mathbf{p_x} \|^2 }
		           {\| \mathbf{p_1} - \mathbf{p_0} \|^2},

		\f]
*		Obviously if \f$\mathbf{p_x}\f$ is on the line segment then
*		\f$d^2 = 0\f$.
*		ALG_DBL_TOLLERANCE as a tolerance for squared distances in this
*		function.
* \param	pX				Test vertex.
* \param	p0				First vertex of line segment.
* \param	p1				Second vertex of line segment.
* \param	dstN				Destination pointer for the
* 						intersection, may be NULL.
*/
int             WlzGeomVtxOnLineSegment3D(WlzDVertex3 pX,
                                        WlzDVertex3 p0, WlzDVertex3 p1,
                                        WlzDVertex3 *dstN)
{
  int           isn = 0;
  double	s,
  		tSq,
  		l10Sq;
  WlzDVertex3	pT,
  		p10,
  		p0X;
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_3_SUB(p10, p1, p0);
  WLZ_VTX_3_SUB(p0X, p0, pX);
  l10Sq = WLZ_VTX_3_SQRLEN(p10);
  if(l10Sq < tol)
  {
    tSq = WLZ_VTX_3_SQRLEN(p0X);
    if(tSq < tol)
    {
      isn = 1;
      if(dstN)
      {
        *dstN = p0;
      }
    }
  }
  else
  {
    WLZ_VTX_3_CROSS(pT, p10, p0X);
    tSq = WLZ_VTX_3_SQRLEN(pT) / l10Sq;
    if(tSq < tol)
    {
      s = WLZ_VTX_3_DOT(p0X, p10) / l10Sq;
      if((1.0 - s) * (1.0 - s) < tol)
      {
        isn = 1;
	if(dstN)
	{
	  *dstN = p1;
	}
      }
      else if((s > 0) && (s < 1.0))
      {
        isn = 2;
	if(dstN)
	{
	  WLZ_VTX_3_SCALE_ADD(*dstN, p10, s, p0);
	}
      }
    }
  }
  return(isn);
}

/*!
* \return	The arc length.
* \ingroup	WlzGeometry
* \brief	Computes the arc length from a to b traveling CCW on a
* 		circle with centre c.
* \param	a			Start point.
* \param	b			End point.
* \param	c			Cirecle centre.
*/
double		WlzGeomArcLength2D(WlzDVertex2 a, WlzDVertex2 b, WlzDVertex2 c)
{
  int		qa,
  		qb;
  double	ang,
  		angA,
  		angB,
		chordSq,
  		r,
		rSq,
		len = 0.0;
  WlzDVertex2	p,
  		t;
  int		quadTbl[4] = {2, 3, 1, 0};    /* o--> x  2|3
                                               * |       -+-
					       * v y     1|0
					       */
  double	quadEndX[4] = { 0.0, -1.0,  0.0,  1.0},
  		quadEndY[4] = { 1.0,  0.0, -1.0,  0.0};
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_2_SUB(t, a, c);
  rSq = WLZ_VTX_2_SQRLEN(t);
  if(rSq > tol)
  {
    r = sqrt(rSq);
    /* Compute the quadrants that contain points a and b. */
    qa = quadTbl[((a.vtY > c.vtY) << 1) | (a.vtX > c.vtX)];
    qb = quadTbl[((b.vtY > c.vtY) << 1) | (b.vtX > c.vtX)];
    if(qa == qb)
    {
      /* Compute angle from point a to b. */
      WLZ_VTX_2_SUB(t, a, b);
      chordSq = WLZ_VTX_2_SQRLEN(t);
      ang = acos(1.0 - (0.5 * chordSq / rSq));
    }
    else
    {
      /* Compute angle at centre from point a to the end of it's quadrant. */
      p.vtX = c.vtX + r * quadEndX[qa];
      p.vtY = c.vtY + r * quadEndY[qa];
      WLZ_VTX_2_SUB(t, a, p);
      chordSq = WLZ_VTX_2_SQRLEN(t);
      angA = acos(1.0 - (0.5 * chordSq / rSq));
      /* Compute angle at centre from point b to the start of it's quadrant. */
      p.vtX = c.vtX + r * quadEndX[(qb + 3) % 4];
      p.vtY = c.vtY + r * quadEndY[(qb + 3) % 4];
      WLZ_VTX_2_SUB(t, b, p);
      chordSq = WLZ_VTX_2_SQRLEN(t);
      angB = acos(1.0 - (0.5 * chordSq / rSq));
      /* Compute total CCW arc length from a to b. */
      ang = angA + angB + ((3 + qb - qa) * ALG_M_PI_2);
      while(ang > 2 * ALG_M_PI)
      {
        ang -= 2 * ALG_M_PI;
      }
    }
    len = r * ang;
  }
  return(len);
}

/*!
* \return	Non-zero if \f$S\f$ and \f$T\f$ are coincident.
* \ingroup	WlzGeometry
* \brief 	Computes the coordinates of vertices that may be used
* 		to draw a wide rectangle.
*
* 		Given two vertices \f$S\f$, \f$T\f$ which define a line
* 		segment and width \f$w\f$ perpendicular to the line segent,
* 		this function computes the coordinates of the vertices
* 		\f$V_i\f$, \f$i\in[0-3]\f$ of the rectangle. The vertices
* 		are sorted such that the rectangle may be drawn using
* 		the line segmants:
* 		\f$(V_0,V_1)\f$, \f$(V_1,V_2)\f$, \f$(V_2,V_3)\f$,
* 		\f$(V_3,V_0)\f$.
* 		Given the line segment \f$(S,T)\f$ it can be shown that
* 		the vertices which share a line segment passing through
* 		\f$S\f$ are:
* 		\f[
 		V_i = S + (-\frac{r (T_x - S_x)}{l},  \frac{r (T_y - S_y)}{l})
		\f]
* 		and
* 		\f[
 		V_j = S + ( \frac{r (T_x - S_x)}{l}, -\frac{r (T_x - S_y)}{l})
		\f]
*		where \f$l = ||T - S||\f$ and \f$r = w / 2\f$.
*		Should  \f$S\f$ and \f$T\f$ be coincident then the
*		destination rectangle vertices are left unmodified..
* \param	s			First vertex of line segment.
* \param	t			Second vertex of line segment.
* \param	w			line width.
* \param	v			Destination pointer for the four
* 					vertices of the rectangle.
*/
int		WlzGeomRectFromWideLine(WlzDVertex2 s, WlzDVertex2 t,
				        double w, WlzDVertex2 *v)
{
  int		coincident = 1;
  double	len;
  WlzDVertex2	del,
  		f;
  const double  tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_2_SUB(del, t, s);
  len = WLZ_VTX_2_LENGTH(del);
  if(len > tol)
  {
    coincident = 0;
    f.vtX =   w * del.vtY / (2.0 * len);
    f.vtY = -(w * del.vtX / (2.0 * len));
    WLZ_VTX_2_ADD(v[0], s, f);
    WLZ_VTX_2_SUB(v[1], s, f);
    WLZ_VTX_2_SUB(v[2], t, f);
    WLZ_VTX_2_ADD(v[3], t, f);
  }
  return(coincident);
}

/*!
* \return	Intersection of line with plane.
* \ingroup	WlzGeometry
* \brief	Computes the intersection of a line with a plane.
* 		
* 		Computes the intersection of line (vector \f$\mathbf{v}\f$)
* 		through an off plane vertex (\f$\mathbf{p_3}\f$) with a plane.
* 		The plane is defined by three on plane vertices
* 		\f$(\mathbf{p_0}, \mathbf{p_1}, \mathbf{p_2})\f$.
* 		\f$q\f$ the vertex at the intersection is found by solving
* 		the vector equations
* 		\f[
		(\mathbf{p_1} - \mathbf{p_0}) \times
		(\mathbf{p_2} - \mathbf{p_0}) .
		(\mathbf{q}   - \mathbf{p_0}) = \mathbf{0},
		\mathbf{n} \times (\mathbf{p_3} - \mathbf{q}) = \mathbf{0}
		\f]
*		These are equivalunt to solving
*		\f[
                \left|\begin{array}{ccc}
                p_{1x} - p_{0x} & p_{1y} - p_{0y} & p_{1z} - p_{0z} \\
                p_{2x} - p_{0x} & p_{2y} - p_{0y} & p_{2z} - p_{0z} \\
                q_x    - p_{0x} & q_y    - p_{0y} & q_z    - p_{0z}
                \end{array}\right| = 0,
                \left|\begin{array}{ccc}
                \mathbf{i}   & \mathbf{j}   & \mathbf{k} \\
                n_x          & n_y          & n_z \\
                p_{3x} - q_x & p_{3y} - q_y & p_{3z} - q_z
                \end{array}\right| = 0
		\f]
* \param	v			Unit vector for ray.
* \param	p0			First point on the plane.
* \param	p1			Second point on the plane.
* \param	p2			Third point on the plane.
* \param	p3			Off plane point that ray passes
* 					through.
* \param	dstPar			Destination value set to 1 if the
* 					vector is parrallel to the plane,
* 					must not be NULL.
*/
WlzDVertex3	WlzGeomLinePlaneIntersection(WlzDVertex3 v,
					    WlzDVertex3 p0, WlzDVertex3 p1,
					    WlzDVertex3 p2, WlzDVertex3 p3,
					    int *dstPar)
{
  double	a,
  		b,
		c,
		d,
		e,
		f,
		g,
		h;
  WlzDVertex3	p1r,
		q;
  const double  tol = ALG_DBL_TOLLERANCE;
  		
  q = p0;
  WLZ_VTX_3_SUB(p1r, p1, p0);
  a = p0.vtY * p1.vtZ - p0.vtZ * p1.vtY +
      p1r.vtY * p2.vtZ - p1r.vtZ * p2.vtY;
  b = p0.vtZ * p1.vtX - p0.vtX * p1.vtZ +
      p1r.vtZ * p2.vtX - p1r.vtX * p2.vtZ;
  c = p0.vtX * p1.vtY - p0.vtY * p1.vtX +
      p1r.vtX * p2.vtY - p1r.vtY * p2.vtX;
  d = -(v.vtX * a + v.vtY * b + v.vtZ * c);
  if(d * d < tol)
  {
    *dstPar = 1;
  }
  else
  {
    *dstPar = 0;
    e = (p0.vtY * p1.vtX - p0.vtX * p1.vtY) * p2.vtZ +
	(p0.vtX * p1.vtZ - p0.vtZ * p1.vtX) * p2.vtY +
	(p0.vtZ * p1.vtY - p0.vtY * p1.vtZ) * p2.vtX;
    f = p0.vtZ * p1.vtY - p0.vtY * p1.vtZ +
	p1r.vtZ * p2.vtY - p1r.vtY * p2.vtZ;
    g = p0.vtY * p1.vtX - p0.vtX * p1.vtY +
	p1r.vtY * p2.vtX - p1r.vtX * p2.vtY;
    h = p0.vtX * p1.vtZ - p0.vtZ * p1.vtX +
	p1r.vtX * p2.vtZ - p1r.vtZ * p2.vtX;
    q.vtX = (v.vtX * (b * p3.vtY + c * p3.vtZ + e) +
	     v.vtY * h * p3.vtX +
	     v.vtZ * g * p3.vtX)/ d;
    q.vtY = (v.vtX * f * p3.vtY +
	     v.vtY * (c * p3.vtZ + a * p3.vtX + e) +
	     v.vtZ * g * p3.vtY)/ d;
    q.vtZ = (v.vtX * f * p3.vtZ +
	     v.vtY * h * p3.vtZ +
	     v.vtZ * (a * p3.vtX + b * p3.vtY + e))/ d;
  }
  return(q);
}

/*!
* \return	Value indicating the position of the vertex with respect
*               to the triangle in 3D:
*		<ul>
*		  <li>0 if the vertex is outside the triangle.</li>
*		  <li>1 if the vertex is on an edge of the triangle
*		        or line passes through triangle in it's plane.</li>
*		  <li>2 if the vertex is inside the triangle.</li>
*		</ul>
* \ingroup	WlzGeometry
* \brief	Tests to set if a line directed from a given origin
* 		intersects a triangle in 3D space. This function is
* 		based on the algorithm: Tomas Moller and Ben Trumbore,
* 		"Fast, Minimum Storage Ray/Triangle Intersection",
* 		Journal of Graphics Tools, 1997(2), pp 25--30.
*
* 		Given a parameterised line
* 		\f[
 		R(t) = O + t D
		\f]
*		and a triangle specified by it's vertices
*		\f$(V_0, V_1, V_2)\f$, the point of intersection may
*		be written in terms of the barycentric coordinates
*		\f$u, v\f$
* 		\f[
                O + t D = T(u, v) = (1 - u - v) V_0 + u V_1 + v V_2
                \f]
*		Solving for \f$t\f$, \f$u\f$ and \f$v\f$ gives
*		\f[
		\left[ \begin{array}{c}
		       t \\
		       u \\
		       v
		       \end{array}
		\right] = 
		\frac{1}{P \cdot E_1}
		\left[ \begin{array}{c}

		       Q \cdot E_2 \\
		       P \cdot t \\
		       Q \cdot D
		       \end{array}
		\right]
		\f]
*		where \f$E_1 = V_1 - V_0\f$, \f$E_2 = V_2 - V_0\f$
*		\f$t = O - V_0\f$, \f$P = D \times E_2\f$ and
*		\f$Q = T \times E_1\f$.
*
*		The \f$t, u\f$ and \f$v\f$ are only set if the line passes
*		through the triangle.
*
* \param	org			Line origin, \f$O\f$.
* \param	dir			Line direction, \f$D\f$ (does not
* 					need to be a unit vector).
* \param	v0			First vertex on triangle.
* \param	v1			Second vertex on the triangle
* \param	v2			Third vertex on the triangle
* \param	dstPar			Destination pointer for flag set to
* 					1 if the vector is parrallel to the
* 					plane, may be NULL.
* \param	dstT			Destination pointer for t parameter,
* 					may be NULL.
* \param	dstU			Destination pointer for u parameter,
* 					may be NULL.
* \param	dstV			Destination pointer for v parameter,
* 					may be NULL.
*/
int		WlzGeomLineTriangleIntersect3D(WlzDVertex3 org, WlzDVertex3 dir,
				WlzDVertex3 v0, WlzDVertex3 v1, WlzDVertex3 v2,
				int *dstPar, double *dstT,
				double *dstU, double *dstV)
{
  int		par,
  		isn = 0;
  double	det,
		u = 0,
		v = 0,
		w;
  WlzDVertex3	e1,
  		e2,
		p,
		q,
		m;
  const double  tol = ALG_DBL_TOLLERANCE;

  /* Find vectors for two edges sharing v0. */
  WLZ_VTX_3_SUB(e1, v1, v0);
  WLZ_VTX_3_SUB(e2, v2, v0);
  /* Compute determinant, if near zero line passes through plane of triangle. */
  WLZ_VTX_3_CROSS(p, dir, e2);
  det = WLZ_VTX_3_DOT(e1, p);
  if(det * det < tol)
  {
    par = 1;
    isn = WlzGeomLineTriangleIntersectEdge3D(org, dir, v0, v1, v2);
  }
  else
  {
    par = 0;
    det = 1.0 / det;
    /* Calculate distance from v0 to org. */
    WLZ_VTX_3_SUB(m, org, v0);
    /* Calculate barycentric u parameter and test bounds. */
    u = WLZ_VTX_3_DOT(m, p) * det;
    if((u > -tol) && (u - 1.0 < tol))
    {
      /* Calculate barycentric v parameter and test bounds. */
      WLZ_VTX_3_CROSS(q, m, e1);
      v = WLZ_VTX_3_DOT(dir, q) * det;
      if((v > -tol) && (v - 1.0 < tol))
      {
        /* Test bounds of the trird barycentric coordinate. */
        w = 1.0 - (u + v);
	if((w > -tol) && (w - 1.0 < tol))
	{
	  if((u < tol) || (v < tol) || (w < tol))
	  {
	    isn = 1;
	  }
	  else
	  {
	    isn = 2;
	  }
	  if(dstU)
	  {
	    *dstU = u;
	  }
	  if(dstV)
	  {
	    *dstV = v;
	  }
	  if(dstT)
	  {
	    *dstT = WLZ_VTX_3_DOT(e2, q) * det;
	  }
	}
      }
    }
  }
  if(dstPar)
  {
    *dstPar = par;
  }
  return(isn);
}

/*!
* \return	1 for intersection 0 for no intersection.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection of the given line with the edges
* 		of a given triangle.
* \param	org			Line origin, \f$O\f$.
* \param	dir			Line direction, \f$D\f$ (does not
* 					need to be a unit vector).
* \param	v0			First vertex on triangle.
* \param	v1			Second vertex on the triangle
* \param	v2			Third vertex on the triangle
*/
static int	WlzGeomLineTriangleIntersectEdge3D(WlzDVertex3 org,
				WlzDVertex3 dir, WlzDVertex3 v0,
				WlzDVertex3 v1, WlzDVertex3 v2)
{
  int		isn;

  isn = (WlzGeomLineLineSegmentIntersect3D(org, dir, v0, v1, NULL) > 0) ||
        (WlzGeomLineLineSegmentIntersect3D(org, dir, v0, v2, NULL) > 0);
  return(isn);
}

/*!
* \return	Integer value which classifies the intersection of
*		the line with the line line segment. Values are:
*		<ul>
*		  <li>0 no intersection.</li>
*		  <li>1 intersection along line segment, all points
*			on the line segment are coincident.</li>
*                 <li>2 intersection at end points.</li>
*                 <li>3 intersection at a single point (not end points).</li>
*		</ul>
* \ingroup	WlzGeometry
* \brief	Tests to see if the two given line segment is intersected
* 		by the given line using the ALG_DBL_TOLLERANCE tollerance value.
* 		The line is a line which passes through the given point
* 		to infinity (on both sides) with the given direction.
*
* 		Given a line parameterised by \f$t\f$:
* 		\f[
		\mathbf{x} = \mathbf{R_0} + t \mathbf{R_d}
		\f]
* 		and a line segment parameterised by \f$s\f$:
* 		\f[
		\mathbf{x} = \mathbf{P_0} + s (\mathbf{P_1} - \mathbf{P_0})
		\f]
*		their intersection at \f$\mathbf{x}\f$ is given by
* 		\f[
 		\mathbf{P_0} + s (\mathbf{P_1} - \mathbf{P_0}) =
		\mathbf{R_0} + t \mathbf{R_d}
 		\f]
* 		giving
* 		\f[
		s (\mathbf{P_1} - \mathbf{P_0}) \times \mathbf{R_d} =
		  (\mathbf{R_0} - \mathbf{P_0} \times \mathbf{R_d}
		t \mathbf{R_d} \times (\mathbf{P_1} - \mathbf{P_0}) =
		  (\mathbf{P_0} - \mathbf{R_0} \times
		  (\mathbf{P_1} - \mathbf{P_0})
 		\f]
* 		\f[
 		s ((\mathbf{P_1} - \mathbf{P_0}) \times \mathbf{R_d}) \cdot
 		  ((\mathbf{P_1} - \mathbf{P_0}) \times \mathbf{R_d}) =
		  ((\mathbf{R_0} - \mathbf{P_0}) \times \mathbf{R_d}) \cdot
		  ((\mathbf{P_1} - \mathbf{P_0}) \times \mathbf{R_d})
 		t (\mathbf{R_d} \times (\mathbf{P_1} - \mathbf{P_0})) \cdot
 		  (\mathbf{R_d} \times (\mathbf{P_1} - \mathbf{P_0})) =
		  ((\mathbf{P_0} - \mathbf{R_0} \times
		   (\mathbf{P_1} - \mathbf{P_0})) \cdot
		  (\mathbf{R_d} \times (\mathbf{P_1} - \mathbf{P_0}))
		\f]
* 		\f[
  		s = \frac{\mbox{det}(\mathbf{R_0} - \mathbf{P_0},
  		                     \mathbf{R_d},
  		                     (\mathbf{P_1} - \mathbf{P_0}) \times
  		                     \mathbf{R_d}}
  		         {\|(\mathbf{P_1} - \mathbf{P_0}) \times
  		            \mathbf{R_d}\|^2}
  		t = \frac{\mbox{det}(\mathbf{R_0} - \mathbf{P_0},
  		                     \mathbf{R_d},
  		                     (\mathbf{P_1} - \mathbf{P_0}) \times
  		                     \mathbf{R_d}}
  		         {\|(\mathbf{P_1} - \mathbf{P_0}) \times
  		            \mathbf{R_d}\|^2}
		\f]
*		If the denominator
*	      \f$\|(\mathbf{P_1} - \mathbf{P_0}) \times \mathbf{R_d}\|^2 = 0\f$
*		then the line and the line segment are parrallel, provided
*		that \f$\mathbf{P_1} \neq \mathbf{P_0}\f$ and
*		\f$\mathbf{R_d} \neq \mathbf{0}\f$.
*		The line segment is intersected by the line if
*		\f$s\f$ is in the range [0-1].
*		The point of intersection \f$\mathbf{x}\f$ is
*		\f[
		\mathbf{x} = \mathbf{P_0} + s (\mathbf{P_1} - \mathbf{P_0})
		\f]
*		Special cases are considered for
*		\f$\mathbf{R_0} = \mathbf{P_0},\mathbf{P_1}\f$,
*		\f$\mathbf{R_0} = 2 \mathbf{P_0} - \mathbf{P_1}\f$ and
*		\f$\mathbf{R_d} = \alpha (\mathbf{P_1} - \mathbf{P_0})\f$.
* \param	r0			A vertex on the line.
* \param	rD			Direction of the line.
* \param	p0			First end vertex of the line segment.
* \param	p1			Second end vertex of the line segment.
* \param	dstN			Destination ptr for intersection
*					vertex, may be NULL. The intersection
*					value will be set if the ray and line
*					segment intersect at a single point.
*/
int		WlzGeomLineLineSegmentIntersect3D(WlzDVertex3 r0,
				WlzDVertex3 rD,
				WlzDVertex3 p0, WlzDVertex3 p1,
				WlzDVertex3 *dstN)
{
  int		isn = 0;
  double	den,
  		lSq,
		s;
  WlzDVertex3	p,
		r,
  		v;
  const double  tol = ALG_DBL_TOLLERANCE;

  /* TODO Test this function. */
  WLZ_VTX_3_SUB(p, p1, p0);
  lSq = WLZ_VTX_3_SQRLEN(p);
  if(lSq < tol)
  {
    /* Given line segment is a single point at p0. */
    if(WlzGeomVtxOnLine3D(p0, r0, rD) != 0)
    {
      isn = 1;
    }
  }
  else
  {
    lSq = WLZ_VTX_3_SQRLEN(rD);
    if(lSq < tol)
    {
      /* Line is a single point at r0. */
      switch(WlzGeomVtxOnLineSegment3D(r0, p0, p1, dstN))
      {
        case 1:
	  isn = 2;
	  break;
	case 2:
	  isn = 3;
	  break;
        default:
	  break;
      }
    }
    else
    {
      WLZ_VTX_3_CROSS(v, p, rD);
      den = WLZ_VTX_3_SQRLEN(v);
      if(den < tol)
      {
	/* Line and line segment are parralel. */
	if(WlzGeomVtxOnLine3D(p0, r0, rD) > 0)
	{
	  isn = 1;
	}
      }
      else
      {
	WLZ_VTX_3_SUB(p, p1, p0);
	WLZ_VTX_3_SUB(r, r0, p0);
	s = (r.vtX * (rD.vtY * v.vtZ - rD.vtZ * v.vtY) +
	     r.vtY * (rD.vtZ * v.vtX - rD.vtX * v.vtZ) +
	     r.vtZ * (rD.vtX * v.vtY - rD.vtY * v.vtX)) / den;
	if((s * s) < tol)
	{
	  isn = 2;
	  if(dstN)
	  {
	    *dstN = p0;
	  }
	}
	else if((1 - s) * (1 - s) < tol)
	{
	  isn = 2;
	  if(dstN)
	  {
	    *dstN = p1;
	  }
	}
	else if((s > 0) && (s < 1.0))
	{
	  isn = 3;
	  if(dstN)
	  {
	    WLZ_VTX_3_SCALE_ADD(*dstN, p, s, p0);
	  }
	}
      }
    }
  }
  return(isn);
}

/*!
* \return	If the vertex is on the ray 1, otherwise 0.
* \ingroup	WlzGeometry
* \brief	Tests whether a vertex is a line.
*
* 		Tests whether the given vertex \f$p_0\f$ is on the
* 		given line
* 		\f$\mathbf{x} = \mathbf{r_0} + \mathbf{r_D}\f$.
* 		If the point on the line and the given vertex are
* 		coincident then obviously the vertex is on the line
* 		otherwise the vertex is on the line if the squared
* 		length of the cross product of the line direction and
* 		the direction of the vertex (with respect to the point
* 		on the line) is less than ALG_DBL_TOLLERANCE, ie if
* 		\f[
		\| (\mathbf{p_0} - \mathbf{r_0}) \times \mathbf{r_D} \|
		< \epsilon
 		\f]
* \param	p0			Given vertex.
* \param	r0			Point on line.
* \param	rD			Direction of line.
*/
int		WlzGeomVtxOnLine3D(WlzDVertex3 p0, WlzDVertex3 r0,
				  WlzDVertex3 rD)
{
  int		isn = 0;
  double	lSq;
  const double	tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_3_SUB(p0, p0, r0);
  lSq = WLZ_VTX_3_SQRLEN(p0);
  if(lSq < tol)
  {
    isn = 1;
  }
  else
  {
    WLZ_VTX_3_CROSS(r0, p0, rD);
    lSq = WLZ_VTX_3_SQRLEN(r0);
    if(lSq < tol)
    {
      isn = 1;
    }
  }
  return(isn);
}

/*!
* \return	Interpolated value.
* \ingroup	WlzGeometry
* \brief	Given the coordinates of the vertices of a 2D triangle
* 		and a set of values at each of these vertices, this
* 		function interpolates the value at the given position
* 		which is either inside or on an edge of the triangle.
* 		This is implimented using barycentric coordinates.
* 		Once the barycentric coordinates (\f$\lambda_0\f$,
* 		\f$\lambda_0\f$, \f$\lambda_2\f$) have been computed
* 		then the interpolated value is given by:
* 		\f$v = \sum_{i=0}^{2}{v_i \lambda_i}\f$.
* 		If the determinant is zero in solving for the barycentric
* 		coordinates then the interpolated value is just the
* 		mean of the given values.
* \param	p0			First vertex of triangle.
* \param	p1			Second vertex of triangle.
* \param	p2			Third vertex of triangle.
* \param	v0			Value at first vertex of triangle.
* \param	v1			Value at second vertex of triangle.
* \param	v2			Value at third vertex of triangle.
* \param	pX			Given position, which is within
* 					(or on) the triangle.
*/
extern double	WlzGeomInterpolateTri2D(WlzDVertex2 p0, WlzDVertex2 p1,
                                        WlzDVertex2 p2,
					double v0, double v1, double v2,
					WlzDVertex2 pX)
{
  double	
		l0,
		l1,
		l2,
		del,
		val;
  WlzDVertex2	q0,
  		q1,
		qX;
  const double	eps = 1.0e-10;

  q0.vtX = p0.vtX - p2.vtX;
  q1.vtX = p1.vtX - p2.vtX;
  q0.vtY = p0.vtY - p2.vtY;
  q1.vtY = p1.vtY - p2.vtY;
  del = (q0.vtX * q1.vtY) - (q1.vtX * q0.vtY);
  if(fabs(del) > eps)
  {
    qX.vtX = pX.vtX - p2.vtX;
    qX.vtY = pX.vtY - p2.vtY;
    l0 = ((q1.vtY * qX.vtX) - (q1.vtX * qX.vtY)) / del;
    l1 = ((q0.vtX * qX.vtY) - (q0.vtY * qX.vtX)) / del;
    l2 = 1.0 - (l0 + l1);
    val = (l0 * v0) + (l1 * v1) + (l2 * v2);
  }
  else
  {
    val = (v0 + v1 + v2) / 3.0;
  }
  return(val);
}

/*!
* \return       Interpolated value.
* \ingroup	WlzGeometry
* \brief	Given the vertex coordinates of an irregular possible
* 		non-convex 2D polygon ordered counter-clockwise and a set
* 		of values at each of these vertices, this function
* 		interpolates the value at the given position which must be
* 		inside the convex hull of the polygon, on an edge of the
* 		convex hull of the polygon or coincident with one of the
* 		vertices of it's convex hull.
* 		This function first calls WlzConvHullClarkson2D() to compute
* 		the convex hull and then WlzGeomInterpolateConvexPoly2D()
* 		to perform the interpolation. If the polygonis known to be
* 		convex then WlzGeomInterpolateConvexPoly2D() should be called
* 		directly since temporary workspaces are allocated.
* \param	n			Number of polygon vertices, which is
* 					the same as the number of values and
* 					generalised barycentric coordinates.
* \param	p			The vertices of the convex polygon
* 					(ordered counter-clockwise).
* \param        v			Values at the corresponding polygon
* 					vertices.
* \param	w			Used to compute and return the
* 					generalised barycentric coordinates
* 					of the given position.
* 					The coordinate values of vertices
* 					outside the polygon's convex hull will
* 					be zero.
* \param	q			Given position, which is must be
* 					within or on the convex hull of the
* 					polygon.
* \param	dstErr			Destination error pointer, may be NULL.
*/
extern double	WlzGeomInterpolatePoly2D(int n, WlzDVertex2 *p,
				         double *v, double *w,
					 WlzDVertex2 q, WlzErrorNum *dstErr)
{
  int		m;
  double	rVal = 0.0;
  int		*idx = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Compute the convex hull of thge given polygon. */
  m = WlzConvHullClarkson2D(p, n, &idx, &errNum);
  if(m == n)
  {
    /* Polygon is convex. */
    rVal = WlzGeomInterpolateConvexPoly2D(n, p, v, w, q);
  }
  else if(m > 0)
  {
    /* Polygon is not convex so use it's convex hull for interpolation. */
    double 	*vC = NULL,
		*wC = NULL;
    WlzDVertex2	*pC = NULL;

    if(((pC = (WlzDVertex2 *)AlcMalloc(m * sizeof(WlzDVertex2))) == NULL) ||
       ((vC = (double *)AlcMalloc(m * sizeof(double))) == NULL) ||
       ((wC = (double *)AlcMalloc(m * sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int	i;

      for(i = 0; i < m; ++i)
      {
        pC[i] = p[idx[i]];
        vC[i] = v[idx[i]];
      }
      rVal = WlzGeomInterpolateConvexPoly2D(m, pC, vC, wC, q);
      for(i = 0; i < m; ++i)
      {
        w[idx[i]] = wC[i];
      }
      for(i = 0; i < m; ++i)
      {
        w[idx[i]] = wC[i];
      }
    }
    AlcFree(pC);
    AlcFree(vC);
    AlcFree(wC);
  }
  AlcFree(idx);
  return(rVal);
}

/*!
* \return       Interpolated value.
* \ingroup	WlzGeometry
* \brief	Given the vertex coordinates of an irregular convex 2D
* 		polygon ordered counter-clockwise and a set of values at
* 		each of these vertices, this function interpolates the value
* 		at the given position which must be inside the polygon, on
* 		an edge of the polygon or coincident with one of it's vertices.
* 		This is implimented using general barycentric coordinates,
* 		see the paper: "Generalized Barycentric Coordinates on
* 		Irregular Polygons" Mark Mayer, etal, Journal of Graphics
* 		Tools 2002. All parameters of this function must be valid.
* \param	n			Number of polygon vertices, which is
* 					the same as the number of values and
* 					generalised barycentric coordinates.
* \param	p			The vertices of the convex polygon
* 					(ordered counter-clockwise).
* \param        v			Values at the corresponding polygon
* 					vertices.
* \param	w			Used to compute and return the
* 					generalised barycentric coordinates
* 					of the given position.
* \param	q			Given position, which is must be
* 					in or on the polygon.
*/
extern double	WlzGeomInterpolateConvexPoly2D(int n, WlzDVertex2 *p,
				        double *v, double *w,
					WlzDVertex2 q)
{
  int		i,
		iPrv,
		iNxt,
  		qOnP = -1,
		qOnE = -1;
  double 	a = 0.0,
  		s = 0.0;
  WlzDVertex2	d[2];
  double	l[3];
  const double	tol = 1.0e-10;

  for(i = 0; i < n; ++i)
  {
    iPrv = (i + n - 1) % n;
    iNxt = (i + 1) % n;
    WLZ_VTX_2_SUB(d[0], q, p[i]);
    l[0] = WLZ_VTX_2_SQRLEN(d[0]);
    if(l[0] < tol)
    {
      /* The given position is coincident with the current vertex of the
       * polygon. */
      qOnP = i;
      break;
    }
    /* Check for given position on edge between p[i] and p[iNxt] using
     * squares of lengths to avoid expensive sqrt(). */
    WLZ_VTX_2_SUB(d[1], p[iNxt], p[i]);
    l[1] = ((d[0].vtX * d[1].vtY) - (d[0].vtY - d[1].vtX));
    l[1] *= l[1];
    l[2] = WLZ_VTX_2_SQRLEN(d[1]);
    if(l[1] < tol * l[2])
    {
      /* Position is the edge from the current to the next vertex of the
       * polygon. */
      qOnE = i;
      break;
    }
    w[i] = (WlzGeomCot2D3(q, p[i], p[iPrv]) +
            WlzGeomCot2D3(q, p[i], p[iNxt]) / l[0]);
    s += w[i];
  }
  if(qOnP >= 0)
  {
    /* The given position is coincident with a vertex of the polygon. */
    for(i = 0; i < n; ++i)
    {
      w[i] = 0.0;
    }
    w[qOnP] = 1.0;
    a = v[qOnP];
  }
  else if(qOnE >= 0)
  {
    /* The given position is coincident with the edge p[i], p[iNxt] of the
     * polygon. */
    for(i = 0; i < n; ++i)
    {
      w[i] = 0.0;
    }
    l[0] = WLZ_VTX_2_LENGTH(d[0]);
    l[1] = WLZ_VTX_2_LENGTH(d[1]);
    s = 1.0 / (l[0] + l[1]);
    w[i] = l[1] * s;
    w[iNxt] = l[0] * s;
    a = (v[i] * w[i]) + (v[iNxt] * w[iNxt]);
  }
  else
  {
    s = 1.0 / s;
    for(i = 0; i < n; ++i)
    {
      w[i] *= s;
      a += v[i] * w[i];
    }
  }
  return(a);
}

/*!
* \return	Cotangent of the angle.
* \ingroup	WlzGeometry
* \brief	Computes the cotangent of the angle at B within the
* 		triable A,B,C.
* \param	a			Position of triangle vertex A.
* \param	b			Position of triangle vertex B.
* \param	c			Position of triangle vertex C.
*/
static double	WlzGeomCot2D3(WlzDVertex2 a, WlzDVertex2 b, WlzDVertex2 c)
{
  double	d,
  		e;
  const double	tol = 1.0e-10;

  WLZ_VTX_2_SUB(a, a, b);
  WLZ_VTX_2_SUB(c, c, b);
  d = WLZ_VTX_2_DOT(a, c);
  if(fabs(d) > tol)
  {
    e = fabs((a.vtX * b.vtY) - (a.vtY * b.vtX));
    d = (e > tol)? d / e: DBL_MAX;
  }
  return(d);
}

/*!
* \return	Interpolated value.
* \ingroup	WlzGeometry
* \brief	Given the coordinates of the vertices of a 3D tetrahedron
* 		and a set of values at each of these vertices, this
* 		function interpolates the value at the given position
* 		which is either inside or on a face/edge of the tetrahedron.
* 		This is implimented using barycentric coordinates.
* 		Once the barycentric coordinates (\f$\lambda_0\f$,
* 		\f$\lambda_0\f$, \f$\lambda_2\f$, \f$\lambda_3\f$)
* 		have been computed then the interpolated value is given by:
* 		\f$v = \sum_{i=0}^{3}{v_i \lambda_i}\f$.
* 		If the determinant is zero in solving for the barycentric
* 		coordinates then the interpolated value is just the
* 		mean of the given values.
* \param	p0			First vertex of tetrahedron.
* \param	p1			Second vertex of tetrahedron.
* \param	p2			Third vertex of tetrahedron.
* \param	p3			Fourth vertex of tetrahedron.
* \param	v0			Value at first vertex of tetrahedron.
* \param	v1			Value at second vertex of tetrahedron.
* \param	v2			Value at third vertex of tetrahedron.
* \param	v3			Value at fourth vertex of tetrahedron.
* \param	pX			Given position, which is within
* 					(or on) the tetrahedron.
*/
extern double	WlzGeomInterpolateTet3D(WlzDVertex3 p0, WlzDVertex3 p1,
                                        WlzDVertex3 p2, WlzDVertex3 p3,
					double v0, double v1,
					double v2, double v3,
					WlzDVertex3 pX)
{
  double	d,
  		vX;
  double	c[3],
  		l[4];
  WlzDVertex3	q0,
  		q1,
		q2,
		qX;
  const double	eps = 1.0e-10;

  WLZ_VTX_3_SUB(q0, p0, p3);
  WLZ_VTX_3_SUB(q1, p1, p3);
  WLZ_VTX_3_SUB(q2, p2, p3);
  WLZ_VTX_3_SUB(qX, pX, p3);
  d = q0.vtX * (q1.vtY * q2.vtZ - q2.vtY * q1.vtZ) +
      q1.vtX * (q2.vtY * q0.vtZ - q0.vtY * q2.vtZ) +
      q2.vtX * (q0.vtY * q1.vtZ - q1.vtY * q0.vtZ);
  if(fabs(d) > eps)
  {
    c[0] = qX.vtX * (q1.vtY * q2.vtZ - q2.vtY * q1.vtZ) +
	   q1.vtX * (q2.vtY * qX.vtZ - qX.vtY * q2.vtZ) +
	   q2.vtX * (qX.vtY * q1.vtZ - q1.vtY * qX.vtZ);
    c[1] = q0.vtX * (qX.vtY * q2.vtZ - q2.vtY * qX.vtZ) +
           qX.vtX * (q2.vtY * q0.vtZ - q0.vtY * q2.vtZ) +
           q2.vtX * (q0.vtY * qX.vtZ - qX.vtY * q0.vtZ);
    c[2] = q0.vtX * (q1.vtY * qX.vtZ - qX.vtY * q1.vtZ) +
           q1.vtX * (qX.vtY * q0.vtZ - q0.vtY * qX.vtZ) +
           qX.vtX * (q0.vtY * q1.vtZ - q1.vtY * q0.vtZ);
    l[0] = c[0] / d;
    l[1] = c[1] / d;
    l[2] = c[2] / d;
    l[3] = 1.0 - (l[0] + l[1] + l[2]);
    vX = l[0] * v0  + l[1] * v1 + l[2] * v2 + l[3] * v3;
  }
  else
  {
    vX = (v0 + v1 + v2 + v3) / 4.0;
  }
  return(vX);
}

/*!
* \ingroup	WlzGeometry
* \brief	Given the three vertices of a triangle in 3D computes the 2D
* 		coordinates within the plane of the thriangle.
*
* 		If the 3D coordinates of the vertices of the triangle are
* 		\f$\mathbf{p_0}\f$, \f$\mathbf{p_1}\f$ and \f$\mathbf{p_2}\f$
* 		then this function computes the coordinates of vertices in
* 		the plane of the triangle, such that:
* 		\f$\mathbf{q_0} = (0,0)\f$,
* 		\f$\mathbf{q_1} = (\|\mathbf{l_1}\|,0)\f$ and
* 		\f$\mathbf{q_2} = (\mathbf{l_2}\cdot\mathbf{u}\f$.
* 		Where:
* 		  \f$\mathbf{l_1} = \mathbf{p_1} - \mathbf{p_0}\f$,
* 		  \f$\mathbf{l_2} = \mathbf{p_2} - \mathbf{p_0}\f$,
* 		  \f$\mathbf{u} = \frac{\mathbf{1}}{\|\mathbf{l_1}\|}
 		                  \mathbf{l_1}\f$,
* 		  \f$\mathbf{n} = \frac{\mathbf{1}}
 		                       {\|\mathbf{l_1} \times \mathbf{l_2}\|}
				  \mathbf{l_1} \times \mathbf{l_2}\f$
* 		and
* 		  \f$\mathbf{v} = \mathbf{n} \times \mathbf{l_1}\f$.
*
*		The first vertex in the plane is not returned because it's
*		coordinates are always (0,0).
* \param	p0			First vertex of triangle (origin
* 					of the 2D plane).
* \param	p1			Second vertex of triangle (on the
* 					\f$\mathbf{u}\f$ axis in the 2D
* 					plane).
* \param	p2			Third vertex of triangle.
* \param	dstQ1			Destination pointer for the second
* 					vertex in the plane, must not be
* 					NULL.
* \param	dstQ2			Destination pointer for the third
* 					vertex in the plane, must not be
* 					NULL.
*/
void		WlzGeomMap3DTriangleTo2D(WlzDVertex3 p0,
				WlzDVertex3 p1, WlzDVertex3 p2,
				WlzDVertex2 *dstQ1, WlzDVertex2 *dstQ2)
{
  double	ln1,
  		ln2,
		t;
  WlzDVertex2	q1,
  		q2;
  WlzDVertex3	l1,
  		l2,
		n,
  		u,
  		v;
  const double	tol = ALG_DBL_TOLLERANCE;

  WLZ_VTX_3_SUB(l1, p1, p0);
  WLZ_VTX_3_SUB(l2, p2, p0);
  ln1 = WLZ_VTX_3_LENGTH(l1);
  ln2 = WLZ_VTX_3_LENGTH(l1);
  if(fabs(ln1) < tol)
  {
    q1.vtX = 0.0;
    q1.vtY = 0.0;
    q2.vtX = 0.0;
    q2.vtY = ln2;
  }
  else if(fabs(ln2) < tol)
  {
    q1.vtX = ln1;
    q1.vtY = 0.0;
    q2.vtX = 0.0;
    q2.vtY = 0.0;
  }
  else
  {
    q1.vtX = ln1;
    q1.vtY = 0.0;
    t = 1.0 / ln1;
    WLZ_VTX_3_SCALE(u, l1, t);
    WLZ_VTX_3_CROSS(n, l1, l2);
    t = WLZ_VTX_3_LENGTH(n);
    if(fabs(t) < tol)
    {
      q2.vtX = ln1 + ln2;
      q2.vtY = 0.0;
    }
    else
    {
      t = 1.0 / t;
      WLZ_VTX_3_SCALE(n, n, t);
      WLZ_VTX_3_CROSS(v, n, u);
      q2.vtX = WLZ_VTX_3_DOT(l2, u);
      q2.vtY = WLZ_VTX_3_DOT(l2, v);
    }
  }
  *dstQ1 = q1;
  *dstQ2 = q2;
}

/*!
* \return	The result of the intersection test: 0 - no intersection,
* 		1 - triangle and box touch or intersect.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the given triangle and
* 		the axis aligned bounding box in 3D using the Separating Axis
* 		Theorem (SAT).
*
* 		Given an axis aligned bounding box and a triangle this
* 		function tests for an intersection using the Separating
* 		Axis Theorem (SAT) which states : Two 3D convex domains do
* 		not intersect iff there exists a line, called a separating
* 		axis, on which projection intervals of the domains do not
* 		intersect.  The minimal set of axes that need to be considered
* 		is formed by the normals to all faces of the (polyhedral)
* 		domains and the cross product of all edge combinations in
* 		which one edge is from each polyhedron. For an axis aligned
* 		bounding box and a triangle in 3D, this is equivalent to
* 		testing for the intersection of the given axis aligned
* 		bounding box with the axis aligned bounding box of the
* 		triangle, testing for intersections on the axes
* 		normal to the face of the triangle and testing
* 		for intersection along the cross product of the axis
* 		aligned bounding box - triangle edges. The mathematics
* 		are simplified by the box being axis aligned.
*
* 		The algorithm may return false positives when the domains
* 		are very close to touching.
* \param	t0			First vertex of the triangle.
* \param	t1			Second vertex of the triangle.
* \param	t2			Third vertex of the triangle.
* \param	b0			Minimum of the bounding box.
* \param	b1			Maximum of the bounding box.
* \param	tst			Determines the actual intersection
* 					tests used:
* 					0 - AABB / triangle.
* 					1 - AABB / AABB(triangle) only.
* 					2 - AABB / triangle omitting the
* 					    AABB / AABB(tetrahedron) test
* 					    this is probably only useful if
* 					    the AABB / AABB(triangle) are
* 					    known to intersect.
*/
int		WlzGeomTriangleAABBIntersect3D(WlzDVertex3 t0, WlzDVertex3 t1,
				    WlzDVertex3 t2, WlzDVertex3 b0,
				    WlzDVertex3 b1, int tst)
{
  int		idx,
  		isn = 1;
  double	l;
  WlzDVertex3	b,
  		c,
		d;
  WlzDVertex3	e[3],
  		t[3],
		x[8];
  const double	tol = ALG_DBL_TOLLERANCE;

  /* Make origin centroid of the AABB. */
  c.vtX = (b0.vtX + b1.vtX) * 0.5;
  c.vtY = (b0.vtY + b1.vtY) * 0.5;
  c.vtZ = (b0.vtZ + b1.vtZ) * 0.5;
  WLZ_VTX_3_SUB(b, b1, c);
  WLZ_VTX_3_SUB(t[0], t0, c);
  WLZ_VTX_3_SUB(t[1], t1, c);
  WLZ_VTX_3_SUB(t[2], t2, c);
  /* Check AABB and the AABB of the triangle intersect. This is equivalent
   * to checking for an intersection using projections of the triangle
   * onto vectors perpendicular to the faces of the AABB. */
  if((tst == 0) || (tst == 1))
  {
    /* Compute the AABB of the triangle. */
    WlzDVertex3	bT[2];

    bT[0] = bT[1] = t[0];
    for(idx = 1; idx <= 2; ++idx)
    {
      if(t[idx].vtX < bT[0].vtX)
      {
	bT[0].vtX = t[idx].vtX;
      }
      else if(t[idx].vtX > bT[1].vtX)
      {
	bT[1].vtX = t[idx].vtX;
      }
      if(t[idx].vtY < bT[0].vtY)
      {
	bT[0].vtY = t[idx].vtY;
      }
      else if(t[idx].vtY > bT[1].vtY)
      {
	bT[1].vtY = t[idx].vtY;
      }
      if(t[idx].vtZ < bT[0].vtZ)
      {
	bT[0].vtZ = t[idx].vtZ;
      }
      else if(t[idx].vtZ > bT[1].vtZ)
      {
	bT[1].vtZ = t[idx].vtZ;
      }
    }
    /* Compare AABB of triangle with given AABB. */
    if((-b.vtX - bT[1].vtX > tol) || ( bT[0].vtX - b.vtX > tol) ||
       (-b.vtY - bT[1].vtY > tol) || ( bT[0].vtY - b.vtY > tol) ||
       (-b.vtZ - bT[1].vtZ > tol) || ( bT[0].vtZ - b.vtZ > tol))
    {
      isn = 0;        /* No intersection of the AABB with AABB(tetrahedron). */
    }
  }
  /* Check for intersection using projections of the AABB onto the vector
   * perpendicular to the faces of the tetrahedron. */
  if((tst == 0) || (tst == 2))
  {
    if(isn != 0)
    {
      /* Compute the 3 edge vectors for the triangle. */
      WLZ_VTX_3_SUB(e[0], t[1], t[0]);
      WLZ_VTX_3_SUB(e[1], t[2], t[1]);
      WLZ_VTX_3_SUB(e[2], t[0], t[2]);
      x[0].vtX = -b.vtX; x[0].vtY = -b.vtY; x[0].vtZ = -b.vtZ;
      x[1].vtX = -b.vtX; x[1].vtY = -b.vtY; x[1].vtZ =  b.vtZ;
      x[2].vtX = -b.vtX; x[2].vtY =  b.vtY; x[2].vtZ = -b.vtZ;
      x[3].vtX = -b.vtX; x[3].vtY =  b.vtY; x[3].vtZ =  b.vtZ;
      x[4].vtX =  b.vtX; x[4].vtY = -b.vtY; x[4].vtZ = -b.vtZ;
      x[5].vtX =  b.vtX; x[5].vtY = -b.vtY; x[5].vtZ =  b.vtZ;
      x[6].vtX =  b.vtX; x[6].vtY =  b.vtY; x[6].vtZ = -b.vtZ;
      x[7].vtX =  b.vtX; x[7].vtY =  b.vtY; x[7].vtZ =  b.vtZ;
      /* Check for an intersection between the normal vector and the AABB. */
      WLZ_VTX_3_CROSS(d, e[0], e[1]); 
      l = WLZ_VTX_3_SQRLEN(d);
      if(l > tol)
      {
	isn = WlzGeomTriAABBIsnDir(d, t, x);
      }
    }
    /* Check for intersection along vectors which are the cross product of
     * all possible edge combinations where one edge is from the AABB and
     * the other is from the triangle. */
    if(isn != 0)
    {
      /* Cross product of tetrahedron edges with the AABB edges can be
       * computed avoiding a full cross product as (0, -z, y), (z, 0, -x),
       * and (-y, x, 0) where (x, y, z) is the tetrahedron edge vector.
       * Edge - edge intersections with the edges parallel and consequent
       * zero length cross product are ignored because the AABB - AABB
       * intersection test above will find these intersections. */
      idx = 0;
      do
      {
	d.vtX =  0.0;
	d.vtY = -e[idx].vtZ;
	d.vtZ =  e[idx].vtY;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTriAABBIsnDir(d, t, x);
	  if(isn == 0)
	  {
	    break;
	  }
	}
	d.vtX =  e[idx].vtZ;
	d.vtY =  0.0;
	d.vtZ = -e[idx].vtX;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTriAABBIsnDir(d, t, x);
	  if(isn == 0)
	  {
	    break;
	  }
	}
	d.vtX = -e[idx].vtY;
	d.vtY =  e[idx].vtX;
	d.vtZ =  0.0;
	l = WLZ_VTX_3_SQRLEN(d);
	if(l > tol)
	{
	  isn = WlzGeomTriAABBIsnDir(d, t, x);
	}
      } while((isn >= 0) && (++idx < 3));
    }
  }
  return(isn);
}

/*!
* \return	Intersection code : 0 - no intersection,
* 		1 - triangle and box touch or intersect.
* \ingroup	WlzGeometry
* \brief	Intersection interval test code for
* 		WlzGeomTriangleAABBIntersect3D().
* \param	d			Direction vector.
* \param	t			Array of three traingle vertices.
* \param	b			Array of eight box vertices.
*/
static int	WlzGeomTriAABBIsnDir(WlzDVertex3 d,
                                     WlzDVertex3 t[], WlzDVertex3 b[])
{
  int		idx;
  double	f;
  double	p[2],
  		q[2];
  int		isn = 1;
  const double  tol = ALG_DBL_TOLLERANCE;

  p[0] = p[1] = WLZ_VTX_3_DOT(d, t[0]);
  for(idx = 1; idx < 3; ++idx)
  {
    f = WLZ_VTX_3_DOT(d, t[idx]);
    if(f < p[0])
    {
      p[0] = f;
    }
    else if(f > p[1])
    {
      p[1] = f;
    }
  }
  q[0] = q[1] = WLZ_VTX_3_DOT(d, b[0]);
  for(idx = 1; idx < 8; ++idx)
  {
    f = WLZ_VTX_3_DOT(d, b[idx]);
    if(f < q[0])
    {
      q[0] = f;
    }
    else if(f > q[1])
    {
      q[1] = f;
    }
  }
  if((p[0] - q[1] > tol) || (q[0] - p[1] > tol))
  {
    isn = 0;                   /* No intersection of the AABB with triangle. */
  }
  return(isn);
}

/*!
* \return	Non zero if the two triangles intersect or touch.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the two triangles in 2D
* 		using the by testing for intersections between the line
* 		segments of the triangles and then for either triangle
* 		being contained within the other.
*
* 		See WlzGeomTriangleTriangleIntersect2DA().
* \param	s0			First vertex of the first triangle.
* \param	s1			Second vertex of the first triangle.
* \param	s2			Third vertex of the first triangle.
* \param	t0			First vertex of the second triangle.
* \param	t1			Second vertex of the second triangle.
* \param	t2			Third vertex of the second triangle.
*/
int		WlzGeomTriangleTriangleIntersect2D(WlzDVertex2 s0,
				WlzDVertex2 s1, WlzDVertex2 s2,
				WlzDVertex2 t0, WlzDVertex2 t1,
				WlzDVertex2 t2)
{
  int		isn;
  WlzDVertex2	s[3],
  		t[3];

  s[0] = s0; s[1] = s1; s[2] = s2;
  t[0] = t0; t[1] = t1; t[2] = t2;
  isn = WlzGeomTriangleTriangleIntersect2DA(s, t);
  return(isn);
}

/*!
* \return	Non zero if the two triangles intersect or touch.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the two triangles in 2D
* 		using the by testing for intersections between the line
* 		segments of the triangles and then for either triangle
* 		being contained within the other.
*
* 		The algorithm may return false positives when the domains
* 		are very close to touching.
* \param	s			Array of 3 vertices in 1st triangle.
* \param	t			Array of 3 vertices in 2nd triangle.
*/
int		WlzGeomTriangleTriangleIntersect2DA(WlzDVertex2 s[],
				WlzDVertex2 t[])
{
  int		i0,
  		i1,
		j0,
		j1,
		isn = 0;

  i0 = 2;
  for(i1 = 0; i1 < 3; ++i1)
  {
    j0 = 2;
    for(j1 = 0; j1 < 3; ++j1)
    {
      if(WlzGeomLineSegmentsIntersect(s[i0], s[i1], t[j0], t[j1], NULL) != 0)
      {
        isn = 1;
	goto RETURN;
      }
      j0 = j1;
    }
    i0 = i1;
  }
  isn = WlzGeomVxInTriangle2D(s[0], s[1], s[2], t[0]) >= 0;
  if(isn == 0)
  {
    isn = WlzGeomVxInTriangle2D(t[0], t[1], t[2], s[0]) >= 0;
  }
RETURN:
  return(isn);
}

/*!
* \return	Non zero if the two triangles intersect or touch.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the two triangles in 3D
* 		using the Separating Axis Theorem (SAT).
*
* 		See WlzGeomTriangleTriangleIntersect3DA().
* \param	s0			First vertex of the first triangle.
* \param	s1			Second vertex of the first triangle.
* \param	s2			Third vertex of the first triangle.
* \param	t0			First vertex of the second triangle.
* \param	t1			Second vertex of the second triangle.
* \param	t2			Third vertex of the second triangle.
*/
int		WlzGeomTriangleTriangleIntersect3D(WlzDVertex3 s0,
				WlzDVertex3 s1, WlzDVertex3 s2,
				WlzDVertex3 t0, WlzDVertex3 t1,
				WlzDVertex3 t2)
{
  int		isn;
  WlzDVertex3	s[3],
  		t[3];

  s[0] = s0; s[1] = s1; s[2] = s2;
  t[0] = t0, t[1] = t1; t[2] = t2;
  isn = WlzGeomTriangleTriangleIntersect3DA(s, t);
  return(isn);
}

/*!
* \return	Non zero if the two triangles intersect or touch.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between the two triangles in 3D
* 		using the Separating Axis Theorem (SAT).
*
* 		Given two triangles this function tests for an intersection
* 		using the Separating Axis Theorem (SAT) which states:
* 		Two 3D convex domains do not intersect iff there exists
* 		a line, called a separating axis, on which projection
* 		intervals of the domains do not intersect.  The minimal
* 		set of axes that need to be considered is formed by the
* 		normals to all faces of the (polyhedral) domains and the
* 		cross product of all edge combinations in which one edge
* 		is from each polyhedron.
*
* 		This code is based on "A Fast Triangle-Triangle Intersection
* 		Test" Tomas Moller, Journal of Graphics, GPU, and Game Tools
* 		1997 2(2) pp25-30.
*
* 		The algorithm may return false positives when the domains
* 		are very close to touching.
* \param	s			Array of vertices in 1st triangle.
* \param	t			Array of vertices in 2nd triangle.
*/
int		WlzGeomTriangleTriangleIntersect3DA(WlzDVertex3 s[],
				WlzDVertex3 t[])
{
  int		idx,
		cop = 0,
  		isn = 0;
  double	tmp,
  		ds0ds1,
		ds0ds2,
  		dt0dt1,
		dt0dt2;
  WlzDVertex3	e,
  		ae;
  double	d[2],
  		ds[3],
		dt[3],
		is[2],
		it[2],
		ps[3],
		pt[3];
  WlzDVertex3	n[2];

  isn = WlzGeomTriTriPlaneTest(dt, n + 0, d + 0, &dt0dt1, &dt0dt2, t, s);
  if(isn != 0)
  {
    isn = WlzGeomTriTriPlaneTest(ds, n + 1, d + 1, &ds0ds1, &ds0ds2, s, t);
    if(isn != 0)
    {
      isn = 0;
      /* Compute direction of the intersection line. */
      WLZ_VTX_3_CROSS(e, n[0], n[1]);
      /* Find the largest component of e and use this to do a simplified
       * projection. */
      WLZ_VTX_3_FABS(ae, e);
      idx = 0;
      tmp = ae.vtX;
      if(ae.vtY > tmp)
      {
        tmp = ae.vtY;
	idx = 1;
      }
      if(ae.vtZ > tmp)
      {
	idx = 2;
      }
      switch(idx)
      {
        case 0: /* 0 - x component largest. */
	  ps[0] = s[0].vtX;
	  ps[1] = s[1].vtX;
	  ps[2] = s[2].vtX;
          break;
	case 1: /* 1 - y component largest. */
	  ps[0] = s[0].vtY;
	  ps[1] = s[1].vtY;
	  ps[2] = s[2].vtY;
          break;
	default: /* 2 - z component largest. */
	  ps[0] = s[0].vtZ;
	  ps[1] = s[1].vtZ;
	  ps[2] = s[2].vtZ;
          break;
      }
      /* Compute intersection for the 1st triangle (s). */
      cop = WlzGeomTriTri3DIsn(s, ps, ds, ds0ds1, ds0ds2, is);
      if(cop)
      {
        isn = WlzGeomTriTri3DCoplanar(n[0], s, t);
      }
      else
      {
	switch(idx)
	{
	  case 0: /* 0 - x component largest. */
	    pt[0] = t[0].vtX;
	    pt[1] = t[1].vtX;
	    pt[2] = t[2].vtX;
	    break;
	  case 1: /* 1 - y component largest. */
	    pt[0] = t[0].vtY;
	    pt[1] = t[1].vtY;
	    pt[2] = t[2].vtY;
	    break;
	  default: /* 2 - z component largest. */
	    pt[0] = t[0].vtZ;
	    pt[1] = t[1].vtZ;
	    pt[2] = t[2].vtZ;
	    break;
	}
	/* Compute intersection for the 2nd triangle (t). */
	(void )WlzGeomTriTri3DIsn(t, pt, dt, dt0dt1, dt0dt2, it);
	if(is[0] > is[1])
	{
	  tmp = is[0]; is[0] = is[1]; is[1] = tmp;
	}
	if(it[0] > it[1])
	{
	  tmp = it[0]; it[0] = it[1]; it[1] = tmp;
	}
	if((is[1] > it[0]) && (it[1] > is[0]))
	{
	  /* Now know the triangles intersect and are not coplanar. */
	  isn = 1;
	}
      }
    }
  }
  return(isn);
}

/*!
* \return	Planarity code:
*		<ul>
*		  <li>0 1st triangle does not intersect plane of 2nd.</li>
*		  <li>1 1st triangle intersects plane of 2nd.</li>
*		  <li>2 1st and 2nd triangle are co-planar.</li>
* \ingroup	WlzGeometry
* \brief	Given the the vertices of two triangles; Computes the
* 		plane equation of the 2nd triangle and then tests for
* 		possible intersection of the 1st triangle with the plane
* 		of the 2nd. Where the equation of the plane is of the
* 		form:  n.x + d = 0.
* 		All destination pointers must be valid.
* \param	ds			Valid array for the three distances
* 					from the vertices of the 1st triangle
* 					to the plane of the 2nd.
* \param	n			Destination pointer for the normal
* 					of the plane equation.
* \param	d			Destination pointer for the distance
* 					of the plane equation.
* \param	d0d1			Destination pointer for d[0] * d[1],
* 					with robust treatment near zero.
* \param	d0d2			Destination pointer for d[0] * d[2],
* 					with robust treatment near zero.
* \param	s[]			Vertices of the 1st triangle.
* \param	t[]			Vertices of the 2nd triangle.
*/
static int	WlzGeomTriTriPlaneTest(double ds[], WlzDVertex3 *n,
                                       double *d, double *d0d1, double *d0d2,
				       WlzDVertex3 s[], WlzDVertex3 t[])
{
  int		isn = 0;
  WlzDVertex3	e[2];
  const double  tol = ALG_DBL_TOLLERANCE;

  /* Compute the plane of the 2nd triangle. */
  WLZ_VTX_3_SUB(e[0], t[1], t[0]);
  WLZ_VTX_3_SUB(e[1], t[2], t[1]);
  WLZ_VTX_3_CROSS(n[0], e[0], e[1]);
  *d = -(WLZ_VTX_3_DOT(n[0], t[0]));
  /* Put the 1st traingle (s) into the plane equation of the 2nd
   * traingle (t) to compute signed distances to the  plane. */
  ds[0] = WLZ_VTX_3_DOT(n[0], s[0]) + *d;
  ds[1] = WLZ_VTX_3_DOT(n[0], s[1]) + *d;
  ds[2] = WLZ_VTX_3_DOT(n[0], s[2]) + *d;
  if(fabs(ds[0]) < tol)
  {
    ++isn;
    ds[0] = 0.0;
  }
  if(fabs(ds[1]) < tol)
  {
    ++isn;
    ds[1] = 0.0;
  }
  if(fabs(ds[2]) < tol)
  {
    ++isn;
    ds[2] = 0.0;
  }
  if(isn == 3)
  {
    /* All distances are zero so the triangles are planar. */
    isn = 2;
  }
  else
  {
    *d0d1 = ds[0] * ds[1];
    *d0d2 = ds[0] * ds[2];
    if((*d0d1 > 0.0) && (*d0d2 > 0.0))
    {
      /* Distances all have the same sign therefore on intersection is
       * not possible. */
      isn = 0;
    }
    else
    {
      /* Distances have different signs therefore on intersection is
       * possible. */
      isn = 1;
    }
  }
  return(isn);
}

/*!
* \return	Non-zero if the two triangles are coplanar.
* \ingroup	WlzGeometry
* \brief	Checks for coplanarity and computes intersection intervals
* 		for WlzGeomTriangleTriangleIntersect3D().
* \param	s			Vertices of 1st triangle.
* \param	sp			Maximal component of each vertex of
* 					1st triangle.
* \param	d			Distances from vertices of the 1st
* 					triangle to the plane of the 2nd.
* \param	d0d1			Product d[0] * d[1], with robust
* 					treatment near zero.
* \param	d0d2			Product d[0] * d[2], with robust
* 					treatment near zero.
* \param	is			Destination pointer for the two
* 					intersection interval values.
*/
static int	WlzGeomTriTri3DIsn(WlzDVertex3 s[], double sp[], double d[],
				   double d0d1, double d0d2, double is[])
{
  int		coplanar = 0;
  double	tmp;
  double	dd[3],
  		pp[3];
  WlzDVertex3	ss[3];
  const double  tol = ALG_DBL_TOLLERANCE;

  if(d0d1 > 0.0)
  {
    /* Know that d0d2 <= 0.0, ie d[0], d[1] are on the same side, d[2] on the
     * other or on the plane. */
    ss[0] = s[2]; pp[0] = sp[2]; dd[0] = d[2];
    ss[1] = s[0]; pp[1] = sp[0]; dd[1] = d[0];
    ss[2] = s[1]; pp[2] = sp[1]; dd[2] = d[1];
  }
  else if(d0d2 > 0.0)
  {
    /* Know that d0d1 <= 0.0. */
    ss[0] = s[1]; pp[0] = sp[1]; dd[0] = d[1];
    ss[1] = s[0]; pp[1] = sp[0]; dd[1] = d[0];
    ss[2] = s[2]; pp[2] = sp[2]; dd[2] = d[2];
  }
  else if((fabs(d[0]) > tol) || (d[1] * d[2] > tol))
  {
    /* Know that d0d1 <= 0.0 or that d[0] != 0.0 */
    ss[0] = s[0]; pp[0] = sp[0]; dd[0] = d[0];
    ss[1] = s[1]; pp[1] = sp[1]; dd[1] = d[1];
    ss[2] = s[2]; pp[2] = sp[2]; dd[2] = d[2];
  }
  else if(fabs(d[1]) > tol)
  {
    ss[0] = s[1]; pp[0] = sp[1]; dd[0] = d[1];
    ss[1] = s[0]; pp[1] = sp[0]; dd[1] = d[0];
    ss[2] = s[2]; pp[2] = sp[2]; dd[2] = d[2];
  }
  else if(fabs(d[2]) > 0)
  {
    ss[0] = s[2]; pp[0] = sp[2]; dd[0] = d[2];
    ss[1] = s[0]; pp[1] = sp[0]; dd[1] = d[0];
    ss[2] = s[1]; pp[2] = sp[1]; dd[2] = d[1];
  }
  else
  {
    /* Triangles are coplanar */
    coplanar = 1;
  }
  if(coplanar == 0)
  {
    tmp = dd[0] / (dd[0] - dd[1]);
    is[0] = pp[0] + (pp[1] - pp[0]) * tmp;
    tmp = dd[0] /(dd[0] - dd[2]);
    is[1] = pp[0] + (pp[2] - pp[0]) * tmp;
  }
  return(coplanar);
}

/*!
* \return	Non zero if the triangles intersect or touch.
* \ingroup	WlzGeometry
* \brief	Tests for an intersection between two 3D triangles that
* 		are known to be coplanar. This is done by projecting the
* 		triangles onto the common plane and then testing for an
* 		intersection in that plane.
* 		No tests are made for a zero length normal or coincident
* 		vertices within a triangle.
* \param	n			Unit normal to plane of the triangles,
* 					must not be 0.0.
* \param	s			Vertices of the 1st triangle which
* 					must not be coincident.
* \param	t			Vertices of the 2nd triangle which
* 					must not be coincident.
*/
static int 	WlzGeomTriTri3DCoplanar(WlzDVertex3 n,
				        WlzDVertex3 s[], WlzDVertex3 t[])
{
  int		i,
  		isn = 0;
  WlzDVertex3	m;
  WlzDVertex2	u[3],
  		v[3];

  /* Project onto an axis aligned plane that maximises the area of the
   * triangles which avoids vector arithmetic including square roots.
   * Compute projected triangles. */
  WLZ_VTX_3_FABS(m, n);
  if(m.vtX > m.vtY)
  {
    if(m.vtX > m.vtZ)
    {
      for(i = 0; i < 3; ++i)
      {
	u[i].vtX = s[i].vtY; u[i].vtY = s[i].vtZ;
	v[i].vtX = t[i].vtY; v[i].vtY = t[i].vtZ;
      }
    }
    else
    {
      for(i = 0; i < 3; ++i)
      {
	u[i].vtX = s[i].vtX; u[i].vtY = s[i].vtY;
	v[i].vtX = t[i].vtX; v[i].vtY = t[i].vtY;
      }
    }
  }
  else   //* m.vtX<=m.vtY *//
  {
    if(m.vtZ > m.vtY)
    {
      for(i = 0; i < 3; ++i)
      {
	u[i].vtX = s[i].vtX; u[i].vtY = s[i].vtY;
	v[i].vtX = t[i].vtX; v[i].vtY = t[i].vtY;
      }
    }
    else
    {
      for(i = 0; i < 3; ++i)
      {
	u[i].vtX = s[i].vtX; u[i].vtY = s[i].vtZ;
	v[i].vtX = t[i].vtX; v[i].vtY = t[i].vtZ;
      }
    }
  }
  isn = WlzGeomTriangleTriangleIntersect2DA(u, v);
  return(isn);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeometry
* \brief	Computes the principle curvatures of a parabolic surface
* 		fitted to the given vertices at the first of these vertices.
* \param	nC			Number of curvatures required:
* 					  1 computes the Gaussian curvature,
* 					  2 computes both principle
* 					    curvatures.
* \param	dstC			Destination pointer for the curvature
* 					value(s), must not be NULL. If nC == 2
* 					then the values atr Gaussian followed
* 					by mean curvature.
* \param	nrm			Normal at the first vertex.
* \param	nV			Number of vertices, must be >= 3.
* \param	vtx			Array of vertex positions, the first of
* 					which must be the vertex at which the
* 					curvature is to be computed. The
* 					array contents are modified by this
* 					function. Must not be NULL.
*/
WlzErrorNum	WlzGeomCurvature(int nC, double *dstC, WlzDVertex3 nrm,
                                 int nV, WlzDVertex3 *vtx)
{
  int		idV;
  double	len;
  double	*bV= NULL;
  double	**aA;
  AlgMatrix	aM;
  WlzDVertex3   t;
  WlzDVertex3   b[3],
  		r[3];
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  aM.core = NULL;
  if((nV < 3) || (nC < 1) || (nC > 2) ||
     ((len = WLZ_VTX_3_LENGTH(nrm)) < ALG_DBL_TOLLERANCE))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(((bV = (double *)AlcMalloc(sizeof(double) * (nV - 1))) == NULL) ||
	  ((aM.rect = AlgMatrixRectNew((nV - 1), 3, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    aA = aM.rect->array;
    /* Make b[2] = nrm  / |nrm|. */
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(b[2], nrm, len);
    /* Shift st vtx[0] is at the origin. */
    for(idV = 1; idV < nV; ++idV)
    {
      WLZ_VTX_3_SUB(vtx[idV], vtx[idV], vtx[0]);
    }
    /* Find a unit vector b[0] that is perpendicular to b[2]. */
    WLZ_VTX_3_SET(b[1], fabs(nrm.vtX), fabs(nrm.vtY), fabs(nrm.vtZ));
    if(b[1].vtX < b[1].vtY)
    {
      if(b[1].vtX < b[1].vtZ)
      {
        WLZ_VTX_3_SET(b[1], 1, 0, 0);
      }
      else
      {
        WLZ_VTX_3_SET(b[1], 0, 0, 1);
      }
    }
    else
    {
      if(b[1].vtY < b[1].vtZ)
      {
        WLZ_VTX_3_SET(b[1], 0, 1, 0);
      }
      else
      {
        WLZ_VTX_3_SET(b[1], 0, 0, 1);
      }
    }
    WLZ_VTX_3_CROSS(b[0], b[1], b[2]);
    len = WLZ_VTX_3_LENGTH(b[0]);
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(b[0], b[0], len);
    /* Complete the three orthogonal vectors with b[1] = b[2] x b[0]. */
    WLZ_VTX_3_CROSS(b[1], b[2], b[0]);
    len = WLZ_VTX_3_LENGTH(b[1]);
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(b[1], b[1], len);
    /* Compute rotation matrix R which transforms b[0], b[1], b[2] to the
     * x, y, z axes, ie R = B^{-1} where B = (b_0, b_1, b_2). Inverting a
     * 3 x 3 matrix can be done using the cross and tripple product. This
     * gives
     *   R = \frac{1}{b_0 . b_1 x b_2} (b_1 x b_2, b_2 * b_0, b_0 x b_1)
     * where b_i x b_j are now row vectors.
     * Because we've choosen b_0 = b_1 * b_2 and b_0, b_1, b_2 are unit
     * vectors then b_0 . b_1 x b_2 = 1.
     * This gives R = (b_1 x b_2, b_2 * b_0, b_0 x b_1). */
    WLZ_VTX_3_CROSS(r[0], b[1], b[2]);
    WLZ_VTX_3_CROSS(r[1], b[2], b[0]);
    WLZ_VTX_3_CROSS(r[2], b[0], b[1]);
    /* Rotate the shifted vertices about the first using the rotation
     * matrix r and fill in the matrices ready to compute the least squares
     * estimate of the parameters a,b,c in z = ax^2 + bxy + cy^2. */
    for(idV = 1; idV < nV; ++idV)
    {
      t.vtX = WLZ_VTX_3_DOT(r[0], vtx[idV]);
      t.vtY = WLZ_VTX_3_DOT(r[1], vtx[idV]);
      t.vtZ = WLZ_VTX_3_DOT(r[2], vtx[idV]);
      aA[idV - 1][0] = t.vtX * t.vtX;
      aA[idV - 1][1] = t.vtX * t.vtY;
      aA[idV - 1][2] = t.vtY * t.vtY;
      bV[idV - 1] = t.vtZ;
    }
    /* Solve for a, b, c ie matrix equation for x. */
    errNum = WlzErrorFromAlg(AlgMatrixSVSolve(aM, bV, 1.0e-06, NULL));
  }
  /* Compute the curvature values, Gaussian = 4ac - b^2, mean = a + c. */
  dstC[0] = 4.0 * bV[0] * bV[2] - bV[1] * bV[1];
  if(nC == 2)
  {
    dstC[1] = bV[0] + bV[2];
  }
  /* Free allocated storage. */
  AlcFree(bV);
  AlgMatrixFree(aM);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzGeometry
* \brief	Computes the plane which is the least squares best fit to
* 		the given vertices, where this is with respect to the
* 		orthogonal distance from the vertices to the plane.
*
* 		The orthogonal distance regression plane is an eigenvector
* 		problem. This solution is based on one from
* 		http://mathforum.org which is credited to Doctor George.
* 		Starting with the distance from a point to a plane we
* 		wish to find \f$(a,b,c,d)\f$ such as to minimise
* 		\f[
                f(a,b,c,d) = \sum{
		             \frac{a x_i + b y_i + c z_i}
			          {a^2 + b^2 + c^2}}
*                \f]
		setting \f$\frac{\partial f}{\partial d} = 0\f$ gives
*		\f[
		d = -(a x_0 + b y_0 + c z_0)
		\f]
*		where \f$(x_0, y_0, z_0)\f$ is the centroid of the vertices.
*		Using the centroid
*		\f[
		f(a,b,c,d) = \sum{
		             \frac{|a (x_i - x_0) +
			            b (y_i - y_0) +
				    c (z_i - z_0) |}
				  {a^2 + b^2 + c^2}}
*		\f]
*		Define \f$v\f$ \f$M\f$ such that:
*		\f[
		v^T = [a \, b \, c]
		\f]
*		\f[
		M = 
		\left[
		\begin{array}{ccc}
		x_1 - x_0 & y_1 - y_0 & z_1 - z_0 \\
		x_2 - x_0 & y_2 - y_0 & z_2 - z_0 \\
		\cdots    & \cdots    & \cdots    \\
		x_n - x_0 & y_n - y_0 & z_n - z_0
		\end{array}
		\right]
		\f]
*		\f[
		f(v) = \frac{(v^T M^T)(M v)}{v^T v}
		\f]
*		\f[
		f(v) = \frac{v^T (M^T M) v}{v^T v}
		\f]
*		Define \f$A\f$
*		\f[
*		A = M^T M
		\f]
*		The Rayleigh Quotient \f$f(v)\f$ is minimised by the
*		eigenvector of \f$A\f$. However there is no need to compute
*		the eigenvectors of \f$A\f$. The SVD of \f$M\f$ is
*		\f[
*		M = U S V^T
		\f]
*		where \f$S\f$ is a diagonal vector containing the singular
*		values of \f$M\f$. The columns of \f$V\f$ are it's singular
*		vectors and \f$U\f$ is an orthogonal matrix.
*		\f[
*		A = M M^T
		\f]
*		\f[
*		A = (U S V^T)^T (U S V^T)
		\f]
*		\f[
*		A = (V S^T U^T) (U S V^T)
		\f]
*		\f[
*		A = V S^2 V^T
		\f]
*		The decomposition of \f$A\f$ diagonalises the matrix and gives
*		an eigenvector decomposition. It means that the eigenvectors of
*		\f$A\f$ are the squares of the singular values of \f$M\f$ and
*		the eigenvectors of \f$A\f$ are the singular vectors of
*		\f$M\f$. \f$M\f$ is the covarience matrix.
* \param	dstNrm			Destination pointer for the
* 					plane normal.
* \param	dstCen			Destination pointer for the centroid
* 					of the vertices which is on the plane.
* \param	nVtx			umber of vertices.
* \param	vtx			Vector of vertices.
*/
WlzErrorNum	WlzGeometryLSqOPlane(WlzDVertex3 *dstNrm, WlzDVertex3 *dstCen,
		     	int nVtx, WlzDVertex3 *vtx)
{
  double	*sV = NULL;
  WlzDVertex3	cen,
  		nrm;
  AlgMatrix	mM,
  		vM;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  mM.core = NULL;
  vM.core = NULL;
  if(nVtx < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if(((sV = AlcMalloc(sizeof(double) * 3)) == NULL) ||
       ((mM.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL))
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;

    AlgMatrixZero(mM);
    cen = WlzCentreOfMassVtx3D(nVtx, vtx);
    for(idN = 0; idN < nVtx; ++idN)
    {
      double	t;
      double	**m;
      WlzDVertex3 p;

      m = mM.rect->array;
      WLZ_VTX_3_SUB(p, vtx[idN], cen);
      m[0][0] += p.vtX * p.vtX;
      m[1][1] += p.vtY * p.vtY;
      m[2][2] += p.vtZ * p.vtZ;
      t = p.vtX * p.vtY; m[0][1] += t; m[1][0] += t;
      t = p.vtX * p.vtZ; m[0][2] += t; m[2][0] += t;
      t = p.vtY * p.vtZ; m[1][2] += t; m[2][1] += t;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(mM, sV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
    		idM;
    double	vMin;
    double	**v;

    /* Find the smallest eigenvalue (ie singular value). */
    idM = 0;
    vMin = sV[idM] * sV[idM];
    for(idN = 1; idN < 3; ++idN)
    {
      double	val;

      if((val = sV[idN] * sV[idN]) < vMin)
      {
        idM = idN;
	vMin = val;
      }
    }
    v = vM.rect->array;        /* It's the transpose so index appropriately. */
    nrm.vtX = v[0][idM];
    nrm.vtY = v[1][idM];
    nrm.vtZ = v[2][idM];
    if(dstNrm != NULL)
    {
      *dstNrm = nrm;
    }
    if(dstCen != NULL)
    {
      *dstCen = cen;
    }
  }
  AlcFree(sV);
  AlgMatrixFree(mM);
  AlgMatrixFree(vM);
  return(errNum);
}

/*!
* \return	Square of the minimum distance from the test vertex to the
* 		triangle.
* \ingroup	WlzGeometry
* \brief	Computes the minimum distance from the test vertex to the
* 		triangle. This algorithm is based on "Distance Between Point
* 		and Triangle in 3D", David Eberly, Geometric Tools, 1999.
* \param	dstPT			Destination pointer for the position
* 					of vertex in the triangle that is
* 					closest to the test vertex .
* \param	dstZT			Destination pointer, the value of which
* 					will be set to a non-zero value if the
* 					triangle has zero area. May be NULL.
* \param	dstIT			Destination pointer, the value of which
* 					will be set to a non-zero value if the
* 					projected test vertex is within the
* 					triangle. May be NULL.
* \param	dstL0			Destination pointer for the first
* 					barycentric coordinates, may be NULL.
* \param	dstL1			Destination pointer for the second
* 					barycentric coordinates, may be NULL.
* \param	dstL2			Destination pointer for the third
* 					barycentric coordinates, may be NULL.
* \param	vT			The position of the test vertex.
* \param	v0			First vertex of the triangle.
* \param	v1			Second vertex of the triangle.
* \param	v2			Third vertex of the triangle.
*/
double	 	WlzGeomTriangleVtxDistSq3D(WlzDVertex3 *dstPT,
					   int *dstZT, int *dstIT,
					   double *dstL0, double *dstL1,
					   double *dstL2,
					   WlzDVertex3 vT, WlzDVertex3 v0,
					   WlzDVertex3 v1, WlzDVertex3 v2)
{
  int		iT = 0,
  		zT = 1;
  double	a,
  		b,
		c,
		d,
		e,
		f,
		det,
		dist = 0.0;
  double	l[2];
  WlzDVertex3	u1,
  		u2,
		uT,
		pT;
  const double	eps = 1.0e-10;

  WLZ_VTX_3_SUB(u1, v1, v0);
  WLZ_VTX_3_SUB(u2, v2, v0);
  WLZ_VTX_3_SUB(uT, v0, vT);
  a = WLZ_VTX_3_DOT(u1, u1);
  b = WLZ_VTX_3_DOT(u1, u2);
  c = WLZ_VTX_3_DOT(u2, u2);
  d = WLZ_VTX_3_DOT(u1, uT);
  e = WLZ_VTX_3_DOT(u2, uT);
  det = (a * c) - (b * b);
  if(fabs(det) > eps)
  {
    zT = 0;
    l[0] = (b * e) - (c * d);
    l[1] = (b * d) - (a * e);
    if((l[0] + l[1]) < det)
    {
      if(l[0] < 0.0)
      {
	if(l[1] < 0.0)
	{
	  if(d < 0.0)
	  {
	    f = -d / a;
	    l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	    l[1] = 0.0;
	  }
	  else
	  {
	    f = -e / c;
	    l[0] = 0.0;
	    l[1] = WLZ_CLAMP(f, 0.0, 1.0);
	  }
	}
	else
	{
	  f = -e / c;
	  l[0] = 0.0;
	  l[1] = WLZ_CLAMP(f, 0.0, 1.0);
	}
      }
      else if(l[1] < 0.0)
      {
	f = -d / a;
	l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	l[1] = 0.0;
      }
      else
      {
	iT = 1;
	f = 1.0 / det;
	l[0] *= f;
	l[1] *= f;
      }
    }
    else
    {
      if(l[0] < 0.0)
      {
	double	g,
		  h;

	g = b + d;
	h = c + e;
	if(h > g)
	{
	  f = (h - g) / (a - (2 * b) + c);
	  l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	  l[1] = 1.0 - l[0];
	}
	else
	{
	  f = -e / c;
	  l[1] = WLZ_CLAMP(f, 0.0, 1.0);
	  l[0] = 0.0;
	}
      }
      else if(l[1] < 0.0)
      {
	if((a + d) > (b + e))
	{
	  f = (c + e - b - d) / (a - (2 * b) + c);
	  l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	  l[1] = 1.0 - l[0];
	}
	else
	{
	  f = -e / c;
	  l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	  l[1] = 0.0;
	}
      }
      else
      {
	f = (c + e - b - d) /(a - (2 * b) + c);
	l[0] = WLZ_CLAMP(f, 0.0, 1.0);
	l[1] = 1.0 - l[0];
      }
    }
    pT.vtX = v0.vtX + (l[0] * u1.vtX) + (l[1] * u2.vtX);
    pT.vtY = v0.vtY + (l[0] * u1.vtY) + (l[1] * u2.vtY);
    pT.vtZ = v0.vtZ + (l[0] * u1.vtZ) + (l[1] * u2.vtZ);
    if(dstIT)
    {
      *dstIT = iT;
    }
    if(dstPT)
    {
      *dstPT = pT;
    }
    if(dstL0)
    {
      *dstL0 = 1.0 - (l[0] + l[1]);
    }
    if(dstL1)
    {
      *dstL1 = l[0];
    }
    if(dstL2)
    {
      *dstL2 = l[1];
    }
    WLZ_VTX_3_SUB(pT, pT, vT);
    dist = WLZ_VTX_3_SQRLEN(pT);
  }
  if(dstZT)
  {
    *dstZT = zT;
  }
  return(dist);
}
