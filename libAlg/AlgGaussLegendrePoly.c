#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgGaussLegendrePoly_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         AlgGaussLegendrePoly.c
* \author       Bill Hill
* \date         July 2020
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2020],
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
* \brief	Functions for Gauss-Legendre polynomials.
* \ingroup	AlgPoly
*/

#include <Alg.h>

/*!
* \return	Gauss-Legendre polynomial weight.
* \ingroup	AlgPoly
* \brief	Returns the Gauss-Legendre polynomial weight for the
* 		given order and point number. The weights are ordered
* 		such that
* 		AlgGaussLegendrePoints[n,i] < AlgGaussLegendrePoints[n,i+1].
* \param	n		Order, must be in the valid range
* 				[1 - ALG_GAUSSLEGENDRE_ORDER_MAX].
* \param	i		Point number, must be in the valid range
* 				[0 - (n - 1)].
*/
double				AlgGaussLegendreWeights(
				  int n,
				  int i)
{
  double	w;
  const double	weights_0[] = {0.0000000000000000},
                weights_1[] = {2.0000000000000000},
  		weights_2[] = {1.0000000000000000,
		               1.0000000000000000},
		weights_3[] = {0.5555555555555556,
		               0.8888888888888888,
			       0.5555555555555556},
  		weights_4[] = {0.3478548451374538,
			       0.6521451548625461,
			       0.6521451548625461,
			       0.3478548451374538},
                weights_5[] = {0.2369268850561891,
		               0.4786286704993665,
			       0.5688888888888889,
			       0.4786286704993665,
			       0.2369268850561891},
                weights_6[] = {0.1713244923791704,
		               0.3607615730481386,
			       0.4679139345726910,
			       0.4679139345726910,
			       0.3607615730481386,
			       0.1713244923791704};
  const double	*weights[] = {weights_0,
			      weights_1,
			      weights_2,
			      weights_3,
			      weights_4,
			      weights_5,
			      weights_6};


  w = weights[n][i];
  return(w);
}

/*!
* \return	Gauss-Legendre polynomial point.
* \ingroup	AlgPoly
* \brief	Returns the Gauss-Legendre polynomial point for the
* 		given order and point number. The points are ordered
* 		such that
* 		AlgGaussLegendrePoints[n,i] < AlgGaussLegendrePoints[n,i+1].
* \param	n		Order, must be in the valid range
* 				[1 - ALG_GAUSSLEGENDRE_ORDER_MAX].
* \param	i		Point number, must be in the valid range
* 				[0 - (n - 1)].
*/
double				AlgGaussLegendrePoints(
				  int n,
				  int i)
{
  double	p;
  const double	points_0[] = { 0.000000000000000},
                points_1[] = { 0.000000000000000},
  		points_2[] = {-0.5773502691896257,
		               0.5773502691896257},
		points_3[] = {-0.7745966692414834,
		               0.0000000000000000,
			       0.7745966692414834},
  		points_4[] = {-0.8611363115940526,
			      -0.3399810435848563,
			       0.3399810435848563,
			       0.8611363115940526},
                points_5[] = {-0.9061798459386640,
		              -0.5384693101056831,
			       0.0000000000000000,
			       0.5384693101056831,
			       0.9061798459386640},
                points_6[] = {-0.9324695142031521,
		              -0.6612093864662645,
			      -0.2386191860831969,
			       0.2386191860831969,
			       0.6612093864662645,
			       0.9324695142031521};
  const double	*points[] = {points_0,
			     points_1,
			     points_2,
			     points_3,
			     points_4,
			     points_5,
			     points_6};


  p = points[n][i];
  return(p);
}
 
