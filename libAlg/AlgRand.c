#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgRand.c
* \author       richard Baldock, Bill Hill
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
* \brief        Provides functions which produce pseudo-random values.
* \ingroup	AlgRand
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <math.h>
#include <Alg.h>

/*!
* \return	void
* \ingroup	AlgRand
* \brief	Seeds the pseudo-random number generators.
* \param	 seed			Given seed value.
*/
void		AlgRandSeed(long seed)
{
  srand((unsigned int )seed);
}

/*!
* \return	Pseudo-random value.
* \ingroup	AlgRand
* \brief	Produces a pseudo-random value from a uniform
*		distribution over the interval [0.0, 1.0].
*/
double		AlgRandUniform(void)
{
  double	value;

  value = ((double) rand()) / RAND_MAX;
  return(value);
}

/*!
* \return	Pseudo-random value.
* \ingroup	AlgRand
* \brief	Produces a pseudo-random value from a normal
*		distribution over the interval [-1.0, 1.0].
* \param	mu			Mean of distribution.
* \param	sigma			Standard deviation of
*					distribution.
*/
double		AlgRandNormal(double mu, double sigma)
{
  double	value;

  value = AlgRandUniform() + AlgRandUniform() + AlgRandUniform() +
  	  AlgRandUniform() + AlgRandUniform() + AlgRandUniform() +
	  AlgRandUniform() + AlgRandUniform() + AlgRandUniform() +
	  AlgRandUniform() + AlgRandUniform() + AlgRandUniform();
  value = ((value - 6.0) * sigma) + mu;
  return(value);
}
