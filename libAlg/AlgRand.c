#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgRand.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides functions which produce pseudo-random values.
* \ingroup	AlgRand
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <math.h>
#include <Alg.h>

/*!
* \return	<void>
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
