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
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgRand
* @{
*/

#include <Alg.h>

#if defined (CYGWIN) || defined (DARWIN)
#define drand48() (((double) rand()) / RAND_MAX)
#define srand48(X) (srand((unsigned int) X))
#define lrand48() ((long) ((((double) rand()) / RAND_MAX) * (1<<31)))
#endif /* CYGWIN || DARWIN */

/*!
* \return	<void>
* \brief	Seeds the pseudo-random number generators.
* \param	 seed			Given seed value.
*/
void		AlgRandSeed(long seed)
{
  srand48(seed);
}

/*!
* \return				Pseudo-random value.
* \brief	Produces a pseudo-random value from a uniform
*		distribution over the interval [0.0, 1.0].
* \param	<void>
*/
double		AlgRandUniform(void)
{
  double	value;

  value = drand48();
  return(value);
}

/*!
* \return				Pseudo-random value.
* \brief	Produces a pseudo-random value from a normal
*		distribution over the interval [-1.0, 1.0].
* \param	mu			Mean of distribution.
* \param	sigma			Standard deviation of
*					distribution.
*/
double		AlgRandNormal(double mu, double sigma)
{
  double	value;

  value = drand48() + drand48() + drand48() + drand48() +
          drand48() + drand48() + drand48() + drand48() +
          drand48() + drand48() + drand48() + drand48();
  value = ((value - 6.0) * sigma) + mu;
  return(value);
}

/*!
* @}
*/
