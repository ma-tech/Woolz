#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgRand.c
* Date:         March 1999
* Author:       Richard Baldock, Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Provides functions which produce pseudo-random values
*		for the	MRC Human Genetics Unit numerical algorithm
*		library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
************************************************************************/
#include <Alg.h>

/************************************************************************
* Function:	AlgRandSeed						*
* Returns:	void							*
* Purpose:	Seeds the pseudo-random number generators.		*
* Global refs:	-							*
* Parameters:	long seed:		Given seed value.		*
************************************************************************/
void		AlgRandSeed(long seed)
{
  srand48(seed);
}

/************************************************************************
* Function:	AlgRandUniform						*
* Returns:	double:			Pseudo-random value.		*
* Purpose:	Produces a pseudo-random value from a uniform		*
*		distribution over the interval [0.0, 1.0].		*
* Global refs:	-							*
* Parameters:	void							*
************************************************************************/
double		AlgRandUniform(void)
{
  double	value;

  value = drand48();
  return(value);
}

/************************************************************************
* Function:	AlgRandNormal						*
* Returns:	double:			Pseudo-random value.		*
* Purpose:	Produces a pseudo-random value from a normal		*
*		distribution over the interval [-1.0, 1.0].		*
* Global refs:	-							*
* Parameters:	double mu:		Mean of distribution.		*
*		double sigma:		Standard deviation of 		*
*					distribution.			*
************************************************************************/
double		AlgRandNormal(double mu, double sigma)
{
  double	value;

  value = drand48() + drand48() + drand48() + drand48() +
          drand48() + drand48() + drand48() + drand48() +
          drand48() + drand48() + drand48() + drand48();
  value = ((value - 6.0) * sigma) + mu;
  return(value);
}
