#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgMixture_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgMixture.c
* \author       Bill Hill
* \date         January 2000
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
* \brief        Provides a function for computing the maximum liklihood
*               parameters of a mixture of distributions which fit the
*               given data.
* \ingroup      AlgMixture
*/

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

#if defined (CYGWIN) || defined (DARWIN) || defined (_WIN32)
#define drand48() (((double) rand()) / RAND_MAX)
#define srand48(X) (srand((unsigned int) X))
#define lrand48() ((long) ((((double) rand()) / RAND_MAX) * (1<<31)))
#endif /* CYGWIN || DARWIN */

#ifdef _WIN32
#define isnan _isnan
#endif
/*!
* \return	Error code.
* \ingroup 	AlgMixture
* \brief	Computes the maximum liklihood estimate of the
*		parameters of a mixture of normal distributions which
*		best fit the give frequencies. These parameters are the
*		miximg proportions, means, and standard deviations.
*		The log-liklihood and the number of iterations are also
*		computed.
*		This algorithm is based on the article:
*		Agha M. and Ibrahim M.T.
*		Maximum Liklihood Estimation of Mixtures of Distributions
*	 	AS203 Applied Statistics 33(3):327-332, 1984.
* \param	nDbn			Number of distributions in the
*					mixture.
* \param	nCls			Number of classes in the
*					frequency distribution.
* \param	samOrg			Position of first sample (origin).
* \param	samItv			Sample interval.
* \param	freq			Given frequency samples.
* \param	alpha			Estimates of nDbn alpha values,
*					both given and returned.
*				 	Alpha values must be in the
*					range (0,1].
* \param	mu			Estimates of nDbn mean values,
*					both given and returned.
* \param	sd			Estimates of nDbn standard
*					deviation values,
*					both given and returned.
* \param	tol			Difference between two
*					consecutive log liklihood values
*					required to terminate iteration.
* \param	sumFreq			Number of observations.
* \param	dstLL			Destination pointer for the
*					log liklihood.
* \param	dstNItn			Destination pointer for the
*					number of iterations taken.
*/
AlgError	AlgMixtureMLG(int nDbn, int nCls,
			      double samOrg, double samItv, double *freq,
			      double *alpha, double *mu, double *sd,
			      double tol, double sumFreq,
			      double *dstLL, int *dstNItn)
{
  int		idC,
  		idD,
		idD1,
		converged = 0;
  double	tD0,
  		tD1,
		tD2,
		oldLL,
		sumAlpha;
  double	*g = NULL,
  		*vt = NULL,
		*dt = NULL,
		*nt = NULL,
		*newalpha = NULL,
		*newmean = NULL,
		*newsd = NULL;
  double	**f = NULL;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters */
  if((nDbn < 2) || (nDbn > (nCls / 2)) || (samItv < DBL_EPSILON))
  {
    algErr = ALG_ERR_FUNC;
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(nCls < (nDbn * 2))
    {
      algErr = ALG_ERR_FUNC;
    }
  }
  if(algErr == ALG_ERR_NONE)
  {
    for(idC = 0; idC < nCls; ++idC)
    {
      if(freq[idC] < 0.0)
      {
	algErr = ALG_ERR_FUNC;
      }
    }
  }
  if(algErr == ALG_ERR_NONE)
  {
    /* Allocate temporary storage. */
    if((AlcDouble2Malloc(&f, nCls, nDbn) != ALC_ER_NONE) ||
       ((newalpha = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((newmean = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((newsd = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((vt = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((dt = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((nt = (double *)AlcMalloc(sizeof(double) * nDbn)) == NULL) ||
       ((g = (double *)AlcMalloc(sizeof(double) * nCls)) == NULL))
    {
      algErr = ALG_ERR_MALLOC;
    }
  }
  if(algErr == ALG_ERR_NONE)
  {
    /* The iteration loop. */
    oldLL = 0.0;
    *dstNItn = 0;
    do
    {
      ++*dstNItn;
      /* Check alpha, mu and standard deviation. */
      for(idD = 0; idD < nDbn; ++idD)
      {
	if(isnan(alpha[idD]) ||
	   (alpha[idD]) > 1.0 || (alpha[idD] < 0.0))
	{
	  algErr = ALG_ERR_CONVERGENCE;
	}
	else if(isnan(mu[idD]) ||
		(mu[idD] >= samOrg + (samItv * (nCls - 1))) ||
		(mu[idD] <= samOrg))
	{
	  algErr = ALG_ERR_CONVERGENCE;
	}
      }
      if(algErr == ALG_ERR_NONE)
      {
	for(idD = 0; idD < nDbn; ++idD)
	{
	  if(isnan(sd[idD]) || (sd[idD] <= 0.0))
	  {
	    algErr = ALG_ERR_CONVERGENCE;
	  }
	}
      }
      if(algErr == ALG_ERR_NONE)
      {
	idD = 0;
	while((idD < nDbn - 1) && (algErr == ALG_ERR_NONE))
	{
	  idD1 = idD + 1;
	  while((idD1 < nDbn) && (algErr == ALG_ERR_NONE))
	  {
	    if(mu[idD] == mu[idD1])
	    {
	      if(sd[idD] == sd[idD1])
	      {
		algErr = ALG_ERR_FUNC;
	      }
	    }
	    ++idD1;
	  }
	  ++idD;
	}
      }
      if(algErr == ALG_ERR_NONE)
      {

	*dstLL = 0.0;
	for(idC = 0; idC < nCls; ++idC)
	{
	  g[idC] = 0.0;
	  for(idD = 0; idD < nDbn; ++idD)
	  {
	    tD0 = (samOrg + (samItv * idC) - mu[idD]) / sd[idD];
	    tD1 = exp(-0.5 * (tD0 * tD0)) / sd[idD];
	    *(*(f + idC) + idD) = tD1;
	    g[idC] += alpha[idD] * tD1;
	  }
	  /* Avoid numerical problems (eg log(0.0)) by adding some noise. */
	  g[idC] += (drand48() * 1.0e-6) +  1.0e-8;
	  *dstLL += freq[idC] * log(g[idC]);
	}
	/* Calculate the probability densities of the sub-populations
	 * which form the mixture, and the log-likelihood function.  */
	converged = 0;
	sumAlpha = 0.0;
	for(idD = 0; idD < nDbn; ++idD)
	{
	  nt[idD] = dt[idD] = vt[idD] = 0.0;
	  for(idC = 0; idC < nCls; ++idC)
	  {
	    tD0 = samOrg + (samItv * idC);
	    tD1 = *(*(f + idC) + idD) * freq[idC] / g[idC];
	    dt[idD] += tD1;
	    nt[idD] += tD0 * tD1;
	    tD2 = tD0 - mu[idD];
	    vt[idD] += tD1 * tD2 * tD2;
	  }
	  /* Calculate denominators and numerators of new estimates. */
	  newmean[idD] = nt[idD] / dt[idD];
	  newalpha[idD] = alpha[idD] * dt[idD] / sumFreq;
	  sumAlpha += newalpha[idD];
	  newsd[idD] = sqrt(vt[idD] / dt[idD]);
	  /* Convergence test. */
	  converged = fabs(oldLL - *dstLL) <= tol;
	  oldLL = *dstLL;
	  alpha[idD] = newalpha[idD];
	  mu[idD] = newmean[idD];
	  sd[idD] = newsd[idD];
	}
	/* Ensure alpha's are normalised. */
	tD0 = 1.0 / sumAlpha;
	for(idD = 0; idD < nDbn; ++idD)
	{
	  alpha[idD] *= tD0;
	}
      }
    }
    while(!converged && (algErr == ALG_ERR_NONE));
  }
  /* Free temporary storage. */
  if(f)
  {
    AlcDouble2Free(f);
  }
  if(newalpha)
  {
    AlcFree(newalpha);
  }
  if(newmean)
  {
    AlcFree(newmean);
  }
  if(newsd)
  {
    AlcFree(newsd);
  }
  if(vt)
  {
    AlcFree(vt);
  }
  if(nt)
  {
    AlcFree(nt);
  }
  if(g)
  {
    AlcFree(g);
  }
  return(algErr);
}

/*!
* \return	Error code.
* \ingroup     AlgMixture
* \brief	Synthesise a mixture of normal distributions.
* \param	nCls			Number of classes in the
*					frequency distribution.
* \param	synFreq			Buffer for synthesised
*					frequencies.
* \param	nObv			Number of observations.
* \param	synOrg			Origin of synFreq buffer.
* \param	synStep			Sample interval for the synFreq
*					buffer.
* \param	nDbn			Number of distributions in the
*					mixture.
* \param	alpha			Estimates of nDbn alpha values.
* \param	mu			Estimates of nDbn mu values.
* \param	sigma			Estimates of nDbn standard
*					deviation values.
*/
AlgError	AlgMixtureSyn(int nCls, int *synFreq, int nObv,
			      double synOrg, double synStep,
			      int nDbn,
			      double *alpha, double *mu, double *sigma)
{
  int		idC,
  		idD;
  double	tD0,
  		synPos,
  		synVal,
		synTot;
  AlgError	algErr = ALG_ERR_NONE;

  if((synFreq == NULL) || 
     (alpha == NULL) || (mu == NULL) || (sigma == NULL) ||
     (nCls <= 0) || (nDbn <= 0))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    synTot = 0.0;
    for(idC = 0; idC < nCls; ++idC)
    {
      synVal = 0.0;
      synPos = synOrg + (synStep * idC);
      for(idD = 0; idD < nDbn; ++idD)
      {
	tD0 = (synPos - *(mu + idD)) / *(sigma + idD);
	synVal +=  *(alpha + idD) *
		   exp(-0.5 * tD0 * tD0) / sqrt(*(sigma + idD));
      }
      synTot += *(synFreq + idC) = ALG_NINT(synVal);
    }
    if(synTot > DBL_EPSILON)
    {
      tD0 = nObv / synTot;
      for(idC = 0; idC < nCls; ++idC)
      { 
	synVal = tD0 * *(synFreq + idC);
	*(synFreq + idC) = ALG_NINT(synVal);
      }
    }
  }
  return(algErr);
}
