#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgMixture.c
* Date:         January 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Provides a function for computing the maximum liklihood
*		parameters of a mixture of distributions which fit the
*		given data for the MRC Human Genetics Unit numerical
*		algorithm library.
*		This algorithm is based on the article:
*		Agha M. and Ibrahim M.T.
*		Maximum Liklihood Estimation of Mixtures of Distributions
*	 	Applied Statistics 33(3):327-332, 1983.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
************************************************************************/
#include <Alg.h>
#include <math.h>
#include <float.h>

/************************************************************************
* Function:	AlgMixture
* Returns:	AlgError:		Error code.
* Purpose:	Computes the maximum liklihood estimate of the
*		parameters of a mixture of normal, exponential,
*		Poisson or binomial distributions which best fits
*		the give frequencies. These parameters are the miximg
*		proportions, means, and (for normal distributions)
*		standard deviations. The log-liklihood and the number
*		of iterations are also computed.
*		The frequencies are sampled at unit intervals.
* Global refs:	-
* Parameters:	AlgDistribution dbnType	Type of distribution.
*		int nDbn:		Number of distributions in the
*					mixture.
*		int nCls:		Number of classes in the
*					frequency distribution.
*		double *x:		Position of given frequencies.
*		int *freq:		Given frequencies.
*		double *alpha:		Estimates of nDbn alpha values,
*					both given and returned.
*				 	Alpha values must be in the
*					range (0,1].
*		double *mean:		Estimates of nDbn mean values,
*					both given and returned.
*		double *sd:		Estimates of nDbn standard
*					deviation values,
*					both given and returned.
*		double tol:		Difference between two
*					consecutive log liklihood values
*					required to terminate iteration.
*		int nObv:		Number of observations.
*		double *dstLL:		Destination pointer for the
*					log liklihood.
*		double *nItn:		Destination pointer for the
*					number of iterations taken.
************************************************************************/
AlgError	AlgMixture(AlgDistribution dbnType, int nDbn, int nCls,
    			   double *x, int *freq, double *alpha, double *mean,
			   double *sd, double tol, int nObv, double *dstLL,
			   int *nItn)
{
  int		idC,
  		idD,
		idD1,
		converged;
  double	tD0,
  		tD1,
		tD2,
		part,
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
  switch(dbnType)
  {
    case ALG_DISTRIBUTION_NORMAL:
    case ALG_DISTRIBUTION_EXP:
    case ALG_DISTRIBUTION_POISSON:
    case ALG_DISTRIBUTION_BINOMIAL:
      break;
    default:
      algErr = ALG_ERR_FUNC;
      break;
  }
  if(algErr == ALG_ERR_NONE)
  {
    for(idC = 1; idC < nCls; ++idC)
    {
      if(x[idC - 1] > x[idC])
      {
	algErr = ALG_ERR_FUNC;
      }
    }
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(nObv < (nCls * 2))
    {
      algErr = ALG_ERR_FUNC;
    }
  }
  if(algErr == ALG_ERR_NONE)
  {
    for(idC = 0; idC < nCls; ++idC)
    {
      if (freq[idC] < 0)
      {
	algErr = ALG_ERR_FUNC;
      }
    }
  }
  if((algErr == ALG_ERR_NONE) && (dbnType != ALG_DISTRIBUTION_NORMAL))
  {
    for(idC = 0; idC < nCls; ++idC)
    {
      if(x[idC] < 0.0)
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
    *nItn = 0;
    do
    {
      ++*nItn;
      /* Check alpha, mean and standard deviation. */
      for(idD = 0; idD < nDbn; ++idD)
      {
	if(isnan(alpha[idD]) ||
	   (alpha[idD]) > 1.0 || (alpha[idD] < 0.0))
	{
	  algErr = ALG_ERR_CONVERGENCE;
	}
	else if(isnan(mean[idD]) ||
		(mean[idD] >= x[nCls - 1]) || (mean[idD] <= x[0]))
	{
	  algErr = ALG_ERR_CONVERGENCE;
	}
      }
      if((algErr == ALG_ERR_NONE) && (dbnType == ALG_DISTRIBUTION_NORMAL))
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
	    if(mean[idD] == mean[idD1])
	    {
	      if(dbnType == ALG_DISTRIBUTION_NORMAL)
	      {
		if(sd[idD] == sd[idD1])
		{
		  algErr = ALG_ERR_FUNC;
		}
	      }
	      else
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
	    switch(dbnType)
	    {
	      case ALG_DISTRIBUTION_NORMAL:
		tD0 = (x[idC] - mean[idD]) / sd[idD];
		tD1 = exp(-0.5 * (tD0 * tD0)) / sd[idD];
		break;
	      case ALG_DISTRIBUTION_EXP:
		tD1 = exp(-x[idC] / mean[idD]) / mean[idD];
		break;
	      case ALG_DISTRIBUTION_POISSON:
		if (x[idC] == x[0])
		{
		  tD1 = exp(-mean[idD]) * pow(mean[idD], x[idC]);
		}
		else 
		{
		  tD1 = *(*(f + idC - 1) + idD) * mean[idD];
		}
		break;
	      case ALG_DISTRIBUTION_BINOMIAL:
		if (x[idC] == x[0])
		{
		  tD1 = pow(1.0 - (mean[idD] / x[nCls - 1]), x[nCls - 1]) *
			pow(mean[idD - 1] / (x[nCls - 1] - mean[idD - 1]),
			    x[idC]);
		}
		else
		{
		  tD1  = *(*(f + idC - 1) + idD) *
		  	 (mean[idD] / (x[nCls - 1] - mean[idD]));
		}
		break;
	    }
	    *(*(f + idC) + idD) = tD1;
	    g[idC] += alpha[idD] * tD1;
	  }
	  /* Avoid numerical problems (eg log(0.0)) by adding some noise. */
	  g[idC] += (drand48() * 1.0e-10) +  1.0e-12;
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
	    part = *(*(f + idC) + idD) * freq[idC] / g[idC];
	    dt[idD] += part;
	    nt[idD] += part * x[idC];
	    if(dbnType == ALG_DISTRIBUTION_NORMAL)
	    {
	      tD0 = x[idC] - mean[idD];
	      vt[idD] += part * tD0 * tD0;
	    }
	  }
	  /* Calculate denominators and numerators of new estimates. */
	  newmean[idD] = nt[idD] / dt[idD];
	  if(idD != nDbn - 1)
	  {
	    newalpha[idD] = alpha[idD] * dt[idD] / nObv;
	    sumAlpha += newalpha[idD];
	  }
	  else
	  {
	    newalpha[nDbn - 1] = 1.0 - sumAlpha;
	  }
	  if(dbnType == ALG_DISTRIBUTION_NORMAL)
	  {
	    newsd[idD] = sqrt(vt[idD] / dt[idD]);
	  }
	  /* Convergence test. */
	  converged = fabs(oldLL - *dstLL) <= tol;
	  oldLL = *dstLL;
	  alpha[idD] = newalpha[idD];
	  mean[idD] = newmean[idD];
	  if(dbnType == ALG_DISTRIBUTION_NORMAL)
	  {
	    sd[idD] = newsd[idD];
	  }
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
