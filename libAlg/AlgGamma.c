#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgGamma.c
* \author       Bill Hill
* \date         May 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for computing gamma and incomplete gamma functions.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgGamma
* @{
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

static double	AlgGammaCF(double a, double x, AlgError *dstErr);
static double	AlgGammaS(double a, double x, AlgError *dstErr);

/*!
* \return				Log of gamma function value.
* \brief	Computes the log gamma function log(Gamma(x)), ie
*		exp(AlgGammaLog(n + 1)) = n!.
*		This function is based on the function gammln():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
* \param   	x			Given value.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
double		AlgGammaLog(double x, AlgError *dstErr)
{
  int 		j;
  double	xx,
		y,
		tmp,
		ser,
		lnGam;
  AlgError	algErr = ALG_ERR_NONE;
  const double	cof[6]={76.18009172947146,-86.50532032941677,
  			24.01409824083091,-1.231739572450155,
  			0.1208650973866179e-2,-0.5395239384953e-5};

  if(x <= 0.0)
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    y = xx = x;
    tmp = xx + 5.5;
    tmp -= (xx + 0.5) * log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++)
    {
      y += 1.0;
      ser += cof[j] / y;
    }
    lnGam = -tmp + log(2.5066282746310005 * ser / xx);
  }
  if(dstErr)
  {
    *dstErr = algErr;
  }
  return(lnGam);
}

/*!
* \return				Incomplete gamma function value.
* \brief	Computes the incomplete gamma function P(a,x), which
*		has the limiting values P(a,0) = 0, P(a, oo) = 1..
*		This function is based on the function gammp():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
* \param	a			Incomplete gamma fn parameter.
* \param	x			Incomplete gamma fn parameter.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
double		AlgGammaP(double a, double x, AlgError *dstErr)
{
  double	gamma = 0.0;
  AlgError	algErr = ALG_ERR_NONE;

  if((x < 0.0) || (a < DBL_EPSILON))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    if(x < (a + 1.0))
    {
      gamma = 1.0 - AlgGammaS(a, x, &algErr);
    }
    else
    {
      gamma = AlgGammaCF(a, x, &algErr);
    }
  }
  if(dstErr)
  {
    *dstErr = algErr;
  }
  return(gamma);
}

/*!
* \return				Incomplete gamma function value.
* \brief	Computes the incomplete gamma function P(a,x), using
*		a series method.
*		This function is based on the function gser():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
* \param	a			Incomplete gamma fn parameter.
* \param	x			Incomplete gamma fn parameter.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
double		AlgGammaS(double a, double x, AlgError *dstErr)
{
  int 		cnt;
  double	ap,
  		del,
		gln,
		sum,
		gamma = 0.0;
  const int	maxIt = 100;
  const double	eps = 1.0e-06;
  AlgError	algErr = ALG_ERR_NONE;

  if((x < 0.0) || (a < DBL_EPSILON))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    if(x > DBL_EPSILON)
    {
      ap = a;
      cnt = maxIt;
      del = sum = 1.0 / a;
      do
      {
	++ap;
	del *= x / ap;
	sum += del;
      } while((fabs(del) > (fabs(sum) * eps)) && (--cnt > 0));
      if(cnt > 0)
      {
    	gln = AlgGammaLog(a, NULL);
	gamma = sum * exp((a * log(x)) - (gln + x));
      }
      else
      {
	algErr = ALG_ERR_CONVERGENCE;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = algErr;
  }
  return(gamma);
}

/*!
* \return				Incomplete gamma function value.
* \brief	Computes the incomplete gamma function P(a,x) using
*		continued fractions.
*		This function is based on the function gcf():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
* \param	a			Incomplete gamma fn parameter.
* \param	x			Incomplete gamma fn parameter.
* \param	dstErr			Destination ptr for error code,
*					may be NULL.
*/
static double	AlgGammaCF(double a, double x, AlgError *dstErr)
{
  int		idx;
  double	an,
  		b,
		c,
		d,
		del,
		h,
		gln,
		gamma = 0.0;
  AlgError	algErr = ALG_ERR_NONE;
  const int	maxIt = 100;
  const double	eps = 1.0e-06,
  		min = 1.0e-30;

  if((x < 0.0) || (a < DBL_EPSILON))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    b = x + 1.0 - a;
    c = 1.0 / min;
    d = 1.0 / b;
    h = d;
    idx = 1;
    do
    {
      an = -idx * (idx - a);
      b += 2.0;
      d = (an * d) + b;
      if(fabs(d) < min)
      {
        d = min;
      }
      c = b + (an / c);
      if(fabs(c) < min)
      {
        c = min;
      }
      d = 1.0 / d;
      del = d * c;
      h *= del;
    } while((fabs(del - 1.0) > eps) && (++idx <= maxIt));
    if(idx <= maxIt)
    {
      gln = AlgGammaLog(a, NULL);
      gamma = h * exp((a * log(x)) - (gln + x));
    }
    else
    {
      algErr = ALG_ERR_CONVERGENCE;
    }
  }
  if(dstErr)
  {
    *dstErr = algErr;
  }
  return(gamma);
}

#ifdef ALG_GAMMA_TEST

int		main(int argc, char **argv)
{
  int		ok = 0;
  double	a = 0.0,
		g,
		gln,
  		x;
  AlgError	algErr = ALG_ERR_NONE;

  switch(argc)
  {
    case 2:
      if((sscanf(*(argv + 1), "%lg", &x) == 1) && (x >= 0.0))
      {
        ok = 1;
      }
      break;
    case 3:
      if((sscanf(*(argv + 1), "%lg", &x) == 1) && (x >= 0.0) &&
         (sscanf(*(argv + 2), "%lg", &a) == 1) && (a >= DBL_EPSILON))
      {
        ok = 1;
      }
      break;
    default:
      ok = 0;
      break;
  }
  if(ok == 0)
  {
    (void )fprintf(stderr,
    "Usage: %s <x> [<a>]\n%s\n",
    *argv,
    "A test for AlgGammaP() and AlgGammaLog() which computes either the\n"
    "gamma function P(x) from a single value, or the incomplete gamma\n"
    "function P(a,x) from a pair of values.\n");
  }
  else
  {
    if(a >= DBL_EPSILON)
    {
      g = AlgGammaP(a, x, &algErr);
    }
    else
    {
      gln = AlgGammaLog(x, &algErr);
      g = exp(gln);
    }
    if(algErr == ALG_ERR_NONE)
    {
      printf("%g\n", g);
    }
    else
    {
      (void )fprintf(stderr, "%s: Failed to compute gamma function.\n", *argv);
    }
  }
  return(!ok);
}

#endif /* ALG_GAMMA_TEST */

/*!
* @}
*/
