#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgGamma_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgGamma.c
* \author       Bill Hill
* \date         May 2000
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
* \brief	Functions for computing gamma and incomplete gamma functions.
* \ingroup      AlgGamma
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

static double	AlgGammaCF(double a, double x, AlgError *dstErr);
static double	AlgGammaS(double a, double x, AlgError *dstErr);

/*!
* \return	Log of gamma function value.
* \ingroup	AlgGamma
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
		lnGam = 1.0;
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
* \return	Incomplete gamma function value.
* \ingroup     AlgGamma
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
* \return	Incomplete gamma function value.
* \ingroup     AlgGamma
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
* \return	Incomplete gamma function value.
* \ingroup     AlgGamma
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
