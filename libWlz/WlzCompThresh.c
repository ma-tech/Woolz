#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzCompThresh.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for computing a threshold value from a
*		histogram object.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static double	WlzCompThreshFoot(WlzHistogramDomain *),
		WlzCompThreshDepth(WlzHistogramDomain *),
		WlzCompThreshGradient(WlzHistogramDomain *);
static double	WlzCompThreshLSqFit(WlzGreyP vec, int idx0, int idx1,
				     WlzObjectType histDomtype,
				     double *dstSlope);

/************************************************************************
* Function:	WlzCompThreshold					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Computes a threshold value from the given histogram	*
*		object using the given method.				*
* Global refs:	-							*
* Parameters:	double *dstThrVal:	Destination pointer for 	*
*					threshold value.		*
*		WlzObject *histObj:	Given histogram object.		*
*		WlzCompThreshType method: Given method.			*
*		double extraFrac:	Extra fraction added onto the	*
*					threshold value, ie:		*
*					thrVal += thrVal * extraFrac.	*
************************************************************************/
WlzErrorNum	WlzCompThreshold(double *dstThrVal, WlzObject *histObj,
				 WlzCompThreshType method,
				 double extraFrac)
{
  double	thrVal;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzCompThreshold FE 0x%lx 0x%lx %d %g\n",
	   (unsigned long )dstThrVal, (unsigned long )histObj,
	   (int )method, extraFrac));
  if(histObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(histObj->type != WLZ_HISTOGRAM)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((histDom = histObj->domain.hist) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((histDom->type != WLZ_HISTOGRAMDOMAIN_INT) &&
          (histDom->type != WLZ_HISTOGRAMDOMAIN_FLOAT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((histDom->nBins < 0) || (histDom->maxBins < 0) ||
          (histDom->binSize <= DBL_EPSILON))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if(dstThrVal == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(method)
    {
      case WLZ_COMPTHRESH_FOOT:
	thrVal = WlzCompThreshFoot(histDom);
        break;
      case WLZ_COMPTHRESH_DEPTH:
	thrVal = WlzCompThreshDepth(histDom);
        break;
      case WLZ_COMPTHRESH_GRADIENT:
	thrVal = WlzCompThreshGradient(histDom);
        break;
      default:
        errNum = WLZ_ERR_COMPTHRESH_TYPE;
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      thrVal += extraFrac * thrVal;
      *dstThrVal = thrVal;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzCompThreshold FX %d\n",
	   (int )errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzCompThreshFoot					*
* Returns:	double:			Threshold value.		*
* Purpose:	Computes a threshold value from the given histogram	*
*		domain using the WLZ_COMPTHRESH_FOOT method, ie:	*
*		The threshold value is intercept of a line fitted to	*
*		the upper slope of the histogram main peak with the	*
*		abscissa.						*
*		* The line is fitted from 90% of peak to 33% of peak.	*
* 		* Zero points in histogram are ignored.			*
*		* The histogram may need to be smoothed for this	*
*		  algorithm to work.					*
* Global refs:	-							*
* Parameters:	WlzHistogramDomain *histDom: Given histogram domain.	*
************************************************************************/
static double	WlzCompThreshFoot(WlzHistogramDomain *histDom)
{
  double	tD0,
		lsqA,
		lsqB,
  		maxVal,
		thrVal;
  int		tI0,
		binCount,
		topIdx,
		botIdx,
  		maxIdx;
  WlzGreyP	binVal;

  maxIdx = WlzHistogramBinMax(histDom);
  maxVal = (histDom->type == WLZ_HISTOGRAMDOMAIN_INT)?
	   *(histDom->binValues.inp + maxIdx):
	   *(histDom->binValues.dbp + maxIdx);
  /* Find the bin at which the histogram first drops to 90% of the maximum.
     If the drop is to 30% then take previous bin.  */
  topIdx = maxIdx + 1;
  binCount = histDom->nBins;
  switch(histDom->type)
  {
    case WLZ_HISTOGRAMDOMAIN_INT:
      tI0 = maxVal * 0.9;
      binVal.inp = histDom->binValues.inp + topIdx;
      while((topIdx < binCount) && (*(binVal.inp) > tI0)) 
      {
        ++topIdx;
	++(binVal.inp);
      }
      if(*(binVal.inp) <= (maxVal * 0.3))
      {
        --topIdx;
      }
      break;
    case WLZ_HISTOGRAMDOMAIN_FLOAT:
      tD0 = maxVal * 0.9;
      binVal.dbp = histDom->binValues.dbp + topIdx;
      while((topIdx < binCount) && (*(binVal.dbp) > tD0)) 
      {
        ++topIdx;
	++(binVal.dbp);
      }
      if((*(binVal.dbp) <= (maxVal * 0.3)) ||
         (topIdx == binCount))
      {
        --topIdx;
      }
      break;
  }
  /* Find the bin at which the histogram first drops to 33% of the maximum. */
  botIdx = topIdx + 1;
  switch(histDom->type)
  {
    case WLZ_HISTOGRAMDOMAIN_INT:
      tI0 = maxVal * 0.33;
      binVal.inp = histDom->binValues.inp + botIdx;
      while((botIdx < binCount) && (*(binVal.inp) > tI0)) 
      {
        ++topIdx;
	++(binVal.inp);
      }
      if(botIdx == binCount)
      {
        --botIdx;
      }
      break;
    case WLZ_HISTOGRAMDOMAIN_FLOAT:
      tD0 = maxVal * 0.33;
      binVal.dbp = histDom->binValues.dbp + botIdx;
      while((botIdx < binCount) && (*(binVal.dbp) > tD0)) 
      {
        ++topIdx;
	++(binVal.dbp);
      }
      if(botIdx == binCount)
      {
        --botIdx;
      }
      break;
  }
  /* Compute least squares fit of line to histogram values between topIdx
     and botIdx */
  if(((botIdx - topIdx) <= 2) ||
     (((lsqA = WlzCompThreshLSqFit(histDom->binValues, topIdx, botIdx,
    			           histDom->type, &lsqB)) <= DBL_EPSILON) &&
      (lsqA >= -(DBL_EPSILON))) ||
     ((lsqB <= DBL_EPSILON) && (lsqB >= -(DBL_EPSILON))))
  {
    thrVal = histDom->origin + (botIdx * histDom->binSize);
  }
  else
  {
    thrVal = histDom->origin - ((lsqA * histDom->binSize) / lsqB);
  }
  return (thrVal) ;
}

/************************************************************************
* Function:	WlzCompThreshDepth					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Computes a threshold value from the given histogram	*
*		domain using the WLZ_COMPTHRESH_DEPTH method, ie:	*
*		The threshold value is that point to the right of the	*
*		histogram peak that is maximally distant from the chord	*
*		joining the peak and the histogram right hand end 	*
*		point.							*
*		For each point from maxIdx to histMaxIdx the depth	*
*		below the chord joining these thwo points is computed	*
*		to find the index of the maximum depth.			*
*		The depth d is given by:				*
*		  d = 2A / L						*
*		Where A is the area of the triangle with the chord at	*
*		it's base and histogram point at its apex, and L is	*
*		the length of the cord. Since the search only needs	*
*		to perform a comparison only the area A needs to be	*
*		computed.						*
*		The area A is given by:					*
*		   |(Px - Hx)(My - Hy) - (Py - Hy)(Mx - Hx)| / 2	*
*		Where 	P is the histogram maximum (peak), 		*
*			H is the point on the histogram and		*
*			M is the point at the histogram maximum.	*
*		* The histogram may need to be smoothed for this	*
*		  algorithm to work.					*
* Global refs:	-							*
* Parameters:	WlzHistogramDomain *histDom: Given histogram domain.	*
************************************************************************/
static double	WlzCompThreshDepth(WlzHistogramDomain *histDom)
{
  double	depth,
  		depthMax,
  		endVal,
		newVal,
  		maxVal,
		thrVal;
  int		tI0,
  		tI1,
		endIdx,
		deepIdx,
  		newIdx,
  		maxIdx;
  WlzGreyP	binVal;

  maxIdx = WlzHistogramBinMax(histDom);
  /* Find the bin at which the histogram value is at the maximum depth below
     cord joining histogram peak and histogram end bin */
  endIdx = histDom->nBins - 1;
  deepIdx = 0;
  depthMax = 0.0;
  if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
  {
    binVal.inp = histDom->binValues.inp + maxIdx;
    maxVal = *(binVal.inp);
    endVal = *(histDom->binValues.inp + endIdx);
  }
  else
  {
    binVal.dbp = histDom->binValues.dbp + maxIdx;
    maxVal = *(binVal.dbp);
    endVal = *(histDom->binValues.dbp + endIdx);
  }
  newIdx = maxIdx + 1;
  tI0 = newIdx - maxIdx;
  tI1 = endIdx - newIdx;
  while(newIdx <= endIdx)
  {
    if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      newVal = *(binVal.inp)++;
    }
    else
    {
      newVal = *(binVal.dbp)++;
    }
    depth = fabs((tI0 * (endVal - newVal)) - (tI1 * (maxVal - newVal)));
    if(depth > depthMax)
    {
      depthMax = depth;
      deepIdx = newIdx;
    }
    ++tI0;
    --tI1;
    ++newIdx;
  }
  thrVal = histDom->origin + (deepIdx * histDom->binSize);
  return(thrVal);
}

/************************************************************************
* Function:	WlzCompThreshGradient					*
* Returns:	double:			Threshold value.		*
* Purpose:	Computes a threshold value from the given histogram	*
*		domain using the WLZ_COMPTHRESH_GRADIENT method, ie:	*
*		The threshold value is the first point to the right	*
*		of the histogram main peak at which the gradient falls	*
*		to zero (cf finding a minimum).				*
*		To find the slope of the histogram at some point a	*
*		straight line is fitted through the point +/- a fixed	*
*		number of points to either side.			*
*		* The histogram may need to be smoothed for this	*
*		  algorithm to work.					*
* Global refs:	-							*
* Parameters:	WlzHistogramDomain *histDom: Given histogram domain.	*
************************************************************************/
static double	WlzCompThreshGradient(WlzHistogramDomain *histDom)
{
  double	tD0,
		slope,
  		maxVal,
		thrVal;
  int		tI0,
		binCount,
		fitEndIdx,
		topIdx,
		endIdx,
  		maxIdx;
  WlzGreyP	binVal;
  WlzGreyP	fitBuf;
  const int	hFitWidth = 3,
  		fitWidth = (hFitWidth * 2) + 1;

  maxIdx = WlzHistogramBinMax(histDom);
  maxVal = (histDom->type == WLZ_HISTOGRAMDOMAIN_INT)?
	   *(histDom->binValues.inp + maxIdx):
	   *(histDom->binValues.dbp + maxIdx);
  /* Find the bin at which the histogram first drops to 60% of the maximum.
     If the drop is to 30% then take previous bin.  */
  topIdx = maxIdx + 1;
  binCount = histDom->nBins;
  switch(histDom->type)
  {
    case WLZ_HISTOGRAMDOMAIN_INT:
      tI0 = maxVal * 0.60;
      binVal.inp = histDom->binValues.inp + topIdx;
      while((topIdx < binCount) && (*(binVal.inp) > tI0)) 
      {
        ++topIdx;
	++(binVal.inp);
      }
      if(*(binVal.inp) <= (maxVal * 0.3))
      {
        --topIdx;
      }
      break;
    case WLZ_HISTOGRAMDOMAIN_FLOAT:
      tD0 = maxVal * 0.60;
      binVal.dbp = histDom->binValues.dbp + topIdx;
      while((topIdx < binCount) && (*(binVal.dbp) > tD0)) 
      {
        ++topIdx;
	++(binVal.dbp);
      }
      if((*(binVal.dbp) <= (maxVal * 0.3)) ||
         (topIdx == binCount))
      {
        --topIdx;
      }
      break;
  }
  /* Find the bin at which the histogram gradient drops to zero. */
  endIdx = histDom->nBins - 1;
  if(topIdx < hFitWidth)
  {
    topIdx = hFitWidth;
  }
  fitEndIdx = topIdx + hFitWidth;
  if(fitEndIdx < endIdx)
  {
    if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      binVal.inp = histDom->binValues.inp + topIdx - hFitWidth;
    }
    else
    {
      binVal.dbp = histDom->binValues.dbp + topIdx - hFitWidth;
    }
    slope = -1.0;
    while((fitEndIdx <= endIdx) && (slope <= DBL_EPSILON))
    {
      if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
        fitBuf.inp = (binVal.inp)++;
      }
      else
      {
        fitBuf.dbp = (binVal.dbp)++;
      }
      (void )WlzCompThreshLSqFit(fitBuf, 0, fitWidth - 1,
				 histDom->type, &slope);

      ++fitEndIdx;
    }
    if(slope >= -(DBL_EPSILON))
    {
      fitEndIdx -= 2;
    }
  }
  thrVal = histDom->origin + ((fitEndIdx - hFitWidth) * histDom->binSize);
  return(thrVal);
}

/************************************************************************
* Function:	WlzCompThreshLSqFit					*
* Returns:	double:			Intercept of fitted line.	*
* Purpose:	Computes the least squares fit of a straight line to	*
*		the integer or double values in the given vector	*
*		between the two given indicies.				*
*		With the straight line of the form			*
*		  y = a*x + b						*
*		then							*
*		  a = ((sxy*sx) - (sxx*sy)) / ((sx^2) - (n*sxx))	*
*		  b = ((sx*sy) - (n*sxy)) / ((sx^2) - (n*sxx))		*
*		Where							*
*		  sx =  Sum(xi), sy = Sum(yi)				*
*		  sxx = Sum(xi*xi), syy = Sum(yi*yi)			*
*		  sxy = Sum(xi*yi)					*
*		If {xi} are known to be integers of the form		*
*		  {x0, x0 + 1, ..., x0 + n - 1}				*
*		then							*
*		  sum(xi) = (n*((2*x0) + n - 1)) / 2			*
*		  sum(xi^2) = (n/6) - (n*x0) - ((n^2)/2) + ((n^3)/3) +	*
*			      (n*(x0^2)) + ((n^2)*x0)			*
*		This gives						*
*		  a = ((6*sxy*(1  - n - (2*x0))) + 			*
*		       (2*sy*(1 - (3*n) + (2*(n^2)) +			*
*			      (6*x0*(-1 + n + x0))))) / 		*
*		      (n*((n^2) - 1))					*
*		  b = ((12*sxy) + (6*sy*(1 - n - (2*x0)))) /		*
*		      (n*((n^2) - 1))					*
* Global refs:	-							*
* Parameters:	WlzGreyP vec:		Vector of values.		*
*		int idx0:		Lowest of the two given search	*
*					range indicies.			*
*		int idx1:		Highest of the two given search	*
*					range indicies.			*
*		WlzObjectType histDomtype: Histogram domain type which	*
*					specifies int or double values.	*
*		double *dstSlope:	Destination pointer for the 	*
*					slope of the fitted line.	*
************************************************************************/
static double	WlzCompThreshLSqFit(WlzGreyP vec, int idx0, int idx1,
				    WlzObjectType histDomtype,
				    double *dstSlope)
{
  int		tI0,
  		idx,
  		nSam,
		nSamSq;
  double	tD1,
		tD2,
  		tD3,
		lsqA,
  		lsqB,
		sXY = 0.0,
  		sY = 0.0;
  int		*tIP0;
  double	*tDP0;

  if(idx1 > idx0)
  {
    nSam = idx1 - idx0 + 1;
    nSamSq = nSam * nSam;
    idx = idx0;
    if(histDomtype == WLZ_HISTOGRAMDOMAIN_INT)
    {
      tIP0 = vec.inp + idx0;
      while(idx <= idx1)
      {
	sY += *tIP0;
	sXY  += idx * *tIP0;
	++tIP0;
	++idx;
      }
    }
    else
    {
      tDP0 = vec.dbp + idx0;
      while(idx <= idx1)
      {
	sY += *tDP0;
	sXY  += idx * *tDP0;
	++tDP0;
	++idx;
      }
    }
    tI0 = 2 * idx0;
    tD1 = 1 - nSam - (2 * tI0); 
    tD2 = 1.0 - (3.0 * nSam) + (2.0 * nSamSq) + 
	  (6.0 * tI0 * (-1 + nSam + tI0));
    tD3 = nSam*((nSamSq) - 1.0);
    lsqA = ((6.0 * sXY * tD1) + (2.0 * sY * tD2)) / tD3;
    lsqB = ((12.0 * sXY) + (6.0 * sY * tD1)) / tD3;
  }
  if(dstSlope)
  {
    *dstSlope = lsqB;
  }
  return(lsqA);
}
