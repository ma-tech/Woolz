#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzCompThresh.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for computing a threshold value from a histogram
* 		object.
* \ingroup	WlzThreshold
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static double			WlzCompThreshDepth(
				  WlzHistogramDomain *hDom);
static double			WlzCompThreshFoot(
				  WlzHistogramDomain *hDom);
static double			WlzCompThreshFracMin(
				  WlzHistogramDomain *hDom,
				  double bFrac,
				  WlzThresholdType *dstThrType);
static double			WlzCompThreshGradient(
				  WlzHistogramDomain *hDom);
static double			WlzCompThreshLSqFit(
				  WlzGreyP vec,
				  int idx0,
				  int idx1,
				  WlzObjectType histDomtype,
				  double *dstSlope);

/*!
* \return	Woolz error code.
* \ingroup	WlzThreshold
* \brief	Computes a threshold value and type using the given
*		histogram and method.
* \param	hObj
* \param	method			Threshold value method.
* \param	param0			First passed parameter:
*					<ul>
*					  <li>
*					  WLZ_COMPTHRESH_FRACMIN = minimum
*					  fraction of values which are
*					  background.
*					</ul>
*					Otherwise not used.
* \param	param1			First passed parameter:
*					  Not used.
* \param	extraFrac		Extra fraction to be added or
*					subtracted from the threshold value.
* \param	dstTV			Destination pointer for threshold
*					value, may be NULL.
* \param	dstTType		Destination pointer for threshold
*					type, may be NULL.
*/
WlzErrorNum	WlzCompThresholdVT(WlzObject *hObj, WlzCompThreshType method,
				   double param0, double param1,
				   double extraFrac, WlzPixelV *dstTV,
				   WlzThresholdType *dstTType)
{
  int		hMaxI;
  WlzPixelV	tV;
  WlzThresholdType tType;
  WlzHistogramDomain *hDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tV.type = WLZ_GREY_DOUBLE;
  if(hObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(hObj->type != WLZ_HISTOGRAM)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((hDom = hObj->domain.hist) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((hDom->type != WLZ_HISTOGRAMDOMAIN_INT) &&
          (hDom->type != WLZ_HISTOGRAMDOMAIN_FLOAT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((hDom->nBins < 0) || (hDom->maxBins < 0) ||
          (hDom->binSize <= DBL_EPSILON))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    switch(method)
    {
      case WLZ_COMPTHRESH_FOOT:  /* FALLTHROUGH */
      case WLZ_COMPTHRESH_DEPTH: /* FALLTHROUGH */
      case WLZ_COMPTHRESH_GRADIENT:
	hMaxI = WlzHistogramBinMax(hDom);
	if(hMaxI > (((hDom->nBins * hDom->binSize) + hDom->origin) / 2))
	{
	  /* These algorithms only work if the background peak is low. */
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	}
	else
	{
	  tType = WLZ_THRESH_HIGH;
	  switch(method)
	  {
	    case WLZ_COMPTHRESH_FOOT:
	      tV.v.dbv = WlzCompThreshFoot(hDom);
	      break;
	    case WLZ_COMPTHRESH_DEPTH:
	      tV.v.dbv = WlzCompThreshDepth(hDom);
	      break;
	    case WLZ_COMPTHRESH_GRADIENT:
	      tV.v.dbv = WlzCompThreshGradient(hDom);
	      break;
	  }
	}
	break;
      case WLZ_COMPTHRESH_FRACMIN:
        tV.v.dbv = WlzCompThreshFracMin(hDom, param0, &tType);
        break;
      default:
        errNum = WLZ_ERR_COMPTHRESH_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(tType == WLZ_THRESH_LOW)
    {
      tV.v.dbv -= extraFrac * tV.v.dbv;
    }
    else
    {
      tV.v.dbv += extraFrac * tV.v.dbv;
    }
    if(hDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      (void )WlzValueConvertPixel(&tV, tV, WLZ_GREY_INT);
    }
    if(dstTV)
    {
      *dstTV = tV;
    }
    if(dstTType)
    {
      *dstTType = tType;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup      WlzThreshold
* \brief	Computes a threshold value from the given histogram
*		object using the given method.
* \param	dstThrVal		Destination pointer for the threshold
*					value.
* \param	histObj			Given histogram object.
* \param	method			Given method.
* \param	extraFrac		Extra fraction added onto the
*					threshold value, ie:
*					thrVal += thrVal * extraFrac.
*/
WlzErrorNum	WlzCompThreshold(double *dstThrVal, WlzObject *histObj,
				 WlzCompThreshType method,
				 double extraFrac)
{
  double	thrVal;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	minBgdFrac = 0.20;

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
      case WLZ_COMPTHRESH_FRACMIN:
	thrVal = WlzCompThreshFracMin(histDom, minBgdFrac, NULL);
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

/*!
* \return	Threshold value.
* \ingroup 	WlzThreshold
* \brief	Computes both a threshold value and type using
*		a fraction and minimum algorithm in which:
*		<ul>
*		<li>
*		  The histogram maximum is found and assumed to
*		  be background. The position of the maximum in 
*		  the histogram determines the threshold type.
*		<li>
*		  The first minimum in th histogram after the
*		  given fraction of background values is the
*		  threshold value.
*		</ul>
* \param	hDom			Histogram domain.
* \param	bFrac			Minimum fraction of values
*					which are background.
* \param	dstThrType		Destination pointer for
*					threshold type, may be NULL.
*/
static double	WlzCompThreshFracMin(WlzHistogramDomain *hDom,
				     double bFrac,
				     WlzThresholdType *dstThrType)
{
  int		idH,
  		hMaxI;
  double	hC,
  		hF,
		hV0,
		hV1,
		thrVal,
  		hMaxV;
  WlzThresholdType thrType;

  /* Find histogram maximum. */
  hMaxI = WlzHistogramBinMax(hDom);
  hMaxV = (hDom->type == WLZ_HISTOGRAMDOMAIN_INT)?
	   *(hDom->binValues.inp + hMaxI):
	   *(hDom->binValues.dbp + hMaxI);
  hF = WlzHistogramBinSum(hDom);
  hF *= bFrac;
  thrType = (hMaxI > (((hDom->nBins * hDom->binSize) + hDom->origin) / 2))?
          WLZ_THRESH_LOW: WLZ_THRESH_HIGH;
  if(thrType == WLZ_THRESH_HIGH)
  {
    idH = 0;
    if(hDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      hC = *(hDom->binValues.inp + 0);
      while(hC < hF)
      {
	hC += *(hDom->binValues.inp + ++idH);
      }
      if(idH < hMaxI)
      {
	idH = hMaxI;
      }
      hV0 = *(hDom->binValues.inp + idH);
      while((++idH < hDom->nBins) &&
            (hV0 > (hV1 = *(hDom->binValues.inp + idH))))
      {
        hV0 = hV1;
      }
    }
    else /* hDom->type == WLZ_HISTOGRAMDOMAIN_FLOAT */
    {
      hC = *(hDom->binValues.dbp + 0);
      while(hC < hF)
      {
	hC += *(hDom->binValues.dbp + ++idH);
      }
      if(idH < hMaxI)
      {
	idH = hMaxI;
      }
      hV0 = *(hDom->binValues.dbp + idH);
      while((++idH < hDom->nBins) &&
            (hV0 > (hV1 = *(hDom->binValues.dbp + idH))))
      {
        hV0 = hV1;
      }
    }
    if(idH >= hDom->nBins)
    {
      idH = hDom->nBins;
    }
  }
  else /* thrType == WLZ_THRESH_LOW */
  {
    idH = hDom->nBins - 1;
    if(hDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      hC = *(hDom->binValues.inp + hDom->nBins - 1);
      while(hC < hF)
      {
	hC += *(hDom->binValues.inp + --idH);
      }
      if(idH > hMaxI)
      {
	idH = hMaxI;
      }
      hV0 = *(hDom->binValues.inp + idH);
      while((--idH >= 0) && (hV0 > (hV1 = *(hDom->binValues.inp + idH))))
      {
        hV0 = hV1;
      }
    }
    else /* hDom->type == WLZ_HISTOGRAMDOMAIN_FLOAT */
    {
      hC = *(hDom->binValues.dbp + hDom->nBins - 1);
      while(hC < hF)
      {
	hC += *(hDom->binValues.dbp + --idH);
      }
      if(idH > hMaxI)
      {
	idH = hMaxI;
      }
      hV0 = *(hDom->binValues.dbp + idH);
      while((--idH >= 0) && (hV0 > (hV1 = *(hDom->binValues.dbp + idH))))
      {
        hV0 = hV1;
      }
    }
    if(idH < 0)
    {
      idH = 0;
    }
  }
  if(dstThrType)
  {
    *dstThrType = thrType;
  }
  thrVal = hDom->origin + (idH * hDom->binSize);
  return(thrVal);
}

/*!
* \return	Threshold value.
* \ingroup 	WlzThreshold
* \brief	Computes a threshold value from the given histogram
*               domain using the WLZ_COMPTHRESH_FOOT method, ie:
*               The threshold value is intercept of a line fitted to
*               the upper slope of the histogram main peak with the
*               abscissa.
*		<ul>
*		<li>
*               The line is fitted from 90% of peak to 33% of peak.
*		<li>
*               Zero points in histogram are ignored.
*		<li>
*               The histogram may need to be smoothed for this
*               algorithm to work.
*		</ul>
* \param	histDom			Given histogram domain.
*/
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

/*!
* \return	Woolz error code.
* \ingroup	WlzThreshold
* \brief	Computes a threshold value from the given histogram
*               domain using the WLZ_COMPTHRESH_DEPTH method, ie:
*               The threshold value is that point to the right of the
*               histogram peak that is maximally distant from the chord
*               joining the peak and the histogram right hand end
*               point.
*               For each point from maxIdx to histMaxIdx the depth
*               below the chord joining these thwo points is computed
*               to find the index of the maximum depth.
*               The depth d is given by:
*                 \f$d = \frac{2A}{L}\f$
*               Where A is the area of the triangle with the chord at
*               it's base and histogram point at its apex, and L is
*               the length of the cord. Since the search only needs
*               to perform a comparison only the area A needs to be
*               computed.
*               The area A is given by:
*		\f[
                  \frac{|{(P_x - H_x)(M_y - H_y) - (P_y - H_y)(M_x - H_x)}|}
		       {2}
		\f]
*               Where   P is the histogram maximum (peak),
*                       H is the point on the histogram and
*                       M is the point at the histogram maximum.
*               The given histogram may need to be smoothed for this
*               algorithm to work.
* \param	histDom			Given histogram domain.
*/
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

/*!
* \return	Threshold value.
* \ingroup	WlzThreshold
* \brief	Computes a threshold value from the given histogram
*               domain using the WLZ_COMPTHRESH_GRADIENT method, ie:
*               The threshold value is the first point to the right
*               of the histogram main peak at which the gradient falls
*               to zero (cf finding a minimum).
*               To find the slope of the histogram at some point a
*               straight line is fitted through the point \f$\pm\f$ a fixed
*               number of points to either side.
*               The histogram may need to be smoothed for this
*               algorithm to work.
* \param	histDom			Given histogram domain.
*/
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

/*!
* \return	Intercept of fitted line.
* \ingroup	WlzThreshold
* \brief	Computes the least squares fit of a straight line to
*               the integer or double values in the given vector
*               between the two given indicies.
*               With the straight line of the form
*		\f[
                  y = a x + b
		\f]
*               then
*		\f[
                  a = \frac{(s_{xy} s_x) - (s_{xx} s_y)}
		           {(s_x^2) - (n s_{xx})}
		\f]
*		\f[
                  b = \frac{(s_x s_y) - (n s_{xy})}{(s_x^2) - (n s_{xx})}
		\f]
*               Where
*               \f[
                  s_x =  \sum_i{x_i}, s_y = \sum_i{y_i}
		\f]
*		\f[
                  s_{xx} = \sum_i{x_i^2}, s_{yy} = \sum_i{y_i^2}
		\f]
*		\f[
                  s_{xy} = \sum_i{x_i y_i}
		\f]
*               If \f$\{x_i\}\f$ are known to be integers of the form
*               \f$\{x_0, x_0 + 1, ..., x_0 + n - 1\}\f$
*               then
*		\f[
                  \sum_i{x_i} = \frac{n (2 x_0 + n - 1)}{2}
		\f]
*		\f[
                  \sum_i{x_i^2} = \frac{n}{6} - n x0 - \frac{n^2}{2} +
                                  \frac{n^3}{3} + n x_0^2 + n^2 x_0
		\f]
*               This gives
*		\f[
                  a = \frac{6 s_{xy} (1  - n - 2 x_0) +
                            2 s_y (1 - 3 n + 2 n^2) +
                            6 x_0 (-1 + n + x_0)}
                           {n (n^2 - 1)}
		\f]
*		\f[
                  b = \frac{12 s_{xy} + 6 s_y (1 - n - 2 x_0)}
                           {n (n^2 - 1)}
		\f]
* \param	vec			Vector of values.
* \param	idx0			Lowest of the two given search range
* 					indicies.
* \param	idx1			Highest of the two given search range
* 					indicies.
* \param	histDomtype		Histogram domain type which specifies
* 					int or double values.
* \param	dstSlope		Destination pointer for the slope
* 					of the fitted line.
*/
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

