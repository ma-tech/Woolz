#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzHistogram_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzHistogram.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Functions for computing and transforming Woolz histograms
* 		domains.
* \ingroup	WlzHistogram
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Checks the given objects domain and values if either is
*               NULL the appropriate error code is returned. Also finds
*               the grey type from the given object.
*               Because this is a static function the object pointer
*               doesn't need to be checked.
* \param	greyType		Destination pointer for grey
*                                       type.
* \param	isTiled			Destination pointer set to non-zero
* 					if the object has tiled values.
* \param	srcObj			Given source object.
*/
static WlzErrorNum WlzHistogramCheckDomainAndValues(WlzGreyType *greyType,
						    int *isTiled,
						    WlzObject *srcObj)
{
  int		planeCount;
  WlzValues	*planeValues;
  WlzGreyType	prvGreyType0,
  		prvGreyType1 = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    *isTiled = 0;
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        break;
      case WLZ_TRANS_OBJ:
        break;
      case WLZ_2D_DOMAINOBJ:
	*isTiled = WlzGreyTableIsTiled(srcObj->values.core->type);
	*greyType = WlzGreyTableTypeToGreyType(srcObj->values.core->type,
					       &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	*isTiled = WlzGreyTableIsTiled(srcObj->values.core->type);
	if(*isTiled)
	{
          *greyType = WlzGreyTypeFromObj(srcObj, &errNum);
	}
	else
	{
	  planeCount = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
	  planeValues = srcObj->values.vox->values;
	  while((planeCount-- > 0) && (errNum == WLZ_ERR_NONE))
	  {
	    if(((*planeValues).core) &&
	       ((*planeValues).core->type != WLZ_EMPTY_OBJ))
	    {
	      prvGreyType0 = WlzGreyTableTypeToGreyType(
						    (*planeValues).core->type,
						    &errNum);
	      if((prvGreyType1 != WLZ_GREY_ERROR) &&
		 (prvGreyType0 != prvGreyType1) &&
		 (errNum == WLZ_ERR_NONE))
	      {
		errNum = WLZ_ERR_VALUES_DATA;
	      }
	      prvGreyType1 = prvGreyType0;
	    }
	    ++planeValues;
	  }
	  if(prvGreyType1 == WLZ_GREY_ERROR)
	  {
	    errNum = WLZ_ERR_VALUES_TYPE;
	  }
	  *greyType = prvGreyType1;
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Checks the given histogram object for basic validity.
* \param	histObj			Given histogram object.
*/
static WlzErrorNum WlzHistogramCheckHistObj(WlzObject *histObj)
{
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
  else if((histDom->nBins < 0) || (histDom->maxBins < 0) ||
          (histDom->binSize <= DBL_EPSILON))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if((histDom->type != WLZ_HISTOGRAMDOMAIN_INT) &&
          (histDom->type != WLZ_HISTOGRAMDOMAIN_FLOAT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup      WlzHistogram
* \brief	Adds the bins of the second histogram domain to those
*               of the first. Both histogram domains are assumed to
*               be the same except for the bin values.
*               Because this is a static function the object pointer
*               doesn't need to be checked.
* \param	dstHistDom		Destination histogram domain.
* \param	srcHistDom		Source histogram domain.
*/
static void	WlzHistogramAddIntBins(WlzHistogramDomain *dstHistDom,
				       WlzHistogramDomain *srcHistDom)
{
  int		tI0;
  int		*tIP0,
		*tIP1;

  tIP0 = dstHistDom->binValues.inp;
  tIP1 = srcHistDom->binValues.inp;
  tI0 = dstHistDom->maxBins;
  while(tI0-- > 0)
  {
    *tIP0++ += *tIP1++;
  }
}

/*!
* \return	Woolz error code.
* \ingroup      WlzHistogram
* \brief	Computes the histogram occupancies of the histogram
*               bins for the given histogram domain and 2D domain
*               object.
*               Because this is a static function there's no need to
*               check the parameters.
* \param	histDom			Histogram domain.
* \param	srcObj			Source 2D domain object.
* \param	pln			Plane for 3D value tables.
*/
static WlzErrorNum WlzHistogramCompute2D(WlzHistogramDomain *histDom,
					 WlzObject *srcObj, int pln)
{
  int		idx,
  		ivCount,
		originI,
		unityBinSize;
  int		*histBin;
  WlzGreyP	objPix;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzHistogramCompute2D FE %p %p %d\n",
	   histDom, srcObj, pln));
  originI = (int )floor(histDom->origin + DBL_EPSILON);
  unityBinSize = (histDom->binSize >= (1.0 - DBL_EPSILON)) &&
                 (histDom->binSize <= (1.0 + DBL_EPSILON));
  histBin = histDom->binValues.inp;
  if((errNum = WlzInitGreyScan(srcObj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
  {
    iWSp.plnpos = pln;
    while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
    {
      ivCount = iWSp.rgtpos - iWSp.lftpos + 1;
      switch(gWSp.pixeltype)
      {
	case WLZ_GREY_INT:
	  objPix.inp = gWSp.u_grintptr.inp;
	  if(unityBinSize)
	  {
	    while(ivCount-- > 0)
	    {
	      idx = *(objPix.inp)++ - originI;
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  else
	  {
	    while(ivCount-- > 0)
	    {
	      idx = (int )floor((*(objPix.inp)++ - histDom->origin) /
	                        histDom->binSize);
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  break;
	case WLZ_GREY_SHORT:
	  objPix.shp = gWSp.u_grintptr.shp;
	  if(unityBinSize)
	  {
	    while(ivCount-- > 0)
	    {
	      idx = *(objPix.shp)++ - originI;
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  else
	  {
	    while(ivCount-- > 0)
	    {
	      idx = (int )floor((*(objPix.shp)++ - histDom->origin) /
	                        histDom->binSize);
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  objPix.ubp = gWSp.u_grintptr.ubp;
	  if(unityBinSize)
	  {
	    if((originI == 0) && (histDom->nBins == 256))
	    {
	      while(ivCount-- > 0)
	      {
		++*(histBin + *(objPix.ubp)++);
	      }
	    }
	    else
	    {
	      idx = *(objPix.ubp)++ - originI;
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  else
	  {
	    while(ivCount-- > 0)
	    {
	      idx = (int )floor((*(objPix.ubp)++ - histDom->origin) /
	                        histDom->binSize);
	      if((idx >= 0) && (idx < histDom->nBins))
	      {
		++*(histBin + idx);
	      }
	    }
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  objPix.flp = gWSp.u_grintptr.flp;
	  while(ivCount-- > 0)
	  {
	    idx = (int )floor((*(objPix.flp)++ - histDom->origin) /
	    		      histDom->binSize);
	    if((idx >= 0) && (idx < histDom->nBins))
	    {
	      ++*(histBin + idx);
	    }
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  objPix.dbp = gWSp.u_grintptr.dbp;
	  while(ivCount-- > 0)
	  {
	    idx = (int )floor((*(objPix.dbp)++ - histDom->origin) / 
	    		      histDom->binSize);
	    if((idx >= 0) && (idx < histDom->nBins))
	    {
	      ++*(histBin + idx);
	    }
	  }
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    (void )WlzEndGreyScan(&iWSp, &gWSp);
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzHistogramCompute2D FX %d\n",
	   (int )errNum));
  return(errNum);
}


/*!
* \return	void
* \ingroup      WlzHistogram
* \brief	Recomputes the given histogram domain for the required
*               number of bins, origin and bin size.
*               This function will only be called for a histogram
*               computed from WlzUByte grey values with the range [0-255]
*               with histDom->maxBins == 256. Also the given histograms
*               binSize is known to be 1.0 and it's origin is known to
*               be 0.0.
*               Because this is a static function there's no need to
*               check the parameters.
* \param	histDom			Histogram domain.
* \param	nBins			Required number of histogram
*                                       bins.
* \param	origin			Required lowest grey value in
*                                       first histogram bin.
* \param	binSize			Required grey value range for
*                                       each histogram bin.
*/
static void	WlzHistogramReBinUbyte(WlzHistogramDomain *histDom,
				       int nBins,
				       double origin, double binSize)
{
  int		srcBinIdx,
  		dstBinIdx;
  int		*histBins;
  double	tD0;
  int		oldBins[256];

  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzHistogramReBinUbyte FE %p %d %g %g\n",
	   histDom, nBins, origin, binSize));
  if((nBins != histDom->nBins) ||
     ((tD0 = (origin - histDom->origin)) > DBL_EPSILON) ||
     (tD0 < DBL_EPSILON) ||
     ((tD0 = (binSize - histDom->binSize)) > DBL_EPSILON) ||
     (tD0 < DBL_EPSILON))
  {
    histBins = histDom->binValues.inp;
    WlzValueCopyIntToInt(oldBins, histBins, 256); /* histDom->maxBins == 256 */
    WlzValueSetInt(histBins, 0, 256);
    if(nBins == 0)
    {
      binSize = 1.0;
      srcBinIdx = 0;
      while((oldBins[srcBinIdx] == 0) && (srcBinIdx < 255))
      {
        ++srcBinIdx;
      }
      origin = srcBinIdx;
      srcBinIdx = 255;
      while((oldBins[srcBinIdx] == 0) && (srcBinIdx > 0))
      {
        --srcBinIdx;
      }
      nBins = srcBinIdx - (int )origin + 1;
    }
    srcBinIdx = 0;
    while(srcBinIdx < 256)
    {
      dstBinIdx = (int )floor((srcBinIdx - origin) / binSize);
      if((dstBinIdx >= 0) && (dstBinIdx < nBins))
      {
        *(histBins + dstBinIdx) += oldBins[srcBinIdx];
      }
      ++srcBinIdx;
    }
    histDom->nBins = nBins;
    histDom->binSize = binSize;
    histDom->origin = origin;

  }
  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzHistogramReBinUbyte FX\n"));
}

/*!
* \return	Woolz error code.
* \ingroup      WlzHistogram
* \brief	Modifies the grey values of the given Woolz 2D domain
*               object so that the resulting object's histogram (given)
*               matches the target histogram.
*               Because this is a static function there's no need to
*               check the parameters.
*               Both the source and target histogram are known to have
*               the same binSize = 1.0, the same number bins and
*               integral type.
* \param	srcObj			Given 2D domain object.
* \param	tgtHD			Target cumulative histogram.
* \param	srcHD			Source object cumulative histogram.
* \param	dither			If non zero then dither mapped
*                                       values to avoid generating a
*                                       histogram with a comb appearence.
*/
static WlzErrorNum WlzHistogramMatchPlane(WlzObject *srcObj,
					  WlzHistogramDomain *tgtHD,
					  WlzHistogramDomain *srcHD,
					  int dither)
{
  int		mapIdx,
  		tgtIdx;
  double	tD0,
		normFac,
  		srcMaxD,
		tgtMaxD;
  int		*mapping;
  WlzObject	*mapHistObj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  

  if(((mapHistObj = WlzMakeHistogram(WLZ_HISTOGRAMDOMAIN_INT,
  				     srcHD->nBins, &errNum)) == NULL) &&
     (errNum == WLZ_ERR_NONE))
  {
   errNum = WLZ_ERR_UNSPECIFIED;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Compute mapping table */
    mapIdx = 0;
    tgtIdx = 0;
    mapHistObj->domain.hist->nBins = srcHD->nBins;
    mapHistObj->domain.hist->origin = srcHD->origin;
    mapHistObj->domain.hist->binSize = srcHD->binSize;
    mapping = mapHistObj->domain.hist->binValues.inp;
    srcMaxD = *(srcHD->binValues.inp + srcHD->nBins - 1);
    tgtMaxD = *(tgtHD->binValues.inp + tgtHD->nBins - 1);
    normFac = srcMaxD / tgtMaxD;
    while(tgtIdx < tgtHD->nBins)
    {
      tD0 = *(tgtHD->binValues.inp + tgtIdx) * normFac;
      while((mapIdx < srcHD->nBins) &&
            (*(srcHD->binValues.inp + mapIdx) < tD0))
      {
	*(mapping + mapIdx++) = tgtIdx;
      }
      ++tgtIdx;
    }
    if(mapIdx < srcHD->nBins)
    {
      /* Need to distribute values better than this! */
      *(mapping + mapIdx) = *(mapping + mapIdx - 1);
    }
    /* Use mapping histogram to map each grey value */
    errNum = WlzHistogramMapValues(srcObj, mapHistObj, dither);
  }
  if(mapHistObj)
  {
    WlzFreeObj(mapHistObj);
  }
  return(errNum);
}

/*!
* \return	Mapped pixel value.
* \ingroup	WlzHistogram
* \brief	Computes an integer pixel value from the given mapping and map
*               index, dithering the mapped value so that it has a
*               uniform probability of falling anywhere between the
*               midpoints of the previous - current and next - current
*               mapping values.
* \param	mapping			Simple integer look up table.
* \param	mapIdx			Index into the look up table.
* \param	nBins			Number of look up table entries.
*/
static int	WlzHistogramMapValuesDitherI(int *mapping, int mapIdx,
					     int nBins)
{
  int		mapP,
  		mapC,
		mapN,
		mappedVal;
  double	mapRnd,
  		mapRngP,
		mapRngN;

  mapC = *(mapping + mapIdx);
  mapP = (mapIdx > 0)? *(mapping + mapIdx - 1): mapC;
  mapN = (mapIdx < nBins - 1)?  *(mapping + mapIdx + 1): mapC;
  if((mapP >= (mapC - 1)) && (mapN <= (mapC + 1)))
  {
    mappedVal = mapC;
  }
  else
  {
    mapRngP = (double )(mapC - mapP) * 0.5;
    mapRngN = (double )(mapN - mapC) * 0.5;
    mapRnd = (AlgRandUniform() * (mapRngP + mapRngN)) - mapRngP;
    mappedVal = mapC + WLZ_NINT(mapRnd);
  }
  return(mappedVal);
}

/*!
* \return	Mapped pixel value.
* \ingroup	WlzHistogram
* \brief	Computes a double pixel value from the given mapping and map
*               index, dithering the mapped value so that it has a
*               uniform probability of falling anywhere between the
*               midpoints of the previous - current and next - current
*               mapping values.
* \param	mapping			Simple integer look up table.
* \param	mapIdx			Index into the look up table.
* \param	nBins			Number of look up table entries.
*/
static double	WlzHistogramMapValuesDitherD(int *mapping, int mapIdx,
					     int nBins)
{
  double	mapP,
  		mapC,
		mapN,
		mapRnd,
  		mapRngP,
		mapRngN,
		mappedVal;

  mapC = *(mapping + mapIdx);
  mapP = (mapIdx > 0)? *(mapping + mapIdx - 1): mapC;
  mapN = (mapIdx < nBins - 1)?  *(mapping + mapIdx + 1): mapC;
  mapRngP = (double )(mapC - mapP) * 0.5;
  mapRngN = (double )(mapN - mapC) * 0.5;
  mapRnd = (AlgRandUniform() * (mapRngP + mapRngN)) - mapRngP;
  mappedVal = mapC + mapRnd;
  return(mappedVal);
}

/*!
* \return	New woolz object with histogram domain.
* \ingroup	WlzHistogram
* \brief	Creates a woolz histogram object (with an int histogram
*               domain) using the given 2D or 3D woolz domain object
*               and histogram parameters.
*               If the requested number of bins is zero then the
*               number of bins, the origin and the bin size will be
*               computed as follows:
*		\verbatim
                  nBins = ceil(max(g) - min(g) + 1.0)
                  binOrigin = min(g)
                  binSize = 1.0
		\endverbatim
*               Where min(g) and max(g) are the minimum and maximum
*               grey values in the source object.
* \param	srcObj			Given source object.
* \param	nBins			Required number of histogram bins.
* \param	binOrigin		Lowest grey value in first histogram
* 					bin.
* \param	binSize			Grey value range for each
*                                       histogram bin.
* \param	dstErrNum		Destination error pointer.
*/
WlzObject	*WlzHistogramObj(WlzObject *srcObj, int nBins,
				 double binOrigin, double binSize,
				 WlzErrorNum *dstErrNum)
{
  int		planeCount,
  		planeIdx,
  		nBins0 = 0,
		isTiled = 0;
  double	binOrigin0 = 0.0,
  		binSize0 = 1.0;
  WlzGreyType	greyType;
  WlzHistogramDomain *histDom,
  		*histDom2D;
  WlzObject	*histObj = NULL,
  		*histObj2D = NULL,
		*srcObj2D = NULL;
  WlzDomain	dummyDom;
  WlzValues	dummyValues;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzPixelV	greyMinV,
  		greyMaxV;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramObj FE %p %d %g %g %p\n",
	   srcObj, nBins, binOrigin, binSize, dstErrNum));
  dummyDom.core = NULL;
  dummyValues.core = NULL;
  if(nBins < 0)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	if(((histObj = WlzMakeEmpty(&errNum)) == NULL) &&
	   (errNum == WLZ_ERR_NONE))
	{
	  errNum = WLZ_ERR_UNSPECIFIED;
	}
        break;
      case WLZ_2D_DOMAINOBJ:
	if((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						      srcObj)) == WLZ_ERR_NONE)
	{
	  switch(greyType)
	  {
	    case WLZ_GREY_UBYTE:
	      nBins0 = 256;    /* Compute range after histogram if WlzUByte */
	      binOrigin0 = 0.0;
	      binSize0 = 1.0;
	      break;
	    case WLZ_GREY_INT:
	    case WLZ_GREY_SHORT:
	    case WLZ_GREY_FLOAT:
	    case WLZ_GREY_DOUBLE:
	      if(nBins == 0)
	      {
	        if((errNum = WlzGreyRange(srcObj, &greyMinV,
					  &greyMaxV)) == WLZ_ERR_NONE)
		{
		  (void )WlzValueConvertPixel(&greyMinV, greyMinV,
					      WLZ_GREY_DOUBLE);
		  (void )WlzValueConvertPixel(&greyMaxV, greyMaxV,
					      WLZ_GREY_DOUBLE);
		  nBins = (int )ceil(greyMaxV.v.dbv - greyMinV.v.dbv + 1.0);
		  binOrigin = greyMinV.v.dbv;
		  nBins0 = nBins;
		  binOrigin0 = binOrigin;
		  binSize0 = binSize;
		}
	      }
	      else
	      {
	        nBins0 = nBins;
		binOrigin0 = binOrigin;
		binSize0 = binSize;
	      }
	      break;
	    case WLZ_GREY_RGBA: /* RGBA to be done RAB */
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(((histObj = WlzMakeHistogram(WLZ_HISTOGRAMDOMAIN_INT,
	  				  nBins0, &errNum)) == NULL) &&
	      (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WLZ_ERR_UNSPECIFIED;
	  }
	  else
	  {
	    histDom = histObj->domain.hist;
	    histDom->nBins = nBins0;
	    histDom->origin = binOrigin0;
	    histDom->binSize = binSize0;
	    errNum = WlzHistogramCompute2D(histDom, srcObj, 0);
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (greyType == WLZ_GREY_UBYTE) &&
	   ((nBins != nBins0) ||
	    (fabs(binOrigin - binOrigin0) > DBL_EPSILON) ||
	    (fabs(binSize - binSize0) > DBL_EPSILON)))
	{
	  WlzHistogramReBinUbyte(histObj->domain.hist, nBins, binOrigin,
	  			 binSize);
	}
        break;
      case WLZ_TRANS_OBJ:
	if((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						      srcObj)) == WLZ_ERR_NONE)
	{
	  histObj = WlzHistogramObj(srcObj->values.obj, nBins,
	  			    binOrigin, binSize, &errNum);
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						      srcObj)) == WLZ_ERR_NONE)
	{
	  switch(greyType)
	  {
	    case WLZ_GREY_UBYTE:
	      nBins0 = 256;
	      binOrigin0 = 0.0;
	      binSize0 = 1.0;
	      break;
	    case WLZ_GREY_INT:
	    case WLZ_GREY_SHORT:
	    case WLZ_GREY_FLOAT:
	    case WLZ_GREY_DOUBLE:
	      if(nBins == 0)
	      {
	        if((errNum = WlzGreyRange(srcObj, &greyMinV,
					  &greyMaxV)) == WLZ_ERR_NONE)
		{
		  (void )WlzValueConvertPixel(&greyMinV, greyMinV,
					      WLZ_GREY_DOUBLE);
		  (void )WlzValueConvertPixel(&greyMaxV, greyMaxV,
					      WLZ_GREY_DOUBLE);
		  nBins = (int )ceil(greyMaxV.v.dbv - greyMinV.v.dbv + 1.0);
		  binOrigin = greyMinV.v.dbv;
		}
	      }
	      else
	      {
	        nBins0 = nBins;
		binOrigin0 = binOrigin;
		binSize0 = binSize;
	      }
	      break;
	    case WLZ_GREY_RGBA: /* RGBA to be done RAB */
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  histObj = WlzMakeHistogram(WLZ_HISTOGRAMDOMAIN_INT, nBins0,
	  			     &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	     histObj2D = WlzMakeHistogram(WLZ_HISTOGRAMDOMAIN_INT, nBins0,
	     				  &errNum);
	  }
	  if(((histObj == NULL) || (histObj2D == NULL)) &&
	     (errNum == WLZ_ERR_NONE))
	  {
	    errNum = WLZ_ERR_UNSPECIFIED;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    histDom = histObj->domain.hist;
	    histDom->nBins = nBins0;
	    histDom->origin = binOrigin0;
	    histDom->binSize = binSize0;
	    histDom2D = histObj2D->domain.hist;
	    histDom2D->nBins = nBins0;
	    histDom2D->origin = binOrigin0;
	    histDom2D->binSize = binSize0;
	    planeCount = srcObj->domain.p->lastpl -
	    		 srcObj->domain.p->plane1 + 1;
	    srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom, dummyValues,
	    			   NULL, NULL, &errNum);
	    if((srcObj2D == NULL) && (errNum == WLZ_ERR_NONE))
	    {
	      errNum = WLZ_ERR_UNSPECIFIED;
	    }
	    else
	    {
	      planeIdx = 0;
	      WlzValueSetInt(histDom->binValues.inp, 0, histDom->maxBins);
	      while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	      {
		srcObj2D->domain = *(srcObj->domain.p->domains + planeIdx);
		if(isTiled)
		{
		  srcObj2D->values = srcObj->values;
		}
		else
		{
		  srcObj2D->values = *(srcObj->values.vox->values + planeIdx);
		}
		if(srcObj2D->domain.core && srcObj2D->values.core &&
		   (srcObj2D->domain.core->type != WLZ_EMPTY_OBJ) &&
		   (srcObj2D->values.core->type != WLZ_EMPTY_OBJ))
		{
		  WlzValueSetInt(histDom2D->binValues.inp, 0,
				 histDom2D->maxBins);
		  if((errNum = WlzHistogramCompute2D(histDom2D, srcObj2D,
		                                    planeIdx)) == WLZ_ERR_NONE)
		  {
		    WlzHistogramAddIntBins(histDom, histDom2D);
		  }
	        }
		++planeIdx;
	      }
	      srcObj2D->domain = dummyDom;
	      srcObj2D->values = dummyValues;
	      WlzFreeObj(srcObj2D);
	    }
	    WlzFreeObj(histObj2D);
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (greyType == WLZ_GREY_UBYTE) &&
	   ((nBins != nBins0) ||
	    (fabs(binOrigin - binOrigin0) > DBL_EPSILON) ||
	    (fabs(binSize - binSize0) > DBL_EPSILON)))
	{
	  WlzHistogramReBinUbyte(histObj->domain.hist, nBins, binOrigin,
	  			 binSize);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(histObj)
    {
      WlzFreeObj(histObj);
      histObj = NULL;
    }
    if(dstErrNum)
    {
      *dstErrNum = errNum;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramObj 01 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramObj FX %p\n",
	   histObj));
  return(histObj);
}

/*!
* \return	Copied Woolz histogram object or NULL on error.
* \ingroup	WlzHistogram
* \brief	Copies the given Woolz histogram object to a new one of
*               the required type.
* \param	srcHistObj		Given histogram object.
* \param	dstType			Required type of histogram domain
* 					values.
* \param	dstErrNum		Destination error pointer,
*                                       may be NULL.
*/
WlzObject	*WlzHistogramCopy(WlzObject *srcHistObj,
				  WlzObjectType dstType,
				  WlzErrorNum *dstErrNum)
{
  WlzHistogramDomain *srcHistDom,
		*dstHistDom;
  WlzObject	*dstHistObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramCopy FE %p %d %p\n",
	   srcHistObj, (int )dstType, dstErrNum));
  if((errNum = WlzHistogramCheckHistObj(srcHistObj)) == WLZ_ERR_NONE)
  {
    srcHistDom = srcHistObj->domain.hist;
    dstHistObj = WlzMakeHistogram(dstType, srcHistDom->maxBins, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstHistDom = dstHistObj->domain.hist;
    dstHistDom->nBins = srcHistDom->nBins;
    dstHistDom->origin = srcHistDom->origin;
    dstHistDom->binSize = srcHistDom->binSize;
    WlzValueCopyGreyToGrey(dstHistDom->binValues, 0,
    			   (dstHistDom->type == WLZ_HISTOGRAMDOMAIN_INT)?
			   (WLZ_GREY_INT): (WLZ_GREY_DOUBLE),
			   srcHistDom->binValues, 0,
			   (srcHistDom->type == WLZ_HISTOGRAMDOMAIN_INT)?
			   (WLZ_GREY_INT): (WLZ_GREY_DOUBLE),
			   srcHistDom->nBins);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramCopy 01 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramCopy FX %p\n",
	   dstHistObj));
  return(dstHistObj);
}

/*!
* \return	Copied Woolz histogram object or NULL on error.
* \ingroup	WlzHistogram
* \brief	Re-bins the given Woolz histogram object to a new one
*               of the required type, number of bins, origin and bin
*               size.
* \param	srcHistObj		Given histogram object.
* \param	dstType			Required type of histogram
*                                       domain values.
* \param	dstMaxBins		Total number of bins to be
*                                       allocated.
* \param	dstNBins		Number of bins to be used in
*                                       new histogram.
* \param	dstOrigin		New histogram origin.
* \param	dstBinSize		New histogram bin size.
* \param	dstErrNum		Destination error pointer,
*                                       may be NULL.
*/
WlzObject	*WlzHistogramRebin(WlzObject *srcHistObj,
				   WlzObjectType dstType,
				   int dstMaxBins,
				   int dstNBins,
				   double dstOrigin,
				   double dstBinSize,
				   WlzErrorNum *dstErrNum)
{
  int		srcIdx,
		dstCount;
  double	srcPos,
  		srcStep;
  int		*dstValIP;
  double	*dstValDP;
  WlzHistogramDomain *srcHistDom,
		*dstHistDom;
  WlzObject	*dstHistObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramRebin FE %p %d %d %d %g %g %p\n",
	   srcHistObj, (int )dstType, dstMaxBins, dstNBins,
	   dstOrigin, dstBinSize, dstErrNum));
  if((errNum = WlzHistogramCheckHistObj(srcHistObj)) == WLZ_ERR_NONE)
  {
    srcHistDom = srcHistObj->domain.hist;
    dstHistObj = WlzMakeHistogram(dstType, dstMaxBins, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstHistDom = dstHistObj->domain.hist;
    dstHistDom->maxBins = dstMaxBins;
    dstHistDom->nBins = dstNBins;
    dstHistDom->origin = dstOrigin;
    dstHistDom->binSize = dstBinSize;
    srcPos = (dstOrigin - srcHistDom->origin) / srcHistDom->binSize;
    srcStep = dstBinSize / srcHistDom->binSize;
    dstCount = dstNBins;
    if(dstType == WLZ_HISTOGRAMDOMAIN_INT)
    {
      dstValIP = dstHistDom->binValues.inp;
      if(srcHistDom->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
	while(dstCount-- > 0)
	{
	  srcIdx = WLZ_NINT(srcPos);
	  *dstValIP++ = ((srcIdx < 0) || (srcIdx >= srcHistDom->nBins))?
	  		0: *(srcHistDom->binValues.inp + srcIdx);
	  srcPos += srcStep;
	}
      }
      else
      {
	while(dstCount-- > 0)
	{
	  srcIdx = WLZ_NINT(srcPos);
	  *dstValIP++ = ((srcIdx < 0) || (srcIdx >= srcHistDom->nBins))?
	  		0:  WLZ_NINT(*(srcHistDom->binValues.dbp + srcIdx));
	  srcPos += srcStep;
        }
      }
    }
    else
    {
      dstValDP = dstHistDom->binValues.dbp;
      if(srcHistDom->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
	while(dstCount-- > 0)
	{
	  srcIdx = WLZ_NINT(srcPos);
	  *dstValDP++ = ((srcIdx < 0) || (srcIdx >= srcHistDom->nBins))?
	                0.0: *(srcHistDom->binValues.inp + srcIdx);
	  srcPos += srcStep;
	}
      }
      else
      {
	while(dstCount-- > 0)
	{
	  srcIdx = WLZ_NINT(srcPos);
	  *dstValDP++ = ((srcIdx < 0) || (srcIdx >= srcHistDom->nBins))?
	  		0.0: *(srcHistDom->binValues.dbp + srcIdx);
	  srcPos += srcStep;
        }
      }
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramRebin 01 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramRebin FX %p\n",
	   dstHistObj));
  return(dstHistObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Applies an infinite impulse response (IIR) filter to the
*		histogram object'c bin values.
* \param	histObj			Given histogram object.
* \param	flt			IIR filter.
*/
WlzErrorNum     WlzHistogramRsvFilter(WlzObject *histObj, WlzRsvFilter *flt)
{
  double	*datBuf = NULL,
  		*tmpBuf0 = NULL,
  		*tmpBuf1 = NULL;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramRsvFilter FE %p %p\n",
	   histObj, flt));
  if(flt == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzHistogramCheckHistObj(histObj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    histDom = histObj->domain.hist;
    if(((datBuf  = (double *)AlcMalloc(sizeof(double) * 
    				       histDom->nBins)) == NULL) ||
       ((tmpBuf0  = (double *)AlcMalloc(sizeof(double) * 
    				        histDom->nBins)) == NULL) ||
       ((tmpBuf1  = (double *)AlcMalloc(sizeof(double) * 
    				        histDom->nBins)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      WlzValueCopyIntToDouble(datBuf, histDom->binValues.inp,
			      histDom->nBins);
    }
    else
    {
      WlzValueCopyDoubleToDouble(datBuf, histDom->binValues.dbp,
				 histDom->nBins);
    }
    errNum = WlzRsvFilterBuffer(flt, datBuf, tmpBuf0, tmpBuf1,
				histDom->nBins);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      WlzValueCopyDoubleToInt(histDom->binValues.inp, datBuf,
      			      histDom->nBins);
    }
    else
    {
      WlzValueCopyDoubleToDouble(histDom->binValues.dbp, datBuf,
      				 histDom->nBins);
    }
  }
  /* Free up the buffers. */
  if(datBuf)
  {
    AlcFree(datBuf);
  }
  if(tmpBuf0)
  {
    AlcFree(tmpBuf0);
  }
  if(tmpBuf1)
  {
    AlcFree(tmpBuf1);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramRsvFilter FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Convolves the given Woolz histogram object's histogram
*               bin values with the given convolution kernel.
* \param	histObj			Given histogram object.
* \param	krnSz			Kernel size, must be odd.
* \param	krn			The convolution kernel.
*/
WlzErrorNum     WlzHistogramConvolve(WlzObject *histObj,
				     int krnSz, double *krn)
{
  double	*cnvBuf = NULL,
  		*datBuf = NULL;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  AlgError 	algErr = ALG_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramConvolve FE %p %d %p\n",
	   histObj, krnSz, krn));
  if(((errNum = WlzHistogramCheckHistObj(histObj)) == WLZ_ERR_NONE) &&
     ((histDom = histObj->domain.hist)->nBins > 0))
  {
    if((krnSz < 1) || ((krnSz % 2) != 1))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    /* Allocate buffer(s). */
    if(errNum == WLZ_ERR_NONE)
    {
      if((cnvBuf = (double *)AlcMalloc(sizeof(double) *
				       histDom->nBins)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else if(histDom->type ==  WLZ_HISTOGRAMDOMAIN_INT)
      {
	if((datBuf = (double *)AlcMalloc(sizeof(double) *
	 				 histDom->nBins)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Fill the histogram bin data buffer if int histogram bins and convolve
       * histogram bin values with the kernel. */
      if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
	WlzValueCopyIntToDouble(datBuf, histDom->binValues.inp,
				histDom->nBins);
        algErr = AlgConvolve(histDom->nBins, cnvBuf, 
			     krnSz, krn,
			     histDom->nBins, datBuf,
			     ALG_PAD_END);
        if(algErr == ALG_ERR_NONE)
	{
	  WlzValueCopyDoubleToInt(histDom->binValues.inp, cnvBuf,
			          histDom->nBins);
	}
	else
	{
	  errNum = WlzErrorFromAlg(algErr);
	}
      }
      else
      {
        algErr = AlgConvolve(histDom->nBins, cnvBuf,
			     histDom->nBins, histDom->binValues.dbp,
			     krnSz, krn, ALG_PAD_END);
        if(algErr == ALG_ERR_NONE)
	{
	  WlzValueCopyDoubleToDouble(histDom->binValues.dbp, cnvBuf,
	  			 histDom->nBins);
	}
	else
	{
	  errNum = WlzErrorFromAlg(algErr);
	}
      }
    }
    /* Free up the buffers. */
    if(cnvBuf)
    {
      AlcFree(cnvBuf);
    }
    if(datBuf)
    {
      AlcFree(datBuf);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramConvolve FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Convolves the given Woolz histogram object's histogram
*               bin values with a Gaussian, 1st derivative of Gaussian
*               or 2nd derivative of Gaussian convolution kernel.
* \param	histObj			Given histogram object.
* \param	sigma			Sigma value of the Gaussian
*                                       which is scaled within this
*                                       function to account for the
*                                       histogram bin size.
* \param	deriv			Derivative must be in range [0-2].
*/
WlzErrorNum	WlzHistogramCnvGauss(WlzObject *histObj,
				     double sigma, int deriv)
{
  int		krnCnt,
  		krnSz,
  		krnSz2;
  double	tD0,
  		tD1,
		tD2,
  		sigma2,
		norm;
  double	*kP0,
  		*kP1,
		*krnBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramCnvGauss FE %p %g %d\n",
	   histObj, sigma, deriv));
  if((errNum = WlzHistogramCheckHistObj(histObj)) == WLZ_ERR_NONE)
  {
    if((sigma < DBL_EPSILON) || (deriv < 0) || (deriv > 2))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Compute the required kernel size and allocate the kernel. */
  if(errNum == WLZ_ERR_NONE)
  {
    /* Scale the given sigma value to reflect the histogram bin size. */
    sigma /= histObj->domain.hist->binSize;
    if(sigma < (1.0 + DBL_EPSILON))
    {
      krnSz2 = 1;
    }
    else
    {
      krnSz2 = (int )ceil(sigma) * (1 + deriv);
    }
    krnSz = (krnSz2 * 2) + 1;
    if((krnBuf = (double *)AlcMalloc(sizeof(double) * krnSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    /* Don't worry about constants because renormalise kernel anyway. */
    norm = 0.0;
    sigma2 = 1.0 / (2.0 * sigma * sigma);
    switch(deriv)
    {
      case 0:
	kP0 = krnBuf;
        kP1 = krnBuf + krnSz - 1;
	krnCnt = krnSz2;
	while(krnCnt > 0)
	{
	  norm += *kP1-- = *kP0++ = exp(-(krnCnt * krnCnt * sigma2));
	  --krnCnt;
	}
	*(krnBuf + krnSz2) = 1.0;
	norm = 1.0 / ((norm * 2.0) + 1.0);
	break;
      case 1:
	kP0 = krnBuf;
	kP1 = krnBuf + krnSz - 1;
	krnCnt = krnSz2;
	while(krnCnt > 0)
	{
	  tD1 = exp(-(krnCnt * krnCnt * sigma2));
	  norm += *kP0++ =  krnCnt * tD1;
	  *kP1-- =  -krnCnt * tD1;
	  --krnCnt;
	}
	*(krnBuf + krnSz2) = 0.0;
	norm = 1.0 / (norm * 2.0);
        break;
      case 2:
	kP0 = krnBuf;
	kP1 = krnBuf + krnSz - 1;
	krnCnt = krnSz2;
	while(krnCnt > 0)
	{
	  tD0 = krnCnt * krnCnt;
	  tD1 = exp(-(tD0 * sigma2));
	  tD2 = ((tD0 * sigma2 * 2.0) - 1.0) *tD1;
	  norm += *kP1-- = *kP0++ =  tD2;
	  --krnCnt;
	}
	*(krnBuf + krnSz2) = -1.0;
	norm = 1.0 / ((norm * 2.0) - 1.0);
        break;
    }
    kP0 = krnBuf;
    krnCnt = krnSz;
    while(krnCnt-- > 0)
    {
      *kP0++ *= norm;
    }
/* #define WLZ_HISTOGRAM_DEBUG_CNVGAUSS */
#ifdef WLZ_HISTOGRAM_DEBUG_CNVGAUSS
    for(krnCnt = 0; krnCnt < krnSz; ++krnCnt)
    {
      printf("%d %g\n", krnCnt - krnSz2, *(krnBuf + krnCnt));
    }
    exit(0);
#endif /* WLZ_HISTOGRAM_DEBUG_CNVGAUSS */
  }
  /* Set the kernel values and do the convolution. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzHistogramConvolve(histObj, krnSz, krnBuf);
  }
  /* Free up the kernel. */
  if(krnBuf)
  {
    AlcFree(krnBuf);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramCnvGauss FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Uses an infinite impulse response (IIR) filter to
*		convolve the given Woolz histogram object's bin values with
*		a Gaussian, 1st derivative of Gaussian or 2nd derivative of
*		Gaussian.
* \param	histObj			Given histogram object.
* \param	sigma			Sigma value of the Gaussian
*                                       which is scaled within this
*                                       function to account for the
*                                       histogram bin size.
* \param	deriv			Derivative must be in range
*                                       [0-2].
*/
WlzErrorNum	WlzHistogramRsvGauss(WlzObject *histObj, double sigma,
				     int deriv)
{
  WlzRsvFilter *ftr = NULL;
  WlzRsvFilterName ftrName;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramRsvGauss FE %p %g %d\n",
	   histObj, sigma, deriv));
  if((errNum = WlzHistogramCheckHistObj(histObj)) == WLZ_ERR_NONE)
  {
    if(sigma < DBL_EPSILON)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else
    {
      switch(deriv)
      {
        case 0:
	  ftrName = WLZ_RSVFILTER_NAME_GAUSS_0;
	  break;
	case 1:
	  ftrName = WLZ_RSVFILTER_NAME_GAUSS_1;
	  break;
	case 2:
	  ftrName = WLZ_RSVFILTER_NAME_GAUSS_2;
	  break;
	default:
	  errNum = WLZ_ERR_PARAM_DATA;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ftr = WlzRsvFilterMakeFilter(ftrName, sigma, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzHistogramRsvFilter(histObj, ftr);
  }
  /* Free up the filter. */
  if(ftr)
  {
    WlzRsvFilterFreeFilter(ftr);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramRsvGauss FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Low pass filters the given Woolz histogram object's
*               histogram bin values by convolving them in the space
*               domain with a Gaussian kernel of the given half height
*               full width.
* \param	histObj			Given histogram object.
* \param	width			Gaussian kernel half height
*                                       full width.
*/
WlzErrorNum	WlzHistogramSmooth(WlzObject *histObj, int width)
{
  const double sigmaFromWidth = 0.424660900; 	      /* 1/(2 * sqrt(ln(4))) */
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramSmooth FE %p %d\n",
	   histObj, width));
  if(width > 0)
  {
    errNum = WlzHistogramRsvGauss(histObj, sigmaFromWidth * width, 0);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramSmooth FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Modifies the given histogram so that its values are
*               cummulative.
* \param	srcHist			Given histogram object.
*/
WlzErrorNum	WlzHistogramCummulative(WlzObject *srcHist)
{
  int		count;
  int		*tIP0,
  		*tIP1;
  double	*tDP0,
  		*tDP1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramCummulative FE %p\n",
	   srcHist));
  if((errNum = WlzHistogramCheckHistObj(srcHist)) == WLZ_ERR_NONE)
  {
    if((count = srcHist->domain.hist->nBins) > 1)
    {
      if(srcHist->domain.hist->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
	tIP0 = srcHist->domain.hist->binValues.inp;
	tIP1 = tIP0 + 1;
	while(--count > 0)
	{
	  *tIP1++ += *tIP0++;
	}
      }
      else
      {
	tDP0 = srcHist->domain.hist->binValues.dbp;
	tDP1 = tDP0 + 1;
	while(--count > 0)
	{
	  *tDP1++ += *tDP0++;
	}
      }
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramCummulative FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Uses the given mapping histogram with integral bin
*               values to remap the grey values of the given 2D or
*               3D domain object.
*               The mapping histogram MUST have integral bin values
*               and bins appropriate for all domain object values.
* \param	srcObj			Given 2D or 3D domain object.
* \param	mapHistObj		Mapping histogram.
* \param	dither			If non zero then dither mapped
*                                       values to avoid generating a
*                                       histogram with a comb appearence.
*/
WlzErrorNum	WlzHistogramMapValues(WlzObject *srcObj,
				      WlzObject *mapHistObj,
				      int dither)
{
  int		tI0,
		ivCount,
		originI,
  		planeIdx,
		planeCount,
		isTiled = 0;
  double	originD;
  int		*mapping;
  WlzGreyType	greyType;
  WlzHistogramDomain *mapHistDom = NULL;
  WlzGreyP	objPix;
  WlzObject	*srcObj2D = NULL;
  WlzDomain	dummyDom;
  WlzValues	dummyValues;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramMapValues FE %p %p %d\n",
	   srcObj, mapHistObj, dither));
  if(((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						 srcObj)) == WLZ_ERR_NONE) &&
     ((errNum = WlzHistogramCheckHistObj(mapHistObj)) == WLZ_ERR_NONE))
  {
    mapHistDom = mapHistObj->domain.hist;
    if((mapHistDom->type != WLZ_HISTOGRAMDOMAIN_INT) ||
       (mapHistDom->binSize <= (1.0 - DBL_EPSILON)) ||
       (mapHistDom->binSize >= (1.0 + DBL_EPSILON)))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if((errNum = WLZ_ERR_NONE) && isTiled)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        break;
      case WLZ_TRANS_OBJ:
	errNum = WlzHistogramMapValues(srcObj->values.obj, mapHistObj, dither);
        break;
      case WLZ_3D_DOMAINOBJ:
	dummyDom.core = NULL;
	dummyValues.core = NULL;
	if(((srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
				    dummyValues, NULL, NULL,
				    &errNum)) == NULL) &&
	   (errNum == WLZ_ERR_NONE))
	{
	  errNum = WLZ_ERR_UNSPECIFIED;
	}
	else
	{
	  planeIdx =  0;
	  planeCount = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
	  while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	  {
	    srcObj2D->domain = *(srcObj->domain.p->domains + planeIdx);
	    srcObj2D->values = *(srcObj->values.vox->values + planeIdx);
	    if((srcObj2D->values.core) &&
	       (srcObj2D->values.core->type != WLZ_EMPTY_OBJ))
	    {
		errNum = WlzHistogramMapValues(srcObj2D, mapHistObj, dither);
	    }
	    ++planeIdx;
	  }
	  srcObj2D->domain.core = NULL;
	  srcObj2D->values.core = NULL;
	  WlzFreeObj(srcObj2D);
	}
        break;
      case WLZ_2D_DOMAINOBJ:
	if((errNum = WlzInitGreyScan(srcObj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
	{
          mapping = mapHistDom->binValues.inp;
	  originD = mapHistDom->origin;
	  originI = (int )floor(originD);
	  while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
	  {
	    ivCount = iWSp.rgtpos - iWSp.lftpos + 1;
	    switch(gWSp.pixeltype)
	    {
	      case WLZ_GREY_INT:
		objPix.inp = gWSp.u_grintptr.inp;
		if(dither == 0)
		{
		  if(originI)
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.inp) = *(mapping + *(objPix.inp) - originI);
		      ++(objPix.inp);
		    }
		  }
		  else
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.inp) = *(mapping + *(objPix.inp));
		      ++(objPix.inp);
		    }
		  }
		}
		else
		{
		  while(ivCount-- > 0)
		  {
		    *(objPix.inp) = WlzHistogramMapValuesDitherI(mapping,
		    			*(objPix.inp) - originI,
					mapHistDom->nBins);
		    ++(objPix.inp);
		  }

		}
		break;
	      case WLZ_GREY_SHORT:
		objPix.shp = gWSp.u_grintptr.shp;
		if(dither == 0)
		{
		  if(originI)
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.shp) = (short )
		                      *(mapping + *(objPix.shp) - originI);
		      ++(objPix.shp);
		    }
		  }
		  else
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.shp) = (short )
		                      *(mapping + *(objPix.shp));
		      ++(objPix.shp);
		    }
		  }
		}
		else
		{
		  while(ivCount-- > 0)
		  {
		    *(objPix.shp) = (short )
		                    WlzHistogramMapValuesDitherI(mapping,
		    			*(objPix.shp) - originI,
					mapHistDom->nBins);
		    ++(objPix.shp);
		  }
		}
		break;
	      case WLZ_GREY_UBYTE:
		objPix.ubp = gWSp.u_grintptr.ubp;
		if(dither == 0)
		{
		  if(originI)
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.ubp) = (WlzUByte )
		                      *(mapping + *(objPix.ubp) - originI);
		      ++(objPix.ubp);
		    }
		  }
		  else
		  {
		    while(ivCount-- > 0)
		    {
		      *(objPix.ubp) = (WlzUByte )*(mapping + *(objPix.ubp));
		      ++(objPix.ubp);
		    }
		  }
		}
		else
		{
		  while(ivCount-- > 0)
		  {
		    *(objPix.ubp) = (WlzUByte )
		                    WlzHistogramMapValuesDitherI(mapping,
		    			*(objPix.ubp) - originI,
					mapHistDom->nBins);
		    ++(objPix.ubp);
		  }
		}
		break;
	      case WLZ_GREY_FLOAT:
		objPix.flp = gWSp.u_grintptr.flp;
		if(dither == 0)
		{
		  while(ivCount-- > 0)
		  {
		    tI0 = (int )floor(*(objPix.flp) - originD);
		    *(objPix.flp) = (float )*(mapping + tI0);
		    ++(objPix.flp);
		  }
		}
		else
		{
		  while(ivCount-- > 0)
		  {
		    *(objPix.flp) = (float )
		    		    WlzHistogramMapValuesDitherD(mapping,
		    			(int )(floor(*(objPix.flp) - originD)),
					mapHistDom->nBins);
		    ++(objPix.flp);
		  }
		}
		break;
	      case WLZ_GREY_DOUBLE:
		objPix.dbp = gWSp.u_grintptr.dbp;
		if(dither == 0)
		{
		  while(ivCount-- > 0)
		  {
		    tI0 = (int )floor(*(objPix.dbp) - originD);
		    *(objPix.dbp) = *(mapping + tI0);
		    ++(objPix.dbp);
		  }
		}
		else
		{
		  while(ivCount-- > 0)
		  {
		    *(objPix.dbp) = WlzHistogramMapValuesDitherD(mapping,
		    			(int )(floor(*(objPix.flp) - originD)),
					mapHistDom->nBins);
		    ++(objPix.dbp);
		  }
		}
		break;
	      default:
	        errNum = WLZ_ERR_GREY_TYPE;
	        break;
	    }
	  }
          (void )WlzEndGreyScan(&iWSp, &gWSp);
	  if(errNum == WLZ_ERR_EOO)	/* Reset error from end of intervals */
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramMapValues FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Find the peaks and/or troughs  of the given histogram.
*               Peaks are found by searching through the 1st derivative
*               (h') of the histogram bin values (h). For a peak
*		\verbatim
                  (h'[i - 1] > 0) && (h'[i + 1] < 0) &&
                  (h[i] > thresh).
		\endverbatim
*               If the first or last histogram bin has the maximum
*               bin occupancy of the histogram it is classified as a
*               peak too.
* \param	histObj			Given histogram object.
* \param	sigma			Gaussian sigma smooth histogram
*                                       and compute 1st derivative.
* \param	thresh			Minimum (smoothed) histogram
*                                       value at peak, not used for
*                                       troughs.
* \param	dstFeatSz		Destination ptr for the number
*                                       of peaks/troughs found.
* \param	dstFeat			Destination ptr for the array of
*                                       peak/trough positions (indicies
*                                       into the histogram bin values).
*                                       May have *dstFeat == NULL on
*                                       return if *dstFeatSz == 0.
* \param	feat			Features to find: peaks, troughs
*                                       or both.
*/
WlzErrorNum     WlzHistogramFindPeaks(WlzObject *histObj,
				      double sigma, double thresh,
				      int *dstFeatSz, int **dstFeat,
				      WlzHistFeature feat)
{
  int		binIdx,
		binCnt1,
		binLst,
		fndFlg,
  		featIdx,
  		featMax,
		threshI,
		maxI,
		maxBinIdx,
		minSep;
  double	maxD;
  int		*dIC,
  		*dIN,
		*sIC,
		*featDat = NULL;
  double	*dDC,
  		*dDN,
		*sDC;
  WlzObject 	*sHist = NULL,
  		*dHist = NULL;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if((errNum = WlzHistogramCheckHistObj(histObj)) == WLZ_ERR_NONE)
  {
    if((dstFeatSz == NULL) || (dstFeat == NULL) ||
       ((feat & ~(WLZ_HIST_FEATURE_PEAK | WLZ_HIST_FEATURE_TROUGH)) != 0))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if((histDom = histObj->domain.hist)->nBins > 2)
  {
    /* Copy the given histogram and then compute the smoothed and
     * derivative histograms using these copies. */
    if(errNum == WLZ_ERR_NONE)
    {
      sHist = WlzCopyObject(histObj, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      dHist = WlzCopyObject(histObj, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Smooth. */
      errNum = WlzHistogramCnvGauss(sHist, sigma, 0);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* 1st derivative. */
      errNum = WlzHistogramCnvGauss(dHist, sigma, 1);
    }
    /* Allocate a feature array. */
    if(errNum == WLZ_ERR_NONE)
    {
      featMax = histDom->nBins / 2;
      if((featDat = (int *)AlcMalloc(sizeof(int) * featMax)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    /* Search through the histogram bins for peak/trough features. */
    if(errNum == WLZ_ERR_NONE)
    {
      minSep = (int )ceil(sigma + 0.1);
      binLst = -(minSep + 1);
      featIdx = 0;
      binIdx = 0;
      binCnt1 = histDom->nBins - 1;
      maxBinIdx = WlzHistogramBinMax(sHist->domain.hist);
      switch(histDom->type)
      {
        case WLZ_HISTOGRAMDOMAIN_INT:
	  threshI = (int )ceil(thresh);
	  dIC = dHist->domain.hist->binValues.inp;
	  dIN = dIC + 1;
	  sIC = sHist->domain.hist->binValues.inp;
	  maxI = *(sIC + maxBinIdx);
	  while(binIdx < binCnt1)
	  {
	    fndFlg = 0;
	    if(feat & WLZ_HIST_FEATURE_PEAK)
	    {
	      if(*sIC > threshI)
	      {
		if(*sIC == maxI)
		{
		  if(((binIdx - binLst) > minSep) || (featIdx == 0))
		  {
		    fndFlg = 1;
		    *(featDat + featIdx++) = binIdx;
		    binLst = binIdx;
		  }
		  else
		  {
		    fndFlg = 1;
		    *(featDat + featIdx - 1) = binIdx;
		    binLst = binIdx;
		  }
		}
		else if(((binIdx - binLst) > minSep) &&
			(*dIC <= 0) && (*dIN > 0))
		{
		  fndFlg = 1;
		  *(featDat + featIdx++) = binIdx;
		  binLst = binIdx;
		}
	      }
	    }
	    if((fndFlg == 0) && (feat & WLZ_HIST_FEATURE_TROUGH))
	    {
	      if(((binIdx - binLst) > minSep) &&
	      	  (*dIC > 0) && (*dIN <= 0))
	      {
	        fndFlg = 1;
		*(featDat + featIdx++) = binIdx;
		binLst = binIdx;
	      }
	    }
	    ++dIC; ++dIN; ++sIC;
	    ++binIdx;
	  }
	  break;
        case WLZ_HISTOGRAMDOMAIN_FLOAT:
	  dDC = dHist->domain.hist->binValues.dbp;
	  dDN = dDC + 1;
	  sDC = sHist->domain.hist->binValues.dbp;
	  maxD = *(sDC + maxBinIdx);
	  while(binIdx < binCnt1)
	  {
	    if(*sDC > thresh)
	    {
	      if(fabs(*sDC - maxD) < DBL_EPSILON)
	      {
		if(((binIdx - binLst) > minSep) || (featIdx == 0))
		{
		  *(featDat + featIdx++) = binIdx;
		  binLst = binIdx;
		}
		else
		{
		  *(featDat + featIdx - 1) = binIdx;
		  binLst = binIdx;
		}
	      }
	      else if(((binIdx - binLst) > minSep) &&
		      (*dDC <= -(DBL_EPSILON)) && (*dDN > DBL_EPSILON))
	      {
		*(featDat + featIdx++) = binIdx;
		binLst = binIdx;
	      }
	    }
	    ++dDC; ++dDN; ++sDC;
	    ++binIdx;
	  }
	  break;
        default:
	  break;
      }
      if(featIdx == 0)
      {
        AlcFree(featDat);
	featDat = NULL;
      }
      *dstFeatSz = featIdx;
      *dstFeat = featDat;
    }
    if(sHist)
    {
      WlzFreeObj(sHist);
    }
    if(dHist)
    {
      WlzFreeObj(dHist);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Fit the given Woolz histogram object's histogram
*               bin values using a mixture of Gaussian distributions.
*
*               First guess at peak positions, the method of guessing
*               depends on the parameter numDbn. If (numDbn > 0) guess
*               numDbn equally spaced Gaussians, else use the function
*               WlzHistogramFindPeaks().
*               After initial guess use the function AlgMixture() to find
*               the maximum liklihood mixture of Gaussian distributions
*               which best fit the histogram
*		\f[
                  h(x) = \sum_{i=0}^k {\alpha_i g_i(x)}
		\f]
*               where
*		\f[
                   \sum_{i=0}^k {\alpha_i} = 1.0
		\f]
*                 \f$g_i(x)\f$ are the gaussians defined by standard deviation
*                 \f$sigma_i\f$ and mean \f$mu_i\f$.
* \param	histObj			Given histogram object.
* \param	numDbn			If greater than zero the number
*                                       of Gaussians to fit. If less
*                                       than or equal to zero then the
*                                       Gaussians are fitted to peaks
*                                       found by WlzHistogramFindPeaks().
* \param	smooth			Gaussian \f$sigma\f$ smooth histogram
*                                       and compute 1st derivative.
* \param	thresh			Mimimum bin count for a peak.
* \param	tol			Log-likihood tolerance value
*                                       for terminating the iteration.
* \param	dstNumMu		Destination ptr for the number
*                                       of Gaussian mean values.
* \param	dstMu			Destination ptr for \f$mu\f$, the
*					Gaussian mean values.
* \param	dstNumSigma		Destination ptr for the number
*                                       of Gaussian \f$sigma\f$ values.
* \param	dstSigma		Destination ptr for the Gaussian
*                                       \f$sigma\f$ values.
* \param	dstNumAlpha		Destination ptr for the Gaussian
*                                       \f$alpha\f$ values.
* \param	dstAlpha		Destination ptr for the
*                                       log-liklihood, may be NULL.
* \param	dstLL
*/
WlzErrorNum     WlzHistogramFitPeaks(WlzObject *histObj, int numDbn,
				     double smooth, double thresh, double tol,
				     int *dstNumMu, double **dstMu,
				     int *dstNumSigma, double **dstSigma,
				     int *dstNumAlpha, double **dstAlpha,
				     double *dstLL)
{
  int		idx,
		fitNI,
  		pkCnt = 0,
		binCnt = 0;
  double	fitLL,
		fitObs = 0.0,
		pkHgt,
		pkHWHM,
  		sumAlpha;
  int		*pkIdx = NULL;
  double	*pkFreq = NULL,
  		*pkApha = NULL,
		*pkMu = NULL,
		*pkSigma = NULL;
  WlzObject	*smHistObj = NULL;
  WlzHistogramDomain *histDom = NULL,
  		*smHistDom = NULL;
  AlgError 	algErr = ALG_ERR_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstNumMu == NULL) || (dstMu == NULL) ||
     (dstNumSigma == NULL) || (dstSigma == NULL) ||
     (dstNumAlpha == NULL) || (dstAlpha == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzHistogramCheckHistObj(histObj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    histDom = histObj->domain.hist;
    if(histDom->nBins < 3)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    smHistObj = WlzCopyObject(histObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzHistogramCnvGauss(smHistObj, smooth, 0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    smHistDom =smHistObj->domain.hist;
    if(numDbn > 0)
    {
      if(numDbn > (smHistDom->nBins / 3))
      {
        errNum = WLZ_ERR_PARAM_DATA;
      }
      else
      {
	pkCnt = numDbn;
	if((pkIdx = (int *)AlcMalloc(sizeof(int) * pkCnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idx = 0; idx < numDbn; ++idx)
	  {
	    *(pkIdx + idx) = (smHistDom->nBins * (idx + 1)) / (numDbn + 1);
	  }
	}
      }
    }
    else
    {
      errNum = WlzHistogramFindPeaks(histObj, smooth, thresh,
  				     &pkCnt, &pkIdx, WLZ_HIST_FEATURE_PEAK);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    binCnt = smHistDom->nBins;
    if(((pkFreq = (double *)AlcMalloc(sizeof(double) * binCnt)) == NULL) ||
       ((pkApha = (double *)AlcMalloc(sizeof(double) * pkCnt)) == NULL) ||
       ((pkMu = (double *)AlcMalloc(sizeof(double) * pkCnt)) == NULL) ||
       ((pkSigma = (double *)AlcMalloc(sizeof(double) * pkCnt)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill in distribution parameters. */
    fitObs = 0.0;
    switch(smHistDom->type)
    {
      case WLZ_HISTOGRAMDOMAIN_INT:
	for(idx = 0; idx < binCnt; ++idx)
	{
	  fitObs += *(pkFreq + idx) = *(smHistDom->binValues.inp + idx);
	}
	break;
      case WLZ_HISTOGRAMDOMAIN_FLOAT:
	for(idx = 0; idx < binCnt; ++idx)
	{
	  fitObs += *(pkFreq + idx) = *(smHistDom->binValues.dbp + idx);
	}
        break;
      default:
        break;
    }
    if(pkCnt == 0)
    {
      *pkMu = smHistDom->origin + (binCnt * smHistDom->binSize * 0.5);
      *pkSigma = smooth * 2.0;
      *pkApha = 1.0;
    }
    else if(pkCnt == 1)
    {
      *pkMu = smHistDom->origin + (*pkIdx * smHistDom->binSize);
      *pkSigma = smooth * 2.0;
      *pkApha = 1.0;
    }
    else
    {
      /* Find sum of peak heights. */
      sumAlpha = 0.0;
      for(idx = 0; idx < pkCnt; ++idx)
      {
        sumAlpha += *(pkFreq + *(pkIdx + idx));
      }
      /* Set initial estimates */
      pkHgt = *(pkFreq + *pkIdx);
      *pkMu = smHistDom->origin + (smHistDom->binSize * *pkIdx);
      if(pkHgt > *(pkFreq + *(pkIdx + 1)) * 0.5)
      {
	/* Find half height half width */
	idx = *pkIdx;
	pkHWHM = 1;
	while((idx < *(pkIdx + 1)) && (*(pkFreq + idx) > pkHgt / 2))
	{
	  ++pkHWHM;
	  ++idx;
	}
	/* Convert half height half width to sigma */
	*pkSigma = 1.2 * pkHWHM;
      }
      else
      {
	*pkSigma = (*(pkIdx + 1) - *pkIdx) *
		   0.25 * smooth * smHistDom->binSize;
      }
      *pkApha = pkHgt / sumAlpha;
      for(idx = 1; idx < (pkCnt - 1); ++idx)
      {
	*(pkMu + idx) = smHistDom->origin +
			(smHistDom->binSize * *(pkIdx + idx));
	*(pkSigma + idx) = (*(pkIdx + idx + 1) - *(pkIdx + idx - 1)) *
			   smHistDom->binSize * 0.25 * smooth;

        *(pkApha + idx) = *(pkFreq + *(pkIdx + idx)) / sumAlpha;
      }
      pkHgt = *(pkFreq + *(pkIdx + pkCnt - 1));

      *(pkMu + pkCnt - 1) = smHistDom->origin +
      			    (smHistDom->binSize * *(pkIdx + pkCnt - 1));

      if(pkHgt > *(pkFreq + *(pkIdx + pkCnt - 2)) * 0.5)
      {
        /* Find half height half width */
	idx = *(pkIdx + pkCnt - 1);
	pkHWHM = 1;
	while((idx > *(pkIdx + pkCnt - 2)) && (*(pkFreq + idx) > pkHgt * 0.5))
	{
	  ++pkHWHM;
	  --idx;
	}
	/* Convert half height half width to sigma */
	*(pkSigma + pkCnt - 1) = 1.2 * pkHWHM;
      }
      else
      {
	 *(pkSigma + pkCnt - 1) = (*(pkIdx + pkCnt - 1) -
	 			   *(pkIdx + pkCnt - 2)) * smHistDom->binSize;
      }
      *(pkApha + pkCnt - 1) = *(pkFreq + *(pkIdx + pkCnt - 1)) / sumAlpha;
    }
  }
  if(pkIdx)
  {
    AlcFree(pkIdx);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fit mixture of Gaussians to the histogram. */
    algErr = AlgMixtureMLG(pkCnt, binCnt,
    			   smHistDom->origin, smHistDom->binSize, pkFreq,
			   pkApha, pkMu, pkSigma,
			   tol, fitObs, &fitLL, &fitNI);
    errNum = WlzErrorFromAlg(algErr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNumMu = pkCnt;
    *dstMu = pkMu;
    *dstNumSigma = pkCnt;
    *dstSigma = pkSigma;
    *dstNumAlpha = pkCnt;
    *dstAlpha = pkApha;
    if(dstLL)
    {
      *dstLL = fitLL;
    }
  }
  else
  {
    if(pkApha)
    {
      AlcFree(pkApha);
    }
    if(pkMu)
    {
      AlcFree(pkMu);
    }
    if(pkSigma)
    {
      AlcFree(pkSigma);
    }
  }
  if(pkFreq)
  {
    AlcFree(pkFreq);
  }
  if(smHistObj)
  {
    WlzFreeObj(smHistObj);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Normalise the given Woolz histogram object's histogram
*               bin values to the given range [0-maxVal]. If the given
*               maximum value is 0.0 then the histogram bin values are
*               normalised so that maxVal is equal to the number of
*               bins.
* \param	histObj			Given histogram object.
* \param	maxVal			Given maximum value.
*/
WlzErrorNum	WlzHistogramNorm(WlzObject *histObj, double maxVal)
{
  int		tI0,
		count,
		oldMaxValI;
  double	tD0,
		oldMaxValD;
  int		*histBinPI;
  double	*histBinPD;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramNorm FE %p %g\n",
	   histObj, maxVal));
  if(((errNum = WlzHistogramCheckHistObj(histObj)) == WLZ_ERR_NONE) &&
     ((histDom = histObj->domain.hist)->nBins > 0))
  {
    if((maxVal >= -(DBL_EPSILON)) && (maxVal <= DBL_EPSILON))
    {
      maxVal = histDom->origin + ((histDom->nBins  - 1) * histDom->binSize);
    }
    if(histDom->type == WLZ_HISTOGRAMDOMAIN_INT)
    {
      count = histDom->nBins;
      histBinPI = histDom->binValues.inp;
      oldMaxValI = *histBinPI;
      while(--count > 0)
      {
	if((tI0 = *++histBinPI) > oldMaxValI)
	{
	  oldMaxValI = tI0;
	}
      }
      tD0 = maxVal / oldMaxValI;
      count = histDom->nBins;
      histBinPI = histDom->binValues.inp;
      while(count-- > 0)
      {
        *histBinPI = WLZ_NINT(*histBinPI * tD0);
	++histBinPI;
      }
    }
    else
    {
      count = histDom->nBins;
      histBinPD = histDom->binValues.dbp;
      oldMaxValD = *histBinPD;
      while(--count > 0)
      {
	if((tD0 = *++histBinPD) > oldMaxValD)
	{
	  oldMaxValD = tD0;
	}
      }
      tD0 = maxVal / oldMaxValD;
      count = histDom->nBins;
      histBinPD = histDom->binValues.dbp;
      while(count-- > 0)
      {
        *histBinPD++ *= tD0;
      }
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramNorm FX %d\n",
	   (int )errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Computes the sum of the histogram domain's bin values.
* \param	histDom			Given histogram domain.
*/
double		WlzHistogramBinSum(WlzHistogramDomain *histDom)
{
  int		binCount;
  double	binSum = 0.0;
  WlzGreyP	binVal;

  if(histDom)
  {
    binCount = histDom->nBins;
    binVal = histDom->binValues;
    switch(histDom->type)
    {
      case WLZ_HISTOGRAMDOMAIN_INT:
	while(binCount-- > 0)
	{
	  binSum += *(binVal.inp)++;
	}
	break;
      case WLZ_HISTOGRAMDOMAIN_FLOAT:
	while(binCount-- > 0)
	{
	  binSum += *(binVal.dbp)++;
	}
	break;
      default:
        break;
    }
  }
  return(binSum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Finds the bin of the histogram domain which has the
*               greatest occupancy.
* \param	histDom			Given histogram domain.
*/
int		WlzHistogramBinMax(WlzHistogramDomain *histDom)
{
  int		binIdx = 0,
		maxIdx = 0,
  		maxValI,
		binCount;
  double	maxValD;
  WlzGreyP	binVal;

  if(histDom && ((binCount = histDom->nBins) > 1))
  {
    binVal = histDom->binValues;
    switch(histDom->type)
    {
      case WLZ_HISTOGRAMDOMAIN_INT:
	maxValI = *(binVal.inp);
	while(++binIdx < binCount)
	{
	  if(*++(binVal.inp) > maxValI)
	  {
	    maxValI = *(binVal.inp);
	    maxIdx = binIdx;
	  }
	}
	break;
      case WLZ_HISTOGRAMDOMAIN_FLOAT:
	maxValD = *(binVal.dbp);
	while(++binIdx < binCount)
	{
	  if(*++(binVal.dbp) > maxValD)
	  {
	    maxValD = *(binVal.dbp);
	    maxIdx = binIdx;
	  }
	}
	break;
      default:
        break;
    }
  }
  return(maxIdx);
}

/*!
* \return	Histogram distance measure, range [0.0 - 1.0].
* \ingroup	WlzHistogram
* \brief	Calculates a distance measure for comparing two
*               histograms. The histogram distance is in the range
*               [0.0 - 1.0] with 1.0 being a perfect match between
*               histograms.
*               The given histograms are not modified.
* \param	histObj0		First histogram object.
* \param	histObj1		Second histogram object.
* \param	dstErrNum		Destination error pointer, may be NULL.
*/
double		WlzHistogramDistance(WlzObject *histObj0,
				     WlzObject *histObj1,
				     WlzErrorNum *dstErrNum)
{
  int		count,
		idx,
		idx0,
		idx1;
  double	tD0,
  		tD1,
		tD2,
		tD3,
		tD4,
		tD5,
		coefA = 0,
		coefB = 0,
		coefC = 0,
		newOrigin,
		newBinSize,
		histDistance = 0.0;
  WlzGreyP	binVal0,
  		binVal1;
  WlzHistogramDomain *histDom0,
  		*histDom1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramDistance FE %p %p %p\n",
	   histObj0, histObj1, dstErrNum));
  if(((errNum = WlzHistogramCheckHistObj(histObj0)) == WLZ_ERR_NONE) &&
     ((errNum = WlzHistogramCheckHistObj(histObj1)) == WLZ_ERR_NONE) &&
     ((histDom0 = histObj0->domain.hist)->nBins > 0) &&
     ((histDom1 = histObj1->domain.hist)->nBins > 0))
  {
    binVal0 = histDom0->binValues;
    binVal1 = histDom1->binValues;
    if((histDom0->type != histDom1->type) ||
       ((tD0 = histDom0->origin - histDom1->origin) < DBL_EPSILON) ||
       (tD0 > DBL_EPSILON) ||
       ((tD1 = histDom0->nBins -  histDom1->nBins) < DBL_EPSILON) ||
       (tD1 > DBL_EPSILON))
    {
      /* Histograms need re-binning when calculating distance */
      newOrigin = WLZ_MIN(histDom0->origin, histDom1->origin);
      newBinSize = WLZ_MIN(histDom0->binSize, histDom1->binSize);
      tD0 = ((histDom0->nBins - 1) * histDom0->binSize) + histDom0->origin;
      tD1 = ((histDom1->nBins - 1) * histDom1->binSize) + histDom1->origin;
      tD2 = newBinSize / histDom0->binSize;
      tD3 = (newOrigin - histDom0->origin) / histDom0->binSize;
      tD4 = newBinSize / histDom1->binSize;
      tD5 = (newOrigin - histDom1->origin) / histDom1->binSize;
      for(idx = 0; idx < newBinSize; ++idx)
      {
        idx0 = WLZ_NINT((idx * tD2) + tD3);
	if((idx0 < 0) || (idx0 >= histDom0->nBins))
	{
	  tD0 = 0.0;
	}
	else
	{
	  tD0 = (histDom0->type == WLZ_HISTOGRAMDOMAIN_INT)?
		*(binVal0.inp + idx0): *(binVal0.dbp + idx0);
	}
	idx1 = (int )((idx * tD4) + tD5);
	if((idx1 < 0) || (idx1 >= histDom1->nBins))
	{
	  tD1 = 0.0;
	}
	else
	{
	  tD1 = (histDom1->type == WLZ_HISTOGRAMDOMAIN_INT)?
		*(binVal1.inp + idx1): *(binVal1.dbp + idx1);
	}
	coefA += tD0 * tD0;
	coefB += tD1 * tD1;
	coefC += tD0 * tD1;
      }
    }
    else
    {
      count = histDom1->nBins;
      if(histDom0->type == WLZ_HISTOGRAMDOMAIN_INT)
      {
	while(count-- > 0)
	{
	  tD0 = *(binVal0.inp)++;
	  tD1= *(binVal1.inp)++;
	  coefA += tD0 * tD0;
	  coefB += tD1 * tD1;
	  coefC += tD0 * tD1;
	}
      }
      else
      {
	while(count-- > 0)
	{
	  tD0 = *(binVal0.dbp)++;
	  tD1= *(binVal1.dbp)++;
	  coefA += tD0 * tD0;
	  coefB += tD1 * tD1;
	  coefC += tD0 * tD1;
	}
      }
      tD0 = coefA + coefB - coefC;
      histDistance = (WLZ_ABS(tD0) > DBL_EPSILON)? (coefC / tD0): (1.0);
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramDistance 01 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramDistance FX %g\n",
	   histDistance));
  return(histDistance);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Modifies the grey values of the given Woolz 2D or 3D
*               domain object so that the resulting object's histogram
*               matches the target histogram object.
*
*               If the given domain object is a 3D object and the
*               independentPlanes flag is set then each of the 3D
*               object's planes is matched independently.
*               A domain object or plane will only be matched to the
*               given histogram if the distance between it's histogram
*               and the given histogram is within the given distance
*               range.
*               If the specified smoothing is non-zero then the source
*               objects histogram is smoothed by WlzHistogramSmooth().
*               The histogram distance is in the range [0.0 - 1.0] with
*               1.0 being a perfect match between histograms. If the
*               given minimum distance is <= 0.0 and the maximum
*               distance is >= 1.0 then matching will always be done.
*               The given histogram is not modified.
* \param	srcObj			Given domain object.
* \param	targetHist		Target histogram object.
* \param	independentPlanes	If srcObj is a 3D domain object
*                                       and independentPlanes is set
*                                       (non-zero) then the domain
*                                       objects planes are matched
*                                       independently.
* \param	smoothing		Gaussian smoothing of the
*                                       source object histogram before
*                                       matching.
* \param	minDist			Minimum distance between the
*                                       histograms.
* \param	maxDist			Maximum distance between the
*                                       histograms.
* \param	dither			If non zero then dither mapped
*                                       values to avoid generating a
*                                       histogram with a comb appearence.
*/
WlzErrorNum	WlzHistogramMatchObj(WlzObject *srcObj, WlzObject *targetHist,
				     int independentPlanes, int smoothing,
				     double minDist, double maxDist,
				     int dither)
{
  int		tI0,
  		planeIdx,
  		planeCount,
		matchFlag = 0,
		isTiled = 0;
  double 	dist;
  WlzGreyType	greyType;
  WlzObject	*targetHistCum = NULL,
  		*targetHistNew = NULL,
  		*srcObjHist = NULL,
  		*srcObj2D = NULL;
  WlzHistogramDomain *histDom;
  WlzDomain	dummyDom;
  WlzValues	dummyValues;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramMatchObj FE %p %p %d %g %g %d\n",
	   srcObj, targetHist,
	   independentPlanes, minDist, maxDist, dither));
  if(((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						srcObj)) == WLZ_ERR_NONE) &&
     ((errNum = WlzHistogramCheckHistObj(targetHist)) == WLZ_ERR_NONE) &&
     ((histDom = targetHist->domain.hist)->nBins > 0))
  {
    if((errNum = WLZ_ERR_NONE) && isTiled)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
    else if((srcObj->type == WLZ_2D_DOMAINOBJ) ||
            (srcObj->type == WLZ_3D_DOMAINOBJ))
    {
      tI0 = (int )ceil((histDom->binSize * histDom->nBins) + histDom->origin);
      if(tI0 < 256)
      {
        tI0 = 256;
      }
      targetHistNew = WlzHistogramRebin(targetHist, WLZ_HISTOGRAMDOMAIN_INT,
					tI0, tI0, 0.0, 1.0, &errNum);
    }
  }
  if((errNum == WLZ_ERR_NONE) && (histDom->nBins > 0))
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        break;
      case WLZ_TRANS_OBJ:
	errNum = WlzHistogramMatchObj(srcObj->values.obj, targetHist,
				      independentPlanes, smoothing,
				      minDist, maxDist, dither);
        break;
      case WLZ_2D_DOMAINOBJ:
	histDom = targetHistNew->domain.hist;
	srcObjHist = WlzHistogramObj(srcObj, histDom->nBins,
				     0.0, 1.0, &errNum);
	if((smoothing > 0) && (errNum == WLZ_ERR_NONE))
	{
	  errNum = WlzHistogramSmooth(srcObjHist, smoothing);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(srcObjHist->domain.hist->nBins <= 0)
	  {
	    matchFlag = 0;
	  }
	  else
	  {
	    if((minDist >= DBL_EPSILON) || (maxDist <= (1.0 - DBL_EPSILON)))
	    {
	      dist = WlzHistogramDistance(targetHistNew, srcObjHist, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		matchFlag = (dist >= minDist) && (dist <= maxDist);
	      }
	    }
	    else
	    {
	      matchFlag = 1;
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && matchFlag)
	{
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (void )WlzHistogramCummulative(targetHistNew);
	    (void )WlzHistogramCummulative(srcObjHist);
	    errNum = WlzHistogramMatchPlane(srcObj, histDom,
	    				    srcObjHist->domain.hist, dither);
	  }
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	histDom = targetHistNew->domain.hist;
	dummyDom.core = NULL;
	dummyValues.core = NULL;
	if(((srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
				    dummyValues, NULL, NULL,
				    &errNum)) == NULL) &&
	   (errNum == WLZ_ERR_NONE))
	{
	  errNum = WLZ_ERR_UNSPECIFIED;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  targetHistCum = WlzHistogramCopy(targetHistNew, histDom->type,
	  				   &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(independentPlanes == 0)
	  {
	    srcObjHist = WlzHistogramObj(srcObj, histDom->nBins, 0.0, 1.0,
					 &errNum);
	    if((smoothing > 0) && (errNum == WLZ_ERR_NONE))
	    {
	      errNum = WlzHistogramSmooth(srcObjHist, smoothing);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(srcObjHist->domain.hist->nBins <= 0)
	      {
		matchFlag = 0;
	      }
	      else
	      {
		if((minDist >= DBL_EPSILON) ||
		   (maxDist <= (1.0 - DBL_EPSILON)))
		{
		  dist = WlzHistogramDistance(targetHistNew, srcObjHist,
		  			      &errNum);
		  matchFlag = (dist >= minDist) && (dist <= maxDist);
		}
		else
		{
		  matchFlag = 1;
		}
	      }
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (independentPlanes || matchFlag))
	{
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (void )WlzHistogramCummulative(targetHistCum);
	    if(independentPlanes == 0)
	    {
	      (void )WlzHistogramCummulative(srcObjHist);
	    }
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (independentPlanes || matchFlag))
	{
	  planeIdx =  0;
	  planeCount = srcObj->domain.p->lastpl - srcObj->domain.p->plane1 + 1;
	  while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	  {
	    srcObj2D->domain = *(srcObj->domain.p->domains + planeIdx);
	    srcObj2D->values = *(srcObj->values.vox->values + planeIdx);
	    if((srcObj2D->values.core) &&
	       (srcObj2D->values.core->type != WLZ_EMPTY_OBJ))
	    {
	      if(independentPlanes)
	      {
		if(srcObjHist)
		{
		  WlzFreeObj(srcObjHist);
		  srcObjHist = NULL;
		}
		srcObjHist = WlzHistogramObj(srcObj2D,
					     histDom->nBins, 0.0, 1.0,
					     &errNum);
		if((smoothing > 0) && (errNum == WLZ_ERR_NONE))
		{
		  errNum = WlzHistogramSmooth(srcObjHist, smoothing);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  if(srcObjHist->domain.hist->nBins <= 0)
		  {
		    matchFlag = 0;
		  }
		  else
		  {
		    if((minDist >= DBL_EPSILON) ||
		       (maxDist <= (1.0 - DBL_EPSILON)))
		    {
		      dist = WlzHistogramDistance(targetHistNew, srcObjHist,
						  &errNum);
		      matchFlag = (dist >= minDist) && (dist <= maxDist);
		    }
		    else
		    {
		      matchFlag = 1;
		    }
		  }
		}
	      }
	      if((errNum == WLZ_ERR_NONE) && matchFlag)
	      {
		if(independentPlanes)
		{
		  (void )WlzHistogramCummulative(srcObjHist);
		}
		errNum = WlzHistogramMatchPlane(srcObj2D,
						targetHistCum->domain.hist,
						srcObjHist->domain.hist,
						dither);
	      }
	    }
	    ++planeIdx;
	  }
	  srcObj2D->domain.core = NULL;
	  srcObj2D->values.core = NULL;
	  WlzFreeObj(srcObj2D);
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
    if(targetHistNew)
    {
      WlzFreeObj(targetHistNew);
    }
    if(targetHistCum)
    {
      WlzFreeObj(targetHistCum);
    }
    if(srcObjHist)
    {
      WlzFreeObj(srcObjHist);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramMatchObj FX %d\n",
	   (int )errNum));
  return(errNum);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzHistogram
* \brief	Modifies the grey values of the given Woolz 2D or 3D
*               domain object so that the resulting object's histogram
*               approximates a uniform histogram.
*               If the specified smoothing is non-zero then the objects
*               histogram is smoothed by WlzHistogramSmooth().
* \param	srcObj			Given domain object.
* \param	smoothing		Gaussian smoothing of the
*                                       histogram before equalisation.
* \param	dither			If non zero then dither mapped
*                                       values to avoid generating a
*                                       histogram with a comb appearence.
*/
WlzErrorNum	WlzHistogramEqualiseObj(WlzObject *srcObj, int smoothing,
					int dither)
{
  int		count,
  		isTiled = 0;
  int		*map;
  double	mass,
  		normFac;
  WlzGreyType	greyType;
  WlzObject	*histObj = NULL;
  WlzHistogramDomain *histDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzHistogramEqualiseObj FE %p %d %d\n",
	   srcObj, smoothing, dither));
  if((errNum = WlzHistogramCheckDomainAndValues(&greyType, &isTiled,
						srcObj)) == WLZ_ERR_NONE)
  {
    if(isTiled)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        break;
      case WLZ_TRANS_OBJ:
	errNum = WlzHistogramEqualiseObj(srcObj->values.obj, smoothing,
					 dither);
        break;
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
        histObj = WlzHistogramObj(srcObj, 0, 0.0, 1.0, &errNum);
	if((errNum == WLZ_ERR_NONE) && (smoothing > 0))
	{
	  errNum = WlzHistogramSmooth(histObj, smoothing);
	}
        if(errNum == WLZ_ERR_NONE)
        {
	  (void )WlzHistogramCummulative(histObj);
	  histDom = histObj->domain.hist;
	  map = histDom->binValues.inp;
	  mass = *(map + histDom->nBins - 1);
	  normFac = (fabs(mass) > DBL_EPSILON)? (255.0 / mass): 0.0;
	  count = histDom->nBins;
	  while(count-- > 0)
	  {
	    *map = WLZ_NINT(*map * normFac);
	    ++map;
	  }
          errNum = WlzHistogramMapValues(srcObj, histObj, dither);
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzHistogramEqualiseObj FX %d\n",
	   (int )errNum));
  return(errNum);
}
