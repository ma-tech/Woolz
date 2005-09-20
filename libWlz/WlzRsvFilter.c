#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzRsvFilter.c
* \author       Bill Hill
* \date         May 1999
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
* \brief	Recursive filters for Woolz, these are also known as
* 		Deriche and infinite impulse responce (IIR) filters.
* \ingroup	WlzValueFilters
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

/* These tests are for debuging only. */
/* #define WLZ_RSVFILTER_TEST_1D */
/* #define WLZ_RSVFILTER_TEST_2D */
/* #define WLZ_RSVFILTER_TEST_3D */

static WlzObject *WlzRsvFilterObj2DX(WlzObject *, WlzRsvFilter *,
				     WlzErrorNum *);
static WlzObject *WlzRsvFilterObj2DY(WlzObject *, WlzRsvFilter *,
				     WlzErrorNum *);
static WlzObject *WlzRsvFilterObj3DXY(WlzObject *, WlzRsvFilter *,
			              int, WlzErrorNum *);
static WlzObject *WlzRsvFilterObj3DZ(WlzObject *, WlzRsvFilter *,
				     WlzErrorNum *);
static void	WlzRsvFilterFilterBufXF(WlzRsvFilter *,
				      double *, double *, double *,
				      int);
static void	WlzRsvFilterFilterBufYF(WlzRsvFilter *,
				      double **, double **, UBYTE **,
				      WlzIVertex2, int, int);
static void	WlzRsvFilterFilterBufZF(WlzRsvFilter *ftr,
				 	double ***wrkBuf,
					double ***srcBuf,
					UBYTE ***itvBuf,
					WlzIVertex3 bufPos,
					int itvLen,
					int goingUp);

/*!
* \return	void
* \ingroup      WlzValueFilters
* \brief	Free a recursive filter.
* \param	ftr			Filter to free.
*/
void		WlzRsvFilterFreeFilter(WlzRsvFilter *ftr)
{
  if(ftr != NULL)
  {
    AlcFree(ftr);
  }
}

/*!
* \return	New recursive filter.
* \ingroup      WlzValueFilters
* \brief	Creates a new recursive filter.
* \param	name			Name of filter.
* \param	prm			Parameter is alpha for Deriche's
*                                       operators and sigma for
*                                       Gaussians, must be \f$\geq\f$ 0.0.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzRsvFilter *WlzRsvFilterMakeFilter(WlzRsvFilterName name,
				       double prm, WlzErrorNum *dstErr)
{
  int		idx;
  double	tD0,
		tAS,
		tWS,
  		tExp,
  		tExp2,
		tNorm,
		tCos,
		tSin;
  WlzRsvFilter *ftr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	cGauss = 2.4726;
  const double	aGauss[3] = { 1.2600,  0.9338,  0.6134},
  		bGauss[3] = { 0.9629, -0.1051,  1.2280},
		gGauss[3] = { 1.9420,  1.8980, -0.2976},
		wGauss[3] = { 0.8448,  0.9459,  1.3320};

  if((ftr = (WlzRsvFilter *)AlcCalloc(sizeof(WlzRsvFilter),
  					  1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(name == WLZ_RSVFILTER_NAME_NONE)
  {
    ftr->name = WLZ_RSVFILTER_NAME_NONE;
  }
  else if(prm < DBL_EPSILON)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(name)
    {
      case WLZ_RSVFILTER_NAME_DERICHE_0: /* FALLTHROUGH */
      case WLZ_RSVFILTER_NAME_DERICHE_1: /* FALLTHROUGH */
      case WLZ_RSVFILTER_NAME_DERICHE_2:
	tExp = exp(-prm);
	tExp2 = tExp * tExp;
	tD0 = 1.0 - tExp;
	tNorm = (tD0 * tD0) / (1.0 + (2.0 * prm * tExp) - tExp2);
	ftr->name = name;
	ftr->b[0] = -2.0 * tExp;
	ftr->b[1] = tExp2;
	switch(name)
	{
	  case WLZ_RSVFILTER_NAME_DERICHE_0:
	    ftr->a[0] = 1.0;
	    ftr->a[1] = (prm - 1.0) * tExp;
	    ftr->a[2] = (prm + 1.0) * tExp;
	    ftr->a[3] = -tExp2;
	    ftr->c = tNorm;
	    break;
	  case WLZ_RSVFILTER_NAME_DERICHE_1:
	    ftr->a[0] = 0.0;
	    tD0 = (prm * prm) * tExp;
	    ftr->a[1] = -tD0;
	    ftr->a[2] = tD0;
	    ftr->a[3] = 0.0;
	    ftr->c = -tNorm;
	    break;
	  case WLZ_RSVFILTER_NAME_DERICHE_2:
	    tD0 = (prm * prm);
	    ftr->a[0] = -tD0;
	    ftr->a[1] = (prm + 1) * tD0 * tExp;
	    ftr->a[2] = (prm - 1) * tD0 * tExp;
	    ftr->a[3] = tD0 * tExp2;
	    ftr->c = tNorm;
	    break;
	}
        break;
      case WLZ_RSVFILTER_NAME_GAUSS_0: /* FALLTHROUGH */
      case WLZ_RSVFILTER_NAME_GAUSS_1: /* FALLTHROUGH */
      case WLZ_RSVFILTER_NAME_GAUSS_2:
	ftr->name = name;
	switch(name)
	{
	  case WLZ_RSVFILTER_NAME_GAUSS_0:
	    idx = 0;
	    break;
	  case WLZ_RSVFILTER_NAME_GAUSS_1:
	    idx = 1;
	    break;
	  case WLZ_RSVFILTER_NAME_GAUSS_2:
	    idx = 2;
	    break;
	}
	tWS = wGauss[idx] / prm;
	tAS = aGauss[idx] / prm;
	tCos = cos(tWS);
	tSin = sin(tWS);
	tExp = exp(-tAS);
	tExp2 = tExp * tExp;
        tNorm = 1.0 / (cGauss * prm);
	ftr->a[0] = bGauss[idx];
	tD0 = ((-bGauss[idx] * tCos) + (gGauss[idx] * tSin)) * tExp;
	ftr->a[1] = tD0;
	ftr->b[0] = -2.0 * tCos * tExp;
	ftr->b[1] = tExp2;
	if(name == WLZ_RSVFILTER_NAME_GAUSS_1) 	        /* Asymetric */
	{
	  ftr->a[3] = -(tD0 + (2.0 * bGauss[idx] * tCos * tExp));
	  ftr->a[2] = bGauss[idx] * tExp2;
	}
	else
	{
	  ftr->a[2] = tD0 + (2.0 * bGauss[idx] * tCos * tExp);
	  ftr->a[3] = -(bGauss[idx] * tExp2);
	}
	ftr->c = tNorm;
        break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ftr);
}

/*!
* \return	The filtered object, or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Applies a recursive filter to the given object.
* \param	srcObj			Given object.
* \param	ftr			Recursive filter.
* \param	actionMsk		Action mask.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzObject	*WlzRsvFilterObj(WlzObject *srcObj, WlzRsvFilter *ftr,
			         int actionMsk, WlzErrorNum *dstErr)
{
  WlzValues	tVal;
  WlzObject	*xObj = NULL,
		*yObj = NULL,
		*xyObj = NULL,
		*zObj = NULL,
  		*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcObj == NULL) || (ftr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	if((actionMsk & (WLZ_RSVFILTER_ACTION_X |
			 WLZ_RSVFILTER_ACTION_Y)) != 0)
	{
	  /* Filter in each required direction. */
	  if((actionMsk & WLZ_RSVFILTER_ACTION_X) != 0)
	  {
	    xObj = WlzRsvFilterObj2DX(srcObj, ftr, &errNum);
	  }
	  if((errNum == WLZ_ERR_NONE) &&
	     ((actionMsk & WLZ_RSVFILTER_ACTION_Y) != 0))
	  {
	    yObj = WlzRsvFilterObj2DY((xObj)? xObj: srcObj,
	    		              ftr, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(yObj)
	    {
	      dstObj = yObj;
	      yObj = NULL;
	    }
	    else if(xObj)
	    {
	      dstObj = xObj;
	      xObj = NULL;
	    }
	  }
	  if(xObj)
	  {
	    WlzFreeObj(xObj);
	  }
	  if(yObj)
	  {
	    WlzFreeObj(yObj);
	  }
	}
	else
	{
	  /* No filtering required. */
	  dstObj = WlzNewGrey(srcObj, &errNum);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	if((actionMsk & (WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y |
			  WLZ_RSVFILTER_ACTION_Z)) != 0)
	{
	  if((actionMsk & (WLZ_RSVFILTER_ACTION_X | 
	  		   WLZ_RSVFILTER_ACTION_Y)) != 0)
	  {
	    xyObj = WlzRsvFilterObj3DXY(srcObj, ftr, actionMsk,
	    				&errNum);
	  }
	  if((errNum == WLZ_ERR_NONE) && 
	     ((actionMsk & WLZ_RSVFILTER_ACTION_Z) != 0))
	  {
	    zObj = WlzRsvFilterObj3DZ((xyObj)? xyObj: srcObj, ftr, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(zObj)
	    {
	      dstObj = zObj;
	      zObj = NULL;
	    }
	    else if(xyObj)
	    {
	      dstObj = xyObj;
	      xyObj = NULL;
	    }
	  }
	  if(xyObj)
	  {
	    WlzFreeObj(xyObj);
	  }
	  if(zObj)
	  {
	    WlzFreeObj(zObj);
	  }
	}
	else
	{
	  /* No filtering required. */
	  tVal = WlzCopyValues(srcObj->type, srcObj->values, srcObj->domain,
	  		       &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj= WlzMakeMain(srcObj->type, srcObj->domain, tVal,
	    			NULL, NULL, &errNum);
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	void
* \ingroup	WlzValueFilters
* \brief	 Filters the given single line data buffer using an
*               IIR filter defined by the filter coefficients and
*               double precision arithmetic.
* \note         Before optimization this and all the other IIR filter
*               code looks like:
*		\verbatim
                  for(i = 2; i < dataSz; ++i)
                  {
                    buf0[i] =   (a0 * data[i + 0]) + (a1 * data[i - 1])
                              - (b0 * buf0[i - 1]) - (b1 * buf0[i - 2]);
                  }
                  for(i = dataSz - 3; i >= 0; --i)
                  {
                    buf1[i] =   (a2 * data[i + 1]) + (a3 * data[i + 2])
                              - (b0 * buf1[i + 1]) - (b1 * buf1[i + 2]);
                  }
                  for(i = 0; i < dataSz; ++i)
                  {
                    data[i] =  c * (buf0[i] + buf1[i]);
                  }
		\endverbatim
* \param	ftr			The filter.
* \param	data			The line of data to be filtered.
* \param	buf0			Buffer with at least dataSz
*                                       elements.
* \param	buf1			Buffer with at least dataSz
*                                       elements.
* \param	dataSz			Number of data.
*/
static void	WlzRsvFilterFilterBufXF(WlzRsvFilter *ftr, double *data,
				      double *buf0, double *buf1,
				      int dataSz)
{
  int		cnt;
  double	a0,
  		a1,
		a2,
		a3,
		b0,
		b1,
		c,
		d0,
		d1,
		d2,
		f1,
		f2;
  double	*dP,
  		*fP0,
		*fP1;

  a0 = ftr->a[0];
  a1 = ftr->a[1];
  a2 = ftr->a[2];
  a3 = ftr->a[3];
  b0 = ftr->b[0];
  b1 = ftr->b[1];
  c = ftr->c;
  dP = data;
  fP0 = buf0;
  d0 = d1 = *dP;
  *fP0 = f1 = f2 = ((a0 + a1) * d0) / (b0 + b1 + 1);
  cnt = dataSz;
  while(--cnt > 0)
  {
    d1 = d0;
    d0 = *++dP;
    f2 = f1;
    f1 = *fP0;
    *++fP0 = (a0 * d0) + (a1 * d1) - (b0 * f1) - (b1 * f2);
  }
  dP = data + dataSz - 1;
  fP1 = buf1 + dataSz - 1;
  d2 = d1 = *dP;
  *fP1 = f1 = f2 = ((a2 + a3) * d1) / (b0 + b1 + 1);
  cnt = dataSz;
  while(--cnt > 0)
  {
    d1 = d2;
    d2 = *dP--;
    f1 = f2;
    f2 = *fP1;
    *--fP1 = (a2 * d2) + (a3 * d1) - (b0 * f2) - (b1 * f1);
  }
  dP = data;
  fP0 = buf0; 
  fP1 = buf1;
  cnt = dataSz;
  while(cnt-- > 0)
  {
    *dP++ = c * (*fP0++ + *fP1++);
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValueFilters
* \brief	Filters a given buffer using an IIR filter defined by
*               the filter coefficients.
* \param	ftr			The filter.
* \param	datBuf			The data buffer to be filtered
*                                       with bufSz elements.
* \param	tmpBuf0			Temporary buffer with at least
*                                       dataSz elements.
* \param	tmpBuf1			Temporary buffer with at least
*                                       dataSz elements.
* \param	bufSz			Number of data in buffer.
*/
WlzErrorNum	WlzRsvFilterBuffer(WlzRsvFilter *ftr, double *datBuf,
				   double *tmpBuf0, double *tmpBuf1,
				   int bufSz)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ftr == NULL) ||
     (datBuf == NULL) || (tmpBuf0 == NULL) || (tmpBuf1 == NULL) ||
     (bufSz <= 0))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    WlzRsvFilterFilterBufXF(ftr, datBuf, tmpBuf0, tmpBuf1, bufSz);
  }
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzValueFilters
* \brief	Filters a single interval using a vertical IIR filter
*               defined by the filter coefficients and double precision
*               arithmetic.
* \note         See the double precision horizontal function for
*               simple code.
* \param	ftr			The filter.
* \param	wrkBuf			Working buffer.
* \param	srcBuf			Source buffer.
* \param	itvBuf			Within interval bit buffer.
* \param	bufPos			Position within the buffer.
* \param	itvLen			Interval length.
* \param	goingUp			Non-zero if going up through
*                                       the lines (ie raster direction
*                                       is WLZ_RASTERDIR_DLIC).
*/
static void	WlzRsvFilterFilterBufYF(WlzRsvFilter *ftr,
				      double **wrkBuf, double **srcBuf,
				      UBYTE **itvBuf,
				      WlzIVertex2 bufPos, int itvLen,
				      int goingUp)
{
  int		cnt0,
		cnt1,
		kol,
		iBM,
		iBO,
		iBS,
  		idL0,
		idL1,
		idL2;
  double	a0,
  		a1,
		a2,
		a3,
		b0,
		b1,
		c,
		d0,
		d1,
		d2,
		f0,
		f1,
		f2;
  UBYTE		*iBP1,
  		*iBP2,
  		*iP1,
  		*iP2;
  double	*dP0,
  		*dP1,
		*dP2,
		*fP0,
		*fP1,
		*fP2;

  if(goingUp)
  {
    idL0 = (bufPos.vtY + 3 + 0) % 3;
    idL1 = (bufPos.vtY + 3 + 1) % 3;
    idL2 = (bufPos.vtY + 3 + 2) % 3;
    a2 = ftr->a[2];
    a3 = ftr->a[3];
  }
  else
  {
    idL0 = (bufPos.vtY + 3 - 0) % 3;
    idL1 = (bufPos.vtY + 3 - 1) % 3;
    idL2 = (bufPos.vtY + 3 - 2) % 3;
    a0 = ftr->a[0];
    a1 = ftr->a[1];
  }
  b0 = ftr->b[0];
  b1 = ftr->b[1];
  c = ftr->c;
  kol = bufPos.vtX;
  iP1 = *(itvBuf + idL1);
  iP2 = *(itvBuf + idL2);
  dP0 = *(srcBuf + idL0) + kol;
  fP0 = *(wrkBuf + idL0) + kol;
  dP1 = *(srcBuf + idL1) + kol;
  fP1 = *(wrkBuf + idL1) + kol;
  dP2 = *(srcBuf + idL2) + kol;
  fP2 = *(wrkBuf + idL2) + kol;
  cnt0 = itvLen;
  if(goingUp)
  {
    while(cnt0 > 0)
    {
      iBO = kol >> 3;
      iBS = kol & 7;
      iBP1 = iP1 + iBO;
      iBP2 = iP2 + iBO;
      if((iBS == 0) && ((*iBP1 & *iBP2) == 0xff) && (cnt0 >= 8))
      {
	cnt1 = 8;
	while(cnt1-- > 0)
	{
	  f0 = *fP0;
	  *fP0 = (a2 * *dP1++) + (a3 * *dP2++) - (b0 * *fP1++) - (b1 * *fP2);
	  *fP2++ = c * (f0 + *fP0++);
	}
	kol += 8;
	cnt0 -= 8;
      }
      else
      {
	d0 = *dP0;
	f0 = *fP0;
	iBM = 1 << iBS;
	if((*iBP1 & iBM == 0) || ((*iBP2 & iBM) == 0))
	{
	  d1 = d0;
	  d2 = d0;
	  f1 = ((a2 * a3) * d0) / (b0 + b1 + 1);
	  f2 = ((a2 * a3) * d2) / (b0 + b1 + 1);
	}
	else
	{
	  d1 = *dP1;
	  d2 = *dP2;
	  f1 = *fP1;
	  f2 = *fP2;
	}
	++dP0;
	++dP1;
	++dP2;
	++fP1;
	*fP0 = (a2 * d1) + (a3 * d2) - (b0 * f1) - (b1 * f2);
	*fP2++ = c * (f0 + *fP0++);
	--cnt0;
	++kol; 
      }
    }
  }
  else
  {
    while(cnt0 > 0)
    {
      iBO = kol >> 3;
      iBS = kol & 7;
      iBP1 = iP1 + iBO;
      iBP2 = iP2 + iBO; 
      if((iBS == 0) && ((*iBP1 & *iBP2) == 0xff) && (cnt0 >= 8))
      {
        cnt1 = 8;
	while(cnt1-- > 0)
	{
	  *fP0++ = (a0 * *dP0++) + (a1 * *dP1++) -
	  	   (b0 * *fP1++) - (b1 * *fP2++);
	}
	kol += 8;
	cnt0 -= 8;
      }
      else
      {
	d0 = *dP0;
	iBM = 1 << iBS;
	if((*iBP1 & iBM == 0) || ((*iBP2 & iBM) == 0))
	{
	  d1 = d0;
	  f1 = ((a0 + a1) * d0) / (b0 + b1 + 1);
	  f2 = ((a0 + a1) * d0) / (b0 + b1 + 1);
	}
	else
	{
	  d1 = *dP1;
	  f1 = *fP1;
	  f2 = *fP2;
	}
	++dP0;
	++dP1;
	++fP1;
	++fP2;
	*fP0++ = (a0 * d0) + (a1 * d1) - (b0 * f1) - (b1 * f2);
	--cnt0;
	++kol;
      }
    }
  }
}

/*!
* \return	void
* \ingroup	WlzValueFilters
* \brief	Filters a single interval using a vertical IIR filter
*               defined by the filter coefficients and double precision
*               arithmetic.
* \note         See the double precision horizontal function for
*               simple code.
* \param	ftr			The filter.
* \param	wrkBuf			Working buffer.
* \param	srcBuf			Source buffer.
* \param	itvBuf			Within interval bit buffer.
* \param	bufPos			Position within the buffer.
* \param	itvLen			Interval length.
* \param	goingUp			Non-zero if going up through
*                                       the planes.
*/
static void	WlzRsvFilterFilterBufZF(WlzRsvFilter *ftr,
				 	double ***wrkBuf,
					double ***srcBuf,
					UBYTE ***itvBuf,
					WlzIVertex3 bufPos,
					int itvLen,
					int goingUp)
{
  int		cnt0,
		cnt1,
		kol,
		lin,
		iBM,
		iBO,
		iBS,
  		idP0,
		idP1,
		idP2;
  double	a0,
  		a1,
		a2,
		a3,
		b0,
		b1,
		c,
		d0,
		d1,
		d2,
		f0,
		f1,
		f2;
  UBYTE		*iBP1,
  		*iBP2,
  		*iP1,
  		*iP2;
  double	*dP0,
  		*dP1,
		*dP2,
		*fP0,
		*fP1,
		*fP2;

  if(goingUp)
  {
    idP0 = (bufPos.vtZ + 3 + 0) % 3;
    idP1 = (bufPos.vtZ + 3 + 1) % 3;
    idP2 = (bufPos.vtZ + 3 + 2) % 3;
    a2 = ftr->a[2];
    a3 = ftr->a[3];
  }
  else
  {
    idP0 = (bufPos.vtZ + 3 - 0) % 3;
    idP1 = (bufPos.vtZ + 3 - 1) % 3;
    idP2 = (bufPos.vtZ + 3 - 2) % 3;
    a0 = ftr->a[0];
    a1 = ftr->a[1];
  }
  b0 = ftr->b[0];
  b1 = ftr->b[1];
  c = ftr->c;
  kol = bufPos.vtX;
  lin = bufPos.vtY;
  iP1 = *(*(itvBuf + idP1) + lin);
  iP2 = *(*(itvBuf + idP2) + lin);
  dP0 = *(*(srcBuf + idP0) + lin) + kol;
  fP0 = *(*(wrkBuf + idP0) + lin) + kol;
  dP1 = *(*(srcBuf + idP1) + lin) + kol;
  fP1 = *(*(wrkBuf + idP1) + lin) + kol;
  dP2 = *(*(srcBuf + idP2) + lin) + kol;
  fP2 = *(*(wrkBuf + idP2) + lin) + kol;
  cnt0 = itvLen;
  if(goingUp)
  {
    while(cnt0 > 0)
    {
      iBO = kol >> 3;
      iBS = kol & 7;
      iBP1 = iP1 + iBO;
      iBP2 = iP2 + iBO;
      if((iBS == 0) && ((*iBP1 & *iBP2) == 0xff) && (cnt0 >= 8))
      {
	cnt1 = 8;
	while(cnt1-- > 0)
	{
	  f0 = *fP0;
	  *fP0 = (a2 * *dP1++) + (a3 * *dP2++) - (b0 * *fP1++) - (b1 * *fP2);
	  *fP2++ = c * (f0 + *fP0++);
	}
	kol += 8;
	cnt0 -= 8;
      }
      else
      {
	d0 = *dP0;
	f0 = *fP0;
	iBM = 1 << iBS;
	if((*iBP1 & iBM == 0) || ((*iBP2 & iBM) == 0))
	{
	  d1 = d0;
	  d2 = d0;
	  f1 = ((a2 * a3) * d0) / (b0 + b1 + 1);
	  f2 = ((a2 * a3) * d2) / (b0 + b1 + 1);
	}
	else
	{
	  f1 = *fP1;
	  f2 = *fP2;
	  d1 = *dP1;
	  d2 = *dP2;
	}
	++dP0;
	++dP1;
	++dP2;
	++fP1;
	*fP0 = (a2 * d1) + (a3 * d2) - (b0 * f1) - (b1 * f2);
	*fP2++ = c * (f0 + *fP0++);
	--cnt0;
	++kol; 
      }
    }
  }
  else
  {
    while(cnt0 > 0)
    {
      iBO = kol >> 3;
      iBS = kol & 7;
      iBP1 = iP1 + iBO;
      iBP2 = iP2 + iBO; 
      if((iBS == 0) && ((*iBP1 & *iBP2) == 0xff) && (cnt0 >= 8))
      {
        cnt1 = 8;
	while(cnt1-- > 0)
	{
	  *fP0++ = (a0 * *dP0++) + (a1 * *dP1++) -
	  	   (b0 * *fP1++) - (b1 * *fP2++);
	}
	kol += 8;
	cnt0 -= 8;
      }
      else
      {
	d0 = *dP0;
	iBM = 1 << iBS;
	if((*iBP1 & iBM == 0) || ((*iBP2 & iBM) == 0))
	{
	  d1 = d0;
	  f1 = ((a0 + a1) * d0) / (b0 + b1 + 1);
	  f2 = ((a0 + a1) * d0) / (b0 + b1 + 1);
	}
	else
	{
	  d1 = *dP1;
	  f1 = *fP1;
	  f2 = *fP2;
	}
	++dP0;
	++dP1;
	++fP1;
	++fP2;
	*fP0++ = (a0 * d0) + (a1 * d1) - (b0 * f1) - (b1 * f2);
	--cnt0;
	++kol;
      }
    }
  }
}

/*!
* \return	The filtered object, or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Applies a recursive filter along the lines of the given
*               2D domain object with grey values using either double
*               precision floating point arithmetic or fixed
*               point arithmetic.
*               It is assumed that the object type has already been
*               checked, the domain and values are non-null.
* \param	srcObj			Given 2D domain object.
* \param	ftr			Recursive filter.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzRsvFilterObj2DX(WlzObject *srcObj, WlzRsvFilter *ftr,
				     WlzErrorNum *dstErr)
{
  int		bufSz,
		bufSpace,
		itvLen;
  WlzGreyType	bufType,
  		srcGType,
  		dstGType;
  WlzObjectType	vType;
  WlzGreyP	bufGP;
  WlzPixelV	bgdPix;
  WlzDomain	srcDom;
  WlzValues	dstVal;
  WlzObject	*dstObj = NULL;
  void		*lnBuf0 = NULL,
  		*lnBuf1 = NULL,
  		*datBuf = NULL;
  WlzIntervalWSpace srcIWSp,
  		dstIWSp;
  WlzGreyWSpace srcGWSp,
  		dstGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Gather information about the source object. */
  if(((srcDom = srcObj->domain).core->type != WLZ_INTERVALDOMAIN_INTVL) &&
     (srcDom.core->type != WLZ_INTERVALDOMAIN_RECT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    srcGType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcGType)
    {
      /* Promote grey type for destination object. */
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT:
        dstGType = WLZ_GREY_SHORT;
	break;
      case WLZ_GREY_INT:
	dstGType = WLZ_GREY_INT;
        break;
      case WLZ_GREY_FLOAT:
        dstGType = WLZ_GREY_FLOAT;
	break;
      case WLZ_GREY_DOUBLE:
        dstGType = WLZ_GREY_DOUBLE;
	break;
      case WLZ_GREY_RGBA: /* RGBA to be done RAB */
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vType = WlzGreyTableTypeToTableType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bgdPix = WlzGetBackground(srcObj, &errNum);
  }
  /* Make destination object with it's own values but a shared domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstVal.v = WlzNewValueTb(srcObj, WlzGreyTableType(vType, dstGType, NULL),
    			     bgdPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj= WlzMakeMain(srcObj->type, srcDom, dstVal, srcObj->plist,
			NULL, &errNum);
  }
  /* Make the working buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz = srcDom.i->lastkl - srcDom.i->kol1 + 1;
    bufSpace = sizeof(double) * bufSz;
    if(((datBuf = AlcMalloc(bufSpace)) == NULL) ||
       ((lnBuf0 = AlcMalloc(bufSpace)) == NULL) ||
       ((lnBuf1 = AlcMalloc(bufSpace)) ==  NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Work down through the object from the first line to last. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufType = WLZ_GREY_DOUBLE;
    bufGP.dbp = (double *)datBuf;
    if(((errNum = WlzInitGreyRasterScan(srcObj, &srcIWSp, &srcGWSp,
    					WLZ_RASTERDIR_ILIC,
					0)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyRasterScan(dstObj, &dstIWSp, &dstGWSp,
       					WLZ_RASTERDIR_ILIC,
					0)) == WLZ_ERR_NONE))
    {
      while((errNum == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&srcIWSp)) == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&dstIWSp)) == WLZ_ERR_NONE))
      {
	/* Copy interval to working buffer. */
	itvLen = srcIWSp.rgtpos - srcIWSp.lftpos + 1;
	WlzValueCopyGreyToGrey(bufGP, 0, bufType,
			       srcGWSp.u_grintptr, 0, srcGWSp.pixeltype,
			       itvLen);
	/* Apply filter. */
	WlzRsvFilterFilterBufXF(ftr, (double *)datBuf,
			        (double *)lnBuf0, (double *)lnBuf1,
				itvLen);
	/* Clamp data from buffer into the dst interval. */
	WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0, dstGWSp.pixeltype,
			          bufGP, 0, bufType, itvLen);
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  /* Free buffers. */
  if(datBuf)
  {
    AlcFree(datBuf);
  }
  if(lnBuf0)
  {
    AlcFree(lnBuf0);
  }
  if(lnBuf1)
  {
    AlcFree(lnBuf1);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      (void )WlzFreeObj(dstObj);
    }
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	The filtered object, or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Applies a recursive filter through the rows and columns
*               of the given 2D domain object with grey values using
*               either double precision floating point arithmetic or
*               fixed point arithmetic.
*               It is assumed that the object type has already been
*               checked, the domain and values are non-null.
* \param	srcObj			Given 2D domain object.
* \param	ftr			Recursive filter.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzRsvFilterObj2DY(WlzObject *srcObj, WlzRsvFilter *ftr,
			             WlzErrorNum *dstErr)
{
  int		idD,
  		idN,
		bufLnIdx,
		dstLnIdx,
  		itvLen,
		itvBufWidth;
  WlzRasterDir	rasDir;
  WlzGreyType	bufType,
  		srcGType,
		dstGType;
  WlzIVertex2	bufPos,
  		bufSz;
  WlzGreyP	dstBufGP,
  		srcBufGP,
  		wrkBufGP;
  WlzObjectType	vType;
  WlzPixelV	bgdPix;
  WlzDomain	srcDom;
  WlzValues	dstVal;
  WlzObject	*dstObj = NULL;
  void		**srcBuf = NULL,
  		**wrkBuf = NULL;
  UBYTE		**itvBuf = NULL;
  WlzIntervalWSpace srcIWSp,
  		dstIWSp;
  WlzGreyWSpace srcGWSp,
  		dstGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Gather information about the source object. */
  if(((srcDom = srcObj->domain).core->type != WLZ_INTERVALDOMAIN_INTVL) &&
     (srcDom.core->type != WLZ_INTERVALDOMAIN_RECT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    srcGType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcGType)
    {
      /* Promote grey type for destination object. */
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT:
        dstGType = WLZ_GREY_SHORT;
	break;
      case WLZ_GREY_INT:
	dstGType = WLZ_GREY_INT;
        break;
      case WLZ_GREY_FLOAT:
        dstGType = WLZ_GREY_FLOAT;
	break;
      case WLZ_GREY_DOUBLE:
        dstGType = WLZ_GREY_DOUBLE;
	break;
      case WLZ_GREY_RGBA: /* RGBA to be done RAB */
      default:
        errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vType = WlzGreyTableTypeToTableType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bgdPix = WlzGetBackground(srcObj, &errNum);
  }
  /* Make destination object with it's own values but a shared domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstVal.v = WlzNewValueTb(srcObj, WlzGreyTableType(vType, dstGType, NULL),
    			     bgdPix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj= WlzMakeMain(srcObj->type, srcDom, dstVal, srcObj->plist,
			NULL, &errNum);
  }
  /* Make the buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtX = srcDom.i->lastkl - srcDom.i->kol1 + 1;
    bufSz.vtY = 3;
    if(AlcBit2Malloc(&itvBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bufType = WLZ_GREY_DOUBLE;
    if((AlcDouble2Malloc((double ***)&srcBuf,
			 bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
       (AlcDouble2Malloc((double ***)&wrkBuf,
			 bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Work down and then back up through the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    idD = 0;
    itvBufWidth = (bufSz.vtX + 7) / 8;
    while((errNum == WLZ_ERR_NONE) && (idD < 2))
    {
      /* Initialise buffers. */
      for(idN = 0; idN < 3; ++idN)
      {
	WlzValueSetUByte(*(itvBuf + idN), 0, itvBufWidth);
      }
      rasDir = (idD == 0)? WLZ_RASTERDIR_ILIC: WLZ_RASTERDIR_DLIC;
      if(((errNum = WlzInitGreyRasterScan(srcObj, &srcIWSp, &srcGWSp, rasDir,
					  0)) == WLZ_ERR_NONE) &&
	 ((errNum = WlzInitGreyRasterScan(dstObj, &dstIWSp, &dstGWSp, rasDir,
					  0)) == WLZ_ERR_NONE))
      {
	while((errNum == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&srcIWSp)) == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&dstIWSp)) == WLZ_ERR_NONE))
	{
	  itvLen = srcIWSp.rgtpos - srcIWSp.lftpos + 1;
	  bufPos.vtX = srcIWSp.lftpos - srcDom.i->kol1;
	  /* Process any lines between this and the last. */
	  if(srcIWSp.nwlpos > 0)
	  {
	    idN = (srcIWSp.nwlpos >= 3)? 3: srcIWSp.nwlpos;
	    while(--idN >= 0)
	    {
	      bufPos.vtY = srcIWSp.linpos - srcDom.i->line1 - idN;
	      bufLnIdx = (3 + bufPos.vtY) % 3;
	      WlzValueSetUByte(*(itvBuf + bufLnIdx), 0, itvBufWidth);
	    }
	    dstLnIdx = (bufPos.vtY + 3 + ((idD)? 2: - 2)) % 3;
	  }
	  /* Copy interval to buffer. */
	  srcBufGP.dbp = *((double **)srcBuf + bufLnIdx);
	  wrkBufGP.dbp = *((double **)wrkBuf + bufLnIdx);
	  dstBufGP.dbp = *((double **)wrkBuf + dstLnIdx);
	  WlzBitLnSetItv(*(itvBuf + bufLnIdx),
	  		 bufPos.vtX, bufPos.vtX + itvLen - 1, bufSz.vtX);
	  WlzValueCopyGreyToGrey(srcBufGP, bufPos.vtX, bufType,
				 srcGWSp.u_grintptr, 0, srcGWSp.pixeltype,
				 itvLen);
	  if(rasDir == WLZ_RASTERDIR_DLIC)
	  {
	    WlzValueCopyGreyToGrey(wrkBufGP, bufPos.vtX, bufType,
	    			   dstGWSp.u_grintptr, 0, dstGWSp.pixeltype,
				   itvLen);
	  }
	  /* Apply filter to this interval. */
	  WlzRsvFilterFilterBufYF(ftr, (double **)wrkBuf, (double **)srcBuf,
				itvBuf, bufPos, itvLen, idD);
	  /* Clamp data buffer into the dst interval. */
	  if(rasDir == WLZ_RASTERDIR_ILIC)
	  {
	    WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0, dstGWSp.pixeltype,
				      wrkBufGP, bufPos.vtX, bufType, itvLen);
	  }
	  else
	  {
	    WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0, dstGWSp.pixeltype,
				      dstBufGP, bufPos.vtX, bufType, itvLen);
	  }
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
      ++idD;
    }
  }
  if(itvBuf)
  {
    Alc2Free((void **)itvBuf);
  }
  if(srcBuf)
  {
    Alc2Free(srcBuf);
  }
  if(wrkBuf)
  {
    Alc2Free(wrkBuf);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      (void )WlzFreeObj(dstObj);
    }
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	The filtered object, or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Applies a recursive filter along the lines and/or
*               through the columns of the given 3D domain
*               object with grey values using either double
*               precision floating point arithmetic or fixed
*               point arithmetic.
*               It is assumed that the object type has already been
*               checked, the domain and values are non-null.
* \param	srcObj			Given object.
* \param	ftr			Recursive filter.
* \param	actionMsk		Action mask.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzRsvFilterObj3DXY(WlzObject *srcObj, WlzRsvFilter *ftr,
			              int actionMsk, WlzErrorNum *dstErr)
{
  int		nPlanes;
  WlzObject	*srcObj2D,
  		*dstObj2D,
		*dstObj = NULL;
  WlzDomain 	*srcDom2D;
  WlzValues	*srcVal2D,
  		*dstVal2D;
  WlzDomain	srcDom;
  WlzValues	dstVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstVal.core = NULL;
  srcDom = srcObj->domain;
  srcDom2D = srcDom.p->domains;
  srcVal2D = srcObj->values.vox->values;
  nPlanes = srcDom.p->lastpl - srcDom.p->plane1 + 1;
  dstVal.vox = WlzMakeVoxelValueTb(srcObj->values.vox->type,
				   srcDom.p->plane1, srcDom.p->lastpl,
				   WlzGetBackground(srcObj, NULL),
				   NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dstVal2D = dstVal.vox->values;
    dstObj = WlzMakeMain(srcObj->type, srcDom, dstVal, NULL, NULL, &errNum);
  }
  while((errNum == WLZ_ERR_NONE) &&
        (nPlanes-- > 0))
  {
    srcObj2D = NULL;
    dstObj2D = NULL;
    if(srcDom2D->core)
    {
      dstObj2D = NULL;
      srcObj2D = WlzAssignObject(
      		 WlzMakeMain(WLZ_2D_DOMAINOBJ, *srcDom2D++, *srcVal2D++,
			     NULL, NULL, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	dstObj2D = WlzAssignObject(
		   WlzRsvFilterObj(srcObj2D, ftr, actionMsk, &errNum), NULL);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	*dstVal2D = WlzAssignValues(dstObj2D->values, NULL);
      }
      if(srcObj2D)
      {
        WlzFreeObj(srcObj2D);
      }
      if(dstObj2D)
      {
        WlzFreeObj(dstObj2D);
      }
    }
    else
    {
      dstVal2D->core = NULL;
    }
    ++dstVal2D;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
    }
    else if(dstVal.core)
    {
      WlzFreeVoxelValueTb(dstVal.vox);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	The filtered object, or NULL on error.
* \ingroup	WlzValueFilters
* \brief	Applies a recursive filter through the planes of the
*               given 3D domain object with grey values using either
*               double precision floating point arithmetic or fixed
*               point arithmetic.
*               It is assumed that the object type has already been
*               checked, the domain and values are non-null.
* \param	srcObj			Given object.
* \param	ftr			Recursive filter.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzRsvFilterObj3DZ(WlzObject *srcObj, WlzRsvFilter *ftr,
				     WlzErrorNum *dstErr)
{
  int		tI0,
  		idD,
		idN,
  		idP,
		itvLen,
		dstPnIdx,
		bufPlIdx,
  		nPlanes,
		itvBufArea;
  WlzIVertex3	bufPos,
  		bufSz;
  WlzDomain	srcDom;
  WlzValues	tVal0,
  		srcVal,
		dstVal;
  WlzObjectType	dstValTbType2D;
  WlzGreyType	tmpGType,
		bufType,
  		srcGType,
  		dstGType;
  WlzGreyP	dstBufGP,
  		srcBufGP,
  		wrkBufGP;
  void		***srcBuf = NULL,
  		***wrkBuf = NULL;
  void		**srcBuf2D,
		**dstBuf2D,
  		**wrkBuf2D;
  UBYTE		***itvBuf = NULL;
  UBYTE		**itvBuf2D;
  WlzDomain	*srcDom2D;
  WlzValues	*srcVal2D,
  		*dstVal2D;
  WlzObject	*srcObj2D,
  		*dstObj2D,
		*dstObj = NULL;
  WlzPixelV	bgdPix;
  WlzIntervalWSpace srcIWSp,
  		dstIWSp;
  WlzGreyWSpace srcGWSp,
  		dstGWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Gather information about the source object. */
  if((srcDom = srcObj->domain).core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((srcVal = srcObj->values).core->type != WLZ_VOXELVALUETABLE_GREY)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bgdPix = WlzGetBackground(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find grey type, checking all planes have same grey type. */
    srcDom2D = srcDom.p->domains;
    srcVal2D = srcVal.vox->values;
    nPlanes = srcDom.p->lastpl - srcDom.p->plane1 + 1;
    if(nPlanes > 0)
    {
      idP = 0;
      srcGType = WLZ_GREY_ERROR;
      while((errNum == WLZ_ERR_NONE) && (idP < nPlanes))
      {
	if(srcDom2D->core && srcVal2D->core)
	{
	  tmpGType = WlzGreyTableTypeToGreyType(srcVal2D->core->type, NULL);
	  if(srcGType == WLZ_GREY_ERROR)
	  {
	    srcGType = tmpGType;
	  }
	  else if(srcGType != tmpGType)
	  {
	    errNum = WLZ_ERR_GREY_TYPE;
	  }
	}
	++srcDom2D;
	++srcVal2D;
        ++idP;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nPlanes < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if(nPlanes == 0)
    {
      dstObj = WlzCopyObject(srcObj, &errNum);
    }
    else
    {
      /* Promote grey type for destination object. */
      switch(srcGType)
      {
	case WLZ_GREY_UBYTE: /* FALLTHROUGH */
	case WLZ_GREY_SHORT:
	  dstGType = WLZ_GREY_SHORT;
	  break;
	case WLZ_GREY_INT:
	  dstGType = WLZ_GREY_INT;
	  break;
	case WLZ_GREY_FLOAT:
	  dstGType = WLZ_GREY_FLOAT;
	  break;
	case WLZ_GREY_DOUBLE:
	  dstGType = WLZ_GREY_DOUBLE;
	  break;
        case WLZ_GREY_RGBA: /* RGBA to be done RAB */
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        dstValTbType2D = WlzGreyTableType(WLZ_GREY_TAB_RAGR, dstGType,
					  &errNum);

      }
      /* Make buffers. */
      if(errNum == WLZ_ERR_NONE)
      {
	bufSz.vtX = srcDom.p->lastkl - srcDom.p->kol1 + 1;
	bufSz.vtY = srcDom.p->lastln - srcDom.p->line1 + 1;
	bufSz.vtZ = 3;
	tI0 = (bufSz.vtX + 7) / 8;
	itvBufArea = tI0 * bufSz.vtY;
	if(AlcBit3Malloc(&itvBuf, bufSz.vtZ, bufSz.vtY,
			 bufSz.vtX) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	bufType = WLZ_GREY_DOUBLE;
	if((AlcDouble3Malloc((double ****)&srcBuf,
			     bufSz.vtZ, bufSz.vtY,
			     bufSz.vtX) != ALC_ER_NONE) ||
	   (AlcDouble3Malloc((double ****)&wrkBuf,
			     bufSz.vtZ, bufSz.vtY,
			     bufSz.vtX) != ALC_ER_NONE))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      /* Make destination object with it's own voxel value table but with a
       * shared domain. The 2D values are created/added later. */
      if(errNum == WLZ_ERR_NONE)
      {
	dstVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					 srcDom.p->plane1, srcDom.p->lastpl,
					 bgdPix, NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	dstObj= WlzMakeMain(srcObj->type, srcDom, dstVal, srcObj->plist,
			    NULL, &errNum);
      }
      /* Work down and then back up through the object's planes. */
      if(errNum == WLZ_ERR_NONE)
      {
	idD = 0;
	while((errNum == WLZ_ERR_NONE) && (idD < 2))
	{
	  /* Initialise buffers. */
	  for(idN = 0; idN < 3; ++idN)
	  {
	    WlzValueSetUByte(**(itvBuf + idN), 0, itvBufArea);
	  }
	  idP = 0;
	  bufPos.vtZ = (idD == 0)? 0: nPlanes - 1;
	  while((errNum == WLZ_ERR_NONE) && (idP < nPlanes))
	  {
	    bufPlIdx = (bufPos.vtZ + 3 + 0) % 3;
	    dstPnIdx = (bufPos.vtZ + 3 + 2) % 3;
	    srcDom2D = srcDom.p->domains + bufPos.vtZ;
	    srcVal2D = srcVal.vox->values + bufPos.vtZ;
	    dstVal2D = dstVal.vox->values + bufPos.vtZ;
	    itvBuf2D = *(itvBuf + bufPlIdx);
	    srcBuf2D = *(srcBuf + bufPlIdx); 
	    wrkBuf2D = *(wrkBuf + bufPlIdx);
	    dstBuf2D = *(wrkBuf + dstPnIdx);
	    /* Clear this plane's interval buffer bit mask. */
	    WlzValueSetUByte(*itvBuf2D, 0, itvBufArea);
	    /* Process non-empty planes. */
	    if(srcDom2D->core)
	    {
	      /* Make a 2D object from current source plane. */
	      srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *srcDom2D, *srcVal2D,
				     NULL, NULL, &errNum);
	      if((errNum == WLZ_ERR_NONE) && (idD == 0))
	      {
		/* If first pass through object (going down through planes)
		 * then make a new 2D value table. */
		tVal0.v  = WlzNewValueTb(srcObj2D, dstValTbType2D,
					 bgdPix, &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  *dstVal2D = WlzAssignValues(tVal0, NULL);
		}
	      }
	      /* Make a 2D object from destination plane. */
	      dstObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *srcDom2D, *dstVal2D,
				     NULL, NULL, &errNum);
	      /* Process each interval of this plane. */
	      if(((errNum = WlzInitGreyScan(srcObj2D, &srcIWSp,
	      				    &srcGWSp)) == WLZ_ERR_NONE) &&
	         ((errNum = WlzInitGreyScan(dstObj2D, &dstIWSp,
	      				    &dstGWSp)) == WLZ_ERR_NONE))
	      {
		while((errNum == WLZ_ERR_NONE) &&
		      ((errNum = WlzNextGreyInterval(
		      			&srcIWSp)) == WLZ_ERR_NONE) &&
		      ((errNum = WlzNextGreyInterval(
		      			&dstIWSp)) == WLZ_ERR_NONE))
		{
		  itvLen = srcIWSp.rgtpos - srcIWSp.lftpos + 1;
		  bufPos.vtX = srcIWSp.lftpos - srcDom.p->kol1;
		  bufPos.vtY = srcIWSp.linpos - srcDom.p->line1;
		  /* Copy interval to buffer. */
		  srcBufGP.dbp = *((double **)srcBuf2D + bufPos.vtY);
		  wrkBufGP.dbp = *((double **)wrkBuf2D + bufPos.vtY);
		  dstBufGP.dbp = *((double **)dstBuf2D + bufPos.vtY);
		  WlzBitLnSetItv(*(itvBuf2D + bufPos.vtY),
		  		 bufPos.vtX, bufPos.vtX + itvLen - 1,
				 bufSz.vtX);
		  WlzValueCopyGreyToGrey(srcBufGP, bufPos.vtX, bufType,
		  			 srcGWSp.u_grintptr, 0,
					 srcGWSp.pixeltype,
					 itvLen);
		  if(idD)
		  {
		    WlzValueCopyGreyToGrey(wrkBufGP, bufPos.vtX, bufType,
					   dstGWSp.u_grintptr, 0,
					   dstGWSp.pixeltype,
					   itvLen);
		  }
		  /* Apply filter to this interval. */
		  WlzRsvFilterFilterBufZF(ftr, (double ***)wrkBuf,
					  (double ***)srcBuf, itvBuf,
					  bufPos, itvLen, idD);
		  /* Clamp data buffer back into the destination plane. */
		  if(idD == 0)
		  {
		    WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0,
					      dstGWSp.pixeltype,
					      wrkBufGP, bufPos.vtX, bufType,
					      itvLen);
		  }
		  else
		  {
		    WlzValueClampGreyIntoGrey(dstGWSp.u_grintptr, 0,
					      dstGWSp.pixeltype,
					      dstBufGP, bufPos.vtX, bufType,
					      itvLen);
		  }
		}
		if(errNum == WLZ_ERR_EOO)
		{
		  errNum = WLZ_ERR_NONE;
		}
	      }
	      if(srcObj2D)
	      {
	        WlzFreeObj(srcObj2D);
	      }
	      if(dstObj2D)
	      {
	        WlzFreeObj(dstObj2D);
	      }
	    }
	    ++idP;
	    bufPos.vtZ -= (idD * 2) - 1; /* ++ for idD == 0, -- for idD == 1 */
	  }
	  ++idD;
	}
      }
    }
  }
  if(itvBuf)
  {
    Alc3Free((void ***)itvBuf);
  }
  if(srcBuf)
  {
    Alc3Free(srcBuf);
  }
  if(wrkBuf)
  {
    Alc3Free(wrkBuf);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      (void )WlzFreeObj(dstObj);
    }
    else if(dstVal.core)
    {
      WlzFreeVoxelValueTb(dstVal.vox);
    }
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

#ifdef WLZ_RSVFILTER_TEST_1D
int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
  		usage = 0,
		order = 1,
  		idN;
  double	param = 1.0,
  		sum;
  WlzRsvFilterName name = WLZ_RSVFILTER_NAME_DERICHE_1;
  WlzRsvFilter *ftr = NULL;
  static double	data[100],
  		buf0[100],
		buf1[100];
  static char  optList[] = "dgp:r:";

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	name = WLZ_RSVFILTER_NAME_DERICHE_1;
	break;
      case 'g':
	name = WLZ_RSVFILTER_NAME_GAUSS_1;
	break;
      case 'p':
	if(sscanf(optarg, "%lg", &param) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'r':
	if((sscanf(optarg, "%d", &order) != 1) || (order < 0) || (order > 2))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    switch(order)
    {
      case 0:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_0: WLZ_RSVFILTER_NAME_GAUSS_0;
        break;
      case 1:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_1: WLZ_RSVFILTER_NAME_GAUSS_1;
        break;
      case 2:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_2: WLZ_RSVFILTER_NAME_GAUSS_2;
        break;
    }
    for(idN = 0; idN < 100; ++idN)
    {
      data[idN] = 0.0;
    }
    for(idN = 50 - 10; idN < 50 + 10; ++idN)
    {
      data[idN] = 100.0;
    }
    for(idN = 0; idN < 100; ++idN)
    {
      printf("%3d %g\n", idN, data[idN]); 
    }
    printf("\n"); 
    ftr = WlzRsvFilterMakeFilter(name, param, NULL);
    WlzRsvFilterFilterBufXF(ftr, data, buf0, buf1, 100);
    sum = 0.0;
    for(idN = 0; idN < 100; ++idN)
    {
      printf("%3d %g\n", idN, data[idN]); 
      sum += data[idN];
    }
    fprintf(stderr, "%g %g %g\n", param, sum, ftr->c);
    WlzRsvFilterFreeFilter(ftr);
  }
  else if(usage)
  {
    ok = 1;
    fprintf(stderr, "Usage: %s [-d] [-g] [-p#] [-r#]\n"
    	    "Test for woolz deriche filter.\n"
	    "Options are:\n"
	    "  -d  Deriche filter\n"
	    "  -g  Gaussian IIR filter\n"
	    "  -p  Parameter (alpha for Deriche, sigma for Gaussian)\n"
	    "  -r  Differential [0-2]\n",
	    argv[0]);

  }
  return(ok? 0: 1);
}
#endif /* WLZ_RSVFILTER_TEST_1D */

#if defined(WLZ_RSVFILTER_TEST_2D) || defined(WLZ_RSVFILTER_TEST_3D)

int		main(int argc, char *argv[])
{
  int		ok = 1,
#ifdef WLZ_RSVFILTER_TEST_2D
  		dir = WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y,
#else /* ! WLZ_RSVFILTER_TEST_2D */
#ifdef WLZ_RSVFILTER_TEST_3D
  		dir = WLZ_RSVFILTER_ACTION_X | WLZ_RSVFILTER_ACTION_Y |
		      WLZ_RSVFILTER_ACTION_Z,
#endif /* WLZ_RSVFILTER_TEST_3D */
#endif /* WLZ_RSVFILTER_TEST_2D */
  		option,
  		usage = 0,
		order = 1,
  		idN;
  double	param = 1.0;
  WlzRsvFilterName name = WLZ_RSVFILTER_NAME_DERICHE_1;
  WlzRsvFilter *dFtr = NULL;
  static char  optList[] = "dgnp:r:";
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	name = WLZ_RSVFILTER_NAME_DERICHE_1;
	break;
      case 'g':
	name = WLZ_RSVFILTER_NAME_GAUSS_1;
	break;
      case 'n':
        dir = WLZ_RSVFILTER_ACTION_NONE;
	break;
      case 'p':
	if(sscanf(optarg, "%lg", &param) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'r':
	if((sscanf(optarg, "%d", &order) != 1) || (order < 0) || (order > 2))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    switch(order)
    {
      case 0:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_0: WLZ_RSVFILTER_NAME_GAUSS_0;
        break;
      case 1:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_1: WLZ_RSVFILTER_NAME_GAUSS_1;
        break;
      case 2:
        name = (name == WLZ_RSVFILTER_NAME_DERICHE_1)?
		WLZ_RSVFILTER_NAME_DERICHE_2: WLZ_RSVFILTER_NAME_GAUSS_2;
        break;
    }
    if(((dFtr = WlzRsvFilterMakeFilter(name, param, &errNum)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      fprintf(stderr, "%s Failed to create filter (%s)\n",
      	      argv[0], WlzStringFromErrorNum(errNum, NULL));
    }
    if(ok)
    {
      if((inObj = WlzAssignObject(WlzReadObj(stdin, &errNum), NULL)) == NULL)
      {
        ok = 0;
	fprintf(stderr, "%s Failed to read input object (%s)\n",
		argv[0], WlzStringFromErrorNum(errNum, NULL));
      }
    }
    if(ok)
    {
      if(((outObj = WlzAssignObject(WlzRsvFilterObj(inObj, dFtr,
				    dir, &errNum), NULL)) == NULL) ||
	 (errNum != WLZ_ERR_NONE))
      {
        ok = 0;
	fprintf(stderr, "%s Failed to filter object (%s)\n",
		argv[0], WlzStringFromErrorNum(errNum, NULL));
      }
    }
    if(ok)
    {
      if(WlzWriteObj(stdout, outObj) != WLZ_ERR_NONE)
      {
        ok = 0;
	fprintf(stderr, "%s Failed to write object (%s)\n",
		argv[0], WlzStringFromErrorNum(errNum, NULL));
      }
    }
    WlzRsvFilterFreeFilter(dFtr);
  }
  else if(usage)
  {
    ok = 1;
    fprintf(stderr, "Usage: %s [-d] [-g] [-n] [-p#] [-r#]\n"
    	    "Test for woolz deriche filter.\n"
	    "Options are:\n"
	    "  -n  No filter\n"
	    "  -d  Deriche filter\n"
	    "  -g  Gaussian filter\n"
	    "  -p  Parameter (alpha for Deriche, sigma for Gaussian)\n"
	    "  -r  Differential [0-2]\n",
	    argv[0]);

  }
  return(ok? 0: 1);
}
#endif /* WLZ_RSVFILTER_TEST_2D || WLZ_RSVFILTER_TEST_3D */
