#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzRsvFilter.c
* Date:         May 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Recursice filters for Woolz, these are also known as
*		Deriche and infinite impulse responce (IIR) filters.
*		The recursive filter code can be tested by defining
*		WLZ_RSVFILTER_TEST_1D, WLZ_RSVFILTER_TEST_2D or
*		WLZ_RSVFILTER_TEST_3D. See code below.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <float.h>
#include <Wlz.h>

/* #define WLZ_RSVFILTER_TEST_1D */
#define WLZ_RSVFILTER_TEST_2D
/* #define WLZ_RSVFILTER_TEST_3D */

#define WLZ_RSVFILTER_IFAC	(1<<10)     /* Fixed point arithmetic factor */


static void	WlzRsvFilterFilterBufXF(WlzRsvFilter *,
				      double *, double *, double *,
				      int);
static void	WlzRsvFilterFilterBufXI(WlzRsvFilter *,
				      int *, int *, int *,
				      int);
static void	WlzRsvFilterFilterBufYF(WlzRsvFilter *,
				      double **, double **, UBYTE **,
				      WlzIVertex2, int, int);
static void	WlzRsvFilterFilterBufYI(WlzRsvFilter *,
				      int **, int **, UBYTE **,
				      WlzIVertex2, int, int);

/************************************************************************
* Function:	WlzRsvFilterFreeFilter
* Returns:	void
* Purpose:	Free a recursive filter.
* Global refs:	-
* Parameters:	WlzRsvFilter *ftr:	Filter to free.
************************************************************************/
void		WlzRsvFilterFreeFilter(WlzRsvFilter *ftr)
{
  if(ftr == NULL)
  {
    AlcFree(ftr);
  }
}

/************************************************************************
* Function:	WlzRsvFilterMakeFilter
* Returns:	WlzRsvFilter *:	New recursive filter.
* Purpose:	Creates a new recursive filter.
* Global refs:	-
* Parameters:	WlzRsvFilterName name:  Name of filter.
*		double prm:		Parameter is alpha for Deriche's
*					operators and sigma for
*					Gaussians, must be >= 0.0.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
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

/************************************************************************
* Function:	WlzRsvFilterObj2DX
* Returns:	WlzObject:		The filtered object, or NULL on
*					error.
* Purpose:	Applies a recursive filter along the lines of the given
*		2D domain object with grey values using either double
*		precision floating point arithmetic or fixed
*		point arithmetic.
*		It is assumed that the object type has already been
*		checked, the domain and values are non-null.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given 2D domain object.
*		WlzRsvFilter *ftr:	Recursive filter.
*		int useFP:		Use floating point if non-zero.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzObject *WlzRsvFilterObj2DX(WlzObject *srcObj, WlzRsvFilter *ftr,
				   int useFP, WlzErrorNum *dstErr)
{
  int		bufSz,
  		bufItemSz,
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
			srcObj, &errNum);
  }
  /* Make the working buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz = srcDom.i->lastkl - srcDom.i->kol1 + 1;
    bufItemSz = useFP? sizeof(double): sizeof(int);
    bufSpace = bufItemSz * bufSz;
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
    if(useFP)
    {
      bufType = WLZ_GREY_DOUBLE;
      bufGP.dbp = (double *)datBuf;
    }
    else
    {
      bufType = WLZ_GREY_INT;
      bufGP.inp = (int *)datBuf;
    }
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
	if(useFP)
	{
	  WlzRsvFilterFilterBufXF(ftr, (double *)datBuf,
	  		        (double *)lnBuf0, (double *)lnBuf1,
				itvLen);
	}
	else
	{
	  WlzRsvFilterFilterBufXI(ftr, (int *)datBuf,
	  			(int *)lnBuf0, (int *)lnBuf1,
				itvLen);
	}
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

/************************************************************************
* Function:	WlzRsvFilterFilterBufXF
* Returns:	void
* Purpose:	Filters the given single line data buffer using an
*		IIR filter defined by the filter coefficients and
*		double precision arithmetic.
*		Before optimization this and all the other IIR filter
*		code looks like:
*		  for(i = 2; i < dataSz; ++i)
*		  {
*		    buf0[i] =   (a0 * data[i + 0]) + (a1 * data[i - 1])
*			      - (b0 * buf0[i - 1]) - (b1 * buf0[i - 2]);
*		  }
*		  for(i = dataSz - 3; i >= 0; --i)
*		  {
*		    buf1[i] =   (a2 * data[i + 1]) + (a3 * data[i + 2])
*			      - (b0 * buf1[i + 1]) - (b1 * buf1[i + 2]);
*		  }
*		  for(i = 0; i < dataSz; ++i)
*		  {
*		    data[i] =  c * (buf0[i] + buf1[i]);
*		  }
* Global refs:	-
* Parameters:	WlzRsvFilter *ftr:	The filter.
*		double *data:		The line of data to be filtered.
*		double *buf0:		Buffer with at least dataSz
*					elements.
*		double *buf1:		Buffer with at least dataSz
*					elements.
*		int dataSz:		Number of data.
************************************************************************/
static void	WlzRsvFilterFilterBufXF(WlzRsvFilter *ftr, double *data,
				      double *buf0, double *buf1,
				      int dataSz)
{
  int		cnt;
  double	a0,
  		a1,
		b0,
		b1,
		c,
		d0,
		d1,
		f1,
		f2;
  double	*dP,
  		*fP0,
		*fP1;

  if(dataSz > 2)
  {
    a0 = ftr->a[0];
    a1 = ftr->a[1];
    b0 = ftr->b[0];
    b1 = ftr->b[1];
    dP = data;
    fP0 = buf0;
    d1 = d0 = *dP;
    f2 = (a0 + a1) * d0;
    f1 = f2 * (1.0 - b0);
    cnt = dataSz;
    while(cnt-- > 0)
    {
      *fP0 = (a0 * d0) + (a1 * d1) - (b0 * f1) - (b1 * f2);
      d1 = d0;
      d0 = *++dP;
      f2 = f1;
      f1 = *fP0++;
    }
    a0 = ftr->a[2];
    a1 = ftr->a[3];
    dP = data + dataSz - 1;
    fP1 = buf1 + dataSz - 1;
    d1 = d0 = *dP;
    f1 = (a0 + a1) * d0;
    *fP1 = f1 * (1.0 - b0);
    cnt = dataSz;
    while(--cnt > 0)
    {
      f2 = f1;
      f1 = *fP1--;
      *fP1 = (a0 * d0) + (a1 * d1) - (b0 * f1) - (b1 * f2);
      d1 = d0;
      d0 = *--dP;
    }
    c = ftr->c;
    dP = data;
    fP0 = buf0;
    fP1 = buf1;
    cnt = dataSz; 
    while(cnt-- > 0)
    {
      *dP++ = c * (*fP0++ + *fP1++);
    }
  }
}

/************************************************************************
* Function:	WlzRsvFilterFilterBufXI
* Returns:	void
* Purpose:	Filters the given single line data buffer using an
*		IIR filter defined by the filter coefficients and
*		fixed point arithmetic.
*		See the double precision function for simple code.
* Global refs:	-
* Parameters:	WlzRsvFilter *ftr:	The filter.
*		int *data:		The line of data to be filtered.
*		int *buf0:		Buffer with at least dataSz
*					elements.
*		int *buf1:		Buffer with at least dataSz
*					elements.
*		int dataSz:		Number of data.
************************************************************************/
static void	WlzRsvFilterFilterBufXI(WlzRsvFilter *ftr, int *data,
				      int *buf0, int *buf1,
				      int dataSz)
{

  int		cnt;
  int		a0,
  		a1,
		b0,
		b1,
		c,
		d0,
		d1,
		f1,
		f2;
  int		*dP,
  		*fP0,
		*fP1;

  if(dataSz > 2)
  {
    a0 = ftr->a[0] * WLZ_RSVFILTER_IFAC;
    a1 = ftr->a[1] * WLZ_RSVFILTER_IFAC;
    b0 = ftr->b[0] * WLZ_RSVFILTER_IFAC;
    b1 = ftr->b[1] * WLZ_RSVFILTER_IFAC;
    dP = data;
    fP0 = buf0;
    d1 = d0 = *dP;
    f2 = ((a0 + a1) * d0) / WLZ_RSVFILTER_IFAC;
    f1 = f2  - ((b0 * f2) / WLZ_RSVFILTER_IFAC);
    cnt = dataSz;
    while(cnt-- > 0)
    {
      *fP0 = ((a0 * d0) / WLZ_RSVFILTER_IFAC) +
             ((a1 * d1) / WLZ_RSVFILTER_IFAC) -
	     ((b0 * f1) / WLZ_RSVFILTER_IFAC) -
	     ((b1 * f2) / WLZ_RSVFILTER_IFAC);
      d1 = d0;
      d0 = *++dP;
      f2 = f1;
      f1 = *fP0++;
    }
    a0 = ftr->a[2] * WLZ_RSVFILTER_IFAC;
    a1 = ftr->a[3] * WLZ_RSVFILTER_IFAC;
    dP = data + dataSz - 1;
    fP1 = buf1 + dataSz - 1;
    d1 = d0 = *dP;
    f1 = ((a0 + a1) * d0) / WLZ_RSVFILTER_IFAC;
    *fP1 = (f1 * (WLZ_RSVFILTER_IFAC - b0)) / WLZ_RSVFILTER_IFAC;
    cnt = dataSz;
    while(--cnt > 0)
    {
      f2 = f1;
      f1 = *fP1--;
      *fP1 = ((a0 * d0) / WLZ_RSVFILTER_IFAC) +
      	     ((a1 * d1) / WLZ_RSVFILTER_IFAC) -
	     ((b0 * f1) / WLZ_RSVFILTER_IFAC) -
	     ((b1 * f2) / WLZ_RSVFILTER_IFAC);
      d1 = d0;
      d0 = *--dP;
    }
    c = ftr->c * WLZ_RSVFILTER_IFAC;
    dP = data;
    fP0 = buf0;
    fP1 = buf1;
    cnt = dataSz; 
    while(cnt-- > 0)
    {
      *dP++ = (c * (*fP0++ + *fP1++)) / WLZ_RSVFILTER_IFAC;
    }
  }
}

/************************************************************************
* Function:	WlzRsvFilterObj2DY
* Returns:	WlzObject:		The filtered object, or NULL on
*					error.
* Purpose:	Applies a recursive filter through the columns of the
*		given 2D domain object with grey values using either
*		double precision floating point arithmetic or fixed
*		point arithmetic.
*		It is assumed that the object type has already been
*		checked, the domain and values are non-null.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given 2D domain object.
*		WlzRsvFilter *ftr:	Recursive filter.
*		int useFP:		Use floating point if non-zero.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
static WlzObject *WlzRsvFilterObj2DY(WlzObject *srcObj, WlzRsvFilter *ftr,
				   int useFP, WlzErrorNum *dstErr)
{
  int		idD,
  		idN,
		bufLnIdx,
  		itvLen;
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
			srcObj, &errNum);
  }
  /* Make the buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtX = srcDom.i->lastkl - srcDom.i->kol1 + 1;
    bufSz.vtY = 3;
    if(AlcUnchar2Malloc(&itvBuf, bufSz.vtY, bufSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(useFP)
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
    else
    {
      bufType = WLZ_GREY_INT;
      if((AlcInt2Malloc((int ***)&srcBuf,
      	       		bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
         (AlcInt2Malloc((int ***)&wrkBuf,
	 	        bufSz.vtY, bufSz.vtX) != ALC_ER_NONE))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  /* Work down and then back up through the object. */
  if(errNum == WLZ_ERR_NONE)
  {
    idD = 0;
    while((errNum == WLZ_ERR_NONE) && (idD < 2))
    {
      /* Initialise buffers. */
      for(idN = 0; idN < 3; ++idN)
      {
	WlzValueSetUByte(*(itvBuf + idN), 0, bufSz.vtX);
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
	      WlzValueSetUByte(*(itvBuf + bufLnIdx), 0, bufSz.vtX);
	    }
	  }
	  /* Copy interval to buffer. */
	  if(useFP)
	  {
	    srcBufGP.dbp = *((double **)srcBuf + bufLnIdx);
	    wrkBufGP.dbp = *((double **)wrkBuf + bufLnIdx);
	  }
	  else
	  {
	    srcBufGP.inp = *((int **)srcBuf + bufLnIdx);
	    wrkBufGP.inp = *((int **)wrkBuf + bufLnIdx);
	  }
	  WlzValueSetUByte(*(itvBuf + bufLnIdx) + bufPos.vtX, 1, itvLen);
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
	  if(useFP)
	  {
	    wrkBufGP.dbp = *((double **)wrkBuf + bufLnIdx);
	    dstBufGP.inp = *((int **)wrkBuf + ((bufPos.vtY + 2) % 3));
	    WlzRsvFilterFilterBufYF(ftr, (double **)wrkBuf, (double **)srcBuf,
				  itvBuf, bufPos, itvLen,
				  idD);
	  }
	  else
	  {
	    wrkBufGP.inp = *((int **)wrkBuf + bufLnIdx);
	    dstBufGP.inp = *((int **)wrkBuf + ((bufPos.vtY + 2) % 3));
	    WlzRsvFilterFilterBufYI(ftr, (int **)wrkBuf, (int **)srcBuf,
				  itvBuf, bufPos, itvLen,
				  idD);
	  }
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

/************************************************************************
* Function:	WlzRsvFilterFilterBufYF
* Returns:	void
* Purpose:	Filters a single interval using a vertical IIR filter
*		defined by the filter coefficients and double precision
*		arithmetic.
*		See the double precision horizontal function for
*		simple code.
* Global refs:	-
* Parameters:	WlzRsvFilter *ftr:	The filter.
*		double **wrkBuf:	Working buffer.
*		double **srcBuf:	Source buffer.
*		UBYTE **itvBuf:		Within interval buffer.
*		WlzIVertex2 bufPos:	Position within the buffer.
*		int itvLen:		Interval length.
*		int goingUp:		Non-zero if going up through
*					the lines.
************************************************************************/
static void	WlzRsvFilterFilterBufYF(WlzRsvFilter *ftr,
				      double **wrkBuf, double **srcBuf,
				      UBYTE **itvBuf,
				      WlzIVertex2 bufPos, int itvLen,
				      int goingUp)
{
  int		cnt,
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
		tD0;
  UBYTE		*iP1,
  		*iP2;
  double	*dP0,
  		*dP1,
		*dP2,
		*fP0,
		*fP1,
		*fP2;

  if(itvLen > 1)
  {
    if(goingUp)
    {
      idL0 = (bufPos.vtY + 0) % 3;
      idL1 = (bufPos.vtY + 1) % 3;
      idL2 = (bufPos.vtY + 2) % 3;
      a0 = ftr->a[2];
      a1 = ftr->a[3];
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
    iP1 = *(itvBuf + idL1) + bufPos.vtX;
    iP2 = *(itvBuf + idL2) + bufPos.vtX;
    dP0 = *(srcBuf + idL0) + bufPos.vtX;
    fP0 = *(wrkBuf + idL0) + bufPos.vtX;
    dP1 = *(srcBuf + idL1) + bufPos.vtX;
    fP1 = *(wrkBuf + idL1) + bufPos.vtX;
    dP2 = *(srcBuf + idL2) + bufPos.vtX;
    fP2 = *(wrkBuf + idL2) + bufPos.vtX;
    cnt = itvLen;
    while(cnt-- > 0)
    {
      if(*iP1 == 0)
      {
	if(goingUp)
	{
	  *dP2 = *dP0;
	}
	*dP1 = *dP0;
	*fP1 = (a0 + a1) * *dP0;
	*fP2 = *fP1 * (1.0 - b0);
      }
      else
      {
	if(*iP2 == 0)
	{
	  if(goingUp)
	  {
	    *dP2 = *dP1;
	  }
	  *fP2 = *fP1 * (1.0 - b0);
	}
      }
      if(goingUp)
      {

        tD0 = *fP0;
	*fP0 = (a0 * *dP1) + (a1 * *dP2) - (b0 * *fP1) - (b1 * *fP2);
	*fP2 = c * (tD0 + *fP0);
      }
      else
      {
	*fP0 =   (a0 * *dP0) + (a1 * *dP1) - (b0 * *fP1) - (b1 * *fP2);
      }
      ++iP1; ++iP2;
      ++dP0; ++fP0;
      ++dP1; ++fP1;
      ++dP2; ++fP2;
    }
  }
}

/************************************************************************
* Function:	WlzRsvFilterFilterBufYI
* Returns:	void
* Purpose:	Filters a single interval using a vertical IIR filter
*		defined by the filter coefficients and double precision
*		arithmetic.
*		See the double precision horizontal function for
*		simple code.
* Global refs:	-
* Parameters:	WlzRsvFilter *ftr:	The filter.
*		int **wrkBuf:		Working buffer.
*		int **srcBuf:		Source buffer.
*		UBYTE **itvBuf:		Within interval buffer.
*		WlzIVertex2 bufPos:	Position within the buffer.
*		int itvLen:		Interval length.
*		int goingUp:		Non-zero if going up through
*					the lines.
************************************************************************/
static void	WlzRsvFilterFilterBufYI(WlzRsvFilter *ftr,
				      int **wrkBuf, int **srcBuf,
				      UBYTE **itvBuf,
				      WlzIVertex2 bufPos, int itvLen,
				      int goingUp)
{
  int		cnt,
  		idL0,
		idL1,
		idL2;
  int		a0,
  		a1,
		a2,
		a3,
		b0,
		b1,
		c,
		tI0;
  UBYTE		*iP1,
  		*iP2;
  int		*dP0,
  		*dP1,
		*dP2,
		*fP0,
		*fP1,
		*fP2;

  if(itvLen > 1)
  {
    if(goingUp)
    {
      idL0 = (bufPos.vtY + 0) % 3;
      idL1 = (bufPos.vtY + 1) % 3;
      idL2 = (bufPos.vtY + 2) % 3;
      a0 = ftr->a[2] * WLZ_RSVFILTER_IFAC;
      a1 = ftr->a[3] * WLZ_RSVFILTER_IFAC;
    }
    else
    {
      idL0 = (bufPos.vtY + 3 - 0) % 3;
      idL1 = (bufPos.vtY + 3 - 1) % 3;
      idL2 = (bufPos.vtY + 3 - 2) % 3;
      a0 = ftr->a[0] * WLZ_RSVFILTER_IFAC;
      a1 = ftr->a[1] * WLZ_RSVFILTER_IFAC;
    }
    b0 = ftr->b[0] * WLZ_RSVFILTER_IFAC;
    b1 = ftr->b[1] * WLZ_RSVFILTER_IFAC;
    c = ftr->c * WLZ_RSVFILTER_IFAC;
    iP1 = *(itvBuf + idL1) + bufPos.vtX;
    iP2 = *(itvBuf + idL2) + bufPos.vtX;
    dP0 = *(srcBuf + idL0) + bufPos.vtX;
    fP0 = *(wrkBuf + idL0) + bufPos.vtX;
    dP1 = *(srcBuf + idL1) + bufPos.vtX;
    fP1 = *(wrkBuf + idL1) + bufPos.vtX;
    dP2 = *(srcBuf + idL2) + bufPos.vtX;
    fP2 = *(wrkBuf + idL2) + bufPos.vtX;
    cnt = itvLen;
    while(cnt-- > 0)
    {
      if(*iP1 == 0)
      {
	if(goingUp)
	{
	  *dP2 = *dP0;
	}
	*dP1 = *dP0;

	*fP1 = ((a0 + a1) * *dP0) / WLZ_RSVFILTER_IFAC;
	*fP2 = (*fP1 * (1.0 - b0)) / WLZ_RSVFILTER_IFAC;
      }
      else
      {
	if(*iP2 == 0)
	{
	  if(goingUp)
	  {
	    *dP2 = *dP1;
	  }
	  *fP2 = (*fP1 * (1.0 - b0)) / WLZ_RSVFILTER_IFAC;
	}
      }
      if(goingUp)
      {

        tI0 = *fP0;
	*fP0 =   ((a0 * *dP1) / WLZ_RSVFILTER_IFAC)
	       + ((a1 * *dP2) / WLZ_RSVFILTER_IFAC)
	       - ((b0 * *fP1) / WLZ_RSVFILTER_IFAC)
	       - ((b1 * *fP2) / WLZ_RSVFILTER_IFAC);
	*fP2 = (c * (tI0 + *fP0)) / WLZ_RSVFILTER_IFAC;
      }
      else
      {
	*fP0 =   ((a0 * *dP0) / WLZ_RSVFILTER_IFAC)
	       + ((a1 * *dP1) / WLZ_RSVFILTER_IFAC)
	       - ((b0 * *fP1) / WLZ_RSVFILTER_IFAC)
	       - ((b1 * *fP2) / WLZ_RSVFILTER_IFAC);
      }
      ++iP1; ++iP2;
      ++dP0; ++fP0;
      ++dP1; ++fP1;
      ++dP2; ++fP2;
    }
  }
}

/************************************************************************
* Function:	WlzRsvFilterObj
* Returns:	WlzObject:		The filtered object, or NULL on
*					error.
* Purpose:	Applies a recursive filter to the given object.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		WlzRsvFilter *ftr:	Recursive filter.
*		int actionMsk:		Action mask.
*		WlzErrorNum *dstErr:	Destination error pointer, may
*					be null.
************************************************************************/
WlzObject	*WlzRsvFilterObj(WlzObject *srcObj, WlzRsvFilter *ftr,
			         int actionMsk, int useFP,
				 WlzErrorNum *dstErr)
{
  WlzObject	*xObj = NULL,
		*yObj = NULL,
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
	/* Filter in each required direction. */
	if((errNum == WLZ_ERR_NONE) &&
	   ((actionMsk & WLZ_RSVFILTER_ACTION_X) != 0))
	{
	  xObj = WlzRsvFilterObj2DX(srcObj, ftr, useFP, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) &&
	   ((actionMsk & WLZ_RSVFILTER_ACTION_Y) != 0))
	{
	  if(xObj)
	  {
	    yObj = WlzRsvFilterObj2DY(xObj, ftr, useFP, &errNum);
	  }
	  else
	  {
	    yObj = WlzRsvFilterObj2DY(srcObj, ftr, useFP, &errNum);
	  }
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
        break;
      case WLZ_3D_DOMAINOBJ:
	/* TODO No 3D objects yet. */
        WLZ_ERR_OBJECT_TYPE;
        break;
      default:
        WLZ_ERR_OBJECT_TYPE;
	break;
    }
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
		useFP = 1,
  		option,
  		usage = 0,
		order = 1,
  		idN;
  double	param = 1.0;
  WlzRsvFilterName name = WLZ_RSVFILTER_NAME_DERICHE_1;
  WlzRsvFilter *dFtr = NULL;
  static char  optList[] = "dgfinp:r:";
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
      case 'f':
        useFP = 1;
	break;
      case 'i':
        useFP = 0;
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
				    dir, useFP,
				    &errNum), NULL)) == NULL) ||
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
    fprintf(stderr, "Usage: %s [-f] [-i] [-d] [-g] [-n] [-p#] [-r#]\n"
    	    "Test for woolz deriche filter.\n"
	    "Options are:\n"
	    "  -f  Use floating point arithmetic\n"
	    "  -i  Use integer arithmetic\n"
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
