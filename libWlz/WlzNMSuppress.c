#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzNMSuppress_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzNMSuppress.c
* \author       Bill Hill
* \date         May 1999
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
* \brief	A non-maximal supression filter, which constructs a new
* 		domain object using a Canny-like non-maximal suppression
* 		algorithm. The domain is the non-maximally suppressed
* 		domain and the values are the encoded gradient direction.
* \ingroup	WlzFeatures
* \todo		WlzNMSuppress3D() has yet to be implemented.
*/

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

static WlzErrorNum 		WlzNMSuppress2DBufI(
				  WlzIntervalDomain *dstIDom,
				  int **grdMBuf,
				  int *grdYBuf,
				  int *grdXBuf,
				  WlzDynItvPool *iPool,
				  WlzUByte *dstBuf,
				  int dstLen,
				  WlzIVertex2 dstPos,
				  WlzIVertex2 orgPos,
				  int minGM);
static WlzErrorNum 		WlzNMSuppress2DBufD(
				  WlzIntervalDomain *dstIDom,
				  double **grdMBuf,
				  double *grdYBuf,
				  double *grdXBuf,
				  WlzDynItvPool *iPool,
				  WlzUByte *dstBuf,
				  int dstLen,
				  WlzIVertex2 dstPos,
				  WlzIVertex2 orgPos,
				  double minGM);
static WlzObject 		*WlzNMSuppress2D(
				  WlzObject *grdM,
				  WlzObject *grdY,
				  WlzObject *grdX,
				  WlzPixelV minThrV,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzNMSuppress3D(
				  WlzObject *grdM,
				  WlzObject *grdZ,
				  WlzObject *grdY,
				  WlzObject *grdX,
				  WlzPixelV minThrV,
				  WlzErrorNum *dstErr);

/*!
* \return	Woolz error code.
* \ingroup      WlzFeatures
* \brief	Performs non-maximal suppression on the given buffers.
* \param	dstIDom			Interval domain to which
*                                       intervals are to be appended.
* \param	grdMBuf			Integer modulus of grey gradient
*                                       buffer.
* \param	grdYBuf			Integer vertical grey gradient
*                                       buffer.
* \param	grdXBuf			buffer.
*               int *grdXBuf:           Integer horizontal grey gradient
*                                       buffer.
* \param	iPool			Interval pool.
* \param	dstBuf			Buffer for direction values.
* \param	dstLen			Buffer (given interval) length.
* \param	dstPos			Position of start of buffer wrt
*                                       the origin.
* \param	orgPos			The origin.
* \param	minGM			Minimum (modulus) gradient value
*                                       to be considered.
*/
static WlzErrorNum WlzNMSuppress2DBufI(WlzIntervalDomain *dstIDom,
				       int **grdMBuf, int *grdYBuf,
				       int *grdXBuf,
				       WlzDynItvPool *iPool,
				       WlzUByte *dstBuf, int dstLen,
				       WlzIVertex2 dstPos, WlzIVertex2 orgPos,
				       int minGM)
{
  int		cnt,
		itvLft,
		itvLen,
		maximal,
  		g0,
  		g1,
		gM,
		gY,
		gX,
		gMLft = 0,
		gMRgt = 0;
  WlzUByte	dCode;
  int		*grdMBufPrv,
  		*grdMBufCur,
		*grdMBufNxt;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzUByte dTable[8] = {3, 2, 0, 1, 4, 5, 7, 6};

  grdMBufPrv = *(grdMBuf + ((dstPos.vtY + 3 - 1) % 3)) + dstPos.vtX;
  grdMBufCur = *(grdMBuf + ((dstPos.vtY + 3 + 0) % 3)) + dstPos.vtX;
  grdMBufNxt = *(grdMBuf + ((dstPos.vtY + 3 + 1) % 3)) + dstPos.vtX;
  itvLen = 0;
  cnt = dstLen - 2;
  *(dstBuf + 0) = 0;
  *(dstBuf + dstLen - 1) = 0;
  itvLft = ++(dstPos.vtX);
  while((errNum == WLZ_ERR_NONE) && (cnt-- > 0))
  {
    ++grdXBuf;
    ++grdYBuf;
    *++dstBuf = 0;
    maximal = 0;
    gM = *(grdMBufCur + 1);
    if(gM && (gM > minGM))
    {
      gX = *grdXBuf;
      gY = *grdYBuf;
      /* Compute direction of maximum gradient. */
      dCode = dTable[((gY >= 0) << 2) |
      	             ((gX >= 0) << 1) |
	             ((gY * gY) >= (gX * gX))];
      /* Interpolate magnitudes. */
      switch(dCode)
      {
	case 0: /* \" */
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufNxt + 0);
	  gMLft = (((gM - g0) * gX) - ((g0 - g1) * gY)) / gM;
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufPrv + 2);
	  gMRgt = (((gM - g0) * gX) - ((g0 - g1) * gY)) / gM;
	  break;
	case 1: /* |\ */
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 0);
	  gMLft = (((g0 - g1) * gX) - ((gM - g0) * gY)) / gM;
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 2);
	  gMRgt = (((g0 - g1) * gX) - ((gM - g0) * gY)) / gM;
	  break;
	case 2: /* /| */
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 2);
	  gMLft = (((g1 - g0) * gX) - ((gM - g0) * gY)) / gM;
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 0);
	  gMRgt = (((g1 - g0) * gX) - ((gM - g0) * gY)) / gM;
	  break;
	case 3: /* "/ */
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufNxt + 2);
	  gMLft = (((g0 - gM) * gX) - ((g0 - g1) * gY)) / gM;
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufPrv + 0);
	  gMRgt = (((g0 - gM) * gX) - ((g0 - g1) * gY)) / gM;
	  break;
	case 4: /* _\ */
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufPrv + 2);
	  gMLft = (((g0 - gM) * gX) - ((g1 - g0) * gY)) / gM;
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufNxt + 0);
	  gMRgt = (((g0 - gM) * gX) - ((g1 - g0) * gY)) / gM;
	  break;
	case 5: /* \| */
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 2);
	  gMLft = (((g1 - g0) * gX) - ((g0 - gM) * gY)) / gM;
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 0);
	  gMRgt = (((g1 - g0) * gX) - ((g0 - gM) * gY)) / gM;
	  break;
	case 6: /* |/ */
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 0);
	  gMLft = (((g0 - g1) * gX) - ((g0 - gM) * gY)) / gM;
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 2);
	  gMRgt = (((g0 - g1) * gX) - ((g0 - gM) * gY)) / gM;
	  break;
	case 7: /* /_ */
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufPrv + 0);
	  gMLft = (((gM - g0) * gX) - ((g1 - g0) * gY)) / gM;
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufNxt + 2);
	  gMRgt = (((gM - g0) * gX) - ((g1 - g0) * gY)) / gM;
	  break;
      }
      /* Determine if maximal */
      maximal = (gMLft > 0) && (gMRgt > 0);
    }
    if(maximal)
    {
      if(itvLen)
      {
        ++itvLen;
      }
      else
      {
	itvLen = 1;
        itvLft = dstPos.vtX;
      }
      *dstBuf = (WlzUByte )(dCode | 0x80);
    }
    if((itvLen > 0) && ((cnt == 0) || (maximal == 0)))
    {
      errNum = WlzDynItvAdd(dstIDom, iPool, dstPos.vtY + orgPos.vtY,
			    itvLft + orgPos.vtX, itvLen);
      itvLen = 0;
    }
    ++grdMBufPrv;
    ++grdMBufCur;
    ++grdMBufNxt;
    ++(dstPos.vtX);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Performs non-maximal suppression on the given buffers.
* \param	dstIDom			Interval domain to which
*                                       intervals are to be appended.
* \param	grdMBuf			Double modulus of grey gradient
*                                       buffer.
* \param	grdYBuf			Double vertical grey gradient
*                                       buffer.
* \param	grdXBuf			Double horizontal grey gradient
*                                       buffer.
* \param	iPool			Interval pool.
* \param	dstBuf			Buffer for direction values.
* \param	dstLen			Buffer (given interval) length.
* \param	dstPos			Position of start of buffer wrt
*                                       the origin.
* \param	orgPos			The origin.
* \param	minGM			Minimum (modulus) gradient value
*                                       to be considered.
*/
static WlzErrorNum WlzNMSuppress2DBufD(WlzIntervalDomain *dstIDom,
				       double **grdMBuf, double *grdYBuf,
				       double *grdXBuf,
				       WlzDynItvPool *iPool,
				       WlzUByte *dstBuf, int dstLen,
				       WlzIVertex2 dstPos, WlzIVertex2 orgPos,
				       double minGM)
{
  int		cnt,
		itvLft,
		itvLen,
  		maximal;
  WlzUByte	dCode;
  double	g0,
  		g1,
		gM,
		gY,
		gX,
		gMLft = 0.0,
		gMRgt = 0.0;
  double	*grdMBufPrv,
  		*grdMBufCur,
		*grdMBufNxt;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const WlzUByte dTable[8] = {3, 2, 0, 1, 4, 5, 7, 6};

  grdMBufPrv = *(grdMBuf + ((dstPos.vtY + 3 - 1) % 3)) + dstPos.vtX;
  grdMBufCur = *(grdMBuf + ((dstPos.vtY + 3 + 0) % 3)) + dstPos.vtX;
  grdMBufNxt = *(grdMBuf + ((dstPos.vtY + 3 + 1) % 3)) + dstPos.vtX;
  itvLen = 0;
  cnt = dstLen - 2;
  *(dstBuf + 0) = 0;
  *(dstBuf + dstLen - 1) = 0;
  itvLft = ++(dstPos.vtX);
  while((errNum == WLZ_ERR_NONE) && (cnt-- > 0))
  {
    ++grdXBuf;
    ++grdYBuf;
    maximal = 0;
    *++dstBuf = 0;
    gM = *(grdMBufCur + 1);
    if(((gM * gM) > DBL_EPSILON) && (gM > minGM))
    {
      gX = *grdXBuf;
      gY = *grdYBuf;
      /* Compute direction of maximum gradient. */
      dCode = dTable[((gY >= 0.0) << 2) |
      	             ((gX >= 0.0) << 1) |
	             ((gY * gY) >= (gX * gX))];
      /* Interpolate magnitudes. */
      switch(dCode)
      {
	case 0: /* \" */
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufNxt + 0);
	  gMLft = (((gM - g0) * gX) - ((g0 - g1) * gY)) / gM;
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufPrv + 2);
	  gMRgt = (((gM - g0) * gX) - ((g0 - g1) * gY)) / gM;
	  break;
	case 1: /* |\ */
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 0);
	  gMLft = (((g0 - g1) * gX) - ((gM - g0) * gY)) / gM;
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 2);
	  gMRgt = (((g0 - g1) * gX) - ((gM - g0) * gY)) / gM;
	  break;
	case 2: /* /| */
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 2);
	  gMLft = (((g1 - g0) * gX) - ((gM - g0) * gY)) / gM;
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 0);
	  gMRgt = (((g1 - g0) * gX) - ((gM - g0) * gY)) / gM;
	  break;
	case 3: /* "/ */
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufNxt + 2);
	  gMLft = (((g0 - gM) * gX) - ((g0 - g1) * gY)) / gM;
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufPrv + 0);
	  gMRgt = (((g0 - gM) * gX) - ((g0 - g1) * gY)) / gM;
	  break;
	case 4: /* _\ */
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufPrv + 2);
	  gMLft = (((g0 - gM) * gX) - ((g1 - g0) * gY)) / gM;
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufNxt + 0);
	  gMRgt = (((g0 - gM) * gX) - ((g1 - g0) * gY)) / gM;
	  break;
	case 5: /* \| */
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 2);
	  gMLft = (((g1 - g0) * gX) - ((g0 - gM) * gY)) / gM;
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 0);
	  gMRgt = (((g1 - g0) * gX) - ((g0 - gM) * gY)) / gM;
	  break;
	case 6: /* |/ */
	  g0 = *(grdMBufPrv + 1);
	  g1 = *(grdMBufPrv + 0);
	  gMLft = (((g0 - g1) * gX) - ((g0 - gM) * gY)) / gM;
	  g0 = *(grdMBufNxt + 1);
	  g1 = *(grdMBufNxt + 2);
	  gMRgt = (((g0 - g1) * gX) - ((g0 - gM) * gY)) / gM;
	  break;
	case 7: /* /_ */
	  g0 = *(grdMBufCur + 0);
	  g1 = *(grdMBufPrv + 0);
	  gMLft = (((gM - g0) * gX) - ((g1 - g0) * gY)) / gM;
	  g0 = *(grdMBufCur + 2);
	  g1 = *(grdMBufNxt + 2);
	  gMRgt = (((gM - g0) * gX) - ((g1 - g0) * gY)) / gM;
	  break;
      }
      /* Determine if maximal */
      maximal = (gMLft > DBL_EPSILON) && (gMRgt > DBL_EPSILON);
    }
    if(maximal)
    {
      *dstBuf = (WlzUByte )(dCode | 0x80);
    }
    if((itvLen > 0) && ((cnt == 0) || (maximal == 0)))
    {
      errNum = WlzDynItvAdd(dstIDom, iPool, dstPos.vtY + orgPos.vtY,
			    itvLft + orgPos.vtX, itvLen);
      itvLen = 0;
    }
    ++grdMBufPrv;
    ++grdMBufCur;
    ++grdMBufNxt;
  }
  return(errNum);
}

/*!
* \return	New Woolz domain object with maximal domain and grey
*		values which encode the gradient's direction or NULL
*		on error.
* \ingroup	WlzFeatures
* \brief	Computes the maximal domain and gradient direction of
*               given Woolz 2D domain object.
* \note         All the objects domains are known to be the same.
* \param	grdM			Gradient magnitude.
* \param	grdY			Gradient (partial derivative)
*                                       through lines.
* \param	grdX			Gradient (partial derivative)
*                                       through columns.
* \param	minThrV			Minimum gradient value to
*                                       consider.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzNMSuppress2D(WlzObject *grdM,
				  WlzObject *grdY, WlzObject *grdX,
				  WlzPixelV minThrV, WlzErrorNum *dstErr)
{
  int		idN,
		inLen,
		outLen,
  		inLnIdx = 0;
  WlzGreyType	gType,
  		bufType;
  WlzIVertex2	bufSz,
  		inPos,
		outPos,
		orgPos;
  WlzValues	tmpVal;
  WlzDomain	dstDom,
  		grdDom;
  WlzIntervalWSpace tmpIWSp = {0},
  		grdMIWSp = {0},
  		grdYIWSp = {0},
		grdXIWSp = {0};
  WlzGreyWSpace tmpGWSp,
  		grdMGWSp,
  		grdYGWSp,
		grdXGWSp;
  WlzPixelV	zeroV;
  WlzGreyP	grdMBufGP,
  		grdYBufGP,
		grdXBufGP;
  WlzDynItvPool pool;
  WlzObject	*dstObj = NULL,
  		*tmpObj = NULL;
  void 		*grdYBuf = NULL,
		*grdXBuf = NULL;
  void		**grdMBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tmpVal.core = NULL;
  pool.itvBlock = NULL;
  dstDom.core = NULL;
  if((grdM->type != WLZ_2D_DOMAINOBJ) ||
     (grdY->type != WLZ_2D_DOMAINOBJ) ||
     (grdX->type != WLZ_2D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((grdM->domain.core == NULL) ||
          (grdY->domain.core == NULL) ||
          (grdX->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((grdM->values.core == NULL) ||
  	  (grdY->values.core == NULL) ||
	  (grdX->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    /* Find required buffer type (WLZ_GREY_DOUBLE or WLZ_GREY_INT). */
    bufType = WLZ_GREY_INT;
    gType = WlzGreyTableTypeToGreyType(grdM->values.core->type, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if((gType == WLZ_GREY_FLOAT) || (gType == WLZ_GREY_DOUBLE))
      {
	bufType = WLZ_GREY_DOUBLE;
      }
      else
      {
	gType = WlzGreyTableTypeToGreyType(grdY->values.core->type, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  if((gType == WLZ_GREY_FLOAT) || (gType == WLZ_GREY_DOUBLE))
	  {
	    bufType = WLZ_GREY_DOUBLE;
	  }
	  else
	  {
	    gType = WlzGreyTableTypeToGreyType(grdX->values.core->type,
	    				       &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if((gType == WLZ_GREY_FLOAT) || (gType == WLZ_GREY_DOUBLE))
	      {
		bufType = WLZ_GREY_DOUBLE;
	      }
	    }
	  }
	}
      }
    }
  }
  /* Convert minimum gradient threshold value. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(bufType == WLZ_GREY_INT)
    {
      errNum = WlzValueConvertPixel(&minThrV, minThrV, WLZ_GREY_INT);
    }
    else /* bufType == WLZ_GREY_DOUBLE */
    {
      errNum = WlzValueConvertPixel(&minThrV, minThrV, WLZ_GREY_DOUBLE);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    grdDom = grdM->domain;
    /* Make destination object with WLZ_GREY_UBYTE greys. */
    zeroV.type = WLZ_GREY_UBYTE;
    zeroV.v.inv = 0;
    tmpVal.v = WlzNewValueTb(grdM, WlzGreyTableType(WLZ_GREY_TAB_RAGR,
						    WLZ_GREY_UBYTE, NULL),
			     zeroV, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      /* Use the input domain while calculating the new maximal domain. */
      tmpObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, grdM->domain, tmpVal,
			   NULL, NULL, &errNum);
    }
  }
  /* Initialize the memory pool with some size of block. Any +ve number
   * greater than the maximum number of intervals in any destination line
   * would work but the fewer allocations then the more efficient the code,
   * hence this attempt to guess the required number of intervals in the
   * destination domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    pool.itvsInBlock = (((grdDom.i->lastkl - grdDom.i->kol1 + 1) *
			(grdDom.i->lastln - grdDom.i->line1 + 1)) / 64) +
		       grdDom.i->lastkl - grdDom.i->kol1 + 1024;
  }
  /* Make gradient buffers. */
  if(errNum == WLZ_ERR_NONE)
  {
    bufSz.vtY = 3;
    bufSz.vtX = grdDom.i->lastkl - grdDom.i->kol1 + 1;
    if(bufType == WLZ_GREY_INT)
    {
      if((AlcInt2Malloc((int ***)&grdMBuf,
			bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
	 ((grdYBuf = AlcMalloc(sizeof(int) * bufSz.vtX)) == NULL) ||
	 ((grdXBuf = AlcMalloc(sizeof(int) * bufSz.vtX)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	grdYBufGP.inp = (int *)grdYBuf;
	grdXBufGP.inp = (int *)grdXBuf;
      }
    }
    else /* bufType == WLZ_GREY_DOUBLE */
    {
      if((AlcDouble2Malloc((double ***)&grdMBuf,
			   bufSz.vtY, bufSz.vtX) != ALC_ER_NONE) ||
	 ((grdYBuf = AlcMalloc(sizeof(double) * bufSz.vtX)) == NULL) ||
	 ((grdXBuf = AlcMalloc(sizeof(double) * bufSz.vtX)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	grdYBufGP.dbp = (double *)grdYBuf;
	grdXBufGP.dbp = (double *)grdXBuf;
      }
    }
  }
  /* Make destination interval domain with interval lines but not intervals. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
    				     grdDom.i->line1, grdDom.i->lastln,
				     grdDom.i->kol1, grdDom.i->lastkl,
				     &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Scan down through the gradient objects. */
    if(((errNum = WlzInitGreyScan(tmpObj,
    				  &tmpIWSp, &tmpGWSp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(grdM,
    				  &grdMIWSp, &grdMGWSp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(grdY,
    				  &grdYIWSp, &grdYGWSp)) == WLZ_ERR_NONE) &&
       ((errNum = WlzInitGreyScan(grdX,
    				  &grdXIWSp, &grdXGWSp)) == WLZ_ERR_NONE))
    {
      orgPos.vtX = grdDom.i->kol1;
      orgPos.vtY = grdDom.i->line1;
      while((errNum == WLZ_ERR_NONE) &&
            ((errNum = WlzNextGreyInterval(&grdMIWSp)) == WLZ_ERR_NONE))
      {
        inLen = grdMIWSp.rgtpos - grdMIWSp.lftpos + 1;
	inPos.vtX = grdMIWSp.lftpos - orgPos.vtX;
	/* Process any lines between this and the last by clearing the
	 * gradient magnitude buffer . */
	if(grdMIWSp.nwlpos > 0)
	{
	  idN = (grdMIWSp.nwlpos >= 3)? 3: grdMIWSp.nwlpos;
	  while(--idN >= 0)
	  {
	    inPos.vtY = grdMIWSp.linpos - orgPos.vtY - idN;
	    inLnIdx = (3 + inPos.vtY) % 3;
	    if(bufType == WLZ_GREY_INT)
	    {
	      WlzValueSetInt(*((int **)grdMBuf + inLnIdx), 0,
	      		     bufSz.vtX);
	    }
	    else /* bufType == WLZ_GREY_DOUBLE */
	    {
	      WlzValueSetDouble(*((double **)grdMBuf + inLnIdx), 0,
	      			bufSz.vtX);
	    }
	  }
	}
	/* Copy intervals to values buffers. */
	if(bufType == WLZ_GREY_INT)
	{
	  grdMBufGP.inp = *((int **)grdMBuf + inLnIdx);
	}
	else /* bufType == WLZ_GREY_DOUBLE */
	{
	  grdMBufGP.dbp = *((double **)grdMBuf + inLnIdx);
	}
	WlzValueCopyGreyToGrey(grdMBufGP, inPos.vtX, bufType,
			       grdMGWSp.u_grintptr, 0, grdMGWSp.pixeltype,
			       inLen);
	if(grdMIWSp.intrmn == 0)
	{
	  while((errNum == WLZ_ERR_NONE) &&
	        (tmpIWSp.linpos < grdMIWSp.linpos))
	  {
	    outPos.vtY = tmpIWSp.linpos - orgPos.vtY;
	    if(outPos.vtY >= 0)
	    {
	      outLen = tmpIWSp.rgtpos - tmpIWSp.lftpos + 1;
	      outPos.vtX = tmpIWSp.lftpos - orgPos.vtX;
	      WlzValueCopyGreyToGrey(grdYBufGP, 0, bufType,
				     grdYGWSp.u_grintptr, 0,
				     grdYGWSp.pixeltype, outLen);
	      WlzValueCopyGreyToGrey(grdXBufGP, 0, bufType,
				     grdXGWSp.u_grintptr, 0,
				     grdXGWSp.pixeltype, outLen);
	      if(bufType == WLZ_GREY_INT)
	      {
		errNum = WlzNMSuppress2DBufI(dstDom.i,
					     (int **)grdMBuf,
					     (int *)grdYBuf,
					     (int *)grdXBuf,
					     &pool,
					     tmpGWSp.u_grintptr.ubp,
					     outLen, outPos, orgPos,
					     minThrV.v.inv);
	      }
	      else /* bufType == WLZ_GREY_DOUBLE */
	      {
		errNum = WlzNMSuppress2DBufD(dstDom.i,
				    	     (double **)grdMBuf,
				    	     (double *)grdYBuf,
					     (double *)grdXBuf,
					     &pool,
					     tmpGWSp.u_grintptr.ubp,
					     outLen, outPos, orgPos,
					     minThrV.v.dbv);
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzNextGreyInterval(&tmpIWSp);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzNextGreyInterval(&grdYIWSp);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzNextGreyInterval(&grdXIWSp);
	    }
	  }
        }
      }
      if(errNum == WLZ_ERR_EOO)
      {
        errNum = WLZ_ERR_NONE;
      }
    }
    if(tmpIWSp.gryptr == &tmpGWSp)
    {
      (void )WlzEndGreyScan(&tmpIWSp, &tmpGWSp);
    }
    if(grdMIWSp.gryptr == &grdMGWSp)
    {
      (void )WlzEndGreyScan(&grdMIWSp, &grdMGWSp);
    }
    if(grdYIWSp.gryptr == &grdYGWSp)
    {
      (void )WlzEndGreyScan(&grdYIWSp, &grdYGWSp);
    }
    if(grdXIWSp.gryptr == &grdXGWSp)
    {
      (void )WlzEndGreyScan(&grdXIWSp, &grdXGWSp);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((errNum = WlzStandardIntervalDomain(dstDom.i)) == WLZ_ERR_NONE)
    {
      dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, tmpVal,
			   NULL, NULL, &errNum);
    }
  }
  if(tmpObj)
  {
    WlzFreeObj(tmpObj);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(tmpObj == NULL)
    {
      if(dstDom.core)
      {
        (void )WlzFreeDomain(dstDom);
      }
      if(tmpVal.core)
      {
        (void )WlzFreeValues(tmpVal);
      }
    }
  }
  if(grdMBuf)
  {
    Alc2Free(grdMBuf);
  }
  if(grdYBuf)
  {
    AlcFree(grdYBuf);
  }
  if(grdXBuf)
  {
    AlcFree(grdXBuf);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New Woolz domain object with maximal domain and grey
*		values which encode the gradient's direction or NULL
*		on error.
* \ingroup	WlzFeatures
* \brief	Computes the maximal domain and gradient direction of
*               given Woolz 3D domain object.
* \note		All the objects domains are known to be the same.
* \param	grdM			Gradient magnitude.
* \param	grdZ			Gradient (partial derivative)
*                                       through planes.
* \param	grdY			Gradient (partial derivative)
*                                       through lines.
* \param	grdX			Gradient (partial derivative)
*                                       through columns.
* \param	minThrV			Minimum gradient value to
*                                       consider.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
static WlzObject *WlzNMSuppress3D(WlzObject *grdM, WlzObject *grdZ,
				  WlzObject *grdY, WlzObject *grdX,
				  WlzPixelV minThrV, WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return
* \brief	Computes the maximal domain and gradient direction of
*               given Woolz domain object.
*
*		Currently only implemented for 2D domain objects.
*               The domain is the maximaly suppressed domain and the values
*               are the encoded gradient direction. The direction is encoding 
*               is from the +ve x-axis counter clockwise in eight steps
*               with a mask of 0x80, ie directions values are in the
*               range 128 -> 128 + 7.
*		\verbatim
                             ^ Y axis (downwards when displayed)
                             |
                  +----------+---------+
                  | \        |        /|
                  |  \128 + 2|128 + 1/ |
                  |   \      |      /  |
                  |    \     |     /   |
                  |     \    |    /    |
                  |      \   |   /     |
                  |       \  |  /      |
                  |128 + 3 \ | /128 + 0|
                  |         \|/        |
                  +----------O---------+--> X axis
                  |         /|\        |
                  |128 + 4 / | \128 + 7|
                  |       /  |  \      |
                  |      /   |   \     |
                  |     /    |    \    |
                  |    /     |     \   |
                  |   /      |      \  |
                  |  /128 + 5|128 + 6\ |
                  | /        |        \|
                  +----------+---------+
		\endverbatim

* \param	grdM			Gradient magnitude.
* \param	grdZ			Gradient (partial derivative)
*                                       through planes.
* \param	grdY			Gradient (partial derivative)
*                                       through lines.
* \param	grdX			Gradient (partial derivative)
*                                       through columns.
* \param	minThrV			Minimum gradient value to
*                                       consider.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzObject	*WlzNMSuppress(WlzObject *grdM, WlzObject *grdZ,
			       WlzObject *grdY, WlzObject *grdX,
			       WlzPixelV minThrV, WlzErrorNum *dstErr)
{
  int		idN;
  WlzObject	*istObj = NULL,
  		*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*tObjs[4];

  if(grdM == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(grdM->type)
    {
      case WLZ_2D_DOMAINOBJ:
	tObjs[0] = grdM;
	tObjs[1] = grdY;
	tObjs[2] = grdX;
	istObj = WlzIntersectN(3, tObjs, 0, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(istObj->type)
	  {
            case WLZ_2D_DOMAINOBJ:
	      tObjs[1] = tObjs[2] = NULL;
	      tObjs[0] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	      			     istObj->domain, grdM->values,
	      		     	     NULL, NULL, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		tObjs[1] = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
				       istObj->domain, grdY->values,
				       NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		tObjs[2] = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
				       istObj->domain, grdX->values,
				       NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		dstObj = WlzNMSuppress2D(tObjs[0], tObjs[1], tObjs[2],
					 minThrV, &errNum);
	      }
	      for(idN = 0; idN < 3; ++idN)
	      {
		if(tObjs[0])
		{
		  WlzFreeObj(tObjs[0]);
		}
	      }
	      break;
	    case WLZ_EMPTY_OBJ:
	      dstObj = WlzMakeEmpty(&errNum);
	      break;
	    default:
	      errNum = WLZ_ERR_OBJECT_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	tObjs[0] = grdM;
	tObjs[1] = grdZ;
	tObjs[2] = grdY;
	tObjs[3] = grdX;
	istObj = WlzIntersectN(4, tObjs, 0, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  switch(istObj->type)
	  {
            case WLZ_2D_DOMAINOBJ:
	      tObjs[1] = tObjs[2] = tObjs[3] = NULL;
	      tObjs[0] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	      			     istObj->domain, grdM->values,
	      		     	     NULL, NULL, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		tObjs[1] = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
				       istObj->domain, grdZ->values,
				       NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		tObjs[2] = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
				       istObj->domain, grdY->values,
				       NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		tObjs[3] = WlzMakeMain(WLZ_2D_DOMAINOBJ, 
				       istObj->domain, grdX->values,
				       NULL, NULL, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		dstObj = WlzNMSuppress3D(tObjs[0], tObjs[1], tObjs[2], tObjs[3],
					 minThrV, &errNum);
	      }
	      for(idN = 0; idN < 4; ++idN)
	      {
		if(tObjs[0])
		{
		  WlzFreeObj(tObjs[0]);
		}
	      }
	      break;
	    case WLZ_EMPTY_OBJ:
	      dstObj = WlzMakeEmpty(&errNum);
	      break;
	    default:
	      errNum = WLZ_ERR_OBJECT_TYPE;
	      break;
	  }
	}
        break;
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
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

/* #define TEST_WLZNMSUPRESS */
#ifdef TEST_WLZNMSUPRESS
static WlzObject *WlzCanny2(WlzObject *srcObj, double alpha,
			   WlzPixelV lMinGrdV, WlzPixelV hMinGrdV,
			   WlzErrorNum *dstErr)
{
  WlzObject     *gMObj = NULL,
		*gXObj = NULL,
		*gYObj = NULL,
		*dstObj = NULL;
  WlzRsvFilter  *ftr = NULL;
  WlzPixelV     zeroBV;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  zeroBV.type = WLZ_GREY_INT;
  zeroBV.v.inv = 0;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_2D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((ftr = WlzRsvFilterMakeFilter(WLZ_RSVFILTER_NAME_DERICHE_1,
					alpha, &errNum)) != NULL)
  {
    ftr->c *= 4;
    gMObj = WlzAssignObject(
	    WlzGreyGradient(NULL, &gYObj, &gXObj, srcObj, ftr,
			    &errNum), NULL);
    WlzRsvFilterFreeFilter(ftr);
    if(gXObj && gYObj && gMObj && (errNum == WLZ_ERR_NONE))
    {
      dstObj = WlzAssignObject(
	       WlzNMSuppress(gMObj, NULL, gYObj, gXObj, lMinGrdV,
	       		     &errNum), NULL);
    }
  }
  if(gMObj)
  {
    WlzFreeObj(gMObj);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstErr)
    {
      *dstErr = errNum;
    }
    dstObj = NULL;
  }
  return(dstObj);
}

int             main(int argc, char *argv[])
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzPixelV	lMinGrdV,
  		hMinGrdV;
  double        alpha;
  WlzObject     *inObj = NULL,
		*outObj = NULL;

  alpha = 1.0;
  lMinGrdV.type = WLZ_GREY_INT;
  lMinGrdV.v.inv = 2;
  hMinGrdV.type = WLZ_GREY_INT;
  lMinGrdV.v.inv = 8;
  if(((inObj = WlzAssignObject(WlzReadObj(stdin, &errNum), NULL)) == NULL) ||
     (errNum != WLZ_ERR_NONE))
  {
    (void )fprintf(stderr,
		   "%s: failed to read object from stdin\n",
		   *argv);
  }
  else if(((outObj = WlzCanny2(inObj, alpha, lMinGrdV, hMinGrdV,
  			       &errNum)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
  {
    (void )fprintf(stderr,
		   "%s: failed to Canny filter object\n",
		   *argv);
  }
  else if(WlzWriteObj(stdout, outObj) != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr,
		   "%s: failed to write object\n",
		   *argv);
  }
  return(errNum);
}
#endif /* TEST_WLZNMSUPRESS */
