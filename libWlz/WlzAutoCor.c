#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzAutoCor_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzAutoCor.c
* \author       Bill Hill
* \date         November 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Functions to compute the autocorrelation of an object.
* \ingroup	WlzRegistration
* \todo         -
* \bug          None known.
*/

#include <float.h>
#include <limits.h>
#include <Wlz.h>

static WlzObject 		*WlzAutoCor2D(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);
static void			WlzAutoCorRearrange2D(
				  double **aAr,
				  WlzIVertex2 aSz,
				  double **wAr,
				  WlzIVertex2 wSz);

/*!
* \return	Autocorrelated object or NULL on error.
* \ingroup	WlzRegistration
* \brief	Computes the autocorrelation of the given object.
*		The autocorrelation object will have double values,
*		an origin of (0,0) and a column and line sizes which
*		are the integer power of two greater than or equal to
*		those of the given object.
*		See AlgAutoCorrelate2D() for the organisation of the
*		autocorrelation data.
* \param	gObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzAutoCor(WlzObject *gObj, WlzErrorNum *dstErr)
{
  WlzObject	*aObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	aObj = WlzAutoCor2D(gObj, &errNum);
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
  return(aObj);
}

/*!
* \return	Autocorrelated object or NULL on error.
* \ingroup	WlzRegistration
* \brief	Computes the autocorrelation of the given 2D object, see
*		WlzAutoCor().
* \param	gObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzAutoCor2D(WlzObject *gObj, WlzErrorNum *dstErr)
{
  WlzIVertex2	aSz,
		wSz,
		aOrg,
  		wOrg;
  WlzIBox2	box;
  double	**wAr = NULL,
  		**aAr = NULL;
  WlzObject     *aObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    box = WlzBoundingBox2I(gObj, &errNum);
    aSz.vtX = box.xMax - box.xMin + 1;
    aSz.vtY = box.yMax - box.yMin + 1;
    /* Make sure aSz is even in x and y. */
    if((aSz.vtX & 1) != 0)
    {
      aSz.vtX += 1;
    }
    if((aSz.vtY & 1) != 0)
    {
      aSz.vtY += 1;
    }
    wOrg.vtX = box.xMin - (aSz.vtX / 2);
    wOrg.vtY = box.yMin - (aSz.vtY / 2);
    wSz.vtX = aSz.vtX * 2;
    wSz.vtY = aSz.vtY * 2;
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(wSz.vtX), wSz.vtX);
    (void )AlgBitNextPowerOfTwo((unsigned int *)&(wSz.vtY), wSz.vtY);
    errNum = WlzToArray2D((void ***)&wAr, gObj, wSz, wOrg, 0, WLZ_GREY_DOUBLE);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )AlgAutoCorrelate2D(wAr, wSz.vtX, wSz.vtY);
    if(AlcDouble2Malloc(&aAr, aSz.vtY, aSz.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzAutoCorRearrange2D(aAr, aSz, wAr, wSz);
    aOrg.vtX = -(aSz.vtX / 2);
    aOrg.vtY = -(aSz.vtY / 2);
    aObj = WlzFromArray2D((void **)aAr, aSz, aOrg,
                          WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
			  0, 0, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(aObj != NULL)
    {
      (void )WlzFreeObj(aObj);
    }
  }
  (void )Alc2Free((void **)wAr);
  (void )Alc2Free((void **)aAr);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(aObj);
}

static void	WlzAutoCorRearrange2D(double **aAr, WlzIVertex2 aSz,
				      double **wAr, WlzIVertex2 wSz)
{
  int		idX,
  		idY;
  double	*botLftA,
		*botRgtA,
        	*topLftA,
		*topRgtA,
  		*botLftW,
		*botRgtW,
  		*topLftW,
		*topRgtW;
  WlzIVertex2	hlfA;

  hlfA.vtX = aSz.vtX / 2;
  hlfA.vtY = aSz.vtY / 2;
  for(idY = 0; idY < hlfA.vtY; ++idY)
  {
    botLftA = *(aAr + idY);
    botRgtA = botLftA + hlfA.vtX;
    topLftA = *(aAr + hlfA.vtY + idY);
    topRgtA = topLftA + hlfA.vtX;

    botLftW = *(wAr + idY);
    botRgtW = botLftW + wSz.vtX - hlfA.vtX;
    topLftW = *(wAr + wSz.vtY - hlfA.vtY + idY);
    topRgtW = topLftW + wSz.vtX - hlfA.vtX;

    for(idX = 0; idX < hlfA.vtX; ++idX)
    {
      *botLftA++ = *topRgtW++;
      *botRgtA++ = *topLftW++;
      *topLftA++ = *botRgtW++;
      *topRgtA++ = *botLftW++;
    }
  }
}
