#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzHyThreshold_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzHyThreshold.c
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
* \brief	A hysteresis threshold filter.
* \ingroup	WlzThreshold
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/*!
* \return	Thresholded object or NULL on error.
* \ingroup      WlzThreshold
* \brief	Hysteresis thresholds the given Woolz object.
*               Values are in the domain of the hysteresis threshold'd
*               object if they are above/below the primary threshold
*               or above/below the secondary threshold and connected
*               to values above/below the primary threshold.
* \param	srcObj			Object to be thresholded.
* \param	pThrV			Primary hysteresis threshold
*                                       value.
* \param	sThrV			Threshold for above or below
*                                       values.
* \param	hilo			Threshold for above or below
*                                       values.
* \param	con			Connectivity to examine for
*                                       hysteresis.
* \param	dstErr			Destination error pointer, may
*                                       be null.
*/
WlzObject	*WlzHyThreshold(WlzObject *srcObj,
			       WlzPixelV pThrV, WlzPixelV sThrV,
			       WlzThresholdType hilo, WlzConnectType con,
			       WlzErrorNum *dstErr)
{
  int		simpleThr = 0;
  WlzPixelV	tmpV;
  WlzObject	*dstObj = NULL,
  		*pThrObj = NULL,
  		*sThrObj = NULL,
		*dThrObj = NULL,
		*iThrObj = NULL,
		*uThrObj;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcObj == NULL)
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
    if(con == WLZ_0_CONNECTED)
    {
      simpleThr = 1;
    }
    else
    {
      if((errNum = WlzValueConvertPixel(&tmpV, sThrV,
      					pThrV.type)) == WLZ_ERR_NONE)
      {
        switch(pThrV.type)
	{
	  case WLZ_GREY_INT:
	    if(tmpV.v.inv == pThrV.v.inv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    if(tmpV.v.shv == pThrV.v.shv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    if(tmpV.v.ubv == pThrV.v.ubv)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    if(fabs(tmpV.v.flv - pThrV.v.flv) <= FLT_EPSILON)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    if(fabs(tmpV.v.dbv - pThrV.v.dbv) <= DBL_EPSILON)
	    {
	      simpleThr = 1;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    if( tmpV.v.rgbv == pThrV.v.rgbv )
	    {
	      simpleThr = 1;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(simpleThr)
    {
      dstObj = WlzThreshold(srcObj, pThrV, hilo, &errNum);
    }
    else
    {
      pThrObj = WlzThreshold(srcObj, pThrV, hilo, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	if(pThrObj->type == WLZ_EMPTY_OBJ)
	{
	  dstObj = pThrObj;
	  pThrObj = NULL;
	}
	else
	{
	  sThrObj = WlzThreshold(srcObj, sThrV, hilo, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(sThrObj->type == WLZ_EMPTY_OBJ)
	    {
	      dstObj = pThrObj;
	      pThrObj = NULL;
	    }
	    else
	    {
	      dThrObj = WlzDilation(pThrObj, con, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		iThrObj = WlzIntersect2(dThrObj, sThrObj, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		uThrObj = WlzUnion2(pThrObj, iThrObj, &errNum);
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		dstObj = WlzMakeMain(uThrObj->type,
				     uThrObj->domain, srcObj->values,
				     srcObj->plist, srcObj, &errNum);
	      }
	      if(dThrObj)
	      {
		WlzFreeObj(dThrObj);
	      }
	      if(iThrObj)
	      {
		WlzFreeObj(iThrObj);
	      }
	      if(uThrObj)
	      {
		WlzFreeObj(uThrObj);
	      }
	    }
	    if(sThrObj)
	    {
	      WlzFreeObj(sThrObj);
	    }
	    if(pThrObj)
	    {
	      WlzFreeObj(pThrObj);
	    }
	  }
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}
