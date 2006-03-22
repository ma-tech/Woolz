#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzPolarSample_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzPolarSample.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	A rectangular to polar image resampling filter.
* \ingroup	WlzValueFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Minimum value.
* \ingroup      WlzValueFilters
* \brief	Given four values, finds and returns the minimum.
* \param	val1			First value.
* \param	val2			Second value.
* \param	val3			Third value.
* \param	val4			Forth value.
*/
static int	WlzPolarSampleMinOf4I(int val1, int val2, int val3, int val4)
{
  if(val1 > val2)
  {
    val1 = val2;
  }
  if(val3 > val4)
  {
    val3 = val4;
  }
  if(val1 > val3)
  {
    val1 = val3;
  }
  return(val1);
}

/*!
* \return	Maximum value.
* \ingroup      WlzValueFilters
* \brief	Given four values, finds and returns the maximum.
* \param	val1			First value.
* \param	val2			Second value.
* \param	val3			Third value.
* \param	val4			Forth value.
*/
static int	WlzPolarSampleMaxOf4I(int val1, int val2, int val3, int val4)
{
  if(val1 < val2)
  {
    val1 = val2;
  }
  if(val3 < val4)
  {
    val3 = val4;
  }
  if(val1 < val3)
  {
    val1 = val3;
  }
  return(val1);
}

/*!
* \return	New polar sampled woolz object.
* \ingroup      WlzValueFilters
* \brief	Polar resamples the given woolz object using the given
*               angle and radial distance increments. All angles are in
*               radians.
* Notes:        The linkcount of the returned object is zero.
* \param	srcObj			Given Woolz object.
* \param	org			Origin in given object about
*                                       which to polar resample.
* \param	angleInc		Angle increment (radians).
* \param	distInc			Radial distance increment.
* \param	nLines			Number of lines (angular
*                                       samples).
* \param	outFlag			If non zero then the outer
*                                       bounding circle of the given
*                                       object is used, else the inner
*                                       bound is used for resampling.
* \param	wlzErr			Destination error pointer, may
*                                       be NULL.
*/
WlzObject	*WlzPolarSample(WlzObject *srcObj, WlzIVertex2 org,
			        double angleInc, double distInc,
			        int nLines, int outFlag,
				WlzErrorNum *wlzErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		tI0,
 		 countX,
		nCols;
  double	angle,
		dist,
		cosAngle,
		sinAngle;
  WlzGreyType	vType;
  WlzPixelV	bkgPix;
  WlzGreyP	dstValP;
  WlzObject	*dstObj = NULL;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzIBox2	box;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzPolarSample FE 0x%lx {%d %d} %g %g %d %d 0x%lx\n",
	   (unsigned long )srcObj, org.vtX, org.vtY, angleInc, distInc,
	   nLines, outFlag, (unsigned long )wlzErr));
  dstValP.ubp = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcObj->values.core == NULL)
	{
	  errNum = WLZ_ERR_VALUES_NULL;
	}
	else if((srcObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
	        (srcObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(distInc <= 0)
	{
	    errNum = WLZ_ERR_PARAM_DATA;
	}
	else
	{
	  box.xMin = org.vtX - srcObj->domain.i->kol1;
	  box.yMin = org.vtY - srcObj->domain.i->line1;
	  box.yMax = srcObj->domain.i->lastln - org.vtY;
	  box.xMax = srcObj->domain.i->lastkl - org.vtX;
	  if(outFlag)				    /* Find longest diagonal */
	  {
	    box.xMax *= box.xMax;
	    box.xMin *= box.xMin;
	    box.yMax *= box.yMax;
	    box.yMin *= box.yMin;
	    tI0 = WlzPolarSampleMaxOf4I(box.xMax + box.yMax,
					box.xMin + box.yMax,
					box.xMin + box.yMin,
					box.xMax  + box.yMin);
	    nCols = WLZ_NINT(sqrt((double )tI0) / distInc) + 1;
	  }
	  else					       /* Find shortest line */
	  {
	    tI0 = WlzPolarSampleMinOf4I(box.xMax, box.xMin,
					box.yMax, box.yMin);
	    nCols = WLZ_NINT(tI0 / distInc) + 1;
	  }
	  if((nCols <= 0) || (nLines <= 0))
	  {
	    errNum = WLZ_ERR_PARAM_DATA;
	  }
	  WLZ_DBG(WLZ_DBG_LVL_2,
		  ("WlzPolarSample 01 %d %d %d\n",
		   nCols, nLines, (int )errNum));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  vType = WlzGreyTableTypeToGreyType(srcObj->values.core->type,
	  				     &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    bkgPix = WlzGetBackground(srcObj, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tI0 = nCols * nLines;
	    switch(vType)
	    {
	      case WLZ_GREY_INT:
		if((dstValP.inp = (int *)AlcMalloc((unsigned long )tI0 *
						   sizeof(int))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_GREY_SHORT:
		if((dstValP.shp = (short *)AlcMalloc((unsigned long )tI0 *
						     sizeof(short))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_GREY_UBYTE:
		if((dstValP.ubp = (UBYTE *)AlcMalloc((unsigned long )tI0 *
						     sizeof(UBYTE))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_GREY_FLOAT:
		if((dstValP.flp = (float *)AlcMalloc((unsigned long )tI0 *
						     sizeof(float))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		if((dstValP.dbp = (double *)AlcMalloc((unsigned long )tI0 *
						      sizeof(double))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_GREY_RGBA:
		if((dstValP.rgbp = (UINT *)AlcMalloc((unsigned long )tI0 *
						      sizeof(UINT))) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      default:
		errNum = WLZ_ERR_VALUES_TYPE;
		break;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzMakeRect(0, nLines - 1, 0, nCols - 1, vType, dstValP.inp,
		      	       bkgPix, NULL, NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj->values.r->freeptr = AlcFreeStackPush(
	  					     dstObj->values.r->freeptr,
	  					     (void *)dstValP.ubp,
						     NULL);
	  angle = 0;
	  while(nLines-- > 0)
	  {
	    dist = 0;
	    countX = nCols;
	    sinAngle = sin(angle);
	    cosAngle = cos(angle);
	    while(countX-- > 0)
	    {
	      WlzGreyValueGet(gVWSp, 0, (dist * sinAngle) + org.vtY,
	      		      (dist * cosAngle) + org.vtX);
  	      switch(vType)
	      {
	        case WLZ_GREY_INT: 
		  *(dstValP.inp)++ = (*(gVWSp->gVal)).inv;
		  break;
		case WLZ_GREY_SHORT:
		  *(dstValP.shp)++ = (*(gVWSp->gVal)).shv;
		  break;
		case WLZ_GREY_UBYTE:
		  *(dstValP.ubp)++ = (*(gVWSp->gVal)).ubv;
		  break;
		case WLZ_GREY_FLOAT:
		  *(dstValP.flp)++ = (*(gVWSp->gVal)).flv;
		  break;
		case WLZ_GREY_DOUBLE:
		  *(dstValP.dbp)++ = (*(gVWSp->gVal)).dbv;
		  break;
		case WLZ_GREY_RGBA:
		  *(dstValP.rgbp)++ = (*(gVWSp->gVal)).rgbv;
		  break;
	      }
	      dist += distInc;
	    }
	    angle += angleInc;
	  }
	}
	WlzGreyValueFreeWSp(gVWSp);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(dstObj)
      {
	(void )WlzFreeObj(dstObj);
      }
      else if(dstValP.ubp)
      {
	AlcFree(dstValP.ubp);
      }
      dstObj = NULL;
    }
  }
  WLZ_DBG(WLZ_DBG_LVL_1,
	  ("WlzPolarSample 01 %d\n",
	   (int )errNum));
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzPolarSample FX 0x%lx\n",
	   (unsigned long )dstObj));
  return (dstObj);
}
