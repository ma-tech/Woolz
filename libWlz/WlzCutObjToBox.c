#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCutObjToBox_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzCutObjToBox.c
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
* \brief	Functions for creating new domain objects with rectangular
* 		value tables.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

static void			WlzCutObjSetRand(
				  WlzGreyP vec,
				  int vecOff,
				  WlzGreyType gType,
				  int count,
				  double mu,
				  double sigma);

/*!
* \return	New object with rectangular value table or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with a rectangular value table from
*		the given woolz object.
* \param	srcObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dstGreyType		Required grey type for the value table.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox2D(WlzObject *srcObj, WlzIBox2 cutBox,
				  WlzGreyType dstGreyType,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox2D(srcObj, cutBox, dstGreyType, NULL,
  			     bgdNoise, bgdMu, bgdSigma,
			     dstErrNum));
}

/*!
* \return	New object with rectangular value table or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with a rectangular value table from
*               the given woolz object and allows access to grey table.
* \param	srcObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dstGreyType		Required grey type for the value table.
* \param	valP			If non-NULL allocated space for grey
* 					values.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToValBox2D(WlzObject *srcObj, WlzIBox2 cutBox,
				  WlzGreyType dstGreyType, void *valP,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  int		tI0,
		dstOff,
		srcOff,
		dstCount,
		dstBkgCount,
  		lastFilledOff = -1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;
  WlzGreyP	dstValP;
  WlzPixelV	dstBkgPix,
  		srcBkgPix;
  WlzIVertex2	cutSz;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox2D FE 0x%lx {%d %d %d %d} %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj,
	   cutBox.xMin, cutBox.yMin, cutBox.xMax, cutBox.yMax,
	   (int )dstGreyType, (unsigned long )valP,
	   (unsigned long )dstErrNum));
  dstBkgPix.type = dstGreyType;
  dstValP.inp = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
    	    ("WlzCutObjToValBox2D 01 %d\n",
	     (int )(srcObj->type)));
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
	else
	{
	  srcBkgPix = WlzGetBackground(srcObj, NULL);
	  errNum = WlzValueConvertPixel(&dstBkgPix, srcBkgPix, dstGreyType);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  cutSz.vtX = cutBox.xMax - cutBox.xMin + 1;
	  cutSz.vtY = cutBox.yMax - cutBox.yMin + 1;
	  if((cutSz.vtX > 0) && (cutSz.vtY > 0))
	  {
	    tI0 = cutSz.vtX * cutSz.vtY;
	    switch(dstGreyType)
	    {
	      case WLZ_GREY_LONG:
		if(valP)
		{
		  dstValP.lnp = (long *)valP;
		}
		else
		{
		  if((dstValP.lnp = (long *)AlcMalloc((unsigned long )tI0 *
						     sizeof(long))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_INT:
		if(valP)
		{
		  dstValP.inp = (int *)valP;
		}
		else
		{
		  if((dstValP.inp = (int *)AlcMalloc((unsigned long )tI0 *
						     sizeof(int))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_SHORT:
		if(valP)
		{
		  dstValP.shp = (short *)valP;
		}
		else
		{
		  if((dstValP.shp = (short *)AlcMalloc((unsigned long )tI0 *
						       sizeof(short))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_UBYTE:
		if(valP)
		{
		  dstValP.ubp = (WlzUByte *)valP;
		}
		else
		{
		  if((dstValP.ubp = (WlzUByte *)
		                    AlcMalloc((unsigned long )tI0 *
					      sizeof(WlzUByte))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_FLOAT:
		if(valP)
		{
		  dstValP.flp = (float *)valP;
		}
		else
		{
		  if((dstValP.flp = (float *)AlcMalloc((unsigned long )tI0 *
						       sizeof(float))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_DOUBLE:
		if(valP)
		{
		  dstValP.dbp = (double *)valP;
		}
		else
		{
		  if((dstValP.dbp = (double *)AlcMalloc((unsigned long )tI0 *
						      sizeof(double))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_RGBA:
		if(valP)
		{
		  dstValP.rgbp = (WlzUInt *)valP;
		}
		else
		{
		  if((dstValP.rgbp = (WlzUInt *)
		                     AlcMalloc((unsigned long )tI0 *
					       sizeof(WlzUInt))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      default:
	        errNum = WLZ_ERR_VALUES_TYPE;
		break;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzInitGreyScan(srcObj, &iWSp, &gWSp);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
	      {
		if((iWSp.linpos >= cutBox.yMin) &&
		   (iWSp.linpos <= cutBox.yMax) &&
		   (iWSp.lftpos <= cutBox.xMax) &&
		   (iWSp.rgtpos >= cutBox.xMin))
		{
		  dstOff = (iWSp.linpos - cutBox.yMin) * cutSz.vtX;
		  if((tI0 = cutBox.xMin - iWSp.lftpos) > 0)
		  {
		    srcOff = tI0;
		  }
		  else
		  {
		    srcOff = 0;
		    dstOff -= tI0;
		  }
		  dstCount = WLZ_MIN(iWSp.rgtpos, cutBox.xMax) -
			     WLZ_MAX(iWSp.lftpos, cutBox.xMin) + 1;;
		  dstBkgCount = dstOff - lastFilledOff + 1;
		  WLZ_DBG((WLZ_DBG_LVL_3),
			  ("WlzCutObjToValBox2D 02 %d %d %d %d\n",
			   lastFilledOff, dstOff, dstCount, dstBkgCount));
		  if(dstBkgCount > 0)
		  {
		    if(bgdNoise)
		    {
		      WlzCutObjSetRand(dstValP, lastFilledOff + 1,
				       dstGreyType, dstBkgCount,
				       bgdMu, bgdSigma);
		    }
		    else
		    {
		      WlzValueSetGrey(dstValP, lastFilledOff + 1, dstBkgPix.v,
				      dstGreyType, dstBkgCount);
		    }
		    lastFilledOff += dstBkgCount;
		  }
		  if(dstCount > 0)
		  {
		    WlzValueCopyGreyToGrey(dstValP, dstOff,
					   dstGreyType,
					   gWSp.u_grintptr, srcOff,
					   gWSp.pixeltype,
					   dstCount);
		    lastFilledOff = dstOff + dstCount - 1;
		  }
		}
	      }
	      if(errNum == WLZ_ERR_EOO) /* Reset error from end of intervals */ 
	      {
		errNum = WLZ_ERR_NONE;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstBkgCount = (cutSz.vtX * cutSz.vtY) - (lastFilledOff + 1);
	      WLZ_DBG((WLZ_DBG_LVL_3),
		      ("WlzCutObjToValBox2D 03 %d %d\n",
		       lastFilledOff, dstBkgCount));
	      if(dstBkgCount > 0)
	      {
		if(bgdNoise)
		{
		  WlzCutObjSetRand(dstValP, lastFilledOff + 1,
				   dstGreyType, dstBkgCount,
				   bgdMu, bgdSigma);
		}
		else
		{
		  WlzValueSetGrey(dstValP, lastFilledOff + 1, dstBkgPix.v,
				      dstGreyType, dstBkgCount);
	        }
	      }
	      dstObj = WlzMakeRect(cutBox.yMin, cutBox.yMax,
				   cutBox.xMin, cutBox.xMax,
				   dstGreyType, dstValP.inp,
				   dstBkgPix, NULL, NULL, &errNum);
	      if( errNum == WLZ_ERR_NONE )
	      {
		dstObj->values.r->freeptr = AlcFreeStackPush(NULL,
	         				(void *)(dstValP.inp), NULL);
	      }
	    }
	  }
	  else   /* No intersection between the cut box and the given domain */
	  {
	    dstObj = WlzMakeEmpty(&errNum);
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(valP == NULL)
    {
      if(dstValP.inp)
      {
	AlcFree(dstValP.inp);
      }
    }
    if(dstErrNum)
    {
      *dstErrNum = errNum;
    }
    if(dstObj)
    {
      WlzFreeObj(dstObj);
      dstObj = NULL;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox2D FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}


/*!
* \return	New object with rectangular value table(s) or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with (a) rectangular value table(s)
*               from the given woolz object.
* \param	srcObj			Given source object.
* \param	cutBox			Cut box.
* \param	dstGreyType		Required grey type for the value table.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox3D(WlzObject *srcObj, WlzIBox3 cutBox,
				  WlzGreyType dstGreyType,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox3D(srcObj, cutBox, dstGreyType, NULL,
  			     bgdNoise, bgdMu, bgdSigma,
			     dstErrNum));
}

/*!
* \return	New object with rectangular value table(s) or NULL on error.
* \ingroup      WlzValuesUtils
* \brief	Cuts a new object with (a) rectangular value table(s) from the
* 		given woolz object.
* \param	srcObj			Given source object.
* \param	cutBox			Cut box.
* \param	dstGreyType		Required grey type for the value table.
* \param	valP			If non-NULL allocated space for grey
* 					values.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
*                                       may be NULL if not required.
*/
WlzObject	*WlzCutObjToValBox3D(WlzObject *srcObj, WlzIBox3 cutBox,
				  WlzGreyType dstGreyType, void *valP,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  int		size2D,
		size3D,
  		dstPlaneIdx,
  		srcPlaneIdx,
		srcPlanePos,
  		planeCount;
  WlzObject	*dstObj = NULL,
  		*srcObj2D = NULL,
		*dstObj2D = NULL,
		*dummyObj2D = NULL;
  WlzDomain	dstDom,
  		srcDom,
		srcDom2D,
		tDom0;
  WlzValues	dstValues,
  		srcValues2D,
		tVal0;
  WlzGreyP	dstValP,
  		dstVal2DP;
  WlzPixelV	dstBkgPix,
  		srcBkgPix;
  WlzIBox2	cutBox2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox3D FE 0x%lx {%d %d %d %d %d %d} %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj,
	   cutBox.xMin, cutBox.yMin, cutBox.zMin, cutBox.xMax,
	   cutBox.yMax, cutBox.zMax,
	   (int )dstGreyType, (unsigned long )valP,
	   (unsigned long )dstErrNum));
  dstDom.core = NULL;
  dstValues.core = NULL;
  dstValP.inp = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
	    ("WlzCutObjToValBox3D 01 %d\n",
	     (int )(srcObj->type)));
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
        cutBox2D.xMin = cutBox.xMin;
        cutBox2D.xMax = cutBox.xMax;
        cutBox2D.yMin = cutBox.yMin;
        cutBox2D.yMax = cutBox.yMax;
	dstObj = WlzCutObjToValBox2D(srcObj, cutBox2D, dstGreyType,
				  valP, bgdNoise, bgdMu, bgdSigma,
				  &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	if((srcDom = srcObj->domain).core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	    srcBkgPix = WlzGetBackground(srcObj, NULL);
	    errNum = WlzValueConvertPixel(&dstBkgPix, srcBkgPix,
					  dstGreyType);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dummyObj2D = WlzMakeRect(0, 1, 0, 1,
				   srcBkgPix.type, &(srcBkgPix.v.inv),
				   srcBkgPix, NULL, NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((cutBox.xMin <= cutBox.xMax) &&
	     (cutBox.yMin <= cutBox.yMax) &&
	     (cutBox.zMin <= cutBox.zMax))
 	  {
	    size2D = (cutBox.xMax - cutBox.xMin + 1) *
	             (cutBox.yMax - cutBox.yMin + 1);
	    size3D = size2D * (cutBox.zMax - cutBox.zMin + 1);
	    switch(dstGreyType)
	    {
	      case WLZ_GREY_LONG:
		if(valP)
		{
		  dstValP.lnp = (long *)valP;
		}
		else
		{
		  if((dstValP.lnp = (long *)AlcMalloc((unsigned long )size3D *
						     sizeof(long))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_INT:
		if(valP)
		{
		  dstValP.inp = (int *)valP;
		}
		else
		{
		  if((dstValP.inp = (int *)AlcMalloc((unsigned long )size3D *
						     sizeof(int))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_SHORT:
		if(valP)
		{
		  dstValP.shp = (short *)valP;
		}
		else
		{
		  if((dstValP.shp = (short *)AlcMalloc((unsigned long )size3D *
						       sizeof(short))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_UBYTE:
		if(valP)
		{
		  dstValP.ubp = (WlzUByte *)valP;
		}
		else
		{
		  if((dstValP.ubp = (WlzUByte *)
		  		    AlcMalloc((unsigned long )size3D *
					       sizeof(WlzUByte))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_FLOAT:
		if(valP)
		{
		  dstValP.flp = (float *)valP;
		}
		else
		{
		  if((dstValP.flp = (float *)AlcMalloc((unsigned long )size3D *
						       sizeof(float))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_DOUBLE:
		if(valP)
		{
		  dstValP.dbp = (double *)valP;
		}
		else
		{
		  if((dstValP.dbp = (double *)AlcMalloc((unsigned long )size3D *
						      sizeof(double))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      case WLZ_GREY_RGBA:
		if(valP)
		{
		  dstValP.rgbp = (WlzUInt *)valP;
		}
		else
		{
		  if((dstValP.rgbp = (WlzUInt *)
		                     AlcMalloc((unsigned long )size3D *
					       sizeof(WlzUInt))) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		break;
	      default:
		errNum = WLZ_ERR_VALUES_TYPE;
		break;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      srcDom2D.core = NULL;
	      srcValues2D.core = NULL;
	      cutBox2D.xMin = cutBox.xMin;
	      cutBox2D.xMax = cutBox.xMax;
	      cutBox2D.yMin = cutBox.yMin;
	      cutBox2D.yMax = cutBox.yMax;
	      dstDom.p = 
		WlzMakePlaneDomain(srcDom.p->type,
				   cutBox.zMin, cutBox.zMax,
				   cutBox.yMin, cutBox.yMax,
				   cutBox.xMin, cutBox.xMax,
				   &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstValues.vox = 
		WlzMakeVoxelValueTb(srcObj->values.vox->type,
				    cutBox.zMin, cutBox.zMax,
				    dstBkgPix, NULL,
				    &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      srcObj2D =
		WlzMakeMain(WLZ_2D_DOMAINOBJ, srcDom2D,
			    srcValues2D, NULL, NULL, &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstPlaneIdx = 0;
	      srcPlaneIdx = cutBox.zMin - srcDom.p->plane1;
	      srcPlanePos = cutBox.zMin;
	      planeCount = cutBox.zMax - cutBox.zMin + 1;
	      dstVal2DP = dstValP;
	      while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	      {
		/*
		 * Use the dummy 2D object to get a background plane if
		 * the plane is outside of the source objects domain or
		 * it is an old-style NULL plane.
		 */
		if((srcPlanePos < srcDom.p->plane1) ||
		   (srcPlanePos > srcDom.p->lastpl) ||
		   ((tDom0 = *(srcObj->domain.p->domains +
		   	       srcPlaneIdx)).core == NULL) ||
		   ((tVal0 = *(srcObj->values.vox->values +
		   	       srcPlaneIdx)).core == NULL))
		{
		  srcObj2D->domain = dummyObj2D->domain;
		  srcObj2D->values = dummyObj2D->values;
		}
		else
		{
		  srcObj2D->domain = tDom0;
		  srcObj2D->values = tVal0;
		}
		if(((dstObj2D = WlzCutObjToValBox2D(srcObj2D, cutBox2D,
						  dstGreyType,
						  (void *)(dstVal2DP.inp),
						  bgdNoise, bgdMu, bgdSigma,
						  &errNum)) != NULL) &&
		   (errNum == WLZ_ERR_NONE))
		{
		  dstObj2D->values.r->freeptr = NULL;
		  *(dstDom.p->domains + dstPlaneIdx) =
		    WlzAssignDomain(dstObj2D->domain, NULL);
		  *(dstValues.vox->values + dstPlaneIdx) =
		    WlzAssignValues(dstObj2D->values, NULL);
		  ++srcPlanePos;
		  ++srcPlaneIdx;
		  ++dstPlaneIdx;
		  switch(dstGreyType)
		  {
		    case WLZ_GREY_LONG:
		      dstVal2DP.lnp += size2D;
		      break;
		    case WLZ_GREY_INT:
		      dstVal2DP.inp += size2D;
		      break;
		    case WLZ_GREY_SHORT:
		      dstVal2DP.shp += size2D;
		      break;
		    case WLZ_GREY_UBYTE:
		      dstVal2DP.ubp += size2D;
		      break;
		    case WLZ_GREY_FLOAT:
		      dstVal2DP.flp += size2D;
		      break;
		    case WLZ_GREY_DOUBLE:
		      dstVal2DP.dbp += size2D;
		      break;
		    case WLZ_GREY_RGBA:
		      dstVal2DP.rgbp += size2D;
		      break;
		    default:
		      errNum = WLZ_ERR_GREY_TYPE;
		      break;
		  }
		  WlzFreeObj(dstObj2D);
		}
	      }
	      srcObj2D->domain.core = NULL;
	      srcObj2D->values.core = NULL;
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dstValues.vox->freeptr = AlcFreeStackPush(NULL,
	      					      (void *)(dstValP.inp),
						      NULL);
	      dstDom.p->voxel_size[0] = srcObj->domain.p->voxel_size[0];
	      dstDom.p->voxel_size[1] = srcObj->domain.p->voxel_size[1];
	      dstDom.p->voxel_size[2] = srcObj->domain.p->voxel_size[2];
	      WlzStandardPlaneDomain(dstDom.p, dstValues.vox);
	      dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues,
				   NULL, NULL, &errNum);

	    }
	  }
	  else
	  {
	    dstObj = WlzMakeEmpty(&errNum);
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(srcObj2D)
  {
    AlcFree(srcObj2D);
  }
  if(dummyObj2D)
  {
    WlzFreeObj(dummyObj2D);
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
      dstObj = NULL;
    }
    else if(dstValP.inp && (valP == NULL))
    {
      AlcFree(dstValP.inp);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzCutObjToValBox3D 03 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzCutObjToValBox3D FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Set random values from a normal distribution.
* \param	vec			Vector of values.
* \param	vecOff			Offset into the vector.
* \param	gType			Grey type of the values in the vector.
* \param	count			Number of values to set.
* \param	mu			Mean of distribution.
* \param	sigma			Standard deviation of distribution.
*/
static void	WlzCutObjSetRand(WlzGreyP vec, int vecOff, WlzGreyType gType,
				 int count, double mu, double sigma)
{
  double	noise;

  switch(gType)
  {
    case WLZ_GREY_LONG:
      vec.lnp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.lnp)++ = WLZ_CLAMP(noise, INT_MIN, INT_MAX);
      }
      break;
    case WLZ_GREY_INT:
      vec.inp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.inp)++ = WLZ_CLAMP(noise, INT_MIN, INT_MAX);
      }
      break;
    case WLZ_GREY_SHORT:
      vec.shp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.shp)++ = WLZ_CLAMP(noise, SHRT_MIN, SHRT_MAX);
      }
      break;
    case WLZ_GREY_UBYTE:
      vec.ubp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.ubp)++ = WLZ_CLAMP(noise, 0, 255);
      }
      break;
    case WLZ_GREY_FLOAT:
      vec.flp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.flp)++ = WLZ_CLAMP(noise, FLT_MIN, FLT_MAX);
      }
      break;
    case WLZ_GREY_DOUBLE:
      vec.dbp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.dbp)++ = noise;
      }
      break;
    case WLZ_GREY_RGBA:
      vec.inp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.rgbp)++ = WLZ_CLAMP(noise, 0, 0xffffff) + 0xff000000;
      }
      break;
    default:
      break;
  }
}
