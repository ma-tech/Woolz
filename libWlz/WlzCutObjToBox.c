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
* \param	sObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dGreyType		Required grey type for the value table.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox2D(WlzObject *sObj, WlzIBox2 cutBox,
				  WlzGreyType dGreyType,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox2D(sObj, cutBox, dGreyType, NULL,
  			     bgdNoise, bgdMu, bgdSigma,
			     dstErrNum));
}

/*!
* \return	New object with rectangular value table or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with a rectangular value table from
*               the given woolz object and allows access to grey table.
*               If the source object has no values then the destination
*               object will always have foreground value 1 and background
*               value 0, ie the background parameters are ignored.
* \param	sObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dGreyType		Required grey type for the value table.
* \param	gValP			If non-NULL allocated space for grey
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
WlzObject	*WlzCutObjToValBox2D(WlzObject *sObj, WlzIBox2 cutBox,
				  WlzGreyType dGreyType, void *gValP,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  int		tI0,
		dOff,
		sOff,
		dCnt,
		dBkgCnt,
		sHasVal = 0,
  		lastFilledOff = -1;
  size_t	tSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dObj = NULL;
  WlzGreyP	dValP;
  WlzPixelV	dFgdPix,
  		dBkgPix,
		sFgdPix,
  		sBkgPix;
  WlzIVertex2	cutSz;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox2D FE %p {%d %d %d %d} %d %p %p\n",
	   sObj, cutBox.xMin, cutBox.yMin, cutBox.xMax, cutBox.yMax,
	   (int )dGreyType, gValP, dstErrNum));
  dBkgPix.type = dGreyType;
  dValP.inp = NULL;
  if(sObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
    	    ("WlzCutObjToValBox2D 01 %d\n",
	     (int )(sObj->type)));
    switch(sObj->type)
    {
      case WLZ_EMPTY_OBJ:
        dObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_2D_DOMAINOBJ:
	if(sObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if((sObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
		(sObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  if(sObj->values.core == NULL)
	  {
	    sHasVal = 0;
	    sFgdPix.type = WLZ_GREY_UBYTE;
	    sBkgPix.type = WLZ_GREY_UBYTE;
	    sBkgPix.v.ubv = 0;
	    sFgdPix.v.ubv = 1;
	    errNum = WlzValueConvertPixel(&dFgdPix, sFgdPix, dGreyType);
	  }
	  else
	  {
	    sHasVal = 1;
	    sBkgPix = WlzGetBackground(sObj, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzValueConvertPixel(&dBkgPix, sBkgPix, dGreyType);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  cutSz.vtX = cutBox.xMax - cutBox.xMin + 1;
	  cutSz.vtY = cutBox.yMax - cutBox.yMin + 1;
	  if((cutSz.vtX > 0) && (cutSz.vtY > 0))
	  {
	    tSz = cutSz.vtX * cutSz.vtY * WlzGreySize(dGreyType);
	    if(tSz == 0)
	    {
	      errNum = WLZ_ERR_GREY_TYPE;
	    }
	    else
	    {
	      if(gValP)
	      {
	        dValP.v = gValP;
	      }
	      else if((dValP.v = AlcMalloc(tSz)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = (sHasVal == 0)?
	               WlzInitRasterScan(sObj, &iWSp, WLZ_RASTERDIR_ILIC):
	               WlzInitGreyScan(sObj, &iWSp, &gWSp);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sHasVal == 0)
	      {
	        while((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE)
		{
		  if((iWSp.linpos >= cutBox.yMin) &&
		     (iWSp.linpos <= cutBox.yMax) &&
		     (iWSp.lftpos <= cutBox.xMax) &&
		     (iWSp.rgtpos >= cutBox.xMin))
		  {
		    dOff = (iWSp.linpos - cutBox.yMin) * cutSz.vtX;
		    if((tI0 = cutBox.xMin - iWSp.lftpos) > 0)
		    {
		      sOff = tI0;
		    }
		    else
		    {
		      sOff = 0;
		      dOff -= tI0;
		    }
		    dCnt = WLZ_MIN(iWSp.rgtpos, cutBox.xMax) -
			       WLZ_MAX(iWSp.lftpos, cutBox.xMin) + 1;;
		    dBkgCnt = dOff - lastFilledOff + 1;
		    WLZ_DBG((WLZ_DBG_LVL_3),
			    ("WlzCutObjToValBox2D 02 %d %d %d %d\n",
			     lastFilledOff, dOff, dCnt, dBkgCnt));
		    if(dBkgCnt > 0)
		    {
		      WlzValueSetGrey(dValP, lastFilledOff + 1,
				      dBkgPix.v, dGreyType, dBkgCnt);
		      lastFilledOff += dBkgCnt;
		    }
		    if(dCnt > 0)
		    {
		      WlzValueSetGrey(dValP, dOff,
				      dFgdPix.v, dGreyType, dCnt);
		      lastFilledOff = dOff + dCnt - 1;
		    }
		  }
		}
	      }
	      else
	      {
		while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
		{
		  if((iWSp.linpos >= cutBox.yMin) &&
		     (iWSp.linpos <= cutBox.yMax) &&
		     (iWSp.lftpos <= cutBox.xMax) &&
		     (iWSp.rgtpos >= cutBox.xMin))
		  {
		    dOff = (iWSp.linpos - cutBox.yMin) * cutSz.vtX;
		    if((tI0 = cutBox.xMin - iWSp.lftpos) > 0)
		    {
		      sOff = tI0;
		    }
		    else
		    {
		      sOff = 0;
		      dOff -= tI0;
		    }
		    dCnt = WLZ_MIN(iWSp.rgtpos, cutBox.xMax) -
			       WLZ_MAX(iWSp.lftpos, cutBox.xMin) + 1;;
		    dBkgCnt = dOff - lastFilledOff + 1;
		    WLZ_DBG((WLZ_DBG_LVL_3),
			    ("WlzCutObjToValBox2D 03 %d %d %d %d\n",
			     lastFilledOff, dOff, dCnt, dBkgCnt));
		    if(dBkgCnt > 0)
		    {
		      if(bgdNoise)
		      {
			WlzCutObjSetRand(dValP, lastFilledOff + 1,
					 dGreyType, dBkgCnt,
					 bgdMu, bgdSigma);
		      }
		      else
		      {
			WlzValueSetGrey(dValP, lastFilledOff + 1,
					dBkgPix.v, dGreyType, dBkgCnt);
		      }
		      lastFilledOff += dBkgCnt;
		    }
		    if(dCnt > 0)
		    {
		      WlzValueCopyGreyToGrey(dValP, dOff,
					     dGreyType,
					     gWSp.u_grintptr, sOff,
					     gWSp.pixeltype,
					     dCnt);
		      lastFilledOff = dOff + dCnt - 1;
		    }
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
	      dBkgCnt = (cutSz.vtX * cutSz.vtY) - (lastFilledOff + 1);
	      WLZ_DBG((WLZ_DBG_LVL_3),
		      ("WlzCutObjToValBox2D 04 %d %d\n",
		       lastFilledOff, dBkgCnt));
	      if(dBkgCnt > 0)
	      {
		if(sHasVal == 1)
		{
		  WlzValueSetGrey(dValP, lastFilledOff + 1, dBkgPix.v,
				  dGreyType, dBkgCnt);
		}
		else if(bgdNoise)
		{
		  WlzCutObjSetRand(dValP, lastFilledOff + 1,
				   dGreyType, dBkgCnt, bgdMu, bgdSigma);
		}
		else
		{
		  WlzValueSetGrey(dValP, lastFilledOff + 1, dBkgPix.v,
				      dGreyType, dBkgCnt);
	        }
	      }
	      dObj = WlzMakeRect(cutBox.yMin, cutBox.yMax,
				 cutBox.xMin, cutBox.xMax,
				 dGreyType, dValP.inp,
				 dBkgPix, NULL, NULL, &errNum);
	      if((errNum == WLZ_ERR_NONE ) &&
	         (gValP == NULL) && (dValP.v != NULL))
	      {
		dObj->values.r->freeptr =
		  AlcFreeStackPush(dObj->values.r->freeptr, dValP.v, NULL);
	      }
	    }
	  }
	  else   /* No intersection between the cut box and the given domain */
	  {
	    dObj = WlzMakeEmpty(&errNum);
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
    if(gValP == NULL)
    {
      if(dValP.inp)
      {
	AlcFree(dValP.inp);
      }
    }
    if(dstErrNum)
    {
      *dstErrNum = errNum;
    }
    if(dObj)
    {
      WlzFreeObj(dObj);
      dObj = NULL;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox2D FX %p\n",
	   dObj));
  return(dObj);
}


/*!
* \return	New object with rectangular value table(s) or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with (a) rectangular value table(s)
*               from the given woolz object.
*               If the source object has no values then the destination
*               object will always have foreground value 1 and background
*               value 0, ie the background parameters are ignored.
* \param	sObj			Given source object.
* \param	cutBox			Cut box.
* \param	dGreyType		Required grey type for the value table.
* \param	bgdNoise		If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgdMu			Mean of background noise.
* \param	bgdSigma		Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox3D(WlzObject *sObj, WlzIBox3 cutBox,
				  WlzGreyType dGreyType,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox3D(sObj, cutBox, dGreyType, NULL,
  			     bgdNoise, bgdMu, bgdSigma,
			     dstErrNum));
}

/*!
* \return	New object with rectangular value table(s) or NULL on error.
* \ingroup      WlzValuesUtils
* \brief	Cuts a new object with (a) rectangular value table(s) from the
* 		given woolz object.
*               If the source object has no values then the destination
*               object will always have foreground value 1 and background
*               value 0, ie the background parameters are ignored.
* \param	sObj			Given source object.
* \param	cutBox			Cut box.
* \param	dGreyType		Required grey type for the value table.
* \param	gValP			If non-NULL allocated space for grey
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
WlzObject	*WlzCutObjToValBox3D(WlzObject *sObj, WlzIBox3 cutBox,
				  WlzGreyType dGreyType, void *gValP,
				  int bgdNoise, double bgdMu, double bgdSigma,
				  WlzErrorNum *dstErrNum)
{
  int		dPlIdx,
  		sPlIdx,
		sPlPos,
  		plCnt,
		sHasVal = 0;
  size_t	sz2D,
  		sz3D,
		tSz;
  WlzObject	*sObj2D,
		*dObj2D,
		*dObj = NULL;
  WlzDomain	dom2D,
  		dDom,
  		sDom,
		sDom2D;
  WlzValues	val2D,
  		dVal,
  		sVal2D;
  WlzGreyP	dValP,
  		dVal2DP;
  WlzPixelV	dBkgPix,
  		sBkgPix;
  WlzIBox2	cutBox2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox3D FE %p {%d %d %d %d %d %d} %d %p %p\n",
	   sObj, cutBox.xMin, cutBox.yMin, cutBox.zMin, cutBox.xMax,
	   cutBox.yMax, cutBox.zMax, (int )dGreyType, gValP, dstErrNum));
  dDom.core = NULL;
  dVal.core = NULL;
  dValP.inp = NULL;
  cutBox2D.xMin = cutBox.xMin;
  cutBox2D.xMax = cutBox.xMax;
  cutBox2D.yMin = cutBox.yMin;
  cutBox2D.yMax = cutBox.yMax;
  if(sObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    WLZ_DBG((WLZ_DBG_LVL_2),
	    ("WlzCutObjToValBox3D 01 %d\n",
	     (int )(sObj->type)));
    switch(sObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dObj = WlzMakeEmpty(&errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
	dObj = WlzCutObjToValBox2D(sObj, cutBox2D, dGreyType,
				  gValP, bgdNoise, bgdMu, bgdSigma,
				  &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
	if((sDom = sObj->domain).core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(sDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  if(sObj->values.core == NULL)
	  {
	    sHasVal = 0;
	    sBkgPix.type = WLZ_GREY_UBYTE;
	    sBkgPix.v.ubv = 0;
	  }
	  else
	  {
	    sHasVal = 1;
	    sBkgPix = WlzGetBackground(sObj, &errNum);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzValueConvertPixel(&dBkgPix, sBkgPix, dGreyType);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((cutBox.xMin <= cutBox.xMax) &&
	     (cutBox.yMin <= cutBox.yMax) &&
	     (cutBox.zMin <= cutBox.zMax))
 	  {
	    sz2D = (cutBox.xMax - cutBox.xMin + 1) *
	           (cutBox.yMax - cutBox.yMin + 1);
	    sz3D = sz2D * (cutBox.zMax - cutBox.zMin + 1);
	    tSz = sz3D * WlzGreySize(dGreyType);
	    if(tSz == 0)
	    {
	      errNum = WLZ_ERR_GREY_TYPE;
	    }
	    else
	    {
	      if(gValP)
	      {
	        dValP.v = gValP;
	      }
	      else if((dValP.v = AlcMalloc(tSz)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      sDom2D.core = NULL;
	      sVal2D.core = NULL;
	      dDom.p = WlzMakePlaneDomain(sDom.p->type,
					  cutBox.zMin, cutBox.zMax,
					  cutBox.yMin, cutBox.yMax,
					  cutBox.xMin, cutBox.xMax,
					  &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					     cutBox.zMin, cutBox.zMax,
					     dBkgPix, NULL, &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dPlIdx = 0;
	      sPlIdx = cutBox.zMin - sDom.p->plane1;
	      sPlPos = cutBox.zMin;
	      plCnt = cutBox.zMax - cutBox.zMin + 1;
	      dVal2DP.v = dValP.v;
	      while((errNum == WLZ_ERR_NONE) && (plCnt-- > 0))
	      {
		sObj2D = NULL;
		if((sPlPos < sDom.p->plane1) ||
		   (sPlPos > sDom.p->lastpl) ||
		   ((sObj->domain.p->domains + sPlIdx)->core == NULL))
		{
		  dom2D.core = NULL;
		  val2D.core = NULL;
		}
		else
		{
		  dom2D.core = (sObj->domain.p->domains +
		                sPlIdx)->core;
		  val2D.core = (sHasVal == 0)?
		               NULL:
			       (sObj->values.vox->values + sPlIdx)->core;
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2D, val2D,
		                       NULL, NULL, &errNum);
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  dObj2D = WlzCutObjToValBox2D(sObj2D, cutBox2D,
					       dGreyType, dVal2DP.v,
					       bgdNoise, bgdMu, bgdSigma,
					       &errNum);
		}
		(void )WlzFreeObj(sObj2D);
		if(errNum == WLZ_ERR_NONE)
		{
		  *(dDom.p->domains + dPlIdx) =
		    WlzAssignDomain(dObj2D->domain, NULL);
		  *(dVal.vox->values + dPlIdx) =
		    WlzAssignValues(dObj2D->values, NULL);
		  (void )WlzFreeObj(dObj2D);
		  ++sPlPos;
		  ++sPlIdx;
		  ++dPlIdx;
		  switch(dGreyType)
		  {
		    case WLZ_GREY_LONG:
		      dVal2DP.lnp += sz2D;
		      break;
		    case WLZ_GREY_INT:
		      dVal2DP.inp += sz2D;
		      break;
		    case WLZ_GREY_SHORT:
		      dVal2DP.shp += sz2D;
		      break;
		    case WLZ_GREY_UBYTE:
		      dVal2DP.ubp += sz2D;
		      break;
		    case WLZ_GREY_FLOAT:
		      dVal2DP.flp += sz2D;
		      break;
		    case WLZ_GREY_DOUBLE:
		      dVal2DP.dbp += sz2D;
		      break;
		    case WLZ_GREY_RGBA:
		      dVal2DP.rgbp += sz2D;
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
	      if((gValP == NULL) && (dValP.v != NULL))
	      {
		dVal.vox->freeptr = AlcFreeStackPush(dVal.vox->freeptr,
						     dValP.v, NULL);
	      }
	      dDom.p->voxel_size[0] = sObj->domain.p->voxel_size[0];
	      dDom.p->voxel_size[1] = sObj->domain.p->voxel_size[1];
	      dDom.p->voxel_size[2] = sObj->domain.p->voxel_size[2];
	      WlzStandardPlaneDomain(dDom.p, dVal.vox);
	      dObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dDom, dVal,
				   NULL, NULL, &errNum);
	    }
	  }
	  else
	  {
	    dObj = WlzMakeEmpty(&errNum);
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dObj)
    {
      WlzFreeObj(dObj);
      dObj = NULL;
    }
    else if((dValP.v != NULL) && (gValP == NULL))
    {
      AlcFree(dValP.v);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_2),
	  ("WlzCutObjToValBox3D 03 %d\n",
	   (int )errNum));
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzCutObjToValBox3D FX %p\n",
	   dObj));
  return(dObj);
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
	*(vec.lnp)++ = (WlzLong )(WLZ_CLAMP(noise, LLONG_MIN, LLONG_MAX));
      }
      break;
    case WLZ_GREY_INT:
      vec.inp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.inp)++ = (int )(WLZ_CLAMP(noise, INT_MIN, INT_MAX));
      }
      break;
    case WLZ_GREY_SHORT:
      vec.shp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.shp)++ = (short )(WLZ_CLAMP(noise, SHRT_MIN, SHRT_MAX));
      }
      break;
    case WLZ_GREY_UBYTE:
      vec.ubp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.ubp)++ = (WlzUByte )(WLZ_CLAMP(noise, 0, 255));
      }
      break;
    case WLZ_GREY_FLOAT:
      vec.flp += vecOff;
      while(count-- > 0)
      {
	noise = AlgRandNormal(mu, sigma);
	*(vec.flp)++ = (float )(WLZ_CLAMP(noise, FLT_MIN, FLT_MAX));
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
	*(vec.rgbp)++ = (WlzUInt )(WLZ_CLAMP(noise, 0, 0xffffff) + 0xff000000);
      }
      break;
    default:
      break;
  }
}
