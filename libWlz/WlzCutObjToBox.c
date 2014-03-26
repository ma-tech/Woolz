#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCutObjToBox_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCutObjToBox.c
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
* \brief	Functions for creating new domain objects with rectangular
* 		value tables.
* \ingroup	WlzValuesUtils
*/

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

static WlzObject 		*WlzCutObjToValBgBox2D(
				  WlzObject *sObj,
				  WlzIBox2 cutBox,
				  WlzGreyType dGreyType,
				  void *gValP,
				  int pln,
				  int bgNoise,
				  double bgMu,
				  double bgSigma,
				  int forceBg,
				  WlzPixelV forcedBgPix,
				  WlzErrorNum *dstErrNum);
static void			WlzCutObjToBoxFillBg(
				  WlzGreyP vec,
				  WlzGreyType gType,
				  int bgNoise,
				  WlzPixelV bgPix,
				  double mu,
				  double sigma,
				  int vecOff,
				  int vecCnt);
static void			WlzCutObjSetRand(
				  WlzGreyP vec,
				  WlzGreyType gType,
				  double mu,
				  double sigma,
				  int vecOff,
				  int count);

/*!
* \return	New object with rectangular value table or NULL on error.
* \ingroup	WlzValuesUtils
* \brief	Cuts a new object with a rectangular value table from
*		the given woolz object.
* \param	sObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dGreyType		Required grey type for the value table.
* \param	bgNoise			If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgMu			Mean of background noise.
* \param	bgSigma			Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox2D(WlzObject *sObj, WlzIBox2 cutBox,
				  WlzGreyType dGreyType,
				  int bgNoise, double bgMu, double bgSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox2D(sObj, cutBox, dGreyType, NULL,
  			     bgNoise, bgMu, bgSigma,
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
* \param	bgNoise			If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgMu			Mean of background noise.
* \param	bgSigma			Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject 	*WlzCutObjToValBox2D(WlzObject *sObj, WlzIBox2 cutBox,
				  WlzGreyType dGreyType, void *gValP,
				  int bgNoise, double bgMu, double bgSigma,
				  WlzErrorNum *dstErrNum)
{
  WlzPixelV	forcedBgPix;

  forcedBgPix.type = WLZ_GREY_ERROR;
  return(WlzCutObjToValBgBox2D(sObj, cutBox, dGreyType, gValP, 0,
                               bgNoise, bgMu, bgSigma, 0, forcedBgPix,
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
*               However if force background is set then the given value is
*               used.
* \param	sObj			Given source object.
* \param	cutBox			Rectangle box to cut from the object.
* \param	dGreyType		Required grey type for the value table.
* \param	gValP			If non-NULL allocated space for grey
* 					values.
* \param	pln			Position of plane, only used for 3D
* 					value tables.
* \param	bgNoise			If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgMu			Mean of background noise.
* \param	bgSigma			Standard deviation of the background
* 					noise.
* \param	forceBg			Force the background value.
* \param	forcedBgPix		Background value of object.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
static WlzObject *WlzCutObjToValBgBox2D(WlzObject *sObj, WlzIBox2 cutBox,
				  WlzGreyType dGreyType, void *gValP, int pln,
				  int bgNoise, double bgMu, double bgSigma,
				  int forceBg, WlzPixelV forcedBgPix,
				  WlzErrorNum *dstErrNum)
{
  int		dOff,
		dCnt,
		dBgCnt,
		sHasDom = 0,
		sHasVal = 0,
  		dLastOff = -1;
  size_t	tSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dObj = NULL;
  WlzGreyP	dValP;
  WlzPixelV	dFgPix,
  		dBgPix,
		sFgPix,
  		sBgPix;
  WlzIVertex2	cutSz;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCutObjToValBox2D FE %p {%d %d %d %d} %d %p %p\n",
	   sObj, cutBox.xMin, cutBox.yMin, cutBox.xMax, cutBox.yMax,
	   (int )dGreyType, gValP, dstErrNum));
  dBgPix.type = dGreyType;
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
	  sHasDom = 0;
	  if(forceBg != 0)
	  {
	    sBgPix = forcedBgPix;
	    errNum = WlzValueConvertPixel(&dBgPix, sBgPix, dGreyType);
	  }
	}
	else if((sObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
		(sObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  sHasDom = 1;
	  if(sObj->values.core == NULL)
	  {
	    sHasVal = 0;
	    sFgPix.type = WLZ_GREY_UBYTE;
	    sFgPix.v.ubv = 1;
	    errNum = WlzValueConvertPixel(&dFgPix, sFgPix, dGreyType);
	    sBgPix.type = WLZ_GREY_UBYTE;
	    sBgPix.v.ubv = 0;
	  }
	  else
	  {
	    sHasVal = 1;
	    if(forceBg == 0)
	    {
	      sBgPix = WlzGetBackground(sObj, &errNum);
	    }
	  }
	  if(forceBg != 0)
	  {
	    sBgPix = forcedBgPix;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzValueConvertPixel(&dBgPix, sBgPix, dGreyType);
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
	    if((errNum == WLZ_ERR_NONE) && (sHasDom != 0))
	    {
	      errNum = (sHasVal == 0)?
	               WlzInitRasterScan(sObj, &iWSp, WLZ_RASTERDIR_ILIC):
	               WlzInitGreyScan(sObj, &iWSp, &gWSp);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      if(sHasDom == 0)
	      {
	        WlzCutObjToBoxFillBg(dValP,
				      dGreyType,
		                      bgNoise, dBgPix, bgMu, bgSigma,
		                      dLastOff + 1, cutSz.vtX * cutSz.vtY);
	      }
	      else if(sHasVal == 0)
	      {
	        while((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE)
		{
		  int	dLn;

		  dLn = iWSp.linpos - cutBox.yMin;
		  if((dLn >= 0) && (dLn < cutSz.vtY) &&
		     (iWSp.lftpos <= cutBox.xMax) &&
		     (iWSp.rgtpos >= cutBox.xMin))
		  {
		    int		lft,
			    	rgt;

		    lft = iWSp.lftpos - cutBox.xMin;
		    rgt = iWSp.rgtpos - cutBox.xMin;
		    if(lft < 0)
		    {
		      lft = 0;
		    }
		    if(rgt >= cutSz.vtX)
		    {
		      rgt = cutSz.vtX - 1;
		    }
		    dCnt = rgt - lft + 1;
		    dOff = dLn * cutSz.vtX;
		    dBgCnt = dOff + lft - dLastOff - 1;;
		    if(dBgCnt > 0)
		    {
		      WlzCutObjToBoxFillBg(dValP,
					    dGreyType,
					    bgNoise, dBgPix, bgMu, bgSigma,
					    dLastOff + 1, dBgCnt);
                    
		    }
		    if(dCnt > 0)
		    {
		      WlzValueSetGrey(dValP, dOff + lft,
		                      dFgPix.v, dGreyType, dCnt);
		    }
		    dLastOff = dOff + rgt;
		  }
		}
		if(errNum == WLZ_ERR_EOO)
		{
		  errNum = WLZ_ERR_NONE;
		  dBgCnt = cutSz.vtX * cutSz.vtY - dLastOff - 1;
		  if(dBgCnt > 0)
		  {
		    WlzCutObjToBoxFillBg(dValP,
					  dGreyType,
					  bgNoise, dBgPix, bgMu, bgSigma,
		                          dLastOff + 1, dBgCnt);
		  }
		}
	      }
	      else /* (sHasDom !== 0) && sHasVal !== 0) */
	      {
		if(gWSp.tvb)
		{
		  iWSp.plnpos = pln;
		}
	        while((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)
		{
		  int	dLn;

		  dLn = iWSp.linpos - cutBox.yMin;
		  if((dLn >= 0) && (dLn < cutSz.vtY) &&
		     (iWSp.lftpos <= cutBox.xMax) &&
		     (iWSp.rgtpos >= cutBox.xMin))
		  {
		    int		lft,
			    	rgt,
				sKlOff;

		    sKlOff = 0;
		    lft = iWSp.lftpos - cutBox.xMin;
		    rgt = iWSp.rgtpos - cutBox.xMin;
		    if(lft < 0)
		    {
		      sKlOff = -lft;
		      lft = 0;
		    }
		    if(rgt >= cutSz.vtX)
		    {
		      rgt = cutSz.vtX - 1;
		    }
		    dCnt = rgt - lft + 1;
		    dOff = dLn * cutSz.vtX;
		    dBgCnt = dOff + lft - dLastOff - 1;;
		    if(dBgCnt > 0)
		    {
		      WlzCutObjToBoxFillBg(dValP,
					    dGreyType,
					    bgNoise, dBgPix, bgMu, bgSigma,
					    dLastOff + 1, dBgCnt);
                    
		    }
		    if(dCnt > 0)
		    {
		      WlzValueCopyGreyToGrey(dValP, dOff + lft,
					     dGreyType,
					     gWSp.u_grintptr, sKlOff,
					     gWSp.pixeltype,
					     dCnt);
		    }
		    dLastOff = dOff + rgt;
		  }
		}
		(void )WlzEndGreyScan(&gWSp);
		if(errNum == WLZ_ERR_EOO)
		{
		  errNum = WLZ_ERR_NONE;
		  dBgCnt = cutSz.vtX * cutSz.vtY - dLastOff - 1;
		  if(dBgCnt > 0)
		  {
		    WlzCutObjToBoxFillBg(dValP,
					  dGreyType,
					  bgNoise, dBgPix, bgMu, bgSigma,
		                          dLastOff + 1, dBgCnt);
		  }
		}
	      }
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      dObj = WlzMakeRect(cutBox.yMin, cutBox.yMax,
				 cutBox.xMin, cutBox.xMax,
				 dGreyType, dValP.inp,
				 dBgPix, NULL, NULL, &errNum);
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
* \param	bgNoise			If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgMu			Mean of background noise.
* \param	bgSigma			Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
* 					may be NULL if not required.
*/
WlzObject	*WlzCutObjToBox3D(WlzObject *sObj, WlzIBox3 cutBox,
				  WlzGreyType dGreyType,
				  int bgNoise, double bgMu, double bgSigma,
				  WlzErrorNum *dstErrNum)
{
  return(WlzCutObjToValBox3D(sObj, cutBox, dGreyType, NULL,
  			     bgNoise, bgMu, bgSigma,
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
* \param	bgNoise			If zero background value is used
* 					otherwise if non zero the background
*					is set using gaussian noise.
* \param	bgMu			Mean of background noise.
* \param	bgSigma			Standard deviation of the background
* 					noise.
* \param	dstErrNum		Destination pointer for error number,
*                                       may be NULL if not required.
*/
WlzObject	*WlzCutObjToValBox3D(WlzObject *sObj, WlzIBox3 cutBox,
				  WlzGreyType dGreyType, void *gValP,
				  int bgNoise, double bgMu, double bgSigma,
				  WlzErrorNum *dstErrNum)
{
  int		sIsTiled = 0,
		sHasVal = 0;
  size_t	sz2D,
  		sz3D,
		tSz;
  WlzObject	*dObj = NULL;
  WlzDomain	dDom,
  		sDom;
  WlzValues	dVal;
  WlzGreyP	dValP;
  WlzPixelV	dBgPix,
  		sBgPix;
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
				  gValP, bgNoise, bgMu, bgSigma,
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
	    sBgPix.type = WLZ_GREY_UBYTE;
	    sBgPix.v.ubv = 0;
	  }
	  else
	  {
	    sHasVal = 1;
	    sIsTiled = WlzGreyTableIsTiled(sObj->values.core->type);
	    sBgPix = WlzGetBackground(sObj, &errNum);
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzValueConvertPixel(&dBgPix, sBgPix, dGreyType);
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
					     dBgPix, NULL, &errNum);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      int 	dPlIdx,
	      		plCnt;

	      dPlIdx = 0;
	      plCnt = cutBox.zMax - cutBox.zMin + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	      for(dPlIdx = 0; dPlIdx < plCnt; ++dPlIdx)
	      {
		if(errNum == WLZ_ERR_NONE)
		{
		  int	sPlIdx,
			sPlPos;
  	          WlzGreyP dVal2DP;
       		  WlzDomain dom2D;
  		  WlzValues val2D;
		  WlzObject *sObj2D = NULL,
			    *dObj2D = NULL;
		  WlzErrorNum errNum2 = WLZ_ERR_NONE;

		  dom2D.core = NULL;
		  val2D.core = NULL;
	          sPlPos = cutBox.zMin + dPlIdx;
		  sPlIdx = cutBox.zMin - sDom.p->plane1 + dPlIdx;
		  dVal2DP = WlzValueSetGreyP(dValP, dGreyType, dPlIdx * sz2D);
		  if((sPlPos >= sDom.p->plane1) &&
		     (sPlPos <= sDom.p->lastpl) &&
		     ((sObj->domain.p->domains + sPlIdx)->core != NULL))
		  {
		    dom2D.core = (sObj->domain.p->domains + sPlIdx)->core;
		    if(sHasVal)
		    {
		      if(sIsTiled)
		      {
			val2D.t = sObj->values.t;
		      }
		      else
		      {
		        val2D.core = (sObj->values.vox->values + sPlIdx)->core;
		      }
		    }
		  }
		  if(errNum2 == WLZ_ERR_NONE)
		  {
		    sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2D, val2D,
					 NULL, NULL, &errNum2);
		  }
		  if(errNum2 == WLZ_ERR_NONE)
		  {
		    dObj2D = WlzCutObjToValBgBox2D(sObj2D, cutBox2D,
		    				   dGreyType, dVal2DP.v,
						   sPlIdx,
						   bgNoise, bgMu, bgSigma,
						   1, sBgPix, &errNum2);
		  }
		  if(errNum2 == WLZ_ERR_NONE)
		  {
		    *(dDom.p->domains + dPlIdx) =
		      WlzAssignDomain(dObj2D->domain, NULL);
		    *(dVal.vox->values + dPlIdx) =
		      WlzAssignValues(dObj2D->values, NULL);
		  }
		  (void )WlzFreeObj(sObj2D);
		  (void )WlzFreeObj(dObj2D);
#ifdef _OPENMP
#pragma omp critical
		  {
#endif
		    if((errNum == WLZ_ERR_NONE) && (errNum2 != WLZ_ERR_NONE))
		    {
		      errNum = errNum2;
		    }
#ifdef _OPENMP
		  }
#endif
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
* \ingroup	WlzValuesUtils
* \brief	Fills contiguous values with the appropriate background
* 		values, these may be either constant or Gaussian noise.
* \param	vec			Vector of values.
* \param	gType			Grey type of the values in the vector.
* \param	bgNoise			Non-zero if noise is required.
* \param	bgPix			Background picel value (used if
*                                       bgNoise == 0) which must be of the
*                                       correct grey type for the destination.
* \param	mu			Gaussian noise mean (only used if
* 					bgNoise != 0).
* \param	sigma			Gaussian noise variance (only used if
* 					bgNoise != 0).
* \param	vecOff			Offset of first value from the start of
* 					the given vector or values.
* \param	vecCnt			Number of values to fill.
*/
static void	WlzCutObjToBoxFillBg(WlzGreyP vec, WlzGreyType gType,
		                      int bgNoise, WlzPixelV bgPix,
				      double mu, double sigma,
		                      int vecOff, int vecCnt)
{
  if(bgNoise)
  {
    WlzCutObjSetRand(vec, gType, mu, sigma, vecOff, vecCnt);
  }
  else
  {
    WlzValueSetGrey(vec, vecOff, bgPix.v, gType, vecCnt);
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Set random values from a normal distribution.
* \param	vec			Vector of values.
* \param	gType			Grey type of the values in the vector.
* \param	mu			Mean of distribution.
* \param	sigma			Standard deviation of distribution.
* \param	vecOff			Offset into the vector.
* \param	count			Number of values to set.
*/
static void	WlzCutObjSetRand(WlzGreyP vec, WlzGreyType gType,
				 double mu, double sigma,
				 int vecOff, int count)
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
