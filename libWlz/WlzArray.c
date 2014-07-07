#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzArray_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzArray.c
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
* \brief	Functions for converting between domain objects
* 		and arrays.
* \ingroup	WlzArray
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>
#include <Wlz.h>

static WlzErrorNum 		WlzToArrayBit2D(
				  WlzUByte ***dstP,
				  WlzObject *srcObj,
				  WlzIVertex2 size,
				  WlzIVertex2 origin);
static WlzErrorNum 		WlzToArrayBit3D(
				  WlzUByte ****dstP,
				  WlzObject *srcObj,
				  WlzIVertex3 size,
				  WlzIVertex3 origin);
static WlzErrorNum 		WlzToArrayGrey2D(
				  void ***dstP,
				  WlzObject *srcObj,
				  WlzIVertex2 size,
				  WlzIVertex2 origin,
				  int noiseFlag,
				  WlzGreyType dstGreyType);
static WlzErrorNum 		WlzToArrayGrey3D(
				  void ****dstP,
				  WlzObject *srcObj,
				  WlzIVertex3 size,
				  WlzIVertex3 origin,
				  int noiseFlag,
				  WlzGreyType dstGreyType);
static WlzObject 		*WlzFromArrayBit2D(
				  WlzUByte**arrayP,
				  WlzIVertex2 arraySize,
				  WlzIVertex2 arrayOrigin,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzFromArrayGrey2D(
				  void **arrayP,
				  WlzIVertex2 arraySize,
				  WlzIVertex2 arrayOrigin,
				  WlzGreyType dstGreyType,
				  WlzGreyType srcGreyType,
				  double valOffset,
				  double valScale,
				  int clampFlag,
				  int noCopyFlag,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzFromArrayBit3D(
				  WlzUByte ***arrayP,
				  WlzIVertex3 arraySize,
				  WlzIVertex3 arrayOrigin,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzFromArrayGrey3D(
				  void ***arrayP,
				  WlzIVertex3 arraySize,
				  WlzIVertex3 arrayOrigin,
				  WlzGreyType dstGreyType,
				  WlzGreyType srcGreyType,
				  double valOffset,
				  double valScale,
				  int clampFlag,
				  int noCopyFlag,
				  WlzErrorNum *dstErr);

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a  Alc bit array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToBArray2D(WlzIVertex2 *dstSizeArrayDat,
			  WlzUByte ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    dstSizeArrayDat->vtX = (size.vtX + 7) / 8;
    dstSizeArrayDat->vtY = size.vtY;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin, noiseFlag, WLZ_GREY_BIT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a Alc int array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToIArray2D(WlzIVertex2 *dstSizeArrayDat, int ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dstArrayDat == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if(dstSizeArrayDat)
    {
      *dstSizeArrayDat = size;
    }
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
			  noiseFlag, WLZ_GREY_INT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a Alc short array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToSArray2D(WlzIVertex2 *dstSizeArrayDat, short ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
    			  noiseFlag, WLZ_GREY_SHORT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a Alc unsigned byte array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToUArray2D(WlzIVertex2 *dstSizeArrayDat,
			  WlzUByte ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
    			  noiseFlag, WLZ_GREY_UBYTE);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts aAlc float array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToFArray2D(WlzIVertex2 *dstSizeArrayDat, float ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
    			  noiseFlag, WLZ_GREY_FLOAT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a Alc double array from any Woolz 2D domain
*		object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToDArray2D(WlzIVertex2 *dstSizeArrayDat, double ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
    			  noiseFlag, WLZ_GREY_DOUBLE);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a Alc unsigned int array for RGBA values from
* 		any Woolz 2D domain object.		
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToRArray2D(WlzIVertex2 *dstSizeArrayDat,
			  unsigned int ***dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex2 origin, WlzIVertex2 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  size, origin,
    			  noiseFlag, WLZ_GREY_RGBA);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an Alc array from any Woolz 2D domain object.
*		If the destination pointer points to a non-NULL 
*		pointer then it is assumed to be a suitable Alc array.
*		The data are assumed to be within the valid range.
* \param	dstP			Destination pointer (assumed 
*					valid if *dstP is non-NULL).
* \param	srcObj			Given Woolz object.
* \param	size			Size of array.
* \param	origin			Array origin wrt given object.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
* \param	dstGreyType		Destination array data type.
*/
WlzErrorNum	WlzToArray2D(void ***dstP, WlzObject *srcObj,
			     WlzIVertex2 size, WlzIVertex2 origin,
			     int noiseFlag, WlzGreyType dstGreyType)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzToArray2D FE %p  %p {%d %d} {%d %d} %d %d\n",
	   dstP, srcObj, size.vtX, size.vtY, origin.vtX, origin.vtY,
	   noiseFlag, (int )dstGreyType));
  if(dstP == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(srcObj == NULL)
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
  else
  {
    switch(dstGreyType)
    {
      case WLZ_GREY_BIT:
	errNum = WlzToArrayBit2D((WlzUByte ***)dstP,  srcObj, size, origin);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE: /* FALLTHROUGH */
      case WLZ_GREY_RGBA:
	errNum = WlzToArrayGrey2D(dstP,  srcObj, size, origin,
				  noiseFlag, dstGreyType);
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzToArray2D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an Alc bit array from any Woolz 2D domain
*		object's domain. If the destination pointer points to a
*		non-NULL pointer then it is assumed to be a suitable Alc
*		array.
* \param	dstP			Destination pointer (assumed 
*					valid if *dstP is non-NULL).
* \param	srcObj			Given Woolz object.
* \param	size			Size of the array.
* \param	origin			Array origin wrt given object.
*/
static WlzErrorNum WlzToArrayBit2D(WlzUByte ***dstP, WlzObject *srcObj,
				   WlzIVertex2 size, WlzIVertex2 origin)
{
  int		ivY,
		lstY,
		bytWidth;
  WlzUByte	*bitLnP = NULL;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayBit2D FE %p  %p {%d %d} {%d %d}\n",
	   dstP, srcObj, size.vtX, size.vtY, origin.vtX, origin.vtY));

  if(*dstP == NULL)
  {
    if(AlcBit2Malloc(dstP, size.vtY, size.vtX) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitRasterScan(srcObj, &iWSp, WLZ_RASTERDIR_ILIC);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    lstY = -1;
    bytWidth = (size.vtX + 7) / 8;
    while((errNum == WLZ_ERR_NONE) &&
    	  ((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE))
    {
      if((ivY = iWSp.linpos - origin.vtY) >= 0)
      {
        if(ivY >= size.vtY)
	{
	  errNum = WLZ_ERR_EOO;
	}
	else
	{
	  while(lstY < ivY)
	  {
	    /* Clear lines from last to this one. */
	    bitLnP = *(*dstP + ++lstY);
	    (void )memset(bitLnP, 0, bytWidth);
	  }
	  WlzBitLnSetItv(bitLnP,
	  		iWSp.lftpos - origin.vtX, iWSp.rgtpos - origin.vtX,
			size.vtX);
        }
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      /* Clear lines from last to end of array. */
      while(++lstY < size.vtY)
      {
        (void )memset(*(*dstP + lstY), 0, bytWidth);
      }
      errNum = WLZ_ERR_NONE;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzToArrayBit2D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an Alc array from any Woolz 2D domain
*		object's domain. If the destination pointer points to a
*		non-NULL pointer then it is assumed to be a suitable Alc
*		The data are assumed to be within the valid range.
*		array.
* \param	dstP			Destination pointer (assumed 
*					valid if *dstP is non-NULL).
* \param	srcObj			Given Woolz object.
* \param	size			Size of the array.
* \param	origin			Array origin wrt given object.
* \param	noiseFlag		Fill background with random 
*					noise with the same mean and
*					std. dev. as the given object
*					if non-zero.		
* \param	dstGreyType		Destination array data type.
*/
static WlzErrorNum WlzToArrayGrey2D(void ***dstP, WlzObject *srcObj,
				    WlzIVertex2 size, WlzIVertex2 origin,
				    int noiseFlag, WlzGreyType dstGreyType)
{
  int		idY;
  double	noiseMu = 0.0,
  		noiseSigma = 0.0;
  void		*tVP0;
  void		**valPP = NULL;
  WlzObject	*cutObj = NULL;
  WlzGreyP	gValP;
  WlzIBox2	cutBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayGrey2D FE %p  %p {%d %d} {%d %d} %d %d\n",
	   dstP, srcObj, size.vtX, size.vtY, origin.vtX, origin.vtY,
	   noiseFlag, (int )dstGreyType));
  gValP.inp = NULL;
  cutBox.xMin = origin.vtX;
  cutBox.yMin = origin.vtY;
  cutBox.xMax = origin.vtX + size.vtX - 1;
  cutBox.yMax = origin.vtY + size.vtY - 1;
  if(noiseFlag)
  {
    (void )WlzGreyStats(srcObj, NULL, NULL, NULL, NULL, NULL,
			&noiseMu, &noiseSigma, &errNum);
			  
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((cutObj = WlzCutObjToValBox2D(srcObj, cutBox, dstGreyType,
				     (*dstP)? **dstP: NULL,
				     noiseFlag, noiseMu, noiseSigma,
				     &errNum)) != NULL)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	if(cutObj->type == WLZ_2D_DOMAINOBJ)
	{
	  cutObj->values.r->freeptr = AlcFreeStackPop(
					   cutObj->values.r->freeptr,
					   &tVP0, NULL);
	  gValP.inp = (int *)tVP0;
	  cutObj->values.r->values.inp = NULL;
	}
      }
      WlzFreeObj(cutObj);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(*dstP == NULL)
    {
      if(gValP.inp)
      {
	if((valPP = (void **)AlcMalloc((size_t )(size.vtY *
				       sizeof(void *)))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  switch(dstGreyType)
	  {
	    case WLZ_GREY_INT:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.inp;
		gValP.inp += size.vtX;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.shp;
		gValP.shp += size.vtX;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.ubp;
		gValP.ubp += size.vtX;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.flp;
		gValP.flp += size.vtX;
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.dbp;
		gValP.dbp += size.vtX;
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.rgbp;
		gValP.rgbp += size.vtX;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
      }
      *dstP = valPP;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(gValP.inp)
    {
      AlcFree(gValP.inp);
    }
    if(valPP)
    {
      AlcFree(valPP);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzToArrayGrey2D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a bit Alc array from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					bit array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToBArray3D(WlzIVertex3 *dstSizeArrayDat,
			  WlzUByte ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    dstSizeArrayDat->vtX = (size.vtX + 7) / 8;
    dstSizeArrayDat->vtY = size.vtY;
    dstSizeArrayDat->vtZ = size.vtZ;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_BIT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a int Alc array from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					int array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToIArray3D(WlzIVertex3 *dstSizeArrayDat, int ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_INT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a short Alc array from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					short array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToSArray3D(WlzIVertex3 *dstSizeArrayDat, short ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_SHORT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a unsigned byte Alc array from any Woolz 3D
*		domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					unsigned byte array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToUArray3D(WlzIVertex3 *dstSizeArrayDat, 
			  WlzUByte ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_UBYTE);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a float Alc array from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					float array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToFArray3D(WlzIVertex3 *dstSizeArrayDat, float ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_FLOAT);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts a double Alc array from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					double array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToDArray3D(WlzIVertex3 *dstSizeArrayDat, double ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_DOUBLE);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an unsigned int Alc array for RGBA values
* 		from any Woolz 3D domain object.
* \param	dstSizeArrayDat		Destination pointer for
*					array size, may be NULL.
* \param	dstArrayDat		Destination pointer for Alc
*					double array.
* \param	srcObj			Given Woolz object.
* \param	origin			Array origin wrt given object.
* \param	size			Required region size.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
*/
WlzErrorNum WlzToRArray3D(WlzIVertex3 *dstSizeArrayDat,
			  unsigned int ****dstArrayDat,
			  WlzObject *srcObj,
			  WlzIVertex3 origin, WlzIVertex3 size,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, size,
    			  origin, noiseFlag, WLZ_GREY_RGBA);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup 	WlzArray
* \brief	Extracts an Alc array from any Woolz 3D domain object.
*		If the destination pointer points to a non-NULL 
*		pointer then it is assumed to be a suitable Alc array.
*		The data are assumed to be within the valid range.
* \param	dstP
* \param	srcObj
* \param	size
* \param	origin
* \param	noiseFlag
* \param	dstGreyType
* \param	dstP			Destination pointer (assumed 
*					valid if *dstP is non-NULL).
* \param	srcObj			Given Woolz object.	
* \param	size			Size of the array.	
* \param	origin			Array origin wrt given object.
* \param	noiseFlag 		Fill background with random 
*					noise with the same mean and
*					std. dev. as the given object
*					if non-zero.		
* \param	dstGreyType		Destination array data type.
*/
WlzErrorNum	WlzToArray3D(void ****dstP, WlzObject *srcObj,
			     WlzIVertex3 size, WlzIVertex3 origin,
			     int noiseFlag, WlzGreyType dstGreyType)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzToArray3D FE %p  %p {%d %d %d} {%d %d %d} %d %d\n",
	   dstP, srcObj,
	   size.vtX, size.vtY, size.vtZ, origin.vtX, origin.vtY, origin.vtZ,
	   noiseFlag, (int )dstGreyType));
  if(dstP == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(dstGreyType)
    {
      case WLZ_GREY_BIT:
	errNum = WlzToArrayBit3D((WlzUByte ****)dstP,  srcObj, size, origin);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:/* FALLTHROUGH */
      case WLZ_GREY_RGBA:
	errNum = WlzToArrayGrey3D(dstP,  srcObj, size, origin,
				  noiseFlag, dstGreyType);
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzToArray3D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an Alc bit array from any Woolz 3D domain
*               object's domain.
*               If the destination pointer points to a non-NULL
*               pointer then it is assumed to be a suitable Alc array.
* \param	dstP			Destination pointer (assumed
*                                       valid if *dstP is non-NULL).
* \param	obj			Given Woolz object.
* \param	sz			Size of the array.
* \param	og			Array origin wrt given object.
*/
static WlzErrorNum WlzToArrayBit3D(WlzUByte ****dstP, WlzObject *obj,
				   WlzIVertex3 sz, WlzIVertex3 og)
{
  int		idP,
  		plnCnt,
		plnSz,
		lastPIdx;
  WlzDomain	*domains;
  WlzDomain	dom;
  WlzValues	nulVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  nulVal.core = NULL;
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayBit3D FE %p  %p {%d %d %d} {%d %d %d}\n",
	   dstP, obj,
	   sz.vtX, sz.vtY, sz.vtZ, og.vtX, og.vtY, og.vtZ));

  if((dom = obj->domain).core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(dom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((domains = dom.p->domains) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(*dstP == NULL)
    {
      if(AlcBit3Malloc(dstP, sz.vtZ, sz.vtY, sz.vtX) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIVertex2	sz2,
		og2;

    sz2.vtX = sz.vtX;
    sz2.vtY = sz.vtY;
    og2.vtX = og.vtX;
    og2.vtY = og.vtY;
    plnSz = sz.vtY * ((sz.vtX + 7) / 8);
    plnCnt = sz.vtZ;
    lastPIdx = dom.p->lastpl - dom.p->plane1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = 0; idP < plnCnt; ++idP)
    {
      if(errNum == WLZ_ERR_NONE)
      {
        int	pIdx;
	WlzDomain dom2;
        WlzUByte ***dstP2;

	pIdx = idP + og.vtZ - dom.p->plane1;
	dstP2 = (*dstP + idP);
	dom2 = *(domains + pIdx);
	if((pIdx < 0) || (pIdx > lastPIdx) || (dom2.core == NULL))
	{
	  (void )memset(**dstP2, 0, plnSz);
	}
	else
	{
	  WlzObject *obj2;
	  WlzErrorNum errNum2 = WLZ_ERR_NONE;

	  obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2, nulVal, NULL, NULL,
			     &errNum2);
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    errNum2 = WlzToArrayBit2D(dstP2, obj2, sz2, og2);
	  }
#ifdef _OPENMP
#pragma omp critical
	  {
#endif
	    if(errNum2 != WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
#ifdef _OPENMP
	  }
	  WlzFreeObj(obj2);
#endif
	}
      }
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzToArrayBit3D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Extracts an Alc array from any Woolz 3D domain object.
*               If the destination pointer points to a non-NULL
*               pointer then it is assumed to be a suitable Alc array.
*               The data are assumed to be within the valid range.
* \param	dstP			Destination pointer (assumed
*                                       valid if *dstP is non-NULL).
* \param	srcObj			Given Woolz object.
* \param	size			Size of the array.
* \param	origin			Array origin wrt given object.
* \param	noiseFlag		Fill background with random
*                                       noise with the same mean and
*                                       std. dev. as the given object
*                                       if non-zero.
* \param	dstGreyType		Destination array data type.
*/
static WlzErrorNum WlzToArrayGrey3D(void ****dstP, WlzObject *srcObj,
				    WlzIVertex3 size, WlzIVertex3 origin,
				    int noiseFlag, WlzGreyType dstGreyType)
{
  int		idY,
		idZ;
  double	noiseMu = 0.0,
  		noiseSigma = 0.0;
  void		*tVP0;
  void		**valPP = NULL;
  void		***valPPP = NULL;
  WlzObject	*cutObj = NULL;
  WlzGreyP	gValP;
  WlzIBox3	cutBox;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayGrey3D FE %p  %p {%d %d %d} {%d %d %d} %d %d\n",
	   dstP, srcObj,
	   size.vtX, size.vtY, size.vtZ, origin.vtX, origin.vtY, origin.vtZ,
	   noiseFlag, (int )dstGreyType));
  gValP.inp = NULL;
  cutBox.xMin = origin.vtX;
  cutBox.yMin = origin.vtY;
  cutBox.zMin = origin.vtZ;
  cutBox.xMax = origin.vtX + size.vtX - 1;
  cutBox.yMax = origin.vtY + size.vtY - 1;
  cutBox.zMax = origin.vtZ + size.vtZ - 1;
  if(noiseFlag)
  {
    (void )WlzGreyStats(srcObj, NULL, NULL, NULL, NULL, NULL,
			&noiseMu, &noiseSigma, &errNum);
			  
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((cutObj = WlzCutObjToValBox3D(srcObj, cutBox, dstGreyType,
				     (*dstP)? (***dstP): NULL,
				     noiseFlag, noiseMu, noiseSigma,
				     &errNum)) != NULL)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	if(cutObj->type == WLZ_3D_DOMAINOBJ)
	{
	  cutObj->values.vox->freeptr = AlcFreeStackPop(
					      cutObj->values.vox->freeptr,
					      &tVP0, NULL);
	  gValP.inp = (int *)tVP0;
	}
      }
      WlzFreeObj(cutObj);
    }
  }
  if(*dstP == NULL)
  {
    if(gValP.inp)
    {
      if(((valPP = (void **)AlcMalloc((size_t )(size.vtZ * size.vtY *
				      sizeof(void *)))) == NULL) ||
	 ((valPPP = (void ***)AlcMalloc((size_t )(size.vtZ *
					sizeof(void **)))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	switch(dstGreyType)
	{
	  case WLZ_GREY_INT:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.inp;
		gValP.inp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.shp;
		gValP.shp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.ubp;
		gValP.ubp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.flp;
		gValP.flp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.dbp;
		gValP.dbp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    for(idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(idY = 0; idY < size.vtY; ++idY)
	      {
		*(valPP + idY) = gValP.rgbp;
		gValP.rgbp += size.vtX;
	      }
	      *(valPPP + idZ) = valPP;
	      valPP += size.vtY;
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
	*dstP = valPPP;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(gValP.inp)
    {
      AlcFree(gValP.inp);
    }
    if(valPP)
    {
      AlcFree(valPP);
    }
    if(valPPP)
    {
      AlcFree(valPPP);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzToArrayGrey3D FX %d\n",
	   errNum));
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzArray
* \brief	Transforms and/or clamps a rectangle of data values using
* 		a given buffer.
* \param	dstValP			Destination grey pointer.
* \param	srcValP			Source grey pointer.
* \param	bufP			Buffer with space for at least
*                                       row of double grey values.
* \param	rectSize		The size of the destination and
*                                       source, also the row size for
*                                       the buffer.
* \param	dstOffset		Offset from destination ptr.
* \param	srcOffset		Offset from source ptr.
* \param	dstGreyType		Destination grey type.
* \param	srcGreyType		Source grey type.
* \param	valOffset		Offset added to each value.
* \param	valScale		Scale factor by which each
*                                       value is multiplied before
*                                       adding the offset.
* \param	clampFlag		Values are clamped to the
*                                       destination type range if the
*                                       clamp flag is non-zero.
* \param	txFlag			Transform flag, if non-zero
*					values are transformed.
*/
static void	WlzArrayTxRectValues(WlzGreyP dstValP, WlzGreyP srcValP,
				     double *bufP,
				     WlzIVertex2 rectSize,
				     size_t dstOff,
				     size_t srcOff,
				     WlzGreyType dstGreyType,
				     WlzGreyType srcGreyType,
				     double valOff, double valScale,
				     int clampFlag, int txFlag)
{
  int		lineCount,
  		rowCount;
  double	*tDP0;
  WlzGreyP	bufValP;

  bufValP.dbp = bufP;
  lineCount = rectSize.vtY;
  while(lineCount-- > 0)
  {
    WlzValueCopyGreyToGrey(bufValP, 0, WLZ_GREY_DOUBLE,
			   srcValP, srcOff, srcGreyType,
			   rectSize.vtX);
    if(txFlag)
    {
      tDP0 = bufP;
      rowCount = rectSize.vtX;
      while(rowCount-- > 0)
      {
	*tDP0 = (*tDP0 * valScale) + valOff;
	++tDP0;
      }
    }
    if(clampFlag)
    {
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  WlzValueClampDoubleToInt(bufP, rectSize.vtX);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueClampDoubleToShort(bufP, rectSize.vtX);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueClampDoubleToUByte(bufP, rectSize.vtX);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueClampDoubleToFloat(bufP, rectSize.vtX);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueClampDoubleToRGBA(bufP, rectSize.vtX);
	  break;
        default:
	  break;
      }
    }
    WlzValueCopyGreyToGrey(dstValP, dstOff, dstGreyType,
			   bufValP, 0, WLZ_GREY_DOUBLE,
			   rectSize.vtX);
    dstOff += (size_t )(rectSize.vtX);
    srcOff += (size_t )(rectSize.vtX);
  }
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromBArray2D(WlzIVertex2 arraySizeDat,
				 WlzUByte **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_BIT, WLZ_GREY_BIT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromIArray2D(WlzIVertex2 arraySizeDat,
				 int **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
			WLZ_GREY_INT, WLZ_GREY_INT, 0.0, 1.0,
			0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromSArray2D(WlzIVertex2 arraySizeDat,
				 short **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_SHORT, WLZ_GREY_SHORT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromUArray2D(WlzIVertex2 arraySizeDat,
				 WlzUByte **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_UBYTE, WLZ_GREY_UBYTE, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromFArray2D(WlzIVertex2 arraySizeDat,
				 float **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_FLOAT, WLZ_GREY_FLOAT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromDArray2D(WlzIVertex2 arraySizeDat,
				 double **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc
*               array.
*               The data are assumed to be within the valid range.
*               If the noCopyFlag is set (non-zero) then the array data
*               space is used for the onjects values without copying.
*               For this to be valid both the source and destination
*               grey type must be the same.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstGreyType		Destination object grey type.
* \param	srcGreyType		Array data type.
* \param	valOffset		Offset added to each value.
* \param	valScale		Scale factor by which each
*                                       value is multiplied before
*                                       adding the offset.
* \param	clampFlag		Values are clamped to the
*                                       destination type range if the
*                                       clamp flag is non-zero.
* \param	noCopyFlag		Use the array data for the
*                                       Woolz object values in-place.
*					Take care when using this option!
* \param	dstErr			stination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromArray2D(void **arrayP,
				WlzIVertex2 arraySize,
				WlzIVertex2 arrayOrigin,
				WlzGreyType dstGreyType,
				WlzGreyType srcGreyType,
				double valOffset, double valScale,
				int clampFlag, int noCopyFlag,
				WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFromArray2D FE %p  {%d %d} {%d %d} "
	   "%d %d %g %g %d %d %p\n",
	   arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, dstErr));
  if((arrayP == NULL) || (*arrayP == NULL) ||
     (arraySize.vtX <= 0) || (arraySize.vtY <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(dstGreyType)
    {
      case WLZ_GREY_BIT:
        dstObj = WlzFromArrayBit2D((WlzUByte **)arrayP,
				   arraySize, arrayOrigin, &errNum);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:/* FALLTHROUGH */
      case WLZ_GREY_RGBA:
	dstObj = WlzFromArrayGrey2D(arrayP,
				    arraySize, arrayOrigin,
				    dstGreyType, srcGreyType,
				    valOffset, valScale,
				    clampFlag, noCopyFlag,
				    &errNum);
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzFromArray2D FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 2D domain object with domain,
*               but no values from the given Alc bitmap array.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
static WlzObject *WlzFromArrayBit2D(WlzUByte **arrayP,
				    WlzIVertex2 arraySize,
				    WlzIVertex2 arrayOrigin,
				    WlzErrorNum *dstErr)
{
  int		idY;
  WlzObject	*dstObj = NULL;
  WlzDynItvPool	iPool;
  WlzDomain	dstDom;
  WlzValues	dstVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	ivPoolMin = 1024,
  		ivPoolTune = 32;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
          ("%p {%d %d} {%d %d} %p\n",
  	   arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   dstErr));
  dstDom.core = NULL;
  dstDom.core = NULL;
  dstVal.core = NULL;
  iPool.itvBlock = NULL;
  /* Set number of intervals to be allocated in each block, Any number greater
   * than 1/2 line width would do, but the more allocations the less efficient
   * the code and too large a block could waste memory. The magic number
   * for ivPoolTune was found by running this code on some domains that
   * I had and keeping the number of mallocs to around 2 or 3. Search
   * on WLZ_DYNITV_TUNE_MALLOC. */
  iPool.itvsInBlock = ivPoolMin + arraySize.vtX +
  		      (arraySize.vtX * arraySize.vtY / ivPoolTune);
  dstDom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
  				   arrayOrigin.vtY,
				   arrayOrigin.vtY + arraySize.vtY,
				   arrayOrigin.vtX,
				   arrayOrigin.vtX + arraySize.vtX,
				   &errNum);
  idY = 0;
  while((errNum == WLZ_ERR_NONE) && (idY < arraySize.vtY))
  {
    errNum = WlzDynItvLnFromBitLn(dstDom.i, *(arrayP + idY),
    				  arrayOrigin.vtY + idY,
				  arraySize.vtX, &iPool);
    ++idY;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardIntervalDomain(dstDom.i);
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((dstDom.i->line1 == dstDom.i->lastln) &&
       (dstDom.i->kol1 == dstDom.i->lastkl) &&
       ((dstDom.i->intvlines == NULL) ||
        (dstDom.i->intvlines->nintvs == 0)))
    {
      dstObj = WlzMakeEmpty(&errNum);
      (void )WlzFreeIntervalDomain(dstDom.i);
      dstDom.core = NULL;
    }
    else
    {
      /* Copy the domain just incase some Woolz functions expect all intervals
       * to be contiguous. */
      if(errNum == WLZ_ERR_NONE)
      {
	dstObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dstDom, dstVal, NULL, NULL,
			     &errNum);
      }
    }
  }
  if((errNum != WLZ_ERR_NONE) && (dstDom.core != NULL))
  {
    (void )WlzFreeIntervalDomain(dstDom.i);
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzFromArrayBit2D FX %p\n",
	  dstObj));
  return(dstObj);
}

/*!
* \return	New Woolz object.
* \ingroup 	WlzArray
* \brief	Creates a Woolz 2D domain object from the given Alc
*               array.
*               The data are assumed to be within the valid range.
*               If the noCopyFlag is set (non-zero) then the array data
*               space is used for the onjects values without copying.
*               For this to be valid both the source and destination
*               grey type must be the same.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstGreyType		Destination object grey type.
* \param	srcGreyType		Array data type.
* \param	valOffset		Offset added to each value.
* \param	valScale		Scale factor by which each
*                                       value is multiplied before
*                                       adding the offset.
* \param	clampFlag		Values are clamped to the
*                                       destination type range if the
*                                       clamp flag is non-zero.
* \param	noCopyFlag		Use the array data for the
*                                       Woolz object values in-place.
*					Take care when using this option!
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
static WlzObject *WlzFromArrayGrey2D(void **arrayP,
				     WlzIVertex2 arraySize,
				     WlzIVertex2 arrayOrigin,
				     WlzGreyType dstGreyType,
				     WlzGreyType srcGreyType,
				     double valOffset, double valScale,
				     int clampFlag, int noCopyFlag,
				     WlzErrorNum *dstErr)
{
  int		txFlag = 0;
  double	*bufP = NULL;
  size_t	tSz;
  WlzGreyP	dstValP,
  		srcValP;
  WlzPixelV	dstBkgPix;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzFromArrayGrey2D FE %p  {%d %d} {%d %d} "
	   "%d %d %g %g %d %d %p\n",
	   arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (int )dstGreyType, (int )srcGreyType, valOffset, valScale,
	   clampFlag, noCopyFlag, dstErr));
  dstValP.inp = NULL;
  dstBkgPix.type = dstGreyType;
  (void )memset(&(dstBkgPix.v), 0, sizeof(WlzGreyV));
  switch(srcGreyType)
  {
    case WLZ_GREY_INT:
      srcValP.inp = *(int **)arrayP;
      break;
    case WLZ_GREY_SHORT:
      srcValP.shp = *(short **)arrayP;
      break;
    case WLZ_GREY_UBYTE:
      srcValP.ubp = *(WlzUByte **)arrayP;
      break;
    case WLZ_GREY_FLOAT:
      srcValP.flp = *(float **)arrayP;
      break;
    case WLZ_GREY_DOUBLE:
      srcValP.dbp = *(double **)arrayP;
      break;
    case WLZ_GREY_RGBA:
      srcValP.rgbp = *(WlzUInt **)arrayP;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(noCopyFlag)
    {
      if(srcGreyType != dstGreyType)
      {
	errNum = WLZ_ERR_GREY_TYPE;
      }
      else
      {
        dstValP = srcValP;
      }
    }
    else
    {
      tSz = (size_t)(arraySize.vtX) * (size_t )(arraySize.vtY);
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  if((dstValP.inp = (int *)AlcMalloc(tSz * sizeof(int))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  if((dstValP.shp = (short *)AlcMalloc(tSz * sizeof(short))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  if((dstValP.ubp = (WlzUByte *)AlcMalloc(tSz *
					          sizeof(WlzUByte))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  if((dstValP.flp = (float *)AlcMalloc(tSz * sizeof(float))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  if((dstValP.dbp = (double *)AlcMalloc(tSz * sizeof(double))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_RGBA:
	  if((dstValP.rgbp = (WlzUInt *)AlcMalloc(tSz *
						  sizeof(WlzUInt))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzMakeRect(arrayOrigin.vtY,
			 arrayOrigin.vtY + arraySize.vtY - 1,
			 arrayOrigin.vtX,
			 arrayOrigin.vtX + arraySize.vtX - 1,
			 dstGreyType, dstValP.inp,
			 dstBkgPix, NULL, NULL, &errNum);
  }
  if( (errNum == WLZ_ERR_NONE) && (noCopyFlag == 0) )
  {
    dstObj->values.r->freeptr = AlcFreeStackPush(dstObj->values.r->freeptr,
						 (void *)(dstValP.inp),
						 NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    txFlag = (valOffset < -(DBL_EPSILON)) ||
             (valOffset > DBL_EPSILON) ||
    	     (valScale < (1.0 - DBL_EPSILON)) ||
	     (valScale > (1.0 + DBL_EPSILON));
    if(txFlag || clampFlag)
    {
      if((bufP = (double *)AlcMalloc((size_t )(arraySize.vtX) *
      				     sizeof(double))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((txFlag == 0) && (clampFlag == 0))
    {
      if(noCopyFlag == 0)
      {
	WlzValueCopyGreyToGrey(dstValP, 0, dstGreyType,
			       srcValP, 0, srcGreyType,
			       arraySize.vtX * arraySize.vtY);
      }
    }
    else
    {
      WlzArrayTxRectValues(dstValP, srcValP, bufP, arraySize, 0, 0,
			   dstGreyType, srcGreyType, valOffset, valScale,
			   clampFlag, txFlag);
    }
  }
  if(bufP)
  {
    AlcFree(bufP);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstObj)
    {
      WlzFreeObj(dstObj);
    }
    else if(dstValP.inp)
    {
      AlcFree(dstValP.inp);
    }
    dstObj = NULL;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
          ("WlzFromArrayGrey2D FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromBArray3D(WlzIVertex3 arraySizeDat,
				 WlzUByte ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_BIT, WLZ_GREY_BIT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromIArray3D(WlzIVertex3 arraySizeDat,
				 int ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_INT, WLZ_GREY_INT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromSArray3D(WlzIVertex3 arraySizeDat,
				 short ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_SHORT, WLZ_GREY_SHORT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromUArray3D(WlzIVertex3 arraySizeDat,
				 WlzUByte ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_UBYTE, WLZ_GREY_UBYTE, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromFArray3D(WlzIVertex3 arraySizeDat,
				 float ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_FLOAT, WLZ_GREY_FLOAT, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc array.
* \param	arraySizeDat		Dimensions of the array.
* \param	arrayDat		Given Alc array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromDArray3D(WlzIVertex3 arraySizeDat,
				 double ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErr)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
	                0, 0, dstErr));
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc
*               array.
*               The data are assumed to be within the valid range.
*               If the noCopyFlag is set (non-zero) then the array data
*               space is used for the onjects values without copying.
*               For this to be valid both the source and destination
*               grey type must be the same.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstGreyType		Destination object grey type.
* \param	srcGreyType		Array data type.
* \param	valOffset		Offset added to each value.
* \param	valScale		Scale factor by which each
*                                       value is multiplied before
*                                       adding the offset.
* \param	clampFlag		Values are clamped to the
*                                       destination type range if the
*                                       clamp flag is non-zero.
* \param	noCopyFlag		Use the array data for the
*                                       Woolz object values in-place.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
WlzObject	*WlzFromArray3D(void ***arrayP,
				WlzIVertex3 arraySize, WlzIVertex3 arrayOrigin,
				WlzGreyType dstGreyType,
				WlzGreyType srcGreyType,
				double valOffset, double valScale,
				int clampFlag, int noCopyFlag,
				WlzErrorNum *dstErr)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFromArray3D FE %p  {%d %d %d} {%d %d %d} "
	   "%d %d %g %g %d %d %p\n",
	   arrayP, arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, dstErr));

  if((arrayP == NULL) || (*arrayP == NULL) || (**arrayP == NULL) ||
     (arraySize.vtX <= 0) || (arraySize.vtY <= 0) || (arraySize.vtZ <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if(srcGreyType == WLZ_GREY_BIT)
    {
      if(dstGreyType != WLZ_GREY_BIT)
      {
	errNum = WLZ_ERR_GREY_TYPE;
      }
      else
      {
	dstObj = WlzFromArrayBit3D((WlzUByte ***)arrayP, arraySize,
				   arrayOrigin, &errNum);
      }
    }
    else
    {
      dstObj = WlzFromArrayGrey3D(arrayP, arraySize, arrayOrigin,
				  dstGreyType, srcGreyType,
				  valOffset, valScale,
				  clampFlag, noCopyFlag, &errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzFromArray3D FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object without grey values
*               from the given Alc byte packed bit array.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstErr			Destination pointer for error
*                                       number, may be NULL.
*/
static WlzObject *WlzFromArrayBit3D(WlzUByte ***arrayP,
				    WlzIVertex3 arraySize,
				    WlzIVertex3 arrayOrigin,
				    WlzErrorNum *dstErr)
{
  int		plnIdx,
  		plnCnt;
  WlzDomain	dstDom;
  WlzValues	dstVal;
  WlzIVertex2	size2D,
  		origin2D;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*obj2D,
  		*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzFromArrayBit3D FE %p  {%d %d %d} {%d %d %d} %p\n",
	   arrayP, arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   dstErr));
  dstDom.core = NULL;
  dstVal.core = NULL;
  dstDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
  				arrayOrigin.vtZ,
				arrayOrigin.vtZ + arraySize.vtZ - 1,
				arrayOrigin.vtY,
				arrayOrigin.vtY + arraySize.vtY - 1,
				arrayOrigin.vtX, 
				arrayOrigin.vtX + arraySize.vtX - 1,
				&errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    size2D.vtX = arraySize.vtX;
    size2D.vtY = arraySize.vtY;
    origin2D.vtX = arrayOrigin.vtX;
    origin2D.vtY = arrayOrigin.vtY;
    plnIdx =  0;
    plnCnt = arraySize.vtZ;
    while((errNum == WLZ_ERR_NONE) && (plnCnt-- > 0))
    {
      obj2D = WlzFromArrayBit2D(*(arrayP + plnIdx), size2D, origin2D,
      				&errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	switch(obj2D->type)
	{
	  case WLZ_2D_DOMAINOBJ:
            *(dstDom.p->domains + plnIdx) = WlzAssignDomain(obj2D->domain,
	    						    NULL);
	    break;
	  case WLZ_EMPTY_OBJ:
	    (dstDom.p->domains + plnIdx)->core = NULL;
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      if(obj2D)
      {
        WlzFreeObj(obj2D);
      }
      ++plnIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(dstDom.p, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstVal,
    			 NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(dstDom.core)
    {
      (void )WlzFreePlaneDomain(dstDom.p);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
          ("WlzFromArrayBit3D FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	New Woolz object.
* \ingroup	WlzArray
* \brief	Creates a Woolz 3D domain object from the given Alc
*               array. The data are assumed to be within the valid range. If
*               the noCopyFlag is set (non-zero) then the array data space is
*               used for the onjects values without copying. For this to be
*               valid both the source and destination grey type must be the
*               same.
* \param	arrayP			Given Alc array.
* \param	arraySize		Dimensions of the array.
* \param	arrayOrigin		Array origin wrt given object.
* \param	dstGreyType		Destination object grey type.
* \param	srcGreyType		Array data type.
* \param	valOffset		Offset added to each value.
* \param	valScale		Scale factor by which each
*					value is multiplied before adding the
*					offset.
* \param	clampFlag		Values are clamped to the destination
* 					type range if the clamp flag is
*					non-zero.
* \param	noCopyFlag		Use the array data for the Woolz
*					object values in-place. Take care with
*					this option!
* \param	dstErr			Destination pointer for error 
*                                       number, may be NULL.
*/
static WlzObject	*WlzFromArrayGrey3D(void ***arrayP,
					    WlzIVertex3 arraySize,
					    WlzIVertex3 arrayOrigin,
					    WlzGreyType dstGreyType,
					    WlzGreyType srcGreyType,
					    double valOffset, double valScale,
					    int clampFlag, int noCopyFlag,
					    WlzErrorNum *dstErr)
{
  int  		planeCnt,
  		planeIdx,
		planePos,
		txFlag = 0;
  size_t	aSz,
  		planeOff;
  int		*tIP0 = NULL;
  double	*bufP = NULL;
  WlzGreyP	dstValP,
  		srcValP;
  WlzDomain	tDom0,
  		dstDom;
  WlzValues	tVal0,
  		dstValues;
  WlzIVertex2	arraySize2D;
  WlzPixelV	dstBkgPix;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzFromArrayGrey3D FE %p  {%d %d %d} {%d %d %d} "
	   "%d %d %g %g %d %d %p\n",
	   arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale, clampFlag, noCopyFlag, dstErr));
  dstDom.core = NULL;
  dstValues.core = NULL;
  dstValP.inp = NULL;
  dstBkgPix.type = dstGreyType;
  (void )memset(&(dstBkgPix.v), 0, sizeof(WlzGreyV));
  arraySize2D.vtX = arraySize.vtX;
  arraySize2D.vtY = arraySize.vtY;
  switch(srcGreyType)
  {
    case WLZ_GREY_INT:
      srcValP.inp = **(int ***)arrayP;
      break;
    case WLZ_GREY_SHORT:
      srcValP.shp = **(short ***)arrayP;
      break;
    case WLZ_GREY_UBYTE:
      srcValP.ubp = **(WlzUByte ***)arrayP;
      break;
    case WLZ_GREY_FLOAT:
      srcValP.flp = **(float ***)arrayP;
      break;
    case WLZ_GREY_DOUBLE:
      srcValP.dbp = **(double ***)arrayP;
      break;
    case WLZ_GREY_RGBA:
      srcValP.rgbp = **(WlzUInt ***)arrayP;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(noCopyFlag)
    {
      if(srcGreyType != dstGreyType)
      {
	errNum = WLZ_ERR_GREY_TYPE;
      }
      else
      {
        dstValP = srcValP;
      }
    }
    else
    {
      aSz = (size_t )(arraySize.vtX) * (size_t)(arraySize.vtY) *
            (size_t )(arraySize.vtZ);
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  if((dstValP.inp = (int *)AlcMalloc(aSz * sizeof(int))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  if((dstValP.shp = (short *)AlcMalloc(aSz * sizeof(short))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  if((dstValP.ubp = (WlzUByte *)AlcMalloc(aSz *
	                                           sizeof(WlzUByte))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  if((dstValP.flp = (float *)AlcMalloc(aSz * sizeof(float))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  if((dstValP.dbp = (double *)AlcMalloc(aSz * sizeof(double))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_RGBA:
	  if((dstValP.rgbp = (WlzUInt *)
	                     AlcMalloc(aSz * sizeof(WlzUInt))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    txFlag = (valOffset < -(DBL_EPSILON)) ||
    	     (valOffset > DBL_EPSILON) ||
	     (valScale < (1.0 - DBL_EPSILON)) ||
	     (valScale > (1.0 + DBL_EPSILON));
    if(txFlag || clampFlag)
    {
      if((bufP = (double *)AlcMalloc((size_t )(arraySize.vtX) *
					       sizeof(double))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstDom.p =  WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				   arrayOrigin.vtZ,
				   arrayOrigin.vtZ + arraySize.vtZ - 1,
				   arrayOrigin.vtY,
				   arrayOrigin.vtY + arraySize.vtY - 1,
				   arrayOrigin.vtX,
				   arrayOrigin.vtX + arraySize.vtX - 1,
				   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues.vox =  WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
					 arrayOrigin.vtZ,
					 (arrayOrigin.vtZ + arraySize.vtZ - 1),
					 dstBkgPix, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(noCopyFlag == 0)
    {
      dstValues.vox->freeptr = AlcFreeStackPush(dstValues.vox->freeptr,
						(void *)(dstValP.inp),
						NULL);
    }
    planeIdx = 0;
    planeCnt = arraySize.vtZ;
    planePos = arrayOrigin.vtZ;
    planeOff = 0;
    if((txFlag == 0) && (clampFlag == 0))
    {
      if(noCopyFlag == 0)
      {
	 aSz = (size_t )(arraySize.vtX) * (size_t)(arraySize.vtY) *
	       (size_t )(arraySize.vtZ);
	 WlzValueCopyGreyToGrey(dstValP, 0, dstGreyType,
				srcValP, 0, srcGreyType, aSz);
      }
    }
    while((errNum == WLZ_ERR_NONE) && (planeCnt-- > 0))
    {
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  tIP0 = dstValP.inp + planeOff;
	  break;
	case WLZ_GREY_SHORT:
	  tIP0 = (int *)(dstValP.shp + planeOff);
	  break;
	case WLZ_GREY_UBYTE:
	  tIP0 = (int *)(dstValP.ubp + planeOff);
	  break;
	case WLZ_GREY_FLOAT:
	  tIP0 = (int *)(dstValP.flp + planeOff);
	  break;
	case WLZ_GREY_DOUBLE:
	  tIP0 = (int *)(dstValP.dbp + planeOff);
	  break;
	case WLZ_GREY_RGBA:
	  tIP0 = (int *)(dstValP.rgbp + planeOff);
	  break;
        default:
	  break;
      }
      tDom0.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				      arrayOrigin.vtY,
				      arrayOrigin.vtY + arraySize.vtY - 1,
				      arrayOrigin.vtX,
				      arrayOrigin.vtX +
				      arraySize.vtX - 1, &errNum);
      if( errNum == WLZ_ERR_NONE )
      {
	tVal0.r = WlzMakeRectValueTb(WlzGreyTableType(WLZ_GREY_TAB_RECT,
						      dstGreyType, NULL),
				     arrayOrigin.vtY,
				     arrayOrigin.vtY + arraySize.vtY - 1,
				     arrayOrigin.vtX,
				     arraySize.vtX,
				     dstBkgPix,
				     tIP0, &errNum);
	if( (errNum != WLZ_ERR_NONE) && tDom0.i )
	{
	  WlzFreeDomain(tDom0);
	}
      }
      if( errNum == WLZ_ERR_NONE )
      {
	tVal0.r->freeptr = NULL;
	if(txFlag || clampFlag)
	{
	  WlzArrayTxRectValues(dstValP, srcValP, bufP, arraySize2D, 
			       planeOff, planeOff,
			       dstGreyType, srcGreyType, valOffset, valScale,
			       clampFlag, txFlag);
	}
	*(dstDom.p->domains + planeIdx) = WlzAssignDomain(tDom0, &errNum);
	*(dstValues.vox->values + planeIdx) = WlzAssignValues(tVal0,
							      &errNum);
	++planeIdx;
	++planePos;
	planeOff += (size_t )(arraySize.vtX) * (size_t )(arraySize.vtY);
      }
    }
  }
  if(bufP)
  {
    AlcFree(bufP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstDom.p->voxel_size[0] = 1.0;
    dstDom.p->voxel_size[1] = 1.0;
    dstDom.p->voxel_size[2] = 1.0;
    WlzStandardPlaneDomain(dstDom.p, dstValues.vox);
    dstObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dstDom, dstValues, NULL, NULL,
			 &errNum);
  }
  else
  {
    if(dstDom.core)
    {
      WlzFreeDomain(dstDom);
    }
    if(dstValues.core)
    {
      WlzFreeValues(dstValues);
    }
    else if(dstValP.inp)
    {
      AlcFree(dstValP.inp);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
          ("WlzFromArrayGrey3D FX %p\n",
	   dstObj));
  return(dstObj);
}

/*!
* \return	Number of data in array.
* \ingroup 	WlzArray
* \brief	Calculates simple statistics for the given Alc array.
* \param	arrayP			Given 3D Alc array.
* \param	arraySize		Dimensions of the array.
* \param	greyType		Array data type.
* \param	dstMin			Destination ptr for minimum value,
* 					may be NULL.
* \param	dstMax			Destination ptr for maximum value,
*					may be NULL.
* \param	dstSum			Destination ptr for sum of values,
*					may be NULL.
* \param	dstSumSq		Destination ptr for sum of squares of
* 					values, may be NULL.
* \param	dstMean			Destination ptr for mean of values,
*					may be NULL.
* \param	dstStdDev		Destination ptr for standard deviation
*					of values, may be NULL.
*/
int		WlzArrayStats3D(void ***arrayP,
				WlzIVertex3 arraySize,
				WlzGreyType greyType,
				double *dstMin, double *dstMax,
				double *dstSum, double *dstSumSq,
				double *dstMean, double *dstStdDev)
{
  int		tI0,
  		arrayCount = 0;
  WlzGreyP	dataP;
  double	tD0,
		prvMin = 0.0,
		prvMax = 0.0,
		prvSum = 0.0,
		prvSumSq = 0.0,
		prvMean = -1.0,
		prvStdDev = -1.0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats3D FE %p {%d %d %d} "
	   "%p %p %p %p %p %p\n",
	   arrayP, arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   dstMin, dstMax, dstSum, dstSumSq, dstMean, dstStdDev));
  if((arrayP == NULL) ||
     (arraySize.vtX < 0) || (arraySize.vtY < 0) || (arraySize.vtZ < 0) ||
     ((arraySize.vtX == 0) && (arraySize.vtY == 0) && (arraySize.vtZ == 0)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    arrayCount = arraySize.vtX * arraySize.vtY * arraySize.vtZ;
    tI0 = arrayCount;
    switch(greyType)
    {
      case WLZ_GREY_INT:
        dataP.inp = **(int ***)arrayP;
	tD0 = *(dataP.inp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.inp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      case WLZ_GREY_SHORT:
        dataP.shp = **(short ***)arrayP;
	tD0 = *(dataP.shp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.shp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      case WLZ_GREY_UBYTE:
        dataP.ubp = **(WlzUByte ***)arrayP;
	tD0 = *(dataP.ubp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.ubp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      case WLZ_GREY_FLOAT:
        dataP.flp = **(float ***)arrayP;
	tD0 = *(dataP.flp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.flp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      case WLZ_GREY_DOUBLE:
        dataP.dbp = **(double ***)arrayP;
	tD0 = *(dataP.dbp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.dbp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      case WLZ_GREY_RGBA:
        dataP.rgbp = **(WlzUInt ***)arrayP;
	tD0 = *(dataP.rgbp)++;
	prvMin = tD0;
	prvMax = tD0;
	prvSum = tD0;
	prvSumSq = tD0 * tD0;
	while(--tI0 > 0)
	{
	  tD0 = *(dataP.rgbp)++;
	  if(tD0 < prvMin)
	  {
	    prvMin = tD0;
	  }
	  if(tD0 > prvMax)
	  {
	    prvMax = tD0;
	  }
	  prvSum += tD0;
	  prvSumSq += tD0 * tD0;
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    prvMean = prvSum / arrayCount;
    if(arrayCount > 1)
    {
      prvStdDev = sqrt((prvSumSq -
		       (prvSum * prvSum / arrayCount)) / (arrayCount - 1));
    }
    else
    {
      prvStdDev = 0.0;
    }
  }
  else
  {
    arrayCount = 0;
  }
  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzArrayStats3D 01 %f %f %f %f %f %f\n",
	   prvMin, prvMax, prvSum, prvSumSq, prvMean, prvStdDev));
  if(dstMin)
  {
    *dstMin = prvMin;
  }
  if(dstMax)
  {
    *dstMax = prvMax;
  }
  if(dstSum)
  {
    *dstSum = prvSum;
  }
  if(dstSumSq)
  {
    *dstSumSq = prvSumSq;
  }
  if(dstMean)
  {
    *dstMean = prvMean;
  }
  if(dstStdDev)
  {
    *dstStdDev = prvStdDev;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats3D FX %d\n",
	   arrayCount));
  return(arrayCount);
}

/*!
* \return	Number of data in array.
* \ingroup	WlzArray
* \brief	Calculates simple statistics for the given Alc array.
* \param	arrayP			Given 2D Alc array.
* \param	arraySize		Dimensions of the array.
* \param	greyType		Array data type.
* \param	dstMin			Destination ptr for minimum value,
*					may be NULL.
* \param	dstMax			Destination ptr for maximum value,
*					may be NULL.
* \param	dstSum			Destination ptr for sum of values,
*					may be NULL.
* \param	dstSumSq		Destination ptr for sum of squares of
* 					values, may be NULL.
* \param	dstMean			Destination ptr for mean of values,
*					may be NULL.
* \param	dstStdDev		Destination ptr for standard deviation
*					of values, may be NULL.
*/
int		WlzArrayStats2D(void **arrayP,
				WlzIVertex2 arraySize,
				WlzGreyType greyType,
				double *dstMin, double *dstMax,
				double *dstSum, double *dstSumSq,
				double *dstMean, double *dstStdDev)
{
  int		arrayCount = 0;
  WlzIVertex3	arraySize3D;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats2D FE %p {%d %d} "
	   "%p %p %p %p %p %p\n",
	   arrayP, arraySize.vtX, arraySize.vtY,
	   dstMin, dstMax, dstSum, dstSumSq, dstMean, dstStdDev));
  arraySize3D.vtX = arraySize.vtX;
  arraySize3D.vtY = arraySize.vtY;
  arraySize3D.vtZ = 1;
  arrayCount = WlzArrayStats3D(&arrayP, arraySize3D, greyType,
  			       dstMin, dstMax, dstSum,  dstSumSq,
			       dstMean, dstStdDev);
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats2D FX %d\n",
	   arrayCount));
  return(arrayCount);
}

/*!
* \return	Number of data in array.
* \ingroup	WlzArray
* \brief	Calculates simple statistics for the given Alc array.
* \param	arrayP			Given 1D Alc array.
* \param	arraySize		Dimension of the array.
* \param	greyType		Array data type.
* \param	dstMin			Destination ptr for minimum value,
*					may be NULL.
* \param	dstMax			Destination ptr for maximum value,
*					may be NULL.
* \param	dstSum			Destination ptr for sum of values,
*					may be NULL.
* \param	dstSumSq		Destination ptr for sum of squares of
* 					values, may be NULL.
* \param	dstMean			Destination ptr for mean of values,
*					may be NULL.
* \param	dstStdDev		Destination ptr for standard deviation
*					of values, may be NULL.
*/
int		WlzArrayStats1D(void *arrayP,
				int arraySize,
				WlzGreyType greyType,
				double *dstMin, double *dstMax,
				double *dstSum, double *dstSumSq,
				double *dstMean, double *dstStdDev)
{
  int		arrayCount = 0;
  WlzIVertex2	arraySize2D;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats1D FE %p %d "
	   "%p %p %p %p %p %p\n",
	   arrayP, arraySize, dstMin, dstMax, dstSum, dstSumSq,
	   dstMean, dstStdDev));
  arraySize2D.vtX = arraySize;
  arraySize2D.vtY = 1;
  arrayCount = WlzArrayStats2D(&arrayP, arraySize2D, greyType,
  			       dstMin, dstMax, dstSum,
			       dstSumSq, dstMean, dstStdDev);
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzArrayStats1D FX %d\n",
	   arrayCount));
  return(arrayCount);
}

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	WlzArray
* \brief	Converts a 1D bit array as generated by the Java interfaces to
* 		a 2D woolz object.
* \param	arraySizeDat		Vertex giving the width and height of
* 					the corresponding domain.
* \param	bitData			The bit data with the bits set
* 					contiguously.
* \param	arrayOrigin		The origin of the Woolz object.
* \param	dstErr			Destination pointer for error code,
*					may be NULL.
*/
WlzObject *WlzFromBArray1D(
  WlzIVertex2 arraySizeDat,
  WlzUByte *bitData,
  WlzIVertex2 arrayOrigin,
  WlzErrorNum *dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum 	errNum=WLZ_ERR_NONE;
  WlzUByte	**arrayData;
  int		i, j, srcOffset;
  WlzUByte	bitmask;

  /* check input data */
  if( bitData == NULL ){
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if( (arraySizeDat.vtX <= 0) || (arraySizeDat.vtY <= 0) ){
    errNum = WLZ_ERR_PARAM_DATA;
  }
  /* convert the 1D domain to 2D */
  if( errNum == WLZ_ERR_NONE ){
    AlcBit2Calloc(&arrayData, arraySizeDat.vtY, arraySizeDat.vtX);
    srcOffset = 0;
    for(i=0; i < arraySizeDat.vtY; i++){
      for(j=0; j < arraySizeDat.vtX; j++){
	bitmask = (WlzUByte )(1 << (srcOffset % 8));
	if( bitData[srcOffset/8] & bitmask ){
/*	  arrayData[i][j/8] |= bitmask;*/
	  arrayData[i][j/8] |= (1 << (j % 8));
	}
	srcOffset++;
      }
    }
  }
  /* use array conversion to build the woolz object */
  if( errNum == WLZ_ERR_NONE ){
    rtnObj = WlzFromBArray2D(arraySizeDat, arrayData, arrayOrigin, &errNum);
    AlcBit2Free(arrayData);
  }
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/*!
* \return	New object.
* \ingroup	WlzArray
* \brief	Creates a new 2D or 3D domain object from the given 1D array.
*		The object created is rectangular/cuboid with the given
*		origin and size. The background value of the new object
*		is set to zero.
*		The returned object can be free'd using WlzFreeObj(),
*		irespective of whether the noCopy flag is set, but if
*		an object created with the noCopy flag set is free'd the
*		data must be free'd independently after the returned object.
*		The no copy option was originaly intended for efficient
*		output of Woolz objects to a file in non-Woolz applications.
* \param	oType			Required object type, must be 
      					WLZ_2D_DOMAINOBJ or WLZ_3D_DOMAINOBJ.
* \param	sz			Size of the required object, ie the
* 					number of columns, lines and planes.
*					The number of planes is ignored for
					WLZ_2D_DOMAINOBJ objects.
* \param	org			The origin of the object. The plane
*					origin is ignored for WLZ_2D_DOMAINOBJ
*					objects.
* \param	gType			The grey type of the given data and the
*					resulting object.
* \param	gDat			The 1D array of data.
* \param	noCopy			If non-zero then the data are not
* 					copied but are used in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzFromArray1D(WlzObjectType oType,
				WlzIVertex3 sz, WlzIVertex3 org,
				WlzGreyType gType, WlzGreyP gDat,
				int noCopy, WlzErrorNum *dstErr)
{
  int		idP,
  		gSz2;
  WlzObject	*obj = NULL;
  WlzObjectType	gTabType;
  WlzDomain	dom;
  WlzValues	val;
  WlzGreyP	cDat;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cDat.v = NULL;
  dom.core = NULL;
  val.core = NULL;
  bgdV.type = gType;
  switch(gType)
  {
    case WLZ_GREY_INT:
      bgdV.v.inv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(int);
      break;
    case WLZ_GREY_SHORT:
      bgdV.v.shv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(short);
      break;
    case WLZ_GREY_UBYTE:
      bgdV.v.ubv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(WlzUByte);
      break;
    case WLZ_GREY_FLOAT:
      bgdV.v.flv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(float);
      break;
    case WLZ_GREY_DOUBLE:
      bgdV.v.dbv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(double);
      break;
    case WLZ_GREY_RGBA:
      bgdV.v.rgbv = 0;
      gSz2 = sz.vtX * sz.vtY * sizeof(WlzUInt);
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gTabType = WlzGreyTableType(WLZ_GREY_TAB_RECT, gType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(oType)
    {
      case WLZ_2D_DOMAINOBJ:
	cDat.v = NULL;
	if((sz.vtX < 1) || (sz.vtY < 1))
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if(noCopy)
	  {
	    cDat.v = gDat.v;
	  }
	  else
	  {
	    if((cDat.v = AlcMalloc(gSz2)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      (void )memcpy(cDat.v, gDat.v, gSz2);
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					org.vtY, org.vtY + sz.vtY - 1,
					org.vtX, org.vtX + sz.vtX - 1,
					&errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  val.r = WlzMakeRectValueTb(gTabType, 
				     org.vtY, org.vtY + sz.vtY - 1,
				     org.vtX, sz.vtX,
				     bgdV, cDat.v,
				     &errNum);
	  if((val.r != NULL) && (noCopy == 0))
	  {
	    val.r->freeptr = cDat.inp;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  obj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val, NULL, NULL,
			    &errNum);
	}
	/* Clean up for 2D on error. */
	if(obj == NULL)
	{
	  if((val.core == NULL) && (noCopy == 0))
	  {
	    AlcFree(cDat.v);
	  }
	  (void )WlzFreeValues(val);
	  (void )WlzFreeDomain(dom);
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if((sz.vtX < 1) || (sz.vtY < 1) || (sz.vtZ < 1))
	{
	  errNum = WLZ_ERR_PARAM_DATA;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
	  			      org.vtZ, org.vtZ + sz.vtZ - 1,
	  			      org.vtY, org.vtY + sz.vtY - 1,
				      org.vtX, org.vtX + sz.vtX - 1,
				      &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  val.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
	  				 org.vtZ, org.vtZ + sz.vtZ - 1,
					 bgdV, NULL, &errNum);
	}
	idP = 0;
	while((errNum == WLZ_ERR_NONE) && (idP < sz.vtZ))
	{
	  (dom.p->domains + idP)->i =
		WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				org.vtY, org.vtY + sz.vtY - 1,
				org.vtX, org.vtX + sz.vtX - 1,
				&errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(noCopy)
	    {
	     cDat.v = gDat.ubp + (gSz2 * idP);
	    }
	    else
	    {
	      if((cDat.v = AlcMalloc(gSz2)) == NULL)
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		(void )memcpy(cDat.v, (void *)(gDat.ubp + (gSz2 * idP)), gSz2);
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    (val.vox->values + idP)->r = WlzMakeRectValueTb(gTabType, 
				org.vtY, org.vtY + sz.vtY - 1,
				org.vtX, sz.vtX,
				bgdV, cDat.inp,
				&errNum);
	    if(((val.vox->values + idP)->r != NULL) && (noCopy == 0))
	    {
	      (val.vox->values + idP)->r->freeptr = cDat.inp;
	      cDat.v = NULL;
	    }
	  }
	  ++idP;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  obj = WlzMakeMain(oType, dom, val, NULL, NULL,
			    &errNum);
	}
	/* Clean up for 3D on error. */
	if(obj == NULL)
	{
	  if(noCopy == 0)
	  {
	    AlcFree(cDat.v);
	  }
	  (void )WlzFreeVoxelValueTb(val.vox);
	  (void )WlzFreePlaneDomain(dom.p);
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
  return(obj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArray
* \brief	Sets values in the given buffer to those of the given
*		object. The given buffer must be large enough to hold
*		the given bounding box of image values.
*		Only values within the domain of the object are set.
* \param	gP			Grey pointer to allocated buffer.
* \param	gType			Grey type of the buffer.
* \param	gBufBox			Buffer bounding box. The z component
*					is ignored if the object is a 2D
*					object.
* \param	gOffset			Offset into the buffer.
* \param	obj			Given object.
*/
WlzErrorNum	WlzToArray1D(WlzGreyP gP, WlzGreyType gType,
			     WlzIBox3 gBufBox, int gOffset, WlzObject *obj)
{
  int		fst,
  		lst,
		idx,
		skp,
		offset,
		bufStp;
  WlzObject	*obj2D;
  WlzGreyWSpace	gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gP.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        bufStp = gBufBox.xMax - gBufBox.xMin + 1;
	if((errNum = WlzInitGreyScan(obj, &iWSp, &gWSp)) == WLZ_ERR_NONE)
	{
	  while((errNum == WLZ_ERR_NONE) &&
	        (WlzNextGreyInterval(&iWSp) == 0))
	  {
	    if((iWSp.linpos >= gBufBox.yMin) &&
	       (iWSp.linpos <= gBufBox.yMax) &&
	       (iWSp.lftpos <= gBufBox.xMax) &&
	       (iWSp.rgtpos >= gBufBox.xMin))
	    {
	      fst = (iWSp.lftpos < gBufBox.xMin)? gBufBox.xMin: iWSp.lftpos;
	      lst = (iWSp.rgtpos > gBufBox.xMax)? gBufBox.xMax: iWSp.rgtpos;
	      offset = gOffset +
	               ((iWSp.linpos - gBufBox.yMin) * bufStp) +
		       fst - gBufBox.xMin;
	      WlzValueCopyGreyToGrey(gP, offset, gType,
				     gWSp.u_grintptr, fst - iWSp.lftpos,
				     gWSp.pixeltype, lst - fst + 1);
	    }
	  }
	  (void )WlzEndGreyScan(&iWSp, &gWSp);
	  if(errNum == WLZ_ERR_EOO)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	if((obj->domain.p->plane1 <= gBufBox.zMax) &&
	   (obj->domain.p->lastpl >= gBufBox.zMin))
	{
	  bufStp = (gBufBox.xMax - gBufBox.xMin + 1) *
	           (gBufBox.yMax - gBufBox.yMin + 1);
	  if(obj->domain.p->plane1 < gBufBox.zMin)
	  {
	    fst = gBufBox.zMin - obj->domain.p->plane1;
	    skp = -fst;
	  }
	  else
	  {
	    fst = 0;
	    skp = obj->domain.p->plane1 - gBufBox.zMin;
	  }
	  lst = (obj->domain.p->lastpl > gBufBox.zMax)?
	        gBufBox.zMax - obj->domain.p->plane1:
		obj->domain.p->lastpl - obj->domain.p->plane1;
	  idx = fst;
	  while((errNum == WLZ_ERR_NONE) && (idx <= lst))
	  {
	    if((*(obj->domain.p->domains + idx)).core &&
	       (*(obj->values.vox->values + idx)).core)
	    {
	      obj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                          *(obj->domain.p->domains + idx),
				  *(obj->values.vox->values + idx),
				  NULL, NULL, &errNum);
	      if(errNum == WLZ_ERR_NONE)
	      {
		offset = (skp + idx) * bufStp;
		errNum = WlzToArray1D(gP, gType, gBufBox,
		                      gOffset + offset, obj2D);
	      }
	      (void )WlzFreeObj(obj2D);
	    }
	    ++idx;
	  }
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

#ifdef WLZ_ARRAY_TEST_1
main(int argc, char *argv[])
{
  WlzObject	*inObj= NULL,
  		*outObj = NULL;
  WlzUByte	**array = NULL;
  WlzIVertex2	org,
  		size;
  WlzIBox2	bBox;
  FILE		*fP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  fP = fopen("in.wlz", "r");
  inObj = WlzReadObj(fP, &errNum);
  fclose(fP);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr, "%s: Error failed to read(%s)\n",
    		   argv[0], WlzStringFromErrorNum(errNum, NULL));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bBox = WlzBoundingBox2I(inObj, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr, "%s: Error failed to find bounding box(%s)\n",
    		   argv[0], WlzStringFromErrorNum(errNum, NULL));
    exit(errNum);
  }
  else
  {
    org.vtX = bBox.xMin;
    org.vtY = bBox.yMin;
    size.vtX = bBox.xMax - bBox.xMin + 1;
    size.vtY = bBox.yMax - bBox.yMin + 1;
    errNum = WlzToArray2D((void ***)&array, inObj, size, org, 0,
			  WLZ_GREY_BIT);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr, "%s: Error failed to make bit array(%s)\n",
    		   argv[0], WlzStringFromErrorNum(errNum, NULL));
    exit(errNum);
  }
  else
  {
    outObj = WlzFromArray2D((void **)array, size, org,
			  WLZ_GREY_BIT, WLZ_GREY_BIT, 0.0, 1.0,
			  0, 0, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr, "%s: Error failed to make Woolz object (%s)\n",
    		   argv[0], WlzStringFromErrorNum(errNum, NULL));
    exit(errNum);
  }
  else
  {
    fP = fopen("out.wlz", "w");
    errNum = WlzWriteObj(fP, outObj);
    fclose(fP);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )fprintf(stderr, "%s: Error failed to write Woolz object(%s)\n",
    		   argv[0], WlzStringFromErrorNum(errNum, NULL));
    exit(errNum);
  }
  return(0);
}
#endif /* WLZ_ARRAY_TEST_1 */

#ifdef WLZ_ARRAY_TEST_2
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int		main(int argc, char *argv[])
{
  int		idA,
		aSz,
  		option,
		noCopy = 0,
  		ok = 1,
  		usage = 0;
  double	nrm,
  		val;
  WlzObject	*obj = NULL;
  WlzGreyP	arrayP;
  WlzObjectType	oType = WLZ_2D_DOMAINOBJ;
  WlzIVertex3	*vP;
  WlzIVertex3	idV,
  		dst,
		org,
  		sz;
  char		*outObjFileStr;
  FILE		*fP;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "23hno:f:s:",
  		outObjFileStrDef[] = "-";

  
  arrayP.v = NULL;
  org.vtX = org.vtY = org.vtZ = 0;
  sz.vtX = sz.vtY = sz.vtZ = 64;
  outObjFileStr = outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        oType = WLZ_2D_DOMAINOBJ;
	break;
      case '3':
        oType = WLZ_3D_DOMAINOBJ;
	break;
      case 'n':
        noCopy = 1;
	break;
      case 'o':
	outObjFileStr = optarg;
	break;
      case 'f':  /* FALLTHROUGH */
      case 's':
	vP = (option == 'f')? &org: &sz;
	if(optarg)
	{
	  (void )sscanf(optarg, "%d,%d,%d",
	                &(vP->vtX), &(vP->vtY), &(vP->vtZ));
	}
        break;
      case 'h':  /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(optind < argc)
  {
    usage = 1;
  }
  ok = usage == 0;
  if(ok)
  {
    aSz = sz.vtX * sz.vtY;
    if(oType == WLZ_3D_DOMAINOBJ)
    {
      aSz *= sz.vtZ;
    }
    arrayP.v = AlcMalloc(aSz * sizeof(WlzUByte));
    if(arrayP.v == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )fprintf(stderr, "%s: failed to allocate test array\n", *argv);
    }
  }
  if(ok)
  {
    /* Set the grey values so that they decay away from the centre of
     * the object. */
    if(oType == WLZ_3D_DOMAINOBJ)
    {
      idA = 0;
      nrm = sz.vtX + sz.vtY + sz.vtZ;
      nrm = 255 * 64 / (nrm * nrm);
      for(idV.vtZ = 0; idV.vtZ < sz.vtZ; ++idV.vtZ)
      {
	dst.vtZ = idV.vtZ - (sz.vtZ / 2);
	dst.vtZ *= dst.vtZ;
	for(idV.vtY = 0; idV.vtY < sz.vtY; ++idV.vtY)
	{
	  dst.vtY = idV.vtY - (sz.vtY / 2);
	  dst.vtY *= dst.vtY;
	  for(idV.vtX = 0; idV.vtX < sz.vtX; ++idV.vtX)
	  {
	    dst.vtX = idV.vtX - (sz.vtX / 2);
	    dst.vtX *= dst.vtX;
	    val = (dst.vtZ + dst.vtY + dst.vtX) * nrm;
	    val = WLZ_CLAMP(val, 0, 255);
	    *(arrayP.ubp + idA++) = 255 - WLZ_NINT(val);
	  }
	}
      }
    }
    else /* oType == WLZ_2D_DOMAINOBJ */
    {
      idA = 0;
      nrm = sz.vtX + sz.vtY;
      nrm = 255 * 16 / (nrm * nrm);
      for(idV.vtY = 0; idV.vtY < sz.vtY; ++idV.vtY)
      {
	dst.vtY = idV.vtY - (sz.vtY / 2);
	dst.vtY *= dst.vtY;
	for(idV.vtX = 0; idV.vtX < sz.vtX; ++idV.vtX)
	{
	  dst.vtX = idV.vtX - (sz.vtX / 2);
	  dst.vtX *= dst.vtX;
	  val = (dst.vtY + dst.vtX) * nrm;
	  val = WLZ_CLAMP(val, 0, 255);
	  *(arrayP.ubp + idA++) = 255 - WLZ_NINT(val);
	}
      }
    }
    obj = WlzFromArray1D(oType, sz, org, WLZ_GREY_UBYTE, arrayP,
                         noCopy, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      (void )fprintf(stderr, "%s: Error failed to make Woolz object (%s)\n",
		     *argv, WlzStringFromErrorNum(errNum, NULL));
      ok = 0;
    }
  }
  (void )WlzFreeObj(obj);
  AlcFree(arrayP.v);
  if(ok)
  {
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
              stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to test object\n",
                     *argv);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s\n",
    *argv,
    " [-h] [-o<out object>]\n"
    "        [-2] [-3] [-n] [-s<x>,<y>,<z>] [-f<x>,<y>,<z>]\n"
    "Test for WlzFromArray1D() which writes out a grey valued 2D or 3D\n"
    "object. The origin and size may be set using the command line flags\n"
    "The grey values of the object are 255 at it's centre and decay to 0\n"
    "at the boundary of the object.\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output object file name.\n"
    "  -2  2D object.\n"
    "  -3  3D object.\n"
    "  -n  Object shares data (no copy).\n"
    "  -f  Object origin (offset).\n"
    "  -s  Object size.\n");
  }
  return((ok)? 0: 1);
}
#endif /* WLZ_ARRAY_TEST_2 */
