#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzArray.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for conversion between Woolz domain objects
*		and arrays.
* \ingroup	WlzArray
* \todo         -
* \bug          None known.
*/
#pragma ident "MRC HGU $Id$"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

/* #define WLZ_ARRAY_TEST */

static WlzErrorNum WlzToArrayBit2D(UBYTE ***dstP, WlzObject *srcObj,
				   WlzIVertex2 size, WlzIVertex2 origin);
static WlzErrorNum WlzToArrayBit3D(UBYTE ****dstP, WlzObject *srcObj,
				   WlzIVertex3 size, WlzIVertex3 origin);
static WlzErrorNum WlzToArrayGrey2D(void ***dstP, WlzObject *srcObj,
				    WlzIVertex2 size, WlzIVertex2 origin,
				    int noiseFlag, WlzGreyType dstGreyType);
static WlzErrorNum WlzToArrayGrey3D(void ****dstP, WlzObject *srcObj,
				    WlzIVertex3 size, WlzIVertex3 origin,
				    int noiseFlag, WlzGreyType dstGreyType);
static WlzObject *WlzFromArrayBit2D(UBYTE**arrayP,
				    WlzIVertex2 arraySize,
				    WlzIVertex2 arrayOrigin,
				    WlzErrorNum *dstErr);
static WlzObject *WlzFromArrayGrey2D(void **arrayP,
				     WlzIVertex2 arraySize,
				     WlzIVertex2 arrayOrigin,
				     WlzGreyType dstGreyType,
				     WlzGreyType srcGreyType,
				     double valOffset, double valScale,
				     int clampFlag, int noCopyFlag,
				     WlzErrorNum *dstErr);
static WlzObject *WlzFromArrayBit3D(UBYTE ***arrayP,
				    WlzIVertex3 arraySize,
				    WlzIVertex3 arrayOrigin,
				    WlzErrorNum *dstErr);
static WlzObject *WlzFromArrayGrey3D(void ***arrayP,
				     WlzIVertex3 arraySize,
				     WlzIVertex3 arrayOrigin,
				     WlzGreyType dstGreyType,
				     WlzGreyType srcGreyType,
				     double valOffset, double valScale,
				     int clampFlag, int noCopyFlag,
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
WlzErrorNum WlzToBArray2D(WlzIVertex2 *dstSizeArrayDat, UBYTE ***dstArrayDat,
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

  if((dstArrayDat == NULL) || (dstSizeArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstSizeArrayDat = size;
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
WlzErrorNum WlzToUArray2D(WlzIVertex2 *dstSizeArrayDat, UBYTE ***dstArrayDat,
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
  	  ("WlzToArray2D FE 0x%lx  0x%lx {%d %d} {%d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, origin.vtX, origin.vtY,
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
	errNum = WlzToArrayBit2D((UBYTE ***)dstP,  srcObj, size, origin);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	errNum = WlzToArrayGrey2D(dstP,  srcObj, size, origin,
				  noiseFlag, dstGreyType);
	break;
      deafult:
	WLZ_ERR_GREY_TYPE;
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
static WlzErrorNum WlzToArrayBit2D(UBYTE ***dstP, WlzObject *srcObj,
				   WlzIVertex2 size, WlzIVertex2 origin)
{
  int		ivY,
		lstY,
		bytWidth;
  UBYTE		*bitLnP;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayBit2D FE 0x%lx  0x%lx {%d %d} {%d %d}\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, origin.vtX, origin.vtY));

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
  	  ("WlzToArrayGrey2D FE 0x%lx  0x%lx {%d %d} {%d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, origin.vtX, origin.vtY,
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
	if((valPP = (void **)AlcMalloc((unsigned long )(size.vtY *
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
WlzErrorNum WlzToBArray3D(WlzIVertex3 *dstSizeArrayDat, UBYTE ****dstArrayDat,
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
WlzErrorNum WlzToUArray3D(WlzIVertex3 *dstSizeArrayDat, UBYTE ****dstArrayDat,
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
  	  ("WlzToArray3D FE 0x%lx  0x%lx {%d %d %d} {%d %d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
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
	errNum = WlzToArrayBit3D((UBYTE ****)dstP,  srcObj, size, origin);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	errNum = WlzToArrayGrey3D(dstP,  srcObj, size, origin,
				  noiseFlag, dstGreyType);
	break;
      deafult:
	WLZ_ERR_GREY_TYPE;
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
* \param	srcObj			Given Woolz object.
* \param	size			Size of the array.
* \param	origin			Array origin wrt given object.
*/
static WlzErrorNum WlzToArrayBit3D(UBYTE ****dstP, WlzObject *srcObj,
				   WlzIVertex3 size, WlzIVertex3 origin)
{
  int		plnIdx,
  		plnCnt,
		plnSz;
  WlzDomain	*srcDomains;
  UBYTE		***dstP2D;
  WlzIVertex2	size2D,
  		origin2D;
  WlzDomain	srcDom,
  		dumDom;
  WlzValues	dumVal;
  WlzObject	*srcObj2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzToArrayBit3D FE 0x%lx  0x%lx {%d %d %d} {%d %d %d}\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, size.vtZ, origin.vtX, origin.vtY, origin.vtZ));

  if((srcDom = srcObj->domain).core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((srcDomains = srcDom.p->domains) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(*dstP == NULL)
    {
      if(AlcBit3Malloc(dstP, size.vtZ, size.vtY, size.vtX) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dumDom.core = NULL;
    dumVal.core = NULL;
    srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dumDom, dumVal,
    			   NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    size2D.vtX = size.vtX;
    size2D.vtY = size.vtY;
    origin2D.vtX = origin.vtX;
    origin2D.vtY = origin.vtY;
    plnIdx =  0;
    plnSz = size.vtY * ((size.vtX + 7) / 8);
    plnCnt = srcDom.p->lastpl - srcDom.p->plane1 + 1;
    while((errNum == WLZ_ERR_NONE) && (plnCnt-- > 0))
    {
      dstP2D = (*dstP + plnIdx);
      srcObj2D->domain = *(srcDomains + plnIdx);
      if(srcObj2D->domain.core == NULL)
      {
        (void )memset(**dstP2D, 0, plnSz);
      }
      else
      {
        errNum = WlzToArrayBit2D(dstP2D, srcObj2D, size2D, origin2D);
      }
      ++plnIdx;
    }
  }
  srcObj2D->domain = dumDom;
  WlzFreeObj(srcObj2D);
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
  	  ("WlzToArrayGrey3D FE 0x%lx  0x%lx {%d %d %d} {%d %d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
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
      if(((valPP = (void **)AlcMalloc((unsigned long )(size.vtZ * size.vtY *
				      sizeof(void *)))) == NULL) ||
	 ((valPPP = (void ***)AlcMalloc((unsigned long )(size.vtZ *
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
* \return	<void>
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
				     int dstOffset,
				     int srcOffset,
				     WlzGreyType dstGreyType,
				     WlzGreyType srcGreyType,
				     double valOffset, double valScale,
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
			   srcValP, srcOffset, srcGreyType,
			   rectSize.vtX);
    if(txFlag)
    {
      tDP0 = bufP;
      rowCount = rectSize.vtX;
      while(rowCount-- > 0)
      {
	*tDP0 = (*tDP0 * valScale) + valOffset;
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
      }
    }
    WlzValueCopyGreyToGrey(dstValP, dstOffset, dstGreyType,
			   bufValP, 0, WLZ_GREY_DOUBLE,
			   rectSize.vtX);
    dstOffset += rectSize.vtX;
    srcOffset += rectSize.vtX;
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
				 UBYTE **arrayDat,
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
				 UBYTE **arrayDat,
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
  	  ("WlzFromArray2D FE 0x%lx  {%d %d} {%d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErr));
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
        dstObj = WlzFromArrayBit2D((UBYTE **)arrayP,
				   arraySize, arrayOrigin, &errNum);
	break;
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_INT:   /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	dstObj = WlzFromArrayGrey2D(arrayP,
				    arraySize, arrayOrigin,
				    dstGreyType, srcGreyType,
				    valOffset, valScale,
				    clampFlag, noCopyFlag,
				    &errNum);
	break;
      deafult:
	WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzFromArray2D FX 0x%lx\n",
	   (unsigned long )dstObj));
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
static WlzObject *WlzFromArrayBit2D(UBYTE **arrayP,
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
          ("0x%lx {%d %d} {%d %d} 0x%lx\n",
  	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (unsigned long )dstErr));
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
				   arrayOrigin.vtY + arraySize.vtY - 1,
				   arrayOrigin.vtX,
				   arrayOrigin.vtX + arraySize.vtX - 1,
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
  	  ("WlzFromArrayBit2D FX 0x%lx\n",
	  (unsigned long )dstObj));
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
  unsigned long tUL0;
  double	*bufP = NULL;
  WlzGreyP	dstValP,
  		srcValP;
  WlzPixelV	dstBkgPix;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
  	  ("WlzFromArrayGrey2D FE 0x%lx  {%d %d} {%d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErr));
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
      srcValP.ubp = *(UBYTE **)arrayP;
      break;
    case WLZ_GREY_FLOAT:
      srcValP.flp = *(float **)arrayP;
      break;
    case WLZ_GREY_DOUBLE:
      srcValP.dbp = *(double **)arrayP;
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
      tUL0 = (unsigned long)(arraySize.vtX * arraySize.vtY);
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  if((dstValP.inp = (int *)AlcMalloc(tUL0 *
					     sizeof(int))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  if((dstValP.shp = (short *)AlcMalloc(tUL0 *
					       sizeof(short))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  if((dstValP.ubp = (UBYTE *)AlcMalloc(tUL0 *
					       sizeof(UBYTE))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  if((dstValP.flp = (float *)AlcMalloc(tUL0 *
					       sizeof(float))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  if((dstValP.dbp = (double *)AlcMalloc(tUL0 *
						sizeof(double))) == NULL)
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
      if((bufP = (double *)AlcMalloc((unsigned long )(arraySize.vtX) *
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
          ("WlzFromArrayGrey2D FX 0x%lx\n",
	   (unsigned long )dstObj));
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
				 UBYTE ***arrayDat,
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
				 UBYTE ***arrayDat,
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
  	  ("WlzFromArray3D FE 0x%lx  {%d %d %d} {%d %d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErr));

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
	dstObj = WlzFromArrayBit3D((UBYTE ***)arrayP, arraySize, arrayOrigin,
				   &errNum);
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
          ("WlzFromArray3D FX 0x%lx\n",
	   (unsigned long )dstObj));
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
static WlzObject *WlzFromArrayBit3D(UBYTE ***arrayP,
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
  	  ("WlzFromArrayBit3D FE 0x%lx  {%d %d %d} {%d %d %d} 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (unsigned long )dstErr));
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
          ("WlzFromArrayBit3D FX 0x%lx\n",
	   (unsigned long )dstObj));
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
  int  		planeCount,
  		planeIdx,
		planeOffset,
		planePos,
		txFlag = 0;
  size_t	aSz;
  int		*tIP0;
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
  	  ("WlzFromArrayGrey3D FE 0x%lx  {%d %d %d} {%d %d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErr));
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
      srcValP.ubp = **(UBYTE ***)arrayP;
      break;
    case WLZ_GREY_FLOAT:
      srcValP.flp = **(float ***)arrayP;
      break;
    case WLZ_GREY_DOUBLE:
      srcValP.dbp = **(double ***)arrayP;
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
      aSz = arraySize.vtX * arraySize.vtY * arraySize.vtZ;
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
	  if((dstValP.ubp = (UBYTE *)AlcMalloc(aSz * sizeof(UBYTE))) == NULL)
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
      if((bufP = (double *)AlcMalloc((unsigned long )(arraySize.vtX *
      				                     sizeof(double)))) == NULL)
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
					 (arrayOrigin.vtZ + arraySize.vtZ
					  - 1),
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
    planeCount = arraySize.vtZ;
    planePos = arrayOrigin.vtZ;
    planeOffset = 0;
    if((txFlag == 0) && (clampFlag == 0))
    {
      if(noCopyFlag == 0)
      {
	 WlzValueCopyGreyToGrey(dstValP, 0, dstGreyType,
				srcValP, 0, srcGreyType,
				arraySize.vtX * arraySize.vtY * arraySize.vtZ);
      }
    }
    while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
    {
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  tIP0 = dstValP.inp + planeOffset;
	  break;
	case WLZ_GREY_SHORT:
	  tIP0 = (int *)(dstValP.shp + planeOffset);
	  break;
	case WLZ_GREY_UBYTE:
	  tIP0 = (int *)(dstValP.ubp + planeOffset);
	  break;
	case WLZ_GREY_FLOAT:
	  tIP0 = (int *)(dstValP.flp + planeOffset);
	  break;
	case WLZ_GREY_DOUBLE:
	  tIP0 = (int *)(dstValP.dbp + planeOffset);
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
			       planeOffset, planeOffset,
			       dstGreyType, srcGreyType, valOffset, valScale,
			       clampFlag, txFlag);
	}
	*(dstDom.p->domains + planeIdx) = WlzAssignDomain(tDom0, &errNum);
	*(dstValues.vox->values + planeIdx) = WlzAssignValues(tVal0,
							      &errNum);
	++planeIdx;
	++planePos;
	planeOffset += arraySize.vtX * arraySize.vtY;
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
          ("WlzFromArrayGrey3D FX 0x%lx\n",
	   (unsigned long )dstObj));
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
	  ("WlzArrayStats3D FE 0x%lx {%d %d %d} "
	   "0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   (unsigned long )dstMin, (unsigned long )dstMax,
	   (unsigned long )dstSum, (unsigned long )dstSumSq,
	   (unsigned long )dstMean, (unsigned long )dstStdDev));
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
        dataP.ubp = **(UBYTE ***)arrayP;
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
	  ("WlzArrayStats2D FE 0x%lx {%d %d} "
	   "0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )arrayP, arraySize.vtX, arraySize.vtY,
	   (unsigned long )dstMin, (unsigned long )dstMax,
	   (unsigned long )dstSum, (unsigned long )dstSumSq,
	   (unsigned long )dstMean, (unsigned long )dstStdDev));
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
	  ("WlzArrayStats1D FE 0x%lx %d "
	   "0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )arrayP, arraySize,
	   (unsigned long )dstMin, (unsigned long )dstMax,
	   (unsigned long )dstSum, (unsigned long )dstSumSq,
	   (unsigned long )dstMean, (unsigned long )dstStdDev));
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
  UBYTE *bitData,
  WlzIVertex2 arrayOrigin,
  WlzErrorNum *dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzErrorNum 	errNum=WLZ_ERR_NONE;
  UBYTE		**arrayData;
  int		i, j, srcOffset;
  UBYTE		bitmask;

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
	bitmask = (1 << (srcOffset%8));
	if( bitData[srcOffset/8] & bitmask ){
	  arrayData[i][j/8] |= bitmask;
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

#ifdef WLZ_ARRAY_TEST
main(int argc, char *argv[])
{
  WlzObject	*inObj= NULL,
  		*outObj = NULL;
  UBYTE		**array = NULL;
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
#endif /* WLZ_ARRAY_TEST */
