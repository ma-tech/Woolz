#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzArray.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for conversion between woolz domain objects
*		and arrays.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzTo[ISUFD]Array2D					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Extracts an int, short, UBYTE, float or double Alc	*
*		array from any Woolz 2D domain object.			*
* Global refs:	-							*
* Parameters:	WlzIVertex2 *dstSizeArrayDat: Source and destination	*
*					pointer for array size.		*
*		<TYPE> ***dstArrayDat:	Destination pointer for array	*
*					of type: int, short, UBYTE,	*
*					float or long.			*
*		WlzObject *srcObj:	Given woolz object.		*
*		WlzIVertex2 origin:	Array origin wrt given object.	*
*		int noiseFlag:		Fill background with random 	*
*					noise with the same mean and	*
*					std. dev. as the given object	*
*					if non-zero.			*
************************************************************************/
WlzErrorNum WlzToIArray2D(WlzIVertex2 *dstSizeArrayDat, int ***dstArrayDat,
			  WlzObject *srcObj, WlzIVertex2 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  *dstSizeArrayDat, origin,
			  noiseFlag, WLZ_GREY_INT);
  }
  return(errNum);
}

WlzErrorNum WlzToSArray2D(WlzIVertex2 *dstSizeArrayDat, short ***dstArrayDat,
			  WlzObject *srcObj, WlzIVertex2 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  *dstSizeArrayDat, origin,
    			  noiseFlag, WLZ_GREY_SHORT);
  }
  return(errNum);
}

WlzErrorNum WlzToUArray2D(WlzIVertex2 *dstSizeArrayDat, UBYTE ***dstArrayDat,
			  WlzObject *srcObj, WlzIVertex2 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  *dstSizeArrayDat, origin,
    			  noiseFlag, WLZ_GREY_UBYTE);
  }
  return(errNum);
}

WlzErrorNum WlzToFArray2D(WlzIVertex2 *dstSizeArrayDat, float ***dstArrayDat,
			  WlzObject *srcObj, WlzIVertex2 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  *dstSizeArrayDat, origin,
    			  noiseFlag, WLZ_GREY_FLOAT);
  }
  return(errNum);
}

WlzErrorNum WlzToDArray2D(WlzIVertex2 *dstSizeArrayDat, double ***dstArrayDat,
			  WlzObject *srcObj, WlzIVertex2 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray2D((void ***)dstArrayDat, srcObj,
    			  *dstSizeArrayDat, origin,
    			  noiseFlag, WLZ_GREY_DOUBLE);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzToArray2D						*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Extracts an Alc array from any Woolz 2D domain object.	*
*		If the destination pointer points to a non-NULL 	*
*		pointer then it is assumed to be a suitable Alc array.	*
*		The data are assumed to be within the valid range.	*
* Global refs:	-							*
* Parameters:	void ***dstP:		Destination pointer (assumed 	*
*					valid if *dstP is non-NULL).	*
*		WlzObject *srcObj:	Given woolz object.		*
*		WlzIVertex2 size:	Size of the array.		*
*		WlzIVertex2 origin:	Array origin wrt given object.	*
*		int noiseFlag:		Fill background with random 	*
*					noise with the same mean and	*
*					std. dev. as the given object	*
*					if non-zero.			*
*		WlzGreyType dstGreyType: Destination array data type.	*
************************************************************************/
WlzErrorNum	WlzToArray2D(void ***dstP, WlzObject *srcObj,
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

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzToArray2D FE 0x%lx  0x%lx {%d %d} {%d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, origin.vtX, origin.vtY,
	   noiseFlag, (int )dstGreyType));
  gValP.inp = NULL;
  if(dstP == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
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
	    cutObj->values.r->freeptr = WlzPopFreePtr(
	    				     cutObj->values.r->freeptr,
					     &tVP0, NULL);
	    gValP.inp = (int *)tVP0;
	    cutObj->values.r->values.inp = NULL;
	  }
	}
	WlzFreeObj(cutObj);
      }
    }
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzToArray2D FX %d\n",
	   errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzTo[ISUFD]Array3D					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Extracts an int, short, UBYTE, float or double Alc	*
*		array from any Woolz 3D domain object.			*
* Global refs:	-							*
* Parameters:	WlzIVertex3 *dstSizeArrayDat: Source and destination	*
*					pointer for array size.		*
*		<TYPE> ****dstArrayDat:	Destination pointer for array	*
*					of type: int, short, UBYTE,	*
*					float or long.			*
*		WlzObject *srcObj:	Given woolz object.		*
*		WlzIVertex3 origin:	Array origin wrt given object.	*
*		int noiseFlag:		Fill background with random 	*
*					noise with the same mean and	*
*					std. dev. as the given object	*
*					if non-zero.			*
************************************************************************/
WlzErrorNum WlzToIArray3D(WlzIVertex3 *dstSizeArrayDat, int ****dstArrayDat,
			  WlzObject *srcObj, WlzIVertex3 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, *dstSizeArrayDat,
    			  origin, noiseFlag, WLZ_GREY_INT);
  }
  return(errNum);
}

WlzErrorNum WlzToSArray3D(WlzIVertex3 *dstSizeArrayDat, short ****dstArrayDat,
			  WlzObject *srcObj, WlzIVertex3 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, *dstSizeArrayDat,
    			  origin, noiseFlag, WLZ_GREY_SHORT);
  }
  return(errNum);
}

WlzErrorNum WlzToUArray3D(WlzIVertex3 *dstSizeArrayDat, UBYTE ****dstArrayDat,
			  WlzObject *srcObj, WlzIVertex3 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, *dstSizeArrayDat,
    			  origin, noiseFlag, WLZ_GREY_UBYTE);
  }
  return(errNum);
}

WlzErrorNum WlzToFArray3D(WlzIVertex3 *dstSizeArrayDat, float ****dstArrayDat,
			  WlzObject *srcObj, WlzIVertex3 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, *dstSizeArrayDat,
    			  origin, noiseFlag, WLZ_GREY_FLOAT);
  }
  return(errNum);
}

WlzErrorNum WlzToDArray3D(WlzIVertex3 *dstSizeArrayDat, double ****dstArrayDat,
			  WlzObject *srcObj, WlzIVertex3 origin,
			  int noiseFlag)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstSizeArrayDat == NULL) || (dstArrayDat == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzToArray3D((void ****)dstArrayDat, srcObj, *dstSizeArrayDat,
    			  origin, noiseFlag, WLZ_GREY_DOUBLE);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzToArray3D						*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Extracts an Alc array from any Woolz 3D domain object.	*
*		If the destination pointer points to a non-NULL 	*
*		pointer then it is assumed to be a suitable Alc array.	*
*		The data are assumed to be within the valid range.	*
* Global refs:	-							*
* Parameters:	void ****dstP:		Destination pointer (assumed 	*
*					valid if *dstP is non-NULL).	*
*		WlzObject *srcObj:	Given woolz object.		*
*		WlzIVertex3 size:	Size of the array.		*
*		WlzIVertex3 origin:	Array origin wrt given object.	*
*		int noiseFlag:		Fill background with random 	*
*					noise with the same mean and	*
*					std. dev. as the given object	*
*					if non-zero.			*
*		WlzGreyType dstGreyType: Destination array data type.	*
************************************************************************/
WlzErrorNum	WlzToArray3D(void ****dstP, WlzObject *srcObj,
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

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzToArray3D FE 0x%lx  0x%lx {%d %d %d} {%d %d %d} %d %d\n",
	   (unsigned long )dstP, (unsigned long )srcObj,
	   size.vtX, size.vtY, size.vtZ, origin.vtX, origin.vtY, origin.vtZ,
	   noiseFlag, (int )dstGreyType));
  gValP.inp = NULL;
  if(dstP == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
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
	    cutObj->values.vox->freeptr = WlzPopFreePtr(
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzToArray3D FX %d\n",
	   errNum));
  return(errNum);
}

/************************************************************************
* Function:	WlzArrayTxRectValues					*
* Returns:	void							*
* Purpose:	Transforms and/or clamps a rectangle of data values	*
*		using a given buffer.					*
* Global refs:	-							*
* Parameters:	WlzGreyP dstValP:	Destination grey pointer.	*
*		WlzGreyP srcValP:	Source grey pointer.		*
*		double *bufP:		Buffer with space for at least	*
*					row of double grey values.	*
*		WlzIVertex2 rectSize:	The size of the destination and	*
*					source, also the row size for	*
*					the buffer.			*
*		int dstOffset:		Offset from destination ptr.	*
*		int srcOffset:		Offset from source ptr.		*
*		WlzGreyType dstGreyType: Destination grey type.		*
*		WlzGreyType srcGreyType: Source grey type.		*
*		double valOffset:	Offset added to each value.	*
*		double valScale:	Scale factor by which each	*
*					value is multiplied before	*
*					adding the offset.		*
*		int clampFlag:		Values are clamped to the 	*
*					destination type range if the	*
*					clamp flag is non-zero.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzFrom[ISUFD]Array2D					*
* Returns:	WlzObject *:		New woolz object.		*
* Purpose:	Creates a woolz 2D domain object from the given Alc	*
*		array.							*
* Global refs:	-							*
* Parameters:	WlzIVertex2 arraySizeDat: Dimensions of the array.	*
*		<TYPE> **arrayDat:	Given Alc array of type: int,	*
*					short, UBYTE, float or double.	*
*		WlzIVertex2 arrayOrigin: Array origin wrt given object.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error 	*
*					number, may be NULL.		*
************************************************************************/
WlzObject	*WlzFromIArray2D(WlzIVertex2 arraySizeDat,
				 int **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
			WLZ_GREY_INT, WLZ_GREY_INT, 0.0, 1.0,
			0, 0, dstErrNum));
}

WlzObject	*WlzFromSArray2D(WlzIVertex2 arraySizeDat,
				 short **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_SHORT, WLZ_GREY_SHORT, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromUArray2D(WlzIVertex2 arraySizeDat,
				 UBYTE **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_UBYTE, WLZ_GREY_UBYTE, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromFArray2D(WlzIVertex2 arraySizeDat,
				 float **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_FLOAT, WLZ_GREY_FLOAT, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromDArray2D(WlzIVertex2 arraySizeDat,
				 double **arrayDat,
				 WlzIVertex2 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray2D((void **)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
	                0, 0, dstErrNum));
}

/************************************************************************
* Function:	WlzFromArray2D						*
* Returns:	WlzObject *:		New woolz object.		*
* Purpose:	Creates a woolz 2D domain object from the given Alc	*
*		array.							*
*		The data are assumed to be within the valid range.	*
*		If the noCopyFlag is set (non-zero) then the array data	*
*		space is used for the onjects values without copying.	*
*		For this to be valid both the source and destination	*
*		grey type must be the same.				*
* Global refs:	-							*
* Parameters:	void **arrayP:		Given Alc array.		*
*		WlzIVertex2 arraySize:	Dimensions of the array.	*
*		WlzIVertex2 arrayOrigin:	Array origin wrt given object.	*
*		WlzGreyType dstGreyType: Destination object grey type.	*
*		WlzGreyType srcGreyType: Array data type.		*
*		double valOffset:	Offset added to each value.	*
*		double valScale:	Scale factor by which each	*
*					value is multiplied before	*
*					adding the offset.		*
*		int clampFlag:		Values are clamped to the 	*
*					destination type range if the	*
*					clamp flag is non-zero.		*
*		int noCopyFlag:		Use the array data for the	*
*					woolz object values in-place.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error 	*
*					number, may be NULL.		*
************************************************************************/
WlzObject	*WlzFromArray2D(void **arrayP,
				WlzIVertex2 arraySize, WlzIVertex2 arrayOrigin,
				WlzGreyType dstGreyType,
				WlzGreyType srcGreyType,
				double valOffset, double valScale,
				int clampFlag, int noCopyFlag,
				WlzErrorNum *dstErrNum)
{
  int		txFlag = 0;
  unsigned long tUL0;
  double	*bufP = NULL;
  WlzGreyP	dstValP,
  		srcValP;
  WlzPixelV	dstBkgPix;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject	*dstObj = NULL;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFromArray2D FE 0x%lx  {%d %d} {%d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arrayOrigin.vtX, arrayOrigin.vtY,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErrNum));
  dstValP.inp = NULL;
  dstBkgPix.type = dstGreyType;
  (void )memset(&(dstBkgPix.v), 0, sizeof(WlzGreyV));
  if((arrayP == NULL) || (*arrayP == NULL) ||
     (arraySize.vtX <= 0) || (arraySize.vtY <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
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
    dstObj->values.r->freeptr = WlzPushFreePtr(dstObj->values.r->freeptr,
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzFromArray2D FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/************************************************************************
* Function:	WlzFrom[ISUFD]Array3D					*
* Returns:	WlzObject *:		New woolz object.		*
* Purpose:	Creates a woolz 3D domain object from the given Alc	*
*		array.							*
* Global refs:	-							*
* Parameters:	WlzIVertex3 arraySizeDat: Dimensions of the array.	*
*		<TYPE> **arrayDat:	Given Alc array of type: int,	*
*					short, UBYTE, float or double.	*
*		WlzIVertex3 arrayOrigin: Array origin wrt given object.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error 	*
*					number, may be NULL.		*
************************************************************************/
WlzObject	*WlzFromIArray3D(WlzIVertex3 arraySizeDat,
				 int ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_INT, WLZ_GREY_INT, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromSArray3D(WlzIVertex3 arraySizeDat,
				 short ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_SHORT, WLZ_GREY_SHORT, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromUArray3D(WlzIVertex3 arraySizeDat,
				 UBYTE ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_UBYTE, WLZ_GREY_UBYTE, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromFArray3D(WlzIVertex3 arraySizeDat,
				 float ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_FLOAT, WLZ_GREY_FLOAT, 0.0, 1.0,
	                0, 0, dstErrNum));
}

WlzObject	*WlzFromDArray3D(WlzIVertex3 arraySizeDat,
				 double ***arrayDat,
				 WlzIVertex3 arrayOrigin,
				 WlzErrorNum *dstErrNum)
{
  return(WlzFromArray3D((void ***)arrayDat, arraySizeDat, arrayOrigin,
  	                WLZ_GREY_DOUBLE, WLZ_GREY_DOUBLE, 0.0, 1.0,
	                0, 0, dstErrNum));
}

/************************************************************************
* Function:	WlzFromArray3D						*
* Returns:	WlzObject *:		New woolz object.		*
* Purpose:	Creates a woolz 3D domain object from the given Alc	*
*		array.							*
*		The data are assumed to be within the valid range.	*
*		If the noCopyFlag is set (non-zero) then the array data	*
*		space is used for the onjects values without copying.	*
*		For this to be valid both the source and destination	*
*		grey type must be the same.				*
* Global refs:	-							*
* Parameters:	void ***arrayP:		Given Alc array.		*
*		WlzIVertex3 arraySize:	Dimensions of the array.	*
*		WlzIVertex3 arrayOrigin: Array origin wrt given object.	*
*		WlzGreyType dstGreyType: Destination object grey type.	*
*		WlzGreyType srcGreyType: Array data type.		*
*		double valOffset:	Offset added to each value.	*
*		double valScale:	Scale factor by which each	*
*					value is multiplied before	*
*					adding the offset.		*
*		int clampFlag:		Values are clamped to the 	*
*					destination type range if the	*
*					clamp flag is non-zero.		*
*		int noCopyFlag:		Use the array data for the	*
*					woolz object values in-place.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error 	*
*					number, may be NULL.		*
************************************************************************/
WlzObject	*WlzFromArray3D(void ***arrayP,
				WlzIVertex3 arraySize, WlzIVertex3 arrayOrigin,
				WlzGreyType dstGreyType,
				WlzGreyType srcGreyType,
				double valOffset, double valScale,
				int clampFlag, int noCopyFlag,
				WlzErrorNum *dstErrNum)
{
  int		tI0,
  		planeCount,
  		planeIdx,
		planeOffset,
		planePos,
		txFlag = 0;
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

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzFromArray3D FE 0x%lx  {%d %d %d} {%d %d %d} "
	   "%d %d %g %g %d %d 0x%lx\n",
	   (unsigned long )arrayP,
	   arraySize.vtX, arraySize.vtY, arraySize.vtZ,
	   arrayOrigin.vtX, arrayOrigin.vtY, arrayOrigin.vtZ,
	   (int )dstGreyType, (int )srcGreyType,
	   valOffset, valScale,
	   clampFlag, noCopyFlag, (unsigned long )dstErrNum));
  dstDom.core = NULL;
  dstValues.core = NULL;
  dstValP.inp = NULL;
  dstBkgPix.type = dstGreyType;
  (void )memset(&(dstBkgPix.v), 0, sizeof(WlzGreyV));
  arraySize2D.vtX = arraySize.vtX;
  arraySize2D.vtY = arraySize.vtY;
  if((arrayP == NULL) || (*arrayP == NULL) || (**arrayP == NULL) ||
     (arraySize.vtX <= 0) || (arraySize.vtY <= 0) || (arraySize.vtZ <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
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
      tI0 = arraySize.vtX * arraySize.vtY * arraySize.vtZ;
      switch(dstGreyType)
      {
	case WLZ_GREY_INT:
	  if((dstValP.inp = (int *)AlcMalloc((unsigned long )(tI0 *
					                sizeof(int)))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  if((dstValP.shp = (short *)AlcMalloc((unsigned long )(tI0 *
					              sizeof(short)))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  if((dstValP.ubp = (UBYTE *)AlcMalloc((unsigned long )(tI0 *
					              sizeof(UBYTE)))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  if((dstValP.flp = (float *)AlcMalloc((unsigned long )(tI0 *
					              sizeof(float)))) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  if((dstValP.dbp = (double *)AlcMalloc((unsigned long )(tI0 *
						     sizeof(double)))) == NULL)
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
      dstValues.vox->freeptr = WlzPushFreePtr(dstValues.vox->freeptr,
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
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
          ("WlzFromArray2D FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}

/************************************************************************
* Function:	WlzArrayStats3D						*
* Returns:	int:			Number of data in array.	*
* Purpose:	Calculates simple statistics for the given Alc array.	*
* Global refs:	-							*
* Parameters:	void ***arrayP:		Given 3D Alc array.		*
*		WlzIVertex3 arraySize:	Dimensions of the array.	*
*		WlzGreyType greyType:	Array data type.		*
*		double *dstMin:		Destination ptr for minimum	*
*					value, may be NULL.		*
*		double *dstMax:		Destination ptr for maximum	*
*					value, may be NULL.		*
*		double *dstSum:		Destination ptr for sum of	*
*					values, may be NULL.		*
*		double *dstSumSq:	Destination ptr for sum of	*
*					squares of values, may be NULL.	*
*		double *dstMean:	Destination ptr for mean of	*
*					values, may be NULL.		*
*		double *dstStdDev:	Destination ptr for std. dev. 	*
*					of values, may be NULL.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzArrayStats2D						*
* Returns:	int:			Number of data in array.	*
* Purpose:	Calculates simple statistics for the given Alc array.	*
* Global refs:	-							*
* Parameters:	void **arrayP:		Given 2D Alc array.		*
*		WlzIVertex2 arraySize:	Dimensions of the array.	*
*		WlzGreyType greyType:	Array data type.		*
*		double *dstMin:		Destination ptr for minimum	*
*					value, may be NULL.		*
*		double *dstMax:		Destination ptr for maximum	*
*					value, may be NULL.		*
*		double *dstSum:		Destination ptr for sum of	*
*					values, may be NULL.		*
*		double *dstSumSq:	Destination ptr for sum of	*
*					squares of values, may be NULL.	*
*		double *dstMean:	Destination ptr for mean of	*
*					values, may be NULL.		*
*		double *dstStdDev:	Destination ptr for std. dev. 	*
*					of values, may be NULL.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzArrayStats1D						*
* Returns:	int:			Number of data in array.	*
* Purpose:	Calculates simple statistics for the given Alc array.	*
* Global refs:	-							*
* Parameters:	void *arrayP:		Given 1D Alc array.		*
*		int arraySize:		Dimension of the array.		*
*		WlzGreyType greyType:	Array data type.		*
*		double *dstMin:		Destination ptr for minimum	*
*					value, may be NULL.		*
*		double *dstMax:		Destination ptr for maximum	*
*					value, may be NULL.		*
*		double *dstSum:		Destination ptr for sum of	*
*					values, may be NULL.		*
*		double *dstSumSq:	Destination ptr for sum of	*
*					squares of values, may be NULL.	*
*		double *dstMean:	Destination ptr for mean of	*
*					values, may be NULL.		*
*		double *dstStdDev:	Destination ptr for std. dev. 	*
*					of values, may be NULL.		*
************************************************************************/
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

