#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzValueUtils.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Implements many functions for setting and copying
*		vectors of simple values.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzValueSet<VecType>					*
* Returns:	void							*
* Purpose:	Sets the elements of the given vector to a given 	*
*		value, where the vector type is any one of int,		*
*		short, UBYTE, float, double, WlzIVertex2, WlzFVertex2 or	*
*		WlzDVertex2.						*
* Global refs:	-							*
* Parameters:	<VecType> *vec:		Vector who's elements are to be	*
*					set.				*
*		<VecType> value:	Value to use when setting the	*
*					vector's elements.		*
*		int count:		Number of vector elements to	*
*					be set.				*
************************************************************************/
void		WlzValueSetInt(int *vec, int value,
			       int count)
{
  if(value)
  {
    while(count-- > 0)
    {
      *vec++ = value;
    }
  }
  else
  {
    (void )memset(vec, 0, count * sizeof(int));
  }
}

void		WlzValueSetShort(short *vec, short value,
				 int count)
{
  if(value)
  {
    while(count-- > 0)
    {
      *vec++ = value;
    }
  }
  else
  {
    (void )memset(vec, 0, count * sizeof(short));
  }
}

void		WlzValueSetUByte(UBYTE *vec, UBYTE value,
				 int count)
{
  (void )memset(vec, value, count);
}

void		WlzValueSetFloat(float *vec, float value,
				 int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

void		WlzValueSetDouble(double *vec, double value,
				  int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

void		WlzValueSetDVertex(WlzDVertex2 *vec, WlzDVertex2 value,
				   int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

void		WlzValueSetFVertex(WlzFVertex2 *vec, WlzFVertex2 value,
				   int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

void		WlzValueSetIVertex(WlzIVertex2 *vec, WlzIVertex2 value,
				   int count)
{
  if(value.vtX || value.vtY)
  {
    while(count-- > 0)
    {
      *vec++ = value;
    }
  }
  else
  {
    (void )memset(vec, 0, count * sizeof(WlzIVertex2));
  }
}

/************************************************************************
* Function:	WlzValueSetGrey						*
* Returns:	void							*
* Purpose:	Sets the elements of the given vector to a given 	*
*		value, where the vector type is any one of int,		*
*		short, UBYTE, float, double.				*
* Global refs:	-							*
* Parameters:	WlzGreyP *vec:		Vector who's elements are to be	*
*					set.				*
*		int vecOff:		Offset from vec.		*
*		WlzGreyV value:		Value to use when setting the	*
*					vector's elements.		*
*		WlzGreyType gType:	Grey type, ie: int, short....	*
*		int count		Number of vector elements to	*
*					be set.				*
************************************************************************/
void		WlzValueSetGrey(WlzGreyP vec, int vecOff, WlzGreyV value,
				WlzGreyType gType, int count)
{
  switch(gType)
  {
    case WLZ_GREY_INT:
      WlzValueSetInt(vec.inp + vecOff, value.inv, count);
      break;
    case WLZ_GREY_SHORT:
      WlzValueSetShort(vec.shp + vecOff, value.shv, count);
      break;
    case WLZ_GREY_UBYTE:
      WlzValueSetUByte(vec.ubp + vecOff, value.ubv, count);
      break;
    case WLZ_GREY_FLOAT:
      WlzValueSetFloat(vec.flp + vecOff, value.flv, count);
      break;
    case WLZ_GREY_DOUBLE:
      WlzValueSetDouble(vec.dbp + vecOff, value.dbv, count);
      break;
  }
}

/************************************************************************
* Function:	WlzValueClamp<VecType>To<ReqType>			*
* Returns:	void							*
* Purpose:	Clamps a vector of <VecType> to the limits of		*
*		<ReqType>, where the source and destination types are	*
*		any combination of int, short, UBYTE, float or double	*
*		in which the range of <VecType> is not contained within	*
*		the range of <ReqType>.					*
* Global refs:	-							*
* Parameters:	<VecType> *vec:		Vector who's elements are to be	*
*					clamped.			*
*		int count:		Number of vector elements.	*
************************************************************************/
void		 WlzValueClampIntToShort(int *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > SHRT_MAX)
    {
      *vec = SHRT_MAX;
    }
    else if(*vec < SHRT_MIN)
    {
      *vec = SHRT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampIntToUByte(int *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > UCHAR_MAX)
    {
      *vec = UCHAR_MAX;
    }
    else if(*vec < 0)
    {
      *vec = 0;
    }
    ++vec;
  }
}

void		 WlzValueClampShortToUByte(short *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > UCHAR_MAX)
    {
      *vec = UCHAR_MAX;
    }
    else if(*vec < 0)
    {
      *vec = 0;
    }
    ++vec;
  }
}

void		 WlzValueClampDoubleToInt(double *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > INT_MAX)
    {
      *vec = INT_MAX;
    }
    else if(*vec < INT_MIN)
    {
      *vec = INT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampDoubleToShort(double *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > SHRT_MAX)
    {
      *vec = SHRT_MAX;
    }
    else if(*vec < SHRT_MIN)
    {
      *vec = SHRT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampDoubleToUByte(double *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > UCHAR_MAX)
    {
      *vec = UCHAR_MAX;
    }
    else if(*vec < 0)
    {
      *vec = 0;
    }
    ++vec;
  }
}

void		 WlzValueClampDoubleToFloat(double *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > FLT_MAX)
    {
      *vec = FLT_MAX;
    }
    else if(*vec < FLT_MIN)
    {
      *vec = FLT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampFloatToInt(float *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > INT_MAX)
    {
      *vec = INT_MAX;
    }
    else if(*vec < INT_MIN)
    {
      *vec = INT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampFloatToShort(float *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > SHRT_MAX)
    {
      *vec = SHRT_MAX;
    }
    else if(*vec < SHRT_MIN)
    {
      *vec = SHRT_MIN;
    }
    ++vec;
  }
}

void		 WlzValueClampFloatToUByte(float *vec, int count)
{
  while(count-- > 0)
  {
    if(*vec > UCHAR_MAX)
    {
      *vec = UCHAR_MAX;
    }
    else if(*vec < 0)
    {
      *vec = 0;
    }
    ++vec;
  }
}

/************************************************************************
* Function:	WlzValueMask<VecType>To<ReqType>			*
* Returns:	void							*
* Purpose:	Masks a vector of <VecType> to the limits of		*
*		<ReqType>, where the source and destination types are	*
*		any combination of int, short or UBYTE in which		*
*		the range of <VecType> is not contained within the	*
*		range of <ReqType>.					*
* Global refs:	-							*
* Parameters:	<VecType> *vec:		Vector who's elements are to be	*
*					clamped.			*
*		int count:		Number of vector elements.	*
************************************************************************/
void		 WlzValueMaskIntToShort(int *vec, int count)
{
  while(count-- > 0)
  {
    *vec++ = *vec & USHRT_MAX;
  }
}

void		 WlzValueMaskIntToUByte(int *vec, int count)
{
  while(count-- > 0)
  {
    *vec++ = *vec & UCHAR_MAX;
  }
}

void		 WlzValueMaskShortToUByte(short *vec, int count)
{
  while(count-- > 0)
  {
    *vec++ = *vec & UCHAR_MAX;
  }
}

/************************************************************************
* Function:	WlzValueCopy<SrcType>To<DstType>			*
* Returns:	void							*
* Purpose:	Copies a vector of <SrcType> to a vector of <DstType>,	*
*		where the source and destination types are any		*
*		combination of int, short, UBYTE, float or double.	*
* Global refs:	-							*
* Parameters:	<SrcType> *dst:		Destination vector.		*
*		<DstType> *src:		Source vector.			*
*		int count:		Number of vector elements to	*
*					be copied.			*
************************************************************************/
void		WlzValueCopyIntToInt(int *dst, int *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(int));
}

void		WlzValueCopyIntToShort(short *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyIntToUByte(UBYTE *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (unsigned )*src++;
  }
}

void		WlzValueCopyIntToFloat(float *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyIntToDouble(double *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyShortToInt(int *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyShortToShort(short *dst, short *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(short));
}

void		WlzValueCopyShortToUByte(UBYTE *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (unsigned )*src++;
  }
}

void		WlzValueCopyShortToFloat(float *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyShortToDouble(double *dst, short *src,
					        int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyUByteToInt(int *dst, UBYTE *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyUByteToShort(short *dst, UBYTE *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyUByteToUByte(UBYTE *dst, UBYTE *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(UBYTE));
}

void		WlzValueCopyUByteToFloat(float *dst, UBYTE *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyUByteToDouble(double *dst, UBYTE *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyFloatToInt(int *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyFloatToShort(short *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyFloatToUByte(UBYTE *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (unsigned )WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyFloatToFloat(float *dst, float *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(float));
}

void		WlzValueCopyFloatToDouble(double *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyDoubleToInt(int *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyDoubleToShort(short *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyDoubleToUByte(UBYTE *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (unsigned )WLZ_NINT(*src);
    ++src;
  }
}

void		WlzValueCopyDoubleToFloat(float *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

void		WlzValueCopyDoubleToDouble(double *dst, double *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(double));
}

/************************************************************************
* Function:	WlzValueCopyGreyToGrey					*
* Returns:	void							*
* Purpose:	Copies a source vector to a destination vector,		*
*		where the source and destination types are any		*
*		combination of int, short, UBYTE, float or double.	*
* Global refs:	-							*
* Parameters:	WlzGreyP dst:		Destination vector.		*
*		int dstOff:		Destination offset.		*
*		WlzGreyType dstType:	Destination type, ie: int, ....	*
*		WlzGreyP src:		Source vector.			*
*		int srcOff:		Source offset.			*
*		WlzGreyType srcType:	Source type, ie: int, ....	*
*		int count:		Number of vector elements to	*
*					be copied.			*
************************************************************************/
void		WlzValueCopyGreyToGrey(WlzGreyP dst, int dstOff,
				       WlzGreyType dstType,
				       WlzGreyP src, int srcOff,
				       WlzGreyType srcType,
				       int count)
{
  switch(dstType)
  {
    case WLZ_GREY_INT:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToInt(dst.inp + dstOff, src.inp + srcOff,
	  		       count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToInt(dst.inp + dstOff, src.shp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToInt(dst.inp + dstOff, src.ubp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToInt(dst.inp + dstOff, src.flp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToInt(dst.inp + dstOff, src.dbp + srcOff,
	  			  count);
	  break;
      }
      break;
    case WLZ_GREY_SHORT:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToShort(dst.shp + dstOff, src.inp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToShort(dst.shp + dstOff, src.shp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToShort(dst.shp + dstOff, src.ubp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToShort(dst.shp + dstOff, src.flp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToShort(dst.shp + dstOff, src.dbp + srcOff,
	  			    count);
	  break;
      }
      break;
    case WLZ_GREY_UBYTE:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToUByte(dst.ubp + dstOff, src.inp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToUByte(dst.ubp + dstOff, src.shp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToUByte(dst.ubp + dstOff, src.ubp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToUByte(dst.ubp + dstOff, src.flp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToUByte(dst.ubp + dstOff, src.dbp + srcOff,
	  			    count);
	  break;
      }
      break;
    case WLZ_GREY_FLOAT:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToFloat(dst.flp + dstOff, src.inp + srcOff,
	  			 count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToFloat(dst.flp + dstOff, src.shp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToFloat(dst.flp + dstOff, src.ubp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToFloat(dst.flp + dstOff, src.flp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToFloat(dst.flp + dstOff, src.dbp + srcOff,
	  			    count);
	  break;
      }
      break;
    case WLZ_GREY_DOUBLE:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToDouble(dst.dbp + dstOff, src.inp + srcOff,
	  			  count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToDouble(dst.dbp + dstOff, src.shp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToDouble(dst.dbp + dstOff, src.ubp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToDouble(dst.dbp + dstOff, src.flp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToDouble(dst.dbp + dstOff, src.dbp + srcOff,
	  			     count);
	  break;
      }
      break;
  }
}

/************************************************************************
* Function:	WlzValueCopy<SrcType>To<DstType>			*
* Returns:	void							*
* Purpose:	Copies a vector of <SrcType> verticies to a vector of	*
*		<DstType> verticies, where the source and destination	*
*		types are any combination of WlzDVertex2, WlzFVertex2	*
*		and WlzIVertex2.						*
* Global refs:	-							*
* Parameters:	<SrcType> *dst:		Destination vector.		*
*		<DstType> *src:		Source vector.			*
*		int count:		Number of vector elements to	*
*					be copied.			*
************************************************************************/
void		WlzValueCopyDVertexToDVertex(WlzDVertex2 *dst,
					     WlzDVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzDVertex2));
}

void		WlzValueCopyDVertexToFVertex(WlzFVertex2 *dst,
					     WlzDVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = WLZ_NINT(src->vtX);
    dst->vtY = WLZ_NINT(src->vtY);
    ++dst;
    ++src;
  }
}

void		WlzValueCopyDVertexToIVertex(WlzIVertex2 *dst,
					     WlzDVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = WLZ_NINT(src->vtX);
    dst->vtY = WLZ_NINT(src->vtY);
    ++dst;
    ++src;
  }
}

void		WlzValueCopyFVertexToDVertex(WlzDVertex2 *dst,
					     WlzFVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = src->vtX;
    dst->vtY = src->vtY;
    ++dst;
    ++src;
  }
}

void		WlzValueCopyFVertexToFVertex(WlzFVertex2 *dst,
					     WlzFVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzFVertex2));
}

void		WlzValueCopyFVertexToIVertex(WlzIVertex2 *dst,
					     WlzFVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = WLZ_NINT(src->vtX);
    dst->vtY = WLZ_NINT(src->vtY);
    ++dst;
    ++src;
  }
}

void		WlzValueCopyIVertexToDVertex(WlzDVertex2 *dst,
					     WlzIVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = src->vtX;
    dst->vtY = src->vtY;
    ++dst;
    ++src;
  }
}

void		WlzValueCopyIVertexToFVertex(WlzFVertex2 *dst,
					     WlzIVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = src->vtX;
    dst->vtY = src->vtY;
    ++dst;
    ++src;
  }
}

void		WlzValueCopyIVertexToIVertex(WlzIVertex2 *dst,
					     WlzIVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzIVertex2));
}

/************************************************************************
* Function:	WlzValueConvertPixel					*
* Returns:	WlzErrorNum:		WLZ_ERR_NONE,			*
*					WLZ_ERR_PARAM_NULL or		*
*					WLZ_ERR_GREY_TYPE.		*
* Purpose:	Converts a single pixel value.				*
*		Source values are clamped to the the destination value	*
*		range.							*
* Global refs:	-							*
* Parameters:	WlzPixelV *dstPix	Destination pointer for pixel.	*
*		WlzPixelV srcPix:	Source pixel.			*
*		WlzGreyType dstType:	Destination type, ie: int, ....	*
************************************************************************/
WlzErrorNum	WlzValueConvertPixel(WlzPixelV *dstPix,
				     WlzPixelV srcPix,
				     WlzGreyType dstType)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(dstPix == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(srcPix.type)
    {
      case WLZ_GREY_INT:
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_SHORT:
	    dstPix->type = dstType;
	    WlzValueClampIntToShort(&(srcPix.v.inv), 1);
	    dstPix->v.shv = srcPix.v.inv;
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampIntToUByte(&(srcPix.v.inv), 1);
	    dstPix->v.ubv = (unsigned )(srcPix.v.inv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = srcPix.v.inv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.inv;
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_GREY_SHORT:
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    dstPix->type = dstType;
	    dstPix->v.inv = srcPix.v.shv;
	    break;
	  case WLZ_GREY_SHORT:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampShortToUByte(&(srcPix.v.shv), 1);
	    dstPix->v.ubv =  (unsigned )(srcPix.v.shv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = srcPix.v.shv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.shv;
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_GREY_UBYTE:
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    dstPix->type = dstType;
	    dstPix->v.inv = srcPix.v.ubv;
	    break;
	  case WLZ_GREY_SHORT:
	    dstPix->type = dstType;
	    dstPix->v.shv = srcPix.v.ubv;
	    break;
	  case WLZ_GREY_UBYTE:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = srcPix.v.ubv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.ubv;
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_GREY_FLOAT:
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    dstPix->type = dstType;
	    WlzValueClampFloatToInt(&(srcPix.v.flv), 1);
	    dstPix->v.inv = WLZ_NINT(srcPix.v.flv);
	    break;
	  case WLZ_GREY_SHORT:
	    dstPix->type = dstType;
	    WlzValueClampFloatToShort(&(srcPix.v.flv), 1);
	    dstPix->v.shv = WLZ_NINT(srcPix.v.flv);
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampFloatToUByte(&(srcPix.v.flv), 1);
	    dstPix->v.ubv = (unsigned )(WLZ_NINT(srcPix.v.flv));
	    break;
	  case WLZ_GREY_FLOAT:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.flv;
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_GREY_DOUBLE:
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToInt(&(srcPix.v.dbv), 1);
	    dstPix->v.inv = WLZ_NINT(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_SHORT:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToShort(&(srcPix.v.dbv), 1);
	    dstPix->v.shv =  WLZ_NINT(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToUByte(&(srcPix.v.dbv), 1);
	    dstPix->v.ubv = (unsigned )WLZ_NINT(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToFloat(&(srcPix.v.dbv), 1);
	    dstPix->v.flv = srcPix.v.dbv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    *dstPix = srcPix;
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      default:
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzValueMedian<VecType>					*
* Returns:	void							*
* Purpose:	Computes the median of the given vector (the vector's	*
*		values are overwritten), where the vector type is 	*
*		either int or double.					*
*		The algorithm is based on min-max elimination as used	*
*		by Peath A W. Median finding on a 3x3 Grid. Graphics	*
*		Gems, 171-175. Academic Press. 1990.			*
*		If the vector has less than 1 element 0 is returned.	*
*		In the trivial case of 1 or 2 elements, the first	*
*		element's value is returned.				*
*		If the vector has more than 3 elements:			*
*		* The given vector is partitioned into two portions:	*
*		  the 1st having ((n + 1) / 2) + 1 elements.		*
*		* While all the elements have not been examined:	*
*		*   Find maximum and minimum values in the 1st		*
*		    partition.						*
*		*   Replace the min value with the first in the		*
*		    partition and the max with the last in the		*
*		    partition.						*
*		*   Replace the first value in the first partition with	*
*		    the next value from the second partition.		*
*		*   Reduce the length of the first partion by 1, ie	*
*		    omit the last value.				*
*		* This leaves 3 values in the 1st partion, their median	*
*		  is found by direct comparison.			*
* Global refs:	-							*
* Parameters:	<VecType> *vec:		Vector who's elements are to be	*
*					used in computing the median.	*
*		int count:		Number of vector elements.	*
************************************************************************/
int		WlzValueMedianInt(int *values, int nVal)
{
  int		tmpV,
  		idP,
		idV,
  		nPar,
		medVal;
  int		*minP,
  		*maxP,
		*valP0,
		*valP1;

  if(nVal <= 0)
  {
    medVal = 0;
  }
  else if(nVal < 3)
  {
    medVal = *values;
  }
  else
  {
    /* Partition the given values: values[0 - (nPar - 1)]. */
    nPar = ((nVal + 1) >> 1) + 1;
    idV = nPar;
    while(idV < nVal)
    {
      /* Find the minimum and maximum values in the partition. */
      minP = maxP = valP0 = values;
      for(idP = 1; idP < nPar; ++idP)
      {
	++valP0;
	if(*valP0 > *maxP)
	{
	  maxP = valP0;
	}
	else if(*valP0 < *minP)
	{
	  minP = valP0;
	}
      }
      /* Put 1st value in place of the minimum value in partition. */
      *minP = *values;
      /* Put last value in partition in place of maximum value in partition. */
      *maxP = *(values + nPar - 1);
      /* Set up the next partition. */
      --nPar;
      *values = *(values + idV);
      ++idV;
    }
    /* Partition now has 3 values: Find the median. */
    valP0 = values;
    valP1 = values + 1;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    valP0 = values + 1;
    valP1 = values + 2;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    valP0 = values;
    valP1 = values + 1;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    medVal = *(values + 1);
  }
  return(medVal);
}

double		WlzValueMedianDouble(double *values, int nVal)
{
  int		idP,
		idV,
  		nPar;
  double	tmpV,
  		medVal;
  double	*minP,
  		*maxP,
		*valP0,
		*valP1;

  if(nVal <= 0)
  {
    medVal = 0.0;
  }
  else if(nVal < 3)
  {
    medVal = *values;
  }
  else
  {
    /* Partition the given values: values[0 - (nPar - 1)]. */
    nPar = ((nVal + 1) >> 1) + 1;
    idV = nPar;
    while(idV < nVal)
    {
      /* Find the minimum and maximum values in the partition. */
      minP = maxP = valP0 = values;
      for(idP = 1; idP < nPar; ++idP)
      {
	++valP0;
	if(*valP0 > *maxP)
	{
	  maxP = valP0;
	}
	else if(*valP0 < *minP)
	{
	  minP = valP0;
	}
      }
      /* Put 1st value in place of the minimum value in partition. */
      *minP = *values;
      /* Put last value in partition in place of maximum value in partition. */
      *maxP = *(values + nPar - 1);
      /* Set up the next partition. */
      --nPar;
      *values = *(values + idV);
      ++idV;
    }
    /* Partition now has 3 values: Find the median. */
    valP0 = values;
    valP1 = values + 1;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    valP0 = values + 1;
    valP1 = values + 2;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    valP0 = values;
    valP1 = values + 1;
    if(*valP0 > *valP1)
    {
      tmpV = *valP0; *valP0 = *valP1; *valP1 = tmpV;
    }
    medVal = *(values + 1);
  }
  return(medVal);
}
