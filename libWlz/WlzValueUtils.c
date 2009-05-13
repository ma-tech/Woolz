#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzValueUtils_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzValueUtils.c
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
* \brief	Many small functions for setting, copying and converting
* 		values.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is int.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is short.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is WlzUByte.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetUByte(WlzUByte *vec, WlzUByte value,
				 int count)
{
  (void )memset(vec, value, count);
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is float.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetFloat(float *vec, float value,
				 int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is double.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetDouble(double *vec, double value,
				  int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is rgb-alpha.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetRGBA(WlzUInt *vec, WlzUInt value,
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
    (void )memset(vec, 0, count * sizeof(WlzUInt));
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is WlzDVertex2.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetDVertex(WlzDVertex2 *vec, WlzDVertex2 value,
				   int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is WlzFVertex2.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
void		WlzValueSetFVertex(WlzFVertex2 *vec, WlzFVertex2 value,
				   int count)
{
  while(count-- > 0)
  {
    *vec++ = value;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given value,
*		where the vector type is WlzIVertex2.
* \param	vec			Vector who's elements are to be set.
* \param	value			Value to use when setting the vector's
*					elements.
* \param	count			Number of vector elements to be set.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Sets the elements of the given vector to a given
*               value, where the vector type is any one of int,
*               short, WlzUByte, float, double.
* \param	vec			Vector who's elements are to be set.
* \param	vecOff			Offset from vec.
* \param	value			Value to use when setting the
*                                       vector's elements.
* \param	gType			Grey type, ie: int, short....
* \param	count			Number of vector elements to
*					be set.
*/
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
    case WLZ_GREY_RGBA:
      WlzValueSetRGBA(vec.rgbp + vecOff, value.rgbv, count);
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of int values to the limits of short.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of int values to the limits of WlzUByte.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of short values to the limits of WlzUByte.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values to the limits of int.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values to the limits of short.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values to the limits of WlzUByte.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values to the limits of RGBA.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
void		 WlzValueClampDoubleToRGBA(double *vec, int count)
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values to the limits of float.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values to the limits of int.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values to the limits of short.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values to the limits of WlzUByte.
* \param	vec			Vector who's elements are to be
*					clamped.
* \param	count			Number of vector elements.
*/
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

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of int values into a vector of short
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		 WlzValueClampIntIntoShort(short *dst, int *src, int count)
{
  while(count-- > 0)
  {
    if(*src > SHRT_MAX)
    {
      *dst = SHRT_MAX;
    }
    else if(*src < SHRT_MIN)
    {
      *dst = SHRT_MIN;
    }
    else
    {
      *dst = (short )*src;
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of nt values into a vector of WlzUByte
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		 WlzValueClampIntIntoUByte(WlzUByte *dst, int *src, int count)
{
  while(count-- > 0)
  {
    if(*src > UCHAR_MAX)
    {
      *dst = UCHAR_MAX;
    }
    else if(*src < 0)
    {
      *dst = 0;
    }
    else
    {
      *dst = (WlzUByte )*src;
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of short values into a vector of WlzUByte
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		 WlzValueClampShortIntoUByte(WlzUByte *dst, short *src,
					     int count)
{
  while(count-- > 0)
  {
    if(*src > UCHAR_MAX)
    {
      *dst = UCHAR_MAX;
    }
    else if(*src < 0)
    {
      *dst = 0;
    }
    else
    {
      *dst = (WlzUByte )*src;
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values into a vector of int
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampFloatIntoInt(int *dst, float *src, int count)
{
  while(count-- > 0)
  {
    if(*src > INT_MAX)
    {
      *dst = INT_MAX;
    }
    else if(*src < INT_MIN)
    {
      *dst = INT_MIN;
    }
    else
    {
      *dst = WLZ_NINT(*src);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values into a vector of short
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampFloatIntoShort(short *dst, float *src, int count)
{
  while(count-- > 0)
  {
    if(*src > SHRT_MAX)
    {
      *dst = SHRT_MAX;
    }
    else if(*src < SHRT_MIN)
    {
      *dst = SHRT_MIN;
    }
    else
    {
      *dst = (short )WLZ_NINT(*src);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values into a vector of WlzUByte
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampFloatIntoUByte(WlzUByte *dst, float *src,
					    int count)
{
  while(count-- > 0)
  {
    if(*src > UCHAR_MAX)
    {
      *dst = UCHAR_MAX;
    }
    else if(*src < 0.0)
    {
      *dst = 0;
    }
    else
    {
      *dst = (WlzUByte )(*src + 0.5);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of into double values into a vector of int
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampDoubleIntoInt(int *dst, double *src, int count)
{
  while(count-- > 0)
  {
    if(*src > INT_MAX)
    {
      *dst = INT_MAX;
    }
    else if(*src < INT_MIN)
    {
      *dst = INT_MIN;
    }
    else
    {
      *dst = WLZ_NINT(*src);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values into a vector of short
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampDoubleIntoShort(short *dst, double *src, int count)
{
  while(count-- > 0)
  {
    if(*src > SHRT_MAX)
    {
      *dst = SHRT_MAX;
    }
    else if(*src < SHRT_MIN)
    {
      *dst = SHRT_MIN;
    }
    else
    {
      *dst = (short )WLZ_NINT(*src);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values into a vector of WlzUByte
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampDoubleIntoUByte(WlzUByte *dst, double *src,
					     int count)
{
  while(count-- > 0)
  {
    if(*src > UCHAR_MAX)
    {
      *dst = UCHAR_MAX;
    }
    else if(*src < 0.0)
    {
      *dst = 0;
    }
    else
    {
      *dst = (WlzUByte )(*src + 0.5);
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values into a vector of float
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampDoubleIntoFloat(float *dst, double *src, int count)
{
  while(count-- > 0)
  {
    if(*src > FLT_MAX)
    {
      *dst = FLT_MAX;
    }
    else if(*src < FLT_MIN)
    {
      *dst = FLT_MIN;
    }
    else
    {
      *dst = (float )*src;
    }
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of int values into a vector of RGBA
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampIntIntoRGBA(WlzUInt *dst, int *src, int count)
{
  WlzUInt val;
  while(count-- > 0)
  {
    if(*src > 255)
    {
      val = 255;
    }
    else if(*src < 0)
    {
      val = 0;
    }
    else
    {
      val = (WlzUInt) *src;
    }
    WLZ_RGBA_RGBA_SET(*dst, val, val, val, 255);
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of short values into a vector of RGBA
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampShortIntoRGBA(WlzUInt *dst, short *src, int count)
{
  WlzUInt val;
  while(count-- > 0)
  {
    if(*src > 255)
    {
      val = 255;
    }
    else if(*src < 0)
    {
      val = 0;
    }
    else
    {
      val = (WlzUInt) *src;
    }
    WLZ_RGBA_RGBA_SET(*dst, val, val, val, 255);
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of float values into a vector of RGBA
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampFloatIntoRGBA(WlzUInt *dst, float *src, int count)
{
  WlzUInt val;
  while(count-- > 0)
  {
    if(*src > 255)
    {
      val = 255;
    }
    else if(*src < 0)
    {
      val = 0;
    }
    else
    {
      val = (WlzUInt) *src;
    }
    WLZ_RGBA_RGBA_SET(*dst, val, val, val, 255);
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValuesUtils
* \brief	Clamps a vector of double values into a vector of RGBA
*		values.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements.
*/
void		WlzValueClampDoubleIntoRGBA(WlzUInt *dst, double *src,
					    int count)
{
  WlzUInt val;
  while(count-- > 0)
  {
    if(*src > 255)
    {
      val = 255;
    }
    else if(*src < 0)
    {
      val = 0;
    }
    else
    {
      val = (WlzUInt) *src;
    }
    WLZ_RGBA_RGBA_SET(*dst, val, val, val, 255);
    ++src;
    ++dst;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Clamps a source vector into a destination vector,
*               where the source and destination types are any
*               combination of int, short, WlzUByte, float or double.
* \param	dst			Destination vector.
* \param	dstOff			Destination offset.
* \param	dstType			Destination type, eg: WLZ_GREY_INT.
* \param	src			Source vector.
* \param	srcOff			Source offset.
* \param	srcType			Source type, eg: WLZ_GREY_SHORT.
* \param	count			Number of vector elements to be
*					clamped and copied.
*/
void		WlzValueClampGreyIntoGrey(WlzGreyP dst, int dstOff,
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
	  WlzValueClampFloatIntoInt(dst.inp + dstOff, src.flp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueClampDoubleIntoInt(dst.inp + dstOff, src.dbp + srcOff,
	  			     count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToInt(dst.inp + dstOff, src.rgbp + srcOff,
	  			     count);
	  break;
	default:
	  break;
      }
      break;
    case WLZ_GREY_SHORT:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueClampIntIntoShort(dst.shp + dstOff, src.inp + srcOff,
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
	  WlzValueClampFloatIntoShort(dst.shp + dstOff, src.flp + srcOff,
	  			      count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueClampDoubleIntoShort(dst.shp + dstOff, src.dbp + srcOff,
	  			       count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToShort(dst.shp + dstOff, src.rgbp + srcOff,
	  			       count);
	  break;
	default:
	  break;
      }
      break;
    case WLZ_GREY_UBYTE:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueClampIntIntoUByte(dst.ubp + dstOff, src.inp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueClampShortIntoUByte(dst.ubp + dstOff, src.shp + srcOff,
	  			      count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToUByte(dst.ubp + dstOff, src.ubp + srcOff,
	  			   count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueClampFloatIntoUByte(dst.ubp + dstOff, src.flp + srcOff,
	  			      count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueClampDoubleIntoUByte(dst.ubp + dstOff, src.dbp + srcOff,
	  			       count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToUByte(dst.ubp + dstOff, src.rgbp + srcOff,
	  			       count);
	  break;
	default:
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
	  WlzValueClampDoubleIntoFloat(dst.flp + dstOff, src.dbp + srcOff,
	  			       count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToFloat(dst.flp + dstOff, src.rgbp + srcOff,
	  			       count);
	  break;
	default:
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToDouble(dst.dbp + dstOff, src.rgbp + srcOff,
	  			     count);
	  break;
	default:
	  break;
      }
      break;
    case WLZ_GREY_RGBA:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueClampIntIntoRGBA(dst.rgbp + dstOff, src.inp + srcOff,
	  			  count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueClampShortIntoRGBA(dst.rgbp + dstOff, src.shp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToRGBA(dst.rgbp + dstOff, src.ubp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueClampFloatIntoRGBA(dst.rgbp + dstOff, src.flp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueClampDoubleIntoRGBA(dst.rgbp + dstOff, src.dbp + srcOff,
	  			     count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToRGBA(dst.rgbp + dstOff, src.rgbp + srcOff,
	  			     count);
	  break;
	default:
	  break;
      }
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of int.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToInt(int *dst, int *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(int));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of short.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToShort(short *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (short )*src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of WlzUByte.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToUByte(WlzUByte *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (WlzUByte )*src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of float.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToFloat(float *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (float )*src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of double.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToDouble(double *dst, int *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of int to a vector of RGBA.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIntToRGBA(WlzUInt *dst, int *src, int count)
{
  WlzValueClampIntIntoRGBA(dst, src, count);
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of int.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToInt(int *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of short.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToShort(short *dst, short *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(short));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of WlzUByte.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToUByte(WlzUByte *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (WlzUByte )*src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of float.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToFloat(float *dst, short *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of double.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToDouble(double *dst, short *src,
					        int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of int.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToInt(int *dst, WlzUByte *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of short to a vector of RGBA.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyShortToRGBA(WlzUInt *dst, short *src, int count)
{
  WlzValueClampShortIntoRGBA(dst, src, count);
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of short.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToShort(short *dst, WlzUByte *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of WlzUByte.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToUByte(WlzUByte *dst, WlzUByte *src,
					 int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzUByte));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of float.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToFloat(float *dst, WlzUByte *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of double.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToDouble(double *dst, WlzUByte *src,
					  int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of WlzUByte to a vector of RGBA.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyUByteToRGBA(WlzUInt *dst, WlzUByte *src, int count)
{
  while(count-- > 0)
  {
    WLZ_RGBA_RGBA_SET(*dst, *src, *src, *src, 255);
    dst++;
    src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of int.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToInt(int *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of short.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToShort(short *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (short )WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of WlzUByte.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToUByte(WlzUByte *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (WlzUByte )WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of float.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToFloat(float *dst, float *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(float));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of double.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToDouble(double *dst, float *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of float to a vector of RGBA.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFloatToRGBA(WlzUInt *dst, float *src, int count)
{
  WlzValueClampFloatIntoRGBA(dst, src, count);
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of int.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToInt(int *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of short.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToShort(short *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (short )WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of WlzUByte.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToUByte(WlzUByte *dst, double *src,
					  int count)
{
  while(count-- > 0)
  {
    *dst++ = (WlzUByte )WLZ_NINT(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of float.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToFloat(float *dst, double *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (float )*src++;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of double.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToDouble(double *dst, double *src, int count)
{
  (void )memcpy(dst, src, count * sizeof(double));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of double to a vector of RGBA.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDoubleToRGBA(WlzUInt *dst, double *src, int count)
{
  WlzValueClampDoubleIntoRGBA(dst, src, count);
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of int - uses modulus
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToInt(int *dst, WlzUInt *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (int )WLZ_RGBA_MODULUS(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of short - uses modulus
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToShort(short *dst, WlzUInt *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (short )WLZ_RGBA_MODULUS(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of WlzUByte - uses modulus
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToUByte(WlzUByte *dst, WlzUInt *src,
				        int count)
{
  int	ival;

  while(count-- > 0)
  {
    ival = (int )WLZ_RGBA_MODULUS(*src);
    *dst++ = (WlzUByte )WLZ_CLAMP(ival, 0, 255);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of float - uses modulus
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToFloat(float *dst, WlzUInt *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = (float )WLZ_RGBA_MODULUS(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of double - uses modulus
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToDouble(double *dst, WlzUInt *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = WLZ_RGBA_MODULUS(*src);
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of RGBA to a vector of RGBA
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyRGBAToRGBA(WlzUInt *dst, WlzUInt *src, int count)
{
  while(count-- > 0)
  {
    *dst++ = *src++;
  }
}

/*!
* \return	void
* \ingroup      WlzValueUtils
* \brief	Copies a source vector to a destination vector,
*               where the source and destination types are any
*               combination of int, short, WlzUByte, float or double.
* \param	dst			Destination vector.
* \param	dstOff			Destination offset.
* \param	dstType			Destination type, eg WLZ_GREY_INT.
* \param	src			Source vector.
* \param	srcOff			Source offset.
* \param	srcType			Source type, eg WLZ_GREY_SHORT.
* \param	count			Number of vector elements to copy.
*/
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToInt(dst.inp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToShort(dst.shp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToUByte(dst.ubp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToFloat(dst.flp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
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
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToDouble(dst.dbp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
	  break;
      }
      break;
    case WLZ_GREY_RGBA:
      switch(srcType)
      {
	case WLZ_GREY_INT:
	  WlzValueCopyIntToRGBA(dst.rgbp + dstOff, src.inp + srcOff,
	  			  count);
	  break;
	case WLZ_GREY_SHORT:
	  WlzValueCopyShortToRGBA(dst.rgbp + dstOff, src.shp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_UBYTE:
	  WlzValueCopyUByteToRGBA(dst.rgbp + dstOff, src.ubp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_FLOAT:
	  WlzValueCopyFloatToRGBA(dst.rgbp + dstOff, src.flp + srcOff,
	  			    count);
	  break;
	case WLZ_GREY_DOUBLE:
	  WlzValueCopyDoubleToRGBA(dst.rgbp + dstOff, src.dbp + srcOff,
	  			     count);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueCopyRGBAToRGBA(dst.rgbp + dstOff, src.rgbp + srcOff,
	  			  count);
	  break;
	default:
	  break;
      }
      break;
    default:
      break;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D double vertices to a vector of
*		2D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDVertexToDVertex(WlzDVertex2 *dst,
					     WlzDVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzDVertex2));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D double vertices to a vector of
*		2D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDVertexToFVertex(WlzFVertex2 *dst,
					     WlzDVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = (float )(src->vtX);
    dst->vtY = (float )(src->vtY);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D double vertices to a vector of
*		2D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
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

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D float vertices to a vector of
*		2D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
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

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D float vertices to a vector of
*		2D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFVertexToFVertex(WlzFVertex2 *dst,
					     WlzFVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzFVertex2));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D float vertices to a vector of
*		2D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
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

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D int vertices to a vector of
*		2D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
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

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D int vertices to a vector of
*		2D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIVertexToFVertex(WlzFVertex2 *dst,
					     WlzIVertex2 *src,
					     int count)
{
  while(count-- > 0)
  {
    dst->vtX = (float )(src->vtX);
    dst->vtY = (float )(src->vtY);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 2D int vertices to a vector of
*		2D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIVertexToIVertex(WlzIVertex2 *dst,
					     WlzIVertex2 *src,
					     int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzIVertex2));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D double vertices to a vector of
*		3D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDVertexToDVertex3(WlzDVertex3 *dst,
					      WlzDVertex3 *src,
					      int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzDVertex3));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D double vertices to a vector of
*		3D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDVertexToFVertex3(WlzFVertex3 *dst,
					      WlzDVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = (float )(src->vtX);
    dst->vtY = (float )(src->vtY);
    dst->vtZ = (float )(src->vtZ);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D double vertices to a vector of
*		3D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyDVertexToIVertex3(WlzIVertex3 *dst,
					      WlzDVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = WLZ_NINT(src->vtX);
    dst->vtY = WLZ_NINT(src->vtY);
    dst->vtZ = WLZ_NINT(src->vtZ);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D float vertices to a vector of
*		3D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFVertexToDVertex3(WlzDVertex3 *dst,
					      WlzFVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = src->vtX;
    dst->vtY = src->vtY;
    dst->vtZ = src->vtZ;
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D float vertices to a vector of
*		3D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFVertexToFVertex3(WlzFVertex3 *dst,
					      WlzFVertex3 *src,
					      int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzFVertex3));
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D float vertices to a vector of
*		3D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyFVertexToIVertex3(WlzIVertex3 *dst,
					      WlzFVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = WLZ_NINT(src->vtX);
    dst->vtY = WLZ_NINT(src->vtY);
    dst->vtZ = WLZ_NINT(src->vtZ);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D int vertices to a vector of
*		3D double vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIVertexToDVertex3(WlzDVertex3 *dst,
					      WlzIVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = src->vtX;
    dst->vtY = src->vtY;
    dst->vtZ = src->vtZ;
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D int vertices to a vector of
*		3D float vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIVertexToFVertex3(WlzFVertex3 *dst,
					      WlzIVertex3 *src,
					      int count)
{
  while(count-- > 0)
  {
    dst->vtX = (float )(src->vtX);
    dst->vtY = (float )(src->vtY);
    dst->vtZ = (float )(src->vtZ);
    ++dst;
    ++src;
  }
}

/*!
* \return	void
* \ingroup	WlzValueUtils
* \brief	Copies a vector of 3D int vertices to a vector of
*		3D int vertices.
* \param	dst			Destination vector.
* \param	src			Source vector.
* \param	count			Number of vector elements to copy.
*/
void		WlzValueCopyIVertexToIVertex3(WlzIVertex3 *dst,
					      WlzIVertex3 *src,
					      int count)
{
  (void )memcpy(dst, src, count * sizeof(WlzIVertex3));
}

/*!
* \return	Woolz error code.
* \ingroup	WlzValueUtils
* \brief	Converts a single pixel value.
*               Source values are clamped to the the destination value
*               range.
* \param	dstPix			Destination pointer for pixel.
* \param	srcPix			Source pixel.
* \param	dstType			Destination type, eg WLZ_GREY_INT.
*/
WlzErrorNum	WlzValueConvertPixel(WlzPixelV *dstPix,
				     WlzPixelV srcPix,
				     WlzGreyType dstType)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzUInt	val;

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
	    dstPix->v.shv = (short )(srcPix.v.inv);
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampIntToUByte(&(srcPix.v.inv), 1);
	    dstPix->v.ubv = (WlzUByte )(srcPix.v.inv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = (float )(srcPix.v.inv);
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.inv;
	    break;
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
	    val = WLZ_CLAMP(srcPix.v.inv, 0, 255);
	    WLZ_RGBA_RGBA_SET(dstPix->v.rgbv, val, val, val, 255);
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
	    dstPix->v.ubv =  (WlzUByte )(srcPix.v.shv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = srcPix.v.shv;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.shv;
	    break;
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
	    val = WLZ_CLAMP(srcPix.v.shv, 0, 255);
	    WLZ_RGBA_RGBA_SET(dstPix->v.rgbv, val, val, val, 255);
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
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
	    val = srcPix.v.ubv;
	    WLZ_RGBA_RGBA_SET(dstPix->v.rgbv, val, val, val, 255);
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
	    dstPix->v.shv = (short )WLZ_NINT(srcPix.v.flv);
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampFloatToUByte(&(srcPix.v.flv), 1);
	    dstPix->v.ubv = (WlzUByte )(WLZ_NINT(srcPix.v.flv));
	    break;
	  case WLZ_GREY_FLOAT:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = srcPix.v.flv;
	    break;
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
	    val = WLZ_CLAMP(srcPix.v.flv, 0, 255);
	    WLZ_RGBA_RGBA_SET(dstPix->v.rgbv, val, val, val, 255);
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
	    dstPix->v.shv =  (short )WLZ_NINT(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToUByte(&(srcPix.v.dbv), 1);
	    dstPix->v.ubv = (WlzUByte )WLZ_NINT(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    WlzValueClampDoubleToFloat(&(srcPix.v.dbv), 1);
	    dstPix->v.flv = (float )(srcPix.v.dbv);
	    break;
	  case WLZ_GREY_DOUBLE:
	    *dstPix = srcPix;
	    break;
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
	    val = (WlzUInt )WLZ_CLAMP(srcPix.v.dbv, 0, 255);
	    WLZ_RGBA_RGBA_SET(dstPix->v.rgbv, val, val, val, 255);
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_GREY_RGBA:
	val = (WlzUInt )WLZ_RGBA_MODULUS(srcPix.v.rgbv);
        switch(dstType)
	{
	  case WLZ_GREY_INT:
	    dstPix->type = dstType;
	    dstPix->v.inv = val;
	    break;
	  case WLZ_GREY_SHORT:
	    dstPix->type = dstType;
	    dstPix->v.shv = (short )val;
	    break;
	  case WLZ_GREY_UBYTE:
	    dstPix->type = dstType;
	    dstPix->v.ubv = (WlzUByte )WLZ_CLAMP(val, 0, 255);
	    break;
	  case WLZ_GREY_FLOAT:
	    dstPix->type = dstType;
	    dstPix->v.flv = (float )val;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dstPix->type = dstType;
	    dstPix->v.dbv = val;
	    break;
	  case WLZ_GREY_RGBA:
	    dstPix->type = dstType;
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

/*!
* \return	Median value.
* \ingroup	WlzValueUtils
* \brief	Computes the median of the given vector (the vector's
*               values are overwritten), where the vector type is of
*               int values.
*               The algorithm is based on min-max elimination as used
*               by Peath A W. Median finding on a 3x3 Grid. Graphics
*               Gems, 171-175. Academic Press. 1990.
*               If the vector has less than 1 element 0 is returned.
*               In the trivial case of 1 or 2 elements, the first
*               element's value is returned.
*               If the vector has more than 3 elements:
*		<ul>
*		  <li>
*                 The given vector is partitioned into two portions:
*                 the 1st having ((n + 1) / 2) + 1 elements.
*		  </li>
*		  <li>
*                 While all the elements have not been examined:
*		  </li>
*		  <ul>
*		    <li>
*                   Find maximum and minimum values in the 1st
*		    </li>
*		    <li>
*                   Replace the min value with the first in the
*                   partition and the max with the last in the
*                   partition.
*		    </li>
*		    <li>
*                   Replace the first value in the first partition with
*                   the next value from the second partition.
*		    </li>
*		    <li>
*                   Reduce the length of the first partion by 1, ie
*                   omit the last value.
*		    <li>
*		  </ul>
*		  </li>
*		  <li>
*                 This leaves 3 values in the 1st partion, their median
*                 is found by direct comparison.
*		  </li>
*		</ul>
* \param	values			Vector of int values.
* \param	nVal			Number of values in vector.
*/
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

/*!
* \return	Median value.
* \ingroup	WlzValueUtils
* \brief	Computes the median of the given vector (the vector's
*               values are overwritten), where the vector type is of
*               double values.
*               The algorithm is based on min-max elimination as used
*               by Peath A W. Median finding on a 3x3 Grid. Graphics
*               Gems, 171-175. Academic Press. 1990.
*               If the vector has less than 1 element 0 is returned.
*               In the trivial case of 1 or 2 elements, the first
*               element's value is returned.
*               If the vector has more than 3 elements:
*		<ul>
*		  <li>
*                 The given vector is partitioned into two portions:
*                 the 1st having ((n + 1) / 2) + 1 elements.
*		  </li>
*		  <li>
*                 While all the elements have not been examined:
*		  </li>
*		  <ul>
*		    <li>
*                   Find maximum and minimum values in the 1st
*		    </li>
*		    <li>
*                   Replace the min value with the first in the
*                   partition and the max with the last in the
*                   partition.
*		    </li>
*		    <li>
*                   Replace the first value in the first partition with
*                   the next value from the second partition.
*		    </li>
*		    <li>
*                   Reduce the length of the first partion by 1, ie
*                   omit the last value.
*		    <li>
*		  </ul>
*		  </li>
*		  <li>
*                 This leaves 3 values in the 1st partion, their median
*                 is found by direct comparison.
*		  </li>
*		</ul>
* \param	values			Vector of double values.
* \param	nVal			Number of values in vector.
*/
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

/*!
* \return	Size of the given grey type.
* \ingroup	WlzValueUtils
* \brief	Computes the size of the given grey type.
* \param	gType			Given grey type.
*/
size_t		WlzValueSize(WlzGreyType gType)
{
  size_t	gSz = 0;

  switch(gType)
  {
    case WLZ_GREY_LONG:
      gSz = sizeof(long);
      break;
    case WLZ_GREY_INT:
      gSz = sizeof(int);
      break;
    case WLZ_GREY_SHORT:
      gSz = sizeof(short);
      break;
    case WLZ_GREY_UBYTE:
      gSz = sizeof(WlzUByte);
      break;
    case WLZ_GREY_FLOAT:
      gSz = sizeof(float);
      break;
    case WLZ_GREY_DOUBLE:
      gSz = sizeof(double);
      break;
    case WLZ_GREY_BIT:
      gSz = sizeof(WlzUByte);
      break;
    case WLZ_GREY_RGBA:
      gSz = sizeof(unsigned int);
      break;
    default:
      break;
  }
  return(gSz);
}
