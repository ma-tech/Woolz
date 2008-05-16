#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRGBAPixelUtils_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzRGBAPixelUtils.c
* \author       Richard Baldock
* \date         July 2003
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
* \brief	Utility functions for pixel and RGBA values.
* \ingroup	WlzValuesUtils
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>

/* function:     WlzRGBAPixelValue    */
/*! 
* \ingroup      WlzValuesUtils
* \brief        Calculate the pixel value for a given channel.
For grey-pixel types the colour channel values are set equal to the
 grey value i.e. the pixel is assumed to be (g,g,g). If the grey-channel
 is requested of a colour pixel the modulus is returned. For colour
 pixels the error return value is -1, for grey pixels the error return
 should be tested since all values are valid (except grey-type WlzUByte).
 Hue and saturation are zero for grey-pixel types.
*
* \return       requested value of pixel
* \param    pixVal	input pixel value structure
* \param    chan        requested pixel channel
* \param    dstErr	error destination
* \par      Source:
*                WlzRGBAPixelUtils.c
*/
double WlzRGBAPixelValue(
  WlzPixelV		pixVal,
  WlzRGBAColorChannel	chan,
  WlzErrorNum		*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  double	val=-1;
  int		col[3];

  switch( pixVal.type ){
  case WLZ_GREY_INT:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      val = 0.0;
      break;

    default:
      val = pixVal.v.inv;
      break;
    }
    break;

  case WLZ_GREY_SHORT:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      val = 0.0;
      break;

    default:
      val = pixVal.v.shv;
      break;
    }
    break;

  case WLZ_GREY_UBYTE:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      val = 0.0;
      break;

    default:
      val = pixVal.v.ubv;
      break;
    }
    break;

  case WLZ_GREY_FLOAT:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      val = 0.0;
      break;

    default:
      val = pixVal.v.flv;
      break;
    }
    break;

  case WLZ_GREY_DOUBLE:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      val = 0.0;
      break;

    default:
      val = pixVal.v.dbv;
      break;
    }
    break;

  case WLZ_GREY_RGBA:
    switch( chan ){
    case WLZ_RGBA_CHANNEL_GREY:
      /* ????? */
      val = WLZ_RGBA_MODULUS(pixVal.v.rgbv);
      break;

    case WLZ_RGBA_CHANNEL_RED:
      val = WLZ_RGBA_RED_GET(pixVal.v.rgbv);
      break;

    case WLZ_RGBA_CHANNEL_GREEN:
      val = WLZ_RGBA_GREEN_GET(pixVal.v.rgbv);
      break;

    case WLZ_RGBA_CHANNEL_BLUE:
      val = WLZ_RGBA_BLUE_GET(pixVal.v.rgbv);
      break;

    case WLZ_RGBA_CHANNEL_HUE:
      col[0] = WLZ_RGBA_RED_GET(pixVal.v.rgbv);
      col[1] = WLZ_RGBA_GREEN_GET(pixVal.v.rgbv);
      col[2] = WLZ_RGBA_BLUE_GET(pixVal.v.rgbv);
      WlzRGBAConvertRGBToHSV_UBYTENormalised(col);
      val = col[0];
      break;

    case WLZ_RGBA_CHANNEL_SATURATION:
      col[0] = WLZ_RGBA_RED_GET(pixVal.v.rgbv);
      col[1] = WLZ_RGBA_GREEN_GET(pixVal.v.rgbv);
      col[2] = WLZ_RGBA_BLUE_GET(pixVal.v.rgbv);
      WlzRGBAConvertRGBToHSV_UBYTENormalised(col);
      val = col[1];
      break;

    case WLZ_RGBA_CHANNEL_BRIGHTNESS:
      col[0] = WLZ_RGBA_RED_GET(pixVal.v.rgbv);
      col[1] = WLZ_RGBA_GREEN_GET(pixVal.v.rgbv);
      col[2] = WLZ_RGBA_BLUE_GET(pixVal.v.rgbv);
      WlzRGBAConvertRGBToHSV_UBYTENormalised(col);
      val = col[2];
      break;

    case WLZ_RGBA_CHANNEL_CYAN:
      val = (WLZ_RGBA_BLUE_GET(pixVal.v.rgbv) +
	     WLZ_RGBA_GREEN_GET(pixVal.v.rgbv)) / 2;
      break;

    case WLZ_RGBA_CHANNEL_MAGENTA:
      val = (WLZ_RGBA_BLUE_GET(pixVal.v.rgbv) +
	     WLZ_RGBA_RED_GET(pixVal.v.rgbv)) / 2;
      break;

    case WLZ_RGBA_CHANNEL_YELLOW:
      val = (WLZ_RGBA_RED_GET(pixVal.v.rgbv) +
	     WLZ_RGBA_GREEN_GET(pixVal.v.rgbv)) / 2;
      break;

    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;

    }
    break;

  case WLZ_GREY_BIT:
  default:
    errNum = WLZ_ERR_GREY_TYPE;
    break;
  }
  if(errNum == WLZ_ERR_NONE){
    if( val < 0 ){
      val = -1.0;
      errNum = WLZ_ERR_GREY_DATA;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return val;
}

void WlzRGBAConvertRGBToHSV_UBYTENormalised(
  int		*col)
{
  int	h, s, b;
  int	max, min;

  /* algorithm from Foley, van Dam, Feiner, Hughes,
     Computer Graphics. Modified so each value is in the range
     [0,255]. If saturation is zero then the hue is undefined
     and in this case set to zero */
  max = WLZ_MAX(col[0],col[1]);
  max = WLZ_MAX(max, col[2]);
  min = WLZ_MIN(col[0],col[1]);
  min = WLZ_MIN(min, col[2]);

  b = max;
  if( max > 0 ){
    s = (max - min) * 255 / max;
  }
  else {
    s = 0;
  }
  if( s == 0 ){
    h = 0;
  }
  else {
    if( col[0] == max ){
      h = (col[1] - col[2]) * 42.5 / (max - min);
    }
    else if( col[1] == max ){
      h = 85 + (col[2] - col[0]) * 42.5 / (max - min);
    }
    else if( col[2] == max ){
      h = 170 + (col[0] - col[1]) * 42.5 / (max - min);
    }
    if( h < 0 ){
      h += 255;
    }
  }
  col[0] = h;
  col[1] = s;
  col[2] = b;

  return;
}
