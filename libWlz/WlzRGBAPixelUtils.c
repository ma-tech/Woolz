#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzRGBAPixelUtils.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri Jul 11 18:10:06 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzValuesUtils
* \brief        Utility functions for pixel and RGBA values.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
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
 should be tested since all values are valid (except grey-type UBYTE).
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

    case WLZ_RGBA_CHANNEL_BRIGHTNESS:
    case WLZ_RGBA_CHANNEL_HUE:
    case WLZ_RGBA_CHANNEL_SATURATION:
      /* need definitions here */
      val = 0.0;
      break;

    }
    break;

  case WLZ_GREY_BIT:
  default:
    errNum = WLZ_ERR_GREY_TYPE;
    break;
  }
  if( val < 0 ){
    val = -1.0;
    errNum = WLZ_ERR_GREY_DATA;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return val;
}
