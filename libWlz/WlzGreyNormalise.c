#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyNormalise_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzGreyNormalise.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Normalises the grey-values of an object to the range
* 		[0 - 255].
* 		Colour images have each channel independently normalised.
* 		For proportional normalisation of colour use WlzGreySetRange
* 		directly. To determine the modulus range use
* 		WlzRGBAModulusRange().
* \ingroup	WlzValuesFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/* function:     WlzGreyNormalise    */
/*! 
* \ingroup      WlzValuesFilters
* \brief        
*
* \return       Woolz error number
* \param    obj	Input object to be normalised i.e the grey values reset to
fill the range 0-255. Colour values are independently reset which will
change colour balance. Use WlzGreySetRange directly to avoid this.
Note grey-values are reset in place and not copied.
* \par      Source:
*                WlzGreyNormalise.c
*/
WlzErrorNum WlzGreyNormalise(
  WlzObject	*obj)
{
  WlzPixelV	min, max, Min, Max;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* get then set the grey-range of the object */
  errNum = WlzGreyRange(obj, &min, &max);
  Min = min;
  Max = max;
  if( errNum == WLZ_ERR_NONE ){
    switch( min.type ){
    case WLZ_GREY_INT:
      Min.v.inv = 0;
      Max.v.inv = 255;
      break;
    case WLZ_GREY_SHORT:
      Min.v.shv = 0;
      Max.v.shv = 255;
      break;
    case WLZ_GREY_UBYTE:
      Min.v.ubv = 0;
      Max.v.ubv = 255;
      break;
    case WLZ_GREY_FLOAT:
      Min.v.flv = 0.0;
      Max.v.flv = 255.0;
      break;
    case WLZ_GREY_DOUBLE:
      Min.v.dbv = 0.0;
      Max.v.dbv = 255.0;
      break;
    case WLZ_GREY_RGBA:
      Min.v.rgbv = 0xff000000;
      Max.v.rgbv = 0xffffffff;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzGreySetRange(obj, min, max, Min, Max);
    }
  }

  return(errNum);
}
