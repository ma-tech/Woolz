#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzGreyNormalise.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Fri May 23 07:45:49 2003
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
* \ingroup      WlzValuesFilters
* \brief         Normalise the grey-values of a woolz object to the
 range 0-255. Colour images have each channel independently normalised.
For proportional normalisation of colour use WlzGreySetRange directly.
To determine the modulus range use WlzRGBAModulusRange().
*
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
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
  WlzErrorNum	wlzErrno=WLZ_ERR_NONE;

  /* get then set the grey-range of the object */
  wlzErrno = WlzGreyRange(obj, &min, &max);
  Min = min;
  Max = max;
  if( wlzErrno == WLZ_ERR_NONE ){
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
    }
    wlzErrno = WlzGreySetRange(obj, min, max, Min, Max);
  }

  return wlzErrno;
}
