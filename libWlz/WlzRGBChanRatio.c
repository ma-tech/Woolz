#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzRGBChanRatio_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzRGBChanRatio.c
* \author       Bill Hill
* \date         February 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Computes log ratio of RGB channels in a RGBA object.
* \ingroup	WlzArithmetic
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <Wlz.h>

static WlzObject 		*WlzRGBChanRatio2D(
				  WlzObject *rgbObj,
				  WlzRGBAColorChannel numC,
				  WlzRGBAColorChannel denC,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzRGBAChanValid(
				  WlzRGBAColorChannel chan);
static WlzUByte			WlzRGBAChanGetValue(
				  WlzUInt rgba,
				  WlzRGBAColorChannel chan);

/*!
* \return	Ratio object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Computes log ratio of RGB channels in a RGBA object for each
* 		pixel using ratio \f$r\f$ with
\f[
r = 46 \log(1 + \frac{n}{1 + d}).
\f]
* 		This results in an object normalised to the range 0-255.
* 		The numerator and denominator channels must be red, green
* 		blue, yellow, magenta or cyan.
* \param	rgbObj			The input RGBA object.
* \param	num			Channel for numerator in ratio.
* \param	den			Channel for denominator in ratio.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzRGBChanRatio(WlzObject *rgbObj,
				WlzRGBAColorChannel num,
				WlzRGBAColorChannel den,
				WlzErrorNum *dstErr)
{
  WlzObject	*ratioObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(rgbObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(rgbObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(rgbObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    errNum = WlzRGBAChanValid(num);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzRGBAChanValid(den);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(rgbObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        ratioObj = WlzRGBChanRatio2D(rgbObj, num, den, &errNum);
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
  return(ratioObj);
}

/*!
* \return	Ratio object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Computes log ratio of RGB channels in a 2D RGBA object, see
* 		WlzRGBChanRatio().
* \param	rgbObj			The input 2D RGBA object.
* \param	num			Channel for numerator in ratio.
* \param	den			Channel for denominator in ratio.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzRGBChanRatio2D(WlzObject *rgbObj,
				WlzRGBAColorChannel numC,
				WlzRGBAColorChannel denC,
				WlzErrorNum *dstErr)
{
  int		cnt;
  double	num,
  		den,
		ratio;
  WlzValues	values;
  WlzPixelV	bgdV;
  WlzObject	*ratioObj = NULL;
  WlzObjectType	vType;
  WlzGreyWSpace	gWSp0,
  		gWSp1;
  WlzIntervalWSpace iWSp0,
  		iWSp1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, NULL);
  values.v = WlzNewValueTb(rgbObj, vType, bgdV, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    ratioObj = WlzMakeMain(rgbObj->type, rgbObj->domain, values,
			   rgbObj->plist, rgbObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(rgbObj, &iWSp0, &gWSp0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(ratioObj, &iWSp1, &gWSp1);
  }
  while((errNum == WLZ_ERR_NONE) &&
	((errNum = WlzNextGreyInterval(&iWSp0)) == WLZ_ERR_NONE))
  {
    (void )WlzNextGreyInterval(&iWSp1);
    switch(gWSp0.pixeltype)
    {
      case WLZ_GREY_RGBA:
	cnt = iWSp0.rgtpos - iWSp0.lftpos + 1;
	while(cnt-- > 0)
	{
	  num = WlzRGBAChanGetValue(*(gWSp0.u_grintptr.rgbp), numC);
	  den = WlzRGBAChanGetValue(*(gWSp0.u_grintptr.rgbp), denC);
	  ratio = 46.0 * log(1.0 + (num / (den + 1.0)));
	  *(gWSp1.u_grintptr.ubp) = WLZ_NINT(ratio);
	  ++(gWSp0.u_grintptr.rgbp);
	  ++(gWSp1.u_grintptr.ubp);
	}
	break;
      default:
	errNum = WLZ_ERR_GREY_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(ratioObj)
    {
      WlzFreeObj(ratioObj);
      ratioObj = NULL;
    }
    else if(values.core != NULL)
    {
      WlzFreeValues(values);
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ratioObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Checks that the given channel is valid for
* 		calling WlzRGBAChanGetValue().
* \param	chan			Given colour channel.
*/
static WlzErrorNum WlzRGBAChanValid(WlzRGBAColorChannel chan)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(chan)
  {
    case WLZ_RGBA_CHANNEL_RED:        /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_GREEN:      /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_BLUE:	      /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_CYAN:       /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_MAGENTA:    /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_YELLOW:
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzArithmetic
* \brief	Gets the given colour channel from the given RGBA value.
* \param	rgba			RGBA value.
* \param	chan			Given colour channel.
*/
static WlzUByte	WlzRGBAChanGetValue(WlzUInt rgba, WlzRGBAColorChannel chan)
{
  WlzUByte	c,
		y,
		m,
		k,
		val = 0;

  switch(chan)
  {
    case WLZ_RGBA_CHANNEL_RED:
      val = WLZ_RGBA_RED_GET(rgba);
      break;
    case WLZ_RGBA_CHANNEL_GREEN:
      val = WLZ_RGBA_GREEN_GET(rgba);
      break;
    case WLZ_RGBA_CHANNEL_BLUE:
      val = WLZ_RGBA_BLUE_GET(rgba);
      break;
    case WLZ_RGBA_CHANNEL_CYAN:     /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_MAGENTA:  /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_YELLOW:
      c = 255 - WLZ_RGBA_RED_GET(rgba);
      m = 255 - WLZ_RGBA_GREEN_GET(rgba);
      y = 255 - WLZ_RGBA_BLUE_GET(rgba);
      k = ALG_MIN3(c, m, y);
      switch(chan)
      {
        case WLZ_RGBA_CHANNEL_CYAN:
	  val = c - k;
	  break;
        case WLZ_RGBA_CHANNEL_MAGENTA:
	  val = m - k;
	  break;
        case WLZ_RGBA_CHANNEL_YELLOW:
	  val = y - k;
	  break;
	default:
	  break;
      }
      break;
    default:
      break;
  }
}
