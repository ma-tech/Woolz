#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzRGBChanRatio_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzRGBChanRatio.c
* \author       Bill Hill
* \date         February 2008
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
* \brief	Computes log ratio of RGB channels in a RGBA object.
* \ingroup	WlzArithmetic
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <Wlz.h>

static WlzObject 		*WlzRGBChanRatio2D(
				  WlzObject *rgbObj,
				  WlzRGBAColorChannel numC,
				  WlzRGBAColorChannel denC,
				  WlzRGBAColorChannel mulC,
				  int useMul,
				  int norm,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzRGBChanRatio3D(
				  WlzObject *rgbObj,
				  WlzRGBAColorChannel numC,
				  WlzRGBAColorChannel denC,
				  WlzRGBAColorChannel mulC,
				  int useMul,
				  int norm,
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
r = m \log(1 + \frac{n}{1 + d}).
\f]
*		where m is the multipler channel value or unity if not
*		used.
* 		This results in either an object with float values or if
* 		the normalise parameter is non-zero an object with unsigned
* 		byte values normalised to the range 0-255.
* 		The numerator and denominator channels must be red, green
* 		blue, yellow, magenta, cyan, hue, staturation, brightness,
* 		or grey (modulus).
* \param	rgbObj			The input RGBA object.
* \param	num			Channel for numerator in ratio.
* \param	den			Channel for denominator in ratio.
* \param	mul			Channel for multiplier value.
* \param	useMul			Multiplier used if non zero.
* \param	norm			Normalise the object values and
* 					return a ubyte object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzRGBChanRatio(WlzObject *rgbObj,
				WlzRGBAColorChannel num,
				WlzRGBAColorChannel den,
				WlzRGBAColorChannel mul,
				int useMul, int norm,
				WlzErrorNum *dstErr)
{
  WlzObject	*ratioObj = NULL;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_UBYTE;
  bgdV.v.ubv = 0;
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
        ratioObj = WlzRGBChanRatio2D(rgbObj, num, den, mul, useMul, norm,
				     &errNum);
	break;
      case WLZ_3D_DOMAINOBJ:
        ratioObj = WlzRGBChanRatio3D(rgbObj, num, den, mul, useMul, norm,
				     &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzSetBackground(ratioObj, bgdV);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(ratioObj);
    ratioObj = NULL;
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
* \brief	Computes log ratio of RGB channels in a RGBA object for each
* 		pixel using ratio \f$r\f$ with
\f[
r = m \log(1 + \frac{n}{1 + d}).
\f]
*		where m is the multipler channel value or unity if not
*		used.
* 		This results in either an object with float values or if
* 		the normalise parameter is non-zero an object with unsigned
* 		byte values normalised to the range 0-255.
* 		The numerator and denominator channels must be red, green
* 		blue, yellow, magenta, cyan, hue, staturation, brightness,
* 		or grey (modulus).
* \param	rgbObj			The input RGBA object.
* \param	num			Channel for numerator in ratio.
* \param	den			Channel for denominator in ratio.
* \param	mul			Channel for multiplier value.
* \param	useMul			Multiplier used if non zero.
* \param	norm			Normalise whole object values and
* 					return a ubyte object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzRGBChanRatio3D(WlzObject *rgbObj,
				WlzRGBAColorChannel num,
				WlzRGBAColorChannel den,
				WlzRGBAColorChannel mul,
				int useMul, int norm,
				WlzErrorNum *dstErr)
{
  int		idP,
  		idQ;
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
		*rtnObj = NULL,
  		*ratioObj = NULL;
  WlzPlaneDomain *pDom;
  WlzVoxelValues *vTab0,
  		 *vTab1;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pDom = rgbObj->domain.p;
  vTab0 = rgbObj->values.vox;
  vTab1 = WlzMakeVoxelValueTb(vTab0->type, vTab0->plane1, vTab0->lastpl,
                              vTab0->bckgrnd, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    dom.p = pDom;
    val.vox = vTab1;
    ratioObj = WlzMakeMain(rgbObj->type, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idP = pDom->plane1;
    while((errNum == WLZ_ERR_NONE) && (idP <= pDom->lastpl))
    {
      idQ = idP - pDom->plane1;;
      tObj0 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			  pDom->domains[idQ], vTab0->values[idQ],
			  NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	tObj1 = WlzRGBChanRatio2D(tObj0, num, den, mul, useMul, 0, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        vTab1->values[idQ] = WlzAssignValues(tObj1->values, NULL);
      }
      (void )WlzFreeObj(tObj0);
      (void )WlzFreeObj(tObj1);
      ++idP;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(norm)
    {
      if((errNum = WlzGreyNormalise(ratioObj, 1)) == WLZ_ERR_NONE)
      {
        rtnObj = WlzConvertPix(ratioObj, WLZ_GREY_UBYTE, &errNum);
      }
    }
    else
    {
      rtnObj = ratioObj;
      ratioObj = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}

/*!
* \return	Ratio object or NULL on error.
* \ingroup	WlzArithmetic
* \brief	Computes log ratio of RGB channels in a 2D RGBA object, see
* 		WlzRGBChanRatio().
* \param	rgbObj			The input 2D RGBA object.
* \param	num			Channel for numerator in ratio.
* \param	den			Channel for denominator in ratio.
* \param	mul			Channel for multiplier value.
* \param	useMul			Multiplier used if non zero.
* \param	norm			Normalise the object values and
* 					return a ubyte object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzRGBChanRatio2D(WlzObject *rgbObj,
				WlzRGBAColorChannel numC,
				WlzRGBAColorChannel denC,
				WlzRGBAColorChannel mulC,
				int useMul,
				int norm,
				WlzErrorNum *dstErr)
{
  int		cnt;
  double	num,
  		den,
		mul,
		ratio;
  WlzUInt	rgb;
  WlzValues	values;
  WlzPixelV	bgdV;
  WlzObject	*rtnObj = NULL,
  		*ratioObj = NULL;
  WlzGreyP	gP0,
  		gP1;
  WlzObjectType	vType;
  WlzGreyWSpace	gWSp0,
  		gWSp1;
  WlzIntervalWSpace iWSp0,
  		iWSp1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_FLOAT;
  bgdV.v.flv = 0.0f;
  vType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_FLOAT, NULL);
  values.v = WlzNewValueTb(rgbObj, vType, bgdV, &errNum);
  (void )WlzAssignValues(values, NULL);
  if(errNum == WLZ_ERR_NONE)
  {
    ratioObj = WlzAssignObject(
               WlzMakeMain(rgbObj->type, rgbObj->domain, values,
			   NULL, NULL, &errNum), NULL);
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
	((errNum = WlzNextGreyInterval(&iWSp0)) == WLZ_ERR_NONE) &&
	((errNum = WlzNextGreyInterval(&iWSp1)) == WLZ_ERR_NONE))
  {
    switch(gWSp0.pixeltype)
    {
      case WLZ_GREY_RGBA:
	gP0.rgbp = gWSp0.u_grintptr.rgbp;
	gP1.rgbp = gWSp1.u_grintptr.rgbp;
	cnt = iWSp0.rgtpos - iWSp0.lftpos + 1;
	while(cnt-- > 0)
	{
	  rgb = *(gP0.rgbp);
	  num = WlzRGBAChanGetValue(rgb, numC);
	  den = WlzRGBAChanGetValue(rgb, denC);
	  ratio = log(1.0 + (num / (den + 1.0)));
	  if(useMul)
	  {
	    mul = WlzRGBAChanGetValue(rgb, mulC);
	    ratio = ratio * mul;
	  }
	  ++(gP0.rgbp);
	  *(gP1.flp)++ = (float )ratio;
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
  if(errNum == WLZ_ERR_NONE)
  {
    if(norm)
    {
      if((errNum = WlzGreyNormalise(ratioObj, 1)) == WLZ_ERR_NONE)
      {
        rtnObj = WlzConvertPix(ratioObj, WLZ_GREY_UBYTE, &errNum);
      }
    }
    else
    {
      rtnObj = ratioObj;
      ratioObj = NULL;
    }
  }
  (void )WlzFreeValues(values);
  (void )WlzFreeObj(ratioObj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rtnObj);
    rtnObj = NULL;
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
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
    case WLZ_RGBA_CHANNEL_GREY:       /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_RED:        /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_GREEN:      /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_BLUE:	      /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_CYAN:       /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_MAGENTA:    /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_YELLOW:     /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_HUE:        /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_SATURATION: /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_BRIGHTNESS:
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
  int		b,
		c,
		h,
		k,
		m,
		s,
		y,
		red,
		green,
		blue,
		del,
		min,
		max;
  WlzUByte	val = 0;

  switch(chan)
  {
    case WLZ_RGBA_CHANNEL_GREY:
      val = (WlzUByte )WLZ_NINT((double )(WLZ_RGBA_MODULUS(rgba)) / sqrt(3.0));
      break;
    case WLZ_RGBA_CHANNEL_RED:
      val = (WlzUByte )WLZ_RGBA_RED_GET(rgba);
      break;
    case WLZ_RGBA_CHANNEL_GREEN:
      val = (WlzUByte )WLZ_RGBA_GREEN_GET(rgba);
      break;
    case WLZ_RGBA_CHANNEL_BLUE:
      val = (WlzUByte )WLZ_RGBA_BLUE_GET(rgba);
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
	  val = (WlzUByte )(c - k);
	  break;
        case WLZ_RGBA_CHANNEL_MAGENTA:
	  val = (WlzUByte )(m - k);
	  break;
        case WLZ_RGBA_CHANNEL_YELLOW:
	  val = (WlzUByte )(y - k);
	  break;
	default:
	  break;
      }
      break;
    case WLZ_RGBA_CHANNEL_HUE:        /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_SATURATION: /* FALLTHROUGH */
    case WLZ_RGBA_CHANNEL_BRIGHTNESS:
      red = WLZ_RGBA_RED_GET(rgba);
      green = WLZ_RGBA_GREEN_GET(rgba);
      blue = WLZ_RGBA_BLUE_GET(rgba);
      min = ALG_MIN3(red, green, blue);
      max = ALG_MAX3(red, green, blue);
      del = max - min;
      b = max;
      s = (max == 0)? 0: (int )(255.0 * del / max);
      switch(chan)
      {
        case WLZ_RGBA_CHANNEL_BRIGHTNESS:
	  val = (WlzUByte )b;
	  break;
        case WLZ_RGBA_CHANNEL_SATURATION:
	  val = (WlzUByte )s;
	  break;
        case WLZ_RGBA_CHANNEL_HUE:
	  if(s == 0)
	  {
	    h = -1;
	  }
	  else
	  {
	    if(max == red)
	    {
	      h = (green - blue) / del;
	    }
	    else if(max == green)
	    {
	      h = 2 + (blue - red) / del;
	    }
	    else
	    {
	      h = 4 + (red - green) / del;
	    }
	  }
	  h = (360 + (h * 60)) % 360;
	  val = (WlzUByte )h;
	  break;
        default:
	  break;
      }
      break;
    default:
      break;
  }
  return(val);
}
