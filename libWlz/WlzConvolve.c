#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzConvolve_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzConvolve.c
* \author       Jim Piper, Bill Hill
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
* \brief	Functions for convolving the values of domain objects.
* \ingroup	WlzValuesFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Convolved pixel value.
* \ingroup	WlzValuesFilters
* \brief	Performs a general space-domain convolution for
*               WlzSeqPar().
* \param	spWSpace		Work space data structure from
* 					WlzSeqPar().
* \param	spData			Data passed on by WlzSeqPar(),
*					which is used to pass the convolution
*					structure.
*/
int 		WlzConvolveSeqParFn(WlzSeqParWSpace *spWSpace,
				    void *spData)
{
  WlzConvolution *conv;
  int		*data,
  		*mask;
  int		idX,
  		idY,
		lx,
		ly,
		convPixVal = 0;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
  	  ("WlzConvolveSeqParFn FE 0x%lx 0x%lx\n",
	   (unsigned long )spWSpace, (unsigned long )spData));
  if(spWSpace && ((conv = (WlzConvolution *)spData) != NULL) && 
      conv->divscale)
  {
    mask = conv->cv;
    lx = (conv->xsize - 1) / 2; 
    ly = (conv->ysize - 1) / 2;
    WLZ_DBG((WLZ_DBG_LVL_3), 
	    ("WlzConvolveSeqParFn 01 %d %d\n",
	     lx, ly));
    for(idY = -ly; idY <= ly; ++idY)
    {
      data = spWSpace->adrptr[idY] - lx;
      for(idX = -lx; idX <= lx; ++idX)
      {
	WLZ_DBG((WLZ_DBG_LVL_3), 
		("WlzConvolveSeqParFn 02 %d %d\n",
		 *data, *mask));
	convPixVal += *data++ * *mask++;
      }
    }
    convPixVal = (convPixVal / conv->divscale) + conv->offset;
    if(conv->modflag)
    {
      convPixVal = abs(convPixVal);
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
  	  ("WlzConvolveSeqParFn FX %d\n",
	   convPixVal));
  return (convPixVal);
}

/*!
* \return	Convolved object or NULL on error.
* \ingroup	WlzValuesFilters
* \brief	Performs a general space-domain convolution using WlzSeqPar().
*		Only objects with WLZ_EMPTY_OBJ and WLZ_2D_DOMAINOBJ
*               types are valid. WLZ_2D_DOMAINOBJ ojects must have
*               non null domain and values fields, and only integral
*               values (ie int, short or UBYTE) are valid.
* \param	inObj			Given object.
* \param	conv			Convolution data structure.
* \param	newObjFlag		If zero the convolution is done
*					in place, else a new object is created.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzConvolveObj(WlzObject *inObj, WlzConvolution *conv,
			        int newObjFlag, WlzErrorNum *dstErr)
{
  int		bkgIntVal,
  		convSize;
  WlzObject 	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzPixelV	bkgVal;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzConvolveObj FE 0x%lx 0x%lx %d\n",
	   (unsigned long )inObj, (unsigned long )conv,
	   newObjFlag));
  if((inObj == NULL) || (conv == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(conv->type != WLZ_CONVOLVE_INT)
  {
    errNum = WLZ_ERR_GREY_DATA;
  }
  else
  {
    switch(inObj->type)
    {
      case WLZ_EMPTY_OBJ:
	outObj = (newObjFlag)? WlzMakeEmpty(&errNum): inObj;
	break;
      case WLZ_2D_DOMAINOBJ:
	bkgVal = WlzGetBackground(inObj, &errNum);
	switch(bkgVal.type)
	{
	  case WLZ_GREY_INT:
	    bkgIntVal = bkgVal.v.inv;
	    break;
	  case WLZ_GREY_SHORT:
	    bkgIntVal = bkgVal.v.shv;
	    break;
	  case WLZ_GREY_UBYTE:
	    bkgIntVal = bkgVal.v.ubv;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_DATA;
	    break;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  convSize = (WLZ_MAX(conv->xsize, conv->ysize) - 1) / 2;
	  outObj = WlzSeqPar(inObj, newObjFlag, 0, WLZ_RASTERDIR_ILIC,
			     convSize, bkgIntVal,
			     (void *)conv, WlzConvolveSeqParFn, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if( dstErr )
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
  	  ("WlzConvolveObj FX 0x%lx\n",
	   (unsigned long )outObj));
  return(outObj);
}

/*!
* \return	Sum of convolution values.
* \ingroup	WlzValuesFilters
* \brief	Calculates the sum of the values in the convolution object.
* \param	conv			Given convolution object.
*/
int		WlzConvolutionSum(WlzConvolution *conv)
{
  int		count,
  		sum = 0;
  int		*values;

  if(conv && ((count = conv->xsize * conv->ysize) > 0) &&
     ((values = conv->cv) != NULL))
  {
    sum = *values;
    while(--count > 0)
    {
      sum += *++values;
    }
  }
  return(sum);
}

/*!
* \return	Convolution sum.
* \ingroup	WlzValuesFilters
* \brief	Normalises a convolution object by adjusting the scaling such
* 		that the components effectively sum to unity.
* \param	conv			Given convolution object.
*/
int		WlzConvolutionNormalise(WlzConvolution *conv)
{
  int		sum = 0;

  if(conv)
  {
    if((sum = WlzConvolutionSum(conv)) == 0)
    {
      conv->divscale = 1;
    }
    else
    {
      conv->divscale = sum;
    }
  }
  return(sum);
}
