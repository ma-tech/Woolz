#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzLaplacian_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzLaplacian.c
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
* \brief	Implements a Laplacian edge enhancement filter.
* \ingroup	WlzValueFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Laplacian filtered object.
* \ingroup      WlzValueFilters
* \brief	Applies a Laplacian edge enhancement filter to the
*               given woolz object.
* \param	srcObj			Given source object.
* \param	kSize			Kernel size, must be 3, 5 or 7.
* \param	newObjFlag		If zero the convolution is done
*                                       in place, else a new object is
*                                       created.
* \param	modFlag			Take the absolute value of the
*                                       convolution if non zero (this
*                                       is always the case for a WlzUByte
*                                       object).
* \param	dstErr			Destination error pointer, may
*                                       be NULL.
*/
WlzObject	*WlzLaplacian(WlzObject *srcObj, int kSize,
			      int newObjFlag, int modFlag,
			      WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	lMask3I[] = 		     /* 3 by 3 with central weight 1 */
  {
    -1, -1, -1,
    -1,  9, -1,
    -1, -1, -1
  },
  		lMask5I[] =		     /* 5 by 5 with central weight 4 */
  {
    -1, -1, -1, -1, -1,
    -1,  0,  0,  0, -1,
    -1,  0, 20,  0, -1,
    -1,  0,  0,  0, -1,
    -1, -1, -1, -1, -1
  },
  		lMask7I[] = 		     /* 7 by 7 with central weight 4 */
  {
     0,  0,  0, -1,  0,  0,  0,
     0, -1,  0,  0,  0, -1,  0,
     0,  0,  0,  1,  0,  0,  0,
    -1,  0,  1,  8,  1,  0, -1,
     0,  0,  0,  1,  0,  0,  0,
     0, -1,  0,  0,  0, -1,  0,
     0,  0,  0, -1,  0,  0,  0
  };
  const double	lMask3D[] = 		     /* 3 by 3 with central weight 1 */
  {
    -1.0, -1.0, -1.0,
    -1.0,  9.0, -1.0,
    -1.0, -1.0, -1.0
  },
  		lMask5D[] =		     /* 5 by 5 with central weight 4 */
  {
    -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0,  0.0,  0.0,  0.0, -1.0,
    -1.0,  0.0, 20.0,  0.0, -1.0,
    -1.0,  0.0,  0.0,  0.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0
  },
  		lMask7D[] = 		     /* 7 by 7 with central weight 4 */
  {
     0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0,
     0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
     0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,
    -1.0,  0.0,  1.0,  8.0,  1.0,  0.0, -1.0,
     0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,
     0.0, -1.0,  0.0,  0.0,  0.0, -1.0,  0.0,
     0.0,  0.0,  0.0, -1.0,  0.0,  0.0,  0.0
  };
  WlzConvolution lConv;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzLaplacian FE %p %d %d %d %p\n",
	   srcObj, kSize, newObjFlag, modFlag, dstErr));
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->type != WLZ_2D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    gType = WlzGreyTableTypeToGreyType(srcObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    lConv.linkcount = 0;
    lConv.xsize = kSize;
    lConv.ysize = kSize;
    lConv.cv = NULL;
    lConv.divscale = 1;
    lConv.offset = 0;
    lConv.modflag = gType == WLZ_GREY_UBYTE? 1: modFlag;
    switch(gType)
    {
      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
        lConv.type = WLZ_CONVOLVE_INT;
	switch(kSize)
	{
	  case 3:
    	    lConv.divscale = 1;
	    lConv.cv = (int *)(&(lMask3I[0]));
	    break;
	  case 5:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask5I[0]));
	    break;
	  case 7:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask7I[0]));
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
        break;
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
        lConv.type = WLZ_CONVOLVE_FLOAT;
	switch(kSize)
	{
	  case 3:
    	    lConv.divscale = 1;
	    lConv.cv = (int *)(&(lMask3D[0]));
	    break;
	  case 5:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask5D[0]));
	    break;
	  case 7:
    	    lConv.divscale = 4;
	    lConv.cv = (int *)(&(lMask7D[0]));
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
	break;
      case WLZ_GREY_RGBA: /* RGBA to be done RAB */
      default:
	errNum = WLZ_ERR_GREY_DATA;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzConvolveObj(srcObj, &lConv, newObjFlag, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzLaplacian FX %p\n",
	   dstObj));
  return(dstObj);
}
