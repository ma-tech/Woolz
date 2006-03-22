#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzSobel_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzSobel.c
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
* \brief	Implements a Sobel edge detection filter.
* \ingroup	WlzValueFilters
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	Sobel filtered object.
* \ingroup	WlzValueFilters
* \brief	Applies a Sobel edge detection filter to the given
*               woolz object.
* \param	srcObj			Given source object.
* \param	hFlag			Apply horizontal edge kernel if
* 					non-zero.
* \param	vFlag			Apply vertical edge kernel if
* 					non-zero.
* \param	wlzErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzSobel(WlzObject *srcObj, int hFlag, int vFlag,
			  WlzErrorNum *wlzErr)
{
  WlzObject	*dstObj = NULL,
  		*objH = NULL,
		*objV = NULL;
  WlzGreyType	gType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	sMaskHI[] =
  {
    -1, -2, -1,
     0,  0,  0,
     1,  2,  1
  },
  		sMaskVI[] =
  {
    -1,  0,  1,
    -2,  0,  2,
    -1,  0,  1
  };
  const double	sMaskHD[] =
  {
    -1.0, -2.0, -1.0,
     0.0,  0.0,  0.0,
     1.0,  2.0,  1.0
  },
  		sMaskVD[] =
  {
    -1.0,  0.0,  1.0,
    -2.0,  0.0,  2.0,
    -1.0,  0.0,  1.0
  };
  WlzConvolution sConv;

  WLZ_DBG((WLZ_DBG_LVL_1),
	  ("WlzSobel FE 0x%lx %d %d 0x%lx\n",
	   (unsigned long )srcObj, hFlag, vFlag, (unsigned long)wlzErr));
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
    sConv.linkcount = 0;
    sConv.xsize = 3;
    sConv.ysize = 3;
    sConv.cv = NULL;
    sConv.divscale = 1;
    sConv.offset = 0;
    sConv.modflag = 1;
    switch(gType)
    {
      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
        sConv.type = WLZ_CONVOLVE_INT;
        break;
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
        sConv.type = WLZ_CONVOLVE_FLOAT;
	break;
      case WLZ_GREY_RGBA: /* RGBA to be done RAB */
      default:
	errNum = WLZ_ERR_GREY_DATA;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(hFlag)
    {
      sConv.cv = (sConv.type == WLZ_CONVOLVE_INT)?
      		 (int *)(&(sMaskHI[0])):
      		 (int *)(&(sMaskHD[0]));
      objH = WlzConvolveObj(srcObj, &sConv, 1, &errNum);
    }
    if((errNum == WLZ_ERR_NONE) && vFlag)
    {
      sConv.cv = (sConv.type == WLZ_CONVOLVE_INT)?
      		 (int *)(&(sMaskVI[0])):
      		 (int *)(&(sMaskVD[0]));
      objV = WlzConvolveObj(srcObj, &sConv, 1, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(objH && objV)
      {
	dstObj = WlzImageArithmetic(objV, objH, WLZ_BO_ADD, 1, &errNum);
	(void )WlzFreeObj(objV);
	objV = NULL;
	(void )WlzFreeObj(objH);
	objH = NULL;
      }
      else if(objH && (vFlag == 0))
      {
	dstObj = objH;
	objH = NULL;
      }
      else if(objV && (hFlag == 0))
      {
	dstObj = objV;
	objV = NULL;
      }
    }
    else
    {
      (void )WlzFreeObj(objH);
      (void )WlzFreeObj(objV);
    }
  }
  if(wlzErr)
  {
    *wlzErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzSobel FX 0x%lx\n",
	   (unsigned long )dstObj));
  return(dstObj);
}
