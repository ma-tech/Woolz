#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFRaw_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFRaw.c
* \author       Richard Baldock
* \date         April 2004
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
* \brief        Read/write raw image data.
* \ingroup	WlzExtFF
*/

#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads a raw 'brick of bytes' object into a Woolz object.
* \param	fp			Given file stream.
* \param	numKols			Number of columns, x size.
* \param	numLines		Number of lines, y size.
* \param	numPlanes		Number of planes, z size.
* \param	bytesPerPixel		Number of bytes per pixel.
* \param	bitsPerChannel		Number of bits per channel.
* \param	numChannels		Number of channels.
* \param	byteOrder		Byte order.
* \param	dstErr			Destination pointer for error code,
					may be NULL.
*/
WlzObject *WlzExtFFReadObjRaw(
  FILE		*fp,
  int		numKols,
  int		numLines,
  int		numPlanes,
  int		bytesPerPixel,
  int		bitsPerChannel,
  int		numChannels,
  int		byteOrder,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj=NULL;
  WlzGreyPtr	*wlzData;
  int		wlzSize;
  char		rawData[8];
  int		numBytes;
  WlzGreyType	greyType;
  WlzObjectType	objType, subObjType;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check inputs */
  if((fp == NULL) ||
     (width <= 0) ||
     (height <= 0) ||
     (bytesPerPixel <= 0) ||
     (bitsPerChannel <= 0) ||
     (bitsPerChannel > 32) ||
     (numChannels <= 0)
    ){
    errNum = WLZ_ERR_PARAM_DATA;
  } 

  /* width, height and bytesPerPixel are sufficient to read the data */
  if( errNum == WLZ_ERR_NONE ){
    if( rawData = (char *) AlcMalloc(numBytes) ){
      if( fread(rawData, 1, numBytes, fp) != numBytes ){
	AlcFree(rawData);
	errNum - WLZ_ERR_READ_INCOMPLETE;
      }
    }
    else {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  /* now work out a woolz grey type (or type for a compound image).
     Note we assume the data is unsigned integer, if it is float we are lost */
  if( errNum ==WLZ_ERR_NONE ){
    switch( numChannels ){
    case 1:
      if( bitsPerChannel <= 8 ){
	greyType = WLZ_GREY_UBYTE;
      }
      else if( bitsPerChannel <= 16 ){
	greyType = WLZ_GREY_SHORT;
      }
      else if( bitsPerChannel <= 32 ){
	greyType = WLZ_GREY_INT;
      }
      objType = (numPlanes > 1)?WLZ_3D_DOMAINOBJ:WLZ_2D_DOMAINOBJ;
      break;

    case 2:
    default:
      if( bitsPerChannel <= 8 ){
	greyType = WLZ_GREY_UBYTE;
      }
      else if( bitsPerChannel <= 16 ){
	greyType = WLZ_GREY_SHORT;
      }
      else if( bitsPerChannel <= 32 ){
	greyType = WLZ_GREY_INT;
      }
      objType = WLZ_COMPOUND_ARR_1;
      subObjType = (numPlanes > 1)?WLZ_3D_DOMAINOBJ:WLZ_2D_DOMAINOBJ;
      break;
      
    case 3: /* this could be RGB but only if 8-bits per channel */
      objType = WLZ_COMPOUND_ARR_1;
      subObjType = (numPlanes > 1)?WLZ_3D_DOMAINOBJ:WLZ_2D_DOMAINOBJ;
      if( bitsPerChannel <= 8 ){
	greyType = WLZ_GREY_RGBA;
	objType = (numPlanes > 1)?WLZ_3D_DOMAINOBJ:WLZ_2D_DOMAINOBJ;
      }
      else if( bitsPerChannel <= 16 ){
	greyType = WLZ_GREY_SHORT;
      }
      else if( bitsPerChannel <= 32 ){
	greyType = WLZ_GREY_INT;
      }
      break;
    }

    /* determine woolz size */
    switch( greyType ){
    case WLZ_GREY_UBYTE:
      wlzSize = sizeof(WlzUByte);
      break;
    case WLZ_GREY_SHORT:
      wlzSize = sizeof(short);
      break;
    case WLZ_GREY_INT:
    case WLZ_GREY_RGBA:
      wlzSize = sizeof(int);
      break;
    }
  }

  /* byte swap if necessary */
#if defined (__sparc) || defined (__mips) || defined (__ppc)
#endif( /* __sparc || __mips || __ppc */

#if defined (__sparc) || defined (__mips)
  if( byteOrder ){
    swapFlg = 1;
  }
#endif( /* __sparc || __mips */
#if defined (__x86) || defined (__alpha)
  if( !byteOrder ){
    swapFlg = 1;
  }
#endif( /* __x86 || __alpha */

  /* allocate space for woolz object data and read data */
  if( errNum == WLZ_ERR_NONE ){
    rawData = (WlzUByte *) AlcMalloc(bytesPerPixel);
    switch( objType ){
    case WLZ_3D_DOMAINOBJ:
    case WLZ_2D_DOMAINOBJ:
      wlzData = (WlzGreyPtr *) AlcMalloc(sizeof(WlzGreyPtr));
      numBytes = numKols * numLines * numPlanes * wlzSize;
      if( wlzData[0].ubp = (unsigned char *) AlcMalloc(numBytes) ){
	if( wlzSize == bytesPerPixel ){
	  if( fread(wlzData[0].ubp, 1, numBytes, fp) != numBytes ){
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	}
	else {
	  for(i=0; i < numKols*numLines*numPlanes; i++){
	    if( fread(rawData, bytesPerPixel, 1, fp) != 1 ){
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    }
	    switch( greyType ){
	      /********** bloop bloop *************/
	    }
	  }
	}
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      AlcFree(wlzData);
      break;

    case WLZ_COMPOUND_ARRAY_1:
	 wlzData = (WlzGreyPtr *) AlcMalloc(sizeof(WlzGreyPtr) * numChannels)
      break;

    default:
      /* something gone wrong here */
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
    AlcFree(rawData);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}
