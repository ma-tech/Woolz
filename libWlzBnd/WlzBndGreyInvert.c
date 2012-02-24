#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndGreyInvert_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndGreyInvert.c
* \author       Guangjie Feng
* \date         August 2003
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
* \brief	Binding for Woolz grey inversion.
* \ingroup	LibWlzBnd
*/

#include <WlzBnd.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzBnd
* \brief	A binding for the WlzGreyRange() and WlzGreyInvertMinMax()
* 		which inverts an objects grey value while keeping the same
*		grey range as in the Woolz binary.
* \param	obj			Given object.
*/
WlzErrorNum	WlzBndGreyInvert(WlzObject *obj)
{
  WlzPixelV	min, max, gmin, gmax;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((errNum = WlzGreyRange(obj, &gmin, &gmax)) == WLZ_ERR_NONE)
  {
    WlzValueConvertPixel(&min, gmin, WLZ_GREY_DOUBLE);
    WlzValueConvertPixel(&max, gmax, WLZ_GREY_DOUBLE);

    errNum = WlzGreyInvertMinMax(obj, min, max);
  }
  return(errNum);
}


