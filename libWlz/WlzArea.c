#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzArea_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzArea.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Computes the area of an object.
* \ingroup	WlzFeatures
*/

#include <Wlz.h>

/*!
* \return	Area of object or a negative value on error.
* \ingroup	WlzFeatures
* \brief	Computes the area of an object.
* \todo		Implement for objects other than 2D domain objects.
* \param	obj			Input object.
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
int WlzArea(
  WlzObject *obj,
  WlzErrorNum *dstErr)
{
  int 			size;
  WlzIntervalWSpace 	iwsp;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == (WlzObject *) NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
    size = -1;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	size = -1;
      }
      break;
    
    case WLZ_EMPTY_OBJ:
      size = 0;
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      size = -1;
      break;

    }
  }

  /* calculate size */
  if( errNum == WLZ_ERR_NONE ){
    size = 0;
    errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC);
    if( errNum == WLZ_ERR_NONE ){
      while( (errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	size += iwsp.colrmn;
      }
      if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return size;
}
