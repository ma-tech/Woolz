#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPoints_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzPoints.c
* \author       Bill Hill
* \date         December 2005
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
* \brief	Functions for handling point domains.
* \ingroup	WlzFeatures
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	New domain object which coresponds to the union of
* 		the given points.
* \ingroup	WlzFeatures
* \brief	Creates a domain object which coresponds to the union of
* 		the given points.
* \param	pnt			Point domain.
* \param	scale			Scale, which if greater than zero
* 					is used as the diameter of a circle
* 					or sphere centred on each of the
* 					points vertices and a multiplier
* 					for the point position.
* \param	dstErr			Destination error poiter, may be NULL.
*/
WlzObject	*WlzPointsToDomObj(WlzPoints *pnt, double scale,
    				   WlzErrorNum *dstErr)
{
  int		idP;
  WlzObjectType	dType;
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
		*tObj2 = NULL,
		*dObj = NULL;
  WlzVertex	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idP = 0;
  if(scale > DBL_EPSILON)
  {
    pos.d3.vtZ = 0.0;
  }
  else
  {
    pos.i3.vtZ = 0;
  }
  if(pnt == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pnt->type)
    {
      case WLZ_POINTS_2I: /* FALLTHROUGH */
      case WLZ_POINTS_2D:
	dType = WLZ_2D_DOMAINOBJ;
	break;
      case WLZ_POINTS_3I: /* FALLTHROUGH */
      case WLZ_POINTS_3D:
	dType = WLZ_3D_DOMAINOBJ;
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj0 = WlzMakeEmpty(&errNum);
  }
  while((errNum == WLZ_ERR_NONE) && (idP < pnt->nPoints))
  {
    if(scale > DBL_EPSILON)
    {
      switch(pnt->type)
      {
	case WLZ_POINTS_2I:
	  pos.d3.vtX = pnt->points.i2[idP].vtX;
	  pos.d3.vtY = pnt->points.i2[idP].vtY;
	  break;
	case WLZ_POINTS_2D:
	  pos.d3.vtX = pnt->points.d2[idP].vtX;
	  pos.d3.vtY = pnt->points.d2[idP].vtY;
	  pos.d3.vtZ = 0.0;
	  break;
	case WLZ_POINTS_3I:
	  pos.d3.vtX = pnt->points.i3[idP].vtX;
	  pos.d3.vtY = pnt->points.i3[idP].vtY;
	  pos.d3.vtZ = pnt->points.i3[idP].vtZ;
	  break;
	case WLZ_POINTS_3D:
	  pos.d3 = pnt->points.d3[idP];
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj1 = WlzMakeSphereObject(dType, scale / 2.0,
				    scale * pos.d3.vtX,
				    scale * pos.d3.vtY,
				    scale * pos.d3.vtZ,
				    &errNum);
      }
    }
    else
    {
      switch(pnt->type)
      {
	case WLZ_POINTS_2I:
	  pos.i3.vtX = pnt->points.i2[idP].vtX;
	  pos.i3.vtY = pnt->points.i2[idP].vtY;
	  break;
	case WLZ_POINTS_2D:
	  pos.i3.vtX = WLZ_NINT(pnt->points.d2[idP].vtX);
	  pos.i3.vtY = WLZ_NINT(pnt->points.d2[idP].vtY);
	  break;
	case WLZ_POINTS_3I:
	  pos.i3 = pnt->points.i3[idP];
	  break;
	case WLZ_POINTS_3D:
	  pos.i3.vtX = WLZ_NINT(pnt->points.d3[idP].vtX);
	  pos.i3.vtY = WLZ_NINT(pnt->points.d3[idP].vtY);
	  pos.i3.vtZ = WLZ_NINT(pnt->points.d3[idP].vtZ);
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj1 = WlzMakeSinglePixelObject(dType,
					 pos.i3.vtX, pos.i3.vtY, pos.i3.vtZ,
				         &errNum);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tObj2 = WlzUnion2(tObj0, tObj1, &errNum);
    }
    (void )WlzFreeObj(tObj0);
    (void )WlzFreeObj(tObj1);
    tObj0 = tObj2;
    ++idP;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dObj = tObj0;
  }
  else
  {
    (void )WlzFreeObj(tObj0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dObj);
}
