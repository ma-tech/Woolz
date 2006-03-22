#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzPrinicipalAngle_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzPrinicipalAngle.c
* \author       Elizabeth Guest, Bill Hill
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
* \brief	Calculates the angle which the long principal axis makes with
* 		the  given object. The functions are based on a combination of
* 		methods by Rees and Hibbard. (D.W.A. Rees, Mechanics of Solids
* 		and Structures).
* \ingroup	WlzRegistration
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	Angle in radians.
* \ingroup	WlzRegistration
* \brief	Calculates the angle  which the long principal axis
*               makes with the x-axis in the given object.
*               \f[
		  \theta = \frac{1}{2}
		  	   \arctan{(2 \frac{I_{xy}}{I_{yy}-I_{xx}})}
		\f]
*               where
*               \f[
		  I_{xx} = \sum_y \sum_x ((y - C_y) (y - C_y) G(x, y))
		\f]
*		\f[
                  I_{yy} = \sum_y \sum_x ((x - C_x) (x - C_x) G(x, y))
		\f]
*		\f[
                  I_{xy} = - \sum_y \sum_x (((x - C x) (y - C_y) G(x, y))
		\f]
*               and \f$(C_x, C_y)\f$ are the coordinates of the centre of
*               mass.
* \param	srcObj			Given object.
* \param	cMass			Center of mass of given object.
* \param	binObjFlag		Given object is binary if non
*                                       zero.
* \param	dstErrNum		Destination pointer for error
*                                       number, may be NULL if not
*                                       required.
*/
double		WlzPrincipalAngle(WlzObject *srcObj, WlzDVertex2 cMass,
				  int binObjFlag, WlzErrorNum *dstErrNum)
{
  int 		tI0,
		tI1,
  		iCount;
  double 	pAngle,
		tD0,
		tD1,
		root0,
		root1,
		ixx = 0.0,
		iyy = 0.0,
		ixy = 0.0;
  WlzIVertex2	delta;
  WlzIntervalWSpace iWsp;
  WlzGreyWSpace	gWsp;
  WlzGreyP	gPix;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzPrincipalAngle FE 0x%lx {%d %d} %d\n",
	   (unsigned long )srcObj, cMass.vtX, cMass.vtY, binObjFlag));
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
  else if((srcObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
          (srcObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if((srcObj->values.core == NULL) ||
       (srcObj->values.core->type == WLZ_EMPTY_OBJ))
    {
      binObjFlag = 1;
    }
    if(binObjFlag)
    {
      errNum = WlzInitRasterScan(srcObj, &iWsp, WLZ_RASTERDIR_ILIC);
    }
    else
    {
      errNum = WlzInitGreyScan(srcObj, &iWsp, &gWsp);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(binObjFlag)
    {
      while((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE)
      {
	delta.vtY = iWsp.linpos - cMass.vtY;
	delta.vtX = iWsp.lftpos - cMass.vtX;
	iCount = iWsp.rgtpos - iWsp.lftpos + 1;
	tI0 = iCount * (iCount - 1);
	ixx += iCount * delta.vtY * delta.vtY;
	iyy += (iCount * delta.vtX * delta.vtX) + (tI0 * delta.vtX) +
	       (tI0 * ((2 * iCount) - 1) / 6);
	ixy -= (iCount * delta.vtX * delta.vtY) + (tI0 * delta.vtY / 2);
      }
      if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
      {
        errNum = WLZ_ERR_NONE;
      }
    }
    else
    {
      while(WlzNextGreyInterval(&iWsp) == 0)
      {
	gPix = gWsp.u_grintptr;
	delta.vtX = iWsp.lftpos - cMass.vtX;
	delta.vtY = iWsp.linpos - cMass.vtY;
	iCount = iWsp.rgtpos - iWsp.lftpos;
	switch(gWsp.pixeltype)
	{
	  case WLZ_GREY_INT:
	    while(iCount-- >= 0)
	    {
	      tI0 = *gPix.inp;
	      tI1 = delta.vtX * tI0;
	      ixx += delta.vtY * delta.vtY * tI0;
	      iyy += delta.vtX * tI1;
	      ixy -= delta.vtY * tI1;
	      ++gPix.inp;
	      ++delta.vtX;
	    }
	    break;
	  case WLZ_GREY_SHORT:
	    while(iCount-- >= 0)
	    {
	      tI0 = *gPix.shp;
	      tI1 = delta.vtX * tI0;
	      ixx += delta.vtY * delta.vtY * tI0;
	      iyy += delta.vtX * tI1;
	      ixy -= delta.vtY * tI1;
	      ++gPix.inp;
	      ++delta.vtX;
	    }
	    break;
	  case WLZ_GREY_UBYTE:
	    while(iCount-- >= 0)
	    {
	      tI0 = *gPix.ubp;
	      tI1 = delta.vtX * tI0;
	      ixx += delta.vtY * delta.vtY * tI0;
	      iyy += delta.vtX * tI1;
	      ixy -= delta.vtY * tI1;
	      ++gPix.ubp;
	      ++delta.vtX;
	    }
	    break;
	  case WLZ_GREY_FLOAT:
	    while(iCount-- >= 0)
	    {
	      tD0 = *gPix.flp;
	      tD1 = delta.vtX * tD0;
	      ixx += delta.vtY * delta.vtY * tD0;
	      iyy += delta.vtX * tD1;
	      ixy -= delta.vtY * tD1;
	      ++gPix.flp;
	      ++delta.vtX;
	    }
	    break;
	  case WLZ_GREY_DOUBLE:
	    while(iCount-- >= 0)
	    {
	      tD0 = *gPix.dbp;
	      tD1 = delta.vtX * tD0;
	      ixx += delta.vtY * delta.vtY * tD0;
	      iyy += delta.vtX * tD1;
	      ixy -= delta.vtY * tD1;
	      ++gPix.dbp;
	      ++delta.vtX;
	    }
	    break;
	  case WLZ_GREY_RGBA:
	    while(iCount-- >= 0)
	    {
	      tI0 = WLZ_RGBA_MODULUS(*gPix.rgbp);
	      tI1 = delta.vtX * tI0;
	      ixx += delta.vtY * delta.vtY * tI0;
	      iyy += delta.vtX * tI1;
	      ixy -= delta.vtY * tI1;
	      ++gPix.rgbp;
	      ++delta.vtX;
	    }
	    break;
	}
      }
      if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
      {
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tD0 = ixx + iyy;
    tD1 = sqrt((tD0 * tD0) - (4.0 * ((ixx * iyy) - (ixy * ixy))));
    root0 = (tD0 - tD1) / 2.0;
    root1 = (tD0 + tD1) / 2.0;
    if(WLZ_ABS(root0) > WLZ_ABS(root1))
    {
      root0 = root1;
    }
    if(WLZ_ABS(root0 - ixx) < DBL_EPSILON)
    {
      pAngle = 0.0;
    }
    else if(WLZ_ABS(ixy) < DBL_EPSILON)
    {
      pAngle = WLZ_M_PI_2;
    }
    else
    {
      pAngle = atan((root0 - ixx) / ixy);
    }
  }
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzPrincipalAngle FX %d\n",
	   pAngle));
  return(pAngle);
}
