#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzPrinicipalAngle.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the angle which the long principal axis
*		makes with the  given object.
*		This code is based on code by lizg@hgu.mrc.ac.uk which
*		implemented a combination of methods by Rees and
*		Hibbard. (D.W.A. Rees, Mechanics of Solids and
*		Structures).
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzPrincipalAngle					*
* Returns:	double:			Angle in radians.		*
* Purpose:	Calculates the angle  which the long principal axis	*
*		makes with the x-axis in the given object.		*
*		theta = arctan(Ixy*2/(Iyy-Ixx)) / 2			*
*		where							*
*		Ixx = SUMySUMx{(y-Cy)*(y-Cy)*G(x,y)}			*
*		Iyy = SUMySUMx{(x-Cx)(x-Cx)*G{x,y)}			*
*		Ixy = - SUMySUMx{(x-Cx)*(y-Cy)*G(x,y)}			*
*		and (Cx,Cy) are the coordinates of the centre of	*
*		mass.							*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given object.			*
*		WlzDVertex2 cMass:	Center of mass of given object.	*
*		int binObjFlag:		Given object is binary if non	*
*					zero.				*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number, may be NULL if not 	*
*					required.			*
************************************************************************/
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
      while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
      {
	delta.vtY = iWsp.linpos - cMass.vtY;
	delta.vtX = iWsp.lftpos - cMass.vtX;
	iCount = iWsp.rgtpos - iWsp.lftpos + 1;
	tI0 = (delta.vtX + iCount - 1);
	ixx += iCount * delta.vtY * delta.vtY;
	iyy += ((tI0 * (tI0 + 1) * ((2 * tI0) + 1)) -
	        ((delta.vtX - 1) * delta.vtX * ((2 * delta.vtX) - 1))) / 6.0;
        ixy -= ((tI0 * (tD0 + 1)) - (delta.vtX * (delta.vtX - 1))) *
	       delta.vtY / 2.0;
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
