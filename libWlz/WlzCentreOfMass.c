#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzCentreOfMass.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes the centre of mass of Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzCentreOfMass2D					*
* Returns:	WlzDvertex:		Coordinates of center of mass.	*
* Purpose:	Calculates the centre of mass of a WLZ_2D_DOMAIN_OBJ.	*
*		If the object has values and the binary object flag is	*
*		not set then the centre of mass is calculated using	*
*		the grey level information.				*
*		Cx = SUMxSUMy{y*G(x,y)} / SUMxSUMy{G(x,y)}		*
*		Cy = SUMxSUMy{x*G(x,y)} / SUMxSUMy{G(x,y)}		*
*		Where (Cx,Cy) are the coordinates of the centre of	*
*		mass.							*
*		If the given object does not have grey values or the	*
*		binary object flag is set (ie non zero) then every	*
*		pixel within the objects domain has the same mass.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given object.			*
*		int binObjFlag:		Binary object flag.		*
*		double *dstMass:	Destination pointer for mass,	*
*					may be NULL if not required.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number, may be NULL if not 	*
*					required.			*
************************************************************************/
WlzDVertex2	WlzCentreOfMass2D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErrNum)
{
  int		iCount;
  double        mass = 0.0,
		tmpD;
  WlzIntervalWSpace iWsp;
  WlzGreyWSpace	gWsp;
  WlzGreyP      gPix;
  WlzIVertex2	pos;
  WlzDVertex2	cMass,
  		sum;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass2D FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj, binObjFlag,
	   (unsigned long )dstMass, (unsigned long )dstErrNum));
  sum.vtX = 0.0;
  sum.vtY = 0.0;
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
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
	iCount = iWsp.rgtpos - iWsp.lftpos + 1;
	mass += iCount;
	sum.vtX += ((iWsp.rgtpos * (iWsp.rgtpos + 1)) -
		    (iWsp.lftpos * (iWsp.lftpos - 1))) / 2.0;
	sum.vtY += iWsp.linpos * iCount;
      }
      if(errNum == WLZ_ERR_EOO)	        /* Reset error from end of intervals */ 
      {
	errNum = WLZ_ERR_NONE;
      }
    }
    else
    {
      if((gWsp.pixeltype != WLZ_GREY_INT) &&
	 (gWsp.pixeltype != WLZ_GREY_SHORT) &&
	 (gWsp.pixeltype != WLZ_GREY_UBYTE) &&
	 (gWsp.pixeltype != WLZ_GREY_FLOAT) &&
	 (gWsp.pixeltype != WLZ_GREY_DOUBLE))
      {
        errNum = WLZ_ERR_GREY_TYPE;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
	{
	  pos.vtX = iWsp.lftpos;
	  pos.vtY = iWsp.linpos;
	  gPix = gWsp.u_grintptr;
	  iCount = iWsp.rgtpos - iWsp.lftpos + 1;
	  switch(gWsp.pixeltype)
	  {
	    case WLZ_GREY_INT:
	      while(iCount-- > 0)
	      {
		tmpD = *(gPix.inp);
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.inp);
		++(pos.vtX);
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      while(iCount-- > 0)
	      {
		tmpD = *(gPix.shp);
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.shp);
		++(pos.vtX);
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      while(iCount-- > 0)
	      {
		tmpD = *(gPix.ubp);
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.ubp);
		++(pos.vtX);
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      while(iCount-- > 0)
	      {
		tmpD = *(gPix.flp);
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.flp);
		++(pos.vtX);
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      while(iCount-- > 0)
	      {
		tmpD = *(gPix.dbp);
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.dbp);
		++(pos.vtX);
	      }
	      break;
	  }
	}
	if(errNum == WLZ_ERR_EOO)        /* Reset error from end of intervals */ 
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((mass > DBL_EPSILON) || (mass < (-(DBL_EPSILON))))
    {
      cMass.vtX = sum.vtX / mass;
      cMass.vtY = sum.vtY / mass;
    }
    if(dstMass)
    {
      *dstMass = mass;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzCentreOfMass2D 01 %d\n",
	   (int )errNum));
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass2D FX {%g %g}\n",
	   cMass.vtX, cMass.vtY));
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMass3D					*
* Returns:	WlzDvertex:		Coordinates of center of mass.	*
* Purpose:	Calculates the centre of mass of a WLZ_2D_DOMAIN_OBJ	*
*		or a WLZ_3D_DOMAIN_OBJ.					*
*		If the object has values and the binary object flag is	*
*		not set then the centre of mass is calculated using	*
*		the grey level information.				*
*		Cx = SUMxSUMySUMz{y*G(x,y,z)} / SUMxSUMySUMz{G(x,y,z)}	*
*		Cy = SUMxSUMySUMz{x*G(x,y,z)} / SUMxSUMySUMz{G(x,y,z)}	*
*		Cz = SUMxSUMySUMz{z*G(x,y,z)} / SUMxSUMySUMz{G(x,y,z)}	*
*		Where (Cx,Cy,Cz) are the coordinates of the centre of	*
*		mass.							*
*		If the given object does not have grey values or the	*
*		binary object flag is set (ie non zero) then every	*
*		pixel within the objects domain has the same mass.	*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Given object.			*
*		int binObjFlag:		Binary object flag.		*
*		double *dstMass:	Destination pointer for mass,	*
*					may be NULL if not required.	*
*		WlzErrorNum *dstErrNum:	Destination pointer for error	*
*					number, may be NULL if not 	*
*					required.			*
************************************************************************/
WlzDVertex3	WlzCentreOfMass3D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErrNum)
{
  int		planeIdx,
  		planeCount;
  double        mass = 0.0,
  		mass2D;
  WlzDVertex2	cMass2D;
  WlzDVertex3	cMass,
  		sum;
  WlzDomain	srcDom,
  		dummyDom;
  WlzDomain	*srcDomains;
  WlzValues	dummyValues;
  WlzValues	*srcValues;
  WlzObject	*srcObj2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass3D FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj, binObjFlag,
	   (unsigned long )dstMass, (unsigned long )dstErrNum));
  sum.vtX = 0.0;
  sum.vtY = 0.0;
  sum.vtZ = 0.0;
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
  cMass.vtZ = 0.0;
  dummyDom.core = NULL;
  dummyValues.core = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	cMass2D = WlzCentreOfMass2D(srcObj, binObjFlag, dstMass, &errNum);
	cMass.vtX = cMass2D.vtX;
	cMass.vtY = cMass2D.vtY;
	break;
      case WLZ_3D_DOMAINOBJ:
	if((srcDom = srcObj->domain).core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcDom.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if((srcDomains = srcDom.p->domains) == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  if((srcValues = srcObj->values.vox->values) == NULL)
	  {
	    binObjFlag = 1;
	  }
	  srcObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, dummyDom,
				 dummyValues, NULL, NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  planeIdx =  0;
	  planeCount = srcDom.p->lastpl - srcDom.p->plane1 + 1;
	  while((errNum == WLZ_ERR_NONE) && (planeCount-- > 0))
	  {
	    srcObj2D->domain = *(srcDomains + planeIdx);
	    if(binObjFlag)
	    {
	      srcObj2D->values.core = NULL;
	    }
	    else
	    {
	      srcObj2D->values = *(srcValues + planeIdx);
	    }
	    cMass2D = WlzCentreOfMass2D(srcObj2D, binObjFlag, &mass2D,
		&errNum);
	    sum.vtX += cMass2D.vtX * mass2D;
	    sum.vtY += cMass2D.vtY * mass2D;
	    sum.vtZ += (srcDom.p->plane1 + planeIdx) * mass2D;
	    mass += mass2D;
	    ++planeIdx;
	  }
	  srcObj2D->domain.core = NULL;
	  srcObj2D->values.core = NULL;
	  WlzFreeObj(srcObj2D);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if((mass > DBL_EPSILON) || (mass < (-(DBL_EPSILON))))
	    {
	      cMass.vtX = sum.vtX / mass;
	      cMass.vtY = sum.vtY / mass;
	      cMass.vtZ = sum.vtZ / mass;
	    }
	    if(dstMass)
	    {
	      *dstMass = mass;
	    }
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzCentreOfMass3D 01 %d\n",
	   (int )errNum));
  if(dstErrNum)
  {
    *dstErrNum = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass3D FX {%g %g %g}\n",
	   cMass.vtX, cMass.vtY, cMass.vtZ));
  return(cMass);
}
