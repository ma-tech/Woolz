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
* 29-08-00 bill	Add/rearange code to include contours, geometric
*		models and transform objects.
************************************************************************/
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static WlzDVertex2 		WlzCentreOfMassDom2D(
				  WlzObject *srcObj,
				  int binObjFlag,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex3 		WlzCentreOfMassDom3D(
				  WlzObject *srcObj,
				  int binObjFlag,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex2 		WlzCentreOfMassCtr2D(
				  WlzObject *srcObj,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex3 		WlzCentreOfMassCtr3D(
				  WlzObject *srcObj,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex3 		WlzCentreOfMassGM(
				  WlzGMModel *model,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex2 		WlzCentreOfMassTrans2D(
				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  int binObjFlag,
				  double *dstMass,
				  WlzErrorNum *dstErr);
static WlzDVertex3 		WlzCentreOfMassTrans3D(
				  WlzObject *srcObj,
				  WlzAffineTransform *trans,
				  int binObjFlag,
				  double *dstMass,
				  WlzErrorNum *dstErr);

/************************************************************************
* Function:	WlzCentreOfMass2D
* Returns:	WlzDvertex2:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a Woolz object.
*		If the given object does not have grey values or the
*		binary object flag is set (ie non zero) then every
*		pixel or vertex within the object's domain has the same
*		mass.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
WlzDVertex2	WlzCentreOfMass2D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzDVertex2	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass2D FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj, binObjFlag,
	   (unsigned long )dstMass, (unsigned long )dstErr));
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_TRANS_OBJ:
	cMass = WlzCentreOfMassTrans2D(srcObj->values.obj,
				       srcObj->domain.t,
				       binObjFlag, &mass, &errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
	cMass = WlzCentreOfMassDom2D(srcObj, binObjFlag, &mass, &errNum);
        break;
      case WLZ_CONTOUR:
	cMass = WlzCentreOfMassCtr2D(srcObj, &mass, &errNum);
	break;
      default:
	WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstMass)
  {
    *dstMass = mass;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzCentreOfMass2D 01 %d\n",
	   (int )errNum));
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass2D FX {%g %g}\n",
	   cMass.vtX, cMass.vtY));
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMass3D
* Returns:	WlzDvertex3:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a Woolz object.
*		If the given object does not have grey values or the
*		binary object flag is set (ie non zero) then every
*		pixel or vertex within the objects domain has the same
*		mass.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
WlzDVertex3	WlzCentreOfMass3D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzDVertex2	cMass2D;
  WlzDVertex3	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass3D FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )srcObj, binObjFlag,
	   (unsigned long )dstMass, (unsigned long )dstErr));
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
  cMass.vtZ = 0.0;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_TRANS_OBJ:
	cMass = WlzCentreOfMassTrans3D(srcObj->values.obj, 
				       srcObj->domain.t,
				       binObjFlag, &mass, &errNum);
        break;
      case WLZ_2D_DOMAINOBJ:
	cMass2D = WlzCentreOfMassDom2D(srcObj, binObjFlag, &mass, &errNum);
	cMass.vtX = cMass2D.vtX;
	cMass.vtY = cMass2D.vtY;
	break;
      case WLZ_3D_DOMAINOBJ:
	cMass = WlzCentreOfMassDom3D(srcObj, binObjFlag, &mass, &errNum);
	break;
      case WLZ_CONTOUR:
	cMass = WlzCentreOfMassCtr3D(srcObj, &mass, &errNum);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && dstMass)
  {
    *dstMass = mass;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_2),
	  ("WlzCentreOfMass3D 01 %d\n",
	   (int )errNum));
  if(dstErr)
  {
    *dstErr = errNum;
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass3D FX {%g %g %g}\n",
	   cMass.vtX, cMass.vtY, cMass.vtZ));
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassDom2D
* Returns:	WlzDvertex2:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a WLZ_2D_DOMAIN_OBJ.
*		If the object has values and the binary object flag is
*		not set then the centre of mass is calculated using
*		the grey level information.
*		Cx = SUMxSUMy{x*G(x,y)} / SUMxSUMy{G(x,y)}
*		Cy = SUMxSUMy{y*G(x,y)} / SUMxSUMy{G(x,y)}
*		Where (Cx,Cy) are the coordinates of the centre of
*		mass.
*		If the given object does not have grey values or the
*		binary object flag is set (ie non zero) then every
*		pixel within the objects domain has the same mass.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex2 WlzCentreOfMassDom2D(WlzObject *srcObj, int binObjFlag,
				        double *dstMass,
					WlzErrorNum *dstErr)
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

  sum.vtX = 0.0;
  sum.vtY = 0.0;
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
  if(srcObj->domain.core == NULL)
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassDom3D
* Returns:	WlzDvertex3:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a WLZ_3D_DOMAIN_OBJ.
*		If the object has values and the binary object flag is
*		not set then the centre of mass is calculated using
*		the grey level information.
*		Cx = SUMxSUMySUMz{y*G(x,y,z)} / SUMxSUMySUMz{G(x,y,z)}
*		Cy = SUMxSUMySUMz{x*G(y,y,z)} / SUMxSUMySUMz{G(x,y,z)}
*		Cz = SUMxSUMySUMz{z*G(x,y,z)} / SUMxSUMySUMz{G(x,y,z)}
*		Where (Cx,Cy,Cz) are the coordinates of the centre of
*		mass.
*		If the given object does not have grey values or the
*		binary object flag is set (ie non zero) then every
*		pixel within the objects domain has the same mass.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex3 WlzCentreOfMassDom3D(WlzObject *srcObj, int binObjFlag,
				        double *dstMass, WlzErrorNum *dstErr)
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

  sum.vtX = 0.0;
  sum.vtY = 0.0;
  sum.vtZ = 0.0;
  cMass.vtX = 0.0;
  cMass.vtY = 0.0;
  cMass.vtZ = 0.0;
  dummyDom.core = NULL;
  dummyValues.core = NULL;
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
    if((srcObj->values.core == NULL) ||
       ((srcValues = srcObj->values.vox->values) == NULL))
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassCtr2D
* Returns:	WlzDvertex2:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a contour.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex2 WlzCentreOfMassCtr2D(WlzObject *srcObj, double *dstMass,
				        WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzContour	*ctr;
  WlzDVertex2	cMass;
  WlzDVertex3	cMass3;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctr = srcObj->domain.ctr) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(ctr->model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(ctr->model->type)
    {
      case WLZ_GMMOD_2I:
      case WLZ_GMMOD_2D:
	cMass3 = WlzCentreOfMassGM(ctr->model, &mass, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  cMass.vtX = cMass3.vtX; cMass.vtY = cMass3.vtY;
	  if(dstMass)
	  {
	    *dstMass = mass;
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassCtr3D
* Returns:	WlzDvertex3:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a contour.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex3 WlzCentreOfMassCtr3D(WlzObject *srcObj, double *dstMass,
				        WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzContour	*ctr;
  WlzDVertex3	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctr = srcObj->domain.ctr) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(ctr->model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(ctr->model->type)
    {
      case WLZ_GMMOD_2I:
      case WLZ_GMMOD_2D:
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
	cMass = WlzCentreOfMassGM(ctr->model, &mass, &errNum);
	if((errNum == WLZ_ERR_NONE) && dstMass)
	{
	  *dstMass = mass;
	}
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassGM
* Returns:	WlzDvertex3:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a geometric model.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex3 WlzCentreOfMassGM(WlzGMModel *model, double *dstMass,
				     WlzErrorNum *dstErr)
{
  int        	iCnt,
  		idI = 0,
  		mass = 0;
  WlzGMVertex	*vertex;
  AlcVector	*vec;
  WlzDVertex3	pos,
  		cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vec = model->res.vertex.vec;
  iCnt = model->res.vertex.numIdx;
  cMass.vtX = cMass.vtY = cMass.vtZ = 0.0;
  while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
  {
    vertex = (WlzGMVertex *)AlcVectorItemGet(vec, idI++);
    if(vertex->idx >= 0)
    {
      ++mass;
      (void )WlzGMVertexGetG3D(vertex, &pos);
      cMass.vtX += pos.vtX; cMass.vtY += pos.vtY; cMass.vtZ += pos.vtZ;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(mass > 0)
    {
      cMass.vtX /= mass; cMass.vtY /= mass; cMass.vtZ /= mass;
    }
    if(dstMass)
    {
      *dstMass = mass;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassTrans2D
* Returns:	WlzDvertex2:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a transformed object.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		WlzAffineTransform *trans: Given transform.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex2 WlzCentreOfMassTrans2D(WlzObject *srcObj,
				     WlzAffineTransform *trans,
				     int binObjFlag,
				     double *dstMass,
				     WlzErrorNum *dstErr)
{
  double	mass;
  WlzDVertex2	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cMass = WlzCentreOfMass2D(srcObj, binObjFlag, &mass, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    cMass = WlzAffineTransformVertexD(trans, cMass, &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && dstMass)
  {
    *dstMass = mass;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/************************************************************************
* Function:	WlzCentreOfMassTrans3D
* Returns:	WlzDvertex3:		Coordinates of center of mass.
* Purpose:	Calculates the centre of mass of a transformed object.
* Global refs:	-
* Parameters:	WlzObject *srcObj:	Given object.
*		WlzAffineTransform *trans: Given transform.
*		int binObjFlag:		Binary object flag.
*		double *dstMass:	Destination pointer for mass,
*					may be NULL if not required.
*		WlzErrorNum *dstErr:	Destination pointer for error
*					number, may be NULL if not
*					required.
************************************************************************/
static WlzDVertex3 WlzCentreOfMassTrans3D(WlzObject *srcObj,
				     WlzAffineTransform *trans,
				     int binObjFlag,
				     double *dstMass,
				     WlzErrorNum *dstErr)
{
  double	mass;
  WlzDVertex3	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cMass = WlzCentreOfMass3D(srcObj, binObjFlag, &mass, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    /* TODO Implement 3D Affine transforms!
     * cMass = WlzAffineTransformVertex3D(trans, cMass, &errNum);
     */
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  if((errNum == WLZ_ERR_NONE) && dstMass)
  {
    *dstMass = mass;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}
