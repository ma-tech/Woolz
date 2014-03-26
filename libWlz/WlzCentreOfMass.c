#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCentreOfMass_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCentreOfMass.c
* \author       Bill Hill
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
* \brief	Computes the centre of mass of objects.
* \ingroup	WlzFeatures
*/

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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a Woolz object.
*		If the given object does not have grey values or the
*		binary object flag is set (ie non zero) then every pixel or
*		vertex within the object's domain has the same mass.
* \param	srcObj			Given object.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL if not required.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
WlzDVertex2	WlzCentreOfMass2D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzDVertex2	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass2D FE %p %d %p %p\n",
	   srcObj, binObjFlag, dstMass, dstErr));
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
	errNum = WLZ_ERR_OBJECT_TYPE;
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a Woolz object.
*		If the given object does not have grey values or the binary
*		object flag is set (ie non zero) then every pixel or vertex
*		within the objects domain has the same mass.
* \param	srcObj			Given object.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL if not required.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
WlzDVertex3	WlzCentreOfMass3D(WlzObject *srcObj, int binObjFlag,
				  double *dstMass, WlzErrorNum *dstErr)
{
  double        mass = 0.0;
  WlzDVertex2	cMass2D;
  WlzDVertex3	cMass;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_1),
	  ("WlzCentreOfMass3D FE %p %d %p %p\n",
	   srcObj, binObjFlag, dstMass, dstErr));
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

/*!
* \return	Coordinates of centre of mass.
* \ingroup	WlzFeatures
* \brief	Computes the centre of mass of a vector of 2D vertices.
* \param	nVtx			Number of vertices.
* \param	vtx			The vertices.
*/
WlzDVertex2	WlzCentreOfMassVtx2D(int nVtx, WlzDVertex2 *vtx)
{
  WlzDVertex2	cen;

  WLZ_VTX_2_ZERO(cen);
  if(nVtx > 0)
  {
    int		idx;

    for(idx = 0; idx < nVtx; ++idx)
    {
      WLZ_VTX_2_ADD(cen, cen, vtx[idx]);
    }
    WLZ_VTX_2_SCALE(cen, cen, 1.0 / nVtx);
  }
  return(cen);
}

/*!
* \return	Coordinates of centre of mass.
* \ingroup	WlzFeatures
* \brief	Computes the centre of mass of a vector of 3D vertices.
* \param	nVtx			Number of vertices.
* \param	vtx			The vertices.
*/
WlzDVertex3	WlzCentreOfMassVtx3D(int nVtx, WlzDVertex3 *vtx)
{
  WlzDVertex3	cen;

  WLZ_VTX_3_ZERO(cen);
  if(nVtx > 0)
  {
    int		idx;

    for(idx = 0; idx < nVtx; ++idx)
    {
      WLZ_VTX_3_ADD(cen, cen, vtx[idx]);
    }
    WLZ_VTX_3_SCALE(cen, cen, 1.0 / nVtx);
  }
  return(cen);
}

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a WLZ_2D_DOMAIN_OBJ.
*               If the object has values and the binary object flag is
*               not set then the centre of mass is calculated using
*               the grey level information.
*		\f[
                C_x = \frac{\sum_x{\sum_y{x G(x,y)}}}
		           {\sum_x{\sum_y{G(x,y)}}} ,
                C_y = \frac{\sum_x{\sum_y{y G(x,y)}}}
		           {\sum_x{\sum_y{G(x,y)}}}
		\f]
*               Where \f$(C_x,C_y)\f$ are the coordinates of the centre of
*               mass.
*               If the given object does not have grey values or the
*               binary object flag is set (ie non zero) then every
*               pixel within the objects domain has the same mass.
* \param	srcObj			Given object.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
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
	 (gWsp.pixeltype != WLZ_GREY_DOUBLE) &&
	 (gWsp.pixeltype != WLZ_GREY_RGBA))
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
	    case WLZ_GREY_RGBA:
	      while(iCount-- > 0)
	      {
		tmpD = WLZ_RGBA_MODULUS(*(gPix.rgbp));
		sum.vtY += pos.vtY * tmpD;
		sum.vtX += pos.vtX * tmpD;
		mass += tmpD;
		++(gPix.rgbp);
		++(pos.vtX);
	      }
	      break;
	    default:
	      break;
	  }
	}
	if(errNum == WLZ_ERR_EOO)        /* Reset error from end of intervals */ 
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      (void )WlzEndGreyScan(&gWsp);
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a WLZ_3D_DOMAIN_OBJ.
*		If the object has values and the binary object flag is not set
*		then the centre of mass is calculated using the grey level
*		information.
*		\f[
                C_x = \frac{\sum_x{\sum_y{\sum_z{x G(x,y,z)}}}}
		           {\sum_x{\sum_y{\sum_z{G(x,y,z)}}}} ,
                C_y = \frac{\sum_x{\sum_y{\sum_z{y G(x,y,z)}}}}
		           {\sum_x{\sum_y{\sum_z{G(x,y,z)}}}}
                C_z = \frac{\sum_x{\sum_y{\sum_z{z G(x,y,z)}}}},
		           {\sum_x{\sum_y{\sum_z{G(x,y,z)}}}}
		\f]
*               Where \f$(C_x,C_y,C_z)\f$ are the coordinates of the centre of
*               mass.
*               If the given object does not have grey values or the
*               binary object flag is set (ie non zero) then every
*               pixel within the objects domain has the same mass.
* \param	srcObj			Given object.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
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
  else if((binObjFlag == 0) && (srcObj->values.core != NULL) &&
          WlzGreyTableIsTiled(srcObj->values.core->type))
  {
    errNum = WLZ_ERR_VOXELVALUES_TYPE;
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
    (void )WlzFreeObj(srcObj2D);
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a contour.
* \param	srcObj			Given object.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
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
      case WLZ_GMMOD_2N:
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a contour.
* \param	srcObj			Given object.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
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
      case WLZ_GMMOD_2N:
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
      case WLZ_GMMOD_3N:
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a geometric model.
* \param	model			Given geometric model.
* \param	dstMass			Destination pointer for mass, may be
* 					NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
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
      *dstMass = (double )mass;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cMass);
}

/*!
* \return	Coordinates of center of mass.
* \brief	Calculates the centre of mass of a transformed object.
* \param	srcObj			Given object.
* \param	trans			Given transform.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass,
*					may be NULL.
* \param	dstErr			Destination pointer for error,
*					may be NULL.
*/
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
    cMass = WlzAffineTransformVertexD2(trans, cMass, &errNum);
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

/*!
* \return	Coordinates of center of mass.
* \ingroup	WlzFeatures
* \brief	Calculates the centre of mass of a transformed object.
* \param	srcObj			Given object.
* \param	trans			Given transform.
* \param	binObjFlag		Binary object flag.
* \param	dstMass			Destination pointer for mass,
*					may be NULL.
* \param	dstErr			Destination pointer for error,
*					may be NULL.
*/
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
    cMass = WlzAffineTransformVertexD3(trans, cMass, &errNum);
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
