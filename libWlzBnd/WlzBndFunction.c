#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndFunction_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndFunction.c
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
* \brief	Binding for access to Woolz objects.
* \ingroup	LibWlzBnd
*/
#include <WlzBnd.h>


WlzErrorNum	WlzSetVoxelSize(WlzObject *obj,
				double x, double y, double z)
{
  WlzErrorNum errNum = WLZ_ERR_NONE;

  if(obj->type == WLZ_3D_DOMAINOBJ)
  {
    obj->domain.p->voxel_size[0] = x;
    obj->domain.p->voxel_size[1] = y;
    obj->domain.p->voxel_size[2] = z;
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  return(errNum);
}


WlzObjectType	WlzGetObjectType(WlzObject *obj)
{
  return(obj->type);
}

WlzErrorNum	WlzExplode(int *dstExpObjCount,
			   WlzObject ***dstExpObjVecP,
		           WlzObject *srcObj)
{
  WlzErrorNum errNum = WLZ_ERR_NONE;
  WlzCompoundArray *cmpObj = NULL;

  if((srcObj->type == WLZ_COMPOUND_ARR_1) ||
     (srcObj->type == WLZ_COMPOUND_ARR_2))
  {
    cmpObj = (WlzCompoundArray *)srcObj;
    *dstExpObjVecP = cmpObj->o;
    *dstExpObjCount = cmpObj->n;
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  return(errNum);
}

const char	*WlzGetPropName(WlzObject *obj)
{
  WlzErrorNum errNum = WLZ_ERR_NONE;
  WlzProperty prop;
  const char *name = NULL;
  char *dst = NULL;

  if(obj->plist)
  {
    prop = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_NAME, &errNum);
    if(prop.core && (errNum == WLZ_ERR_NONE))
    {
      switch(prop.core->type)
      {
	case WLZ_PROPERTY_NAME:
	  dst = prop.name->name;
	  break;
	case WLZ_PROPERTY_GREY:
	  dst = prop.greyV->name;
	  break;
	default:
	  errNum = WLZ_ERR_PROPERTY_TYPE;
	  break;
      }
    }
  }
  name = AlcStrDup((const char *) dst);
  return(name);
}

WlzObject	*WlzGetContourObj(WlzObject *inObj){
  int		flip = 1,
		nrm = 0,
		nItr = 10,
		nonMan = 0,
		filterGeom = 1,
		setbackVz = 0;
  double 	lambda = 0,
		mu = 0,
		ctrVal = 100,
		ctrWth = 1.0,
		filterPB = 0.1,
		filterSB = 1.1,
		xsize = 1,
		ysize = 1,
		zsize = 1;
  WlzDomain	ctrDom;
  WlzValues	dumVal;
  WlzContourMethod ctrMtd = WLZ_CONTOUR_MTD_BND;
  const double	filterDPB = 0.25,
		filterDSB = 0.10;
  const char	*errMsgStr;
  WlzObject	*outObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  ctrDom.core = NULL;
  dumVal.core = NULL;
  if(inObj && (inObj->type == WLZ_3D_DOMAINOBJ) && (inObj->domain.core) &&
     (inObj->domain.core->type == WLZ_2D_DOMAINOBJ))
  {
    xsize = inObj->domain.p->voxel_size[0];
    ysize = inObj->domain.p->voxel_size[1];
    zsize = inObj->domain.p->voxel_size[2];
    inObj->domain.p->voxel_size[0] = 1.0;
    inObj->domain.p->voxel_size[1] = 1.0;
    inObj->domain.p->voxel_size[2] = 1.0;
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ctrDom.ctr = WlzContourObj(inObj, ctrMtd,ctrVal, ctrWth, nrm, &errNum);
    if((errNum != WLZ_ERR_NONE) && filterGeom)
    {
      errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &nItr, 
	                              filterPB, filterSB, filterDPB,
				      filterDSB);
      if(errNum != WLZ_ERR_NONE)
      {
	errNum = WlzGMFilterGeomLPLM(ctrDom.ctr->model,
	                             lambda, mu, nItr, nonMan);
	if(errNum != WLZ_ERR_NONE)
	{
	  outObj = WlzMakeMain(WLZ_CONTOUR, ctrDom, 
	                       dumVal, NULL, NULL, &errNum);
	  if(setbackVz)
	  {
	    inObj->domain.p->voxel_size[0] = xsize;
	    inObj->domain.p->voxel_size[1] = ysize;
	    inObj->domain.p->voxel_size[2] = zsize;
	  }
	  if((errNum != WLZ_ERR_NONE) &&
	     flip && ctrDom.core && ctrDom.ctr->model)
	  {
	    errNum  = WlzGMFilterFlipOrient(ctrDom.ctr->model);
	  }
	}
      }
    }
  }
  (void)WlzStringFromErrorNum(errNum, &errMsgStr);
  return(outObj);
}

int		WlzDestroyObj(WlzObject *obj)
{
  int success = 0;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  errNum = WlzFreeObj(obj) ;
  if(errNum == WLZ_ERR_NONE)
  {
    success = 1;
  }
  return(success);
}
