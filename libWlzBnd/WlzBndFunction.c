#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndFunction_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndFunction.c
* \author       Guangjie Feng, Bill Hill
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


/*!
* \ingroup	LibWlzBnd
* \brief	Simple wrapper to get version.
* \param	dstStr			Destination pointer for the version
* 					string.
*/
void		WlzGetVersion(const char **dstStr)
{
  if(dstStr)
  {
    *dstStr = WlzVersion();
  }
}

/*!
* \return	Voxel size.
* \ingroup	LibWlzBnd
* \brief	Gets the object's voxelsize. It is and error if the
* 		object is NULL, the object is not a 3D spatial domain
* 		object or the object's domain is NULL.
* \param	obj			Given object.
* \param	x			Voxel x size.
* \param	y			Voxel y size.
* \param	z			Voxel z size.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDVertex3	WlzGetVoxelSize(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzDVertex3	voxSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    voxSz.vtX = obj->domain.p->voxel_size[0];
    voxSz.vtY = obj->domain.p->voxel_size[1];
    voxSz.vtZ = obj->domain.p->voxel_size[2];
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(voxSz);
}

/*!
* \return	Woolz error code.
* \ingroup	LibWlzBnd
* \brief	Sets the object's voxelsize. It is and error if the
* 		object is NULL, the object is not a 3D spatial domain
* 		object or the object's domain is NULL.
* \param	obj			Given object.
* \param	x			Voxel x size.
* \param	y			Voxel y size.
* \param	z			Voxel z size.
*/
WlzErrorNum	WlzSetVoxelSize(WlzObject *obj,
				double x, double y, double z)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    obj->domain.p->voxel_size[0] = x;
    obj->domain.p->voxel_size[1] = y;
    obj->domain.p->voxel_size[2] = z;
  }
  return(errNum);
}


/*!
* \return	Non zero if the given object is NULL.
* \ingroup	LibWlzBnd
* \brief	Tests if an object is NULL.
* \param	obj			Given object.
*/
int		WlzObjectIsNull(WlzObject *obj)
{
  int		t;

  t = (obj == NULL);
  return(t);
}

/*!
* \return	Non zero if the given object's domain is NULL.
* \ingroup	LibWlzBnd
* \brief	Tests if an object's domain is NULL.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzObjectDomainIsNull(WlzObject *obj, WlzErrorNum *dstErr)
{
  int		t = 1;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    t = (obj->domain.core == NULL);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(t);
}

/*!
* \return	Non zero if the given object's values is NULL.
* \ingroup	LibWlzBnd
* \brief	Tests if an object's values is NULL.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzObjectValuesIsNull(WlzObject *obj, WlzErrorNum *dstErr)
{
  int		t = 1;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;


  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    t = (obj->values.core == NULL);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(t);
}

/*!
* \return	Woolz object type.
* \ingroup	LibWlzBnd
* \brief	Get's the given obejct's type.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzGetObjectType(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzObjectType	t = WLZ_NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    t = obj->type;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(t);
}

/*!
* \return	Woolz object's domain type.
* \ingroup	LibWlzBnd
* \brief	Get's the given obejct's domain type.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzGetObjectDomainType(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzObjectType t = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    t = obj->domain.core->type;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(t);
}

/*!
* \return	Woolz object's values type.
* \ingroup	LibWlzBnd
* \brief	Get's the given obejct's values type.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObjectType	WlzGetObjectValuesType(WlzObject *obj, WlzErrorNum *dstErr)
{
  WlzObjectType t = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    t = obj->values.core->type;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(t);
}

/*!
* \return	Woolz error code.
* \ingroup	LibWlzBnd
* \brief	Explodes the given compound array object into an array.
* \param	dstExpObjCount		Destination pointer for number of
* 					objects in array, must not be NULL.
* \param	dstExpObjVecP		Destination pointer for the object
* 					array, must not be NULL.
* \param	obj			Given object.
*/
WlzErrorNum	WlzExplode(int *dstExpObjCount,
			   WlzObject ***dstExpObjVecP,
		           WlzObject *obj)
{
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((obj->type == WLZ_COMPOUND_ARR_1) ||
          (obj->type == WLZ_COMPOUND_ARR_2))
  {
    WlzCompoundArray *cmpObj;

    cmpObj = (WlzCompoundArray *)obj;
    *dstExpObjVecP = cmpObj->o;
    *dstExpObjCount = cmpObj->n;
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  return(errNum);
}

/*!
* \return	Property name string.
* \ingroup	LibWlzBnd
* \brief	Get's an object's property name string.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
const char	*WlzGetPropName(WlzObject *obj, WlzErrorNum *dstErr)
{
  const char *name = NULL;
  char *dst = NULL;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->plist)
  {
    WlzProperty prop;

    prop = WlzGetProperty(obj->plist->list, WLZ_PROPERTY_NAME, &errNum);
    if((errNum == WLZ_ERR_NONE) && prop.core)
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
      name = AlcStrDup((const char *)dst);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(name);
}

/*!
* \return	Contour object.
* \ingroup	LibWlzBnd
* \brief	Extracts a smoothed iso-surface contour with iso-value 1.0
* 		from the given 3D domain object with values.
* \param	inObj			Given domain object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzGetContourObj(WlzObject *inObj, WlzErrorNum *dstErr)
{
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(outObj);
}

/*!
* \return	Non zero on success.
* \ingroup	LibWlzBnd
* \brief	Frees a Woolz object when called.
* \param	obj			Given object.
*/
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

/*!
* \return	New Woolz object or NULL on error.
* \ingroup	LibWlzBnd
* \brief	Wrapper avoiding pixel type with WlzThreshold().
* \param	obj			Object to be thresholded.
* \param	thr			Int threshold pixel value.
* \param	hilo			Mode parameter with possible values:
*					<ul>
*					<li> WLZ_THRESH_LOW - thresholded
*					object is of values < given value.
*					</li>
*					<li> WLZ_THRESH_HIGH - thresholded
*					object is of values >= given value.
*					</li>
*					<li> WLZ_THRESH_EQUAL - thresholded
*					object is of values == given value.
*					</li>
*					</ul>
* \param	dstErr			Destination pointer for error number,
*					may be NULL.
*/
WlzObject		*WlzThresholdI(WlzObject *obj, int thr,
				       WlzThresholdType hilo,
			               WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzPixelV	thrV;

  thrV.type = WLZ_GREY_INT;
  thrV.v.inv = thr;
  tObj = WlzThreshold(obj, thrV, hilo, dstErr);
  return(tObj);
}

/*! 
* \return	New 2D domain object with values or NULL on error.
* \ingroup      LibWlzBnd
* \brief        Wrapper avoiding pixel type for WlzMakeRect() which
* 		creates a 2D domain object with values, for which the
* 		domain is a rectangle and the values are of the given
* 		type and initialized to have given integer background value.
* 		See also WlzMakeRectD().
* \param        line1                   First line.
* \param        lastln                  Last line.
* \param        kol1                    First column.
* \param        lastkl                  Last column.
* \param        pixType                 Pixel type for the grey values. If
*                                       WLZ_GREY_ERROR is given then no values
*                                       are created.
* \param	val			Grey value pointer with allocated
* 					grey values.
* \param        bgd                     Integer background value.
* \param        plist                   Property list to be attached.
* \param        assocObj                Associated object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject 			*WlzMakeRectI(
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixeltype,
				  void *val,
				  int bgd,
				  WlzPropertyList *plist,
				  WlzObject *assoc_obj,
				  WlzErrorNum *dstErr)
{
  WlzObject     *obj = NULL;
  WlzPixelV     bgdV;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = bgd;
  obj = WlzMakeRect(line1, lastln, kol1, lastkl, pixeltype, val, bgdV,
                    plist, assoc_obj, dstErr);
  return(obj);
}

/*! 
* \return	New 2D domain object with values or NULL on error.
* \ingroup      LibWlzBnd
* \brief        Wrapper avoiding pixel type for WlzMakeRect() which
* 		creates a 2D domain object with values, for which the
* 		domain is a rectangle and the values are of the given
* 		type and initialized to have given double background value.
* 		See also WlzMakeRectI().
* \param        line1                   First line.
* \param        lastln                  Last line.
* \param        kol1                    First column.
* \param        lastkl                  Last column.
* \param        pixType                 Pixel type for the grey values. If
*                                       WLZ_GREY_ERROR is given then no values
*                                       are created.
* \param	val			Grey value pointer with allocated
* 					grey values.
* \param        bgd                     Double background value.
* \param        plist                   Property list to be attached.
* \param        assocObj                Associated object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject 			*WlzMakeRectD(
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixeltype,
				  void *val,
				  double bgd,
				  WlzPropertyList *plist,
				  WlzObject *assoc_obj,
				  WlzErrorNum *dstErr)
{
  WlzObject     *obj = NULL;
  WlzPixelV     bgdV;

  bgdV.type = WLZ_GREY_DOUBLE;
  bgdV.v.dbv = bgd;
  obj = WlzMakeRect(line1, lastln, kol1, lastkl, pixeltype, val, bgdV,
                    plist, assoc_obj, dstErr);
  return(obj);
}

/*!
* \return       New 3D domain object with values or NULL on error.
* \ingroup      LibWlzBnd
* \brief        Wrapper avoiding pixel type for WlzMakeCuboid() which
* 		creates a 3D domain object with values, for which the
*               domain is a cuboid and the values are of the given
*               type and initialized to have given integer background value.
*               See also WlzMakeCuboidD().
* \param        plane1                  First plane.
* \param        lastpl                  Last plane.
* \param        line1                   First line.
* \param        lastln                  Last line.
* \param        kol1                    First column.
* \param        lastkl                  Last column.
* \param        pixType                 Pixel type for the grey values. If
*                                       WLZ_GREY_ERROR is given then no values
*                                       are created.
* \param        bgd                     Integer background value.
* \param        plist                   Property list to be attached.
* \param        assocObj                Associated object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject			*WlzMakeCuboidI(
				  int plane1,
				  int lastpl,
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixType,
				  int bgd,
				  WlzPropertyList *plist,
				  WlzObject *assocObj,
				  WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzPixelV	bgdV;

  bgdV.type = WLZ_GREY_INT;
  bgdV.v.inv = bgd;
  obj = WlzMakeCuboid(plane1, lastpl, line1, lastln, kol1, lastkl,
                      pixType, bgdV, plist, assocObj, dstErr);
  return(obj);
}

/*!
* \return       New 3D domain object with values or NULL on error.
* \ingroup      LibWlzBnd
* \brief        Wrapper avoiding pixel type for WlzMakeCuboid() which
* 		creates a 3D domain object with values, for which the
*               domain is a cuboid and the values are of the given
*               type and initialized to have given double background value.
*               See also WlzMakeCuboidD().
* \param        plane1                  First plane.
* \param        lastpl                  Last plane.
* \param        line1                   First line.
* \param        lastln                  Last line.
* \param        kol1                    First column.
* \param        lastkl                  Last column.
* \param        pixType                 Pixel type for the grey values. If
*                                       WLZ_GREY_ERROR is given then no values
*                                       are created.
* \param        bgd                     Double background value.
* \param        plist                   Property list to be attached.
* \param        assocObj                Associated object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject			*WlzMakeCuboidD(
				  int plane1,
				  int lastpl,
				  int line1,
				  int lastln,
				  int kol1,
				  int lastkl,
				  WlzGreyType pixType,
				  double bgd,
				  WlzPropertyList *plist,
				  WlzObject *assocObj,
				  WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzPixelV	bgdV;

  bgdV.type = WLZ_GREY_DOUBLE;
  bgdV.v.dbv = bgd;
  obj = WlzMakeCuboid(plane1, lastpl, line1, lastln, kol1, lastkl,
                      pixType, bgdV, plist, assocObj, dstErr);
  return(obj);
}
