#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DProjection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/Wlz3DProjection.c
* \author       Bill Hill, Richard Baldock
* \date         June 2012
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
* \brief	Generates the projection or back-projection of a
* 		domain object from 3D to 2D and visa versa.
* \ingroup	WlzTransform
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

#define	WLZ_FAST_CODE

#ifdef WLZ_DEBUG_PROJECT3D_TIME
#include <sys/time.h>
#endif /* WLZ_DEBUG_PROJECT3D_TIME */

#ifdef _OPENMP
#include <omp.h>
#endif

static void			WlzProjectObjLine(
				  WlzUByte **ary,
				  WlzIVertex2 p0,
				  WlzIVertex2 p1);
static void			WlzProjectObjLineDom(
				  int **ary,
				  int vpp,
				  WlzIVertex2 p0,
				  WlzIVertex2 p1);
static void			WlzProjectObjLineVal(
				  int **ary,
				  WlzUByte lut[256],
				  WlzUByte *gP,
				  WlzGreyValueWSpace *gVWSp,
				  double **vMat,
				  double vMZX,
				  double vMZY,
				  WlzIVertex3 p0,
				  WlzIVertex3 p1);

/*! 
* \return       projection object
* \ingroup      WlzTransform
* \brief        Use the view transform to define a projection from
*		3D to 2D. Currently only the domain is projected as
*		an opaque shadow.
*		This is old code temporarily kept for compatibility.
* \param    obj	source 3D object
* \param    viewStr	view structure defining the projection
* \param    intFunc	grey-value summation function
* \param    intFuncData data to be passed to the integration function
* \param    dstErr	error return
*/
WlzObject *WlzGetProjectionFromObject(
  WlzObject		*obj,
  WlzThreeDViewStruct 	*viewStr,
  Wlz3DProjectionIntFn 	intFunc,
  void			*intFuncData,
  WlzErrorNum		*dstErr)
{
  WlzObject		*rtnObj=NULL,
  			*obj1;
  WlzThreeDViewStruct	*viewStr1=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzGreyType		srcGType = WLZ_GREY_UBYTE,
  			dstGType = WLZ_GREY_UBYTE;
  WlzPixelV		pixval;
  WlzPixelP		pixptr;
  WlzIntervalWSpace	iwsp;
  WlzGreyWSpace		gwsp;
  WlzGreyValueWSpace	*gVWSp = NULL;
  WlzDVertex3		vtx, vtx1;
  double		x, y, z;
  double		*s_to_x=NULL;
  double		*s_to_y=NULL;
  double		*s_to_z=NULL;
  int			k, xp, yp, s, sp;
  int			length = 0, size = 0, occupiedFlg;
  WlzErrorNum 	errNum=WLZ_ERR_NONE;

  /* check inputs */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if( obj->type != WLZ_3D_DOMAINOBJ ){
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if( obj->domain.core == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if( obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  if( (errNum == WLZ_ERR_NONE) && (viewStr == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* create new view transform */
  if( errNum == WLZ_ERR_NONE ){
    if((viewStr1 = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) != NULL){
      /* need to worry about fixed line mode here sometime */
      viewStr1->fixed = viewStr->fixed;
      viewStr1->theta = viewStr->theta;
      viewStr1->phi = viewStr->phi;
      viewStr1->zeta = viewStr->zeta;
      viewStr1->dist = viewStr->dist;
      viewStr1->scale = viewStr->scale;
      viewStr1->voxelSize[0] = viewStr->voxelSize[0];
      viewStr1->voxelSize[1] = viewStr->voxelSize[1];
      viewStr1->voxelSize[2] = viewStr->voxelSize[2];
      viewStr1->voxelRescaleFlg = viewStr->voxelRescaleFlg;
      viewStr1->interp = viewStr->interp;
      viewStr1->view_mode = viewStr->view_mode;
      viewStr1->up = viewStr->up;

      /* now intialize it */
      /* could optimise by setting fixed point to object centre */
      if( (errNum = WlzInit3DViewStruct(viewStr1, obj)) != WLZ_ERR_NONE ){
	WlzFree3DViewStruct(viewStr1);
	viewStr1 = NULL;
      }
    }
  }

  /* set up orthogonal line parameters & luts */
  if( errNum == WLZ_ERR_NONE ){
    length = WLZ_NINT(viewStr1->maxvals.vtZ) -
      WLZ_NINT(viewStr1->minvals.vtZ) + 1;
    s_to_x = (double *) AlcMalloc(sizeof(double) * length );
    s_to_y = (double *) AlcMalloc(sizeof(double) * length );
    s_to_z = (double *) AlcMalloc(sizeof(double) * length );

    /* transform a perpendicular vector */
    vtx.vtX = 0.0;
    vtx.vtY = 0.0;
    vtx.vtZ = 1.0;
    Wlz3DSectionTransformInvVtx(&vtx, viewStr1);
    vtx1.vtX = 0.0;
    vtx1.vtY = 0.0;
    vtx1.vtZ = 0.0;
    Wlz3DSectionTransformInvVtx(&vtx1, viewStr1);
    vtx.vtX -= vtx1.vtX;
    vtx.vtY -= vtx1.vtY;
    vtx.vtZ -= vtx1.vtZ;

    /* assign lut values */
    s = (int )(WLZ_NINT(viewStr1->minvals.vtZ) - viewStr1->dist);
    for(sp=0; sp < length; sp++, s++){
      s_to_x[sp] = s * vtx.vtX;
      s_to_y[sp] = s * vtx.vtY;
      s_to_z[sp] = s * vtx.vtZ;
    }
  }

  /* if there is an integration function then allocate space for
     the grey-level array */
  if( (errNum == WLZ_ERR_NONE) && (intFunc) ){
    srcGType = WlzGreyTypeFromObj(obj, &errNum);
    switch( srcGType ){
    case WLZ_GREY_LONG:
      size = sizeof(WlzLong)*length;
      break;
    case WLZ_GREY_INT:
      size = sizeof(int)*length;
      break;
    case WLZ_GREY_SHORT:
      size = sizeof(short)*length;
      break;
    case WLZ_GREY_UBYTE:
      size = sizeof(WlzUByte)*length;
      break;
    case WLZ_GREY_FLOAT:
      size = sizeof(float)*length;
      break;
    case WLZ_GREY_DOUBLE:
      size = sizeof(double)*length;
      break;
    case WLZ_GREY_RGBA:
      size = sizeof(int)*length;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
    }
    if( (pixptr.p.inp = (int *) AlcMalloc(size)) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    pixptr.type = srcGType;

    /* set up the grey-value workspace for random access */
    gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
  }

  /* create rectangular projection image */
  if( errNum == WLZ_ERR_NONE ){
    if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
				     WLZ_NINT(viewStr1->minvals.vtY),
				     WLZ_NINT(viewStr1->maxvals.vtY),
				     WLZ_NINT(viewStr1->minvals.vtX),
				     WLZ_NINT(viewStr1->maxvals.vtX),
					 &errNum)) != NULL){
      values.core = NULL;
      if((rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values, NULL, NULL,
			       &errNum)) != NULL){
	/* note the grey-values required are determined by the integration
	   function. Here we use WlzUByte and reset later if needed */
	dstGType = WLZ_GREY_UBYTE;
	pixval.type = WLZ_GREY_UBYTE;
	pixval.v.ubv = (WlzUByte )0;
	if((values.v = WlzNewValueTb(rtnObj,
				     WlzGreyTableType(WLZ_GREY_TAB_RECT,
						      dstGType, NULL),
				     pixval, &errNum)) != NULL){
	  rtnObj->values = WlzAssignValues(values, &errNum);
	}
	else {
	  WlzFreeObj(rtnObj);
	  rtnObj = NULL;
	}
      }
      else {
	WlzFreeDomain(domain);
	domain.core = NULL;
      }
    }
  }

  /* scan image setting values */
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitGreyScan(rtnObj, &iwsp, &gwsp);
  }
  if( errNum == WLZ_ERR_NONE ){
    while( (errNum = WlzNextGreyInterval(&iwsp)) == WLZ_ERR_NONE ){
      yp = iwsp.linpos - WLZ_NINT(viewStr1->minvals.vtY);
      for(k=iwsp.lftpos; k <= iwsp.rgtpos; k++){
	xp = k - WLZ_NINT(viewStr1->minvals.vtX);
	vtx.vtX = viewStr1->xp_to_x[xp] + viewStr1->yp_to_x[yp];
	vtx.vtY = viewStr1->xp_to_y[xp] + viewStr1->yp_to_y[yp];
	vtx.vtZ = viewStr1->xp_to_z[xp] + viewStr1->yp_to_z[yp];

	/* get the projection values */
	/* if no function then just check for occupancy */
	if( intFunc == NULL ){
	  occupiedFlg = 0;
	  sp = (int )(viewStr1->dist - WLZ_NINT(viewStr1->minvals.vtZ));
	  for(; !occupiedFlg && (sp < length); sp++){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    if( WlzInsideDomain(obj, z, y, x, &errNum) ){
	      occupiedFlg = 1;
	    }
	  }
	  sp = (int )(viewStr1->dist - WLZ_NINT(viewStr1->minvals.vtZ) - 1);
	  for(; !occupiedFlg && (sp >= 0); sp--){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    if( WlzInsideDomain(obj, z, y, x, &errNum) ){
	      occupiedFlg = 1;
	    }
	  }

	  /* set the integrated value - only WlzUByte at the moment */
	  *(gwsp.u_grintptr.ubp) = (WlzUByte )occupiedFlg;
	  gwsp.u_grintptr.ubp++;
	}
	/* use integration function */
	else {
	  /* set array of pixel values */
	  for(sp=0; sp < length; sp++){
	    x = vtx.vtX + s_to_x[sp];
	    y = vtx.vtY + s_to_y[sp];
	    z = vtx.vtZ + s_to_z[sp];
	    WlzGreyValueGet(gVWSp, WLZ_NINT(z), WLZ_NINT(y),
			    WLZ_NINT(x));
	    switch( srcGType ){
	    case WLZ_GREY_LONG:
	      pixptr.p.lnp[sp] = gVWSp->gVal[0].lnv;
	      break;
	    case WLZ_GREY_INT:
	      pixptr.p.inp[sp] = gVWSp->gVal[0].inv;
	      break;
	    case WLZ_GREY_SHORT:
	      pixptr.p.shp[sp] = gVWSp->gVal[0].shv;
	      break;
	    case WLZ_GREY_UBYTE:
	      pixptr.p.ubp[sp] = gVWSp->gVal[0].ubv;
	      break;
	    case WLZ_GREY_FLOAT:
	      pixptr.p.flp[sp] = gVWSp->gVal[0].flv;
	      break;
	    case WLZ_GREY_DOUBLE:
	      pixptr.p.dbp[sp] = gVWSp->gVal[0].dbv;
	      break;
	    case WLZ_GREY_RGBA:
	      pixptr.p.rgbp[sp] = gVWSp->gVal[0].rgbv;
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	    }
	  }
	  /* call integration function and seet value */
	  intFunc(pixptr, length, (int )(viewStr1->dist -
		                         WLZ_NINT(viewStr1->minvals.vtZ)),
		  intFuncData, &errNum);
	}
      }
    }
    (void )WlzEndGreyScan(&iwsp, &gwsp);
    if(errNum == WLZ_ERR_EOO)	   /* Reset error from end of intervals */ 
    {
      errNum = WLZ_ERR_NONE;
    }

    /* if no integration function then threshold - binary only */
    if( intFunc == NULL ){
      pixval.v.ubv = 1;
      rtnObj = WlzAssignObject(rtnObj, NULL);
      if((obj1 = WlzThreshold(rtnObj, pixval, WLZ_THRESH_HIGH,
                              &errNum)) != NULL){
	WlzFreeObj(rtnObj);
	values.core = NULL;
	rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, obj1->domain, values, NULL, NULL,
			     &errNum);
	WlzFreeObj(obj1);
      }
      else {
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }

  /* clear space */
  if( viewStr1 ){
    errNum = WlzFree3DViewStruct(viewStr1);
  }
  if( s_to_x ){
    AlcFree( s_to_x );
  }
  if( s_to_y ){
    AlcFree( s_to_y );
  }
  if( s_to_z ){
    AlcFree( s_to_z );
  }

  /* check error and return */
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/*! 
* \return       New object with the rojection.
* \ingroup      WlzTransform
* \brief        Use the view transform to define a projection from
*		3D to 2D and then project the object onto this plane.
*		The object supplied to this function must be a 3D
*		spatial domain object (WLZ_3D_DOMAINOBJ) with either
*		no values or for integration WLZ_GREY_UBYTE values.
*		Integration will assign each output pixel the sum of
*		all input voxels mapped via either the domain density
*		or the voxel density.
*		The integration is controled by the integrate parameter
*		with valid values:
*		WLZ_PROJECT_INT_MODE_NONE - a "shadow domain" without values
*               is computed,
*		WLZ_PROJECT_INT_MODE_DOMAIN - the voxels of the domain are
*		integrated using
*		\f[
		p = \frac{1}{255} n d
                \f]
*		WLZ_PROJECT_INT_MODE_VALUES - the voxel values are integrated
*		using
*		\f[
		p = \frac{1}{255} \sum{l\left[v\right]}.
		\f]
*		Where
*		  \f$p\f$ is the projected image value,
*		  \f$n\f$ is the number of voxels projected for \f$p\f$,
*		  \f$d\f$ is the density of domain voxels,
*		  \f$l\f$ is the voxel value density look up table and
*		  \f$v\f$ is a voxel value.
* \param	obj			The given object.
* \param	vStr			Given view structure defining the
* 					projection plane.
* \param	intMod			This may take three values:
* 					WLZ_PROJECT_INT_MODE_NONE,
* 					WLZ_PROJECT_INT_MODE_DOMAIN or
* 					WLZ_PROJECT_INT_MODE_VALUES.
* \param	denDom			Density of domain voxels this value
* 					is not used unless the integration
* 					mode is WLZ_PROJECT_INT_MODE_DOMAIN.
* \param	denVal			Density look up table for object
* 					voxel density values which must be
* 					an array of 256 values. This may be
* 					NULL if the integration mode is not
* 					WLZ_PROJECT_INT_MODE_VALUES.
* \param	depth			If greater than zero, the projection
* 					depth perpendicular to the viewing
* 					plane.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzProjectObjToPlane(WlzObject *obj,
  				WlzThreeDViewStruct *vStr,
  				WlzProjectIntMode intMod,
				WlzUByte denDom, WlzUByte *denVal,
				double depth, WlzErrorNum *dstErr)
{
  int		nThr = 1,
		itvVal = 0;
  WlzIVertex2	prjSz;
  WlzIBox2	prjBox;
  double	pln[4];
  WlzObject	*bufObj = NULL,
  		*prjObj = NULL;
  WlzThreeDViewStruct *vStr1 = NULL;
  double	**vMat;
  WlzValues	nullVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzAffineTransform *rescaleTr = NULL;
  WlzGreyValueWSpace **gVWSp = NULL;
  void		***prjAry = NULL;
  const double	eps = 0.000001;
#ifdef WLZ_DEBUG_PROJECT3D_TIME
struct timeval	times[3];
#endif /* WLZ_DEBUG_PROJECT3D_TIME */

  nullVal.core = NULL;
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
  else if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(vStr == NULL)
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else if((intMod == WLZ_PROJECT_INT_MODE_VALUES) &&
          (obj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
#ifdef WLZ_DEBUG_PROJECT3D_TIME
  gettimeofday(times + 0, NULL);
#endif /* WLZ_DEBUG_PROJECT3D_TIME */
  /* Create new view transform without voxel scaling. The voxel scaling
   * is done after the projection. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((vStr1 = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) != NULL)
    {
      vStr1->fixed = vStr->fixed;
      vStr1->theta = vStr->theta;
      vStr1->phi = vStr->phi;
      vStr1->zeta = vStr->zeta;
      vStr1->dist = vStr->dist;
      vStr1->scale = vStr->scale;
      vStr1->voxelSize[0] = 1.0;
      vStr1->voxelSize[1] = 1.0;
      vStr1->voxelSize[2] = 1.0;
      vStr1->voxelRescaleFlg = 0;
      vStr1->interp = vStr->interp;
      vStr1->view_mode = vStr->view_mode;
      vStr1->up = vStr->up;
      vStr1->initialised = WLZ_3DVIEWSTRUCT_INIT_NONE;
      vMat = vStr1->trans->mat;
      errNum = WlzInit3DViewStructAffineTransform(vStr1);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = Wlz3DViewStructTransformBB(obj, vStr1);
      }
      if(errNum != WLZ_ERR_NONE)
      {
	WlzFree3DViewStruct(vStr1);
	vStr1 = NULL;
      }
    }
  }
  /* Compute bounding box of the projection. */
  if(errNum == WLZ_ERR_NONE)
  {
    prjBox.xMin = WLZ_NINT(vStr1->minvals.vtX);
    prjBox.yMin = WLZ_NINT(vStr1->minvals.vtY);
    prjBox.xMax = WLZ_NINT(vStr1->maxvals.vtX);
    prjBox.yMax = WLZ_NINT(vStr1->maxvals.vtY);
    prjSz.vtX = prjBox.xMax - prjBox.xMin + 1;
    prjSz.vtY = prjBox.yMax - prjBox.yMin + 1;
  }
  /* Compute post projection scaling. */
  if((errNum == WLZ_ERR_NONE) && (vStr->voxelRescaleFlg != 0))
  {
    WlzIBox2	sBox;
    WlzIVertex2 sSz;
    WlzThreeDViewStruct *vStr2;

    vStr2 = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      vStr2->fixed = vStr->fixed;
      vStr2->theta = vStr->theta;
      vStr2->phi = vStr->phi;
      vStr2->zeta = vStr->zeta;
      vStr2->dist = vStr->dist;
      vStr2->scale = vStr->scale;
      vStr2->voxelSize[0] = vStr->voxelSize[0];
      vStr2->voxelSize[1] = vStr->voxelSize[1];
      vStr2->voxelSize[2] = vStr->voxelSize[2];
      vStr2->voxelRescaleFlg = vStr->voxelRescaleFlg;
      vStr2->interp = vStr->interp;
      vStr2->view_mode = vStr->view_mode;
      vStr2->up = vStr->up;
      vStr2->initialised = WLZ_3DVIEWSTRUCT_INIT_NONE;
      errNum = WlzInit3DViewStructAffineTransform(vStr2);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = Wlz3DViewStructTransformBB(obj, vStr2);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	sBox.xMin = WLZ_NINT(vStr2->minvals.vtX);
	sBox.yMin = WLZ_NINT(vStr2->minvals.vtY);
	sBox.xMax = WLZ_NINT(vStr2->maxvals.vtX);
	sBox.yMax = WLZ_NINT(vStr2->maxvals.vtY);
	sSz.vtX = sBox.xMax - sBox.xMin + 1;
	sSz.vtY = sBox.yMax - sBox.yMin + 1;
        rescaleTr = WlzMakeAffineTransform(WLZ_TRANSFORM_2D_AFFINE, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        double	**m;

	m = rescaleTr->mat;
	m[0][0] = (sSz.vtX * eps) / (prjSz.vtX * eps);
	m[1][1] = (sSz.vtY * eps) / (prjSz.vtY * eps);
	m[0][2] = sBox.xMin - WLZ_NINT(m[0][0] * prjBox.xMin);
	m[1][2] = sBox.yMin - WLZ_NINT(m[1][1] * prjBox.yMin);
      }
      (void )WlzFree3DViewStruct(vStr2);
    }
  }
  /* Compute plane equation, used to clip intervals if depth was given. */
  if((errNum == WLZ_ERR_NONE) && (depth > eps))
  {
    Wlz3DViewGetPlaneEqn(vStr1, pln + 0, pln + 1, pln + 2, pln + 3);
  }
  /* Create rectangular projection array buffers, one for each thread,
   * also if integrating values create a grey value workspace per thread. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idB;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        nThr = omp_get_num_threads();
      }
    }
#endif
    if((prjAry = (void ***)AlcCalloc(nThr, sizeof(void **))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      if(intMod == WLZ_PROJECT_INT_MODE_NONE)
      {
	for(idB = 0; idB < nThr; ++idB)
	{
	  if(AlcUnchar2Calloc((WlzUByte ***)&(prjAry[idB]),
			      prjSz.vtY, prjSz.vtX) != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	}
      }
      else
      {
	for(idB = 0; idB < nThr; ++idB)
	{
	  if(AlcInt2Calloc((int ***)&(prjAry[idB]),
			   prjSz.vtY, prjSz.vtX) != ALC_ER_NONE)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	}
      }
    }
    if((errNum == WLZ_ERR_NONE) &&
       (intMod == WLZ_PROJECT_INT_MODE_VALUES))
    {
      itvVal = (WlzGreyTableIsTiled(obj->values.core->type) == 0);
      if(itvVal == 0)
      {
	if((gVWSp = AlcCalloc(nThr, sizeof(WlzGreyValueWSpace *))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idB = 0; idB < nThr; ++idB) 
	  {
	    gVWSp[idB] = WlzGreyValueMakeWSp(obj, &errNum);
	    if(gVWSp[idB]->gType != WLZ_GREY_UBYTE)
	    {
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	    }
	  }
	}
      }
    }
  }
  /* Scan through the 3D domain setting value in the projection array. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		pIdx,
    		pCnt;
    WlzDomain	*doms;
    WlzValues	*vals = NULL;

    doms = obj->domain.p->domains;
    if(itvVal)
    {
      vals = obj->values.vox->values;
    }
    pCnt = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(pIdx =  0; pIdx < pCnt; ++pIdx)
    {
      int	thrId = 0;

      if((errNum == WLZ_ERR_NONE) && (doms[pIdx].core != NULL))
      {
	WlzObject   *obj2;
	WlzGreyWSpace gWSp;
	WlzIntervalWSpace iWSp;
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

#ifdef _OPENMP
	thrId = omp_get_thread_num();
#endif
        obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, doms[pIdx],
	                   (vals)? vals[pIdx]: nullVal,
	                   NULL, NULL, &errNum2);
        if(errNum2 == WLZ_ERR_NONE)
	{
	  if(itvVal)
	  {
	    errNum2 = WlzInitGreyScan(obj2, &iWSp, &gWSp);
	  }
	  else
	  {
	    errNum2 = WlzInitRasterScan(obj2, &iWSp, WLZ_RASTERDIR_ILIC);
	  }
	}
        if(errNum2 == WLZ_ERR_NONE)
	{
	  double      plnZ,
	  	      vMZX,
	  	      vMZY;
	  WlzIVertex3 p0,
	    	      p1;

	  p0.vtZ = p1.vtZ = obj->domain.p->plane1 + pIdx;
	  vMZX = (vMat[0][2] * p0.vtZ) + vMat[0][3] - prjBox.xMin;
	  vMZY = (vMat[1][2] * p0.vtZ) + vMat[1][3] - prjBox.yMin;
	  plnZ = (pln[2] * p0.vtZ) + pln[3];
	  while(((itvVal == 0) &&
	        ((errNum2 = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE)) ||
	        ((itvVal != 0) &&
	        ((errNum2 = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE)))
	  {
	    int		skip = 0;
	    WlzDVertex2 q0,
	    		q1;

            p0.vtX = iWSp.lftpos;
	    p1.vtX = iWSp.rgtpos;
	    p0.vtY = p1.vtY = iWSp.linpos;
	    if(depth > eps)
	    {
	      int	c;
	      double	d0,
	      		d1,
			plnYZ;

	      /* Clip the 3D line segment p0,p1 using the plane equation. */
	      plnYZ = (pln[1] * p0.vtY) + plnZ;
	      d0 = (pln[0] * p0.vtX) + plnYZ;
	      d1 = (pln[0] * p1.vtX) + plnYZ;
	      c = ((d1 >  depth) << 3) || ((d0 >  depth) << 2) ||
		  ((d1 < -depth) << 1) ||  (d0 < -depth);
	      if(c)
	      {
		if((c == 3) || (c == 12)) /* 00-- or ++00 */
		{
		  /* Both out of range, so don't render. */
		  skip = 1;
		}
		else
		{
		  if(fabs(pln[0]) > eps)
		  {
		    double	plnX;

		    plnX = -1.0 / pln[0];
		    if((c &  1) != 0)      /* x0x- */
		    {
		      p0.vtX = plnX * (plnYZ + depth);
		    }
		    else if((c &  4) != 0) /* x+x0 */
		    {
		      p0.vtX = plnX * (plnYZ - depth);
		    }
		    if((c &  2) != 0)      /* 0x-x */
		    {
		      p1.vtX = plnX * (plnYZ + depth);
		    }
		    else if((c &  8) != 0) /* +x0x */
		    {
		      p1.vtX = plnX * (plnYZ - depth);
		    }
		  }
		}
	      }
	    }
	    if(skip == 0)
	    {
	      q0.vtX = (vMat[0][0] * p0.vtX) + (vMat[0][1] * p0.vtY) + vMZX;
	      q0.vtY = (vMat[1][0] * p0.vtX) + (vMat[1][1] * p0.vtY) + vMZY;
	      q1.vtX = (vMat[0][0] * p1.vtX) + (vMat[0][1] * p1.vtY) + vMZX;
	      q1.vtY = (vMat[1][0] * p1.vtX) + (vMat[1][1] * p1.vtY) + vMZY;
	      switch(intMod)
	      {
		case WLZ_PROJECT_INT_MODE_NONE:
		  {
		    WlzIVertex2 u0,
				u1;

		    WLZ_VTX_2_NINT(u0, q0);
		    WLZ_VTX_2_NINT(u1, q1);
		    WlzProjectObjLine((WlzUByte **)(prjAry[thrId]), u0, u1);
		  }
		  break;
		case WLZ_PROJECT_INT_MODE_DOMAIN:
		  {
		    int	        np,
				nq;
		    WlzDVertex3 dq;
		    WlzIVertex2 u0,
				u1;

		    WLZ_VTX_2_NINT(u0, q0);
		    WLZ_VTX_2_NINT(u1, q1);
		    WLZ_VTX_2_SUB(dq, q0, q1);
		    np = denDom * (iWSp.rgtpos - iWSp.lftpos + 1);
		    nq = (int )ceil(WLZ_VTX_2_LENGTH(dq) + eps);
		    WlzProjectObjLineDom((int **)(prjAry[thrId]), np / nq,
					 u0, u1);
		  }
		  break;
		case WLZ_PROJECT_INT_MODE_VALUES:
		  if(itvVal)
		  {
		    WlzProjectObjLineVal((int **)(prjAry[thrId]), denVal,
					 gWSp.u_grintptr.ubp, NULL,
					 vMat, vMZX, vMZY, p0, p1);
		  }
		  else
		  {
		    WlzProjectObjLineVal((int **)(prjAry[thrId]), denVal,
					 NULL, gVWSp[thrId],
					 vMat, vMZX, vMZY, p0, p1);
		  }
		  break;
	      }
	    }
	  }
	  (void )WlzEndGreyScan(&iWSp, &gWSp);
	  if(errNum2 == WLZ_ERR_EOO)
	  {
	    errNum2 = WLZ_ERR_NONE;
	  }
	}
	(void )WlzFreeObj(obj2);
	if(errNum2 != WLZ_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp critical
	  {
#endif
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
#ifdef _OPENMP
	  }
#endif
	}
      }
    }
  }
  /* Free grey value workspaces if they were created. */
  if(gVWSp)
  {
    int		idB;

    for(idB = 0; idB < nThr; ++idB) 
    {
      WlzGreyValueFreeWSp(gVWSp[idB]);
    }
    AlcFree(gVWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idB;
    size_t	idC,
    		bufSz;
    WlzGreyP	buf0,
      		buf1;
    WlzIVertex2	prjOrg;

    prjOrg.vtX = prjBox.xMin;
    prjOrg.vtY = prjBox.yMin;
    bufSz = prjSz.vtX * prjSz.vtY;
    for(idB = 1; idB < nThr; ++idB)
    {
      if(intMod == WLZ_PROJECT_INT_MODE_NONE)
      {
        buf0.ubp = ((WlzUByte ***)(prjAry))[0][0],
	buf1.ubp = ((WlzUByte ***)(prjAry))[idB][0];
	for(idC = 0; idC < bufSz; ++idC)
	{
	  buf0.ubp[idC] += buf1.ubp[idC];
	}
      }
      else
      {
        buf0.inp = ((int ***)(prjAry))[0][0],
	buf1.inp = ((int ***)(prjAry))[idB][0];
	for(idC = 0; idC < bufSz; ++idC)
	{
	  buf0.inp[idC] += buf1.inp[idC];
	}
      }
    }
    switch(intMod != WLZ_PROJECT_INT_MODE_NONE)
    {
      buf0.inp = ((int ***)(prjAry))[0][0];
      for(idC = 0; idC < bufSz; ++idC)
      {
	buf0.inp[idC] /= 256;
      }
    }
    if(intMod == WLZ_PROJECT_INT_MODE_NONE)
    {
      bufObj = WlzAssignObject(
	       WlzFromArray2D((void **)(prjAry[0]), prjSz, prjOrg,
			      WLZ_GREY_UBYTE, WLZ_GREY_UBYTE,
			      0.0, 1.0, 1, 0, &errNum), NULL);
    }
    else
    {
      bufObj = WlzAssignObject(
	       WlzFromArray2D((void **)(prjAry[0]), prjSz, prjOrg,
			      WLZ_GREY_INT, WLZ_GREY_INT,
			      0.0, 1.0, 1, 0, &errNum), NULL);
    }
  }
  /* Free the projection array(s). */
  if(prjAry)
  {
    int		idB;

    for(idB = 0; idB < nThr; ++idB)
    {
      (void )Alc2Free((prjAry[idB]));
    }
    AlcFree(prjAry);
  }
  /* Make return object using threshold. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzPixelV	tV;
    WlzObject	*tObj = NULL;

    tV.type = WLZ_GREY_UBYTE;
    tV.v.ubv = 1;
    tObj = WlzAssignObject(
	   WlzThreshold(bufObj, tV, WLZ_THRESH_HIGH, &errNum), NULL);
    if(tObj)
    {
      if(intMod == WLZ_PROJECT_INT_MODE_NONE)
      {
	prjObj = WlzMakeMain(tObj->type, tObj->domain, nullVal,
			     NULL, NULL, &errNum);
      }
      else
      {
	prjObj = WlzMakeMain(tObj->type, tObj->domain, tObj->values,
			     NULL, NULL, &errNum);
      }
    }
    (void )WlzFreeObj(tObj);
  }
  (void )WlzFreeObj(bufObj);
  (void )WlzFree3DViewStruct(vStr1);
  /* Scale image. */
  if(rescaleTr != NULL)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      WlzObject	*tObj = NULL;

      tObj = WlzAffineTransformObj(prjObj, rescaleTr, 
				   WLZ_INTERPOLATION_NEAREST, &errNum);
      
      (void )WlzFreeObj(prjObj);
      prjObj = tObj;
    }
    (void )WlzFreeAffineTransform(rescaleTr);
  }
#ifdef WLZ_DEBUG_PROJECT3D_TIME
  gettimeofday(times + 1, NULL);
  ALC_TIMERSUB(times + 1, times + 0, times + 2);
  (void )fprintf(stderr, "WlzGetProjectionFromObject: Elapsed time = %g\n",
                 times[2].tv_sec + (0.000001 * times[2].tv_usec));
#endif /* WLZ_DEBUG_PROJECT3D_TIME */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(prjObj);
}

/*!
* \ingroup	WlzTransform
* \brief	Sets values in the array to 1 on the straight line segment
* 		between p0 and p1. This is just an implementation of
* 		Bresenham's line algorithm.
* 		in 2D.
* \param	ary			Rectangular array of bytes to set
* 					projected line values in.
* \param	p0			Projected segment start.
* \param	p1			Projected segment end.
*/
static void	WlzProjectObjLine(WlzUByte **ary,
				  WlzIVertex2 p0, WlzIVertex2 p1)
{
  int		e;
  WlzIVertex2	d,
  		s;
  WlzUByte	*a;

  WLZ_VTX_2_SUB(d, p1, p0);
  s.vtX = ((d.vtX > 0) << 1) - 1;
  s.vtY = ((d.vtY > 0) << 1) - 1;
  WLZ_VTX_2_ABS(d,d);
  a = *(ary + p0.vtY);
  e = ((d.vtX > d.vtY)? d.vtX: -d.vtY) / 2;
  while(1)
  {
    int		e2;

    *(a + p0.vtX) = 1;
    if((p0.vtX == p1.vtX) && (p0.vtY == p1.vtY))
    {
      break;
    }
    e2 = e;
    if(e2 > -d.vtX)
    {
      e      -= d.vtY;
      p0.vtX += s.vtX;
    }
    if(e2 < d.vtY)
    {
      e      += d.vtX;
      p0.vtY += s.vtY;
      a = *(ary + p0.vtY);
    }
  }
}

/*!
* \ingroup	WlzTransform
* \brief	Increments values in the array on the straight line segment
* 		between p0 and p1. This is just an implementation of
* 		Bresenham's line algorithm in 2D.
* \param	ary			Rectangular array of bytes to set
* 					projected line values in.
* \param	vpp			Voxels per pixel, the increment.
* \param	p0			Projected segment start.
* \param	p1			Projected segment end.
*/
static void	WlzProjectObjLineDom(int **ary, int vpp,
				     WlzIVertex2 p0, WlzIVertex2 p1)
{
  int		e;
  WlzIVertex2	d,
  		s;
  int		*a;

  WLZ_VTX_2_SUB(d, p1, p0);
  s.vtX = ((d.vtX > 0) << 1) - 1;
  s.vtY = ((d.vtY > 0) << 1) - 1;
  WLZ_VTX_2_ABS(d,d);
  e = ((d.vtX > d.vtY)? d.vtX: -d.vtY) / 2;
  a = *(ary + p0.vtY);
  while(1)
  {
    int		e2;

    *(a + p0.vtX) += vpp;
    if((p0.vtX == p1.vtX) && (p0.vtY == p1.vtY))
    {
      break;
    }
    e2 = e;
    if(e2 > -d.vtX)
    {
      e      -= d.vtY;
      p0.vtX += s.vtX;
    }
    if(e2 < d.vtY)
    {
      e      += d.vtX;
      p0.vtY += s.vtY;
      a = *(ary + p0.vtY);
    }
  }
}

/*!
* \ingroup	WlzTransform
* \brief	Increments values in the array on the straight line segment
* 		between (x0, y0) and (x1, y1) using the given look up table
* 		and the corresponding grey value as an index into this look
* 		up table.
* \param	ary			Rectangular array of ints to set
* 					projected line values in.
* \param	lut			Value density look up table with
* 					256 entries.
* \param	gP			Grey pointer for access to the
* 					object's grey values.
* \param	gVWSp			Grey value workspace for access to
* 					the object grey values, not used if
* 					the grey pointer is non null.
* \param	vMat			View transform matrix.
* \param	vMZX			Precomputed matrix z -> x.
* \param	vMZY			Precomputed matrix z -> y.
* \param	p0			Line segment start.
* \param	p1			Line segment end.
*/
static void	WlzProjectObjLineVal(int **ary, WlzUByte lut[256],
  				     WlzUByte *gP,
				     WlzGreyValueWSpace *gVWSp,
				     double **vMat, double vMZX, double vMZY,
				     WlzIVertex3 p0, WlzIVertex3 p1)
{
#ifdef WLZ_FAST_CODE
  int		v00,
  		v01,
		vMZXI,
		vMZYI;
#else
  WlzDVertex2	q;
#endif
  WlzIVertex2	u;

  vMZX += vMat[0][1] * p0.vtY;
  vMZY += vMat[1][1] * p0.vtY;
#ifdef WLZ_FAST_CODE
  /*
   * Use integer arithmetic for evaluating the transform instead of floating
   * point and then rounding to the nearest integer, which saves about 20% of
   * the execution time for the grey pointer code (mainly be avoiding nint()).
   * Integer multiplication and division are used with a factor of 2^10 which
   * the compiller should recognise and be able to do an arithmetic shift
   * right for division, we can't use >> as it's behaviour for signed ints
   * is platfrom dependant.
   */
  vMZXI = WLZ_NINT(vMZX);
  vMZYI = WLZ_NINT(vMZY);
  v00 = WLZ_NINT(vMat[0][0] * 1024.0);
  v01 = WLZ_NINT(vMat[1][0] * 1024.0);
#endif
  if(gP)
  {
    /* Use the interval of values given by the grey pointer. */
    if(lut)
    {
      while(1)
      {
#ifdef WLZ_FAST_CODE
	u.vtX = ((v00 * p0.vtX) / 1024) + vMZXI;
	u.vtY = ((v01 * p0.vtX) / 1024) + vMZYI;
#else
	q.vtX = (vMat[0][0] * p0.vtX) + vMZX;
	q.vtY = (vMat[1][0] * p0.vtY) + vMZY;
	WLZ_VTX_2_NINT(u, q);
#endif
	*(*(ary + u.vtY) + u.vtX) += lut[*gP++];
	if(++(p0.vtX) > p1.vtX)
	{
	  break;
	}
      }
    }
    else
    {
      while(1)
      {
#ifdef WLZ_FAST_CODE
	u.vtX = ((v00 * p0.vtX) / 1024) + vMZXI;
	u.vtY = ((v01 * p0.vtX) / 1024) + vMZYI;
#else
	q.vtX = (vMat[0][0] * p0.vtX) + vMZX;
	q.vtY = (vMat[1][0] * p0.vtY) + vMZY;
	WLZ_VTX_2_NINT(u, q);
#endif
	*(*(ary + u.vtY) + u.vtX) += 255 - *gP++;
	if(++(p0.vtX) > p1.vtX)
	{
	  break;
	}
      }
    }
  }
  else
  {
    /* Need to use random access to the grey values since no grey pointer. */
    if(lut)
    {
      while(1)
      {
	WlzGreyValueGet(gVWSp, p0.vtZ, p0.vtY, p0.vtX);
#ifdef WLZ_FAST_CODE
	u.vtX = ((v00 * p0.vtX) / 1024) + vMZXI;
	u.vtY = ((v01 * p0.vtX) / 1024) + vMZYI;
#else
	q.vtX = (vMat[0][0] * p0.vtX) + vMZX;
	q.vtY = (vMat[1][0] * p0.vtY) + vMZY;
	WLZ_VTX_2_NINT(u, q);
#endif
	*(*(ary + u.vtY) + u.vtX) += lut[gVWSp->gVal[0].ubv];
	if(++(p0.vtX) > p1.vtX)
	{
	  break;
	}
      }
    }
    else
    {
      while(1)
      {
	WlzGreyValueGet(gVWSp, p0.vtZ, p0.vtY, p0.vtX);
#ifdef WLZ_FAST_CODE
	u.vtX = ((v00 * p0.vtX) / 1024) + vMZXI;
	u.vtY = ((v01 * p0.vtX) / 1024) + vMZYI;
#else
	q.vtX = (vMat[0][0] * p0.vtX) + vMZX;
	q.vtY = (vMat[1][0] * p0.vtY) + vMZY;
	WLZ_VTX_2_NINT(u, q);
#endif
	*(*(ary + u.vtY) + u.vtX) += 255 - gVWSp->gVal[0].ubv;
	if(++(p0.vtX) > p1.vtX)
	{
	  break;
	}
      }
    }
  }
}
