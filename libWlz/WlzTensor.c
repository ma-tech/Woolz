#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTensor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTensor.c
* \author       Bill Hill
* \date         January 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Functions which derive and manipulate tensor quantities.
* \ingroup	WlzFeatures
*/
#include <stdio.h>
#include <Wlz.h>

#ifdef _OPENMP
#include <omp.h>

#define WLZ_TENSOR_OMP_CHUNKSZ 4096   /* To avoid parallelising small loops. */
#endif

static WlzErrorNum		WlzTensorGetComponentValues2D(
				  WlzObject *rObj,
				  WlzObject *tObj,
				  int cpt,
				  int pln);
static WlzErrorNum		WlzTensorSetComponentValues2D(
				  WlzObject *fObj,
				  WlzObject *tObj,
				  int cpt,
				  int pln);
static WlzErrorNum		WlzTensorComponentValues3D(
				  WlzObject *rObj,
				  WlzObject *tObj,
				  int cpt,
				  int set);
static void			WlzDGTensorSDFeatEigenVecItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt);
static void			WlzDGTensorSDFeatDetJacItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt);
static void			WlzDGTensorSDFeatEigenValItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt);
static void			WlzDGTensorSDFeatDetJacPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP);
static void			WlzDGTensorSDFeatEigenVecPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP);
static void			WlzDGTensorSDFeatEigenValPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP);
static void			WlzDGTensorToRCauchyGreen(
				  AlgMatrix cg,
				  WlzGreyP gP,
				  WlzGreyType gT);
static WlzObject		*WlzCMeshDGTensor3D(
				  WlzObject *cObj,
				  int invert,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshDGTensorAtPts3D(
				  WlzObject *cObj,
				  int invert,
				  WlzDVertex3 sd,
				  int dither,
				  WlzErrorNum *dstErr);
static void	 		WlzCMeshDGToStrainTensor3D(
				  WlzObject *tObj);
static void 			WlzPtsDGToStrainTensor3D(
				  WlzObject *tObj);
static void			WlzCMeshElmSetDGTensor3D(
				  WlzCMeshElm3D *elm,
				  int invert,
				  WlzIndexedValues *cIxv,
				  double *ten);

/*!
* \return	New compound object containing feature objects or NULL on
* 		error.
* \ingroup	WlzFeatures
* \brief	Computes features of a deformation gradient tensor field
* 		and creates a compound object, each element of which is
* 		a required feature of the given deformation gradient field
* 		object. The feature object may be either a points object
* 		or a spatial domain object as required.
* 		See WlzDGTensorSDFeature() which this function calls to
* 		compute the feature objects.
* \param	mObj		Given deformation gradient field object. This
* 				must be a 3D object in which the (tiled) voxel
* 				values are each a nine element Jacobian
* 				deformation gradient tensor.
* \param	features	Mask with required features.
* \param	points		Generate a points object rather than a field.
* \param	dMin		If generating a points the points will have
* 				this minimum seperation distance (if less
* 				than 1.0 then will be set to 1.0).
* \param	dither		If generating a points object dither the points
* 				within this range (no dithering if values are
* 				zero). The dithering is always done in unit
* 				voxel space.
* \param	smooth		If smoothing values are > 0.0 apply a Gaussian
* 				smoothing filter with these sigma values to
* 				spatial field feature values. There is no
* 				smoothing of point feature values.
* \param	voxScaling	Use voxel size for points objects if non zero.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzDGTensorFeatures(
				  WlzObject *mObj,
				  unsigned int features,
				  int points,
				  double dMin,
				  WlzDVertex3 dither,
                                  WlzDVertex3 smooth,
				  int voxScaling,
				  WlzErrorNum *dstErr)
{
  WlzDomain	fDom;
  WlzObject	*rObj = NULL;
  WlzObjectType	rObjType = WLZ_NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  fDom.core = NULL;
  /* Check that the input object is valid. */
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(mObj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(!WlzGreyTableIsTiled(mObj->values.core->type))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    WlzTiledValues	*tv;

    tv = mObj->values.t;
    if(tv->dim != 3)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
    else if (tv->vRank == 1)
    {
      /* Allow 1x9 values per voxel. */
      if(tv->vDim[0] != 9)
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
    }
    else if (tv->vRank == 2)
    {
      /* Allow 3x3 values per voxel. */
      if((tv->vDim[0] != 3) || (tv->vDim[1] != 3))
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
    }
    else
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  /* Find domain within which to compute the features. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(points)
    {
      double	ditherSqLen;
      const double minDitherSqLn = 1.0 - 1.0e-06;

      rObjType = WLZ_POINTS;
      if(dMin < 1.0)
      {
	dMin = 1.0;
      }
      /* Create points with floating point positions. Voxel scaling
       * is done latter if required. */
      fDom.pts = WlzPointsFromDomObj(mObj, dMin, 1, 0,
                                     0, 0.0, 0.0, 0.0, &errNum);
      ditherSqLen = WLZ_VTX_3_SQRLEN(dither);
      if((errNum == WLZ_ERR_NONE) && (ditherSqLen > minDitherSqLn))
      {
        WlzPoints *pts;             

	pts = WlzPointsDither(fDom.pts, dither, mObj, &errNum);
	(void )WlzFreeDomain(fDom);
	fDom.pts = pts;
      }
    }
    else
    {
      rObjType = WLZ_3D_DOMAINOBJ;
      fDom = mObj->domain;
    }
    (void )WlzAssignDomain(fDom, NULL);
  }
  /* Create return compound object. */
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = (WlzObject *)
           WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1,
	                        1, WLZ_DGTENSOR_FEATURE_LIMIT, NULL,
				rObjType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx,
    		feat;
    WlzCompoundArray *cpd;

    idx = 0;
    cpd = (WlzCompoundArray *)rObj;
    for(feat = 1; feat < WLZ_DGTENSOR_FEATURE_LIMIT; ++feat)
    {
      if(features & WLZ_DGTENSOR_FEATURE_MASK(feat))
      {
	if(rObjType == WLZ_3D_DOMAINOBJ)
	{
	  cpd->o[idx] = WlzAssignObject(
                WlzDGTensorSDFeature(mObj, fDom, feat, smooth, &errNum),
		NULL);
	}
	else /* rObjType == WLZ_POINTS */
	{
	  cpd->o[idx] = WlzAssignObject(
                WlzDGTensorPDFeature(mObj, fDom, feat, &errNum),
		NULL);
	}
	++idx;
      }
    }
    cpd->n = idx;
  }
  /* Scale point locations if required. */
  if((errNum == WLZ_ERR_NONE) && (rObjType != WLZ_3D_DOMAINOBJ) && voxScaling)
  {
    int		idP,
    		nPoints;
    WlzDVertex3 voxSz;
    WlzDVertex3 *pts;

    pts = fDom.pts->points.d3;
    nPoints = fDom.pts->nPoints;
    voxSz.vtX = mObj->domain.p->voxel_size[0];
    voxSz.vtY = mObj->domain.p->voxel_size[1];
    voxSz.vtZ = mObj->domain.p->voxel_size[2];
#ifdef _OPENMP
#pragma omp parallel for schedule(static, WLZ_TENSOR_OMP_CHUNKSZ)
#endif
    for(idP = 0; idP < nPoints; ++idP)
    {
      pts[idP].vtX *= voxSz.vtX;
      pts[idP].vtY *= voxSz.vtY;
      pts[idP].vtZ *= voxSz.vtZ;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(rObj);
    rObj = NULL;
  }
  (void )WlzFreeDomain(fDom);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object containing features or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes a feature of a deformation gradient tensor field
* 		throughout the given 3D spatial domain.
* 		This function assumes all input parameters to be valid,
* 		when this is not certain use WlzDGTensorFeatures().
* 		The returned object will have it's name property set to
* 		an appropriate name for the required feature:
*		WLZ_DGTENSOR_FEATURE_DETJAC ("jacobian"),
*		WLZ_DGTENSOR_FEATURE_EIGENVEC ("eigen vectors") or
*		WLZ_DGTENSOR_FEATURE_EIGENVAL ("eigen values").
* \param	mObj		Given deformation gradient field object. This
* 				must be a 3D object in which the (tiled) voxel
* 				values are each a nine element Jacobian
* 				deformation gradient tensor.
* \param	fDom		Given 3D spatial domain (WlzPlaneDomain).
* \param	feat		Required feature.
* \param	smooth		If values are > 0.0 apply a Gaussian smoothing
* 				filter with these sigma values to the feature
* 				values.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzDGTensorSDFeature(
				  WlzObject *mObj,
				  WlzDomain fDom,
				  WlzDGTensorFeatureType feat,
				  WlzDVertex3 smooth,
				  WlzErrorNum *dstErr)
{
  int		plMin,
  		plMax;
  unsigned int	vRank;
  unsigned int  vDim[2];
  char		*pName = NULL;
  WlzObject	*rObj = NULL;
  WlzPlaneDomain *dPDom,
  		 *sPDom;
  WlzValues	nullVal;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-06;
  const size_t  tlSz = WLZ_TILEDVALUES_TILE_SIZE;

  nullVal.core = NULL;
  /* Create values for the feature. */
  switch(feat)
  {
    case WLZ_DGTENSOR_FEATURE_DETJAC:
      vRank = 0;
      vDim[0] = 0;
      pName = "jacobian";
      break;
    case WLZ_DGTENSOR_FEATURE_EIGENVEC:
      vRank = 2;
      vDim[0] = 3;
      vDim[1] = 3;
      pName = "eigen vectors";
      break;
    case WLZ_DGTENSOR_FEATURE_EIGENVAL:
      vRank = 1;
      vDim[0] = 3;
      pName = "eigen values";
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*fObj = NULL;

    fObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, fDom, nullVal, NULL, NULL, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzPixelV	bgdV;

      bgdV.type = WLZ_GREY_DOUBLE;
      bgdV.v.dbv = 0.0;
      rObj = WlzMakeTiledValuesObj3D(fObj, tlSz, 0, WLZ_GREY_DOUBLE,
				     vRank, vDim, bgdV, &errNum);
    }
    (void )WlzFreeObj(fObj);
  }
  /* Compute feature values throughout the domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP;

    dPDom = fDom.p;
    sPDom = mObj->domain.p;
    plMin = ALG_MAX(dPDom->plane1, sPDom->plane1);
    plMax = ALG_MIN(dPDom->lastpl, sPDom->lastpl);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idP = plMin; idP <= plMax; ++idP)
    {
      if(errNum == WLZ_ERR_NONE)
      {
        WlzDomain *dDom2D,
	          *sDom2D;
	AlgMatrix cgt;
        WlzErrorNum errNum2 = WLZ_ERR_NONE;

	cgt.core = NULL;
	if((cgt.rect = AlgMatrixRectNew(3, 3, NULL)) == NULL)
	{
	  errNum2 = WLZ_ERR_MEM_ALLOC;
	}
	if((errNum2 == WLZ_ERR_NONE) &&
	   ((dDom2D = dPDom->domains + idP - dPDom->plane1) != NULL) &&
	   ((sDom2D = sPDom->domains + idP - sPDom->plane1) != NULL) &&
	   ((*dDom2D).core != NULL) &&
	   ((*sDom2D).core != NULL))
        {
	  WlzObject *dObj2D = NULL,
	            *sObj2D = NULL;
          WlzGreyValueWSpace *dGVWSp = NULL,    /* 3D grey value workspaces. */
	                     *sGVWSp = NULL;

	  if(((dObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *dDom2D, nullVal,
	                            NULL, NULL, &errNum2)) != NULL) &&
	     ((sObj2D = WlzMakeMain(WLZ_2D_DOMAINOBJ, *sDom2D, nullVal,
	                            NULL, NULL, &errNum2)) != NULL) &&
	     ((dGVWSp = WlzGreyValueMakeWSp(rObj, &errNum2)) != NULL) &&
	     ((sGVWSp = WlzGreyValueMakeWSp(mObj, &errNum2)) != NULL))
	  {
	    int		lnMin,
			lnMax;
	    WlzIntervalWSpace dIWSp,
			      sIWSp;

	    lnMin = ALG_MAX(dObj2D->domain.i->line1,
	                    sObj2D->domain.i->line1);
	    lnMax = ALG_MIN(dObj2D->domain.i->lastln,
	                    sObj2D->domain.i->lastln);
	    if(lnMin <= lnMax)
	    {
	      if((errNum2 = WlzInitRasterScan(dObj2D, &dIWSp,
				WLZ_RASTERDIR_ILIC) == WLZ_ERR_NONE) &&
		 (errNum2 = WlzInitRasterScan(sObj2D, &sIWSp,
				WLZ_RASTERDIR_ILIC) == WLZ_ERR_NONE))
	      {
		/* Process intervals of the source into destination. */
		if((errNum2 = WlzNextInterval(&sIWSp)) == WLZ_ERR_NONE)
		{
		  do
		  {
		    errNum2 = WlzNextInterval(&dIWSp);
		  } while((errNum2 == WLZ_ERR_NONE) && (dIWSp.linpos < lnMin));
		}
		while((errNum2 == WLZ_ERR_NONE) && (dIWSp.linpos <= lnMax))
	        {
		  while((errNum2 == WLZ_ERR_NONE) &&
		        (sIWSp.linpos == dIWSp.linpos))
		  {
		    int		isS;
		    WlzInterval	itv;

		    isS = WlzIWSpIntersection(&itv, &dIWSp, &sIWSp, NULL);
		    switch(feat)
		    {
		      case WLZ_DGTENSOR_FEATURE_DETJAC:
		        WlzDGTensorSDFeatDetJacItv(dGVWSp, sGVWSp, cgt,
				idP, dIWSp.linpos, itv.ileft, itv.iright);
		        break;
		      case WLZ_DGTENSOR_FEATURE_EIGENVEC:
		        WlzDGTensorSDFeatEigenVecItv(dGVWSp, sGVWSp, cgt,
				idP, dIWSp.linpos, itv.ileft, itv.iright);
		        break;
		      case WLZ_DGTENSOR_FEATURE_EIGENVAL:
		        WlzDGTensorSDFeatEigenValItv(dGVWSp, sGVWSp, cgt,
				idP, dIWSp.linpos, itv.ileft, itv.iright);
		        break;
		      default:
		        /* Already checked for. */
			break;

		    }
		    errNum2 = WlzNextInterval((isS)? &sIWSp: &dIWSp);
		  }
		  if(errNum2 == WLZ_ERR_NONE)
		  {
		    if(sIWSp.linpos < dIWSp.linpos)
		    {
		      errNum2 = WlzNextInterval(&sIWSp);
		    }
		    else if(dIWSp.linpos < sIWSp.linpos)
		    {
		      errNum2 = WlzNextInterval(&dIWSp);
		    }
		  }
		}
		if(errNum2 == WLZ_ERR_EOO)
		{
		  errNum2 = WLZ_ERR_NONE;
		}
	      }
	    }
	  }
	  (void )WlzFreeObj(dObj2D);
	  (void )WlzFreeObj(sObj2D);
	  WlzGreyValueFreeWSp(dGVWSp);
	  WlzGreyValueFreeWSp(sGVWSp);
	}
        if(errNum2 != WLZ_ERR_NONE)
	{
#ifdef _OPENMP
#pragma omp critical (WlzDGTensorSDFeature)
#endif
	  {
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
	  }
	}
        AlgMatrixFree(cgt);
      }
    }
  }
  /* If smoothing is needed, apply a Gaussian smoothing kernel. */
  if(errNum == WLZ_ERR_NONE &&
     ((smooth.vtX > eps) || (smooth.vtY > eps) || (smooth.vtZ > eps)))
  {
    errNum = WlzTensorSmooth(rObj, smooth);
  }
  /* Set object name property. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzProperty      prop = {0};
    WlzPropertyList *pList = NULL;

    if((pList = WlzMakePropertyList(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      prop.name = WlzMakeNameProperty(pName, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
                               WlzFreePropertyListEntry) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rObj->plist = WlzAssignPropertyList(pList, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(pList->list == NULL)
      {
	WlzFreeProperty(prop);
      }
      else
      {
	(void )WlzFreeProperty(prop);
	(void )WlzFreePropertyList(pList);
      }
    }
  }
  /* Clear up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object containing features or NULL on error.
* \ingroup	WlzFeatures
* \brief	Computes a feature of a deformation gradient tensor field
* 		throughout the given 3D points domain.
* 		This function assumes all input parameters to be valid,
* 		when this is not certain use WlzDGTensorFeatures().
* \param	mObj		Given deformation gradient field object. This
* 				must be a 3D object in which the (tiled) voxel
* 				values are each a nine element Jacobian
* 				deformation gradient tensor.
* \param	fDom		Given 3D points domain (WlzPoints).
* \param	feat		Required feature.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzDGTensorPDFeature(
				  WlzObject *mObj,
				  WlzDomain fDom,
				  WlzDGTensorFeatureType feat,
				  WlzErrorNum *dstErr)
{
  int		vRank,
  		nThr = 1;
  int  		vDim[2];
  char		*pName = NULL;
  WlzObject	*rObj = NULL;
  AlgMatrix 	*cgtAry = NULL;
  WlzGreyValueWSpace **sGVWSpAry = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Create values for the feature. */
  switch(feat)
  {
    case WLZ_DGTENSOR_FEATURE_DETJAC:
      vRank = 0;
      vDim[0] = 0;
      pName = "jacobian";
      break;
    case WLZ_DGTENSOR_FEATURE_EIGENVEC:
      vRank = 2;
      vDim[0] = 3;
      vDim[1] = 3;
      pName = "eigen vectors";
      break;
    case WLZ_DGTENSOR_FEATURE_EIGENVAL:
      vRank = 1;
      vDim[0] = 3;
      pName = "eigen values";
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
	nThr = omp_get_num_threads();
      }
    }
#endif
    if(((cgtAry = (AlgMatrix *)
                  AlcCalloc(nThr, sizeof(AlgMatrix))) == NULL) ||
       ((sGVWSpAry = (WlzGreyValueWSpace **)
                     AlcCalloc(nThr, sizeof(WlzGreyValueWSpace *))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int 	idT;

      for(idT = 0; idT < nThr; ++idT)
      {
        if(((sGVWSpAry[idT] = WlzGreyValueMakeWSp(mObj, NULL)) == NULL) ||
	   ((cgtAry[idT].rect = AlgMatrixRectNew(3, 3, NULL)) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues fVal;

    fVal.pts = WlzMakePointValues(fDom.pts->nPoints,
				  vRank, vDim, WLZ_GREY_DOUBLE,
				  &errNum);
    if(errNum == WLZ_ERR_NONE) 
    {
      rObj = WlzMakeMain(WLZ_POINTS, fDom, fVal, NULL, NULL, &errNum);
    }
  }
  /* Compute feature values throughout the domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP;
    WlzPointValues *dPV;

    dPV = rObj->values.pts;
    if(errNum == WLZ_ERR_NONE)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(idP = 0; idP < fDom.pts->maxPoints; ++idP)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  int	      idT = 0;
	  WlzDVertex3 p;
	  AlgMatrix   cgt;
	  WlzGreyValueWSpace *gVWSp;
	  WlzErrorNum errNum2 = WLZ_ERR_NONE;

#ifdef _OPENMP
	  idT = omp_get_thread_num();
#endif
	  cgt = cgtAry[idT];
	  gVWSp = sGVWSpAry[idT];
	  /* Get point location. */
	  p = fDom.pts->points.d3[idP];
	  /* Set grey pointer for point location. */
          WlzGreyValueGet(gVWSp, p.vtZ, p.vtY, p.vtX);
	  if(errNum2 == WLZ_ERR_NONE)
	  {
	    switch(feat)
	    {
	      case WLZ_DGTENSOR_FEATURE_DETJAC:
		WlzDGTensorSDFeatDetJacPnt(dPV, gVWSp->gPtr[0],
					   dPV->vType, gVWSp->gType,
					   cgt, idP);
		break;
	      case WLZ_DGTENSOR_FEATURE_EIGENVEC:
		WlzDGTensorSDFeatEigenVecPnt(dPV, gVWSp->gPtr[0],
					     dPV->vType, gVWSp->gType,
					     cgt, idP);
		break;
	      case WLZ_DGTENSOR_FEATURE_EIGENVAL:
		WlzDGTensorSDFeatEigenValPnt(dPV, gVWSp->gPtr[0], 
					     dPV->vType, gVWSp->gType,
					     cgt, idP);
		break;
	      default:
		/* Already checked for. */
		break;

	    }
	  }
	  if(errNum2 != WLZ_ERR_NONE)
	  {
#ifdef _OPENMP
#pragma omp critical (WlzDGTensorPDFeature)
#endif
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
	  }
	}
      }
    }
  }
  /* Set object name property. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzProperty      prop = {0};
    WlzPropertyList *pList = NULL;

    if((pList = WlzMakePropertyList(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      prop.name = WlzMakeNameProperty(pName, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
                               WlzFreePropertyListEntry) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rObj->plist = WlzAssignPropertyList(pList, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(pList->list == NULL)
      {
	WlzFreeProperty(prop);
      }
      else
      {
	(void )WlzFreeProperty(prop);
	(void )WlzFreePropertyList(pList);
      }
    }
  }
  if(sGVWSpAry || cgtAry)
  {
    int       idT;

    if(sGVWSpAry)
    {
      for(idT = 0; idT < nThr; ++idT)
      {
        WlzGreyValueFreeWSp(sGVWSpAry[idT]);
      }
    }
    if(cgtAry)
    {
      int	idT;

      for(idT = 0; idT < nThr; ++idT)
      {
        AlgMatrixFree(cgtAry[idT]);
      }
    }
    AlcFree(cgtAry);
    AlcFree(sGVWSpAry);
  }
  /* Clear up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Smooths the (possibly) non-scalar features of the given
* 		object in place by applying a Gaussian filter with the
* 		given sigma values (sigma value <~ 0.0 implies no
* 		filtering in the corresponding direction.
* \param	obj		Given object with values to smooth.
* \param	smooth		Gaussian filter sigma values.
*/
WlzErrorNum			WlzTensorSmooth(
				  WlzObject *obj,
				  WlzDVertex3 smooth)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-06;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((obj->type != WLZ_2D_DOMAINOBJ) &&
          (obj->type != WLZ_3D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((smooth.vtX > eps) && (smooth.vtY > eps) && (smooth.vtZ > eps))
  {
    WlzGreyType	gType;


    gType = WlzGreyTypeFromObj(obj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      int	idV,
      		vpe = 1;

      if(WlzGreyTableIsTiled(obj->values.core->type))
      {
        vpe = obj->values.t->vpe;
      }
      for(idV = 0; idV < vpe; ++idV)
      {
	WlzObject *rObj0 = NULL,
		  *rObj1 = NULL;
	WlzIVertex3 order = {0, 0, 0},
		    direction = {1, 1, 1};

        rObj0 = WlzTensorGetComponent(obj, idV, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  rObj1 = WlzGaussFilter(rObj0, smooth, order, direction, gType,
				 ALG_PAD_END, 0.0, 0, &errNum);
	}
	(void )WlzFreeObj(rObj0);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzTensorSetComponent(obj, rObj1, idV);
	}
	(void )WlzFreeObj(rObj1);
      }
    }
  }
  return(errNum);
}

/*!
* \return	New object containing a single component of the given object.
* \ingroup	WlzFeatures
* \brief	Extracts a single value component from a (possibly) non-scalar
* 		object and creates an object with a non-tiled value table
* 		which is returned. This function may also be used to convert
* 		a tiled value object to a non-tiled value object.
* \param	tObj		Given object (possibly) with non-scalar values.
* \param	cpt		The component index which must be in the range
* 				[0-(v - 1)] where v is the values per element.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzObject			*WlzTensorGetComponent(
				  WlzObject *tObj,
				  int cpt,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzPixelV 	bgdV = {0};
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzGreyTableType gTType = 0;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(tObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((tObj->type != WLZ_2D_DOMAINOBJ) &&
          (tObj->type != WLZ_3D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(tObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(tObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    int		tIsTiled;

    tIsTiled = WlzGreyTableIsTiled(tObj->values.core->type);
    if((cpt < 0) || ((cpt > 0) && (!tIsTiled || (cpt > tObj->values.t->vpe))))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bgdV = WlzGetBackground(tObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTableTypeToGreyType(tObj->values.core->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gTType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzNewObjectValues(tObj, gTType, bgdV, 0, bgdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(tObj->type == WLZ_2D_DOMAINOBJ)
    {
      errNum = WlzTensorGetComponentValues2D(rObj, tObj, cpt, 0);
    }
    else /* tObj->type == WLZ_3D_DOMAINOBJ */
    {
      errNum = WlzTensorComponentValues3D(rObj, tObj, cpt, 0);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Sets a single value component in a (possibly) non-scalar
* 		object using scalar values of the second given object.
* \param	tObj		Given object (possibly) with non-scalar values
* 				which is to have component values set.
* \param	fObj            Given scalar valued object with values to
* 				be transfered to the first given object.
* \param	cpt		The component index which must be in the range
* 				[0-(v - 1)] where v is the values per element
* 				of the first object.
*/
WlzErrorNum			WlzTensorSetComponent(
				  WlzObject *tObj,
				  WlzObject *fObj,
				  int cpt)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((tObj == NULL) || (fObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((tObj->type != fObj->type) ||
          ((tObj->type != WLZ_2D_DOMAINOBJ) &&
           (tObj->type != WLZ_3D_DOMAINOBJ)))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((tObj->domain.core == NULL) ||
          (fObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((tObj->values.core == NULL) ||
          (fObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    int		tIsTiled;

    tIsTiled = WlzGreyTableIsTiled(tObj->values.core->type);
    if((cpt < 0) || ((cpt > 0) && (!tIsTiled || (cpt > tObj->values.t->vpe))))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(tObj->type == WLZ_2D_DOMAINOBJ)
    {
      errNum = WlzTensorSetComponentValues2D(fObj, tObj, cpt, 0);
    }
    else /* tObj->type == WLZ_3D_DOMAINOBJ */
    {
      errNum = WlzTensorComponentValues3D(fObj, tObj, cpt, 1);
    }
  }
  return(errNum);
}

/*!
* \return	Conforming mesh object with tensor values attached to
* 		elements or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor for each of it's valid
* 		elements.
* 		Given displacement \f$\vec{u}(\vec{r})\f$ with position
* 		vector \f$\vec{r}\f$ which maps a point from a space
* 		\f$\vec{x}\f$ to a space \f$\vec{u}\f$ the displacement
* 		gradient tensor is defined in 3D as
  		\f[
		{
		\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}
		u_{i,j} =
		\left [
		\begin{array}{ccc}
		  \pd{u_0}{x_0} & \pd{u_0}{x_1} & \pd{u_0}{x_2} \\
		  \pd{u_1}{x_0} & \pd{u_1}{x_1} & \pd{u_1}{x_2} \\
		  \pd{u_2}{x_0} & \pd{u_2}{x_1} & \pd{u_2}{x_2}
		\end{array}
		\right]
		}
  		\f]
*		with
  		\f[
		\Delta r_i = u_{ij} r_j
  		\f]
*		where \f$\vec{u} = \left[u_0, u_1, u_2\right]^T\f$
*		and \f$\Delta \vec{r} = \left[x_0, x_1, x_2\right]^T\f$.
*		The displacement gradient tensor matrix is just the
*		rotation and independent scaling part of the  affine
*		transform that displaces the element.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDGTensor(WlzObject *cObj, int invert,
				  WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	tObj = WlzCMeshDGTensor3D(cObj, invert, &errNum);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Points object with tensor values at the point locations
* 		or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		The tensor values at the sample points are computed at
* 		each point by computing an iverse distance weighted
* 		least squares general affine transform for the ring of
* 		nodes surrounding the closest node.
* 		See WlzCMeshDGTensor() for the description of the tensor.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDGTensorAtPts(WlzObject *cObj, int invert,
				       WlzDVertex3 sd, int dither,
				       WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((sd.vtX < WLZ_MESH_TOLERANCE) ||
          (sd.vtY < WLZ_MESH_TOLERANCE) ||
          (sd.vtZ < WLZ_MESH_TOLERANCE))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	tObj = WlzCMeshDGTensorAtPts3D(cObj, invert, sd, dither, &errNum);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Conforming mesh object with tensor values attached to
* 		elements or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the strain tensor for each of it's valid elements.
* 		This function uses WlzCMeshDGTensor() to compute the
* 		displacement gradient tensor and the derives the strain
* 		tensor \f$e_{ij}\f$ from this using:
		\f[
		e_{ij} = \frac{1}{2} (u_{ij} + u_{ji})
		\f]
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshStrainTensor(WlzObject *cObj, int invert,
 				      WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum  = WLZ_ERR_NONE;

  tObj = WlzCMeshDGTensor(cObj, invert, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	WlzCMeshDGToStrainTensor3D(tObj);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Points object with tensor values or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		The tensor values at the sample points are computed using
* 		WlzCMeshDGTensorAtPts(). The strain tensor is then
* 		computed from the displacement gradient tensor as for
* 		WlzCMeshStrainTensor().
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshStrainTensorAtPts(WlzObject *cObj, int invert,
					   WlzDVertex3 sd,
					   int dither,
					   WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum  = WLZ_ERR_NONE;

  tObj = WlzCMeshDGTensorAtPts(cObj, invert, sd, dither, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	WlzPtsDGToStrainTensor3D(tObj);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Extracts the required component value from a (possibly)
* 		tensor valued (possibly) 2 or 3D tiled value object into
* 		a non-tiled value object. This function assumes all given
* 		parameters are valid and that the domains of the two
* 		objects are identical.
* \param	rObj		Non-tiled value 2D object, already allocated
* 				with appropriate grey value type and background
* 				set.
* \param	tObj		Given (possibly) tensor valued source object.
* 				It is acceptable for this object to either have
* 				either 2D values of a 3D tiled value table.
* \param	cpt		Component index, known to be in range.
* \param	pln		If the source object is 3D then this is
* 				the current plane value, otherwise it is
* 				not used.
*/
static WlzErrorNum		WlzTensorGetComponentValues2D(
				  WlzObject *rObj,
				  WlzObject *tObj,
				  int cpt,
				  int pln)
{
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	    gWSp;
  WlzGreyValueWSpace *tGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((tGVWSp = WlzGreyValueMakeWSp(tObj, &errNum)) != NULL) &&
     ((errNum = WlzInitGreyScan(rObj, &iWSp, &gWSp)) == WLZ_ERR_NONE))
  {
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
    
    {
      int	idK;
      WlzGreyP	rP;

      rP = gWSp.u_grintptr;
      for(idK = iWSp.lftpos; idK <= iWSp.rgtpos; ++idK)
      {
	int	 idI;
	WlzGreyP tP;

        idI = idK - iWSp.lftpos;
        WlzGreyValueGet(tGVWSp, pln, iWSp.linpos, idK);
	tP = tGVWSp->gPtr[0];
	switch(gWSp.pixeltype)
	{
	  case WLZ_GREY_LONG:
	    rP.lnp[idI] = tP.lnp[cpt];
	    break;
	  case WLZ_GREY_INT:
	    rP.inp[idI] = tP.inp[cpt];
	    break;
	  case WLZ_GREY_SHORT:
	    rP.shp[idI] = tP.shp[cpt];
	    break;
	  case WLZ_GREY_UBYTE:
	    rP.ubp[idI] = tP.ubp[cpt];
	    break;
	  case WLZ_GREY_FLOAT:
	    rP.flp[idI] = tP.flp[cpt];
	    break;
	  case WLZ_GREY_DOUBLE:
	    rP.dbp[idI] = tP.dbp[cpt];
	    break;
	  case WLZ_GREY_RGBA:
	    rP.rgbp[idI] = tP.rgbp[cpt];
	    break;
	  default:
	    break;
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WlzGreyValueFreeWSp(tGVWSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Sets values of the required component in a (possibly)
* 		tensor valued (possibly) 2 or 3D tiled value object
* 		from a non-tiled value object. This function assumes that
* 		all the given parameters are valid and that the domains
* 		of the two objects are identical.
* \param	rObj		Non-tiled value 2D object with values to be
* 				transfered to the (possibly) tensor valued
* 				(possibly) tiled object.
* \param	tObj		Given (possibly) tensor valued source object.
* 				It is acceptable for this object to either have
* 				either 2D values of a 3D tiled value table.
* \param	cpt		Component index, known to be in range.
* \param	pln		If the source object is 3D then this is
* 				the current plane value, otherwise it is
* 				not used.
*/
static WlzErrorNum		WlzTensorSetComponentValues2D(
				  WlzObject *rObj,
				  WlzObject *tObj,
				  int cpt,
				  int pln)
{
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	    gWSp;
  WlzGreyValueWSpace *tGVWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((tGVWSp = WlzGreyValueMakeWSp(tObj, &errNum)) != NULL) &&
     ((errNum = WlzInitGreyScan(rObj, &iWSp, &gWSp)) == WLZ_ERR_NONE))
  {
    WlzGreyType rGType,
		tGType;

    rGType = gWSp.pixeltype;
    tGType = tGVWSp->gType;
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
    
    {
      int	idK;
      WlzGreyP	rP;

      rP = gWSp.u_grintptr;
      for(idK = iWSp.lftpos; idK <= iWSp.rgtpos; ++idK)
      {
	int	 idI;
	WlzGreyP tP;

        idI = idK - iWSp.lftpos;
        WlzGreyValueGet(tGVWSp, pln, iWSp.linpos, idK);
	tP = tGVWSp->gPtr[0];
	WlzValueCopyGreyToGrey(tP, cpt, tGType, rP, idI, rGType, 1);
      }
    }
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WlzGreyValueFreeWSp(tGVWSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Extracts or sets the required component values from/to a
* 		(possibly) tensor valued 3D (possibly) tiled value object
* 		into/from a non-tiled value object. This function assumes
* 		that all the given parameters are valid and that the domains
* 		of the two objects are identical.
* \param	rObj		Non-tiled value 3D object, already allocated
* 				with appropriate grey value type and background
* 				set.
* \param	tObj		Given (possibly) tensor valued object. It
* 				is acceptable for this object to either have
* 				either a 3D voxel or 3D tiled value table.
* \param	cpt		Component index, known to be in range.
* \param	set		Set tiled values rather than get them if
* 				non-zero.
*/
static WlzErrorNum		WlzTensorComponentValues3D(
				  WlzObject *rObj,
				  WlzObject *tObj,
				  int cpt,
				  int set)
{
  int		pln,
  		tiled;
  WlzPlaneDomain *pDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  pDom = rObj->domain.p;
  tiled = WlzGreyTableIsTiled(tObj->values.core->type);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(pln = pDom->plane1; pln <= pDom->lastpl; ++pln)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      int	idP;
      WlzDomain	*dom2;
      WlzValues	*rVal2;
      WlzObject	*rObj2 = NULL,
		*tObj2 = NULL;
      WlzErrorNum errNum2 = WLZ_ERR_NONE;

      idP = pln - pDom->plane1;
      dom2 = pDom->domains + idP;
      rVal2 = rObj->values.vox->values + idP;
      rObj2 = WlzAssignObject(
	      WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2, *rVal2,
			  NULL, NULL, &errNum2), NULL);
      if(errNum2 == WLZ_ERR_NONE)
      {
	if(tiled)
	{
	  tObj2 = WlzAssignObject(
		  WlzMakeMain(WLZ_3D_DOMAINOBJ, tObj->domain, tObj->values,
			      NULL, NULL, &errNum2), NULL);
	}
	else
	{
	  WlzValues	*tVal2;

	  tVal2 = tObj->values.vox->values;
	  tObj2 = WlzAssignObject(
		  WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2, *tVal2,
			      NULL, NULL, &errNum2), NULL);
	}
      }
      if(errNum2 == WLZ_ERR_NONE)
      {
	errNum2 = (set)?
	          WlzTensorSetComponentValues2D(rObj2, tObj2, cpt, pln):
		  WlzTensorGetComponentValues2D(rObj2, tObj2, cpt, pln);
      }
      (void )WlzFreeObj(rObj2);
      (void )WlzFreeObj(tObj2);
      if(errNum2 != WLZ_ERR_NONE)
      {
#ifdef _OPENMP
#pragma omp critical (WlzTensorComponentValues3D)
#endif
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = errNum2;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \ingroup	WlzFeatures
* \brief	Computes the Right Cauchy-Green tensor from the deformation
* 		tensor.
* \param	cg			Matrix in which to place the Right
* 					Cauchy-Green tensor values.
* \param	gP 			Grey pointer for the deformation tensor
* 					values.
* \param	gT			Grey type o the values.
*/
static void			WlzDGTensorToRCauchyGreen(
				  AlgMatrix cg,
				  WlzGreyP gP,
				  WlzGreyType gT)
{
  WlzGreyP	m;
  double	f[9];

  /* Get deformation gradient matrix values into f. */
  m.dbp = f;
  WlzValueCopyGreyToGrey(m, 0, WLZ_GREY_DOUBLE, gP, 0, gT, 9);
  /* Use deformation gradient values to compute the right Cauchy-Green
   * tensor ( C = F^T F ). */
  cg.rect->array[0][0] = f[0] * f[0] + f[3] * f[3] + f[6] * f[6];
  cg.rect->array[0][1] = f[0] * f[1] + f[3] * f[4] + f[6] * f[7];
  cg.rect->array[0][2] = f[0] * f[2] + f[3] * f[5] + f[6] * f[8];
  cg.rect->array[1][0] = f[1] * f[0] + f[4] * f[3] + f[7] * f[6];
  cg.rect->array[1][1] = f[1] * f[1] + f[4] * f[4] + f[7] * f[7];
  cg.rect->array[1][2] = f[1] * f[2] + f[4] * f[5] + f[7] * f[8];
  cg.rect->array[2][0] = f[2] * f[0] + f[5] * f[3] + f[8] * f[6];
  cg.rect->array[2][1] = f[2] * f[1] + f[5] * f[4] + f[8] * f[7];
  cg.rect->array[2][2] = f[2] * f[2] + f[5] * f[5] + f[8] * f[8];
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the determinant of the Jacobian matrix values in
* 		the source, the determinants are placed in the destination.
* 		The returned Woolz error code will be WLZ_ERR_NONE unless
* 		the grey type is not one of WLZ_GREY_INT, WLZ_GREY_SHORT,
* 		WLZ_GREY_UBYTE, WLZ_GREY_FLOAT or  WLZ_GREY_DOUBLE.
* \param	dGVWSp			Destination grey value workspace.
* \param	sGVWSp			Source grey value workspace.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	pln			The current plane.
* \param	lin			The current line.
* \param        lft			Left start of interval.
* \param	rgt			Right end of interval.
*/
static void			WlzDGTensorSDFeatDetJacItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt)
{
  int		idK;

  for(idK = lft; idK <= rgt; ++idK)
  {
    double	v = 0.0;
    WlzGreyP	d;

    WlzGreyValueGet(dGVWSp, pln, lin, idK);
    WlzGreyValueGet(sGVWSp, pln, lin, idK);
    WlzDGTensorToRCauchyGreen(cgt, sGVWSp->gPtr[0], sGVWSp->gType);
    (void )AlgMatrixLUDetermRaw(cgt.rect->array, 3, &v);
    d = dGVWSp->gPtr[0];
    v = (v < 0.0)? 0.0: sqrt(v);
    WLZ_CLAMP_DOUBLE_TO_GREYP(d, 0, v, dGVWSp->gType);
  }
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the determinant of the Jacobian matrix value in
* 		the source, the determinant is placed in the destination.
* 		The returned Woolz error code will be WLZ_ERR_NONE unless
* 		the grey type is not one of WLZ_GREY_INT, WLZ_GREY_SHORT,
* 		WLZ_GREY_UBYTE, WLZ_GREY_FLOAT or  WLZ_GREY_DOUBLE.
* \param	dPV			Destination point value.
* \param	sGP			Deformation gradient values.
* \param	dGType			Destination grey type.
* \param	sGType			Source grey type.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	idP			Current point index.
*/
static void			WlzDGTensorSDFeatDetJacPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP)
{
  double	v = 0.0;

  WlzDGTensorToRCauchyGreen(cgt, sGP, sGType);
  (void )AlgMatrixLUDetermRaw(cgt.rect->array, 3, &v);
  v = (v < 0.0)? 0.0: sqrt(v);
  WLZ_CLAMP_DOUBLE_TO_GREYP(dPV->values, idP, v, dGType);
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the principle eigen vectors of the Jacobian
* 		matrix values in the source, the vectors are placed in
* 		the destination as a matrix of three column vectors each
* 		with three values.
* \param	dGVWSp			Destination grey value workspace.
* \param	sGVWSp			Source grey value workspace.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	pln			The current plane.
* \param	lin			The current line.
* \param        lft			Left start of interval.
* \param	rgt			Right end of interval.
*/
static void			WlzDGTensorSDFeatEigenVecItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt)
{
  int		idK;
  double	v[3];

  for(idK = lft; idK <= rgt; ++idK)
  {
    int		i,
		j,
		k;
    WlzGreyP	d;

    WlzGreyValueGet(sGVWSp, pln, lin, idK);
    WlzGreyValueGet(dGVWSp, pln, lin, idK);
    d = dGVWSp->gPtr[0];
    WlzDGTensorToRCauchyGreen(cgt, sGVWSp->gPtr[0], sGVWSp->gType);
    (void )AlgMatrixRSEigen(cgt, v, 1);
    switch(dGVWSp->gType)
    {
      case WLZ_GREY_FLOAT:
	for(j = k = 0; j < 3; ++j)
	{
	  for(i = 0; i < 3; ++i)
	  {
	    d.flp[k++] = WLZ_CLAMP(cgt.rect->array[j][i], -FLT_MAX, FLT_MAX);
	  }
	}
	break;
      case WLZ_GREY_DOUBLE:
	for(j = k = 0; j < 3; ++j)
	{
	  for(i = 0; i < 3; ++i)
	  {
	    d.dbp[k++] = cgt.rect->array[j][i];
	  }
	}
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the principle eigen vectors of the Jacobian
* 		matrix values in the source, the vectors are placed in
* 		the destination as a matrix of three column vectors each
* 		with three values.
* \param	dPV			Destination point value.
* \param	sGP			Deformation gradient values.
* \param	dGType			Destination grey type.
* \param	sGType			Source grey type.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	idP			Current point index.
*/
static void			WlzDGTensorSDFeatEigenVecPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP)
{
  int		i,
	      	j,
	      	k;
  double	v[3];
  WlzGreyP	d;

  WlzDGTensorToRCauchyGreen(cgt, sGP, sGType);
  (void )AlgMatrixRSEigen(cgt, v, 1);
  switch(dPV->vType)
  {
    case WLZ_GREY_FLOAT:
      d.flp = dPV->values.flp + (9 * idP);
      for(j = k = 0; j < 3; ++j)
      {
	for(i = 0; i < 3; ++i)
	{
	  d.flp[k++] = WLZ_CLAMP(cgt.rect->array[j][i], -FLT_MAX, FLT_MAX);
	}
      }
      break;
    case WLZ_GREY_DOUBLE:
      d.dbp = dPV->values.dbp + (9 * idP);
      for(j = k = 0; j < 3; ++j)
      {
	for(i = 0; i < 3; ++i)
	{
	  d.dbp[k++] = cgt.rect->array[j][i];
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the eigen values of the Jacobian matrix values in
* 		the source, the eigen values are placed in the destination
* 		as a vector of three values.
* 		Integral values (int, short or ubyte) are not allowed and are
* 		considered an error.
* \param	dGVWSp			Destination grey value workspace.
* \param	sGVWSp			Source grey value workspace.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	pln			The current plane.
* \param	lin			The current line.
* \param        lft			Left start of interval.
* \param	rgt			Right end of interval.
*/
static void			WlzDGTensorSDFeatEigenValItv(
				  WlzGreyValueWSpace *dGVWSp,
				  WlzGreyValueWSpace *sGVWSp,
				  AlgMatrix cgt,
				  int pln,
				  int lin,
				  int lft,
				  int rgt)
{
  int		idK;
  double	v[3];

  for(idK = lft; idK <= rgt; ++idK)
  {
    WlzGreyP	d;

    WlzGreyValueGet(sGVWSp, pln, lin, idK);
    WlzGreyValueGet(dGVWSp, pln, lin, idK);
    d = dGVWSp->gPtr[0];
    WlzDGTensorToRCauchyGreen(cgt, sGVWSp->gPtr[0], sGVWSp->gType);
    (void )AlgMatrixRSEigen(cgt, v, 0);
    switch(dGVWSp->gType)
    {
      case WLZ_GREY_FLOAT:
	d.flp[0] = WLZ_CLAMP(v[0], -FLT_MAX, FLT_MAX);
	d.flp[1] = WLZ_CLAMP(v[1], -FLT_MAX, FLT_MAX);
	d.flp[2] = WLZ_CLAMP(v[2], -FLT_MAX, FLT_MAX);
	break;
      case WLZ_GREY_DOUBLE:
	d.dbp[0] = v[0];
	d.dbp[1] = v[1];
	d.dbp[2] = v[2];
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup      WlzFeatures
* \brief	Computes the eigen values of the Jacobian matrix values in
* 		the source, the eigen values are placed in the destination
* 		as a vector of three values.
* 		Integral values (int, short or ubyte) are not allowed and are
* 		considered an error.
* \param	dPV			Destination point value.
* \param	sGP			Deformation gradient values.
* \param	dGType			Destination grey type.
* \param	sGType			Source grey type.
* \param	cgt			Working matrix for right Cauchy-Green
* 					tensor.
* \param	idP			Current point index.
*/
static void			WlzDGTensorSDFeatEigenValPnt(
				  WlzPointValues *dPV,
				  WlzGreyP sGP,
				  WlzGreyType dGType,
				  WlzGreyType sGType,
				  AlgMatrix cgt,
				  int idP)
{
  double	v[3];
  WlzGreyP	d;

  WlzDGTensorToRCauchyGreen(cgt, sGP, sGType);
  (void )AlgMatrixRSEigen(cgt, v, 1);
  switch(dPV->vType)
  {
    case WLZ_GREY_FLOAT:
      d.flp = dPV->values.flp + (3 * idP);
      d.flp[0] = WLZ_CLAMP(v[0], -FLT_MAX, FLT_MAX);
      d.flp[1] = WLZ_CLAMP(v[1], -FLT_MAX, FLT_MAX);
      d.flp[2] = WLZ_CLAMP(v[2], -FLT_MAX, FLT_MAX);
      break;
    case WLZ_GREY_DOUBLE:
      d.dbp = dPV->values.dbp + (3 * idP);
      d.dbp[0] = v[0];
      d.dbp[1] = v[1];
      d.dbp[2] = v[2];
      break;
    default:
      break;
  }
}

/*!
* \return	New Woolz point object or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform this function computes
* 		the displacement gradient tensor for each of it's valid
* 		elements. See WlzCMeshDGTensor().
* \param	cObj			Given 3D conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshDGTensor3D(WlzObject *cObj, int invert,
				     WlzErrorNum *dstErr)
{
  WlzObject	*iObj = NULL,
  		*tObj = NULL;
  WlzCMesh3D	*mesh = NULL;
  WlzIndexedValues *cIxv = NULL,
  		   *tIxv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mesh = cObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((cIxv = cObj->values.x)->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((cIxv->rank < 1) || (cIxv->dim[0] < 3) ||
          (cIxv->vType != WLZ_GREY_DOUBLE) ||
	  (cIxv->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    if(invert)
    {
      iObj = WlzCMeshTransformInvert(cObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	cObj = iObj;
        mesh = cObj->domain.cm3;
	cIxv = cObj->values.x;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim[2];

    dim[0] = dim[1] = 3;
    tIxv = WlzMakeIndexedValues(cObj, 2, dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_ELM, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues	tVal;

      tVal.x = tIxv;
      tObj = WlzMakeMain(WLZ_CMESH_3D, cObj->domain, tVal, NULL, NULL,
                         &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        (void )WlzFreeIndexedValues(tIxv);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idE,
    		maxElm;
    AlcVector	*elmVec;

    elmVec = mesh->res.elm.vec;
    maxElm = mesh->res.elm.maxEnt;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idE = 0; idE < maxElm; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
      if(elm->idx >= 0)
      {
	double	*ten;

	ten = (double *)WlzIndexedValueGet(tIxv, elm->idx);
	WlzCMeshElmSetDGTensor3D(elm, 0, cIxv, ten);
      }
    }
  }
  if(invert)
  {
    WlzFreeObj(iObj);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Points object with tensor values at the point locations
* 		or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGTensorAtPts() with all checks having been done.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshDGTensorAtPts3D(WlzObject *cObj, int invert,
					  WlzDVertex3 sd, int dither,
					  WlzErrorNum *dstErr)
{
  WlzDomain	dom;
  WlzValues	val;
  WlzDVertex3	sOrg;
  WlzIVertex3	sNum;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if((mesh = cObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv = cObj->values.x)->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv->rank < 1) || (ixv->dim[0] < 3) ||
          (ixv->vType != WLZ_GREY_DOUBLE) ||
	  (ixv->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    WlzCMeshUpdateBBox3D(mesh);
    sOrg.vtX = ceil(mesh->bBox.xMin / sd.vtX) * sd.vtX;
    sOrg.vtY = ceil(mesh->bBox.yMin / sd.vtY) * sd.vtY;
    sOrg.vtZ = ceil(mesh->bBox.zMin / sd.vtZ) * sd.vtZ;
    sNum.vtX = floor((mesh->bBox.xMax - sOrg.vtX) / sd.vtX) + 1;
    sNum.vtY = floor((mesh->bBox.yMax - sOrg.vtY) / sd.vtY) + 1;
    sNum.vtZ = floor((mesh->bBox.zMax - sOrg.vtZ) / sd.vtZ) + 1;
    if((sNum.vtX <= 0) || (sNum.vtY <= 0) || (sNum.vtZ <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      int	maxPts;
      WlzVertexP dum;

      dum.v = NULL;
      maxPts = sNum.vtX * sNum.vtY * sNum.vtZ;
      dom.pts = WlzMakePoints(WLZ_POINTS_3D, 0, dum, maxPts, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int	dim[2];

	dim[0] = dim[1] = 3;
        val.pts = WlzMakePointValues(maxPts, 2, dim, WLZ_GREY_DOUBLE,
	                             &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj = WlzMakeMain(WLZ_POINTS, dom, val, NULL, NULL, &errNum);
      }
      else
      {
        (void )WlzFreeDomain(dom);
        (void )WlzFreePointValues(val.pts);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		pIdx;
    WlzIVertex3 sIdx;
    WlzDVertex3 sPos;

    pIdx = 0;
    AlgRandSeed(0);
    for(sIdx.vtZ = 0; sIdx.vtZ <sNum.vtZ; ++(sIdx.vtZ))
    {
      sPos.vtZ = sOrg.vtZ + (sIdx.vtZ * sd.vtZ);
      for(sIdx.vtY = 0; sIdx.vtY <sNum.vtY; ++(sIdx.vtY))
      {
        sPos.vtY = sOrg.vtY + (sIdx.vtY * sd.vtY);
	for(sIdx.vtX = 0; sIdx.vtX <sNum.vtX; ++(sIdx.vtX))
	{
	  int	eIdx;
	  WlzDVertex3 dPos;

          sPos.vtX = sOrg.vtX + (sIdx.vtX * sd.vtX);
	  dPos = sPos;
	  if(dither)
	  {
	    dPos.vtX = AlgRandZigNormal(sPos.vtX, sd.vtX * 0.2);
	    dPos.vtY = AlgRandZigNormal(sPos.vtY, sd.vtY * 0.2);
	    dPos.vtZ = AlgRandZigNormal(sPos.vtZ, sd.vtZ * 0.2);
	  }
	  /* If dPos is in an element of the mesh compute it's displacement
	   * gradient tensor. */
          eIdx = WlzCMeshElmEnclosingPos3D(mesh, -1,
	                                   dPos.vtX, dPos.vtY, dPos.vtZ,
					   0, NULL);
	  if(eIdx >= 0)
	  {
	    double	*ten;
	    WlzCMeshElm3D *elm;

	    elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, eIdx);
	    if(elm->idx >= 0)
	    {
	      dom.pts->points.d3[pIdx] = dPos;
	      ten = WlzPointValueGet(val.pts, pIdx);
	      WlzCMeshElmSetDGTensor3D(elm, invert, ixv, ten);
	      ++pIdx;
	    }
	  }
	}
      }
    }
    dom.pts->nPoints = pIdx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform with displacement gradient
* 		tensor values this function converts these values to strain
* 		tensor values.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGToStrainTensor() with all checks having been done.
* \param	cObj			Given 3D conforming mesh object.
*/
static void 	WlzCMeshDGToStrainTensor3D(WlzObject *tObj)
{
  int 		idE,
    		maxElm;
  AlcVector	*elmVec;
  WlzCMesh3D	*mesh = NULL;
  WlzIndexedValues *tIxv = NULL;

  mesh = tObj->domain.cm3;
  tIxv = tObj->values.x;
  elmVec = mesh->res.elm.vec;
  maxElm = mesh->res.elm.maxEnt;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(idE = 0; idE < maxElm; ++idE)
  {
    WlzCMeshElm3D *elm;

    elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
    if(elm->idx >= 0)
    {
      double	*ten;

      ten = (double *)WlzIndexedValueGet(tIxv, elm->idx);
      ten[1] = ten[3] = 0.5 * (ten[1] + ten[3]);
      ten[2] = ten[6] = 0.5 * (ten[2] + ten[6]);
      ten[5] = ten[7] = 0.5 * (ten[5] + ten[7]);
    }
  }
}

/*!
* \ingroup 	WlzFeatures
* \brief	Given a points object with 3D vertices and displacement
* 		gradient tensor values this function converts these values
* 		to strain tensor values.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGToStrainTensor() with all checks having been done.
* \param	cObj			Given 3D conforming mesh object.
*/
static void 	WlzPtsDGToStrainTensor3D(WlzObject *tObj)
{
  int		idx,
  		nPts;
  WlzPoints	*pd;
  WlzPointValues *pv;

  pd = tObj->domain.pts;
  pv = tObj->values.pts;
  nPts = pd->nPoints;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(idx = 0; idx < nPts; ++idx)
  {
    double	*ten;

    ten = (double *)WlzPointValueGet(pv, idx);
    ten[1] = ten[3] = 0.5 * (ten[1] + ten[3]);
    ten[2] = ten[6] = 0.5 * (ten[2] + ten[6]);
    ten[5] = ten[7] = 0.5 * (ten[5] + ten[7]);
  }
}

/*!
* \ingroup	WlzFeatures
* \brief	Given a 3D conforming mesh node, an indexed values with
* 		displacements for the mesh and a tensor pointer; this
* 		function computes the displacement gradient tensor
* 		value.
* \param	elm			Given mesh element.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	cIxv			Mesh displacment indexed values.
* \param	ten			Tensor pointer.
*/
static void	WlzCMeshElmSetDGTensor3D(WlzCMeshElm3D *elm,
					 int invert,
					 WlzIndexedValues *cIxv,
					 double *ten)
{
  int		idN;
  double	tr[16];
  WlzCMeshNod3D *nod[4];
  WlzDVertex3 sVx[4],
	      dVx[4];

  nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
  nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
  nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
  nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
  for(idN = 0; idN < 4; ++idN)
  {
    double	*dsp;

    sVx[idN] = nod[idN]->pos;
    dsp = (double *)WlzIndexedValueGet(cIxv, nod[idN]->idx);
    dVx[idN].vtX = sVx[idN].vtX + dsp[0];
    dVx[idN].vtY = sVx[idN].vtY + dsp[1];
    dVx[idN].vtZ = sVx[idN].vtZ + dsp[2];
  }
  if(invert)
  {
    WlzGeomTetraAffineSolve(tr, dVx, sVx, WLZ_MESH_TOLERANCE_SQ);
  }
  else
  {
    WlzGeomTetraAffineSolve(tr, sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
  }
  ten[0] = tr[0]; ten[1] = tr[1]; ten[2] = tr[2]; 
  ten[3] = tr[4]; ten[4] = tr[5]; ten[5] = tr[6]; 
  ten[6] = tr[8]; ten[7] = tr[9]; ten[8] = tr[10]; 
}
