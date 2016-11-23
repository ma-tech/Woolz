#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCMeshTransform.c
* \author       Bill Hill
* \date         October 2004
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
* \brief	Functions for creating and applying 2D and 3D conforming
* 		mesh transforms.
* \ingroup	WlzTransform
*/

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Wlz.h>

#define WLZ_CMESH_POS_DTOI(X) ((int )floor(X))

/*!
* \enum		_WlzCMeshScanElmFlags
* \ingroup	WlzTransform
* \brief	Flags for the conforming mesh scannning elements.
*/
typedef enum _WlzCMeshScanElmFlags
{
  WLZ_CMESH_SCANELM_NONE 	= (0),	/*! Clear - no flags set for the
  					    element. */
  WLZ_CMESH_SCANELM_SQUASH	= (1),  /*! The transformed element has
					    negligable area/volume and
					    only the transform coefficients
					    for translation are non zero. */
  WLZ_CMESH_SCANELM_FWD 	= (2),	/*! The transform is forward, ie
  					    transforming source to
					    destination. */
  WLZ_CMESH_SCANELM_REV 	= (4)	/*! The transform is reverse, ie
  					    transforming destination to
					    source. */
} WlzCMeshScanElmFlags;

/*!
* \struct	_WlzCMeshScanElm2D
* \ingroup	WlzTransform
* \brief	Conforming mesh scanning element for 2D mesh.
*/
typedef struct _WlzCMeshScanElm2D
{
  int		idx;		/*! Index of the current mesh element. */
  unsigned	flags;		/*! Flags set using ::WlzCMeshScanElmFlags. */
  double	trX[3];		/*! Within element affine transform
  				    coefficients for the x component
				    \f$
  				    t_x = trX_0 s_x + trX_1 s_y + trX_2
				    \f$ */
  double	trY[3];		/*! Within element affine transform
  				    coefficients for the y component
				    \f$
  				    t_y = trY_0 s_x + trY_1 s_y + trY_2
				    \f$ */
} WlzCMeshScanElm2D;

/*!
* \struct	_WlzCMeshScanElm3D
* \ingroup	WlzTransform
* \brief	Conforming mesh scanning element for 2D5 or 3D mesh.
*/
typedef struct _WlzCMeshScanElm3D
{
  int		idx;		/*! Index of the current mesh element. */
  unsigned	flags;		/*! Flags set using ::WlzCMeshScanElmFlags. */
  double	tr[16];		/*! Within element affine transform
  				    coefficients with
			    \f{eqnarray*}
			    t_x = tr_0 s_x + tr_1 s_y + tr_2 s_z + tr_3 \\
			    t_y = tr_4 s_x + tr_5 s_y + tr_6 s_z + tr_7 \\
			    t_z = tr_8 s_x + tr_9 s_y + tr_{10} s_z + tr_{11}
			    \f} */
} WlzCMeshScanElm3D;

/*!
* \struct       _WlzCMeshScanItv2D
* \ingroup      WlzTransform
* \brief        Scan interval within a 2D conforming mesh element.
*/
typedef struct _WlzCMeshScanItv2D
{
  int           elmIdx;                 /*! Element index. */
  int           line;                   /*! Line of interval. */
  int           lftI;                   /*! Start of interval. */
  int           rgtI;                   /*! End of interval. */
} WlzCMeshScanItv2D;

/*!
* \struct       _WlzCMeshScanItv3D
* \ingroup      WlzTransform
* \brief        Scan interval within a 3D conforming mesh element.
*/
typedef struct _WlzCMeshScanItv3D
{
  int           elmIdx;                 /*! Element index. */
  int           line;                   /*! Line of interval. */
  int           plane;                   /*! Line of interval. */
  int           lftI;                   /*! Start of interval. */
  int           rgtI;                   /*! End of interval. */
} WlzCMeshScanItv3D;

/*!
* \struct       _WlzCMeshScanWSp2D
* \ingroup      WlzTransform
* \brief        Conforming mesh scanning workspace for a 2D mesh.
*/
typedef struct _WlzCMeshScanWSp2D
{
  WlzObject	*mTr;			/*! The conforming mesh transform. */
  int		nItvs;			/*! Number of element intervals. */
  WlzCMeshScanItv2D *itvs;		/*! Element intervals sorted by line
  					    then left column. */
  WlzCMeshScanElm2D *dElm;		/*! Destination mesh element data. */
} WlzCMeshScanWSp2D;

/*!
* \struct       _WlzCMeshScanWSp3D
* \ingroup      WlzTransform
* \brief        Conforming mesh scanning workspace for a 3D mesh.
*/
typedef struct _WlzCMeshScanWSp3D
{
  WlzObject	*mTr;			/*! The conforming mesh transform. */
  int		nItvs;			/*! Number of element intervals. */
  WlzCMeshScanItv3D *itvs;		/*! Element intervals sorted by line
  					    then left column. */
  WlzCMeshScanElm3D *dElm;		/*! Destination mesh element data. */
  WlzIBox3	dBox;			/*! Bounding box of the displaced
  					    mesh. */
} WlzCMeshScanWSp3D;

static void 			WlzCMeshUpdateScanElm2D(
				  WlzObject *mObj,
				  WlzCMeshScanElm2D *sElm,
				  int fwd);
static void 			WlzCMeshUpdateScanElm2D5(
				  WlzObject *mObj,
				  WlzCMeshScanElm3D *sElm,
				  int fwd);
static void 			WlzCMeshUpdateScanElm3D(
				  WlzObject *mObj,
				  WlzCMeshScanElm3D *sElm,
				  int fwd);
static void			WlzCMeshScanWSpFree2D(
				  WlzCMeshScanWSp2D *mSWSp);
static void			WlzCMeshScanWSpFree3D(
				  WlzCMeshScanWSp3D *mSWSp);
static void			WlzCMeshScanClearOlpBuf(
				  WlzGreyP olpBuf,
				  int *olpCnt,
  				  WlzGreyType gType,
				  int bufidth,
				  int clrWidth);
static void			WlzCMeshSqzRedundantItv3D(
				  WlzCMeshScanWSp3D *mSWSp);
static WlzErrorNum 		WlzCMeshInterpolateNod2DKrig(
				  WlzGreyP dst,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateElm2DNearest(
				  WlzGreyP dst,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateElm2DLinear(
				  WlzGreyP dst,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateElm3DLinear(
				  WlzGreyP dst,
				  int pl,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateNod2DNearest(
				  WlzGreyP dst,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateNod2DLinear(
				  WlzGreyP dst,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateNod3DNearest(
				  WlzGreyP dst,
				  int pl,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateNod3DLinear(
				  WlzGreyP dst,
				  int pl,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static WlzErrorNum 		WlzCMeshInterpolateElm3DNearest(
				  WlzGreyP dst,
				  int pl,
				  int ln,
				  int kolL,
				  int kolR,
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  int ixi);
static int			WlzCMeshScanTriElm2D(
				  WlzCMeshScanWSp2D *mSWSp,
				  int trans,
				  WlzCMeshElm2D *elm,
				  int iIdx);
static int			WlzCMeshItv2Cmp(
				  const void *cmp0,
				  const void *cmp1);
static int			WlzCMeshItv3Cmp(
				  const void *cmp0,
				  const void *cmp1);
static int			WlzCMeshDVertex3Cmp(
				  const void *cmp0,
				  const void *cmp1);
static WlzObject 		*WlzCMeshTransformInvert2D(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshTransformInvert2D5(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshTransformInvert3D(
				  WlzObject *gObj,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshProduct2D(
				  WlzObject *tr0,
				  WlzObject *tr1,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshProduct3D(
				  WlzObject *tr0,
				  WlzObject *tr1,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzCMeshTransformValues2D(
				  WlzObject *dstObj,
				  WlzObject *srcObj,
				  WlzObject *mObj,
				  WlzInterpolationType interp);
static WlzErrorNum 		WlzCMeshTetElmItv3D(
				  AlcVector *itvVec,
				  int *idI,
				  int elmIdx,
				  WlzDVertex3 *vtx);
static WlzErrorNum 		WlzCMeshTriElmItv3D(AlcVector *itvVec,
				  int *idI,
				  int elmIdx,
				  WlzDVertex3 *vtx);
static WlzErrorNum 		WlzCMeshScanMakeOlpBufs(
				  WlzObject *obj,
				  WlzGreyType gType,
				  WlzGreyP *dstOlpBuf,
				  int **dstOlpCnt,
				  int bufWidth);
static WlzErrorNum 		WlzCMeshAddItv3D(
				  AlcVector *itvVec,
				  int *idI,
				  int elmIdx,
				  WlzDVertex3 *vtx);
static WlzErrorNum 		WlzCMeshScanObjValues3D(
				  WlzObject *dstObj,
				  WlzObject *srcObj,
				  WlzCMeshScanWSp3D *mSWSp,
				  WlzInterpolationType interp);
static WlzErrorNum 		WlzCMeshScanFlushOlpBuf(
				  WlzGreyP dGP,
				  WlzGreyP olpBuf,
				  int *olpCnt,
				  int bufWidth,
				  WlzPixelV bgdV,
				  int iLft,
				  int iRgt,
				  WlzInterpolationType interp,
				  WlzGreyType gType);
static WlzErrorNum 		WlzCMeshAffineProduct2D(
				  WlzObject *trM,
				  WlzAffineTransform *trA,
				  int order);
static WlzErrorNum 		WlzCMeshAffineProduct2D5(
				  WlzObject *trM,
				  WlzAffineTransform *trA,
				  int order);
static WlzErrorNum 		WlzCMeshAffineProduct3D(
				  WlzObject *trM,
				  WlzAffineTransform *trA,
				  int order);
static WlzObject 		*WlzCMeshTransformObjPDomain3D(
				  WlzObject *srcObj,
				  WlzObject *mObj,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshTransformObjV3D(
				  WlzObject *srcObj,
				  WlzObject *mObj,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzContour		*WlzCMeshTransformContour(
				  WlzContour *srcCtr,
				  WlzObject *mObj,
				  int newModFlg,
				  WlzErrorNum *dstErr);
static WlzGMModel		*WlzCMeshTransformGMModel(
				  WlzGMModel *srcM,
				  WlzObject *mObj,
				  int newModFlg,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshScanObjPDomain3D(
				  WlzObject *srcObj,
				  WlzCMeshScanWSp3D *mSWSp,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshToDomObj2D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshToDomObj3D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshToDomObjValues2D(
				  WlzObject *dObj,
				  WlzObject *mObj,
                                  WlzInterpolationType itp,
				  int ixi,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshToDomObjValues3D(
				  WlzObject *dObj,
				  WlzObject *mObj,
                                  WlzInterpolationType itp,
				  int ixi,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshExpansion2D(
				  WlzObject *cObj,
				  int inverse,
				  int method,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshExpansion3D(
				  WlzObject *cObj,
				  int inverse,
				  int method,
				  WlzErrorNum *dstErr);
static WlzPolygonDomain 	*WlzCMeshTransformPoly(
				  WlzPolygonDomain *srcPoly,
				  WlzObject *mObj,
				  WlzErrorNum *dstErr);
static WlzBoundList 		*WlzCMeshTransformBoundList(
				  WlzBoundList *srcBound,
				  WlzObject *mObj,
				  WlzErrorNum *dstErr);
static WlzCMesh2D 		*WlzCMeshTransformCMesh2D(
				  WlzCMesh2D *sMesh,
				  WlzObject *mObj,
				  int newMesh,
				  WlzErrorNum *dstErr);
static WlzCMesh2D5 		*WlzCMeshTransformCMesh2D5(
				  WlzCMesh2D5 *sMesh,
				  WlzObject *mObj,
				  int newMesh,
				  WlzErrorNum *dstErr);
static WlzCMesh3D 		*WlzCMeshTransformCMesh3D(
				  WlzCMesh3D *sMesh,
				  WlzObject *mObj,
				  int newMesh,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp2D 	*WlzCMeshScanWSpInit2D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp3D 	*WlzCMeshMakeScanWSp3D(
				  WlzObject *mObj,
				  int nItv,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp3D 	*WlzCMeshScanWSpInit3D(
				  WlzObject *mObj,
				  int trans,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzScaleIndexedVal(
				  double scale,
                                  WlzIndexedValues *ixv,
                                  int idV,
				  WlzIndexedValues *ixcSrc);
static WlzErrorNum		WlzScaleCMeshValueNodOrElem(
                                  WlzObject *obj,
				  double scale,
				  WlzIndexedValues *ixcSrc);

#ifdef WLZ_CMESHTRANSFORM_DEBUG
static WlzErrorNum 		WlzCMeshVerifyWSp3D(
				  WlzObject *srcObj,
				  WlzCMeshScanWSp3D *mSWSp);
#endif

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Inverts the given constrained mesh transform in place.
* \param	gObj			Given constrained mesh object
* 					with indexed values for transform.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshTransformInvert(WlzObject *gObj, WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(gObj->values.core->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_CMESH_2D:
	rObj = WlzCMeshTransformInvert2D(gObj, &errNum);
        break;
      case WLZ_CMESH_2D5:
	rObj = WlzCMeshTransformInvert2D5(gObj, &errNum);
        break;
      case WLZ_CMESH_3D:
	rObj = WlzCMeshTransformInvert3D(gObj, &errNum);
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
  return(rObj);
}

/*!
* \return	Inverted constrained mesh transform.
* \ingroup	WlzTransform
* \brief	Inverts the given 2D constrained mesh transform. Deleted
* 		entities are squeezed out.
* \param	gObj			Given constrained mesh transform object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformInvert2D(WlzObject *gObj,
				            WlzErrorNum *dstErr)
{
  int		*nTbl = NULL;
  WlzObject	*rObj = NULL;
  WlzIndexedValues *gIxv = NULL,
  		*rIxv = NULL;
  WlzCMesh2D	*gMesh = NULL,
  		*rMesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((gMesh = gObj->domain.cm2)->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    gIxv = gObj->values.x;
    if((gIxv->attach != WLZ_VALUE_ATTACH_NOD) ||
       (gIxv->rank != 1) ||
       (gIxv->dim[0] < 2) ||
       (gIxv->vType != WLZ_GREY_DOUBLE))
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(((errNum == WLZ_ERR_NONE) &&
     gMesh->res.nod.numEnt > 0) && (gMesh->res.elm.numEnt > 0))
  {
    /* Allocate new data structures. */
    rMesh = WlzCMeshNew2D(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(((nTbl = (int *)
                  AlcMalloc(sizeof(int) * gMesh->res.nod.numEnt)) == NULL) ||
         (AlcVectorExtend(rMesh->res.nod.vec,
                          gMesh->res.nod.numEnt) != ALC_ER_NONE) ||
	 (AlcVectorExtend(rMesh->res.elm.vec,
	                  gMesh->res.elm.numEnt) != ALC_ER_NONE))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        errNum = WlzCMeshReassignGridCells2D(rMesh, gMesh->res.nod.numEnt);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain dom;
      WlzValues val;

      dom.cm2 = rMesh;
      val.core = NULL;
      rObj = WlzMakeMain(gObj->type, dom, val, NULL, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rIxv = WlzMakeIndexedValues(rObj, gIxv->rank, gIxv->dim, gIxv->vType,
				  gIxv->attach, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues val;

      val.x = rIxv;
      rObj->values = WlzAssignValues(val, NULL);
    }
    /* Set node positions and displacements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idN;
      double	*gDsp,
      		*rDsp;
      AlcVector	*gVec;
      WlzCMeshNod2D *gNod,
      		    *rNod;

      gVec = gMesh->res.nod.vec;
      for(idN = 0; idN < gMesh->res.nod.maxEnt; ++idN)
      {
        gNod = (WlzCMeshNod2D *)AlcVectorItemGet(gVec, idN);
	if(gNod->idx >= 0)
	{
	  WlzDVertex2 pos;

	  gDsp = (double *)WlzIndexedValueGet(gIxv, idN);
	  pos.vtX = gNod->pos.vtX + gDsp[0];
	  pos.vtY = gNod->pos.vtY + gDsp[1];
	  rNod = WlzCMeshNewNod2D(rMesh, pos, NULL);
	  rDsp = (double *)WlzIndexedValueGet(rIxv, rNod->idx);
	  nTbl[gNod->idx] = rNod->idx;
	  rDsp[0] = -(gDsp[0]);
	  rDsp[1] = -(gDsp[1]);
	}
      }
      WlzCMeshUpdateBBox2D(rMesh);
      errNum = WlzCMeshReassignGridCells2D(rMesh, 0);
    }
    /* Create the elements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idE;
      AlcVector	*gVec,
      		*rVec;
      WlzCMeshElm2D *gElm;

      gVec = gMesh->res.elm.vec;
      rVec = rMesh->res.nod.vec;
      for(idE = 0; idE < gMesh->res.elm.maxEnt; ++idE)
      {
        gElm = (WlzCMeshElm2D *)AlcVectorItemGet(gVec, idE);
	if(gElm->idx >= 0)
	{
	  WlzCMeshNod2D *gNod;
	  WlzCMeshNod2D *rNodes[3];

	  gNod = WLZ_CMESH_ELM2D_GET_NODE_0(gElm);
	  rNodes[0] = (WlzCMeshNod2D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM2D_GET_NODE_1(gElm);
	  rNodes[1] = (WlzCMeshNod2D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM2D_GET_NODE_2(gElm);
	  rNodes[2] = (WlzCMeshNod2D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  (void )WlzCMeshNewElm2D(rMesh, rNodes[0], rNodes[1], rNodes[2],
				  1, &errNum);
          if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  AlcFree(nTbl);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      (void )WlzCMeshFree2D(rMesh);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Inverted constrained mesh transform.
* \ingroup	WlzTransform
* \brief	Inverts the given 2D5 constrained mesh transform.
* \param	mObj			Given constrained mesh transform object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformInvert2D5(WlzObject *gObj,
				             WlzErrorNum *dstErr)
{
  int		*nTbl = NULL;
  WlzObject	*rObj = NULL;
  WlzIndexedValues *gIxv = NULL,
  		*rIxv = NULL;
  WlzCMesh2D5	*gMesh = NULL,
  		*rMesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((gMesh = gObj->domain.cm2d5)->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    gIxv = gObj->values.x;
    if((gIxv->attach != WLZ_VALUE_ATTACH_NOD) ||
       (gIxv->rank != 1) ||
       (gIxv->dim[0] < 3) ||
       (gIxv->vType != WLZ_GREY_DOUBLE))
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(((errNum == WLZ_ERR_NONE) &&
     gMesh->res.nod.numEnt > 0) && (gMesh->res.elm.numEnt > 0))
  {
    /* Allocate new data structures. */
    rMesh = WlzCMeshNew2D5(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(((nTbl = (int *)
                  AlcMalloc(sizeof(int) * gMesh->res.nod.numEnt)) == NULL) ||
         (AlcVectorExtend(rMesh->res.nod.vec,
                          gMesh->res.nod.numEnt) != ALC_ER_NONE) ||
	 (AlcVectorExtend(rMesh->res.elm.vec,
	                  gMesh->res.elm.numEnt) != ALC_ER_NONE))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        errNum = WlzCMeshReassignGridCells2D5(rMesh, gMesh->res.nod.numEnt);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain dom;
      WlzValues val;

      dom.cm2d5 = rMesh;
      val.core = NULL;
      rObj = WlzMakeMain(gObj->type, dom, val, NULL, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rIxv = WlzMakeIndexedValues(rObj, gIxv->rank, gIxv->dim, gIxv->vType,
				  gIxv->attach, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues val;

      val.x = rIxv;
      rObj->values = WlzAssignValues(val, NULL);
    }
    /* Set node positions and displacements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idN;
      double	*gDsp,
      		*rDsp;
      AlcVector	*gVec;
      WlzCMeshNod2D5 *gNod,
      		    *rNod;

      gVec = gMesh->res.nod.vec;
      for(idN = 0; idN < gMesh->res.nod.maxEnt; ++idN)
      {
        gNod = (WlzCMeshNod2D5 *)AlcVectorItemGet(gVec, idN);
	if(gNod->idx >= 0)
	{
	  WlzDVertex3 pos;

	  gDsp = (double *)WlzIndexedValueGet(gIxv, idN);
	  pos.vtX = gNod->pos.vtX + gDsp[0];
	  pos.vtY = gNod->pos.vtY + gDsp[1];
	  pos.vtZ = gNod->pos.vtZ + gDsp[2];
	  rNod = WlzCMeshNewNod2D5(rMesh, pos, NULL);
	  rDsp = (double *)WlzIndexedValueGet(rIxv, rNod->idx);
	  nTbl[gNod->idx] = rNod->idx;
	  rDsp[0] = -(gDsp[0]);
	  rDsp[1] = -(gDsp[1]);
	  rDsp[2] = -(gDsp[2]);
	}
      }
      WlzCMeshUpdateBBox2D5(rMesh);
      errNum = WlzCMeshReassignGridCells2D5(rMesh, 0);
    }
    /* Create the elements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idE;
      AlcVector	*gVec,
      		*rVec;
      WlzCMeshElm2D5 *gElm;

      gVec = gMesh->res.elm.vec;
      rVec = rMesh->res.nod.vec;
      for(idE = 0; idE < gMesh->res.elm.maxEnt; ++idE)
      {
        gElm = (WlzCMeshElm2D5 *)AlcVectorItemGet(gVec, idE);
	if(gElm->idx >= 0)
	{
	  WlzCMeshNod2D5 *gNod;
	  WlzCMeshNod2D5 *rNodes[3];

	  gNod = WLZ_CMESH_ELM2D5_GET_NODE_0(gElm);
	  rNodes[0] = (WlzCMeshNod2D5 *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM2D5_GET_NODE_1(gElm);
	  rNodes[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM2D5_GET_NODE_2(gElm);
	  rNodes[2] = (WlzCMeshNod2D5 *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  (void )WlzCMeshNewElm2D5(rMesh, rNodes[0], rNodes[1], rNodes[2],
				   1, &errNum);
          if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  AlcFree(nTbl);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      (void )WlzCMeshFree2D5(rMesh);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Inverted constrained mesh transform.
* \ingroup	WlzTransform
* \brief	Inverts the given 3D constrained mesh transform.
* \param	mObj			Given mesh transform object.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformInvert3D(WlzObject *gObj,
				             WlzErrorNum *dstErr)
{
  int		*nTbl = NULL;
  WlzObject	*rObj = NULL;
  WlzIndexedValues *gIxv = NULL,
  		*rIxv = NULL;
  WlzCMesh3D	*gMesh = NULL,
  		*rMesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((gMesh = gObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    gIxv = gObj->values.x;
    if((gIxv->attach != WLZ_VALUE_ATTACH_NOD) ||
       (gIxv->rank != 1) ||
       (gIxv->dim[0] < 3) ||
       (gIxv->vType != WLZ_GREY_DOUBLE))
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(((errNum == WLZ_ERR_NONE) &&
     gMesh->res.nod.numEnt > 0) && (gMesh->res.elm.numEnt > 0))
  {
    /* Allocate new data structures. */
    rMesh = WlzCMeshNew3D(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(((nTbl = (int *)
                  AlcMalloc(sizeof(int) * gMesh->res.nod.numEnt)) == NULL) ||
         (AlcVectorExtend(rMesh->res.nod.vec,
                          gMesh->res.nod.numEnt) != ALC_ER_NONE) ||
	 (AlcVectorExtend(rMesh->res.elm.vec,
	                  gMesh->res.elm.numEnt) != ALC_ER_NONE))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        errNum = WlzCMeshReassignGridCells3D(rMesh, gMesh->res.nod.numEnt);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain dom;
      WlzValues val;

      dom.cm3 = rMesh;
      val.core = NULL;
      rObj = WlzMakeMain(gObj->type, dom, val, NULL, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rIxv = WlzMakeIndexedValues(rObj, gIxv->rank, gIxv->dim, gIxv->vType,
				  gIxv->attach, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues val;

      val.x = rIxv;
      rObj->values = WlzAssignValues(val, NULL);
    }
    /* Set node positions and displacements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idN;
      double	*gDsp,
      		*rDsp;
      AlcVector	*gVec;
      WlzCMeshNod3D *gNod,
      		    *rNod;

      gVec = gMesh->res.nod.vec;
      for(idN = 0; idN < gMesh->res.nod.maxEnt; ++idN)
      {
        gNod = (WlzCMeshNod3D *)AlcVectorItemGet(gVec, idN);
	if(gNod->idx >= 0)
	{
	  WlzDVertex3 pos;

	  gDsp = (double *)WlzIndexedValueGet(gIxv, idN);
	  pos.vtX = gNod->pos.vtX + gDsp[0];
	  pos.vtY = gNod->pos.vtY + gDsp[1];
	  pos.vtZ = gNod->pos.vtZ + gDsp[2];
	  rNod = WlzCMeshNewNod3D(rMesh, pos, NULL);
	  rDsp = (double *)WlzIndexedValueGet(rIxv, rNod->idx);
	  nTbl[gNod->idx] = rNod->idx;
	  rDsp[0] = -(gDsp[0]);
	  rDsp[1] = -(gDsp[1]);
	  rDsp[2] = -(gDsp[2]);
	}
      }
      WlzCMeshUpdateBBox3D(rMesh);
      errNum = WlzCMeshReassignGridCells3D(rMesh, 0);
    }
    /* Create the elements. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	idE;
      AlcVector	*gVec,
      		*rVec;
      WlzCMeshElm3D *gElm;

      gVec = gMesh->res.elm.vec;
      rVec = rMesh->res.nod.vec;
      for(idE = 0; idE < gMesh->res.elm.maxEnt; ++idE)
      {
        gElm = (WlzCMeshElm3D *)AlcVectorItemGet(gVec, idE);
	if(gElm->idx >= 0)
	{
	  WlzCMeshNod3D *gNod;
	  WlzCMeshNod3D *rNodes[4];

	  gNod = WLZ_CMESH_ELM3D_GET_NODE_0(gElm);
	  rNodes[0] = (WlzCMeshNod3D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM3D_GET_NODE_1(gElm);
	  rNodes[1] = (WlzCMeshNod3D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM3D_GET_NODE_2(gElm);
	  rNodes[2] = (WlzCMeshNod3D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  gNod = WLZ_CMESH_ELM3D_GET_NODE_3(gElm);
	  rNodes[3] = (WlzCMeshNod3D *)AlcVectorItemGet(rVec, nTbl[gNod->idx]);
	  (void )WlzCMeshNewElm3D(rMesh, rNodes[0], rNodes[1], rNodes[2],
				  rNodes[3], 1, &errNum);
          if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  AlcFree(nTbl);
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else
    {
      (void )WlzCMeshFree3D(rMesh);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return       New conforming mesh transform.
* \ingroup      WlzTransform
* \brief        Creates a conforming mesh transform for the given object
*		with all mesh displacements allocated but set to zero.
*		A mesh transform is a conforming mesh object with indexed
*		values such that the values are double precision displacements
*		that are ordered x, y[, z].
* \param        srcObj                  The given object.
* \param        method                  Mesh generation method to use.
* \param        minDist                 Minimum distance between mesh vertices.
* \param        maxDist                 Maximum distance between mesh vertices.
* \param	dstDilObj		Destination pointer for the dilated
*					object used to build the mesh, may
*					be NULL.
* \param	delOut			Delete all elements with nodes
*					outside the object if non-zero.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshTransformFromObj(WlzObject *srcObj,
				 WlzMeshGenMethod method,
                                 double minDist, double maxDist,
				 WlzObject **dstDilObj, int delOut,
                                 WlzErrorNum *dstErr)
{
  int		dim = 0;
  double	dsp[3];
  WlzDomain	domain;
  WlzValues	values;
  WlzObject 	*mObj = NULL;
  WlzObjectType	mObjType = WLZ_NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  domain.core = NULL;
  values.core = NULL;
  dsp[0] = dsp[1] = dsp[2] = 0.0;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(srcObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(method != WLZ_MESH_GENMETHOD_CONFORM)
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	dim = 2;
	domain.cm2 = WlzCMeshFromObj2D(srcObj, minDist, maxDist, dstDilObj,
				       delOut, &errNum);
	mObjType = WLZ_CMESH_2D;
        break;
      case WLZ_3D_DOMAINOBJ:
	dim = 3;
	domain.cm3 = WlzCMeshFromObj3D(srcObj, minDist, maxDist, dstDilObj,
				       delOut, &errNum);
	mObjType = WLZ_CMESH_3D;
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzMakeMain(mObjType, domain, values, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    values.x = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
                                    WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj->values = WlzAssignValues(values, NULL);
    errNum = WlzIndexedValuesSet(mObj, dim * sizeof(double), (void *)dsp);
  }
  else
  {
    if(mObj)
    {
      (void )WlzFreeObj(mObj);
    }
    else
    {
      WlzFreeDomain(domain);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given integer vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2I(WlzObject *mObj,
					 int nVtx, WlzIVertex2 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  double	*dsp;
  WlzDVertex2	tVtx;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm2;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mesh, lastElmIdx,
			(double )(vtx[idN].vtX), (double )(vtx[idN].vtY),
		        0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
    }
    vtx[idN].vtX = WLZ_NINT(tVtx.vtX);
    vtx[idN].vtY = WLZ_NINT(tVtx.vtY);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given integer vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry3I(WlzObject *mObj,
					 int nVtx, WlzIVertex3 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  double	*dsp;
  WlzDVertex3	tVtx;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm3D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm3;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos3D(mesh, lastElmIdx,
			(double )(vtx[idN].vtX), (double )(vtx[idN].vtY),
		        (double )(vtx[idN].vtZ), 0,
			&nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm3D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.tr[ 0] * vtx[idN].vtX) + (sE.tr[ 1] * vtx[idN].vtY) +
		 (sE.tr[ 2] * vtx[idN].vtZ) +  sE.tr[ 3];
      tVtx.vtY = (sE.tr[ 4] * vtx[idN].vtX) + (sE.tr[ 5] * vtx[idN].vtY) +
		 (sE.tr[ 6] * vtx[idN].vtZ) +  sE.tr[ 7];
      tVtx.vtZ = (sE.tr[ 8] * vtx[idN].vtX) + (sE.tr[ 9] * vtx[idN].vtY) +
		 (sE.tr[10] * vtx[idN].vtZ) +  sE.tr[11];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
      tVtx.vtZ = vtx[idN].vtY + dsp[2];
    }
    vtx[idN].vtX = WLZ_NINT(tVtx.vtX);
    vtx[idN].vtY = WLZ_NINT(tVtx.vtY);
    vtx[idN].vtZ = WLZ_NINT(tVtx.vtZ);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given float vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* 					vertices in the array after
* 					transformation. Must not be NULL.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2F(WlzObject *mObj,
					 int nVtx, WlzFVertex2 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  double	*dsp;
  WlzDVertex2	tVtx;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm2;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mesh, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY,
			0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
    }
    vtx[idN].vtX = (float )(tVtx.vtX);
    vtx[idN].vtY = (float )(tVtx.vtY);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given float vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* 					vertices in the array after
* 					transformation. Must not be NULL.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry3F(WlzObject *mObj,
					 int nVtx, WlzFVertex3 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  double	*dsp;
  WlzDVertex3	tVtx;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm3D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm3;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos3D(mesh, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY, vtx[idN].vtZ,
			0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm3D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.tr[ 0] * vtx[idN].vtX) + (sE.tr[ 1] * vtx[idN].vtY) +
		 (sE.tr[ 2] * vtx[idN].vtZ) +  sE.tr[ 3];
      tVtx.vtY = (sE.tr[ 4] * vtx[idN].vtX) + (sE.tr[ 5] * vtx[idN].vtY) +
		 (sE.tr[ 6] * vtx[idN].vtZ) +  sE.tr[ 7];
      tVtx.vtZ = (sE.tr[ 8] * vtx[idN].vtX) + (sE.tr[ 9] * vtx[idN].vtY) +
		 (sE.tr[10] * vtx[idN].vtZ) +  sE.tr[11];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
      tVtx.vtZ = vtx[idN].vtZ + dsp[2];
    }
    vtx[idN].vtX = (float )(tVtx.vtX);
    vtx[idN].vtY = (float )(tVtx.vtY);
    vtx[idN].vtZ = (float )(tVtx.vtZ);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2D(WlzObject *mObj,
					 int nVtx, WlzDVertex2 *vtx)
{
  int		idN,
		nearNod,
  		lastElmIdx;
  double	*dsp;
  WlzDVertex2	tVtx;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm2;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mesh, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY,
			0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
    }
    vtx[idN] = tVtx;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2D5(WlzObject *mObj,
					   int nVtx, WlzDVertex3 *vtx)
{
  int		idN,
		nearNod,
  		lastElmIdx;
  double	*dsp;
  WlzDVertex3	tVtx;
  WlzCMesh2D5	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm3D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm2d5;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D5(mesh, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY, vtx[idN].vtZ,
			0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D5(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.tr[ 0] * vtx[idN].vtX) + (sE.tr[ 1] * vtx[idN].vtY) +
		 (sE.tr[ 2] * vtx[idN].vtZ) +  sE.tr[ 3];
      tVtx.vtY = (sE.tr[ 4] * vtx[idN].vtX) + (sE.tr[ 5] * vtx[idN].vtY) +
		 (sE.tr[ 6] * vtx[idN].vtZ) +  sE.tr[ 7];
      tVtx.vtZ = (sE.tr[ 8] * vtx[idN].vtX) + (sE.tr[ 9] * vtx[idN].vtY) +
		 (sE.tr[10] * vtx[idN].vtZ) +  sE.tr[11];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
      tVtx.vtZ = vtx[idN].vtZ + dsp[2];
    }
    vtx[idN] = tVtx;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mObj			The mesh transform object.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry3D(WlzObject *mObj,
					 int nVtx, WlzDVertex3 *vtx)
{
  int		idN,
		nearNod,
  		lastElmIdx;
  double	*dsp;
  WlzDVertex3	tVtx;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshScanElm3D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  mesh = mObj->domain.cm3;
  ixv = mObj->values.x;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos3D(mesh, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY, vtx[idN].vtZ,
			0, &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm3D(mObj, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.tr[ 0] * vtx[idN].vtX) + (sE.tr[ 1] * vtx[idN].vtY) +
		 (sE.tr[ 2] * vtx[idN].vtZ) +  sE.tr[ 3];
      tVtx.vtY = (sE.tr[ 4] * vtx[idN].vtX) + (sE.tr[ 5] * vtx[idN].vtY) +
		 (sE.tr[ 6] * vtx[idN].vtZ) +  sE.tr[ 7];
      tVtx.vtZ = (sE.tr[ 8] * vtx[idN].vtX) + (sE.tr[ 9] * vtx[idN].vtY) +
		 (sE.tr[10] * vtx[idN].vtZ) +  sE.tr[11];
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nearNod);
      tVtx.vtX = vtx[idN].vtX + dsp[0];
      tVtx.vtY = vtx[idN].vtY + dsp[1];
      tVtx.vtZ = vtx[idN].vtZ + dsp[2];
    }
    vtx[idN] = tVtx;
  }
  return(errNum);
}

/*!
* \return       New domain object or NULL on error.
* \ingroup      WlzMesh
* \brief        Computes a new domain object, the domain of which corresponds
*               to the region of space enclosed by the mesh, such that the
*               domain object is covered by the given mesh object.
* \param        mObj                    Given mesh object.
* \param 	trans			If non zero domain corresponds to the
* 					transformed (instead of the
* 					untransformed mesh). For a transform,
* 					the objects values must be appropriate.
* \param	scale			Additional scale factor from mesh
* 					to spatial domain.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToDomObj(WlzObject *mObj, int trans, double scale,
				  WlzErrorNum *dstErr)
{
  WlzObject     *dObj = NULL,
  		*tObj = NULL;
  WlzAffineTransform *tr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const double	eps = 0.000001;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(fabs(scale) < eps)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(mObj->domain.core->type)
    {
      case WLZ_CMESH_2D:
        dObj = WlzCMeshToDomObj2D(mObj, trans, &errNum);
	if((errNum == WLZ_ERR_NONE) && (fabs(scale - 1.0) > eps))
	{
	  tr = WlzAffineTransformFromScale(WLZ_TRANSFORM_2D_AFFINE,
	  			scale, scale, 1.0, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj = WlzAffineTransformObj(dObj, tr, WLZ_INTERPOLATION_NEAREST,
					 &errNum);
	  }
	  (void )WlzFreeAffineTransform(tr);
	  (void )WlzFreeObj(dObj);
	  dObj = tObj;
	}
        break;
      case WLZ_CMESH_3D:
        dObj = WlzCMeshToDomObj3D(mObj, trans, &errNum);
	if((errNum == WLZ_ERR_NONE) && (fabs(scale - 1.0) > eps))
	{
	  tr = WlzAffineTransformFromScale(WLZ_TRANSFORM_3D_AFFINE,
	  			scale, scale, scale, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj = WlzAffineTransformObj(dObj, tr, WLZ_INTERPOLATION_NEAREST,
					 &errNum);
	  }
	  (void )WlzFreeAffineTransform(tr);
	  (void )WlzFreeObj(dObj);
	  dObj = tObj;
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dObj);
}

/*!
* \return	New domain object (image) with values interpolated from mesh
* 		values.
* \ingroup      WlzMesh
* \brief	Given a domain object and a mesh object, this function
* 		creates a new domain object using the domain of the given
* 		object and a new value table. The function then interpolates
* 		the mesh values throughout the domain object.
* \param	dObj			Given domain object.
* \param	mObj			Given mesh object.
* \param	itp			Interpolation method.
* \param	ixi			Index into the indexed values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshToDomObjValues(WlzObject *dObj, WlzObject *mObj,
                                        WlzInterpolationType itp, int ixi,
					WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dObj == NULL) || (mObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((dObj->domain.core == NULL) || (mObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;;
  }
  else if(mObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((mObj->values.core->type != (WlzObjectType )WLZ_INDEXED_VALUES) ||
          (mObj->values.x->rank < 0))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    int		i,
    		nv = 1;

    for(i = 0; i < mObj->values.x->rank; ++i)
    {
      nv *= mObj->values.x->dim[i];
    }
    if(ixi >= nv)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if((dObj->type == WLZ_2D_DOMAINOBJ) && (mObj->type == WLZ_CMESH_2D))
  {
    rObj = WlzCMeshToDomObjValues2D(dObj, mObj, itp, ixi, &errNum);
  }
  else if((dObj->type == WLZ_3D_DOMAINOBJ) && (mObj->type == WLZ_CMESH_3D))
  {
    rObj = WlzCMeshToDomObjValues3D(dObj, mObj, itp, ixi, &errNum);
  }
  else
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return       New 2D domain object, empty object if the mesh has no elements
*               or NULL on error.
* \ingroup      WlzMesh
* \brief        Computes a new 2D domain object, the domain of which
*               corresponds to the region of space enclosed by the 2D mesh,
*               such that the domain object is covered by the given mesh.
* \param        mTr                    	Given 2D mesh object.
* \param	trans			If non zero domain corresponds to the
* 					transformed (instead of the
* 					untransformed mesh). For a transform,
* 					the objects values must be appropriate.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshToDomObj2D(WlzObject *mObj, int trans,
				     WlzErrorNum *dstErr)
{
  int           idI = 0,
                line = 0,
                itvLnCnt = 0,
                itvLnWidth = 0,
                itvLnByteWidth = 0;
  WlzIBox2      iBox;
  WlzDBox2      dBox;
  WlzCMesh2D	*mesh;
  WlzCMeshScanWSp2D *mSWSp = NULL;
  WlzCMeshScanItv2D *curItv = NULL,
  		*prvItv = NULL;
  WlzUByte	*lnMsk = NULL;
  WlzObject     *dObj = NULL;
  WlzDomain     dom;
  WlzValues     val;
  WlzDynItvPool itvPool;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  iBox.xMin = iBox.yMin = iBox.xMax = iBox.yMax = 0;
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((mesh = mObj->domain.cm2) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((trans != 0) && (mObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(mesh->res.elm.numEnt == 0)
  {
    dObj = WlzMakeEmpty(&errNum);
  }
  else
  {
    mSWSp = WlzCMeshScanWSpInit2D(mObj, trans, &errNum);
    /* Create a single line bit mask buffer. */
    if(errNum == WLZ_ERR_NONE)
    {
      dBox = WlzCMeshTransformGetBBox2D(mObj, trans, &errNum);
      iBox.xMin = WLZ_CMESH_POS_DTOI(dBox.xMin);
      iBox.xMax = WLZ_CMESH_POS_DTOI(dBox.xMax);
      iBox.yMin = WLZ_CMESH_POS_DTOI(dBox.yMin);
      iBox.yMax = WLZ_CMESH_POS_DTOI(dBox.yMax);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      itvLnWidth = iBox.xMax - iBox.xMin + 1;
      itvLnByteWidth = (itvLnWidth + 7) / 8;
      if(AlcBit1Calloc(&lnMsk, itvLnWidth) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    /* Create the interval domain. */
    if(errNum == WLZ_ERR_NONE)
    {
      dom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				    iBox.yMin, iBox.yMax,
				    iBox.xMin, iBox.xMax,
				    &errNum);
    }
    /* Set up mesh workspace intervals. */
    if(errNum == WLZ_ERR_NONE)
    {
      idI = 0;
      itvLnCnt = 0;
      itvPool.offset = 0;
      itvPool.itvBlock = NULL;
      itvPool.itvsInBlock = 1024; /* This is just the number of intervals to
      				   *  allocate at once - an efficiency
				   *  parameter. */
      prvItv = NULL;
      curItv = mSWSp->itvs;
      while((idI < mSWSp->nItvs) && (curItv->line < iBox.yMin))
      {
	++curItv;
	++idI;
      }
    }
    /* Scan through mesh workspace intervals adding them to the interval
     * domain. */
    while((errNum == WLZ_ERR_NONE) && (idI < mSWSp->nItvs))
    {
      line = curItv->line;
      while((idI < mSWSp->nItvs) && (curItv->line == line))
      {
	WlzBitLnSetItv(lnMsk,
		       curItv->lftI - iBox.xMin, curItv->rgtI - iBox.xMin,
		       itvLnWidth);
	prvItv = curItv;
	++idI;
	++curItv;
	++itvLnCnt;
      }
      /* Add line to interval domain. */
      if(itvLnCnt > 0)
      {
	errNum = WlzDynItvLnFromBitLn(dom.i, lnMsk, prvItv->line,
				      itvLnWidth, &itvPool);
	memset(lnMsk, 0, itvLnByteWidth);
	itvLnCnt = 0;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzStandardIntervalDomain(dom.i);
      dObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
    }
  }
  AlcFree(lnMsk);
  WlzCMeshScanWSpFree2D(mSWSp);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeDomain(dom);
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dObj);
}

/*!
* \return       New 3D domain object, empty object if the mesh has no elements
*               or NULL on error.
* \ingroup      WlzMesh
* \brief        Computes a new 3D domain object, the domain of which
*               corresponds to the region of space enclosed by the 3D mesh,
*               such that the domain object is covered by the given mesh.
* \param        mObj                   	Given 2D mesh object.
* \param	trans			If non zero domain corresponds to the
* 					transformed (instead of the
* 					untransformed mesh). For a transform,
* 					the objects values must be appropriate.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshToDomObj3D(WlzObject *mObj, int trans,
				     WlzErrorNum *dstErr)
{
  WlzObject     *dobj = NULL;
  WlzCMesh3D	*mesh;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((mesh = mObj->domain.cm3) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((trans != 0) && (mObj->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(mesh->res.elm.numEnt == 0)
  {
    dobj = WlzMakeEmpty(&errNum);
  }
  else
  {
    /* Make workspace intervals for the elements in the displaced
     * mesh, with intervals sorted by plane, line and then column. */
    if(errNum == WLZ_ERR_NONE)
    {
      mSWSp = WlzCMeshScanWSpInit3D(mObj, trans, &errNum);
    }
    /* Scan through the sorted intervals creating domains as required. */
    if(errNum == WLZ_ERR_NONE)
    {
      dobj = WlzCMeshScanObjPDomain3D(NULL, mSWSp, &errNum); 
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dobj);
}

/*!
* \return	New domain object (image) with values interpolated from mesh
* 		values.
* \ingroup      WlzMesh
* \brief	Given a domain object and a mesh object, this function
* 		creates a new domain object using the domain of the given
* 		object and a new value table. The function then interpolates
* 		the mesh values throughout the domain object.
* \param	dObj			Given domain object, must be a valid
* 					WLZ_2D_DOMAINOBJ.
* \param	mObj			Given mesh object, must be a valid
* 					WLZ_CMESH_2D object with valid indexed
* 					values..
* \param	itp			Interpolation method.
* \param	ixi			Index into the indexed values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshToDomObjValues2D(WlzObject *dObj, WlzObject *mObj,
                                        WlzInterpolationType itp, int ixi,
					WlzErrorNum *dstErr)
{
  WlzValues	rVal;
  WlzObjectType	rVTT;
  WlzPixelV	bgd;
  WlzIndexedValues *ixv;
  WlzCMesh2D	*mesh;
  WlzGreyWSpace gWsp;
  WlzIntervalWSpace iWsp;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rVal.core = NULL;
  ixv = mObj->values.x;
  mesh = mObj->domain.cm2;
  bgd.type = ixv->vType;
  switch(bgd.type)
  {
    case WLZ_GREY_LONG:
      bgd.v.lnv = 0;
      break;
    case WLZ_GREY_INT:
      bgd.v.inv = 0;
      break;
    case WLZ_GREY_SHORT:
      bgd.v.shv = 0;
      break;
    case WLZ_GREY_UBYTE:
      bgd.v.ubv = 0;
      break;
    case WLZ_GREY_FLOAT:
      bgd.v.flv = 0.0f;
      break;
    case WLZ_GREY_DOUBLE:
      bgd.v.dbv = 0.0;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rVTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, bgd.type, NULL);
    rVal.v = WlzNewValueTb(dObj, rVTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ,
                       dObj->domain, rVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzInitGreyScan(rObj, &iWsp, &gWsp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    while((errNum == WLZ_ERR_NONE) &&
          ((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE))
    {
      switch(ixv->attach)
      {
	case WLZ_VALUE_ATTACH_NOD:
	  switch(itp)
	  {
	    case WLZ_INTERPOLATION_NEAREST:
	      errNum = WlzCMeshInterpolateNod2DNearest(gWsp.u_grintptr,
	      				iWsp.linpos, iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
	      break;
	    case WLZ_INTERPOLATION_LINEAR: /* FALLTHROUGH */
	    case WLZ_INTERPOLATION_BARYCENTRIC:
	      errNum = WlzCMeshInterpolateNod2DLinear(gWsp.u_grintptr, 
				     iWsp.linpos, iWsp.lftpos, iWsp.rgtpos,
				     mesh, ixv, ixi);
	      break;
	    case WLZ_INTERPOLATION_KRIG:
	      errNum = WlzCMeshInterpolateNod2DKrig(gWsp.u_grintptr, 
				    iWsp.linpos, iWsp.lftpos, iWsp.rgtpos,
				    mesh, ixv, ixi);
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	  break;
	case WLZ_VALUE_ATTACH_ELM:
	  switch(itp)
	  {
	    case WLZ_INTERPOLATION_NEAREST:
	      errNum = WlzCMeshInterpolateElm2DNearest(gWsp.u_grintptr,
	      				iWsp.linpos, iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
	      break;
	    case WLZ_INTERPOLATION_LINEAR:
	      errNum = WlzCMeshInterpolateElm2DLinear(gWsp.u_grintptr,
	      				iWsp.linpos, iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	  break;
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	  break;
	default:
	  errNum = WLZ_ERR_VALUES_TYPE;
	  break;
      }
    }
    (void )WlzEndGreyScan(&iWsp, &gWsp);
    if(errNum == WLZ_ERR_EOO)
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj != NULL)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else if(rVal.core != NULL)
    {
      (void )WlzFreeValueTb(rVal.v);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New domain object (image) with values interpolated from mesh
* 		values.
* \ingroup      WlzMesh
* \brief	Given a domain object and a mesh object, this function
* 		creates a new domain object using the domain of the given
* 		object and a new value table. The function then interpolates
* 		the mesh values throughout the domain object.
* \param	dObj			Given domain object, must be a valid
* 					WLZ_3D_DOMAINOBJ.
* \param	mObj			Given mesh object, must be a valid
* 					WLZ_CMESH_3D object with valid indexed
* 					values..
* \param	itp			Interpolation method.
* \param	ixi			Index into the indexed values.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshToDomObjValues3D(WlzObject *dObj, WlzObject *mObj,
                                        WlzInterpolationType itp, int ixi,
					WlzErrorNum *dstErr)
{
  WlzValues	rVal;
  WlzObjectType	rVTT;
  WlzPixelV	bgd;
  WlzIndexedValues *ixv;
  WlzCMesh3D	*mesh;
  WlzGreyWSpace gWsp;
  WlzIntervalWSpace iWsp;
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  rVal.core = NULL;
  ixv = mObj->values.x;
  mesh = mObj->domain.cm3;
  bgd.type = ixv->vType;
  switch(bgd.type)
  {
    case WLZ_GREY_LONG:
      bgd.v.lnv = 0;
      break;
    case WLZ_GREY_INT:
      bgd.v.inv = 0;
      break;
    case WLZ_GREY_SHORT:
      bgd.v.shv = 0;
      break;
    case WLZ_GREY_UBYTE:
      bgd.v.ubv = 0;
      break;
    case WLZ_GREY_FLOAT:
      bgd.v.flv = 0.0f;
      break;
    case WLZ_GREY_DOUBLE:
      bgd.v.dbv = 0.0;
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rVTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, bgd.type, NULL);
    rVal.vox = WlzNewValuesVox(dObj, rVTT, bgd, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ,
                       dObj->domain, rVal, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idP,
    		pCnt;
    
    pCnt = rObj->domain.p->lastpl - rObj->domain.p->plane1 + 1;
    for(idP = 0; idP < pCnt; ++idP)
    {
      int	plnPos;
      WlzObject *objT = NULL;

      objT = WlzMakeMain(WLZ_2D_DOMAINOBJ,
                         *(rObj->domain.p->domains + idP),
			 *(rObj->values.vox->values + idP),
			 NULL, NULL, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        plnPos = rObj->domain.p->plane1 + idP;
	errNum = WlzInitGreyScan(objT, &iWsp, &gWsp);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while((errNum == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE))
	{
	  switch(ixv->attach)
	  {
	    case WLZ_VALUE_ATTACH_NOD:
	      switch(itp)
	      {
	        case WLZ_INTERPOLATION_NEAREST:
		  errNum = WlzCMeshInterpolateNod3DNearest(gWsp.u_grintptr,
		  			plnPos, iWsp.linpos,
					iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
		  break;
	        case WLZ_INTERPOLATION_LINEAR: /* FALLTHROUGH */
	        case WLZ_INTERPOLATION_BARYCENTRIC:
		  errNum = WlzCMeshInterpolateNod3DLinear(gWsp.u_grintptr, 
					plnPos, iWsp.linpos,
					iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
		  break;
		default:
		  errNum = WLZ_ERR_PARAM_TYPE;
		  break;
	      }
	      break;
	    case WLZ_VALUE_ATTACH_ELM:
	      switch(itp)
	      {
	        case WLZ_INTERPOLATION_NEAREST:
		  errNum = WlzCMeshInterpolateElm3DNearest(gWsp.u_grintptr, 
					plnPos, iWsp.linpos,
					iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
		  break;
	        case WLZ_INTERPOLATION_LINEAR:
                  errNum = WlzCMeshInterpolateElm3DLinear(gWsp.u_grintptr,
		  			plnPos, iWsp.linpos,
					iWsp.lftpos, iWsp.rgtpos,
					mesh, ixv, ixi);
		  break;
		default:
		  errNum = WLZ_ERR_PARAM_TYPE;
		  break;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_VALUES_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      (void )WlzFreeObj(objT);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj != NULL)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else if(rVal.core != NULL)
    {
      (void )WlzFreeValueTb(rVal.v);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Interpolates values along the pixels of a single interval
* 		from the given mesh and it's indexed values. Uses kriging
* 		for interpolation within each mesh element.
* 		The mesh must have a valid maxSqEdgLen before calling
* 		this function.
* \param	dst				The interval values.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi			Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateNod2DKrig(WlzGreyP dst,
                                          int ln, int kolL, int kolR,
					  WlzCMesh2D *mesh,
					  WlzIndexedValues *ixv, int ixi)
{
  int		kl,
		idI,
  		idE0,
		idE1,
		idN0,
		idN1,
		nNbr0 = 0,
		nNbr1 = 0,
		nbrChange = 0,
		maxKrigBuf = 0,
		maxNbrIdxBuf = 0;
  double	dRange;
  int		*wSp = NULL,
  		*nbrIdxBuf = NULL;
  double	*posSV = NULL;
  AlgMatrix	modelSV;
  WlzDVertex2	pos;
  WlzDVertex2	*nbrPosBuf = NULL;
  WlzKrigModelFn modelFn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idE0 = -1;
  idN0 = -1;
  modelSV.core = NULL;
  dRange = sqrt(mesh->maxSqEdgLen);
  pos.vtY = ln;
  for(idI = 0, kl = kolL; kl <= kolR; ++idI, ++kl)
  {
    pos.vtX = kl;
    /* Find element enclosing the voxel or if that doesn't exist the
     * nearest node to the voxel. Then find the fing of nodes that
     * include and surround either this element or node. */
    idE1 = WlzCMeshElmEnclosingPos2D(mesh, idE0, pos.vtX, pos.vtY, 0, &idN1);
    if(idE1 >= 0)
    {
      if(idE1 != idE0)
      {
	WlzCMeshElm2D	*elm;

	elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE1);
	nNbr1 = WlzCMeshElmRingNodIndices2D(elm, &maxNbrIdxBuf, &nbrIdxBuf,
					    &errNum);
	nbrChange = 1;
        idE0 = idE1;
	idN0 = -1;
      }
    }
    else if(idN1 >= 0)
    {
      if(idN1 != idN0)
      {
        WlzCMeshNod2D	*nod;

	nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN1);
	nNbr1 = WlzCMeshNodRingNodIndices2D(nod, &maxNbrIdxBuf, &nbrIdxBuf,
					    &errNum);
	nbrChange = 1;
        idN0 = idN1;
	idE0 = -1;
      }
    }
    /* The neighbouring node ring has changed so the position buffers and 
     * it's semi-variogram need updating. */
    if(nNbr1 > 0)
    {
      if(nbrChange)
      {
	/* Check for an increased number of neighbours and reallocate
	 * the buffers if needed. */
	if(nNbr1 > nNbr0)
	{
	  errNum = WlzKrigReallocBuffers2D(&nbrPosBuf, &posSV, &wSp, &modelSV,
	  			           &maxKrigBuf, nNbr1, nNbr0);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  int	i;

	  for(i = 0; i < nNbr1; ++i)
	  {
	    WlzCMeshNod2D *nod;
	    nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec,
	                                            nbrIdxBuf[i]);
	    nbrPosBuf[i] = nod->pos;
	  }
	  WlzKrigSetModelFn(&modelFn, WLZ_KRIG_MODELFN_LINEAR,
	  		    0.0, 0.1, 2.0 * dRange);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzKrigOSetModelSV2D(modelSV, &modelFn, nNbr1, nbrPosBuf,
	  				wSp);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzKrigOSetPosSV2D(posSV, &modelFn, nNbr1, nbrPosBuf, pos);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzKrigOWeightsSolve(modelSV, posSV, wSp, WLZ_MESH_TOLERANCE);
	/* posSV now contains the weights. */
	WlzIndexedValueBufWeight(dst, idI, ixv, posSV, nNbr1, nbrIdxBuf, ixi);
      }
      nNbr0 = nNbr1;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  AlgMatrixFree(modelSV);
  AlcFree(posSV);
  AlcFree(wSp);
  AlcFree(nbrIdxBuf);
  AlcFree(nbrPosBuf);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Interpolates values along the pixels of a single interval
* 		from the given mesh and it's indexed values. Uses nearest
* 		node interpolation within each mesh element.
* \param	dst				The interval values.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateNod2DNearest(WlzGreyP dst,
					int ln, int kolL, int kolR,
					WlzCMesh2D *mesh,
					WlzIndexedValues *ixv, int ixi)
{
  int		kl,
  		idE,
		idI,
		idN;
  WlzDVertex2	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    WlzCMeshNod2D *nNod = NULL;

    pos.vtX = kl;
    /* Find nearest node within the element. */
    if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE, pos.vtX, pos.vtY,
    					0, &idN)) >= 0)
    {
      int	idN;
      double    dMin;
      WlzDVertex2 del;
      WlzCMeshElm2D *elm;
      WlzCMeshNod2D *nod[3];

      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      nod[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm);
      nod[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm);
      nod[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm);
      nNod = nod[0];
      WLZ_VTX_2_SUB(del, pos, nod[0]->pos);
      dMin = WLZ_VTX_2_SQRLEN(del);
      for(idN = 1; idN < 3; ++idN)
      {
	double	d;

        WLZ_VTX_2_SUB(del, pos, nod[idN]->pos);
	d = WLZ_VTX_2_SQRLEN(del);
	if(d < dMin)
	{
	  nNod = nod[idN];
	  dMin = d;
	}
      }
    }
    else if(idN >= 0)
    {
      nNod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
    }
    /* Set value. */
    if(nNod)
    {
      WlzGreyP v;

      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  v.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.lnp + idI) = v.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  v.inp = (int *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.inp + idI) = v.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  v.shp = (short *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.shp + idI) = v.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  v.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.ubp + idI) = v.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  v.flp = (float *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.flp + idI) = v.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  v.dbp = (double *)WlzIndexedValueGet(ixv, nNod->idx);
	  *(dst.dbp + idI) = v.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Interpolates values along the pixels of a single interval
* 		from the given mesh and it's indexed values. Uses linear
* 		interpolation within each mesh element.
* \param	dst				The interval values.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateNod2DLinear(WlzGreyP dst,
					int ln, int kolL, int kolR,
					WlzCMesh2D *mesh,
					WlzIndexedValues *ixv, int ixi)
{
  int		kl,
  		idE,
		idI,
		idN;
  WlzDVertex2	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE, pos.vtX, pos.vtY, 0, 
                                        &idN)) >= 0)
    {
      int	idC;
      double    d[4];
      WlzCMeshElm2D *elm;
      WlzCMeshNod2D *nod[3];

      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      nod[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm);
      nod[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm);
      nod[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm);
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    WlzLong *x;

	    x = (WlzLong *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[3] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  dst.lnp[idI] = WLZ_NINT(d[3]);
	  break;
	case WLZ_GREY_INT:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    int	 *x;

	    x = (int *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[3] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  dst.inp[idI] = WLZ_NINT(d[3]);
	  break;
	case WLZ_GREY_SHORT:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    short *x;

	    x = (short *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[3] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  dst.shp[idI] = WLZ_NINT(d[3]);
	  break;
	case WLZ_GREY_UBYTE:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    WlzUByte *x;

	    x = (WlzUByte *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[3] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  dst.ubp[idI] = WLZ_NINT(d[3]);
	  break;
	case WLZ_GREY_FLOAT:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    float *x;

	    x = (float *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[3] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  dst.flp[idI] = d[3];
	  break;
	case WLZ_GREY_DOUBLE:
	  for(idC = 0; idC < 3; ++idC)
	  {
	    double *x;

	    x = (double *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  dst.dbp[idI] = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, d[0], d[1], d[2], pos);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    else if(idN >= 0)
    {
      WlzGreyP x;
      WlzCMeshNod2D *nod;

      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.lnp + idI) = x.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  x.inp = (int *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.inp + idI) = x.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  x.shp = (short *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.shp + idI) = x.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.ubp + idI) = x.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  x.flp = (float *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.flp + idI) = x.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  x.dbp = (double *)WlzIndexedValueGet(ixv, nod->idx);
	  *(dst.dbp + idI) = x.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Sets values along the pixels of a single interval
* 		from the given mesh and it's indexed values. The pixels
* 		are set to the first value attached to the mesh element.
* \param	dst				The interval values.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateElm2DNearest(WlzGreyP dst,
					int ln, int kolL, int kolR,
					WlzCMesh2D *mesh,
					WlzIndexedValues *ixv, int ixi)
{
  int		kl,
  		idE,
		idI;
  WlzDVertex2	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE, pos.vtX, pos.vtY,
    					0, NULL)) >= 0)
    {
      WlzGreyP	x;
      WlzCMeshElm2D *elm;

      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.lnp + idI) = x.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  x.inp = (int *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.inp + idI) = x.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  x.shp = (short *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.ubp + idI) = x.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.ubp + idI) = x.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  x.flp = (float *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.flp + idI) = x.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  x.dbp = (double *)WlzIndexedValueGet(ixv, elm->idx);
	  *(dst.dbp + idI) = x.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Sets values along the pixels of a single interval
* 		from the given mesh and it's indexed values. The pixels
* 		are set using linear interpolation.
* \param	dst				The interval values.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateElm2DLinear(WlzGreyP dst,
					int ln, int kolL, int kolR,
					WlzCMesh2D *mesh,
					WlzIndexedValues *ixv, int ixi)
{
  int		kl,
  		idE,
		idI,
		maxNbr = 0;
  int		*nbrIdxBuf = NULL;
  double	*nbrDstBuf = NULL;
  WlzDVertex2	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    int		nE = 0;

    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE, pos.vtX, pos.vtY,
    					0, NULL)) >= 0)
    {
      WlzCMeshElm2D *elm;
      int	    lastMaxNbr;

      lastMaxNbr = maxNbr;
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      /* Find elements connected to the given element. */
      nE = WlzCMeshElmRingElmIndices2D(elm, &maxNbr, &nbrIdxBuf, &errNum);
      /* Reallocate the distance buffer if required. */
      if((errNum == WLZ_ERR_NONE) && (lastMaxNbr < maxNbr))
      {
	nbrDstBuf = AlcRealloc(nbrDstBuf, maxNbr * sizeof(double));
        if(nbrDstBuf == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	int	setElmIdx = -1;
	double	v = 0.0;

        /* Compute inverse of distances to centroids of these elements from
	 * position. */
	for(idE = 0; idE < nE; ++idE)
	{
	  double	d;
	  WlzDVertex2	cen,
	  		del;
	  WlzCMeshElm2D *elm0;
	  WlzDVertex2	*nPos[3];

	  elm0 = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec,
	                                           nbrIdxBuf[idE]);
	  nPos[0] = &(elm0->edu[0].nod->pos);
	  nPos[1] = &(elm0->edu[1].nod->pos);
	  nPos[2] = &(elm0->edu[2].nod->pos);
	  cen.vtX = (nPos[0]->vtX + nPos[1]->vtX + nPos[2]->vtX) / 3.0;
	  cen.vtY = (nPos[0]->vtY + nPos[1]->vtY + nPos[2]->vtY) / 3.0;
	  WLZ_VTX_2_SUB(del, cen, pos);
	  d = WLZ_VTX_2_LENGTH(del);
	  if(d < WLZ_MESH_TOLERANCE)
	  {
	    setElmIdx = nbrIdxBuf[idE];
	    break;
	  }
	  else
	  {
	    nbrDstBuf[idE] = 1.0 / d;
	  }
	}
        /* Interpolate values. */
	if(setElmIdx >= 0)
	{
	  WlzGreyP 	x;

	  switch(ixv->vType)
	  {
	    case WLZ_GREY_LONG:
	      x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.lnp[ixi];
	      break;
	    case WLZ_GREY_INT:
	      x.inp = (int *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.inp[ixi];
	      break;
	    case WLZ_GREY_SHORT:
	      x.shp = (short *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.shp[ixi];
	      break;
	    case WLZ_GREY_UBYTE:
	      x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.ubp[ixi];
	      break;
	    case WLZ_GREY_FLOAT:
	      x.flp = (float *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.flp[ixi];
	      break;
	    case WLZ_GREY_DOUBLE:
	      x.dbp = (double *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.dbp[ixi];
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	else
	{
	  double  sum0 = 0.0,
		  sum1 = 0.0;

	  for(idE = 0; idE < nE; ++idE)
	  {
	    WlzGreyP 	x;

	    switch(ixv->vType)
	    {
	      case WLZ_GREY_LONG:
		x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.lnp[ixi];
		break;
	      case WLZ_GREY_INT:
		x.inp = (int *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.inp[ixi];
		break;
	      case WLZ_GREY_SHORT:
		x.shp = (short *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.shp[ixi];
		break;
	      case WLZ_GREY_UBYTE:
		x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.ubp[ixi];
		break;
	      case WLZ_GREY_FLOAT:
		x.flp = (float *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.flp[ixi];
		break;
	      case WLZ_GREY_DOUBLE:
		x.dbp = (double *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.dbp[ixi];
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	    sum0 += v * nbrDstBuf[idE];
	    sum1 += nbrDstBuf[idE];
	  }
	  v = sum0 / sum1;
	}
	switch(ixv->vType)
	{
	  case WLZ_GREY_LONG:
	    dst.lnp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_INT:
	    dst.inp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_SHORT:
	    dst.ubp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_UBYTE:
	    dst.ubp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_FLOAT:
	    dst.flp[idI] = v;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dst.dbp[idI] = v;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
    ++idI;
  }
  AlcFree(nbrIdxBuf);
  AlcFree(nbrDstBuf);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Sets values along the pixels of a single interval
* 		from the given mesh and it's indexed values. The pixels
* 		are set using linear interpolation.
* \param	dst				The interval values.
* \param	pl				Plane coordinate of the
* 						interval.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateElm3DLinear(WlzGreyP dst,
					int pl, int ln, int kolL, int kolR,
					WlzCMesh3D *mesh,
					WlzIndexedValues *ixv,
					int ixi)
{
  int		kl,
  		idE,
		idI,
		maxNbr = 0;
  int		*nbrIdxBuf = NULL;
  double	*nbrDstBuf = NULL;
  WlzDVertex3	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  pos.vtZ = pl;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    int		nE = 0;

    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos3D(mesh, idE, pos.vtX, pos.vtY, pos.vtZ,
    					0, NULL)) >= 0)
    {
      WlzCMeshElm3D *elm;
      int	    lastMaxNbr;

      lastMaxNbr = maxNbr;
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      /* Find elements connected to the given element. */
      nE = WlzCMeshElmRingElmIndices3D(elm, &maxNbr, &nbrIdxBuf, &errNum);
      /* Reallocate the distance buffer if required. */
      if((errNum == WLZ_ERR_NONE) && (lastMaxNbr < maxNbr))
      {
	nbrDstBuf = AlcRealloc(nbrDstBuf, maxNbr * sizeof(double));
        if(nbrDstBuf == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	int	setElmIdx = -1;
	double	v = 0.0;

        /* Compute inverse of distances to centroids of these elements from
	 * position. */
	for(idE = 0; idE < nE; ++idE)
	{
	  double	d;
	  WlzDVertex3	cen,
	  		del;
	  WlzCMeshNod3D *nod0;
	  WlzCMeshElm3D *elm0;
	  WlzDVertex3	*nPos[4];

	  elm0 = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec,
	                                           nbrIdxBuf[idE]);
	  nod0 = WLZ_CMESH_ELM3D_GET_NODE_0(elm0); nPos[0] = &(nod0->pos);
	  nod0 = WLZ_CMESH_ELM3D_GET_NODE_1(elm0); nPos[1] = &(nod0->pos);
	  nod0 = WLZ_CMESH_ELM3D_GET_NODE_2(elm0); nPos[2] = &(nod0->pos);
	  nod0 = WLZ_CMESH_ELM3D_GET_NODE_3(elm0); nPos[3] = &(nod0->pos);
	  cen.vtX = (nPos[0]->vtX + nPos[1]->vtX + 
	             nPos[2]->vtX + nPos[3]->vtX) / 4.0;
	  cen.vtY = (nPos[0]->vtY + nPos[1]->vtY +
	             nPos[2]->vtY + nPos[3]->vtY) / 4.0;
	  cen.vtZ = (nPos[0]->vtZ + nPos[1]->vtZ +
	             nPos[2]->vtZ + nPos[3]->vtZ) / 4.0;
	  WLZ_VTX_3_SUB(del, cen, pos);
	  d = WLZ_VTX_3_LENGTH(del);
	  if(d < WLZ_MESH_TOLERANCE)
	  {
	    setElmIdx = nbrIdxBuf[idE];
	    break;
	  }
	  else
	  {
	    nbrDstBuf[idE] = 1.0 / d;
	  }
	}
        /* Interpolate values. */
	if(setElmIdx >= 0)
	{
	  WlzGreyP x;

	  switch(ixv->vType)
	  {
	    case WLZ_GREY_LONG:
	      x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.lnp[ixi];
	      break;
	    case WLZ_GREY_INT:
	      x.inp = (int *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.inp[ixi];
	      break;
	    case WLZ_GREY_SHORT:
	      x.shp = (short *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.shp[ixi];
	      break;
	    case WLZ_GREY_UBYTE:
	      x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.ubp[ixi];
	      break;
	    case WLZ_GREY_FLOAT:
	      x.flp = (float *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.flp[ixi];
	      break;
	    case WLZ_GREY_DOUBLE:
	      x.dbp = (double *)WlzIndexedValueGet(ixv, setElmIdx);
	      v = x.dbp[ixi];
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
	else
	{
	  double  sum0 = 0.0,
		  sum1 = 0.0;

	  for(idE = 0; idE < nE; ++idE)
	  {
	    WlzGreyP x;

	    switch(ixv->vType)
	    {
	      case WLZ_GREY_LONG:
		x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.lnp[ixi];
		break;
	      case WLZ_GREY_INT:
		x.inp = (int *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.inp[ixi];
		break;
	      case WLZ_GREY_SHORT:
		x.shp = (short *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.shp[ixi];
		break;
	      case WLZ_GREY_UBYTE:
		x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.ubp[ixi];
		break;
	      case WLZ_GREY_FLOAT:
		x.flp = (float *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.flp[ixi];
		break;
	      case WLZ_GREY_DOUBLE:
		x.dbp = (double *)WlzIndexedValueGet(ixv, nbrIdxBuf[idE]);
		v = x.dbp[ixi];
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	    sum0 += v * nbrDstBuf[idE];
	    sum1 += nbrDstBuf[idE];
	  }
	  v = sum0 / sum1;
	}
	switch(ixv->vType)
	{
	  case WLZ_GREY_LONG:
	    dst.lnp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_INT:
	    dst.inp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_SHORT:
	    dst.ubp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_UBYTE:
	    dst.ubp[idI] = WLZ_NINT(v);
	    break;
	  case WLZ_GREY_FLOAT:
	    dst.flp[idI] = v;
	    break;
	  case WLZ_GREY_DOUBLE:
	    dst.dbp[idI] = v;
	    break;
	  default:
	    errNum = WLZ_ERR_GREY_TYPE;
	    break;
	}
      }
    }
    ++idI;
  }
  AlcFree(nbrIdxBuf);
  AlcFree(nbrDstBuf);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Interpolates values along the pixels of a single interval
* 		from the given mesh and it's indexed values. Uses nearest
* 		node interpolation within each mesh element.
* \param	dst				The interval values.
* \param	pl				Plane coordinate of the
* 						interval.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateNod3DNearest(WlzGreyP dst,
					int pl, int ln, int kolL, int kolR,
					WlzCMesh3D *mesh,
					WlzIndexedValues *ixv,
					int ixi)
{
  int		kl,
  		idE,
		idI,
		idN;
  WlzDVertex3	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  pos.vtZ = pl;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    WlzCMeshNod3D *nNod = NULL;

    pos.vtX = kl;
    /* Find nearest node within the element. */
    if((idE = WlzCMeshElmEnclosingPos3D(mesh, idE, pos.vtX, pos.vtY, pos.vtZ,
    					0, &idN)) >= 0)
    {
      int	idN;
      double    dMin;
      WlzDVertex3 del;
      WlzCMeshElm3D *elm;
      WlzCMeshNod3D *nod[4];

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
      nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
      nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
      nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
      nNod = nod[0];
      WLZ_VTX_3_SUB(del, pos, nod[0]->pos);
      dMin = WLZ_VTX_3_SQRLEN(del);
      for(idN = 1; idN < 4; ++idN)
      {
	double	d;

        WLZ_VTX_3_SUB(del, pos, nod[idN]->pos);
	d = WLZ_VTX_3_SQRLEN(del);
	if(d < dMin)
	{
	  nNod = nod[idN];
	  dMin = d;
	}
      }
    }
    else if(idN >= 0)
    {
      nNod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
    }
    /* Set value. */
    if(nNod)
    {
      WlzGreyP	x;

      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.lnp[idI] = x.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  x.inp = (int *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.inp[idI] = x.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  x.shp = (short *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.shp[idI] = x.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.ubp[idI] = x.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  x.flp = (float *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.flp[idI] = x.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  x.dbp = (double *)WlzIndexedValueGet(ixv, nNod->idx);
	  dst.dbp[idI] = x.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Interpolates values along the pixels of a single interval
* 		from the given mesh and it's indexed values. Uses linear
* 		interpolation within each mesh element.
* \param	dst				The interval values.
* \param	pl				Plane coordinate of the
* 						interval.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateNod3DLinear(WlzGreyP dst,
					int pl, int ln, int kolL, int kolR,
					WlzCMesh3D *mesh,
					WlzIndexedValues *ixv,
					int ixi)
{
  int		kl,
  		idE,
		idI,
		idN;
  WlzDVertex3	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  pos.vtZ = pl;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos3D(mesh, idE, pos.vtX, pos.vtY, pos.vtZ,
    					0, &idN)) >= 0)
    {
      int 	idC;
      double 	d[5];
      WlzCMeshElm3D *elm;
      WlzCMeshNod3D *nod[4];

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
      nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
      nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
      nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  for(idC = 0; idC < 4; ++idC)
	  {
            WlzLong 	*x;

	    x = (WlzLong *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[4] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, nod[3]->pos,
			d[0], d[1], d[2], d[3], pos);
	  dst.lnp[idI] = WLZ_NINT(d[4]);
	  break;
	case WLZ_GREY_INT:
	  for(idC = 0; idC < 4; ++idC)
	  {
            int 	*x;

	    x = (int *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[4] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, nod[3]->pos,
			d[0], d[1], d[2], d[3], pos);
	  dst.inp[idI] = WLZ_NINT(d[4]);
	  break;
	case WLZ_GREY_SHORT:
	  for(idC = 0; idC < 4; ++idC)
	  {
            short 	*x;

	    x = (short *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[4] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, nod[3]->pos,
			d[0], d[1], d[2], d[3], pos);
	  dst.shp[idI] = WLZ_NINT(d[4]);
	  break;
	case WLZ_GREY_UBYTE:
	  for(idC = 0; idC < 4; ++idC)
	  {
            WlzUByte 	*x;

	    x = (WlzUByte *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[4] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, nod[3]->pos,
			d[0], d[1], d[2], d[3], pos);
	  dst.ubp[idI] = WLZ_NINT(d[4]);
	  break;
	case WLZ_GREY_FLOAT:
	  for(idC = 0; idC < 4; ++idC)
	  {
            float 	*x;

	    x = (float *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  d[4] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			nod[2]->pos, nod[3]->pos,
			d[0], d[1], d[2], d[3], pos);
	  dst.flp[idI] = d[4];
	  break;
	case WLZ_GREY_DOUBLE:
	  for(idC = 0; idC < 4; ++idC)
	  {
            double 	*x;

	    x = (double *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	    d[idC] = x[ixi];
	  }
	  dst.dbp[idI] = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
			     nod[2]->pos, nod[3]->pos,
			     d[0], d[1], d[2], d[3], pos);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    else if(idN >= 0)
    {
      WlzCMeshNod3D *nod;

      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      switch(ixv->vType)
      {
	WlzGreyP x;

	case WLZ_GREY_LONG:
	  x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.lnp[idI] = x.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  x.inp = (int *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.inp[idI] = x.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  x.shp = (short *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.shp[idI] = x.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.ubp[idI] = x.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  x.flp = (float *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.flp[idI] = x.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  x.dbp = (double *)WlzIndexedValueGet(ixv, nod->idx);
	  dst.dbp[idI] = x.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Sets values along the pixels of a single interval
* 		from the given mesh and it's indexed values. The pixels
* 		are set to the first value attached to the mesh element.
* \param	dst				The interval values.
* \param	pl				Plane coordinate of the
* 						interval.
* \param	ln				Line coordinate of the
* 						interval.
* \param	kolL				leftmost column coordinate of
* 						the interval.
* \param	kolR				Rightmost column coordinate of
* 						the interval.
* \param	mesh				The mesh.
* \param	ixv				The indexed values.
* \param	ixi				Index into the indexed values.
*/
static WlzErrorNum WlzCMeshInterpolateElm3DNearest(WlzGreyP dst,
					int pl, int ln, int kolL, int kolR,
					WlzCMesh3D *mesh,
					WlzIndexedValues *ixv,
					int ixi)
{
  int		kl,
  		idE,
		idI;
  WlzDVertex3	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idI = 0;
  idE = -1;
  pos.vtY = ln;
  pos.vtZ = pl;
  for(kl = kolL; kl <= kolR; ++kl)
  {
    pos.vtX = kl;
    if((idE = WlzCMeshElmEnclosingPos3D(mesh, idE, pos.vtX, pos.vtY, pos.vtZ,
    					0, NULL)) >= 0)
    {
      WlzGreyP	x;
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  x.lnp = (WlzLong *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.lnp[idI] = x.lnp[ixi];
	  break;
	case WLZ_GREY_INT:
	  x.inp = (int *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.inp[idI] = x.inp[ixi];
	  break;
	case WLZ_GREY_SHORT:
	  x.shp = (short *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.shp[idI] = x.shp[ixi];
	  break;
	case WLZ_GREY_UBYTE:
	  x.ubp = (WlzUByte *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.ubp[idI] = x.ubp[ixi];
	  break;
	case WLZ_GREY_FLOAT:
	  x.flp = (float *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.flp[idI] = x.flp[ixi];
	  break;
	case WLZ_GREY_DOUBLE:
	  x.dbp = (double *)WlzIndexedValueGet(ixv, elm->idx);
	  dst.dbp[idI] = x.dbp[ixi];
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
    ++idI;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the transform coefficients for the given conforming
*		2D mesh scan element which must have a valid element index.
* \param	mObj			Conforming mesh transform object.
* \param	sE			Mesh scan element.
* \param	fwd			Non-zero for forward transform
*					mapping source to destination,
* 					otherwise inverse transform.
*/
static void 	WlzCMeshUpdateScanElm2D(WlzObject *mObj,
				        WlzCMeshScanElm2D *sE,
					int fwd)
{
  int		idN,
  		squash;
  double	areaSn2;
  double	*dsp;
  WlzDVertex2	dVx[3],
  		sVx[3];
  WlzCMesh2D	*mesh;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D	*elm;
  AlcVector	*vec;
  WlzIndexedValues *ixv;

  sE->flags = WLZ_CMESH_SCANELM_NONE;
  mesh = mObj->domain.cm2;
  vec = mesh->res.elm.vec;
  elm = (WlzCMeshElm2D *)AlcVectorItemGet(vec, sE->idx);
  ixv = mObj->values.x;
  if(ixv == NULL)
  {
    sE->trX[0] = 1.0;
    sE->trX[1] = 0.0;
    sE->trX[2] = 0.0;
    sE->trY[0] = 0.0;
    sE->trY[1] = 1.0;
    sE->trY[2] = 0.0;
  }
  else
  {
    for(idN = 0; idN < 3; ++idN)
    {
      nod = elm->edu[idN].nod;
      sVx[idN] = nod->pos;
      dsp = (double *)WlzIndexedValueGet(ixv, nod->idx);
      dVx[idN].vtX = dsp[0];
      dVx[idN].vtY = dsp[1];
      WLZ_VTX_2_ADD(dVx[idN], dVx[idN], sVx[idN]);
    }
    if(fwd)
    {
      sE->flags |= WLZ_CMESH_SCANELM_FWD;
      sE->flags &= ~WLZ_CMESH_SCANELM_REV;
      areaSn2 = WlzGeomTriangleSnArea2(sVx[0], sVx[1], sVx[2]);
      squash = WlzGeomTriangleAffineSolve(sE->trX, sE->trY, areaSn2,
					  sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
    }
    else
    {
      sE->flags |= WLZ_CMESH_SCANELM_REV;
      sE->flags &= ~WLZ_CMESH_SCANELM_FWD;
      areaSn2 = WlzGeomTriangleSnArea2(dVx[0], dVx[1], dVx[2]);
      squash = WlzGeomTriangleAffineSolve(sE->trX, sE->trY, areaSn2,
					  dVx, sVx, WLZ_MESH_TOLERANCE_SQ);
    }
    if(squash)
    {
      sE->flags |= WLZ_CMESH_SCANELM_SQUASH;
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the transform coefficients for the given conforming
*		2D5 mesh scan element which must have a valid element index.
* \param	mObj			Conforming mesh transform object.
* \param	sE			Mesh scan element.
* \param	fwd			Non-zero for forward transform
*					mapping source to destination,
* 					otherwise inverse transform.
*/
static void 	WlzCMeshUpdateScanElm2D5(WlzObject *mObj,
				         WlzCMeshScanElm3D *sE,
					 int fwd)
{
  int		idN,
  		squash;
  double	*dsp;
  WlzDVertex3	dVx[4],
  		sVx[4],
		tVx[3];
  WlzCMesh2D5	*mesh;
  WlzCMeshElm2D5 *elm;
  AlcVector	*vec;
  WlzIndexedValues *ixv;
  WlzCMeshNod2D5 *nod[3];

  sE->flags = WLZ_CMESH_SCANELM_NONE;
  mesh = mObj->domain.cm2d5;
  vec = mesh->res.elm.vec;
  elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(vec, sE->idx);
  ixv = mObj->values.x;
  if(ixv == NULL)
  {
    sE->tr[0] = 1.0;
    sE->tr[1] = 0.0;
    sE->tr[2] = 0.0;
    sE->tr[3] = 0.0;
    sE->tr[4] = 0.0;
    sE->tr[5] = 1.0;
    sE->tr[6] = 0.0;
    sE->tr[7] = 0.0;
    sE->tr[8] = 0.0;
    sE->tr[9] = 0.0;
    sE->tr[10] = 1.0;
    sE->tr[11] = 0.0;
  }
  else
  {
    nod[0] = WLZ_CMESH_ELM2D5_GET_NODE_0(elm);
    nod[1] = WLZ_CMESH_ELM2D5_GET_NODE_1(elm);
    nod[2] = WLZ_CMESH_ELM2D5_GET_NODE_2(elm);
    if(ixv == NULL)
    {
      for(idN = 0; idN < 3; ++idN)
      {
	sVx[idN] = dVx[idN] = nod[idN]->pos;
      }
    }
    else
    {
      for(idN = 0; idN < 3; ++idN)
      {
	sVx[idN] = nod[idN]->pos;
	dsp = (double *)WlzIndexedValueGet(ixv, nod[idN]->idx);
	dVx[idN].vtX = dsp[0];
	dVx[idN].vtY = dsp[1];
	dVx[idN].vtZ = dsp[2];
	WLZ_VTX_3_ADD(dVx[idN], dVx[idN], sVx[idN]);
      }
    }
    /* For an affine transform in 3D we need a tetrahedron rather than a
     * triangle; form one by taking the cross product of two edge vectors:
     * \f$v_i\f$ = n[i]->pos \f$\forall i\f$
     * \f$v_3 = ((v_1 - v_0) \times (v_2 - v_0)) + v_0 
     * n[3]->pos = v_3
     */
    WLZ_VTX_3_SUB(tVx[0], dVx[1], dVx[0]);
    WLZ_VTX_3_SUB(tVx[1], dVx[2], dVx[0]);
    WLZ_VTX_3_CROSS(tVx[2], tVx[0], tVx[1]);
    WLZ_VTX_3_ADD(dVx[3], tVx[2], dVx[0]);
    WLZ_VTX_3_SUB(tVx[0], sVx[1], sVx[0]);
    WLZ_VTX_3_SUB(tVx[1], sVx[2], sVx[0]);
    WLZ_VTX_3_CROSS(tVx[2], tVx[0], tVx[1]);
    WLZ_VTX_3_ADD(sVx[3], tVx[2], sVx[0]);
    if(fwd)
    {
      sE->flags |= WLZ_CMESH_SCANELM_FWD;
      sE->flags &= ~WLZ_CMESH_SCANELM_REV;
      squash = WlzGeomTetraAffineSolve(sE->tr, sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
    }
    else
    {
      sE->flags |= WLZ_CMESH_SCANELM_REV;
      sE->flags &= ~WLZ_CMESH_SCANELM_FWD;
      squash = WlzGeomTetraAffineSolve(sE->tr, dVx, sVx, WLZ_MESH_TOLERANCE_SQ);
    }
    if(squash)
    {
      sE->flags |= WLZ_CMESH_SCANELM_SQUASH;
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the transform coefficients for the given conforming
*		3D mesh scan element which must have a valid element index.
* \param	mObj			Conforming mesh transform object.
* \param	sE			Mesh scan element.
* \param	fwd			Non-zero for forward transform
*					mapping source to destination,
* 					otherwise inverse transform.
*/
static void 	WlzCMeshUpdateScanElm3D(WlzObject *mObj,
				        WlzCMeshScanElm3D *sE,
					int fwd)
{
  int		idN,
  		squash;
  double	*dsp;
  WlzDVertex3	dVx[4],
  		sVx[4];
  WlzCMesh3D	*mesh;
  WlzCMeshElm3D	*elm;
  AlcVector	*vec;
  WlzIndexedValues *ixv;
  WlzCMeshNod3D	*nod[4];

  sE->flags = WLZ_CMESH_SCANELM_NONE;
  mesh = mObj->domain.cm3;
  vec = mesh->res.elm.vec;
  elm = (WlzCMeshElm3D *)AlcVectorItemGet(vec, sE->idx);
  ixv = mObj->values.x;
  if(ixv == NULL)
  {
    sE->tr[0] = 1.0;
    sE->tr[1] = 0.0;
    sE->tr[2] = 0.0;
    sE->tr[3] = 0.0;
    sE->tr[4] = 0.0;
    sE->tr[5] = 1.0;
    sE->tr[6] = 0.0;
    sE->tr[7] = 0.0;
    sE->tr[8] = 0.0;
    sE->tr[9] = 0.0;
    sE->tr[10] = 1.0;
    sE->tr[11] = 0.0;
  }
  else
  {
    nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
    nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
    nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
    nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
    if(ixv == NULL)
    {
      for(idN = 0; idN < 4; ++idN)
      {
	sVx[idN] = dVx[idN] = nod[idN]->pos;
      }
    }
    else
    {
      for(idN = 0; idN < 4; ++idN)
      {
	sVx[idN] = nod[idN]->pos;
	dsp = (double *)WlzIndexedValueGet(ixv, nod[idN]->idx);
	dVx[idN].vtX = dsp[0];
	dVx[idN].vtY = dsp[1];
	dVx[idN].vtZ = dsp[2];
	WLZ_VTX_3_ADD(dVx[idN], dVx[idN], sVx[idN]);
      }
    }
    if(fwd)
    {
      sE->flags |= WLZ_CMESH_SCANELM_FWD;
      sE->flags &= ~WLZ_CMESH_SCANELM_REV;
      squash = WlzGeomTetraAffineSolve(sE->tr, sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
    }
    else
    {
      sE->flags |= WLZ_CMESH_SCANELM_REV;
      sE->flags &= ~WLZ_CMESH_SCANELM_FWD;
      squash = WlzGeomTetraAffineSolve(sE->tr, dVx, sVx, WLZ_MESH_TOLERANCE_SQ);
    }
    if(squash)
    {
      sE->flags |= WLZ_CMESH_SCANELM_SQUASH;
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Sets values in the destination object transforming those
*		of the source object.
* \param	dstObj			2D destination object with transformed
* 					domain but no values.
* \param	srcObj			2D source object.
* \param	mObj			Conforming mesh transform object.
* \param	interp			Level of interpolation.
*/
static WlzErrorNum WlzCMeshTransformValues2D(WlzObject *dstObj,
					WlzObject *srcObj,
					WlzObject *mObj,
					WlzInterpolationType interp)
{
  int		idP,
  		idX,
		iLft,
		iRgt,
		mItvIdx0,
		mItvIdx1,
		bufWidth = 0,
		itvWidth = 0;
  double	tD0 ,
  		tD1,
		tD2,
		tD3,
		tD4,
		trXX,
		trXYC,
		trYX,
		trYYC;
  int		*olpCnt = NULL;
  WlzGreyP	dGP,
  		olpBuf;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzPixelV	bgdV;
  WlzDVertex2	sPosD;
  WlzCMeshScanItv2D *mItv0 = NULL,
  		*mItv1 = NULL,
		*mItv2 = NULL;
  WlzCMeshScanWSp2D *mSWSp = NULL;
  WlzCMeshScanElm2D *sElm;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzGreyWSpace gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  olpBuf.inp = NULL;
  bgdV = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTableTypeToGreyType(srcObj->values.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&bgdV, bgdV, gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bufWidth = dstObj->domain.i->lastkl - dstObj->domain.i->kol1 + 1;
    errNum = WlzCMeshScanMakeOlpBufs(dstObj, gType,
                                     &olpBuf, &olpCnt, bufWidth);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItvIdx0 = 0;
    mSWSp = WlzCMeshScanWSpInit2D(mObj, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItv0 = mSWSp->itvs;
    errNum = WlzInitGreyScan(dstObj, &iWSp, &gWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
    while((errNum == WLZ_ERR_NONE) &&
	  (WlzNextGreyInterval(&iWSp) == 0))
    {
      dGP = gWSp.u_grintptr;
      itvWidth = iWSp.rgtpos - iWSp.lftpos + 1;
      switch(gType)
      {
	case WLZ_GREY_INT:   /* FALLTHROUGH */
	case WLZ_GREY_SHORT: /* FALLTHROUGH */
	case WLZ_GREY_UBYTE:
	  WlzValueSetInt(olpBuf.inp, 0, itvWidth);
	  break;
	case WLZ_GREY_FLOAT: /* FALLTHROUGH */
	case WLZ_GREY_DOUBLE:
	  WlzValueSetDouble(olpBuf.dbp, 0.0, itvWidth);
	  break;
	case WLZ_GREY_RGBA:
	  WlzValueSetInt(olpBuf.inp, 0, itvWidth);
	  WlzValueSetInt(olpBuf.inp + bufWidth, 0, itvWidth);
	  WlzValueSetInt(olpBuf.inp + (2 * bufWidth), 0, itvWidth);
	  WlzValueSetInt(olpBuf.inp + (3 * bufWidth), 0, itvWidth);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      WlzValueSetInt(olpCnt, 0, itvWidth);
      /* Update the mesh interval pointer so that it points to the first
       * mesh interval on the which intersects the current grey interval. */
      while((mItv0->line < iWSp.linpos) && (mItvIdx0 < mSWSp->nItvs))
      {
	++mItvIdx0;
	++mItv0;
      }
      while((mItv0->line <= iWSp.linpos) &&
	    (mItv0->rgtI < iWSp.lftpos) &&
	    (mItvIdx0 < mSWSp->nItvs))
      {
	++mItvIdx0;
	++mItv0;
      }
      if((mItv0->line == iWSp.linpos) &&
	 (iWSp.lftpos <= mItv0->rgtI) &&
	 (iWSp.rgtpos >= mItv0->lftI))
      {
	/* Mesh interval mItv0 intersects the current grey interval find
	 * the last mesh interval mItv1 which also intersects the current grey
	 * interval. */
	mItv1 = mItv0;
	mItvIdx1 = mItvIdx0;
	while((mItv1->line == iWSp.linpos) &&
	      (mItv1->lftI <= iWSp.rgtpos) &&
	      (mItvIdx1 < mSWSp->nItvs))
	{
	  ++mItvIdx1;
	  ++mItv1;
	}
	mItv2 = mItv1 - 1;
	mItv1 = mItv0;
	/* For each mesh interval which intersects the current grey interval. */
	while(mItv1 <= mItv2)
	{
	  /* Update mesh scanning. */
	  sElm = mSWSp->dElm + mItv1->elmIdx;
	  WlzCMeshUpdateScanElm2D(mSWSp->mTr, sElm, 0);
	  trXX = sElm->trX[0];
	  trXYC = (sElm->trX[1] * iWSp.linpos) + sElm->trX[2];
	  trYX = sElm->trY[0];
	  trYYC = (sElm->trY[1] * iWSp.linpos) + sElm->trY[2];
	  /* Find length of intersection and set the grey pointer. */
	  iLft = ALG_MAX(mItv1->lftI, iWSp.lftpos);
	  iRgt = ALG_MIN(mItv1->rgtI, iWSp.rgtpos);
	  idX = iLft;
	  switch(interp)
	  {
	    case WLZ_INTERPOLATION_NEAREST:
	      switch(gWSp.pixeltype)
	      {
		case WLZ_GREY_INT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += gVWSp->gVal[0].inv;
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_SHORT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += gVWSp->gVal[0].shv;
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += gVWSp->gVal[0].ubv;
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_FLOAT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += gVWSp->gVal[0].flv;
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_DOUBLE:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += gVWSp->gVal[0].dbv;
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_RGBA:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGet(gVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(gVWSp->bkdFlag == 0)
		    {
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += WLZ_RGBA_RED_GET(
					     gVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + bufWidth + idP) += WLZ_RGBA_GREEN_GET(
						    gVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + (2 * bufWidth) + idP) += WLZ_RGBA_BLUE_GET(
						    gVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + (3 * bufWidth) + idP) += WLZ_RGBA_ALPHA_GET(
						    gVWSp->gVal[0].rgbv);
		    }
		    ++idX;
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	      break;
	    case WLZ_INTERPOLATION_LINEAR:
	      switch(gWSp.pixeltype)
	      {
		case WLZ_GREY_INT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD0 = ((gVWSp->gVal[0]).inv * tD2 * tD3) +
			    ((gVWSp->gVal[1]).inv * tD0 * tD3) +
			    ((gVWSp->gVal[2]).inv * tD2 * tD1) +
			    ((gVWSp->gVal[3]).inv * tD0 * tD1);
		      tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.inp + idP) += gVWSp->gVal[0].inv;
		      }
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_SHORT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD0 = ((gVWSp->gVal[0]).shv * tD2 * tD3) +
			    ((gVWSp->gVal[1]).shv * tD0 * tD3) +
			    ((gVWSp->gVal[2]).shv * tD2 * tD1) +
			    ((gVWSp->gVal[3]).shv * tD0 * tD1);
		      tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.inp + idP) += gVWSp->gVal[0].shv;
		      }
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_UBYTE:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD0 = ((gVWSp->gVal[0]).ubv * tD2 * tD3) +
			    ((gVWSp->gVal[1]).ubv * tD0 * tD3) +
			    ((gVWSp->gVal[2]).ubv * tD2 * tD1) +
			    ((gVWSp->gVal[3]).ubv * tD0 * tD1);
		      tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.inp + idP) += gVWSp->gVal[0].ubv;
		      }
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_FLOAT:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD0 = ((gVWSp->gVal[0]).flv * tD2 * tD3) +
			    ((gVWSp->gVal[1]).flv * tD0 * tD3) +
			    ((gVWSp->gVal[2]).flv * tD2 * tD1) +
			    ((gVWSp->gVal[3]).flv * tD0 * tD1);
		      tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += WLZ_NINT(tD0);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.dbp + idP) += gVWSp->gVal[0].flv;
		      }
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_DOUBLE:
		  while(idX <= iRgt)
		  {
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD0 = ((gVWSp->gVal[0]).dbv * tD2 * tD3) +
			    ((gVWSp->gVal[1]).dbv * tD0 * tD3) +
			    ((gVWSp->gVal[2]).dbv * tD2 * tD1) +
			    ((gVWSp->gVal[3]).dbv * tD0 * tD1);
		      tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
		      idP = idX - iWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += WLZ_NINT(tD0);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.dbp + idP) += gVWSp->gVal[0].dbv;
		      }
		    }
		    ++idX;
		  }
		  break;
		case WLZ_GREY_RGBA:
		  while(idX <= iRgt)
		  {
		    idP = idX - iWSp.lftpos;
		    sPosD.vtX = (trXX * idX) + trXYC;
		    sPosD.vtY = (trYX * idX) + trYYC;
		    WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
		    if(gVWSp->bkdFlag == 0)
		    {
		      tD0 = sPosD.vtX - floor(sPosD.vtX);
		      tD1 = sPosD.vtY - floor(sPosD.vtY);
		      tD2 = 1.0 - tD0;
		      tD3 = 1.0 - tD1;
		      tD4 = (WLZ_RGBA_RED_GET((gVWSp->gVal[0]).rgbv) *
			     tD2 * tD3) +
			    (WLZ_RGBA_RED_GET((gVWSp->gVal[1]).rgbv) *
			     tD0 * tD3) +
			    (WLZ_RGBA_RED_GET((gVWSp->gVal[2]).rgbv) *
			     tD2 * tD1) +
			    (WLZ_RGBA_RED_GET((gVWSp->gVal[3]).rgbv) *
			     tD0 * tD1);
		      tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += WLZ_NINT(tD4);
		      tD4 = (WLZ_RGBA_GREEN_GET((gVWSp->gVal[0]).rgbv) *
			     tD2 * tD3) +
			    (WLZ_RGBA_GREEN_GET((gVWSp->gVal[1]).rgbv) *
			     tD0 * tD3) +
			    (WLZ_RGBA_GREEN_GET((gVWSp->gVal[2]).rgbv) *
			     tD2 * tD1) +
			    (WLZ_RGBA_GREEN_GET((gVWSp->gVal[3]).rgbv) *
			     tD0 * tD1);
		      tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
		      *(olpBuf.inp + bufWidth + idP) += WLZ_NINT(tD4);
		      tD4 = (WLZ_RGBA_BLUE_GET((gVWSp->gVal[0]).rgbv) *
			     tD2 * tD3) +
			    (WLZ_RGBA_BLUE_GET((gVWSp->gVal[1]).rgbv) *
			     tD0 * tD3) +
			    (WLZ_RGBA_BLUE_GET((gVWSp->gVal[2]).rgbv) *
			     tD2 * tD1) +
			    (WLZ_RGBA_BLUE_GET((gVWSp->gVal[3]).rgbv) *
			     tD0 * tD1);
		      tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
		      *(olpBuf.inp + (2 * bufWidth) + idP) += WLZ_NINT(tD4);
		      tD4 = (WLZ_RGBA_ALPHA_GET((gVWSp->gVal[0]).rgbv) *
			     tD2 * tD3) +
			    (WLZ_RGBA_ALPHA_GET((gVWSp->gVal[1]).rgbv) *
			     tD0 * tD3) +
			    (WLZ_RGBA_ALPHA_GET((gVWSp->gVal[2]).rgbv) *
			     tD2 * tD1) +
			    (WLZ_RGBA_ALPHA_GET((gVWSp->gVal[3]).rgbv) *
			     tD0 * tD1);
		      tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
		      *(olpBuf.inp + (3 * bufWidth) + idP) += WLZ_NINT(tD4);
		    }
		    else
		    {
		      WlzGreyValueGet(gVWSp, 0,
				      WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		      if(gVWSp->bkdFlag == 0)
		      {
			idP = idX - iWSp.lftpos;
			++*(olpCnt + idP);
			*(olpBuf.inp + idP) +=
				  WLZ_RGBA_RED_GET(gVWSp->gVal[0].rgbv);
			*(olpBuf.inp + bufWidth + idP) +=
				  WLZ_RGBA_GREEN_GET(gVWSp->gVal[0].rgbv);
			*(olpBuf.inp + (2 * bufWidth) + idP) +=
				  WLZ_RGBA_BLUE_GET(gVWSp->gVal[0].rgbv);
			*(olpBuf.inp + (3 * bufWidth) + idP) +=
				  WLZ_RGBA_ALPHA_GET(gVWSp->gVal[0].rgbv);
		      }
		    }
		    ++idX;
		  }
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	      break;
	    case WLZ_INTERPOLATION_CLASSIFY_1:     /* FALLTHROUGH */
	      errNum = WLZ_ERR_UNIMPLEMENTED;
	      break;
	    default:
	      errNum = WLZ_ERR_INTERPOLATION_TYPE;
	      break;
	  }
	  ++mItv1;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzCMeshScanFlushOlpBuf(dGP, olpBuf, olpCnt, bufWidth, bgdV,
					 iWSp.lftpos, iWSp.rgtpos,
					 interp, gType);
      }
    }
    (void )WlzEndGreyScan(&iWSp, &gWSp);
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  AlcFree(olpBuf.inp);
  AlcFree(olpCnt);
  WlzCMeshScanWSpFree2D(mSWSp);
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

/*!
* \return	New conforming mesh scan workspace.
* \ingroup	WlzTransform
* \brief	Allocate and initialise a 2D conforming mesh scan workspace.
* \param	mObj			Conforming mesh transform object.
* \param	trans			Build workspace for transformed
* 					mesh using the transform in the 
* 					mesh indexed values if non-zero.
* \param	dstErr			Destination error pointer.
*/
static WlzCMeshScanWSp2D *WlzCMeshScanWSpInit2D(WlzObject *mObj,
						int trans,
				    	        WlzErrorNum *dstErr)
{
  int		iIdx;
  unsigned int 	eIdx;
  double	ndLn,
  		eLnMin,
		eLnMax;
  double	*dsp;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D	*elm;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv = NULL;
  WlzCMeshEntRes *elmRes;
  WlzCMeshScanElm2D *dElm;
  WlzCMeshScanWSp2D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Check that the mesh transform object is appropriate to that it doesn't
   * need testing later. */
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((mesh = mObj->domain.cm2) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(trans != 0)
  {
    if((ixv = mObj->values.x) != NULL)
    {
      if(ixv->type != (WlzObjectType )WLZ_INDEXED_VALUES)
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if((ixv->rank != 1) || (ixv->dim[0] < 2) ||
	      (ixv->vType != WLZ_GREY_DOUBLE) ||
	      (ixv->attach != WLZ_VALUE_ATTACH_NOD))
      {
        errNum = WLZ_ERR_VALUES_DATA;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    elmRes = &(mesh->res.elm);
    if((mSWSp = (WlzCMeshScanWSp2D *)
		AlcCalloc(1, sizeof(WlzCMeshScanWSp2D))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp->mTr = mObj;
    /* Compute the maximum number of intervals in the displaced mesh. */
    eIdx = 0;
    for(eIdx = 0; eIdx < elmRes->maxEnt; ++eIdx)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(elmRes->vec, (size_t )eIdx);
      if(elm->idx >= 0)
      {
	nod = elm->edu[0].nod;
	if(ixv == NULL)
	{
	  eLnMin = eLnMax = nod->pos.vtY;
	  nod = elm->edu[1].nod;
	  ndLn = nod->pos.vtY;
	}
	else
	{
	  dsp = WlzIndexedValueGet(ixv, nod->idx);
	  eLnMin = eLnMax = nod->pos.vtY + dsp[1];
	  nod = elm->edu[1].nod;
	  dsp = WlzIndexedValueGet(ixv, nod->idx);
	  ndLn = nod->pos.vtY + dsp[1];
	}
	if(ndLn < eLnMin)
	{
	  eLnMin = ndLn;
	}
	else if(ndLn > eLnMax)
	{
	  eLnMax = ndLn;
	}
	nod = elm->edu[2].nod;
	if(ixv == NULL)
	{
	  ndLn = nod->pos.vtY;
	}
	else
	{
	  dsp = WlzIndexedValueGet(ixv, nod->idx);
	  ndLn = nod->pos.vtY + dsp[1];
	}
	if(ndLn < eLnMin)
	{
	  eLnMin = ndLn;
	}
	else if(ndLn > eLnMax)
	{
	  eLnMax = ndLn;
	}
	mSWSp->nItvs += WLZ_CMESH_POS_DTOI(eLnMax) -
	                WLZ_CMESH_POS_DTOI(eLnMin) + 1;
      }
    }
    if(((mSWSp->itvs = (WlzCMeshScanItv2D *)
    		       AlcMalloc(sizeof(WlzCMeshScanItv2D) *
				 mSWSp->nItvs)) == NULL) ||
       ((mSWSp->dElm = (WlzCMeshScanElm2D *)
		AlcCalloc(elmRes->maxEnt, sizeof(WlzCMeshScanElm2D))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill in the mesh scan intervals */
    eIdx = 0;
    iIdx = 0;
    dElm = mSWSp->dElm;
    while(eIdx < elmRes->maxEnt)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(elmRes->vec, (size_t )eIdx);
      dElm->idx = elm->idx;
      if(elm->idx >= 0)
      {
        iIdx += WlzCMeshScanTriElm2D(mSWSp, trans, elm, iIdx);
      }
      ++eIdx;
      ++dElm;
    }
    mSWSp->nItvs = iIdx;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Sort the conforming mesh scan intervals by line and then left
     * column */
    qsort(mSWSp->itvs, mSWSp->nItvs, sizeof(WlzCMeshScanItv2D),
          WlzCMeshItv2Cmp);
#ifdef WLZ_CMESHTRANSFORM_DEBUG
    for(iIdx = 0; iIdx < mSWSp->nItvs; ++iIdx)
    {
      (void )fprintf(stderr,
      		     "WlzCMeshScanWSpInit2D %d %d %d %d\n",
		     (mSWSp->itvs + iIdx)->elmIdx,
		     (mSWSp->itvs + iIdx)->lftI, (mSWSp->itvs + iIdx)->rgtI,
		     (mSWSp->itvs + iIdx)->line);
    }
#endif
  }
  else
  {
    WlzCMeshScanWSpFree2D(mSWSp);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mSWSp);
}

/*!
* \return	New conforming mesh scan workspace.
* \ingroup	WlzTransform
* \brief	Allocate and initialise a 3D conforming mesh scan workspace.
* \param	mObj			Conforming mesh transform object.
* \param	trans			Build workspace for transformed
* 					mesh using the transform in the 
* 					mesh indexed values if non-zero.
* \param	dstErr			Destination error pointer.
*/
static WlzCMeshScanWSp3D *WlzCMeshScanWSpInit3D(WlzObject *mObj, int trans,
				    	WlzErrorNum *dstErr)
{
  int		idE,
		idI,
		idN,
  		elmCnt,
		fstNod;
  double	*dsp;
  WlzDVertex3	dspP;
  WlzDVertex3	dspPos[4];
  WlzDBox3	dBox;
  AlcVector	*itvVec = NULL;
  AlcVector	*elmVec;
  WlzCMeshNod3D	*nodBuf[4];
  WlzCMeshElm3D	*elm;
  WlzIndexedValues *ixv = NULL;
  WlzCMesh3D	*mesh;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzCMeshScanElm3D *dElm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dBox.xMin = dBox.yMin = dBox.zMin = dBox.xMax = dBox.yMax = dBox.zMax = 0.0;
  /* Check that the mesh transform object is appropriate to that it doesn't
   * need testing later. */
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((mesh = mObj->domain.cm3) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(trans != 0)
  {
    if((ixv = mObj->values.x) != NULL)
    {
      if(ixv->type != (WlzObjectType )WLZ_INDEXED_VALUES)
      {
        errNum = WLZ_ERR_VALUES_TYPE;
      }
      else if((ixv->rank != 1) || (ixv->dim[0] < 3) ||
	      (ixv->vType != WLZ_GREY_DOUBLE) ||
	      (ixv->attach != WLZ_VALUE_ATTACH_NOD))
      {
        errNum = WLZ_ERR_VALUES_DATA;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    elmVec = mesh->res.elm.vec;
    elmCnt = mesh->res.elm.maxEnt;
    /* Collect the intervals in the displaced mesh. */
    /* Create temporary vector in which to accumulate the intervals. */
    if((itvVec = AlcVectorNew(1, sizeof(WlzCMeshScanItv3D),
				   mesh->res.elm.vec->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idE = 0;
    idI = 0;
    fstNod = 1;
    while((errNum == WLZ_ERR_NONE) && (idE < elmCnt))
    {
      /* Compute the displaced nodes and collect the intervals. */
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, (size_t )idE);
      if(elm->idx >= 0)
      {
	nodBuf[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
	nodBuf[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
	nodBuf[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	nodBuf[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
        for(idN = 0; idN < 4; ++idN)
	{
	  if(ixv == NULL)
	  {
	    dspP = nodBuf[idN]->pos;
	  }
	  else
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, nodBuf[idN]->idx);
	    dspP.vtX = nodBuf[idN]->pos.vtX + dsp[0];
	    dspP.vtY = nodBuf[idN]->pos.vtY + dsp[1];
	    dspP.vtZ = nodBuf[idN]->pos.vtZ + dsp[2];
	  }
	  dspPos[idN] = dspP;
	  if(fstNod)
	  {
	    dBox.xMin = dBox.xMax = dspP.vtX;
	    dBox.yMin = dBox.yMax = dspP.vtY;
	    dBox.zMin = dBox.zMax = dspP.vtZ;
	    fstNod = 0;
	  }
	  else
	  {
	    if(dspP.vtX < dBox.xMin)
	    {
	      dBox.xMin = dspP.vtX;
	    }
	    else if(dspP.vtX > dBox.xMax)
	    {
	      dBox.xMax = dspP.vtX;
	    }
	    if(dspP.vtY < dBox.yMin)
	    {
	      dBox.yMin = dspP.vtY;
	    }
	    else if(dspP.vtY > dBox.yMax)
	    {
	      dBox.yMax = dspP.vtY;
	    }
	    if(dspP.vtZ < dBox.zMin)
	    {
	      dBox.zMin = dspP.vtZ;
	    }
	    else if(dspP.vtZ > dBox.zMax)
	    {
	      dBox.zMax = dspP.vtZ;
	    }
	  }
	}
	errNum = WlzCMeshTetElmItv3D(itvVec, &idI, elm->idx, dspPos);
      }
      ++idE;
    }
  }
  /* Create a mesh scan workspace using the collected intervals. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp = WlzCMeshMakeScanWSp3D(mObj, idI, &errNum);
  }
  /* Copy the mesh scan intervals and sort them by plane, then line and
   * then column and then squeeze out the redundant intervals. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp->dBox.xMin = WLZ_CMESH_POS_DTOI(dBox.xMin) - 1;
    mSWSp->dBox.yMin = WLZ_CMESH_POS_DTOI(dBox.yMin) - 1;
    mSWSp->dBox.zMin = WLZ_CMESH_POS_DTOI(dBox.zMin) - 1;
    mSWSp->dBox.xMax = WLZ_CMESH_POS_DTOI(dBox.xMax) + 1;
    mSWSp->dBox.yMax = WLZ_CMESH_POS_DTOI(dBox.yMax) + 1;
    mSWSp->dBox.zMax = WLZ_CMESH_POS_DTOI(dBox.zMax) + 1;
    for(idI = 0; idI < mSWSp->nItvs; ++idI)
    {
      *(mSWSp->itvs + idI) = *(WlzCMeshScanItv3D *)
                             AlcVectorItemGet(itvVec, idI);
    }
    qsort(mSWSp->itvs, mSWSp->nItvs, sizeof(WlzCMeshScanItv3D),
          WlzCMeshItv3Cmp);
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, (size_t )idE);
      dElm = mSWSp->dElm + idE;
      memset(dElm, 0, sizeof(WlzCMeshScanElm3D));
      dElm->idx = elm->idx;
    }
#ifdef WLZ_CMESHTRANSFORM_DEBUG
    for(idI = 0; idI < mSWSp->nItvs; ++idI)
    {
      (void )fprintf(stderr,
      		     "WlzCMeshScanWSpInit3D %d %d %d %d %d\n",
		     (mSWSp->itvs + idI)->elmIdx,
		     (mSWSp->itvs + idI)->lftI, (mSWSp->itvs + idI)->rgtI,
		     (mSWSp->itvs + idI)->line, (mSWSp->itvs + idI)->plane);
    }
#endif
    WlzCMeshSqzRedundantItv3D(mSWSp);
#ifdef WLZ_CMESHTRANSFORM_DEBUG
    for(idI = 0; idI < mSWSp->nItvs; ++idI)
    {
      (void )fprintf(stderr,
      		     "WlzCMeshScanWSpInit3D %d %d %d %d %d\n",
		     (mSWSp->itvs + idI)->elmIdx,
		     (mSWSp->itvs + idI)->lftI, (mSWSp->itvs + idI)->rgtI,
		     (mSWSp->itvs + idI)->line, (mSWSp->itvs + idI)->plane);
    }
#endif
  }
  AlcVectorFree(itvVec);
  if(errNum != WLZ_ERR_NONE)
  {
    WlzCMeshScanWSpFree3D(mSWSp);
    mSWSp = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mSWSp);
}

/*!
* \ingroup	WlzTransform
* \brief	Squeezes out redundant intervals from the 3D conforming
*		mesh scan workspace. This replaces overlaping intervals
*		with a single interval.
* \param	mSWSp			Mesh scan workspace.
*/
static void	WlzCMeshSqzRedundantItv3D(WlzCMeshScanWSp3D *mSWSp)
{
  int		idx0,
  		idx1;
  WlzCMeshScanItv3D *itv0,
  		*itv1;

#ifdef WLZ_CMESHTRANSFORM_DEBUG
  (void )fprintf(stderr, "WlzCMeshSqzRedundantItv3D(E) nItvs = %d\n",
  		 mSWSp->nItvs);
#endif
  idx0 = idx1 = 0;
  itv0 = itv1 = mSWSp->itvs;
  while(idx1 < mSWSp->nItvs)
  {
    if((itv0->elmIdx == itv1->elmIdx) &&
       (itv0->plane == itv1->plane) &&
       (itv0->line == itv1->line) &&
       (itv0->rgtI >= itv1->lftI))
    {
      if(itv0->rgtI < itv1->rgtI)
      {
        itv0->rgtI = itv1->rgtI;
      }
    }
    else
    {
      ++idx0;
      *++itv0 = *itv1;
    }
    ++itv1;
    ++idx1;
  }
  mSWSp->nItvs = idx0 + 1;
#ifdef WLZ_CMESHTRANSFORM_DEBUG
  (void )fprintf(stderr, "WlzCMeshSqzRedundantItv3D(X) nItvs = %d\n",
  		 mSWSp->nItvs);
#endif
}

/*!
* \return	New uninitialised 3D conforming mesh scan workspace.
* \ingroup	WlzTransform
* \brief 	Creates a new uninitialised 3D conforming mesh scan workspace.
* \param	mObj			Conforming mesh transform object.
* \param	nItv			Number of mesh scan intervals.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshScanWSp3D *WlzCMeshMakeScanWSp3D(WlzObject *mObj,
						int nItv,
						WlzErrorNum *dstErr)
{
  WlzCMesh3D	*mesh;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mSWSp = (WlzCMeshScanWSp3D *)
  	      AlcCalloc(1, sizeof(WlzCMeshScanWSp3D))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    mSWSp->mTr = mObj;
    mSWSp->nItvs = nItv;
    mesh = mObj->domain.cm3;
    if(((mSWSp->itvs = (WlzCMeshScanItv3D *)
		       AlcMalloc(sizeof(WlzCMeshScanItv3D) *
			         mSWSp->nItvs)) == NULL) ||
       ((mSWSp->dElm = (WlzCMeshScanElm3D *)
		       AlcCalloc(mesh->res.elm.maxEnt,
			         sizeof(WlzCMeshScanElm3D))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
      WlzCMeshScanWSpFree3D(mSWSp);
      mSWSp = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mSWSp);
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Free's a conforming 2D mesh scan workspace.
* \param	mSnWSp			Conforming mesh scan workspace.
*/
static void	WlzCMeshScanWSpFree2D(WlzCMeshScanWSp2D *mSWSp)
{
  if(mSWSp)
  {
    if(mSWSp->itvs)
    {
      AlcFree(mSWSp->itvs);
    }
    if(mSWSp->dElm)
    {
      AlcFree(mSWSp->dElm);
    }
    AlcFree(mSWSp);
  }
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Free's a conforming 3D mesh scan workspace.
* \param	mSnWSp			Conforming mesh scan workspace.
*/
static void	WlzCMeshScanWSpFree3D(WlzCMeshScanWSp3D *mSWSp)
{
  if(mSWSp)
  {
    if(mSWSp->itvs)
    {
      AlcFree(mSWSp->itvs);
    }
    if(mSWSp->dElm)
    {
      AlcFree(mSWSp->dElm);
    }
    AlcFree(mSWSp);
  }
}

/*!
* \return	Number of intervals added from the conforming mesh
*		scan element.
* \ingroup	WlzTransform
* \brief	Scans a single triangular conforming mesh scan element
*		into conforming mesh intervals.
* \param	mSnWSp			Conforming mesh scan workspace.
* \param	trans			Use transformed element if non zero.
* \param	elm			Mesh element which is assumed to be
*					valid.
* \param	iIdx			Conforming mesh element interval index.
*/
static int	WlzCMeshScanTriElm2D(WlzCMeshScanWSp2D *mSWSp, int trans,
				     WlzCMeshElm2D *elm, int iIdx)
{
  int		kolI0,
		kolI1,
  		lineF,
  		lineI,
		lineL,
		ndIdx0,
  		ndIdx1,
		ndIdx2,
		iCnt = 0;
  double	x0,
  		x1,
		tD0,
  		tD1;
  double	*dsp;
  double	inc[3];
  WlzIVertex2	dNd[3],
  	 	sNd[3];
  WlzCMeshNod2D	*nod;
  WlzDVertex2	dVx0;
  WlzIndexedValues *ixv;
  WlzCMeshScanItv2D *itv;

  ixv = (trans != 0)? mSWSp->mTr->values.x: NULL;
  /* Compute the integer displaced nodes of the element. */
  for(ndIdx0 = 0; ndIdx0 < 3; ++ndIdx0)
  {
    nod = elm->edu[ndIdx0].nod;
    dVx0 = nod->pos;
    if(ixv == NULL)
    {
      tD0 = dVx0.vtX;
      tD1 = dVx0.vtY;
    }
    else
    {
      dsp = (double *)WlzIndexedValueGet(ixv, nod->idx);
      tD0 = dVx0.vtX + dsp[0];
      tD1 = dVx0.vtY + dsp[1];
    }
    dNd[ndIdx0].vtX = WLZ_CMESH_POS_DTOI(tD0);
    dNd[ndIdx0].vtY = WLZ_CMESH_POS_DTOI(tD1);
  }
  /* Sort nodes by line coordinate, min == 0, mid == 1, max == 2. */
  if(dNd[0].vtY < dNd[1].vtY)
  {
    ndIdx0 = (dNd[0].vtY < dNd[2].vtY)? 0: 2;
  }
  else
  {
    ndIdx0 = (dNd[1].vtY < dNd[2].vtY)? 1: 2;
  }
  ndIdx1 = (ndIdx0 + 1) % 3;
  ndIdx2 = (ndIdx0 + 2) % 3;
  if(dNd[ndIdx2].vtY < dNd[ndIdx1].vtY)
  {
    ndIdx1 = ndIdx2;
    ndIdx2 = (ndIdx0 + 1) % 3;
  }
  sNd[0] = dNd[ndIdx0];
  sNd[1] = dNd[ndIdx1];
  sNd[2] = dNd[ndIdx2];
  /* Compute deltas. */
  WLZ_VTX_2_SUB(dNd[0], sNd[0], sNd[1]);
  WLZ_VTX_2_SUB(dNd[1], sNd[1], sNd[2]);
  WLZ_VTX_2_SUB(dNd[2], sNd[2], sNd[0]);
  /* Classify the triangle and then scan convert it. */
  itv = mSWSp->itvs + iIdx;
  if(fabs(dNd[2].vtY) == 0)
  {
    /* Single horizontal interval
     * *-*-* */
    if(sNd[0].vtX <= sNd[1].vtX)
    {
      if(sNd[2].vtX <= sNd[0].vtX)
      {
	kolI0 = sNd[2].vtX;
	kolI1 = WLZ_MAX(sNd[0].vtX, sNd[1].vtX);
      }
      else /* sNd[2].vtX > Nd[0].vtX */
      {
	kolI0 = sNd[0].vtX;
	kolI1 = WLZ_MAX(sNd[1].vtX, sNd[2].vtX);
      }
    }
    else /* sNd[0].vtX > sNd[1].vtX */
    {
      if(sNd[2].vtX <= sNd[1].vtX)
      {
	kolI0 = sNd[2].vtX;
	kolI1 = WLZ_MAX(sNd[0].vtX, sNd[1].vtX);
      }
      else
      {
	kolI0 = sNd[1].vtX;
	kolI1 = WLZ_MAX(sNd[2].vtX, sNd[0].vtX);
      }
    }
    itv->elmIdx = elm->idx;
    itv->line = sNd[0].vtY;
    itv->lftI = kolI0;
    itv->rgtI = kolI1;
    iCnt = 1;
  }
  else if((dNd[0].vtX == 0) && (dNd[1].vtX == 0))
  {
    /* Many possible single column intervals
     * 0
     * |
     * 1
     * |
     * 2 */
    kolI0 = sNd[0].vtX;
    lineI = WLZ_CMESH_POS_DTOI(sNd[0].vtY);
    lineL = WLZ_CMESH_POS_DTOI(sNd[2].vtY);
    while(lineI <= lineL)
    {
      itv->elmIdx = elm->idx;
      itv->line = lineI;
      itv->lftI = itv->rgtI = kolI0;
      ++itv;
      ++lineI;
      ++iCnt;
    }
  }
  else
  {
    /* General case for triangles with non-zero area. */
    itv = mSWSp->itvs + iIdx;
    lineF = sNd[0].vtY;
    lineL = sNd[2].vtY;
    inc[0] = (dNd[0].vtY != 0)?
             (double )(dNd[0].vtX) / (double )(dNd[0].vtY): 0.0;
    inc[1] = (dNd[1].vtY != 0)?
             (double )(dNd[1].vtX) / (double )(dNd[1].vtY): 0.0;
    inc[2] = (dNd[2].vtY != 0)?
             (double )(dNd[2].vtX) / (double )(dNd[2].vtY): 0.0;
    lineI = lineF;
    while(lineI <= lineL)
    {
      if(lineI == lineF)
      {
        x0 = sNd[0].vtX;
	x1 = (lineF + 1 > sNd[1].vtY)?  sNd[1].vtX: sNd[0].vtX;
      }
      else
      {
	x0 = sNd[0].vtX + inc[2] * (lineI - sNd[0].vtY);
        x1 = (lineI >= sNd[1].vtY)? sNd[1].vtX + inc[1] * (lineI - sNd[1].vtY):
	                            sNd[0].vtX + inc[0] * (lineI - sNd[0].vtY);
      }
      if(x0 > x1)
      {
        tD0 = x0; x0 = x1; x1 = tD0;
      }
      kolI0 = WLZ_CMESH_POS_DTOI(x0 + WLZ_MESH_TOLERANCE);
      kolI1 = WLZ_CMESH_POS_DTOI(x1 + WLZ_MESH_TOLERANCE);
      if(kolI0 <= kolI1)
      {
	itv->elmIdx = elm->idx;
        itv->line = lineI;
	itv->lftI = kolI0;
	itv->rgtI = kolI1;
	++itv;
	++iCnt;
      }
      ++lineI;
    }
  }
  return(iCnt);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the intervals which are intersected by a
*		tetrahedron with the given vertex positions.
*		The vertices are sorted by z, then y, then x and sweept
*		through.
* \param	itvVec			Vector in which to accumulate the
* 					intervals.
* \param	idI			On entry this is the current
*					vector index and on return it is
*					the updated index.
* \param	elmIdx			Index of the element corresponding
*					to the tetrahedron.
* \param	vtx			Array of four vertex positions,
*					these are sorted in place by this
*					function.
*/
static WlzErrorNum WlzCMeshTetElmItv3D(AlcVector *itvVec, int *idI,
				       int elmIdx, WlzDVertex3 *vtx)
{
  int		isnCnt;
  double	a,
  		pl,
		plL,
		plU;
  WlzDVertex3	del10,
                del20,
		del30,
		del21,
		del31,
		del32;
  WlzDVertex3	isn[5];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-10;

  /* TODO Optimize when tested. */
  /* Reorder the vertices so that they are sortied by z then y then x. */
  qsort(vtx, 4, sizeof(WlzDVertex3), WlzCMeshDVertex3Cmp);
#ifdef WLZ_CMESHTRANSFORM_DEBUG
  (void )fprintf(stderr,
                 "WlzCMeshTetElmItv3D %d "
		 "{%g,%g,%g},{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}\n",
		 elmIdx, 
		 vtx[0].vtX, vtx[0].vtY, vtx[0].vtZ,
		 vtx[1].vtX, vtx[1].vtY, vtx[1].vtZ,
		 vtx[2].vtX, vtx[2].vtY, vtx[2].vtZ,
		 vtx[3].vtX, vtx[3].vtY, vtx[3].vtZ);
#endif
  /* Sweep through the tetrahedron. */
  pl = WLZ_CMESH_POS_DTOI(vtx[0].vtZ);
  if(pl < vtx[3].vtZ + tol)
  {
    WLZ_VTX_3_SUB(del10, vtx[1], vtx[0]);
    WLZ_VTX_3_SUB(del20, vtx[2], vtx[0]);
    WLZ_VTX_3_SUB(del30, vtx[3], vtx[0]);
    WLZ_VTX_3_SUB(del21, vtx[2], vtx[1]);
    WLZ_VTX_3_SUB(del31, vtx[3], vtx[1]);
    WLZ_VTX_3_SUB(del32, vtx[3], vtx[2]);
    do
    {
      isnCnt = 0;
      plL = pl - tol;
      plU = pl + tol;
      if(plL < vtx[0].vtZ)
      {
	isn[isnCnt++] = vtx[0];
      }
      if(plL < vtx[1].vtZ)
      {
	if(plU > vtx[1].vtZ)
	{
	  isn[isnCnt++] = vtx[1];
	}
	else
	{
	  if((del10.vtZ > tol) &&
	     ((a = (pl - vtx[0].vtZ) / del10.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[0].vtX + (a * del10.vtX);
	    isn[isnCnt].vtY = vtx[0].vtY + (a * del10.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	}
      }
      if(plL < vtx[2].vtZ)
      {
	if(plU > vtx[2].vtZ)
	{
	  isn[isnCnt++] = vtx[2];
	}
	else
	{
	  if((del20.vtZ > tol) &&
	     ((a = (pl - vtx[0].vtZ) / del20.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[0].vtX + (a * del20.vtX);
	    isn[isnCnt].vtY = vtx[0].vtY + (a * del20.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del21.vtZ > tol) &&
	     ((a = (pl - vtx[1].vtZ) / del21.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[1].vtX + (a * del21.vtX);
	    isn[isnCnt].vtY = vtx[1].vtY + (a * del21.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	}
      }
      if(plL < vtx[3].vtZ)
      {
	if(plU > vtx[3].vtZ)
	{
	  isn[isnCnt++] = vtx[3];
	}
	else
	{
	  if((del30.vtZ > tol) &&
	     ((a = (pl - vtx[0].vtZ) / del30.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[0].vtX + (a * del30.vtX);
	    isn[isnCnt].vtY = vtx[0].vtY + (a * del30.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del31.vtZ > tol) &&
	     ((a = (pl - vtx[1].vtZ) / del31.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[1].vtX + (a * del31.vtX);
	    isn[isnCnt].vtY = vtx[1].vtY + (a * del31.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del32.vtZ > tol) &&
	     ((a = (pl - vtx[2].vtZ) / del32.vtZ) > tol) &&
	     (a < 1.0 - tol))
	  {
	    isn[isnCnt].vtX = vtx[2].vtX + (a * del32.vtX);
	    isn[isnCnt].vtY = vtx[2].vtY + (a * del32.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	}
      }
      switch(isnCnt)
      {
        case 1:
	  isn[1] = isn[0];
	  errNum = WlzCMeshAddItv3D(itvVec, idI, elmIdx, isn);
	  break;
        case 2:
	  errNum = WlzCMeshAddItv3D(itvVec, idI, elmIdx, isn);
	  break;
        case 3:
	  errNum = WlzCMeshTriElmItv3D(itvVec, idI, elmIdx, isn);
	  break;
        case 4:
	  errNum = WlzCMeshTriElmItv3D(itvVec, idI, elmIdx, isn);
	  errNum = WlzCMeshTriElmItv3D(itvVec, idI, elmIdx, isn + 1);
	  break;
      }
      pl += 1.0;
    } while((errNum == WLZ_ERR_NONE) &&
            (pl < vtx[3].vtZ + tol));
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the intervals which are intersected by a
*		triangle in 3D space which lies on a plane parallel
*		to the x-y axis with the given vertex positions.
*		The vertices are sorted by y and sweept through.
*		these regions are used to control the sweep.
* \param	itvVec			Vector in which to accumulate the
* 					intervals.
* \param	idI			On entry this is the current
*					vector index and on return it is
*					the updated index.
* \param	elmIdx			Index of the element containing
*					to the triangle.
* \param	vtx			Array of three vertex positions,
*					these are sorted in place by this
*					function.
*/
static WlzErrorNum WlzCMeshTriElmItv3D(AlcVector *itvVec, int *idI,
				       int elmIdx, WlzDVertex3 *vtx)
{
  int		cnt,
		lnCnt,
		idI0,
  		idV0,
		idV1,
		idV2,
		klI,
		lnI,
		plI;
  double	inc,
  		klD;
  int		dX[3],
  		dY[3];
  WlzIVertex2	sVtx[3];
  WlzCMeshScanItv3D *itv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef WLZ_CMESHTRANSFORM_DEBUG
  (void )fprintf(stderr,
                 "WlzCMeshTriElmItv3D %d "
		 "{%g,%g,%g},{%g,%g,%g},{%g,%g,%g}\n",
		 elmIdx, 
		 vtx[0].vtX, vtx[0].vtY, vtx[0].vtZ,
		 vtx[1].vtX, vtx[1].vtY, vtx[1].vtZ,
		 vtx[2].vtX, vtx[2].vtY, vtx[2].vtZ);
#endif
  /* Sort vertices by line coordinate (min == 0, mid == 1, max == 2)
   * and convert to integer. */
  if(vtx[0].vtY < vtx[1].vtY)
  {
    idV0 = (vtx[0].vtY < vtx[2].vtY)? 0: 2;
  }
  else
  {
    idV0 = (vtx[1].vtY < vtx[2].vtY)? 1: 2;
  }
  idV1 = (idV0 + 1) % 3;
  idV2 = (idV0 + 2) % 3;
  if(vtx[idV2].vtY < vtx[idV1].vtY)
  {
    idV1 = idV2;
    idV2 = (idV0 + 1) % 3;
  }
  sVtx[0].vtX = WLZ_CMESH_POS_DTOI(vtx[idV0].vtX);
  sVtx[0].vtY = WLZ_CMESH_POS_DTOI(vtx[idV0].vtY);
  sVtx[1].vtX = WLZ_CMESH_POS_DTOI(vtx[idV1].vtX);
  sVtx[1].vtY = WLZ_CMESH_POS_DTOI(vtx[idV1].vtY);
  sVtx[2].vtX = WLZ_CMESH_POS_DTOI(vtx[idV2].vtX);
  sVtx[2].vtY = WLZ_CMESH_POS_DTOI(vtx[idV2].vtY);
  /* Compute deltas. */
  dX[0] = sVtx[0].vtX - sVtx[1].vtX;
  dX[1] = sVtx[1].vtX - sVtx[2].vtX;
  dY[2] = sVtx[2].vtY - sVtx[0].vtY;
  /* If the triangles vertices are not coincident scan convert it. */
  if(dY[2] && (dX[0] || dX[1]))
  {
    plI = WLZ_CMESH_POS_DTOI(vtx[0].vtZ);
    dY[0] = sVtx[0].vtY - sVtx[1].vtY;
    dY[1] = sVtx[1].vtY - sVtx[2].vtY;
    dX[2] = sVtx[2].vtX - sVtx[0].vtX;
    /* Nodes: min -> max */
    idI0 = *idI;
    klD = sVtx[0].vtX;
    lnI = sVtx[0].vtY;
    inc = (double )(dX[2]) / (double )(dY[2]);
    lnCnt = sVtx[2].vtY - sVtx[0].vtY + 1;
    if(AlcVectorExtend(itvVec, *idI + lnCnt) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      cnt = lnCnt;
      while(cnt-- > 0)
      {
	itv = (WlzCMeshScanItv3D *)AlcVectorItemGet(itvVec, idI0);
	itv->elmIdx = elmIdx;
	itv->plane = plI;
	itv->line = lnI++;
	itv->lftI = itv->rgtI = WLZ_CMESH_POS_DTOI(klD);
	klD += inc;
	++idI0;
      }
      if(dY[0] != 0)
      {
	/* Nodes: mid -> min */
	idI0 = *idI;
	klD = sVtx[0].vtX;
	inc = (double )(dX[0]) / (double )(dY[0]);
	cnt = sVtx[1].vtY - sVtx[0].vtY + 1;
	while(cnt-- > 0)
	{
	  klI = WLZ_CMESH_POS_DTOI(klD);
	  itv = (WlzCMeshScanItv3D *)AlcVectorItemGet(itvVec, idI0);
	  if(klI > itv->lftI)
	  {
	    itv->rgtI = klI;
	  }
	  else
	  {
	    itv->lftI = klI;
	  }
	  klD += inc;
	  ++idI0;
	}
      }
      if(dY[1])
      {
	/* Nodes: max -> mid */
	idI0 = *idI + lnCnt - 1;
	klD = sVtx[2].vtX;
	inc = (double )(dX[1]) / (double )(dY[1]);
	cnt = sVtx[2].vtY - sVtx[1].vtY + 1;
	while(cnt-- > 0)
	{
	  klI = WLZ_CMESH_POS_DTOI(klD);
	  itv = (WlzCMeshScanItv3D *)AlcVectorItemGet(itvVec, idI0);
	  if(klI < itv->lftI)
	  {
	    itv->lftI = klI;
	  }
	  else if(klI > itv->rgtI)
	  {
	    itv->rgtI = klI;
	  }
	  klD -= inc;
	  --idI0;
	}
      }
#ifdef WLZ_CMESHTRANSFORM_DEBUG
      for(idI0 = *idI; idI0 < *idI + lnCnt; ++idI0)
      {
        itv = (WlzCMeshScanItv3D *)AlcVectorExtendAndGet(itvVec, idI0);
	(void )fprintf(stderr,
		       "WlzCMeshTriElmItv3D %d %d %d %d %d\n",
		       itv->elmIdx, itv->lftI, itv->rgtI,
		       itv->line, itv->plane);
      }
#endif
      *idI += lnCnt;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds a 3D scan intervals to the existing interval vector.
* \param	itvVec			The interval vector.
* \param	idI			On entry this is the current
*					vector index and on return it is
*					the updated index.
* \param	elmIdx			Index of the element corresponding
*					to the interval.
* \param	vtx			Array of two vertex positions
*					which define the start and end of
*					interval.
*/
static WlzErrorNum WlzCMeshAddItv3D(AlcVector *itvVec, int *idI,
				    int elmIdx, WlzDVertex3 *vtx)
{
  WlzCMeshScanItv3D *itv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO Optimize when tested. */
  if((itv = (WlzCMeshScanItv3D *)AlcVectorExtendAndGet(itvVec, *idI)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    itv->elmIdx = elmIdx;
    itv->line = WLZ_CMESH_POS_DTOI(vtx[0].vtY);
    itv->plane = WLZ_CMESH_POS_DTOI(vtx[0].vtZ);
    if(vtx[0].vtX < vtx[1].vtX)
    {
      itv->lftI = WLZ_CMESH_POS_DTOI(vtx[0].vtX);
      itv->rgtI = WLZ_CMESH_POS_DTOI(vtx[1].vtX);
    }
    else
    {
      itv->lftI = WLZ_CMESH_POS_DTOI(vtx[1].vtX);
      itv->rgtI = WLZ_CMESH_POS_DTOI(vtx[0].vtX);
    }
    ++*idI;
  }
#ifdef WLZ_CMESHTRANSFORM_DEBUG
  (void )fprintf(stderr,
                 "WlzCMeshAddItv3D %d %d %d %d %d\n",
		 itv->elmIdx, itv->lftI, itv->rgtI, itv->line, itv->plane);
#endif
  return(errNum);
}

/*!
* \return	Sorting value for qsort.
* \ingroup	WlzTransform
* \brief	Callback function for qsort(3) to sort 3D double
*		vertices by z then y then x.
* \param	cmp0			Used to pass first mesh interval.
* \param	cmp1			Used to pass second mesh interval.
*/
static int	WlzCMeshDVertex3Cmp(const void *cmp0, const void *cmp1)
{
  int		rtn = 0;
  double	tst;
  WlzDVertex3	*v0,
  		*v1;

  v0 = (WlzDVertex3 *)cmp0;
  v1 = (WlzDVertex3 *)cmp1;
  if((tst = v0->vtZ - v1->vtZ) < 0.0)
  {
    rtn = -1;
  }
  else if(tst > 0.0)
  {
    rtn = 1;
  }
  else
  {
    if((tst = v0->vtY - v1->vtY) < 0.0)
    {
      rtn = -1;
    }
    else if(tst > 0.0)
    {
      rtn = 1;
    }
    else
    {
      if((tst = v0->vtX - v1->vtX) < 0.0)
      {
	rtn = -1;
      }
      else if(tst > 0.0)
      {
	rtn = 1;
      }
    }
  }
  return(rtn);
}

/*!
* \return	Sorting value for qsort.
* \ingroup	WlzTransform
* \brief	Callback function for qsort(3) to sort 2D conforming mesh
*		element intervals by line and then left column.
* \param	cmp0			Used to pass first mesh interval.
* \param	cmp1			Used to pass second mesh interval.
*/
static int	WlzCMeshItv2Cmp(const void *cmp0, const void *cmp1)
{
  int		rtn;
  WlzCMeshScanItv2D *itv0,
  		 *itv1;

  itv0 = (WlzCMeshScanItv2D *)cmp0;
  itv1 = (WlzCMeshScanItv2D *)cmp1;
  if((rtn = (itv0->line - itv1->line)) == 0)
  {
    if((rtn = itv0->lftI - itv1->lftI) == 0)
    {
      rtn = itv1->rgtI - itv0->rgtI;
    }
  }
  return(rtn);
}

/*!
* \return	Sorting value for qsort.
* \ingroup	WlzTransform
* \brief	Callback function for qsort(3) to sort 3D conforming mesh
*		element intervals by plane, line and then left column.
* \param	cmp0			Used to pass first mesh interval.
* \param	cmp1			Used to pass second mesh interval.
*/
static int	WlzCMeshItv3Cmp(const void *cmp0, const void *cmp1)
{
  int		rtn;
  WlzCMeshScanItv3D *itv0,
  		 *itv1;

  itv0 = (WlzCMeshScanItv3D *)cmp0;
  itv1 = (WlzCMeshScanItv3D *)cmp1;
  if((rtn = (itv0->plane - itv1->plane)) == 0)
  {
    if((rtn = (itv0->line - itv1->line)) == 0)
    {
      if((rtn = itv0->lftI - itv1->lftI) == 0)
      {
	rtn = itv0->rgtI - itv1->rgtI;
      }
    }
  }
  return(rtn);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a conforming mesh transform to the given source
*		object.
* \param	srcObj			Object to be transformed.
* \param	mObj			Conforming mesh transform object.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzCMeshTransformObj(WlzObject *srcObj,
				     WlzObject *mObj,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzPixelV	bgdV;
  WlzDomain	dstDom;
  WlzValues	srcValues,
  		dstValues;
  WlzIndexedValues *mIxv;
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
  dstDom.core = NULL;
  dstValues.core = NULL;
  srcValues.core = NULL;
  if((srcObj == NULL) || (mObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((srcObj->domain.core == NULL) || (mObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mIxv = mObj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(mIxv->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((mIxv->rank != 1) || (mIxv->vType != WLZ_GREY_DOUBLE) ||
	  (mIxv->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(mObj->type)
    {
      case WLZ_CMESH_2D:
	if(mObj->domain.core->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(mIxv->dim[0] < 2)
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
        break;
      case WLZ_CMESH_2D5:
	if(mObj->domain.core->type != WLZ_CMESH_2D5)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(mIxv->dim[0] < 3)
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
        break;
      case WLZ_CMESH_3D:
	if(mObj->domain.core->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(mIxv->dim[0] < 3)
	{
	  errNum = WLZ_ERR_VALUES_DATA;
	}
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
	dstObj = WlzMakeEmpty(&errNum);
	break;
      case WLZ_CMESH_2D:  /* FALLTHROUGH */
      case WLZ_CMESH_2D5: /* FALLTHROUGH */
      case WLZ_CMESH_3D:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(srcObj->type)
	  {
	    case WLZ_CMESH_2D:
	      dstDom.cm2 = WlzCMeshTransformCMesh2D(srcObj->domain.cm2,
	      					    mObj, 1, &errNum);
	      break;
	    case WLZ_CMESH_2D5:
	      dstDom.cm2d5 = WlzCMeshTransformCMesh2D5(srcObj->domain.cm2d5,
	      					       mObj, 1, &errNum);
	      break;
	    case WLZ_CMESH_3D:
	      dstDom.cm3 = WlzCMeshTransformCMesh3D(srcObj->domain.cm3,
	      					    mObj, 1, &errNum);
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzMakeMain(srcObj->type, dstDom, srcValues,
				 NULL, NULL, &errNum);
	  }
	  if((errNum != WLZ_ERR_NONE) && dstDom.core)
	  {
	    (void )WlzFreeDomain(dstDom);
	  }
	}
	break;
      case WLZ_CONTOUR:
        if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  dstDom.ctr = WlzCMeshTransformContour(srcObj->domain.ctr,
	                                         mObj, 1, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzMakeMain(srcObj->type, dstDom,  srcValues,
				 NULL, NULL, &errNum);
	                       
	}
	if((errNum != WLZ_ERR_NONE) && dstDom.core)
	{
	  (void )WlzFreeDomain(dstDom);
	}
	break;
      case WLZ_2D_POLYGON: /* FALLTHROUGH */
      case WLZ_BOUNDLIST:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(srcObj->type)
	  {
	    case WLZ_2D_POLYGON:
	      dstDom.poly = WlzCMeshTransformPoly(srcObj->domain.poly,
						  mObj, &errNum);
	      break;
	    case WLZ_BOUNDLIST:
	      dstDom.b = WlzCMeshTransformBoundList(srcObj->domain.b,
						    mObj, &errNum);
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj = WlzMakeMain(srcObj->type, dstDom, srcValues,
				 NULL, NULL, &errNum);
	  }
	  if((errNum != WLZ_ERR_NONE) && dstDom.core)
	  {
	    (void )WlzFreeDomain(dstDom);
	  }
	}
	break;
      case WLZ_2D_DOMAINOBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
          dstObj = WlzCMeshToDomObj2D(mObj, 1, &errNum);
	}
	if((errNum == WLZ_ERR_NONE) &&
	   (srcObj->values.core))
	{
	  bgdV = WlzGetBackground(srcObj, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstValues.v = WlzNewValueTb(dstObj, srcObj->values.v->type,
					bgdV, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dstObj->values = WlzAssignValues(dstValues, NULL);
	    errNum = WlzCMeshTransformValues2D(dstObj, srcObj, mObj, interp);
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if(srcObj->values.core == NULL)
	{
	  dstObj = WlzCMeshTransformObjPDomain3D(srcObj, mObj, &errNum);
	}
	else if(WlzGreyTableIsTiled(srcObj->values.core->type))
	{
	  errNum = WLZ_ERR_VALUES_TYPE;
	}
	else
	{
	  dstObj = WlzCMeshTransformObjV3D(srcObj, mObj, interp, &errNum);
	}
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
  return(dstObj);
}

/*!
* \return	Compound array with transformed object domains, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a conforming mesh transform to all the domains of the
* 		given compound array, in a single pass, by creating an grey
* 		valued index object from the domains in the compound array;
* 		transforming the grey valued domain object and then extracting
* 		the domains from the grey valued index object.
* 		Because a conforming mesh transform can be comparatively slow
* 		this can reduce the time taken to apply the same transform to
* 		multiple domains. However it is assumed that none of the
* 		domains intersect. If they do then the domains with higher
* 		indices (later in the compound array) will overwrite those
* 		with lower indices (earlier in the compound array). In many
*		cases this assumption may be valid, eg in the case of
*		exclusive anatomy domains or gene expression strength
*		domains where each is known to be a subset of the preceding
*		strength domains. Even so the results of this function may
*		not be identical to seperate calls at the boundaries of the
*		domains. To distinguish the domains from the background value
*		it may be useful to make the first object of the given compound
*		array an empty object.
* \param	srcObj			Object to be transformed. This must
* 					be a compound array with at least
* 					one non empty object.
* \param	mObj			Conforming mesh transform object.
* \param	interp			Type of interpolation, which should
* 					probably either be nearest neighbours
* 					or classify. Other interpolation
* 					methods may be allowed.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCompoundArray  *WlzCMeshTransformManyObjAsIdx(WlzCompoundArray *srcObj,
				     WlzObject *mObj,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzObject  	*srcGObj = NULL,
  		*dstGObj = NULL;
  WlzCompoundArray *dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
  if((srcObj == NULL) || (mObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(((mObj->type != WLZ_CMESH_2D) &&
           (mObj->type != WLZ_CMESH_2D5) &&
           (mObj->type != WLZ_CMESH_3D)) ||
          ((srcObj->type != WLZ_COMPOUND_ARR_1) && 
           (srcObj->type != WLZ_COMPOUND_ARR_2))||
	  (srcObj->n < 1))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    srcGObj = WlzAssignObject(
              WlzIndexObjFromCompound(srcObj, &errNum), NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstGObj = WlzAssignObject(
              WlzCMeshTransformObj(srcGObj, mObj, interp, &errNum), NULL);
  }
  (void )WlzFreeObj(srcGObj);
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzIndexObjToCompound(dstGObj, &errNum);
  }
  (void )WlzFreeObj(dstGObj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj((WlzObject *)dstObj);
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	Transformed 2D conforming mesh.
* \ingroup	WlzTransform
* \brief	Transforms the given 2D conforming mesh using the given
* 		conforming mesh transform.
* \param	sMesh			Given 2D conforming mesh to be
* 					transformed.
* \param	mTrObj			Mesh transform object to apply.
* \param	newMesh			Make a new mesh if non-zero rather than
* 					transform the source mesh in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh2D *WlzCMeshTransformCMesh2D(WlzCMesh2D *sMesh,
					    WlzObject *mTrObj,
					    int newMesh,
					    WlzErrorNum *dstErr)
{
  WlzCMesh2D	*dMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(newMesh != 0)
  {
    dMesh = WlzCMeshCopy2D(sMesh, 0, 0, NULL, NULL, &errNum);
  }
  else
  {
    dMesh = sMesh;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
		dMaxNod,
	        tLastElmIdx = -1;
    AlcVector	*dNV;
    WlzCMesh2D	*tMesh;
    WlzIndexedValues *tIxv;

    tMesh = mTrObj->domain.cm2;
    tIxv = mTrObj->values.x;
    dNV = dMesh->res.nod.vec;
    dMaxNod = dMesh->res.nod.maxEnt;
    for(idN = 0; idN < dMaxNod; ++idN)
    {
      int	tNearNod;
      double	*dsp;
      WlzDVertex2 dVtx;
      WlzCMeshNod2D *dNod;
      WlzCMeshScanElm2D sE;
      
      tNearNod = -1;
      dNod = (WlzCMeshNod2D *)AlcVectorItemGet(dNV, idN);
      if((dNod != NULL) && (dNod->idx >= 0))
      {
	if(((sE.idx = WlzCMeshElmEnclosingPos2D(tMesh, tLastElmIdx,
			    dNod->pos.vtX, dNod->pos.vtY,
			    0, &tNearNod)) < 0) && (tNearNod < 0))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
	if(sE.idx >= 0)
	{
	  if((sE.idx != tLastElmIdx) ||
	     ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
	  {
	    WlzCMeshUpdateScanElm2D(mTrObj, &sE, 1);
	    tLastElmIdx = sE.idx;
	  }
	  dVtx.vtX = (sE.trX[0] * dNod->pos.vtX) +
		     (sE.trX[1] * dNod->pos.vtY) + sE.trX[2];
	  dVtx.vtY = (sE.trY[0] * dNod->pos.vtX) +
		     (sE.trY[1] * dNod->pos.vtY) + sE.trY[2];
	}
	else
	{
	  dsp = (double *)WlzIndexedValueGet(tIxv, tNearNod);
	  dVtx.vtX = dNod->pos.vtX + dsp[0];
	  dVtx.vtY = dNod->pos.vtY + dsp[1];
	}
	dNod->pos = dVtx;
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dMesh);
}

/*!
* \return	Transformed 2D5 conforming mesh.
* \ingroup	WlzTransform
* \brief	Transforms the given 2D5 conforming mesh using the given
* 		conforming mesh transform.
* \param	sMesh			Given 2D5 conforming mesh to be
* 					transformed.
* \param	mObj			Mesh transform object to apply.
* \param	newMesh			Make a new mesh if non-zero rather than
* 					transform the source mesh in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh2D5 *WlzCMeshTransformCMesh2D5(WlzCMesh2D5 *sMesh,
					      WlzObject *mObj,
					      int newMesh,
					      WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*dMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(newMesh != 0)
  {
    dMesh = WlzCMeshCopy2D5(sMesh, 0, 0, NULL, NULL, &errNum);
  }
  else
  {
    dMesh = sMesh;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
		dMaxNod,
	        tLastElmIdx = -1;
    AlcVector	*dNV;
    WlzCMesh2D5	*tMesh;
    WlzIndexedValues *tIxv;

    tMesh = mObj->domain.cm2d5;
    tIxv = mObj->values.x;
    dNV = dMesh->res.nod.vec;
    dMaxNod = dMesh->res.nod.maxEnt;
    for(idN = 0; idN < dMaxNod; ++idN)
    {
      int	tNearNod;
      double	*dsp;
      WlzDVertex3 dVtx;
      WlzCMeshNod2D5 *dNod;
      WlzCMeshScanElm3D sE;
      
      tNearNod = -1;
      dNod = (WlzCMeshNod2D5 *)AlcVectorItemGet(dNV, idN);
      if((dNod != NULL) && (dNod->idx >= 0))
      {
	if(((sE.idx = WlzCMeshElmEnclosingPos2D5(tMesh, tLastElmIdx,
			    dNod->pos.vtX, dNod->pos.vtY, dNod->pos.vtZ,
			    0, &tNearNod)) < 0) && (tNearNod < 0))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
	if(sE.idx >= 0)
	{
	  if((sE.idx != tLastElmIdx) ||
	     ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
	  {
	    WlzCMeshUpdateScanElm2D5(mObj, &sE, 1);
	    tLastElmIdx = sE.idx;
	  }
	  dVtx.vtX = (sE.tr[ 0] * dNod->pos.vtX) +
		     (sE.tr[ 1] * dNod->pos.vtY) +
		     (sE.tr[ 2] * dNod->pos.vtZ) +  sE.tr[ 3];
	  dVtx.vtY = (sE.tr[ 4] * dNod->pos.vtX) +
		     (sE.tr[ 5] * dNod->pos.vtY) +
		     (sE.tr[ 6] * dNod->pos.vtZ) +  sE.tr[ 7];
	  dVtx.vtZ = (sE.tr[ 8] * dNod->pos.vtX) +
		     (sE.tr[ 9] * dNod->pos.vtY) +
		     (sE.tr[10] * dNod->pos.vtZ) +  sE.tr[11];
	}
	else
	{
	  dsp = (double *)WlzIndexedValueGet(tIxv, tNearNod);
	  dVtx.vtX = dNod->pos.vtX + dsp[0];
	  dVtx.vtY = dNod->pos.vtY + dsp[1];
	  dVtx.vtZ = dNod->pos.vtZ + dsp[2];
	}
	dNod->pos = dVtx;
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dMesh);
}

/*!
* \return	Transformed 3D conforming mesh.
* \ingroup	WlzTransform
* \brief	Transforms the given 3D conforming mesh using the given
* 		conforming mesh transform.
* \param	sMesh			Given 3D conforming mesh to be
* 					transformed.
* \param	mTrObj			Mesh transform object to apply.
* \param	newMesh			Make a new mesh if non-zero rather than
* 					transform the source mesh in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMesh3D *WlzCMeshTransformCMesh3D(WlzCMesh3D *sMesh,
					    WlzObject *mTrObj,
					    int newMesh,
					    WlzErrorNum *dstErr)
{
  WlzCMesh3D	*dMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(newMesh != 0)
  {
    dMesh = WlzCMeshCopy3D(sMesh, 0, 0, NULL, NULL, &errNum);
  }
  else
  {
    dMesh = sMesh;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
		dMaxNod,
	        tLastElmIdx = -1;
    AlcVector	*dNV;
    WlzCMesh3D	*tMesh;
    WlzIndexedValues *tIxv;

    tMesh = mTrObj->domain.cm3;
    tIxv = mTrObj->values.x;
    dNV = dMesh->res.nod.vec;
    dMaxNod = dMesh->res.nod.maxEnt;
    for(idN = 0; idN < dMaxNod; ++idN)
    {
      int	tNearNod;
      double	*dsp;
      WlzDVertex3 dVtx;
      WlzCMeshNod3D *dNod;
      WlzCMeshScanElm3D sE;
      
      tNearNod = -1;
      dNod = (WlzCMeshNod3D *)AlcVectorItemGet(dNV, idN);
      if((dNod != NULL) && (dNod->idx >= 0))
      {
	if(((sE.idx = WlzCMeshElmEnclosingPos3D(tMesh, tLastElmIdx,
			    dNod->pos.vtX, dNod->pos.vtY, dNod->pos.vtZ,
			    0, &tNearNod)) < 0) && (tNearNod < 0))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
	if(sE.idx >= 0)
	{
	  if((sE.idx != tLastElmIdx) ||
	     ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
	  {
	    WlzCMeshUpdateScanElm3D(mTrObj, &sE, 1);
	    tLastElmIdx = sE.idx;
	  }
	  dVtx.vtX = (sE.tr[ 0] * dNod->pos.vtX) +
		     (sE.tr[ 1] * dNod->pos.vtY) +
		     (sE.tr[ 2] * dNod->pos.vtZ) +  sE.tr[ 3];
	  dVtx.vtY = (sE.tr[ 4] * dNod->pos.vtX) +
		     (sE.tr[ 5] * dNod->pos.vtY) +
		     (sE.tr[ 6] * dNod->pos.vtZ) +  sE.tr[ 7];
	  dVtx.vtZ = (sE.tr[ 8] * dNod->pos.vtX) +
		     (sE.tr[ 9] * dNod->pos.vtY) +
		     (sE.tr[10] * dNod->pos.vtZ) +  sE.tr[11];
	}
	else
	{
	  dsp = (double *)WlzIndexedValueGet(tIxv, tNearNod);
	  dVtx.vtX = dNod->pos.vtX + dsp[0];
	  dVtx.vtY = dNod->pos.vtY + dsp[1];
	  dVtx.vtZ = dNod->pos.vtZ + dsp[2];
	}
	dNod->pos = dVtx;
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(dMesh);
}

/*!
* \return	Transformed boundary list or NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms the given boundary list using the given conforming
*		mesh transform.
* \param	srcBound		Given boundary list.
* \param	mesh			Mesh transform object to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzBoundList *WlzCMeshTransformBoundList(WlzBoundList *srcBound,
						WlzObject *mObj,
						WlzErrorNum *dstErr)
{
  WlzDomain	dumDom;
  WlzBoundList	*dstBnd = NULL;
  WlzObject	*plyObj = NULL;
  WlzPolygonDomain *plyDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dstBnd = (WlzBoundList *)AlcCalloc(sizeof(WlzBoundList), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* Wrap set to 1 for closed lines by WlzPolyDecimate(). */
    dstBnd->type = srcBound->type;
    dstBnd->wrap = (srcBound->wrap)? 1: 0;
    /* Transform the polygon. */
    /* Don't decimate the poly becauase it may make poly vertices lie outside
     * the mesh. */
    if((plyObj = WlzPolyTo8Polygon(srcBound->poly,
    				   srcBound->wrap, &errNum)) != NULL)
    {
      if((plyDom = WlzCMeshTransformPoly(plyObj->domain.poly,
				         mObj, &errNum)) != NULL)
      {
	  dstBnd->poly = WlzAssignPolygonDomain(plyDom, NULL);
      }
      else
      {
	dstBnd->poly = NULL;
      }
      (void )WlzFreeObj(plyObj);
    }
  }
  /* Transform next boundlist. */
  if((errNum == WLZ_ERR_NONE) && (srcBound->next != NULL))
  {
    if((dumDom.b = WlzCMeshTransformBoundList(srcBound->next, mObj,
					      &errNum)) != NULL)
    {
      (void )WlzAssignDomain(dumDom, &errNum);
      dstBnd->next = dumDom.b;
    }
  }
  /* Transform down boundlist. */
  if((errNum == WLZ_ERR_NONE) && (srcBound->down != NULL))
  {
    if((dumDom.b = WlzCMeshTransformBoundList(srcBound->down, mObj,
					      &errNum)) != NULL)
    {
      (void )WlzAssignDomain(dumDom, &errNum);
      dstBnd->down = dumDom.b;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    AlcFree(dstBnd); dstBnd = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBnd);
}

/*!
* \return	Transformed polygon domain or NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms the given polygon domain using the given conforming
*		mesh transform.
* \param	srcPoly			Given polygon domain.
* \param	mesh			Mesh transform object to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPolygonDomain *WlzCMeshTransformPoly(WlzPolygonDomain *srcPoly,
					       WlzObject *mObj,
					       WlzErrorNum *dstErr)
{
  WlzPolygonDomain *dstPoly = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((srcPoly->type != WLZ_POLYGON_INT) &&
     (srcPoly->type != WLZ_POLYGON_FLOAT) &&
     (srcPoly->type != WLZ_POLYGON_DOUBLE))
  {
    errNum = WLZ_ERR_POLYGON_TYPE;
  }
  else
  {
    dstPoly = WlzMakePolygonDomain(srcPoly->type,
    				   srcPoly->nvertices, srcPoly->vtx,
				   srcPoly->nvertices, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(srcPoly->type)
    {
      case WLZ_POLYGON_INT:
        errNum = WlzCMeshTransformVtxAry2I(mObj, dstPoly->nvertices,
					   dstPoly->vtx);
        break;
      case WLZ_POLYGON_FLOAT:
        errNum = WlzCMeshTransformVtxAry2F(mObj, dstPoly->nvertices,
					   (WlzFVertex2 *)(dstPoly->vtx));
        break;
      case WLZ_POLYGON_DOUBLE:
        errNum = WlzCMeshTransformVtxAry2D(mObj, dstPoly->nvertices,
					   (WlzDVertex2 *)(dstPoly->vtx));
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePolyDmn(dstPoly);
    dstPoly = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPoly);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a 3D conforming mesh transform to the given source
*		object which must be a 3D domain object. The new object
*		will not have values attached.
* \param	srcObj			Object to be transformed.
* \param	mObj			Conforming mesh transform object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformObjPDomain3D(WlzObject *srcObj,
						WlzObject *mObj,
						WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Make workspace intervals for the elements in the displaced
   * mesh, with intervals sorted by plane, line and then column. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp = WlzCMeshScanWSpInit3D(mObj, 1, &errNum);
  }
  /* Scan through the sorted intervals creating domains as required. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzCMeshScanObjPDomain3D(srcObj, mSWSp, &errNum); 
  }
#ifdef WLZ_CMESHTRANSFORM_DEBUG
  if(errNum  == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshVerifyWSp3D(srcObj, mSWSp);
  }
#endif
  /* Free workspace. */
  WlzCMeshScanWSpFree3D(mSWSp);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a 3D conforming mesh transform to the given contour.
* \param	srcCtr			Object to be transformed.
* \param	mObj			Conforming mesh transform object.
* * \param	newModFlg		Make a new model if non-zero,
* 					otherwise transform the given
* 					model in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzContour	*WlzCMeshTransformContour(WlzContour *srcCtr,
					WlzObject *mObj, int newModFlg,
					WlzErrorNum *dstErr)
{
  WlzContour    *dstCtr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dstCtr = WlzMakeContour(&errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if(srcCtr->model == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      dstCtr->model = WlzAssignGMModel(
	  	      WlzCMeshTransformGMModel(srcCtr->model, mObj, newModFlg,
	                                       &errNum), NULL);
    }
  }
  if((errNum != WLZ_ERR_NONE) && dstCtr)
  {
    (void )WlzFreeContour(dstCtr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstCtr);
}

/*!
* \return	Transformed model or NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a conforming mesh transform to the given 2
* 		dimensional geometric model.
* \param	srcM			Given geometric model.
* \param	trObj			Given conforming mesh transform.
* \param	newModFlg		Make a new model if non-zero,
* 					otherwise transform the given
* 					model in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzGMModel	*WlzCMeshTransformGMModel(WlzGMModel *srcM,
					WlzObject *trObj, int newModFlg,
					WlzErrorNum *dstErr)
{
  WlzGMModel	*dstM = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcM == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    dstM = (newModFlg)? WlzGMModelCopy(srcM, &errNum):
	   WlzAssignGMModel(srcM, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idx,
    		cnt;
    AlcVector	*vec;

    vec = dstM->res.vertexG.vec;
    cnt = (int )(dstM->res.vertexG.numIdx);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idx = 0; idx < cnt; ++idx)
    {
      if(errNum == WLZ_ERR_NONE)
      {
        WlzGMElemP      elmP;

        elmP.core = (WlzGMCore *)AlcVectorItemGet(vec, (size_t )idx);
        if(elmP.core && (elmP.core->idx >= 0))
        {
          WlzErrorNum   errNum2 = WLZ_ERR_NONE;

          switch(dstM->type)
          {
            case WLZ_GMMOD_2I:
	      errNum2 = WlzCMeshTransformVtxAry2I(trObj, 1,
	                                          &(elmP.vertexG2I->vtx));
              break;
            case WLZ_GMMOD_2D:
	      errNum2 = WlzCMeshTransformVtxAry2D(trObj, 1,
	                                          &(elmP.vertexG2D->vtx));
              break;
            case WLZ_GMMOD_3I:
	      errNum2 = WlzCMeshTransformVtxAry3I(trObj, 1,
	                                          &(elmP.vertexG3I->vtx));
              break;
            case WLZ_GMMOD_3D:
	      errNum2 = WlzCMeshTransformVtxAry3D(trObj, 1,
	                                          &(elmP.vertexG3D->vtx));
              break;
            default:
              errNum2 = WLZ_ERR_DOMAIN_TYPE;
              break;
          }
#ifdef _OPENMP
#pragma omp critical
          {
#endif
            if((errNum == WLZ_ERR_NONE) && (errNum2 != WLZ_ERR_NONE))
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
  if(errNum == WLZ_ERR_NONE)
  {
    /* Rehash the geometric model. */
    errNum = WlzGMModelRehashVHT(dstM, 0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstM);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a 3D conforming mesh transform to the given source
*		object which must be a 3D domain object with values.
* \param	srcObj			Object to be transformed.
* \param	mObj			Conforming mesh transform object.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformObjV3D(WlzObject *srcObj,
				     WlzObject *mObj,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzPixelV	bgdV;
  WlzGreyType	gType;
  WlzObjectType gTType;
  WlzObject	*dstObj = NULL;
  WlzValues	dstValues;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gType = WlzGreyTypeFromObj(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    bgdV = WlzGetBackground(srcObj, &errNum);
  }
  /* Make workspace intervals for the elements in the displaced
   * mesh, with intervals sorted by plane, line and then column. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp = WlzCMeshScanWSpInit3D(mObj, 1, &errNum);
  }
  /* Scan through the sorted intervals creating domains as required. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzCMeshScanObjPDomain3D(srcObj, mSWSp, &errNum); 
  }
  /* Make a voxel value table for the new domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    gTType = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, gType, NULL);
    dstValues.vox = WlzNewValuesVox(dstObj, gTType, bgdV, &errNum);
  }
  /* Scan through the sorted intervals again setting object values. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj->values = WlzAssignValues(dstValues, NULL);
    errNum = WlzCMeshScanObjValues3D(dstObj, srcObj, mSWSp, interp);
  }
#ifdef WLZ_CMESHTRANSFORM_DEBUG
  if(errNum  == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshVerifyWSp3D(srcObj, mSWSp);
  }
#endif
  /* Free workspace. */
  WlzCMeshScanWSpFree3D(mSWSp);
  /* Clean up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(dstObj);
    dstObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

#ifdef WLZ_CMESHTRANSFORM_DEBUG
/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Verifies that all voxels of the source object are
*		within mesh scan intervals. Outputs a message if
*		a voxel is not in any scan intervals.
*		This function is only intended for use in debuging
*		the C mesh transform code.
* \param	srcObj			Source object with valid planedomain.
* \param	mSWSp			Conforming mesh scan workspace.
*/
static WlzErrorNum WlzCMeshVerifyWSp3D(WlzObject *srcObj,
				       WlzCMeshScanWSp3D *mSWSp)
{
  int		idI,
		idP,
		idN,
		cnt;
  int		*inp;
  WlzPixelV     bgdV;
  WlzValues     tstValues;
  WlzIVertex3	dPos,
  		sPos;
  WlzDVertex3	tV;
  WlzObject	*obj2 = NULL,
  		*tstObj = NULL;
  WlzDomain	dom2;
  WlzCMeshScanItv3D *curItv;
  WlzCMeshScanElm3D *sE;
  WlzIntervalWSpace iWSp;
  WlzGreyWSpace	gWSp;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  bgdV.v.inv = 0;
  bgdV.type = WLZ_GREY_INT;
  /* Make new test object with same domain as the source object and with
   * all values initialized to zero. */
  tstValues.vox = WlzNewValuesVox(srcObj, WLZ_GREY_INT, bgdV, &errNum);
  if(errNum  == WLZ_ERR_NONE)
  {
    tstObj = WlzMakeMain(srcObj->type, srcObj->domain, tstValues,
			 NULL, NULL, &errNum);
  }
  if(errNum  == WLZ_ERR_NONE)
  {
    errNum = WlzGreySetValue(tstObj, bgdV);
  }
  /* Set all test object grey values covered by mesh scan intervals to non
   * zero. */
  if(errNum  == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(tstObj, &errNum);
  }
  if(errNum  == WLZ_ERR_NONE)
  {
    idI = 0;
    curItv = mSWSp->itvs;
  }
  while((errNum == WLZ_ERR_NONE) && (idI < mSWSp->nItvs))
  {
    sE = mSWSp->dElm + curItv->elmIdx;
    dPos.vtY = curItv->line;
    dPos.vtZ = curItv->plane;
    for(dPos.vtX = curItv->lftI; dPos.vtX <= curItv->rgtI; ++dPos.vtX)
    {
      if((sE->flags & WLZ_CMESH_SCANELM_REV) == 0)
      {
	WlzCMeshUpdateScanElm3D(mSWSp->mTr, sE, 0);
      }
      tV.vtX = (sE->tr[ 0] * dPos.vtX) + (sE->tr[ 1] * dPos.vtY) +
	       (sE->tr[ 2] * dPos.vtZ) +  sE->tr[ 3];
      tV.vtY = (sE->tr[ 4] * dPos.vtX) + (sE->tr[ 5] * dPos.vtY) +
	       (sE->tr[ 6] * dPos.vtZ) +  sE->tr[ 7];
      tV.vtZ = (sE->tr[ 8] * dPos.vtX) + (sE->tr[ 9] * dPos.vtY) +
	       (sE->tr[10] * dPos.vtZ) +  sE->tr[11];
      sPos.vtX = WLZ_CMESH_POS_DTOI(tV.vtX);
      sPos.vtY = WLZ_CMESH_POS_DTOI(tV.vtY);
      sPos.vtZ = WLZ_CMESH_POS_DTOI(tV.vtZ);
      if(WlzInsideDomain(srcObj, sPos.vtZ, sPos.vtY,
			 sPos.vtX, NULL) != 0)
      {
	WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
	*(gVWSp->gPtr[0].inp) = sE->idx + 1;
      }
    }
    ++idI;
    ++curItv;
  }
  /* Check the test object for values still zero. Such voxels are not covered
   * by the scan intervals - report them. */
  if(errNum == WLZ_ERR_NONE)
  {
    idP = 0;
    sPos.vtZ = tstObj->domain.p->plane1;
    while((errNum == WLZ_ERR_NONE) && (sPos.vtZ < tstObj->domain.p->lastpl))
    {
      if(((dom2 = *(tstObj->domain.p->domains + idP)).core != NULL) &&
         (dom2.core->type != WLZ_EMPTY_DOMAIN))
      {
	(void )WlzFreeObj(obj2);
        obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                      *(tstObj->domain.p->domains + idP),
			      *(tstObj->values.vox->values + idP),
			      NULL, NULL, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(obj2, &iWSp, &gWSp);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  while((errNum == WLZ_ERR_NONE) &&
		((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
	  {
	    cnt = iWSp.rgtpos - iWSp.lftpos + 1;
	    inp = gWSp.u_grintptr.inp;
	    sPos.vtY = iWSp.linpos;
	    for(idI = 0; idI < cnt; ++idI)
	    {
	      idN = 0;
	      if(*inp == 0)
	      {
		/* Found voxel not covered by the scan intervals. */
		sPos.vtX = iWSp.lftpos + idI;
		/* Find element C mesh which encloses the voxel. */
		idN = WlzCMeshElmEnclosingPos3D(mSWSp->mTr->mesh.m3, -1,
						sPos.vtX, sPos.vtY, sPos.vtZ,
						0, NULL);
		/* Output message. */
		(void )fprintf(stderr, "WlzCMeshVerifyWSp3D() %d %d %d n %d\n",
			       sPos.vtX, sPos.vtY, sPos.vtZ, idN - 1);
	      }
	      ++inp;
	    }
	  }
	  (void )WlzEndGreyScan(&iWSp, &gWSp);
	}
	(void )WlzFreeObj(obj2);
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      ++idP;
      ++sPos.vtZ;
    }
  }

  WlzGreyValueFreeWSp(gVWSp);
  (void )WlzFreeObj(tstObj);
  return(errNum);
}
#endif

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a 3D conforming mesh transform to the given source
*		object (which must be a 3D domain object) using the already
*		initialized mesh transform workspace.
*		If the given source object is NULL a domain will be created
*		corresponding to the given mesh transform.
* \param	srcObj			Object to be transformed, may be NULL.
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshScanObjPDomain3D(WlzObject *srcObj,
					   WlzCMeshScanWSp3D *mSWSp,
					   WlzErrorNum *dstErr) 
{
  int		idI,
		kol,
		itvLnCnt = 0,
		itvPlCnt = 0,
		itvLnWidth,
		itvLnByteWidth;
  WlzIVertex3	dPos,
		sPos;
  WlzDVertex3	tV;
  WlzDynItvPool	itvPool;
  WlzCMeshScanItv3D *curItv = NULL,
  		*prvItv = NULL;
  WlzObjectType	dstObjType;
  WlzCMeshScanElm3D *sE;
  WlzUByte	*lnMsk = NULL;
  WlzDomain	dom2,
  		dom3;
  WlzValues	nullVal;
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minDynItv = 1024; /* This is the number of intervals that are
  				     allocated in a single block. It is a
				     tuning parameter, see use below. */

  dom2.core = NULL;
  dom3.core = NULL;
  nullVal.core = NULL;
  itvPool.offset = 0;
  itvPool.itvBlock = NULL;
  if(mSWSp->nItvs < 1)
  {
    dstObj = WlzMakeEmpty(&errNum);
  }
  else
  {
    dstObjType = (srcObj == NULL)? WLZ_3D_DOMAINOBJ: srcObj->type;
    itvLnWidth = mSWSp->dBox.xMax - mSWSp->dBox.xMin + 1;
    itvLnByteWidth = (itvLnWidth + 7) / 8;
    /* Initialize a dynamic interval pool. Any size greater than the maximum
     * number of intervals per line will do, but for efficiency it shouldn't
     * be too small as this will cause loads of memory allocations. */
    itvPool.itvsInBlock = (itvLnWidth < minDynItv)? minDynItv: itvLnWidth;
    /* Create a new plane domain using the bounding box of the displaced
     * mesh. This is corrected latter. */
    dom3.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				mSWSp->dBox.zMin, mSWSp->dBox.zMax,
				mSWSp->dBox.yMin, mSWSp->dBox.yMax,
				mSWSp->dBox.xMin, mSWSp->dBox.xMax,
				&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if(AlcBit1Calloc(&lnMsk, itvLnWidth) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      idI = 0;
      itvLnCnt = 0;
      itvPlCnt = 0;
      itvPool.offset = 0;
      itvPool.itvBlock = NULL;
      prvItv = NULL;
      curItv = mSWSp->itvs;
    }
    while((errNum == WLZ_ERR_NONE) && (idI < mSWSp->nItvs))
    {
      sE = mSWSp->dElm + curItv->elmIdx;
      if(dom2.core == NULL)
      {
	dom2.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				       mSWSp->dBox.yMin, mSWSp->dBox.yMax,
				       mSWSp->dBox.xMin, mSWSp->dBox.xMax,
				       &errNum);
      }
      dPos.vtY = curItv->line;
      dPos.vtZ = curItv->plane;
      for(kol = curItv->lftI; kol <= curItv->rgtI; ++kol)
      {
        dPos.vtX = kol;
	if((sE->flags & WLZ_CMESH_SCANELM_REV) == 0)
	{
	  WlzCMeshUpdateScanElm3D(mSWSp->mTr, sE, 0);
	}
	tV.vtX = (sE->tr[ 0] * dPos.vtX) + (sE->tr[ 1] * dPos.vtY) +
		 (sE->tr[ 2] * dPos.vtZ) +  sE->tr[ 3];
	tV.vtY = (sE->tr[ 4] * dPos.vtX) + (sE->tr[ 5] * dPos.vtY) +
		 (sE->tr[ 6] * dPos.vtZ) +  sE->tr[ 7];
	tV.vtZ = (sE->tr[ 8] * dPos.vtX) + (sE->tr[ 9] * dPos.vtY) +
		 (sE->tr[10] * dPos.vtZ) +  sE->tr[11];
	sPos.vtX = WLZ_CMESH_POS_DTOI(tV.vtX);
	sPos.vtY = WLZ_CMESH_POS_DTOI(tV.vtY);
	sPos.vtZ = WLZ_CMESH_POS_DTOI(tV.vtZ);
	if((srcObj == NULL) ||
	   (WlzInsideDomain(srcObj, sPos.vtZ, sPos.vtY, sPos.vtX, NULL) != 0))
	{
	  ++itvLnCnt;
	  WlzBitLnSetItv(lnMsk,
	                 kol - mSWSp->dBox.xMin, kol - mSWSp->dBox.xMin,
			 itvLnWidth);
	}
      }
      prvItv = curItv;
      ++idI;
      ++curItv;
      if((errNum == WLZ_ERR_NONE) &&
         (itvLnCnt > 0) &&
	 ((curItv->plane != prvItv->plane) || (curItv->line != prvItv->line)))
      {
	/* Add previous line to interval domain. */
	errNum = WlzDynItvLnFromBitLn(dom2.i, lnMsk, prvItv->line, itvLnWidth,
	                              &itvPool);
	memset(lnMsk, 0, itvLnByteWidth);
	itvPlCnt += itvLnCnt;
	itvLnCnt = 0;
      }
      if((errNum == WLZ_ERR_NONE) &&
         (itvPlCnt > 0) &&
	 (curItv->plane != prvItv->plane))
      {
        /* Add previous plane to plane domain. */
	*(dom3.p->domains + prvItv->plane - dom3.p->plane1) =
				 WlzAssignDomain(dom2, NULL);
        itvPlCnt = 0;
	itvPool.offset = 0;
	itvPool.itvBlock = NULL;
        dom2.core = NULL;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Standardize the domains to account for the correct bounding boxes. */
      errNum = WlzStandardPlaneDomain(dom3.p, NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Set voxel size. */
      if(srcObj != NULL)
      {
	dom3.p->voxel_size[0] = srcObj->domain.p->voxel_size[0];
	dom3.p->voxel_size[1] = srcObj->domain.p->voxel_size[1];
	dom3.p->voxel_size[2] = srcObj->domain.p->voxel_size[2];
      }
      /* Create new object from the transformed plane domain. */
      dstObj = WlzMakeMain(dstObjType, dom3, nullVal, NULL, NULL, &errNum);
    }
  }
  /* Clear up. */
  AlcFree(lnMsk);
  (void )WlzFreeDomain(dom2);
  /* Clear up on error. */
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePlaneDomain(dom3.p);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Fills in the destination object's values from the source
*		object, using the mesh scan workspace.
* \param	dstObj			Destination object with values to be
*					set.
* \param	srcObj			Source object.
* \param	mSWSp			Mesh scan workspace which was used to
*					compute the destination object's
*					domain.
* \param	interp			Interpolation type.
*/
static WlzErrorNum WlzCMeshScanObjValues3D(WlzObject *dstObj,
					WlzObject *srcObj,
					WlzCMeshScanWSp3D *mSWSp,
					WlzInterpolationType interp)
{
  int		idP,
  		idI,
  		iLft,
		iRgt,
		mItvIdx0,
  		mItvIdx1,
		bufWidth,
  		itvWidth;
  double	tD0,
  		tD1,
		tD2,
		tD3,
		tD4;
  int		*olpCnt = NULL;
  WlzGreyP	dGP,
  		olpBuf;
  WlzGreyType	gType;
  WlzPixelV	bgdV;
  WlzIVertex3	dPos,
  		sPos;
  WlzDVertex3	tV,
  		sPosD;
  WlzCMeshScanElm3D *sE;
  WlzCMeshScanItv3D *mItv0,
  		*mItv1,
		*mItv2;
  WlzDomain	dom2,
  		dom3;
  WlzValues	val3;
  WlzObject	*obj2 = NULL;
  WlzGreyWSpace gWSp;
  WlzIntervalWSpace iWSp;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  olpBuf.inp = NULL;
  bgdV = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    gType = WlzGreyTypeFromObj(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&bgdV, bgdV, gType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bufWidth = dstObj->domain.p->lastkl - dstObj->domain.p->kol1 + 1;
    errNum = WlzCMeshScanMakeOlpBufs(dstObj, gType,
                                     &olpBuf, &olpCnt, bufWidth);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItvIdx0 = 0;
    mItv0 = mSWSp->itvs;
    gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idP = 0;
    dom3.p = dstObj->domain.p;
    val3.vox = dstObj->values.vox;
    dPos.vtZ = dom3.p->plane1;
    while((errNum == WLZ_ERR_NONE) && (dPos.vtZ < dom3.p->lastpl))
    {
      if(((dom2 = *(dom3.p->domains + idP)).core != NULL) &&
         (dom2.core->type != WLZ_EMPTY_DOMAIN))
      {
	(void )WlzFreeObj(obj2);
        obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                      *(dom3.p->domains + idP),
			      *(val3.vox->values + idP),
			      NULL, NULL, &errNum);
        if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzInitGreyScan(obj2, &iWSp, &gWSp);
	}
	while((errNum == WLZ_ERR_NONE) &&
	      ((errNum = WlzNextGreyInterval(&iWSp)) == WLZ_ERR_NONE))
        {
          itvWidth = iWSp.rgtpos - iWSp.lftpos + 1;
	  WlzCMeshScanClearOlpBuf(olpBuf, olpCnt, gType, bufWidth, itvWidth);
	  dGP = gWSp.u_grintptr;
	  /* Update the mesh interval pointer so that it points to the
	   * first mesh interval on the which intersects the current grey
	   * interval. */
	  while((mItv0->plane < dPos.vtZ) && (mItvIdx0 < mSWSp->nItvs))
	  {
	    ++mItvIdx0;
	    ++mItv0;
	  }
	  while((mItv0->line < iWSp.linpos) && (mItvIdx0 < mSWSp->nItvs))
	  {
	    ++mItvIdx0;
	    ++mItv0;
	  }
	  while((mItv0->line <= iWSp.linpos) &&
		(mItv0->rgtI < iWSp.lftpos) &&
		(mItvIdx0 < mSWSp->nItvs))
	  {
	    ++mItvIdx0;
	    ++mItv0;
	  }
	  if((mItv0->line == iWSp.linpos) &&
	     (iWSp.lftpos <= mItv0->rgtI) &&
	     (iWSp.rgtpos >= mItv0->lftI))
	  {
	    /* Mesh interval mItv0 intersects the current grey interval find
	     * the last mesh interval mItv1 which also intersects the current
	     * grey interval. */
	    mItv1 = mItv0;
	    mItvIdx1 = mItvIdx0;
	    while((mItv1->line == iWSp.linpos) &&
		  (mItv1->lftI <= iWSp.rgtpos) &&
		  (mItvIdx1 < mSWSp->nItvs))
	    {
	      ++mItvIdx1;
	      ++mItv1;
	    }
	    mItv2 = mItv1 - 1;
	    mItv1 = mItv0;
	    dPos.vtY = mItv0->line;
	    /* For each mesh interval which intersects the current grey
	       interval. */
	    while(mItv1 <= mItv2)
	    {
#ifdef WLZ_CMESHTRANSFORM_DEBUG
      (void )fprintf(stderr,
      		     "WlzCMeshScanObjValues3D %d %d %d %d %d\n",
		     mItv1->elmIdx,
		     mItv1->lftI, mItv1->rgtI, mItv1->line, mItv1->plane);
#endif
	      /* Update mesh scanning. */
	      sE = mSWSp->dElm + mItv1->elmIdx;
	      if((sE->flags & WLZ_CMESH_SCANELM_REV) == 0)
	      {
	         WlzCMeshUpdateScanElm3D(mSWSp->mTr, sE, 0);
	      }
	      tV.vtX = (sE->tr[ 1] * dPos.vtY) + (sE->tr[ 2] * dPos.vtZ) +
		       sE->tr[ 3];
	      tV.vtY = (sE->tr[ 5] * dPos.vtY) + (sE->tr[ 6] * dPos.vtZ) +
		       sE->tr[ 7];
	      tV.vtZ = (sE->tr[ 9] * dPos.vtY) + (sE->tr[10] * dPos.vtZ) +
		       sE->tr[11];
	      /* Find length of intersection and set the grey pointer. */
	      iLft = ALG_MAX(mItv1->lftI, iWSp.lftpos);
	      iRgt = ALG_MIN(mItv1->rgtI, iWSp.rgtpos);
	      dPos.vtX = iLft;
	      switch(interp)
	      {
		case WLZ_INTERPOLATION_NEAREST:
		  switch(gType)
		  {
		    case WLZ_GREY_INT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += gVWSp->gVal[0].inv;
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_SHORT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += gVWSp->gVal[0].shv;
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_UBYTE:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += gVWSp->gVal[0].ubv;
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_FLOAT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.dbp + idI) += gVWSp->gVal[0].flv;
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_DOUBLE:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.dbp + idI) += gVWSp->gVal[0].dbv;
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_RGBA:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += WLZ_RGBA_RED_GET(
						 gVWSp->gVal[0].rgbv);
			  *(olpBuf.inp + bufWidth + idI) +=
			      WLZ_RGBA_GREEN_GET(gVWSp->gVal[0].rgbv);
			  *(olpBuf.inp + (2 * bufWidth) + idI) +=
			      WLZ_RGBA_BLUE_GET(gVWSp->gVal[0].rgbv);
			  *(olpBuf.inp + (3 * bufWidth) + idI) +=
			      WLZ_RGBA_ALPHA_GET(gVWSp->gVal[0].rgbv);
			}
			++dPos.vtX;
		      }
		      break;
		    default:
		      errNum = WLZ_ERR_GREY_TYPE;
		      break;
		  }
		  break;
		case WLZ_INTERPOLATION_LINEAR:
		  switch(gType)
		  {
		    case WLZ_GREY_INT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD0 = ((gVWSp->gVal[0]).inv * tD2 * tD3) +
				((gVWSp->gVal[1]).inv * tD0 * tD3) +
				((gVWSp->gVal[2]).inv * tD2 * tD1) +
				((gVWSp->gVal[3]).inv * tD0 * tD1);
			  tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += WLZ_NINT(tD0);
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY,
					  sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.inp + idI) += gVWSp->gVal[0].inv;
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_SHORT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD0 = ((gVWSp->gVal[0]).shv * tD2 * tD3) +
				((gVWSp->gVal[1]).shv * tD0 * tD3) +
				((gVWSp->gVal[2]).shv * tD2 * tD1) +
				((gVWSp->gVal[3]).shv * tD0 * tD1);
			  tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += WLZ_NINT(tD0);
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY, sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.inp + idI) += gVWSp->gVal[0].shv;
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_UBYTE:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD0 = ((gVWSp->gVal[0]).ubv * tD2 * tD3) +
				((gVWSp->gVal[1]).ubv * tD0 * tD3) +
				((gVWSp->gVal[2]).ubv * tD2 * tD1) +
				((gVWSp->gVal[3]).ubv * tD0 * tD1);
			  tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += WLZ_NINT(tD0);
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY,
					  sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.inp + idI) += gVWSp->gVal[0].ubv;
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_FLOAT:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD0 = ((gVWSp->gVal[0]).flv * tD2 * tD3) +
				((gVWSp->gVal[1]).flv * tD0 * tD3) +
				((gVWSp->gVal[2]).flv * tD2 * tD1) +
				((gVWSp->gVal[3]).flv * tD0 * tD1);
			  tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.dbp + idI) += tD0;
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY,
					  sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.dbp + idI) += gVWSp->gVal[0].flv;
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_DOUBLE:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD0 = ((gVWSp->gVal[0]).dbv * tD2 * tD3) +
				((gVWSp->gVal[1]).dbv * tD0 * tD3) +
				((gVWSp->gVal[2]).dbv * tD2 * tD1) +
				((gVWSp->gVal[3]).dbv * tD0 * tD1);
			  tD0 = WLZ_CLAMP(tD0, 0.0, 255.0);
			  idI = dPos.vtX - iWSp.lftpos;
			  ++*(olpCnt + idI);
			  *(olpBuf.dbp + idI) += tD0;
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY,
					  sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.dbp + idI) += gVWSp->gVal[0].dbv;
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    case WLZ_GREY_RGBA:
		      while(dPos.vtX <= iRgt)
		      {
			idI = dPos.vtX - iWSp.lftpos;
			sPosD.vtX = (sE->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sE->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sE->tr[ 8] * dPos.vtX) + tV.vtZ;
			WlzGreyValueGetCon(gVWSp, sPosD.vtZ, sPosD.vtY,
					   sPosD.vtX);
			if(gVWSp->bkdFlag == 0)
			{
			  tD0 = sPosD.vtX - floor(sPosD.vtX);
			  tD1 = sPosD.vtY - floor(sPosD.vtY);
			  tD2 = 1.0 - tD0;
			  tD3 = 1.0 - tD1;
			  tD4 = (WLZ_RGBA_RED_GET((gVWSp->gVal[0]).rgbv) *
				 tD2 * tD3) +
				(WLZ_RGBA_RED_GET((gVWSp->gVal[1]).rgbv) *
				 tD0 * tD3) +
				(WLZ_RGBA_RED_GET((gVWSp->gVal[2]).rgbv) *
				 tD2 * tD1) +
				(WLZ_RGBA_RED_GET((gVWSp->gVal[3]).rgbv) *
				 tD0 * tD1);
			  tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
			  ++*(olpCnt + idI);
			  *(olpBuf.inp + idI) += WLZ_NINT(tD4);
			  tD4 = (WLZ_RGBA_GREEN_GET((gVWSp->gVal[0]).rgbv) *
				 tD2 * tD3) +
				(WLZ_RGBA_GREEN_GET((gVWSp->gVal[1]).rgbv) *
				 tD0 * tD3) +
				(WLZ_RGBA_GREEN_GET((gVWSp->gVal[2]).rgbv) *
				 tD2 * tD1) +
				(WLZ_RGBA_GREEN_GET((gVWSp->gVal[3]).rgbv) *
				 tD0 * tD1);
			  tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
			  *(olpBuf.inp + bufWidth + idI) += WLZ_NINT(tD4);
			  tD4 = (WLZ_RGBA_BLUE_GET((gVWSp->gVal[0]).rgbv) *
				 tD2 * tD3) +
				(WLZ_RGBA_BLUE_GET((gVWSp->gVal[1]).rgbv) *
				 tD0 * tD3) +
				(WLZ_RGBA_BLUE_GET((gVWSp->gVal[2]).rgbv) *
				 tD2 * tD1) +
				(WLZ_RGBA_BLUE_GET((gVWSp->gVal[3]).rgbv) *
				 tD0 * tD1);
			  tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
			  *(olpBuf.inp + (2 * bufWidth) + idI) +=
			      WLZ_NINT(tD4);
			  tD4 = (WLZ_RGBA_ALPHA_GET((gVWSp->gVal[0]).rgbv) *
				 tD2 * tD3) +
				(WLZ_RGBA_ALPHA_GET((gVWSp->gVal[1]).rgbv) *
				 tD0 * tD3) +
				(WLZ_RGBA_ALPHA_GET((gVWSp->gVal[2]).rgbv) *
				 tD2 * tD1) +
				(WLZ_RGBA_ALPHA_GET((gVWSp->gVal[3]).rgbv) *
				 tD0 * tD1);
			  tD4 = WLZ_CLAMP(tD4, 0.0, 255.0);
			  *(olpBuf.inp + (3 * bufWidth) + idI) +=
			      WLZ_NINT(tD4);
			}
			else
			{
			  sPos.vtX = WLZ_CMESH_POS_DTOI(sPosD.vtX);
			  sPos.vtY = WLZ_CMESH_POS_DTOI(sPosD.vtY);
			  sPos.vtZ = WLZ_CMESH_POS_DTOI(sPosD.vtZ);
			  WlzGreyValueGet(gVWSp, sPos.vtZ, sPos.vtY,
					  sPos.vtX);
			  if(gVWSp->bkdFlag == 0)
			  {
			    idI = dPos.vtX - iWSp.lftpos;
			    ++*(olpCnt + idI);
			    *(olpBuf.inp + idI) +=
				      WLZ_RGBA_RED_GET(gVWSp->gVal[0].rgbv);
			    *(olpBuf.inp + bufWidth + idI) +=
				      WLZ_RGBA_GREEN_GET(gVWSp->gVal[0].rgbv);
			    *(olpBuf.inp + (2 * bufWidth) + idI) +=
				      WLZ_RGBA_BLUE_GET(gVWSp->gVal[0].rgbv);
			    *(olpBuf.inp + (3 * bufWidth) + idI) +=
				      WLZ_RGBA_ALPHA_GET(gVWSp->gVal[0].rgbv);
			  }
			}
			++dPos.vtX;
		      }
		      break;
		    default:
		      errNum = WLZ_ERR_GREY_TYPE;
		      break;
		  }
		  break;
		case WLZ_INTERPOLATION_CLASSIFY_1:     /* FALLTHROUGH */
		  errNum = WLZ_ERR_UNIMPLEMENTED;
		  break;
		default:
		  errNum = WLZ_ERR_INTERPOLATION_TYPE;
		  break;
	      }
	      ++mItv1;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzCMeshScanFlushOlpBuf(dGP, olpBuf, olpCnt, bufWidth,
					     bgdV, iWSp.lftpos, iWSp.rgtpos,
					     interp, gType);
	  }
        }
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
      }
      ++idP;
      ++dPos.vtZ;
    }
  }
  (void )WlzFreeObj(obj2);
  AlcFree(olpBuf.inp);
  AlcFree(olpCnt);
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzTransform
* \brief	Clears the overlap buffers for a given region.
* \param	olpBuf			Overlap grey data buffer.
* \param	olpCnt			Overlap counter.
* \param	gType			Grey type.
* \param	bufWidth		Width of the buffers.
* \param	clrWidth		Width to be cleared.
*/
static void	WlzCMeshScanClearOlpBuf(WlzGreyP olpBuf, int *olpCnt,
  					WlzGreyType gType,
					int bufWidth, int clrWidth)
{
  switch(gType)
  {
    case WLZ_GREY_INT:   /* FALLTHROUGH */
    case WLZ_GREY_SHORT: /* FALLTHROUGH */
    case WLZ_GREY_UBYTE:
      WlzValueSetInt(olpBuf.inp, 0, clrWidth);
      break;
    case WLZ_GREY_FLOAT: /* FALLTHROUGH */
    case WLZ_GREY_DOUBLE:
      WlzValueSetDouble(olpBuf.dbp, 0.0, clrWidth);
      break;
    case WLZ_GREY_RGBA:
      WlzValueSetInt(olpBuf.inp, 0, clrWidth);
      WlzValueSetInt(olpBuf.inp + bufWidth, 0, clrWidth);
      WlzValueSetInt(olpBuf.inp + (2 * bufWidth), 0, clrWidth);
      WlzValueSetInt(olpBuf.inp + (3 * bufWidth), 0, clrWidth);
      break;
    default:
      break;
  }
  WlzValueSetInt(olpCnt, 0, clrWidth);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Allocates overlap buffers for the given target object.
* \param	obj			Target object.
* \param	gType			Grey type.
* \param	dstOlpBuf		Destination pointer for the
*					overlap grey pointer.
* \param	dstOlpCnt		Destination pointer for the
*					overlap count.
* \param	bufWidth		Overlap buffer width.
*/
static WlzErrorNum WlzCMeshScanMakeOlpBufs(WlzObject *obj, WlzGreyType gType,
					WlzGreyP *dstOlpBuf,
					int **dstOlpCnt,
					int bufWidth)
{
  int		*olpCnt = NULL;
  WlzGreyP	olpBuf;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  switch(gType)
  {
    case WLZ_GREY_INT:   /* FALLTHROUGH */
    case WLZ_GREY_SHORT: /* FALLTHROUGH */
    case WLZ_GREY_UBYTE:
      if((olpBuf.inp = (int *)AlcCalloc(bufWidth, sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
    case WLZ_GREY_FLOAT: /* FALLTHROUGH */
    case WLZ_GREY_DOUBLE:
      if((olpBuf.dbp = (double *)AlcCalloc(bufWidth,
					   sizeof(double))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
    case WLZ_GREY_RGBA:
      if((olpBuf.inp = (int *)AlcCalloc(bufWidth * 4, sizeof(int))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
    default:
      errNum = WLZ_ERR_GREY_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((olpCnt = (int *)AlcCalloc(bufWidth, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstOlpBuf = olpBuf;
    *dstOlpCnt = olpCnt;
  }
  else
  {
    AlcFree(olpBuf.v);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Flushes the overlap buffers into the given grey interval.
* \param	dGP			Grey interval pointer.
* \param	olpBuf			Overlap buffer with grey data.
* \param	olpCnt			Overlap buffer with count values.
* \param	bufWidth		Size of the buffers.
* \param	bgdV			Background value.
* \param	iLft			Interval left column.
* \param	iRgt			Interval right column.
* \param	interp			Interpolation.
* \param	gType			Grey type.
*/
static WlzErrorNum WlzCMeshScanFlushOlpBuf(WlzGreyP dGP,
				WlzGreyP olpBuf, int *olpCnt, int bufWidth,
				WlzPixelV bgdV, int iLft,  int iRgt,
                                WlzInterpolationType interp,
				WlzGreyType gType)
{
  int		idX,
  		cnt,
		gVR,
		gVG,
		gVB,
		gVA;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cnt = iRgt - iLft + 1;
  switch(interp)
  {
    case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
    case WLZ_INTERPOLATION_LINEAR:
      switch(gType)
      {
	case WLZ_GREY_INT:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.inp + idX) = (*(olpCnt + idX) >= 1)?
			       *(olpBuf.inp + idX) / *(olpCnt + idX):
			       bgdV.v.inv;
	  }
	  break;
	case WLZ_GREY_SHORT:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.shp + idX) = (*(olpCnt + idX) >= 1)?
			       (short )(*(olpBuf.inp + idX) / *(olpCnt + idX)):
			       bgdV.v.shv;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.ubp + idX) = (*(olpCnt + idX) >= 1)?
			       (WlzUByte )(*(olpBuf.inp + idX) /
			                   *(olpCnt + idX)):
			       bgdV.v.ubv;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.flp + idX) = (*(olpCnt + idX) >= 1)?
			       (float )(*(olpBuf.dbp + idX) / *(olpCnt + idX)):
			       bgdV.v.flv;
	  }
	  break;
	case WLZ_GREY_DOUBLE:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.dbp + idX) = (*(olpCnt + idX) >= 1)?
			       *(olpBuf.dbp + idX) / *(olpCnt + idX):
			       bgdV.v.dbv;
	  }
	  break;
	case WLZ_GREY_RGBA:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    if(*(olpCnt + idX) >= 1)
	    {
	      gVR = *(olpBuf.inp + idX) / *(olpCnt + idX);
	      gVG = *(olpBuf.inp + bufWidth + idX) / *(olpCnt + idX);
	      gVB = *(olpBuf.inp + (2 * bufWidth) + idX) / *(olpCnt + idX);
	      gVA = *(olpBuf.inp + (3 * bufWidth) + idX) / *(olpCnt + idX);
	      WLZ_RGBA_RGBA_SET(*(dGP.rgbp + idX), gVR, gVG, gVB, gVA);
	    }
	    else
	    {
	      *(dGP.rgbp + idX) = bgdV.v.rgbv;
	    }
	  }
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
      break;
    case WLZ_INTERPOLATION_CLASSIFY_1:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
    default:
      errNum = WLZ_ERR_INTERPOLATION_TYPE;
      break;
  }
  return(errNum);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzTransform
* \brief        Gets the nodes, node displacements and edges of the mesh.
*               The nodes and node displacements are returned as simple
*               arrays of vertices. The edges are returned as a degenerate
*               list of triples, with each triple being the indices of the
*               nodes of an element.
* \param        mObj                  	Given mesh object.
* \param        dstNNod                 Destination pointer for the number of
*                                       mesh nodes.
* \param        dstNod                  Destination pointer for the mesh nodes.
* \param        dstNDsp                 Destination pointer for the number of
*                                       mesh node displacements.
* \param        dstDsp                  Destination pointer for the mesh node
*                                       displacement.
* \param        dstNEdg                 Destination pointer for the number of
*                                       edge indices.
* \param        dstEdg                  Destination pointer for the edge
*                                       indices.
*/
WlzErrorNum     WlzCMeshGetNodesAndEdges(WlzObject *mObj,
                                int *dstNNod, WlzDVertex2 **dstNod,
                                int *dstNDsp, WlzDVertex2 **dstDsp,
                                int *dstNEdg, int **dstEdg)
{
  int           nId0,
  		nId1,
		eId0,
		eId1,
		tblSz,
		nNod = 0,
                nEdg = 0;
  WlzDVertex2   *nod = NULL,
                *dsp = NULL;
  int		*tbl = NULL,
  		*edg = NULL;
  WlzCMesh2D	*mesh;
  WlzCMeshNod2D	*mNod;
  WlzCMeshElm2D *mElm;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mObj->domain.core->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((dstNNod == NULL) || (dstNod == NULL) ||
          (dstNDsp == NULL) || (dstDsp == NULL) ||
	  (dstNEdg == NULL) || (dstEdg == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    mesh = mObj->domain.cm2;
    if((mesh->res.nod.maxEnt > 0) && (mesh->res.elm.maxEnt > 0))
    {
      /* Make a table for skipping out invalid node and element id's
       * along with arrays for the nodes, displacements and the edge
       * indices. */
      tblSz = (mesh->res.nod.maxEnt > mesh->res.elm.maxEnt)?
              mesh->res.nod.maxEnt: mesh->res.elm.maxEnt;
      if(((tbl = (int *)
                 AlcMalloc(tblSz * sizeof(int))) == NULL) ||
	 ((nod = (WlzDVertex2 *)
		 AlcMalloc(mesh->res.nod.maxEnt *
		           sizeof(WlzDVertex2))) == NULL) ||
	 ((dsp = (WlzDVertex2 *)
		 AlcMalloc(mesh->res.nod.maxEnt *
		           sizeof(WlzDVertex2))) == NULL) ||
	 ((edg = (int *)
		 AlcMalloc(mesh->res.elm.maxEnt * 3 * sizeof(int))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        /* Skip deleted nodes while building the node table and inserting
	 * the nodes and displacements into their arrays. */
	nId0 = 0;
	while(nId0 < mesh->res.nod.maxEnt)
	{
	  mNod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, nId0);
	  if(mNod->idx >= 0)
	  {
	    tbl[nId0] = nId0;
	    nod[nId0] = mNod->pos;
            dsp[nId0] = *(WlzDVertex2 *)(mNod->prop);
	    ++nId0;
	  }
	  else
	  {
	    break;
	  }
	}
	nId1 = nId0 + 1;
	while(nId1 < mesh->res.nod.maxEnt)
	{
	  tbl[nId1] = nId0;
	  mNod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, nId1);
	  if(mNod->idx >= 0)
	  {
	    nod[nId0] = mNod->pos;
	    dsp[nId0] = *(WlzDVertex2 *)(mNod->prop);
	    ++nId0;
	  }
	  ++nId1;
	}
	nNod = nId0;
	eId0 = 0;
	eId1 = 0;
	while(eId0 < mesh->res.elm.maxEnt)
	{
          mElm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, eId0);
	  if(mElm->idx >= 0)
	  {
	    edg[eId1++] = tbl[mElm->edu[0].nod->idx];
	    edg[eId1++] = tbl[mElm->edu[1].nod->idx];
	    edg[eId1++] = tbl[mElm->edu[2].nod->idx];
	  }
	  ++eId0;
	}
	nEdg = eId1;
      }
      AlcFree(tbl);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNNod = nNod;
    *dstNod = nod;
    *dstNDsp = nNod;
    *dstDsp = dsp;
    *dstNEdg = nEdg;
    *dstEdg = edg;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Scales index values pointed by ixv.
                If ixcSrc is NULL then scales the values
                in place, pointed by ixv. If it not NULL,
                then copies and scales values of ixvSrc
* \author	Zsolt Husz
* \param	scale                   Scale value
* \param	ixv                     Pointer to the array of values
* \param	idV                     Lentgh of index valye array
* \param	ixvSrc                  Pointer to the array of source values.
*/
static WlzErrorNum WlzScaleIndexedVal(double scale,
                                      WlzIndexedValues *ixv,
                                      int idV, WlzIndexedValues *ixcSrc)
{
  int           idX,
                vCnt;
  WlzGreyP      gP;
  WlzGreyP      gPOld;
  WlzUByte      rgba[4];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(ixv == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    vCnt = 1;
    if(ixv->rank > 0)
    {
      for(idX = 0; idX < ixv->rank; ++idX)
      {
	vCnt *= ixv->dim[idX];
      }
    }
    gP.v = WlzIndexedValueGet(ixv, idV);
    if (ixcSrc)
    {
      gPOld.v = WlzIndexedValueGet(ixcSrc, idV);
    }
    else
    {
      gPOld = gP;
    }
    for(idX = 0; idX < vCnt; ++idX)
    {
      switch(ixv->vType)
      {
	case WLZ_GREY_LONG:
	  gP.lnp[idX] = round(gPOld.lnp[idX] * scale);
	  break;
	case WLZ_GREY_INT:
	  gP.inp[idX] = round(gPOld.inp[idX] * scale);
	  break;
	case WLZ_GREY_SHORT:
	  gP.shp[idX] = round(gPOld.shp[idX] * scale);
	  break;
	case WLZ_GREY_UBYTE:
	  gP.ubp[idX] = round(gPOld.ubp[idX] * scale);
	  break;
	case WLZ_GREY_FLOAT:
	  gP.flp[idX] = (float)(gPOld.flp[idX] * scale);
	  break;
	case WLZ_GREY_DOUBLE:
	  gP.dbp[idX] = gPOld.dbp[idX] * scale;
	  break;
	case WLZ_GREY_RGBA:
	  rgba[0] = WLZ_RGBA_RED_GET(gPOld.rgbp[idX]);
	  rgba[1] = WLZ_RGBA_GREEN_GET(gPOld.rgbp[idX]);
	  rgba[2] = WLZ_RGBA_BLUE_GET(gPOld.rgbp[idX]);
	  rgba[3] = WLZ_RGBA_ALPHA_GET(gPOld.rgbp[idX]);
	  rgba[0] = round(rgba[0] * scale);
	  rgba[1] = round(rgba[1] * scale);
	  rgba[2] = round(rgba[2] * scale);
	  rgba[3] = round(rgba[3] * scale);
	  WLZ_RGBA_RED_SET(gP.rgbp[idX] ,rgba[0]);
	  WLZ_RGBA_GREEN_SET(gP.rgbp[idX] ,rgba[1]);
	  WLZ_RGBA_BLUE_SET(gP.rgbp[idX] ,rgba[2]);
	  WLZ_RGBA_ALPHA_SET(gP.rgbp[idX] ,rgba[3]);
	  break;
	default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  return(errNum);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Scales index values pointed of a WoolzObject.
                If ixcSrc is NULL then scales the values
                in place. If it not NULL,
                then copies and scales values of ixvSrc
* \author	Zsolt Husz
* \param	obj                     Woolz object
* \param	scale                   Scale value
* \param	ixvSrc                  Pointer to the array of source values
*/
static WlzErrorNum WlzScaleCMeshValueNodOrElem(WlzObject *obj, double scale,
                                               WlzIndexedValues *ixcSrc)
{
  WlzCMeshP     mesh;
  WlzCMeshEntP  ent;
  WlzIndexedValues *ixv;
  int           idx;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

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
  else
  {
    switch(obj->domain.core->type)
    {
      case WLZ_CMESH_2D:
	mesh.m2 = obj->domain.cm2;
	ixv = obj->values.x;
	switch(ixv->attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    for(idx = 0; idx < mesh.m2->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m2->res.nod.vec, idx);
	      if(ent.n2->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    for(idx = 0; idx < mesh.m2->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m2->res.elm.vec, idx);
	      if(ent.e2->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_CMESH_2D5:
	mesh.m2d5 = obj->domain.cm2d5;
	ixv = obj->values.x;
	switch(ixv->attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    for(idx = 0; idx < mesh.m2d5->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m2d5->res.nod.vec, idx);
	      if(ent.n2d5->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    for(idx = 0; idx < mesh.m2d5->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m2d5->res.elm.vec, idx);
	      if(ent.e2d5->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      case WLZ_CMESH_3D:
	mesh.m3 = obj->domain.cm3;
	ixv = obj->values.x;
	switch(ixv->attach)
	{
	  case WLZ_VALUE_ATTACH_NOD:
	    for(idx = 0; idx < mesh.m3->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m3->res.nod.vec, idx);
	      if(ent.n3->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  case WLZ_VALUE_ATTACH_ELM:
	    for(idx = 0; idx < mesh.m3->res.nod.maxEnt; ++idx)
	    {
	      ent.v = AlcVectorItemGet(mesh.m3->res.elm.vec, idx);
	      if(ent.e3->idx >= 0)
	      {
		if((errNum = WlzScaleIndexedVal(scale, ixv,
			                        idx, ixcSrc)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_VALUES_TYPE;
	    break;
	}
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	The bounding box of the mesh in the mesh transform.
* \ingroup	WlzTransform
* \brief	Computes the bounding box of the mesh in the given mesh
* 		transform, with or without applying the displacements
* 		according to the value of trans.
* \param	mObj			Given mesh transform object.
* \param	trans			Displacements applied if non-zero.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox2	WlzCMeshTransformGetBBox2D(WlzObject *mObj,
					int trans, WlzErrorNum *dstErr)
{
  int		idN,
  		first = 1;
  WlzDVertex2	pos;
  double	*dsp;
  WlzDBox2	bBox;
  WlzCMesh2D	*mesh;
  WlzCMeshNod2D	*nod;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mObj->domain.core->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    mesh = mObj->domain.cm2;
    if(mesh->res.elm.numEnt <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ixv = mObj->values.x;
    if(ixv == NULL)
    {
      trans = 0;
    }
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	pos = nod->pos;
	if(trans)
	{
	  dsp = (double *)WlzIndexedValueGet(ixv, idN);
	  pos.vtX += dsp[0];
	  pos.vtY += dsp[1];
	}
	if(first)
	{
	  bBox.xMin = bBox.xMax = pos.vtX;
	  bBox.yMin = bBox.yMax = pos.vtY;
	  first = 0;
	}
	else
	{
	  if(pos.vtX < bBox.xMin)
	  {
	    bBox.xMin = pos.vtX;
	  }
	  else if(pos.vtX > bBox.xMax)
	  {
	    bBox.xMax = pos.vtX;
	  }
	  if(pos.vtY < bBox.yMin)
	  {
	    bBox.yMin = pos.vtY;
	  }
	  else if(pos.vtY > bBox.yMax)
	  {
	    bBox.yMax = pos.vtY;
	  }
	}
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(bBox);
}

/*!
* \return	The bounding box of the mesh in the mesh transform.
* \ingroup	WlzTransform
* \brief	Computes the bounding box of the mesh in the given mesh
* 		transform, with or without applying the displacements
* 		according to the value of trans.
* \param	mObj			Given mesh transform object.
* \param	trans			Displacements applied if non-zero.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzDBox3	WlzCMeshTransformGetBBox3D(WlzObject *mObj,
					int trans, WlzErrorNum *dstErr)
{
  int		idN,
  		first = 1;
  double	*dsp;
  WlzDVertex3	pos;
  WlzDBox3	bBox;
  WlzCMesh3D	*mesh;
  WlzCMeshNod3D	*nod;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mObj->domain.core->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    mesh = mObj->domain.cm3;
    if(mesh->res.elm.numEnt <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    ixv = mObj->values.x;
    if(ixv == NULL)
    {
      trans = 0;
    }
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	pos = nod->pos;
	if(trans)
	{
	  dsp = (double *)WlzIndexedValueGet(ixv, idN);
	  pos.vtX += dsp[0];
	  pos.vtY += dsp[1];
	  pos.vtZ += dsp[2];
	}
	if(first != 0)
	{
	  bBox.xMin = bBox.xMax = pos.vtX;
	  bBox.yMin = bBox.yMax = pos.vtY;
	  bBox.zMin = bBox.zMax = pos.vtZ;
	  first = 0;
	}
	else
	{
	  if(pos.vtX < bBox.xMin)
	  {
	    bBox.xMin = pos.vtX;
	  }
	  else if(pos.vtX > bBox.xMax)
	  {
	    bBox.xMax = pos.vtX;
	  }
	  if(pos.vtY < bBox.yMin)
	  {
	    bBox.yMin = pos.vtY;
	  }
	  else if(pos.vtY > bBox.yMax)
	  {
	    bBox.yMax = pos.vtY;
	  }
	  if(pos.vtZ < bBox.zMin)
	  {
	    bBox.zMin = pos.vtZ;
	  }
	  else if(pos.vtZ > bBox.zMax)
	  {
	    bBox.zMax = pos.vtZ;
	  }
	}
      }
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(bBox);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Scales index values pointed of a WoolzObject
                in place.
* \author	Zsolt Husz
* \param	scale                   Scale value
* \param	obj                     Woolz object
*/
WlzErrorNum WlzScaleCMeshValue(double scale, WlzObject *obj)
{
  return(WlzScaleCMeshValueNodOrElem(obj, scale, NULL));
}

/*!
* \return	Copied source object.
* \ingroup	WlzTransform
* \brief	Creates an woolz CMesh transform object indentical
                with the input object, but with value scaled.
* \author	Zsolt Husz
* \param	scale                   Scale value
* \param	obj                     Woolz object
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzCopyScaleCMeshValue(double scale, WlzObject *obj,
				        WlzErrorNum *dstErr)
{
  WlzIndexedValues *ixcSrc = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  WlzObject     *scaledObj = NULL;
  WlzValues     values;


  values.x = NULL;
  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((ixcSrc = obj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    scaledObj = WlzMakeMain(obj->type, obj->domain, values, NULL, NULL,
    			    &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    values.x = WlzMakeIndexedValues(scaledObj, ixcSrc->rank, ixcSrc->dim,
                                    ixcSrc->vType, ixcSrc->attach, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    scaledObj->values = WlzAssignValues(values, NULL);
    errNum = WlzScaleCMeshValueNodOrElem(scaledObj, scale, ixcSrc);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(scaledObj);
    scaledObj = NULL;
  }
  if (dstErr)
  {
    *dstErr = errNum;
  }
  return(scaledObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the product of the given affine and conforming mesh
* 		transforms in place, ie the mesh transform has it's
* 		displacements overwritten.
*
* 		Given a conforming mesh transform \f$\mathbf{M}\f$ and an
* 		affine transform \f$\mathbf{A}\f$. The product
* 		\f$\mathbf{P}\f$ can be evaluated as:
*               \f[
                \begin{array}{ll}
                \mathbf{P} = \mathbf{A}\mathbf{M}, & o = 0 \\
                \mathbf{P} = \mathbf{M}\mathbf{A}, & o = 1
                \end{array}
                \f]
*               where \f$o\f$ is the order parameter and
*               \f[
                \mathbf{T_0}\mathbf{T_1}\mathbf{x} =
                \mathbf{T_0}(\mathbf{T_1}\mathbf{x})
                \f]
*               The product with \f$o = 1\f$ applies the mesh displacement
*               at location \f$\mathbf{x}\f$ and not at
*               \f$\mathbf{A}(\mathbf{x})\f$ as might be expected.
* \param	trM			Given conforming mesh transform.
* \param	trA			Given affine transform.
* \param	order			Order of evaluation.
*/
WlzErrorNum	WlzCMeshAffineProduct(WlzObject *trM, WlzAffineTransform *trA,
				      int order)
{
  WlzTransformType aType = WLZ_TRANSFORM_EMPTY;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((trM == NULL) || (trA == NULL) ||
     (trM->domain.core == NULL) || (trM->values.core == NULL))
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else if(trM->values.core->type != (WlzObjectType )WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    switch(trA->type)
    {
      case WLZ_TRANSFORM_2D_AFFINE:
      case WLZ_TRANSFORM_2D_REG:
      case WLZ_TRANSFORM_2D_TRANS:
      case WLZ_TRANSFORM_2D_NOSHEAR:
        aType = WLZ_TRANSFORM_2D_AFFINE;
	break;
      case WLZ_TRANSFORM_3D_AFFINE:
      case WLZ_TRANSFORM_3D_REG:
      case WLZ_TRANSFORM_3D_TRANS:
      case WLZ_TRANSFORM_3D_NOSHEAR:
        aType = WLZ_TRANSFORM_3D_AFFINE;
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(trM->type)
    {
      case WLZ_TRANSFORM_2D_CMESH:
	if(aType != WLZ_TRANSFORM_2D_AFFINE)
	{
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	}
	else
	{
	  errNum = WlzCMeshAffineProduct2D(trM, trA, order);
	}
        break;
      case WLZ_TRANSFORM_2D5_CMESH:
	if(aType != WLZ_TRANSFORM_3D_AFFINE)
	{
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	}
	else
	{
	  errNum = WlzCMeshAffineProduct2D5(trM, trA, order);
	}
        break;
      case WLZ_TRANSFORM_3D_CMESH:
	if(aType != WLZ_TRANSFORM_3D_AFFINE)
	{
	  errNum = WLZ_ERR_TRANSFORM_TYPE;
	}
	else
	{
	  errNum = WlzCMeshAffineProduct3D(trM, trA, order);
	}
        break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the product of the given affine and conforming mesh
* 		transforms in place. See WlzCMeshAffineProduct().
* \param	trM			Given conforming mesh transform.
* \param	trA			Given affine transform.
* \param	order			Order of evaluation.
*/
static WlzErrorNum WlzCMeshAffineProduct2D(WlzObject *trM,
					   WlzAffineTransform *trA, int order)
{
  WlzIndexedValues *ixv;
  WlzCMesh2D	*mesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  ixv = trM->values.x;
  if((ixv->attach != WLZ_VALUE_ATTACH_NOD) ||
     (ixv->rank != 1) ||
     (ixv->dim[0] < 2) ||
     (ixv->vType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    int		maxNod;
    AlcVector	*nVec;

    mesh = trM->domain.cm2;
    nVec = mesh->res.nod.vec;
    maxNod = mesh->res.nod.maxEnt;
    if(order == 0) /* P(x) = A(M(x))*/
    {
      int         idN;

      for(idN = 0; idN < maxNod; ++idN)
      {
	WlzCMeshNod2D *nod;

	nod = (WlzCMeshNod2D *)AlcVectorItemGet(nVec, idN);
	if(nod->idx >= 0)
	{
	  double      *dsp;
	  WlzDVertex2 dPos;

	  dsp = (double *)WlzIndexedValueGet(ixv, idN);
	  dPos.vtX = nod->pos.vtX + dsp[0];
	  dPos.vtY = nod->pos.vtY + dsp[1];
	  dPos = WlzAffineTransformVertexD2(trA, dPos, &errNum);
	  dsp[0] = dPos.vtX - nod->pos.vtX;
	  dsp[1] = dPos.vtY - nod->pos.vtY;
	}
      }
    }
    else /* P(x) = M(A(x)) */
    {
      WlzAffineTransform *trI;

      trI = WlzAffineTransformInverse(trA, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int         idN;

	for(idN = 0; idN < maxNod; ++idN)
	{
	  WlzCMeshNod2D *nod;

	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(nVec, idN);
	  if(nod->idx >= 0)
	  {
	    double	*dsp;
	    WlzDVertex2	iPos;

	    iPos = WlzAffineTransformVertexD2(trI, nod->pos, NULL);
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dsp[0] += nod->pos.vtX - iPos.vtX;
	    dsp[1] += nod->pos.vtY - iPos.vtY;
	    nod->pos = iPos;
	  }
	}
      }
      (void )WlzFreeAffineTransform(trI);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshUpdateBBox2D(mesh);
	WlzCMeshUpdateMaxSqEdgLen2D(mesh);
	errNum = WlzCMeshReassignGridCells2D(mesh, 0);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the product of the given affine and conforming mesh
* 		transforms in place. See WlzCMeshAffineProduct().
* \param	trM			Given conforming mesh transform.
* \param	trA			Given affine transform.
* \param	order			Order of evaluation.
*/
static WlzErrorNum WlzCMeshAffineProduct2D5(WlzObject *trM,
					   WlzAffineTransform *trA, int order)
{
  WlzIndexedValues *ixv;
  WlzCMesh2D5	*mesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  ixv = trM->values.x;
  if((ixv->attach != WLZ_VALUE_ATTACH_NOD) ||
     (ixv->rank != 1) ||
     (ixv->dim[0] < 3) ||
     (ixv->vType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    int		maxNod;
    AlcVector	*nVec;

    mesh = trM->domain.cm2d5;
    nVec = mesh->res.nod.vec;
    maxNod = mesh->res.nod.maxEnt;
    if(order == 0) /* P(x) = A(M(x))*/
    {
      int         idN;

      for(idN = 0; idN < maxNod; ++idN)
      {
	WlzCMeshNod2D5 *nod;

	nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nVec, idN);
	if(nod->idx >= 0)
	{
	  double      *dsp;
	  WlzDVertex3 dPos;

	  dsp = (double *)WlzIndexedValueGet(ixv, idN);
	  dPos.vtX = nod->pos.vtX + dsp[0];
	  dPos.vtY = nod->pos.vtY + dsp[1];
	  dPos.vtZ = nod->pos.vtZ + dsp[2];
	  dPos = WlzAffineTransformVertexD3(trA, dPos, &errNum);
	  dsp[0] = dPos.vtX - nod->pos.vtX;
	  dsp[1] = dPos.vtY - nod->pos.vtY;
	  dsp[2] = dPos.vtZ - nod->pos.vtZ;
	}
      }
    }
    else /* P(x) = M(A(x)) */
    {
      WlzAffineTransform *trI;

      trI = WlzAffineTransformInverse(trA, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int         idN;

	for(idN = 0; idN < maxNod; ++idN)
	{
	  WlzCMeshNod2D5 *nod;

	  nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nVec, idN);
	  if(nod->idx >= 0)
	  {
	    double	*dsp;
	    WlzDVertex3	iPos;

	    iPos = WlzAffineTransformVertexD3(trI, nod->pos, NULL);
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dsp[0] += nod->pos.vtX - iPos.vtX;
	    dsp[1] += nod->pos.vtY - iPos.vtY;
	    dsp[2] += nod->pos.vtZ - iPos.vtZ;
	    nod->pos = iPos;
	  }
	}
      }
      (void )WlzFreeAffineTransform(trI);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshUpdateBBox2D5(mesh);
	WlzCMeshUpdateMaxSqEdgLen2D5(mesh);
	errNum = WlzCMeshReassignGridCells2D5(mesh, 0);
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the product of the given affine and conforming mesh
* 		transforms in place. See WlzCMeshAffineProduct().
* \param	trM			Given conforming mesh transform.
* \param	trA			Given affine transform.
* \param	order			Order of evaluation.
*/
static WlzErrorNum WlzCMeshAffineProduct3D(WlzObject *trM,
					   WlzAffineTransform *trA, int order)
{
  WlzIndexedValues *ixv;
  WlzCMesh3D	*mesh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  ixv = trM->values.x;
  if((ixv->attach != WLZ_VALUE_ATTACH_NOD) ||
     (ixv->rank != 1) ||
     (ixv->dim[0] < 3) ||
     (ixv->vType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    int		maxNod;
    AlcVector	*nVec;

    mesh = trM->domain.cm3;
    nVec = mesh->res.nod.vec;
    maxNod = mesh->res.nod.maxEnt;
    if(order == 0) /* P(x) = A(M(x))*/
    {
      int         idN;

      for(idN = 0; idN < maxNod; ++idN)
      {
	WlzCMeshNod3D *nod;

	nod = (WlzCMeshNod3D *)AlcVectorItemGet(nVec, idN);
	if(nod->idx >= 0)
	{
	  double      *dsp;
	  WlzDVertex3 dPos;

	  dsp = (double *)WlzIndexedValueGet(ixv, idN);
	  dPos.vtX = nod->pos.vtX + dsp[0];
	  dPos.vtY = nod->pos.vtY + dsp[1];
	  dPos.vtZ = nod->pos.vtZ + dsp[2];
	  dPos = WlzAffineTransformVertexD3(trA, dPos, &errNum);
	  dsp[0] = dPos.vtX - nod->pos.vtX;
	  dsp[1] = dPos.vtY - nod->pos.vtY;
	  dsp[2] = dPos.vtZ - nod->pos.vtZ;
	}
      }
    }
    else /* P(x) = M(A(x)) */
    {
      WlzAffineTransform *trI;

      trI = WlzAffineTransformInverse(trA, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int         idN;

	for(idN = 0; idN < maxNod; ++idN)
	{
	  WlzCMeshNod3D *nod;

	  nod = (WlzCMeshNod3D *)AlcVectorItemGet(nVec, idN);
	  if(nod->idx >= 0)
	  {
	    double	*dsp;
	    WlzDVertex3	iPos;

	    iPos = WlzAffineTransformVertexD3(trI, nod->pos, NULL);
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dsp[0] += nod->pos.vtX - iPos.vtX;
	    dsp[1] += nod->pos.vtY - iPos.vtY;
	    dsp[2] += nod->pos.vtZ - iPos.vtZ;
	    nod->pos = iPos;
	  }
	}
      }
      (void )WlzFreeAffineTransform(trI);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzCMeshUpdateBBox3D(mesh);
	WlzCMeshUpdateMaxSqEdgLen3D(mesh);
	errNum = WlzCMeshReassignGridCells3D(mesh, 0);
      }
    }
  }
  return(errNum);
}

/*!
* \return	New conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Computes the product of the two given (convex) mesh
* 		transforms. This is computed within intersection of
* 		the two mesh transforms resulting in a conforming
* 		mesh transform.
* 		\f[
 		\mathbf{T_R}(\mathbf{x}) =
		   \mathbf{T_0}(\mathbf{T_1}(\mathbf{x}))
 		\f]
* 		Where possible the node positions of the second mesh
* 		\f$\mathbf{T_1}\f$ are preserved in the output conforming
* 		mesh \f$\mathbf{T_R}\f$.
* 		The displacements are given in the output conforming
* 		transform are given by
* 		\f[
		\mathbf{d_R}(\mathbf{x}) = \mathbf{d_1}(\mathbf{x}) +
		  \mathbf{d_0}(\mathbf{d_1}(\mathbf{x})) - \mathbf{x}
  		\f]
* \param	tr0			First (convex) mesh transform.
* \param	tr1			Second (convex) mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshMeshMeshProduct(WlzMeshTransform *tr0,
					 WlzMeshTransform *tr1,
					 WlzErrorNum *dstErr)
{
  int		nNbr0 = 0,
  		nNbr1 = 0,
		maxKrigBuf = 0,
		maxNbrIdxBuf = 0;
  double	dRange;
  int		*wSp = NULL,
     		*nbrIdxBuf = NULL;
  double	*posSV = NULL;
  WlzIndexedValues *ixvR = NULL;
  WlzCMesh2D	*meshR = NULL;
  WlzObject	*trR = NULL;
  int		*nodTab = NULL;
  WlzKrigModelFn modelFn;
  AlgMatrix	modelSV;
  WlzDVertex2	*nbrPosBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  modelSV.core = NULL;
  if((tr0 == NULL) || (tr1 == NULL))
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else if((tr0->type != WLZ_TRANSFORM_2D_MESH) || (tr0->type != tr1->type))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    dRange = WlzMeshMaxEdgeLenSq(tr0, &errNum);
    if((errNum == WLZ_ERR_NONE) && (dRange < 1.0))
    {
      errNum = WLZ_ERR_TRANSFORM_DATA;
    }
    else
    {
      dRange = sqrt(dRange);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create a new conforming mesh with nodes and elements of corresponding
     * to those of tr1 and where the displaced nodes of tr1 fall within
     * elements of tr0. */
    meshR = WlzCMeshIntersect2Mesh2D(tr1, tr0, 1, &nodTab, &errNum);
  }
  /* Create a constrained mesh transform from the mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain dom;
    WlzValues val;

    dom.cm2 = meshR;
    val.core = NULL;
    trR = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 2;

    ixvR = WlzMakeIndexedValues(trR, 1, &dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  /* Set displacements for the new constrained mesh transform. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;

    for(idN = 0; idN < tr1->maxNodes; ++idN)
    {
      int	nIdx,
      		lastElm = -1;

      if((nIdx = nodTab[idN]) >= 0)
      {
	int	      idE0;
	WlzDVertex2   v0,
		      v1;
        WlzMeshNode   *nod1;
	double	      *dspR;
        WlzMeshNode   *nod0[3];
        WlzCMeshNod2D *nodR;
	WlzMeshElem   *elm0;

	nod1 = &(tr1->nodes[idN]);
        nodR = (WlzCMeshNod2D *)AlcVectorItemGet(meshR->res.nod.vec, nIdx);
	WLZ_VTX_2_ADD(v0, nod1->position, nod1->displacement);
	/* Find element in tr0 which encloses a vertex at the position of
	 * the new node. */
        idE0 = WlzMeshElemFindVx(tr0, v0, lastElm, &lastElm, NULL, &errNum);
	if(idE0 >= 0)
	{
	  /* Interpolate displacement at the new node position using
	   * barycentric interpolation. */
	  elm0 = &(tr0->elements[idE0]);
	  nod0[0] = &(tr0->nodes[elm0->nodes[0]]);
	  nod0[1] = &(tr0->nodes[elm0->nodes[1]]);
	  nod0[2] = &(tr0->nodes[elm0->nodes[2]]);
	  v1.vtX = WlzGeomInterpolateTri2D(nod0[0]->position,
					   nod0[1]->position,
					   nod0[2]->position,
					   nod0[0]->displacement.vtX,
					   nod0[1]->displacement.vtX,
					   nod0[2]->displacement.vtX,
					   v0);
	  v1.vtY = WlzGeomInterpolateTri2D(nod0[0]->position,
					   nod0[1]->position,
					   nod0[2]->position,
					   nod0[0]->displacement.vtY,
					   nod0[1]->displacement.vtY,
					   nod0[2]->displacement.vtY,
					   v0);
	}
	else
	{
	  /* Node is not in an element of tr0, but the last element visited
	   * should be closest for a convex mesh. Find nearby nodes using
	   * this element. */
          elm0 = &(tr0->elements[lastElm]);
	  if(maxNbrIdxBuf < 6)
	  {
	    maxNbrIdxBuf = 6;
	    if((nbrIdxBuf = (int *)
	                    AlcRealloc(nbrIdxBuf,
			               sizeof(int) * maxNbrIdxBuf)) == NULL)
            {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int 	id1,
	    		id2;
	    double	dSq,
	    		dSqMin;
	    WlzDVertex2 p,
	      		del;
            int		nbrFlgTbl[3];

	    id2 = 0;
	    nbrFlgTbl[0] = WLZ_MESH_ELEM_FLAGS_NBR_0;
	    nbrFlgTbl[1] = WLZ_MESH_ELEM_FLAGS_NBR_1;
	    nbrFlgTbl[2] = WLZ_MESH_ELEM_FLAGS_NBR_2;
	    /* Put closest node of element into the neighbour buffer first. */
	    p = tr0->nodes[elm0->nodes[0]].position;
	    WLZ_VTX_2_SUB(del, p, v0);
	    dSqMin = WLZ_VTX_2_SQRLEN(del);
	    for(id1 = 1; id1 < 3; ++id1)
	    {
	      p = tr0->nodes[elm0->nodes[id1]].position;
	      WLZ_VTX_2_SUB(del, p, v0);
	      dSq = WLZ_VTX_2_SQRLEN(del);
	      if(dSq < dSqMin)
	      {
	        id2 = id1;
	      }
	    }
	    nNbr1 = 3;
	    nbrIdxBuf[0] = elm0->nodes[id2];
	    /* Add remaining nodes of the closest element to the neighbour
	     * buffer. */
	    nbrIdxBuf[1] = elm0->nodes[(id2 + 1) % 3];
	    nbrIdxBuf[2] = elm0->nodes[(id2 + 2) % 3];
	    /* Add nodes of the elements neighbours where the neighbours
	     * exist and their nodes are not already in the neighbours
	     * buffer. */
	    for(id1 = 0; id1 < 3; ++id1)
	    {
	      if((elm0->flags & nbrFlgTbl[id1]) != 0)
	      {
	        WlzMeshElem *nElm;

		nElm = &(tr0->elements[elm0->neighbours[id1]]);
		for(id2 = 0; id2 < 3; ++id2)
		{
		  int	id3,
		  	hit = 0;

		  for(id3 = 0; id3 < nNbr1; ++id3)
		  {
		    if(nElm->nodes[id2] == nbrIdxBuf[id3])
		    {
		      hit = 1;
		      break;
		    }
		  }
		  if(hit == 0)
		  {
		    nbrIdxBuf[nNbr1++] = nElm->nodes[id2];
		  }
		}
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    /* Reallocate buffers if required. */
	    errNum = WlzKrigReallocBuffers2D(&nbrPosBuf, &posSV, &wSp,
					     &modelSV, &maxKrigBuf,
					     nNbr1, nNbr0);
	    nNbr0 = nNbr1;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    for(i = 0; i < nNbr1; ++i)
	    {
	      nbrPosBuf[i] = tr0->nodes[nbrIdxBuf[i]].position;
	    }
	    WlzKrigSetModelFn(&modelFn, WLZ_KRIG_MODELFN_LINEAR,
			      0.0, 0.1, 2.0 * dRange);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetModelSV2D(modelSV, &modelFn, nNbr1, nbrPosBuf,
					  wSp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetPosSV2D(posSV, &modelFn, nNbr1, nbrPosBuf, v0);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    v1.vtX = 0.0;
	    v1.vtY = 0.0;
	    WlzKrigOWeightsSolve(modelSV, posSV, wSp, WLZ_MESH_TOLERANCE);
	    /* posSV now contains the weights. */
	    for(i = 0; i < nNbr1; ++i) 
	    {
	      WlzDVertex2 dsp0;

	      dsp0 = tr0->nodes[nbrIdxBuf[i]].displacement;
	      v1.vtX += posSV[i] * dsp0.vtX;
	      v1.vtY += posSV[i] * dsp0.vtY;
	    }
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WLZ_VTX_2_ADD(v1, v1, v0);
	  WLZ_VTX_2_SUB(v1, v1, nodR->pos);
	  dspR = (double *)WlzIndexedValueGet(ixvR, nodR->idx);
	  dspR[0] = v1.vtX;
	  dspR[1] = v1.vtY;
	}
      }
    }
  }
  AlgMatrixFree(modelSV);
  AlcFree(wSp);
  AlcFree(posSV);
  AlcFree(nbrIdxBuf);
  AlcFree(nbrPosBuf);
  AlcFree(nodTab);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trR);
}

/*!
* \return	New Woolz object containing the conforming mesh transform
* 		product or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes the product of the two given (convex) mesh
* 		transforms. This is computed within intersection of
* 		the two mesh transforms resulting in a conforming
* 		mesh transform.
* 		\f[
 		\mathbf{T_R}(\mathbf{x}) =
		   \mathbf{T_1}(\mathbf{T_0}(\mathbf{x}))
 		\f]
* 		Where possible the node positions of the second mesh
* 		\f$\mathbf{T_1}\f$ are preserved in the output mesh
* 		\f$\mathbf{T_R}\f$.
* 		The displacements in the output transform are given by
* 		\f[
		\mathbf{d_R}(\mathbf{x}) = \mathbf{d_0}(\mathbf{x}) +
		  \mathbf{d_1}(\mathbf{d_0}(\mathbf{x})) - \mathbf{x}
  		\f]
* \param	tr0			First (convex) mesh transform.
* \param	tr1			Second (convex) mesh transform.
* \param	order			Order of evaluation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshMeshProduct(WlzObject *tr0, WlzMeshTransform *tr1,
				     int order, WlzErrorNum *dstErr)
{
  WlzObject	*trR = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trR);
}

/*!
* \return	New Woolz object containing the conforming mesh transform
* 		product or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes the product of the two given (conforming) mesh
* 		transforms. This is computed within intersection of
* 		the two mesh transforms resulting in another conforming
* 		mesh transform.
* 		\f[
 		\mathbf{T_R}(\mathbf{x}) =
		   \mathbf{T_1}(\mathbf{T_0}(\mathbf{x}))
 		\f]
* 		Where possible the node positions of the second mesh
* 		\f$\mathbf{T_1}\f$ are preserved in the output mesh
* 		\f$\mathbf{T_R}\f$.
* 		The displacements in the output transform are given by
* 		\f[
		\mathbf{d_R}(\mathbf{x}) = \mathbf{d_0}(\mathbf{x}) +
		  \mathbf{d_1}(\mathbf{d_0}(\mathbf{x})) - \mathbf{x}
  		\f]
* \param	tr0			First (conforming) mesh transform.
* \param	tr1			Second (conforming) mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshProduct(WlzObject *tr0, WlzObject *tr1,
				 WlzErrorNum *dstErr)
{
  WlzObject	*trR = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((tr0 == NULL) || (tr1 == NULL))
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else if(tr0->type != tr1->type)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((tr0->domain.core == NULL) || (tr1->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((tr0->values.core == NULL) || (tr1->values.core == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(tr0->domain.core->type != tr1->domain.core->type)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(tr0->values.core->type != tr1->values.core->type)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    switch(tr0->type)
    {
      case WLZ_CMESH_2D:
	if(tr0->domain.cm2->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          trR = WlzCMeshProduct2D(tr0, tr1, &errNum);
	}
	break;
      case WLZ_CMESH_3D:
	if(tr0->domain.cm2->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          trR = WlzCMeshProduct3D(tr0, tr1, &errNum);
	}
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
  return(trR);
}

/*!
* \return	New Woolz object containing the conforming mesh transform
* 		product or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes the product of the two given (conforming) mesh
* 		transforms. See WlzCMeshProduct().
* \param	tr0			First (conforming) mesh transform.
* \param	tr1			Second (conforming) mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshProduct2D(WlzObject *tr0, WlzObject *tr1,
				    WlzErrorNum *dstErr)
{
  int		nNbr0 = 0,
  		nNbr1 = 0,
		maxKrigBuf = 0,
		maxNbrIdxBuf = 0;
  double	dRange;
  WlzObject	*trR = NULL;
  WlzCMesh2D	*mesh0 = NULL,
  		*mesh1 = NULL,
		*meshR = NULL;
  WlzIndexedValues *ixv0 = NULL,
  		   *ixv1 = NULL,
		   *ixvR = NULL;
  int		*wSp = NULL,
  		*nodTab = NULL,
  		*nbrIdxBuf = NULL;
  double	*posSV = NULL;
  WlzKrigModelFn modelFn;
  AlgMatrix	modelSV;
  WlzDVertex2	*nbrPosBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  modelSV.core = NULL;
  if(((mesh0 = tr0->domain.cm2) == NULL) ||
     ((mesh1 = tr1->domain.cm2) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(((ixv0 = tr0->values.x) == NULL) ||
          ((ixv1 = tr1->values.x) == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((mesh0->type != WLZ_CMESH_2D) || ( mesh0->type != mesh1->type))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv0->type != (WlzObjectType )WLZ_INDEXED_VALUES) ||
          ( ixv0->type != ixv1->type))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv0->rank != 1) ||
          (ixv1->rank != 1) ||
          (ixv0->dim[0] < 2) ||
	  (ixv1->dim[0] < 2) ||
          (ixv0->vType != WLZ_GREY_DOUBLE) ||
	  (ixv1->vType != WLZ_GREY_DOUBLE) ||
	  (ixv0->attach != WLZ_VALUE_ATTACH_NOD) ||
	  (ixv1->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    /* Create a new conforming mesh with nodes and elements corresponding
     * to those of tr0 and where the displaced nodes of tr0 fall within
     * elements of tr1. */
    trR = WlzCMeshIntersect(tr0, tr1, 1, &nodTab, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 2;

    meshR = trR->domain.cm2;
    ixvR = WlzMakeIndexedValues(trR, 1, &dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;
    WlzValues	val;

    val.x = ixvR;
    dRange = sqrt(mesh1->maxSqEdgLen);
    trR->values = WlzAssignValues(val, &errNum);
    /* Set displacements for the new constrained mesh transform. */
    for(idN = 0; idN < mesh0->res.nod.maxEnt; ++idN)
    {
      int	nIdx;
      WlzCMeshNod2D *nod1,
      		    *nodR;
      double	*dsp1;
      WlzDVertex2 v0,
      		  v1;

      WLZ_VTX_2_ZERO(v1);
      if((nIdx = nodTab[idN]) >= 0)
      {
	int	idE0;
        double  *dspR;

	nod1 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh0->res.nod.vec, idN);
        nodR = (WlzCMeshNod2D *)AlcVectorItemGet(meshR->res.nod.vec, nIdx);
	dsp1 = (double *)WlzIndexedValueGet(ixv0, idN);
	v0.vtX = nod1->pos.vtX + dsp1[0];
	v0.vtY = nod1->pos.vtY + dsp1[1];
	/* Find element in tr1 which encloses a vertex at the position of
	 * the new node. */
        idE0 = WlzCMeshElmEnclosingPos2D(mesh1, -1, v0.vtX, v0.vtY, 0, NULL);
	if(idE0 >= 0)
	{
	  double	*dsp0[3];
	  WlzCMeshNod2D *nod0[3];
	  WlzCMeshElm2D *elm0;

	  /* Interpolate displacement at the new node position using
	   * barycentric interpolation. */
	  elm0 = (WlzCMeshElm2D *)AlcVectorItemGet(mesh1->res.elm.vec, idE0);
	  nod0[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm0);
	  nod0[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm0);
	  nod0[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm0);
	  dsp0[0] = (double *)WlzIndexedValueGet(ixv1, nod0[0]->idx);
	  dsp0[1] = (double *)WlzIndexedValueGet(ixv1, nod0[1]->idx);
	  dsp0[2] = (double *)WlzIndexedValueGet(ixv1, nod0[2]->idx);
	  v1.vtX = WlzGeomInterpolateTri2D(nod0[0]->pos,
					   nod0[1]->pos,
					   nod0[2]->pos,
					   dsp0[0][0],
					   dsp0[1][0],
					   dsp0[2][0],
					   v0);
	  v1.vtY = WlzGeomInterpolateTri2D(nod0[0]->pos,
					   nod0[1]->pos,
					   nod0[2]->pos,
					   dsp0[0][1],
					   dsp0[1][1],
					   dsp0[2][1],
					   v0);
	}
	else /* idE0 <0, the vertex is not in the mesh. Find the closest
	      * node to the vertex in the mesh and then use kriging. */
	{
	  int		idN0;
          WlzCMeshNod2D	*nod0,
	  		*nodN;

	  idN0 = WlzCMeshClosestNod2D(mesh1, v0);
          nod0 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh1->res.nod.vec, idN0);
	  nNbr1 = WlzCMeshNodRingNodIndices2D(nod0, &maxNbrIdxBuf,
	  			              &nbrIdxBuf, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    /* Reallocate buffers if required. */
	    errNum = WlzKrigReallocBuffers2D(&nbrPosBuf, &posSV, &wSp,
					     &modelSV, &maxKrigBuf,
					     nNbr1, nNbr0);
	    nNbr0 = nNbr1;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    for(i = 0; i < nNbr1; ++i)
	    {
	      WlzCMeshNod2D *nod;
	      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh1->res.nod.vec,
						      nbrIdxBuf[i]);
	      nbrPosBuf[i] = nod->pos;
	    }
	    WlzKrigSetModelFn(&modelFn, WLZ_KRIG_MODELFN_LINEAR,
			      0.0, 0.1, 2.0 * dRange);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetModelSV2D(modelSV, &modelFn, nNbr1, nbrPosBuf,
					  wSp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetPosSV2D(posSV, &modelFn, nNbr1, nbrPosBuf, v0);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    v1.vtX = 0.0;
	    v1.vtY = 0.0;
	    WlzKrigOWeightsSolve(modelSV, posSV, wSp, WLZ_MESH_TOLERANCE);
	    /* posSV now contains the weights. */
	    for(i = 0; i < nNbr1; ++i) 
	    {
	      double *dsp0;

	      nodN = (WlzCMeshNod2D *)AlcVectorItemGet(mesh1->res.nod.vec,
	                                               nbrIdxBuf[i]);
	      dsp0 = (double *)WlzIndexedValueGet(ixv1, nodN->idx);
	      v1.vtX += posSV[i] * dsp0[0];
	      v1.vtY += posSV[i] * dsp0[1];
	    }
	  }
	}
	WLZ_VTX_2_ADD(v1, v1, v0);
	WLZ_VTX_2_SUB(v1, v1, nodR->pos);
	dspR = (double *)WlzIndexedValueGet(ixvR, nodR->idx);
	dspR[0] = v1.vtX;
	dspR[1] = v1.vtY;
      }
    }
  }
  AlcFree(wSp);
  AlcFree(posSV);
  AlcFree(nodTab);
  AlcFree(nbrIdxBuf);
  AlcFree(nbrPosBuf);
  AlgMatrixFree(modelSV);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trR);
}

/*!
* \return	New Woolz object containing the conforming mesh transform
* 		product or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes the product of the two given (conforming) mesh
* 		transforms. See WlzCMeshProduct().
* \param	tr0			First (conforming) mesh transform.
* \param	tr1			Second (conforming) mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshProduct3D(WlzObject *tr0, WlzObject *tr1,
				    WlzErrorNum *dstErr)
{
  int		nNbr0 = 0,
  		nNbr1 = 0,
		maxKrigBuf = 0,
		maxNbrIdxBuf = 0;
  double	dRange;
  WlzObject	*trR = NULL;
  WlzCMesh3D	*mesh0 = NULL,
  		*mesh1 = NULL,
		*meshR = NULL;
  WlzIndexedValues *ixv0 = NULL,
  		   *ixv1 = NULL,
		   *ixvR = NULL;
  int		*wSp = NULL,
  		*nodTab = NULL,
  		*nbrIdxBuf = NULL;
  double	*posSV = NULL;
  WlzKrigModelFn modelFn;
  AlgMatrix	modelSV;
  WlzDVertex3	*nbrPosBuf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  modelSV.core = NULL;
  if(((mesh0 = tr0->domain.cm3) == NULL) ||
     ((mesh1 = tr1->domain.cm3) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(((ixv0 = tr0->values.x) == NULL) ||
          ((ixv1 = tr1->values.x) == NULL))
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((mesh0->type != WLZ_CMESH_3D) || ( mesh0->type != mesh1->type))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv0->type != (WlzObjectType )WLZ_INDEXED_VALUES) ||
          ( ixv0->type != ixv1->type))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv0->rank != 1) ||
          (ixv1->rank != 1) ||
          (ixv0->dim[0] < 3) ||
	  (ixv1->dim[0] < 3) ||
          (ixv0->vType != WLZ_GREY_DOUBLE) ||
	  (ixv1->vType != WLZ_GREY_DOUBLE) ||
	  (ixv0->attach != WLZ_VALUE_ATTACH_NOD) ||
	  (ixv1->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    /* Create a new conforming mesh with nodes and elements corresponding
     * to those of tr0 and where the displaced nodes of tr0 fall within
     * elements of tr1. */
    trR = WlzCMeshIntersect(tr0, tr1, 1, &nodTab, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 3;

    meshR = trR->domain.cm3;
    ixvR = WlzMakeIndexedValues(trR, 1, &dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;
    WlzValues	val;

    val.x = ixvR;
    dRange = sqrt(mesh1->maxSqEdgLen);
    trR->values = WlzAssignValues(val, &errNum);
    /* Set displacements for the new constrained mesh transform. */
    for(idN = 0; idN < mesh0->res.nod.maxEnt; ++idN)
    {
      int	nIdx;
      WlzCMeshNod3D *nod1,
      		    *nodR;
      double	*dsp1;
      WlzDVertex3 v0,
      		  v1;

      WLZ_VTX_3_ZERO(v1);
      if((nIdx = nodTab[idN]) >= 0)
      {
	int	idE0;
        double  *dspR;

	nod1 = (WlzCMeshNod3D *)AlcVectorItemGet(mesh0->res.nod.vec, idN);
        nodR = (WlzCMeshNod3D *)AlcVectorItemGet(meshR->res.nod.vec, nIdx);
	dsp1 = (double *)WlzIndexedValueGet(ixv0, idN);
	v0.vtX = nod1->pos.vtX + dsp1[0];
	v0.vtY = nod1->pos.vtY + dsp1[1];
	v0.vtZ = nod1->pos.vtZ + dsp1[2];
	/* Find element in tr1 which encloses a vertex at the position of
	 * the new node. */
        idE0 = WlzCMeshElmEnclosingPos3D(mesh1, -1, v0.vtX, v0.vtY, v0.vtZ,
					 0, NULL);
	if(idE0 >= 0)
	{
	  double	*dsp0[4];
	  WlzCMeshNod3D *nod0[4];
	  WlzCMeshElm3D *elm0;

	  /* Interpolate displacement at the new node position using
	   * barycentric interpolation. */
	  elm0 = (WlzCMeshElm3D *)AlcVectorItemGet(mesh1->res.elm.vec, idE0);
	  nod0[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm0);
	  nod0[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm0);
	  nod0[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm0);
	  nod0[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm0);
	  dsp0[0] = (double *)WlzIndexedValueGet(ixv1, nod0[0]->idx);
	  dsp0[1] = (double *)WlzIndexedValueGet(ixv1, nod0[1]->idx);
	  dsp0[2] = (double *)WlzIndexedValueGet(ixv1, nod0[2]->idx);
	  dsp0[3] = (double *)WlzIndexedValueGet(ixv1, nod0[3]->idx);
	  v1.vtX = WlzGeomInterpolateTet3D(nod0[0]->pos, nod0[1]->pos,
					   nod0[2]->pos, nod0[3]->pos,
					   dsp0[0][0], dsp0[1][0],
					   dsp0[2][0], dsp0[3][0],
					   v0);
	  v1.vtY = WlzGeomInterpolateTet3D(nod0[0]->pos, nod0[1]->pos,
					   nod0[2]->pos, nod0[3]->pos,
					   dsp0[0][1], dsp0[1][1],
					   dsp0[2][1], dsp0[3][1],
					   v0);
	  v1.vtZ = WlzGeomInterpolateTet3D(nod0[0]->pos, nod0[1]->pos,
					   nod0[2]->pos, nod0[3]->pos,
					   dsp0[0][2], dsp0[1][2],
					   dsp0[2][2], dsp0[3][2],
					   v0);
	}
	else /* idE0 <0, the vertex is not in the mesh. Find the closest
	      * node to the vertex in the mesh and then use kriging. */
	{
	  int		idN0;
          WlzCMeshNod3D	*nod0,
	  		*nodN;

	  idN0 = WlzCMeshClosestNod3D(mesh1, v0);
          nod0 = (WlzCMeshNod3D *)AlcVectorItemGet(mesh1->res.nod.vec, idN0);
	  nNbr1 = WlzCMeshNodRingNodIndices3D(nod0, &maxNbrIdxBuf,
	  			              &nbrIdxBuf, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    /* Reallocate buffers if required. */
	    errNum = WlzKrigReallocBuffers3D(&nbrPosBuf, &posSV, &wSp,
					     &modelSV, &maxKrigBuf,
					     nNbr1, nNbr0);
	    nNbr0 = nNbr1;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    for(i = 0; i < nNbr1; ++i)
	    {
	      WlzCMeshNod3D *nod;
	      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh1->res.nod.vec,
						      nbrIdxBuf[i]);
	      nbrPosBuf[i] = nod->pos;
	    }
	    WlzKrigSetModelFn(&modelFn, WLZ_KRIG_MODELFN_LINEAR,
			      0.0, 0.1, 2.0 * dRange);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetModelSV3D(modelSV, &modelFn, nNbr1, nbrPosBuf,
					  wSp);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzKrigOSetPosSV3D(posSV, &modelFn, nNbr1, nbrPosBuf, v0);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	i;

	    v1.vtX = 0.0;
	    v1.vtY = 0.0;
	    v1.vtZ = 0.0;
	    WlzKrigOWeightsSolve(modelSV, posSV, wSp, WLZ_MESH_TOLERANCE);
	    /* posSV now contains the weights. */
	    for(i = 0; i < nNbr1; ++i) 
	    {
	      double *dsp0;

	      nodN = (WlzCMeshNod3D *)AlcVectorItemGet(mesh1->res.nod.vec,
	                                               nbrIdxBuf[i]);
	      dsp0 = (double *)WlzIndexedValueGet(ixv1, nodN->idx);
	      v1.vtX += posSV[i] * dsp0[0];
	      v1.vtY += posSV[i] * dsp0[1];
	      v1.vtZ += posSV[i] * dsp0[2];
	    }
	  }
	}
	WLZ_VTX_3_ADD(v1, v1, v0);
	WLZ_VTX_3_SUB(v1, v1, nodR->pos);
	dspR = (double *)WlzIndexedValueGet(ixvR, nodR->idx);
	dspR[0] = v1.vtX;
	dspR[1] = v1.vtY;
	dspR[2] = v1.vtZ;
      }
    }
  }
  AlcFree(wSp);
  AlcFree(posSV);
  AlcFree(nodTab);
  AlcFree(nbrIdxBuf);
  AlcFree(nbrPosBuf);
  AlgMatrixFree(modelSV);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(trR);
}

/*!
* \return	Woolz object with the given conforming mesh and scalar
* 		expansion factor values attached the the mesh elements,
* 		or NULL on error.
* \ingroup	WlzTransform
* \brief	Compute the scalar expansion factors for the mesh elements
* 		of the given conforming mesh transform. The expansion
* 		factors are attached to the elements of the returned
* 		conforming mesh object.
* 		The expansion factor is defined to be the trace of the
* 		strain tensor but in many cases the maximum eigenvalue
* 		is a more sensitive feature.
* \param	cObj			Given conforming mesh transform.
* \param	inverse			Use inverse of transform if non zero.
* \param	method			Method used:
* 					* 0 - average expansion from strain
* 					      tensor trace
* 					* 1 - maximum expansion from maximum
* 					      eigenvalue of the strain tensor
* 					* 2 - ratio of element volumes,
* 					expansion will be DBL_MAX for
* 					infinite expansion and zero for
* 					complete collapse.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshExpansion(WlzObject *cObj,
				   int inverse, int method,
				   WlzErrorNum *dstErr)
{
  WlzObject	*eObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	eObj = WlzCMeshExpansion2D(cObj, inverse, method, &errNum);
	break;
      case WLZ_CMESH_3D:
	eObj = WlzCMeshExpansion3D(cObj, inverse, method, &errNum);
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
  return(eObj);
}

/*!
* \return	Woolz object with the given conforming mesh and scalar
* 		expansion factor values attached the the mesh elements,
* 		or NULL on error.
* \ingroup	WlzTransform
* \brief	Compute the scalar expansion factors for the mesh elements
* 		of the given 2D conforming mesh transform. The expansion
* 		factors are attached to the elements of the returned
* 		conforming mesh object.
* 		See WlzCMeshExpansion().
* \param	cObj			Given conforming mesh transform.
* \param	inverse			Use inverse of transform if non zero.
* \param	method			Method used, see WlzCMeshExpansion().
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshExpansion2D(WlzObject *cObj,
				      int inverse, int method,
				      WlzErrorNum *dstErr)
{
  WlzObject	*eObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(eObj);
}

/*!
* \return	Woolz object with the given conforming mesh and scalar
* 		expansion factor values attached the the mesh elements,
* 		or NULL on error.
* \ingroup	WlzTransform
* \brief	Compute the scalar expansion factors for the mesh elements
* 		of the given 3D conforming mesh transform. The expansion
* 		factors are attached to the elements of the returned
* 		conforming mesh object.
* 		See WlzCMeshExpansion().
* \param	tObj			Given conforming mesh transform.
* \param	inverse			Use inverse of transform if non zero.
* \param	method			Method used, see WlzCMeshExpansion().
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshExpansion3D(WlzObject *cObj,
				      int inverse, int method,
				      WlzErrorNum *dstErr)
{
  int		nThr = 1,
  		useTensor = 0;
  WlzObject	*eObj = NULL,
  		*tObj = NULL;
  WlzCMesh3D	*mesh;
  AlgMatrixRect	**matTbl = NULL;
  WlzIndexedValues *ixvD,
  		   *ixvT,
  		   *ixvE = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  ixvD = cObj->values.x;
  if((method == 0) || (method == 1))
  {
    useTensor = 1;
    tObj = WlzCMeshStrainTensor(cObj, inverse, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(useTensor)
    {
      ixvT = tObj->values.x;
    }
    mesh = cObj->domain.cm3;
    ixvE = WlzMakeIndexedValues(cObj, 0, NULL, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_ELM, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzDomain	dom;
      WlzValues val;

      dom.cm3 = mesh;
      val.x = ixvE;
      eObj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreeIndexedValues(ixvE);
      ixvE = NULL;
    }
  }
  if((errNum == WLZ_ERR_NONE) && useTensor)
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
    if(method == 1) 				        /* Eigenvalue method */
    {
      if((matTbl = (AlgMatrixRect **)
                   AlcCalloc(nThr, sizeof(AlgMatrixRect *))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	int	idN;

	for(idN = 0; idN < nThr; ++idN)
	{
	  if((matTbl[idN] = AlgMatrixRectNew(3, 3, NULL)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE,
    		maxElm;
    AlcVector	*eVec;

    eVec = mesh->res.elm.vec;
    maxElm = mesh->res.elm.maxEnt;
#ifdef _OPENMP
#pragma omp parallel for private(idE) num_threads(nThr)
#endif
    for(idE = 0; idE < maxElm; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(eVec, idE);
      if(elm->idx >= 0)
      {
        double 	*fac,
		*ten;

	fac = (double *)WlzIndexedValueGet(ixvE, elm->idx);
	if(useTensor)
	{
          ten = (double *)WlzIndexedValueGet(ixvT, elm->idx);
	}
	switch(method)
	{
	  case 0:                                 /* Trace of strain tensor. */
	    *fac = ten[0] + ten[4] + ten[8];
	    break;
	  case 1:                                      /* Eigenvalue method. */
	    {
	      int 		idU,
			    idV,
			    idT,
			    thrId = 0;
	      double	val[3];
	      AlgMatrix 	mat;

	      idT = 0;
#ifdef _OPENMP
	      thrId = omp_get_thread_num();
#endif
	      mat.rect = matTbl[thrId];
	      val[0] = val[1] = val[2] = 0.0;
	      for(idV = 0; idV < 3; ++idV)
	      {
		for(idU = 0; idU < 3; ++idU)
		{
		  mat.rect->array[idV][idU] = ten[idT++];
		}
	      }
	      (void )AlgMatrixRSEigen(mat, val, 0);
	      *fac = ALG_MAX3(val[0], val[1], val[2]);
	    }
	    break;
	  case 2:
	    {
	      int	  idN;
	      double	  sVol,
	      		  dVol;
	      WlzCMeshNod3D *nod[4];
	      WlzDVertex3 sPos[4],
	      		  dPos[4];
              const double eps = 0.000001;

	      nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
	      nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
	      nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	      nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
	      for(idN = 0; idN < 4; ++idN)
	      {
		double	*dsp;

		sPos[idN] = nod[idN]->pos;
		dsp = (double *)WlzIndexedValueGet(ixvD, nod[idN]->idx);
		dPos[idN].vtX = sPos[idN].vtX + dsp[0];
		dPos[idN].vtY = sPos[idN].vtY + dsp[1];
		dPos[idN].vtZ = sPos[idN].vtZ + dsp[2];
	      }
	      sVol = fabs(WlzGeomTetraSnVolume6(sPos[0], sPos[1],
	                                        sPos[2], sPos[3]));
	      dVol = fabs(WlzGeomTetraSnVolume6(dPos[0], dPos[1],
	                                        dPos[2], dPos[3]));
	      if(inverse)
	      {
	        *fac = (dVol < eps)? DBL_MAX: sVol / dVol;
	      }
	      else
	      {
	        *fac = (sVol < eps)? DBL_MAX: dVol / sVol;
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_PARAM_DATA;
	    break;
	}
      }
    }
  }
  if(matTbl)
  {
    int		idN;

    for(idN = 0; idN < nThr; ++idN)
    {
      AlgMatrixRectFree(matTbl[idN]);
    }
  }
  (void )WlzFreeObj(tObj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(eObj);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(eObj);
}

