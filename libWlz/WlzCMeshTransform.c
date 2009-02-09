#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshTransform_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzCMeshTransform.c
* \author       Bill Hill
* \date         October 2004
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
* \brief	Functions for creating and applying 2D and 3D conforming
* 		mesh transforms.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

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
* \brief	Conforming mesh scanning element for 3D mesh.
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
  WlzCMeshTransform *mTr;		/*! The conforming mesh transform. */
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
  WlzCMeshTransform *mTr;		/*! The conforming mesh transform. */
  int		nItvs;			/*! Number of element intervals. */
  WlzCMeshScanItv3D *itvs;		/*! Element intervals sorted by line
  					    then left column. */
  WlzCMeshScanElm3D *dElm;		/*! Destination mesh element data. */
  WlzIBox3	dBox;			/*! Bounding box of the displaced
  					    mesh. */
} WlzCMeshScanWSp3D;

static void 			WlzCMeshUpdateScanElm2D(
				  WlzCMeshTransform *mTr,
				  WlzCMeshScanElm2D *sElm,
				  int fwd);
static void 			WlzCMeshUpdateScanElm3D(
				  WlzCMeshTransform *mTr,
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
static int			WlzCMeshScanTriElm2D(
				  WlzCMeshScanWSp2D *mSWSp,
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
static WlzErrorNum 		WlzCMeshTransMakeDispCb2D(
				  void *meshP,
				  void *nodP, 
				  void *mTrP);
static WlzErrorNum 		WlzCMeshTransMakeDispCb3D(
				  void *meshP,
				  void *nodP,
				  void *mTrP);
static WlzErrorNum 		WlzCMeshTransMakeDisp2D(
				  WlzCMeshTransform *mTr,
				  WlzCMesh2D *mesh,
				  WlzCMeshNod2D *nod,
				  int idx);
static WlzErrorNum 		WlzCMeshTransMakeDisp3D(
				  WlzCMeshTransform *mTr,
				  WlzCMesh3D *mesh,
				  WlzCMeshNod3D *nod,
				  int idx);
static WlzErrorNum 		WlzCMeshTransformValues2D(
				  WlzObject *dstObj,
				  WlzObject *srcObj,
				  WlzCMeshTransform *mTr,
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
static WlzObject 		*WlzCMeshTransformObjPDomain3D(
				  WlzObject *srcObj,
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshTransformObjV3D(
				  WlzObject *srcObj,
				  WlzCMeshTransform *mTr,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshScanObjPDomain3D(
				  WlzObject *srcObj,
				  WlzCMeshScanWSp3D *mSWSp,
				  WlzErrorNum *dstErr);
static WlzPolygonDomain 	*WlzCMeshTransformPoly(
				  WlzPolygonDomain *srcPoly,
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzBoundList 		*WlzCMeshTransformBoundList(
				  WlzBoundList *srcBound,
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp2D 	*WlzCMeshScanWSpInit2D(
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp3D 	*WlzCMeshMakeScanWSp3D(
				  WlzCMeshTransform *mTr,
				  int nItv,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp3D 	*WlzCMeshScanWSpInit3D(
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
#ifdef WLZ_CMESHTRANSFORM_DEBUG
static WlzErrorNum 		WlzCMeshVerifyWSp3D(
				  WlzObject *srcObj,
				  WlzCMeshScanWSp3D *mSWSp);
#endif

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Free's the given conforming mesh transform.
* \param	mTr			Given conforming mesh transform.
*/
WlzErrorNum	WlzFreeCMeshTransform(WlzCMeshTransform *mTr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mTr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mTr->type)
    {
      case WLZ_TRANSFORM_2D_CMESH: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_CMESH:
        (void )AlcVectorFree(mTr->dspVec);
	errNum = WlzCMeshFree(mTr->mesh);
        break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      AlcFree(mTr);
    }
  }
  return(errNum);
}

/*!
* \return	New 2D conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a new conforming mesh transform.
* \param	type			Required conforming mesh type.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzMakeCMeshTransform(WlzTransformType type,
					 WlzErrorNum *dstErr)
{
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(type)
  {
    case WLZ_TRANSFORM_2D_CMESH: /* FALLTHROUGH */
    case WLZ_TRANSFORM_3D_CMESH:
      if((mTr = (WlzCMeshTransform *)
      	        AlcCalloc(1, sizeof(WlzCMeshTransform))) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mTr->type = type;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mTr);
}

/*!
* \return	New 2D conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a new 2D conforming mesh transform from the
*		given 2D conforming mesh.
* \param	mesh			Given 2D conforming mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzMakeCMeshTransform2D(WlzCMesh2D *mesh,
					WlzErrorNum *dstErr)
{
  unsigned int	idN;
  WlzDomain	dom;
  WlzCMeshNod2D	*nod;
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_2D_CMESH, &errNum);
  }
  /* Create vector for the displacements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((mTr->dspVec = AlcVectorNew(1, sizeof(WlzDVertex2),
    				   mesh->res.nod.vec->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Assign the mesh to the transform. */
    dom.cm2 = mesh;
    (void )WlzAssignDomain(dom, NULL);
    mTr->mesh.m2 = mesh;
    /* Set the mesh node properties to be displacements for all existing
     * mesh nodes. */
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < mesh->res.nod.maxEnt))
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, (size_t )idN);
      errNum = WlzCMeshTransMakeDisp2D(mTr, mesh, nod, idN);
      ++idN;
    }
  }
  /* Set a mesh new node callback to set the node property to be a
   * displacement for all new nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshAddNewNodCb2D(mesh, WlzCMeshTransMakeDispCb2D,
    				   mTr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mTr);
}

/*!
* \return	New 3D conforming mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a new 3D conforming mesh transform from the
*		given 3D conforming mesh.
* \param	mesh			Given 3D conforming mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzMakeCMeshTransform3D(WlzCMesh3D *mesh,
					WlzErrorNum *dstErr)
{
  unsigned int	idN;
  WlzDomain	dom;
  WlzCMeshNod3D	*nod;
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_3D_CMESH, &errNum);
  }
  /* Create vector for the displacements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((mTr->dspVec = AlcVectorNew(1, sizeof(WlzDVertex3),
    				   mesh->res.nod.vec->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Assign the mesh to the transform. */
    dom.cm3 = mesh;
    (void )WlzAssignDomain(dom, NULL);
    mTr->mesh.m3 = mesh;
    /* Set the mesh node properties to be displacements for all existing
     * mesh nodes. */
    idN = 0;
    while((errNum == WLZ_ERR_NONE) && (idN < mesh->res.nod.maxEnt))
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, (size_t )idN);
      errNum = WlzCMeshTransMakeDisp3D(mTr, mesh, nod, idN);
      ++idN;
    }
  }
  /* Set a mesh new node callback to set the node property to be a
   * displacement for all new nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshAddNewNodCb3D(mesh, WlzCMeshTransMakeDispCb3D,
    				   mTr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mTr);
}

/*!
* \return       New conforming mesh transform.
* \ingroup      WlzTransform
* \brief        Creates a conforming mesh transform for the given object
*		with all mesh displacements zero.
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
WlzCMeshTransform *WlzCMeshTransformFromObj(WlzObject *srcObj,
				 WlzMeshGenMethod method,
                                 double minDist, double maxDist,
				 WlzObject **dstDilObj, int delOut,
                                 WlzErrorNum *dstErr)
{
  WlzCMeshP	mesh;
  WlzCMeshTransform *mTr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  mesh.v = NULL;
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
	mesh.m2 = WlzCMeshFromObj2D(srcObj, minDist, maxDist, dstDilObj,
				    delOut, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          mTr = WlzMakeCMeshTransform2D(mesh.m2, &errNum);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	mesh.m3 = WlzCMeshFromObj3D(srcObj, minDist, maxDist, dstDilObj,
				    delOut, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          mTr = WlzMakeCMeshTransform3D(mesh.m3, &errNum);
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
  return(mTr);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given integer vertex
*		array in place and using the given conforming mesh
*		transform. If a vertex is outside the mest it is
*		displaced using the displacement of the closest node
*		in the mesh.
* \param	mTr			The mesh transform.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2I(WlzCMeshTransform *mTr,
					 int nVtx, WlzIVertex2 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  WlzDVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  AlcVector     *vec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  vec = mTr->dspVec;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
			(double )(vtx[idN].vtX), (double )(vtx[idN].vtY),
		        &nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) ||
         ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      tVtx = *(WlzDVertex2 *)AlcVectorItemGet(vec, nearNod);
      WLZ_VTX_2_ADD(tVtx, tVtx, vtx[idN]);
    }
    vtx[idN].vtX = WLZ_NINT(tVtx.vtX);
    vtx[idN].vtY = WLZ_NINT(tVtx.vtY);
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
* \param	mTr			The mesh transform.
* 					vertices in the array after
* 					transformation. Must not be NULL.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2F(WlzCMeshTransform *mTr,
					 int nVtx, WlzFVertex2 *vtx)
{
  int		idN,
  		lastElmIdx,
		nearNod;
  WlzDVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  AlcVector     *vec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  vec = mTr->dspVec;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY,
			&nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) ||
         ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      tVtx = *(WlzDVertex2 *)AlcVectorItemGet(vec, nearNod);
      WLZ_VTX_2_ADD(tVtx, tVtx, vtx[idN]);
    }
    vtx[idN].vtX = tVtx.vtX;
    vtx[idN].vtY = tVtx.vtY;
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
* \param	mTr			The mesh transform.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2D(WlzCMeshTransform *mTr,
					 int nVtx,
					 WlzDVertex2 *vtx)
{
  int		idN,
		nearNod,
  		lastElmIdx;
  WlzDVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  AlcVector     *vec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  nearNod = -1;
  lastElmIdx = -1;
  vec = mTr->dspVec;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if(((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
    			vtx[idN].vtX, vtx[idN].vtY,
			&nearNod)) < 0) && (nearNod < 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if(sE.idx > 0)
    {
      if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
      {
	WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
	lastElmIdx = sE.idx;
      }
      tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
		 (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
      tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
		 (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    }
    else
    {
      tVtx = *(WlzDVertex2 *)AlcVectorItemGet(vec, nearNod);
      WLZ_VTX_2_ADD(tVtx, tVtx, vtx[idN]);
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
*               domain object is covered by the given mesh.
* \param        mesh                    Given mesh.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToDomObj(WlzCMeshP mesh, WlzErrorNum *dstErr)
{
  WlzObject     *domObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        domObj = WlzCMeshToDomObj2D(mesh.m2, &errNum);
        break;
      case WLZ_CMESH_TET3D:
        domObj = WlzCMeshToDomObj3D(mesh.m3, &errNum);
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
  return(domObj);
}

/*!
* \return       New 2D domain object, empty object if the mesh has no elements
*               or NULL on error.
* \ingroup      WlzMesh
* \brief        Computes a new 2D domain object, the domain of which
*               corresponds to the region of space enclosed by the 2D mesh,
*               such that the domain object is covered by the given mesh.
* \param        mesh                    Given 2D mesh.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToDomObj2D(WlzCMesh2D *mesh, WlzErrorNum *dstErr)
{
  int           idI,
                line,
                itvLnCnt,
                itvLnWidth,
                itvLnByteWidth;
  WlzIBox2      iBox;
  WlzCMeshScanWSp2D *mSWSp = NULL;
  WlzCMeshTransform *mTr = NULL;
  WlzCMeshScanItv2D *curItv,
  		*prvItv;
  WlzUByte	*lnMsk = NULL;
  WlzObject     *domObj = NULL;
  WlzDomain     dom;
  WlzValues     val;
  WlzDynItvPool itvPool;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(mesh->res.elm.numEnt < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if(mesh->res.elm.numEnt == 0)
  {
    domObj = WlzMakeEmpty(&errNum);
  }
  else
  {
    /* Make and initialise a constrained mesh scan workspace, but with
     * NULL displacements. */
    mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_2D_CMESH, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      mTr->mesh.m2 = mesh;
      mSWSp = WlzCMeshScanWSpInit2D(mTr, &errNum);
    }
    /* Create a single line bit mask buffer. */
    if(errNum == WLZ_ERR_NONE)
    {
      iBox.xMin = (int )floor(mesh->bBox.xMin);
      iBox.xMax = (int )ceil(mesh->bBox.xMax);
      iBox.yMin = (int )floor(mesh->bBox.yMin);
      iBox.yMax = (int )ceil(mesh->bBox.yMax);
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
      domObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
    }
  }
  AlcFree(lnMsk);
  if(mTr)
  {
    mTr->mesh.v = NULL;
    (void )WlzFreeCMeshTransform(mTr);
  }
  WlzCMeshScanWSpFree2D(mSWSp);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeDomain(dom);
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(domObj);
}

/*!
* \return       New 3D domain object, empty object if the mesh has no elements
*               or NULL on error.
* \ingroup      WlzMesh
* \brief        Computes a new 3D domain object, the domain of which
*               corresponds to the region of space enclosed by the 3D mesh,
*               such that the domain object is covered by the given mesh.
* \param        mesh                    Given 3D mesh.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToDomObj3D(WlzCMesh3D *mesh, WlzErrorNum *dstErr)
{
  WlzObject     *domObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(mesh->res.elm.numEnt < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    if(mesh->res.elm.numEnt == 0)
    {
      domObj = WlzMakeEmpty(&errNum);
    }
    else
    {
      errNum = WLZ_ERR_UNIMPLEMENTED; /* TODO */
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(domObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the transform coefficients for the given conforming
*		2D mesh scan element which must have a valid element index.
* \param	mTr			Conforming mesh transform.
* \param	sE			Mesh scan element.
* \param	fwd			Non-zero for forward transform
*					mapping source to destination,
* 					otherwise inverse transform.
*/
static void 	WlzCMeshUpdateScanElm2D(WlzCMeshTransform *mTr,
				        WlzCMeshScanElm2D *sE,
					int fwd)
{
  int		idN,
  		squash;
  double	areaSn2;
  WlzDVertex2	dVx[3],
  		sVx[3];
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D	*elm;
  AlcVector	*vec;

  sE->flags = WLZ_CMESH_SCANELM_NONE;
  vec = mTr->mesh.m2->res.elm.vec;
  elm = (WlzCMeshElm2D *)AlcVectorItemGet(vec, sE->idx);
  vec = mTr->dspVec;
  for(idN = 0; idN < 3; ++idN)
  {
    nod = elm->edu[idN].nod;
    sVx[idN] = nod->pos;
    dVx[idN] = *(WlzDVertex2 *)AlcVectorItemGet(vec, nod->idx);
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

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes the transform coefficients for the given conforming
*		3D mesh scan element which must have a valid element index.
* \param	mTr			Conforming mesh transform.
* \param	sE			Mesh scan element.
* \param	fwd			Non-zero for forward transform
*					mapping source to destination,
* 					otherwise inverse transform.
*/
static void 	WlzCMeshUpdateScanElm3D(WlzCMeshTransform *mTr,
				        WlzCMeshScanElm3D *sE,
					int fwd)
{
  int		idN,
  		squash;
  WlzDVertex3	dVx[4],
  		sVx[4];
  WlzCMeshElm3D	*elm;
  WlzCMeshNod3D	*nod[4];
  AlcVector	*vec;

  sE->flags = WLZ_CMESH_SCANELM_NONE;
  vec = mTr->mesh.m3->res.elm.vec;
  elm = (WlzCMeshElm3D *)AlcVectorItemGet(vec, sE->idx);
  WlzCMeshElmGetNodes3D(elm, nod + 0, nod + 1, nod + 2, nod + 3);
  vec = mTr->dspVec;
  for(idN = 0; idN < 4; ++idN)
  {
    sVx[idN] = nod[idN]->pos;
    dVx[idN] = *(WlzDVertex3 *)AlcVectorItemGet(vec, nod[idN]->idx);
    WLZ_VTX_3_ADD(dVx[idN], dVx[idN], sVx[idN]);
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

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Sets values in the destination object transforming those
*		of the source object.
* \param	dstObj			2D destination object with transformed
* 					domain but no values.
* \param	srcObj			2D source object.
* \param	mTr			Conforming mesh transform.
* \param	interp			Level of interpolation.
*/
static WlzErrorNum WlzCMeshTransformValues2D(WlzObject *dstObj,
					WlzObject *srcObj,
					WlzCMeshTransform *mTr,
					WlzInterpolationType interp)
{
  int		idP,
  		idX,
		iLft,
		iRgt,
		bufWidth,
		itvWidth,
		mItvIdx0,
		mItvIdx1;
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
  WlzGreyType	gType;
  WlzPixelV	bgdV;
  WlzDVertex2	sPosD;
  WlzCMeshScanItv2D *mItv0,
  		*mItv1,
		*mItv2;
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
    mSWSp = WlzCMeshScanWSpInit2D(mTr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItv0 = mSWSp->itvs;
    errNum = WlzInitGreyScan(dstObj, &iWSp, &gWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  }
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
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer.
*/
static WlzCMeshScanWSp2D *WlzCMeshScanWSpInit2D(WlzCMeshTransform *mTr,
				    	        WlzErrorNum *dstErr)
{
  int		iIdx;
  unsigned int 	eIdx;
  double	ndLn,
  		eLnMin,
		eLnMax;
  WlzDVertex2	*dsp;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D	*elm;
  WlzCMeshEntRes *elmRes;
  WlzCMeshScanElm2D *dElm;
  WlzCMeshScanWSp2D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmRes = &(mTr->mesh.m2->res.elm);
  if((mSWSp = (WlzCMeshScanWSp2D *)
  	      AlcCalloc(1, sizeof(WlzCMeshScanWSp2D))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp->mTr = mTr;
    /* Compute the maximum number of intervals in the displaced mesh. */
    eIdx = 0;
    for(eIdx = 0; eIdx < elmRes->maxEnt; ++eIdx)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(elmRes->vec, (size_t )eIdx);
      if(elm->idx >= 0)
      {
	nod = elm->edu[0].nod;
	if(mTr->dspVec == NULL)
	{
	  eLnMin = eLnMax = nod->pos.vtY;
	  nod = elm->edu[1].nod;
	  ndLn = nod->pos.vtY;
	}
	else
	{
	  dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	  eLnMin = eLnMax = nod->pos.vtY + dsp->vtY;
	  nod = elm->edu[1].nod;
	  dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	  ndLn = nod->pos.vtY + dsp->vtY;
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
	if(mTr->dspVec == NULL)
	{
	  ndLn = nod->pos.vtY;
	}
	else
	{
	  dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	  ndLn = nod->pos.vtY + dsp->vtY;
	}
	if(ndLn < eLnMin)
	{
	  eLnMin = ndLn;
	}
	else if(ndLn > eLnMax)
	{
	  eLnMax = ndLn;
	}
	mSWSp->nItvs += WLZ_NINT(eLnMax) - WLZ_NINT(eLnMin) + 1;
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
        iIdx += WlzCMeshScanTriElm2D(mSWSp, elm, iIdx);
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
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer.
*/
static WlzCMeshScanWSp3D *WlzCMeshScanWSpInit3D(WlzCMeshTransform *mTr,
				    	WlzErrorNum *dstErr)
{
  int		idE,
		idI,
		idN,
  		elmCnt,
		fstNod;
  WlzDVertex3	*dsp;
  WlzDVertex3	dspP;
  WlzDVertex3	dspPos[4];
  WlzDBox3	dBox;
  AlcVector	*itvVec = NULL;
  AlcVector	*elmVec;
  WlzCMeshNod3D	*nodBuf[4];
  WlzCMeshElm3D	*elm;
  WlzCMesh3D	*mesh;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzCMeshScanElm3D *dElm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  mesh = mTr->mesh.m3;
  elmVec = mesh->res.elm.vec;
  elmCnt = mesh->res.elm.maxEnt;
  /* Collect the intervals in the displaced mesh. */
  /* Create temporary vector in which to accumulate the intervals. */
  if((itvVec = AlcVectorNew(1, sizeof(WlzCMeshScanItv3D),
				 mesh->res.elm.vec->blkSz, NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
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
        WlzCMeshElmGetNodes3D(elm, nodBuf + 0, nodBuf + 1,
	                      nodBuf + 2, nodBuf + 3);
        for(idN = 0; idN < 4; ++idN)
	{
	  dsp = (WlzDVertex3 *)AlcVectorItemGet(mTr->dspVec, nodBuf[idN]->idx);
	  WLZ_VTX_3_ADD(dspP, nodBuf[idN]->pos, *dsp);
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
    mSWSp = WlzCMeshMakeScanWSp3D(mTr, idI, &errNum);
  }
  /* Copy the mesh scan intervals and sort them by plane, then line and
   * then column and then squeeze out the redundant intervals. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp->dBox.xMin = (int )floor(dBox.xMin) - 1;
    mSWSp->dBox.yMin = (int )floor(dBox.yMin) - 1;
    mSWSp->dBox.zMin = (int )floor(dBox.zMin) - 1;
    mSWSp->dBox.xMax = (int )ceil(dBox.xMax) + 1;
    mSWSp->dBox.yMax = (int )ceil(dBox.yMax) + 1;
    mSWSp->dBox.zMax = (int )ceil(dBox.zMax) + 1;
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
* \param	mTr			Conforming mesh transform.
* \param	nItv			Number of mesh scan intervals.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshScanWSp3D *WlzCMeshMakeScanWSp3D(WlzCMeshTransform *mTr,
						int nItv,
						WlzErrorNum *dstErr)
{
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mSWSp = (WlzCMeshScanWSp3D *)
  	      AlcCalloc(1, sizeof(WlzCMeshScanWSp3D))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    mSWSp->mTr = mTr;
    mSWSp->nItvs = nItv;
    if(((mSWSp->itvs = (WlzCMeshScanItv3D *)
		       AlcMalloc(sizeof(WlzCMeshScanItv3D) *
			         mSWSp->nItvs)) == NULL) ||
       ((mSWSp->dElm = (WlzCMeshScanElm3D *)
		       AlcCalloc(mTr->mesh.m3->res.elm.maxEnt,
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
* \param	elm			Mesh element which is assumed to be
*					valid.
* \param	iIdx			Conforming mesh element interval index.
*/
static int	WlzCMeshScanTriElm2D(WlzCMeshScanWSp2D *mSWSp,
				     WlzCMeshElm2D *elm, int iIdx)
{
  int		kolI0,
		kolI1,
  		lineI,
		ndIdx0,
  		ndIdx1,
		ndIdx2,
		iCnt = 0;
  double	x0,
  		x1,
		tD0,
  		tD1;
  double	inc[3];
  WlzDVertex2	dNd[3],
  		sNd[3];
  WlzCMeshNod2D	*nod;
  WlzDVertex2	dVx0,
  		dVx1;
  WlzCMeshScanItv2D *itv;

  /* Compute the integer displaced nodes of the element. */
  for(ndIdx0 = 0; ndIdx0 < 3; ++ndIdx0)
  {
    nod = elm->edu[ndIdx0].nod;
    dVx0 = nod->pos;
    if(mSWSp->mTr->dspVec == NULL)
    {
      tD0 = dVx0.vtX;
      tD1 = dVx0.vtY;
    }
    else
    {
      dVx1 = *(WlzDVertex2 *)AlcVectorItemGet(mSWSp->mTr->dspVec, nod->idx);
      tD0 = dVx0.vtX + dVx1.vtX;
      tD1 = dVx0.vtY + dVx1.vtY;
    }
    dNd[ndIdx0].vtX = tD0;
    dNd[ndIdx0].vtY = tD1;
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
  if(fabs(dNd[2].vtY) < 1.5)
  {
    /* Possible single horizontal interval
     * *-*-* */
    lineI = WLZ_NINT(sNd[0].vtY);
    if(fabs((double )lineI - sNd[0].vtY) <= DBL_EPSILON)
    {
      if(sNd[0].vtX <= sNd[1].vtX)
      {
        if(sNd[2].vtX <= sNd[0].vtX)
	{
	  tD0 = sNd[2].vtX;
	  tD1 = WLZ_MAX(sNd[0].vtX, sNd[1].vtX);
	}
	else /* sNd[2].vtX > Nd[0].vtX */
	{
	  tD0 = sNd[0].vtX;
	  tD1 = WLZ_MAX(sNd[1].vtX, sNd[2].vtX);
	}
      }
      else /* sNd[0].vtX > sNd[1].vtX */
      {
        if(sNd[2].vtX <= sNd[1].vtX)
	{
	  tD0 = sNd[2].vtX;
	  tD1 = WLZ_MAX(sNd[0].vtX, sNd[1].vtX);
	}
	else
	{
	  tD0 = sNd[1].vtX;
	  tD1 = WLZ_MAX(sNd[2].vtX, sNd[0].vtX);
	}
      }
      kolI0 = (int )floor(tD0 + DBL_EPSILON);
      if(kolI0 < tD1)
      {
	itv->elmIdx = elm->idx;
	itv->line = lineI;
	itv->lftI = itv->rgtI = kolI0;
	itv->rgtI = (int )floor(tD1 + DBL_EPSILON);
	iCnt = 1;
      }
    }
  }
  else if((fabs(dNd[0].vtX) < 1.5) &&
          (fabs(dNd[1].vtX) < 1.5))
  {
    /* Many possible single column intervals
     * 0
     * |
     * 1
     * |
     * 2 */
    kolI0 = WLZ_NINT(sNd[0].vtX);
    if(fabs(kolI0 - sNd[0].vtX) <= DBL_EPSILON)
    {
      lineI = (int )ceil(sNd[0].vtY - DBL_EPSILON);
      while(lineI < sNd[2].vtY)
      {
	itv->elmIdx = elm->idx;
	itv->line = lineI;
	itv->lftI = itv->rgtI = kolI0;
	++itv;
	++lineI;
	++iCnt;
      }
    }

  }
  else
  {
    /* General case for triangles with non-zero area. */
    itv = mSWSp->itvs + iIdx;
    lineI = (int )ceil(sNd[0].vtY - DBL_EPSILON);
    inc[0] = (fabs(dNd[0].vtY) > DBL_EPSILON)? dNd[0].vtX / dNd[0].vtY: 0.0;
    inc[1] = (fabs(dNd[1].vtY) > DBL_EPSILON)? dNd[1].vtX / dNd[1].vtY: 0.0;
    inc[2] = (fabs(dNd[2].vtY) > DBL_EPSILON)? dNd[2].vtX / dNd[2].vtY: 0.0;
    while(lineI < sNd[2].vtY)
    {
      x0 = sNd[0].vtX + inc[2] * (lineI - sNd[0].vtY);
      x1 = (lineI >= sNd[1].vtY)? sNd[1].vtX + inc[1] * (lineI - sNd[1].vtY):
                                  sNd[0].vtX + inc[0] * (lineI - sNd[0].vtY);
      if(x0 > x1)
      {
        tD0 = x0; x0 = x1; x1 = tD0;
      }
      kolI0 = floor(x0 + DBL_EPSILON);
      kolI1 = floor(x1 + DBL_EPSILON);
      if((kolI0 < x1) && (kolI1 > x0))
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
  pl = ceil(vtx[0].vtZ);
  if(pl < vtx[3].vtZ + DBL_EPSILON)
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
	     ((a = (pl - vtx[0].vtZ) / del10.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
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
	     ((a = (pl - vtx[0].vtZ) / del20.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
	  {
	    isn[isnCnt].vtX = vtx[0].vtX + (a * del20.vtX);
	    isn[isnCnt].vtY = vtx[0].vtY + (a * del20.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del21.vtZ > tol) &&
	     ((a = (pl - vtx[1].vtZ) / del21.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
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
	     ((a = (pl - vtx[0].vtZ) / del30.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
	  {
	    isn[isnCnt].vtX = vtx[0].vtX + (a * del30.vtX);
	    isn[isnCnt].vtY = vtx[0].vtY + (a * del30.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del31.vtZ > tol) &&
	     ((a = (pl - vtx[1].vtZ) / del31.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
	  {
	    isn[isnCnt].vtX = vtx[1].vtX + (a * del31.vtX);
	    isn[isnCnt].vtY = vtx[1].vtY + (a * del31.vtY);
	    isn[isnCnt++].vtZ = pl;
	  }
	  if((del32.vtZ > tol) &&
	     ((a = (pl - vtx[2].vtZ) / del32.vtZ) > DBL_EPSILON) &&
	     (a < 1.0 - DBL_EPSILON))
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
    } while((errNum == WLZ_ERR_NONE) && (pl < vtx[3].vtZ + DBL_EPSILON));
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
  sVtx[0].vtX = WLZ_NINT(vtx[idV0].vtX);
  sVtx[0].vtY = WLZ_NINT(vtx[idV0].vtY);
  sVtx[1].vtX = WLZ_NINT(vtx[idV1].vtX);
  sVtx[1].vtY = WLZ_NINT(vtx[idV1].vtY);
  sVtx[2].vtX = WLZ_NINT(vtx[idV2].vtX);
  sVtx[2].vtY = WLZ_NINT(vtx[idV2].vtY);
  /* Compute deltas. */
  dX[0] = sVtx[0].vtX - sVtx[1].vtX;
  dX[1] = sVtx[1].vtX - sVtx[2].vtX;
  dY[2] = sVtx[2].vtY - sVtx[0].vtY;
  /* If the triangles vertices are not coincident scan convert it. */
  if(dY[2] && (dX[0] || dX[1]))
  {
    plI = WLZ_NINT(vtx[0].vtZ);
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
	itv->lftI = itv->rgtI = WLZ_NINT(klD);
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
	  klI = WLZ_NINT(klD);
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
	  klI = WLZ_NINT(klD);
	  itv = (WlzCMeshScanItv3D *)AlcVectorItemGet(itvVec, idI0);
	  if(klI > itv->lftI)
	  {
	    itv->rgtI = klI;
	  }
	  else
	  {
	    itv->lftI = klI;
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
    itv->line = (int )floor(vtx[0].vtY);
    itv->plane = (int )floor(vtx[0].vtZ);
    if(vtx[0].vtX < vtx[1].vtX)
    {
      itv->lftI = (int )floor(vtx[0].vtX);
      itv->rgtI = (int )floor(vtx[1].vtX);
    }
    else
    {
      itv->lftI = (int )floor(vtx[1].vtX);
      itv->rgtI = (int )floor(vtx[0].vtX);
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
* \param	mTr			Conforming mesh transform.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzCMeshTransformObj(WlzObject *srcObj,
				     WlzCMeshTransform *mTr,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzDomain	dstDom;
  WlzValues	srcValues,
  		dstValues;
  WlzObject	*tObj0 = NULL,
  		*tObj1 = NULL,
		*dstObj = NULL;
  WlzPixelV	bgdV;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
  dstDom.core = NULL;
  srcValues.core = NULL;
  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mTr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mTr->type)
    {
      case WLZ_TRANSFORM_2D_CMESH: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_CMESH:
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
						      mTr, &errNum);
	      break;
	    case WLZ_BOUNDLIST:
	      dstDom.b = WlzCMeshTransformBoundList(srcObj->domain.b,
						    mTr, &errNum);
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
	  tObj1 = NULL;
	  tObj0 = WlzObjToBoundary(srcObj, 1, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj1 = WlzCMeshTransformObj(tObj0, mTr, interp, &errNum);
	  }
	}
	(void )WlzFreeObj(tObj0); tObj0 = NULL;
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzBoundToObj(tObj1->domain.b,
				 WLZ_SIMPLE_FILL, &errNum);
	}
	(void )WlzFreeObj(tObj1); tObj1 = NULL;
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
	    errNum = WlzCMeshTransformValues2D(dstObj, srcObj, mTr, interp);
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
	  dstObj = WlzCMeshTransformObjPDomain3D(srcObj, mTr, &errNum);
	}
	else
	{
	  dstObj = WlzCMeshTransformObjV3D(srcObj, mTr, interp, &errNum);
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
* \return	Transformed boundary list or NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms the given boundary list using the given conforming
*		mesh transform.
* \param	srcBound		Given boundary list.
* \param	mesh			Mesh transform to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzBoundList *WlzCMeshTransformBoundList(WlzBoundList *srcBound,
					WlzCMeshTransform *mTr,
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
				        mTr, &errNum)) != NULL)
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
    if((dumDom.b = WlzCMeshTransformBoundList(srcBound->next, mTr,
					      &errNum)) != NULL)
    {
      (void )WlzAssignDomain(dumDom, &errNum);
      dstBnd->next = dumDom.b;
    }
  }
  /* Transform down boundlist. */
  if((errNum == WLZ_ERR_NONE) && (srcBound->down != NULL))
  {
    if((dumDom.b = WlzCMeshTransformBoundList(srcBound->down, mTr,
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
* \param	mesh			Mesh transform to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPolygonDomain *WlzCMeshTransformPoly(WlzPolygonDomain *srcPoly,
						  WlzCMeshTransform *mTr,
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
        errNum = WlzCMeshTransformVtxAry2I(mTr, dstPoly->nvertices,
					   dstPoly->vtx);
        break;
      case WLZ_POLYGON_FLOAT:
        errNum = WlzCMeshTransformVtxAry2F(mTr, dstPoly->nvertices,
					   (WlzFVertex2 *)(dstPoly->vtx));
        break;
      case WLZ_POLYGON_DOUBLE:
        errNum = WlzCMeshTransformVtxAry2D(mTr, dstPoly->nvertices,
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
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformObjPDomain3D(WlzObject *srcObj,
				     WlzCMeshTransform *mTr,
				     WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzCMeshScanWSp3D *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Make workspace intervals for the elements in the displaced
   * mesh, with intervals sorted by plane, line and then column. */
  if(errNum == WLZ_ERR_NONE)
  {
    mSWSp = WlzCMeshScanWSpInit3D(mTr, &errNum);
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
* \brief	Applies a 3D conforming mesh transform to the given source
*		object which must be a 3D domain object with values.
* \param	srcObj			Object to be transformed.
* \param	mTr			Conforming mesh transform.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformObjV3D(WlzObject *srcObj,
				     WlzCMeshTransform *mTr,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzGreyType	gType;
  WlzPixelV	bgdV;
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
    mSWSp = WlzCMeshScanWSpInit3D(mTr, &errNum);
  }
  /* Scan through the sorted intervals creating domains as required. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzCMeshScanObjPDomain3D(srcObj, mSWSp, &errNum); 
  }
  /* Make a voxel value table for the new domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    dstValues.vox = WlzNewValuesVox(dstObj, gType, bgdV, &errNum);
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
      sPos.vtX = WLZ_NINT(tV.vtX);
      sPos.vtY = WLZ_NINT(tV.vtY);
      sPos.vtZ = WLZ_NINT(tV.vtZ);
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
					      NULL);
	      /* Output message. */
	      (void )fprintf(stderr, "WlzCMeshVerifyWSp3D() %d %d %d n %d\n",
	                     sPos.vtX, sPos.vtY, sPos.vtZ, idN - 1);
	    }
	    ++inp;
	  }
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
* \param	srcObj			Object to be transformed.
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshScanObjPDomain3D(WlzObject *srcObj,
					WlzCMeshScanWSp3D *mSWSp,
					WlzErrorNum *dstErr) 
{
  int		idI,
		kol,
		itvLnCnt,
		itvPlCnt,
		itvLnWidth,
		itvLnByteWidth;
  WlzIVertex3	dPos,
		sPos;
  WlzDVertex3	tV;
  WlzDynItvPool	itvPool;
  WlzCMeshScanItv3D *curItv,
  		*prvItv;
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
	sPos.vtX = WLZ_NINT(tV.vtX);
	sPos.vtY = WLZ_NINT(tV.vtY);
	sPos.vtZ = WLZ_NINT(tV.vtZ);
	if(WlzInsideDomain(srcObj, sPos.vtZ, sPos.vtY,
	                   sPos.vtX, NULL) != 0)
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
      dom3.p->voxel_size[0] = srcObj->domain.p->voxel_size[0];
      dom3.p->voxel_size[1] = srcObj->domain.p->voxel_size[1];
      dom3.p->voxel_size[2] = srcObj->domain.p->voxel_size[2];
      /* Create new object from the transformed plane domain. */
      dstObj = WlzMakeMain(srcObj->type, dom3, nullVal, NULL, NULL, &errNum);
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
  WlzCMeshScanElm3D *sElm;
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
	      sElm = mSWSp->dElm + mItv1->elmIdx;
	      tV.vtX = (sElm->tr[ 1] * dPos.vtY) + (sElm->tr[ 2] * dPos.vtZ) +
		       sElm->tr[ 3];
	      tV.vtY = (sElm->tr[ 5] * dPos.vtY) + (sElm->tr[ 6] * dPos.vtZ) +
		       sElm->tr[ 7];
	      tV.vtZ = (sElm->tr[ 9] * dPos.vtY) + (sElm->tr[10] * dPos.vtZ) +
		       sElm->tr[11];
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
			sPos.vtX = WLZ_NINT(sPosD.vtX);
			sPos.vtY = WLZ_NINT(sPosD.vtY);
			sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  *(olpBuf.dbp + idI) += WLZ_NINT(tD0);
			}
			else
			{
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  *(olpBuf.dbp + idI) += WLZ_NINT(tD0);
			}
			else
			{
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			sPosD.vtX = (sElm->tr[ 0] * dPos.vtX) + tV.vtX;
			sPosD.vtY = (sElm->tr[ 4] * dPos.vtX) + tV.vtY;
			sPosD.vtZ = (sElm->tr[ 8] * dPos.vtX) + tV.vtZ;
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
			  sPos.vtX = WLZ_NINT(sPosD.vtX);
			  sPos.vtY = WLZ_NINT(sPosD.vtY);
			  sPos.vtZ = WLZ_NINT(sPosD.vtZ);
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
			       *(olpBuf.inp + idX) / *(olpCnt + idX):
			       bgdV.v.shv;
	  }
	  break;
	case WLZ_GREY_UBYTE:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.ubp + idX) = (*(olpCnt + idX) >= 1)?
			       *(olpBuf.inp + idX) / *(olpCnt + idX):
			       bgdV.v.ubv;
	  }
	  break;
	case WLZ_GREY_FLOAT:
	  for(idX = 0; idX < cnt; ++idX)
	  {
	    *(dGP.flp + idX) = (*(olpCnt + idX) >= 1)?
			       *(olpBuf.dbp + idX) / *(olpCnt + idX):
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
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Callback function for a 2D mesh which creates a displacement
*		property for a node.
* \param	meshP			Used to pass the 2D mesh.
* \param	nodP			Used to pass the 2D node.
* \param	mTrP			Used to pass the mesh transform.
*/
static WlzErrorNum WlzCMeshTransMakeDispCb2D(void *meshP,
					void *nodP, void *mTrP)
{
  WlzCMesh2D	*mesh;
  WlzCMeshNod2D	*nod;
  WlzCMeshTransform *mTr;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  mesh = (WlzCMesh2D *)meshP;
  nod = (WlzCMeshNod2D *)nodP;
  mTr = (WlzCMeshTransform *)mTrP;
  if(mesh && nod && mTr && (nod->idx >= 0))
  {
    errNum = WlzCMeshTransMakeDisp2D(mTr, mesh, nod, nod->idx);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Callback function for a 3D mesh which creates a displacement
*		property for a node.
* \param	meshP			Used to pass the 3D mesh.
* \param	nodP			Used to pass the 3D node.
* \param	mTrP			Used to pass the mesh transform.
*/
static WlzErrorNum WlzCMeshTransMakeDispCb3D(void *meshP,
					void *nodP, void *mTrP)
{
  WlzCMesh3D	*mesh;
  WlzCMeshNod3D	*nod;
  WlzCMeshTransform *mTr;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  mesh = (WlzCMesh3D *)meshP;
  nod = (WlzCMeshNod3D *)nodP;
  mTr = (WlzCMeshTransform *)mTrP;
  if(mesh && nod && mTr && (nod->idx >= 0))
  {
    errNum = WlzCMeshTransMakeDisp3D(mTr, mesh, nod, nod->idx);
  }
  return(errNum);
}

/*!
* \return
* \ingroup      WlzTransform
* \brief	Creates a displacement for the given 2D conforming mesh node.
* \param	mesh			Given 2D conforming mesh.
* \param	vec			Vector from which to allocate the
* 	 				displacement.
* \param	nod			Node to have displacement.
* \param	idx			Index of the node in it's vector,
*					must be valid.
*/
static WlzErrorNum WlzCMeshTransMakeDisp2D(WlzCMeshTransform *mTr,
					   WlzCMesh2D *mesh,
					   WlzCMeshNod2D *nod, int idx)
{
  WlzDVertex2	*dsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dsp = (WlzDVertex2 *)AlcVectorExtendAndGet(mTr->dspVec,
  						 idx)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    nod->prop = (void *)dsp;
    dsp->vtX = dsp->vtY = 0.0;
  }
  return(errNum);
}

/*!
* \return
* \ingroup      WlzTransform
* \brief	Creates a displacement for the given 3D conforming mesh node.
* \param	mesh			Given 3D conforming mesh.
* \param	vec			Vector from which to allocate the
* 	 				displacement.
* \param	nod			Node to have displacement.
* \param	idx			Index of the node in it's vector,
*					must be valid.
*/
static WlzErrorNum WlzCMeshTransMakeDisp3D(WlzCMeshTransform *mTr,
					   WlzCMesh3D *mesh,
					   WlzCMeshNod3D *nod, int idx)
{
  WlzDVertex3	*dsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((dsp = (WlzDVertex3 *)AlcVectorExtendAndGet(mTr->dspVec,
  						 idx)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    nod->prop = (void *)dsp;
    dsp->vtX = dsp->vtY = dsp->vtZ = 0.0;
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
* \param        meshTr                  Given mesh transform.
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
WlzErrorNum     WlzCMeshGetNodesAndEdges(WlzCMeshTransform *meshTr,
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

  if(meshTr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstNNod == NULL) || (dstNod == NULL) ||
          (dstNDsp == NULL) || (dstDsp == NULL) ||
	  (dstNEdg == NULL) || (dstEdg == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    mesh = meshTr->mesh.m2;
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
