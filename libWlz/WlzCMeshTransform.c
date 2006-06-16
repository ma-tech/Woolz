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
					    negligable area and only the
					    transform coefficients for
					    translation are non zero. */
  WLZ_CMESH_SCANELM_FWD 	= (2)	/*! The transform is forward, ie
  					    transforming source to
					    destination. If not set then
					    the transform is inverse, ie
  					    transforming destination to
					    source. */
} WlzCMeshScanElmFlags;

/*!
* \struct	_WlzCMeshScanElm2D
* \ingroup	WlzTransform
* \brief	Conforming mesh scanning element.
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
* \struct       _WlzCMeshScanItv
* \ingroup      WlzTransform
* \brief        Scan interval within a 2D conforming mesh element.
*/
typedef struct _WlzCMeshScanItv
{
  int           elmIdx;                 /*! Element index. */
  int           line;                   /*! Line of interval. */
  int           lftI;                   /*! Start of interval. */
  int           rgtI;                   /*! End of interval. */
} WlzCMeshScanItv;

/*!
* \struct       _WlzCMeshScanWSp
* \ingroup      WlzTransform
* \brief        Conforming mesh scanning workspace.
*/
typedef struct _WlzCMeshScanWSp
{
  WlzCMeshTransform *mTr;		/*! The conforming mesh transform. */
  int		nItvs;			/*! Number of element intervals. */
  WlzCMeshScanItv *itvs;		/*! Element intervals sorted by line
  					    then left column. */
  WlzCMeshScanElm2D *dElm;		/*! Destination mesh element data. */
} WlzCMeshScanWSp;

static void 			WlzCMeshUpdateScanElm2D(
				  WlzCMeshTransform *mTr,
				  WlzCMeshScanElm2D *sE,
				  int fwd);
static void			WlzCMeshScanWSpFree(
				  WlzCMeshScanWSp *mSWSp);
static int			WlzCMeshScanTriElm(
				  WlzCMeshScanWSp *mSWSp,
				  WlzCMeshElm2D *elm,
				  int iIdx);
static int			WlzCMeshItvCmp(
				  const void *cmp0,
				  const void *cmp1);
static WlzErrorNum 		WlzCMeshTransMakeDispCb2D(
				  void *meshP,
				  void *nodP, 
				  void *mTrP);
static WlzErrorNum 		WlzCMeshTransMakeDisp2D(
				  WlzCMeshTransform *mTr,
				  WlzCMesh2D *mesh,
				  WlzCMeshNod2D *nod,
				  int idx);
static WlzErrorNum 		WlzCMeshTransformValues2D(
				  WlzObject *dstObj,
				  WlzObject *srcObj,
				  WlzCMeshTransform *mTr,
				  WlzInterpolationType interp);
static WlzObject 		*WlzCMeshTransformObj2D(
				  WlzObject *srcObj,
				  WlzCMeshTransform *mTr,
				  WlzInterpolationType interp,
				  WlzErrorNum *dstErr);
static WlzPolygonDomain 	*WlzCMeshTransformPoly(
				  WlzPolygonDomain *srcPoly,
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzBoundList 		*WlzCMeshTransformBoundList(
				  WlzBoundList *srcBound,
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);
static WlzCMeshScanWSp 		*WlzCMeshScanWSpInit(
				  WlzCMeshTransform *mTr,
				  WlzErrorNum *dstErr);

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
      case WLZ_TRANSFORM_2D_CMESH:
        (void )AlcVectorFree(mTr->dspVec);
	errNum = WlzCMeshFree2D(mTr->mesh.m2);
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
    case WLZ_TRANSFORM_2D_CMESH:
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
* \return       New conforming mesh transform.
* \ingroup      WlzTransform
* \brief        Creates a conforming mesh transform for the given object
*		with all mesh displacements zero.
* \param        srcObj                  The given object.
* \param        method                  Mesh generation method to use.
* \param        minDist                 Minimum distance between mesh vertices.
* \param        maxDist                 Maximum distance between mesh vertices.
* \param	dstDilObj		Destination pointer for the dilated
*					object used to build the mesh.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzCMeshTransform *WlzCMeshTransformFromObj(WlzObject *srcObj,
				 WlzMeshGenMethod method,
                                 double minDist, double maxDist,
				 WlzObject **dstDilObj,
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
  else
  {
    switch(srcObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	mesh.m2 = WlzCMeshFromObj2D(srcObj, minDist, maxDist, dstDilObj,
				    &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          mTr = WlzMakeCMeshTransform2D(mesh.m2, &errNum);
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
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a conforming mesh transform to the given source
*		object.
* \param	srcObj			Object to be transformed.
* \param	mTr			Conforming mesh transform.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshTransformObj(WlzObject *srcObj,
				     WlzCMeshTransform *mTr,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
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
      case WLZ_TRANSFORM_2D_CMESH:
	dstObj = WlzCMeshTransformObj2D(srcObj, mTr, interp, &errNum);
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
  return(dstObj);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given integer vertex
*		array in place and using the given conforming mesh
*		transform. It is an error if any vertex is not contained
*		within an element of the conforming mesh.
* \param	mTr			The mesh transform.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2I(WlzCMeshTransform *mTr,
					 int nVtx, WlzIVertex2 *vtx)
{
  int		idN,
  		lastElmIdx;
  WlzDVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  lastElmIdx = -1;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
		     (double )(vtx[idN].vtX), (double )(vtx[idN].vtY))) < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
    {
      WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
      lastElmIdx = sE.idx;
    }
    tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
	       (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
    tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
	       (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
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
*		transform. It is an error if any vertex is not contained
*		within an element of the conforming mesh.
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
  		lastElmIdx;
  WlzFVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  lastElmIdx = -1;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
    					   vtx[idN].vtX, vtx[idN].vtY)) < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
    {
      WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
      lastElmIdx = sE.idx;
    }
    tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
	       (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
    tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
	       (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    vtx[idN] = tVtx;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex
*		array in place and using the given conforming mesh
*		transform. It is an error if any vertex is not contained
*		within an element of the conforming mesh.
* \param	mTr			The mesh transform.
* \param	nVtx			Number of vertices in the array.
* \param	vtx			Array of vertices.
*/
WlzErrorNum	WlzCMeshTransformVtxAry2D(WlzCMeshTransform *mTr,
					 int nVtx,
					 WlzDVertex2 *vtx)
{
  int		idN,
  		lastElmIdx;
  WlzDVertex2	tVtx;
  WlzCMeshScanElm2D sE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  lastElmIdx = -1;
  for(idN = 0; idN < nVtx; ++idN)
  {
    if((sE.idx = WlzCMeshElmEnclosingPos2D(mTr->mesh.m2, lastElmIdx,
    					   vtx[idN].vtX, vtx[idN].vtY)) < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
      break;
    }
    if((sE.idx != lastElmIdx) || ((sE.flags & WLZ_CMESH_SCANELM_FWD) == 0))
    {
      WlzCMeshUpdateScanElm2D(mTr, &sE, 1);
      lastElmIdx = sE.idx;
    }
    tVtx.vtX = (sE.trX[0] * vtx[idN].vtX) +
               (sE.trX[1] * vtx[idN].vtY) + sE.trX[2];
    tVtx.vtY = (sE.trY[0] * vtx[idN].vtX) +
               (sE.trY[1] * vtx[idN].vtY) + sE.trY[2];
    vtx[idN] = tVtx;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computs the transform coefficients for the given conforming
*		mesh scan element which must have a valid element index.
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
    nod = elm->edg[idN].nod;
    sVx[idN] = nod->pos;
    dVx[idN] = *(WlzDVertex2 *)AlcVectorItemGet(vec, nod->idx);
    WLZ_VTX_2_ADD(dVx[idN], dVx[idN], sVx[idN]);
  }
  if(fwd)
  {
    sE->flags |= WLZ_CMESH_SCANELM_FWD;
    areaSn2 = WlzGeomTriangleSnArea2(sVx[0], sVx[1], sVx[2]);
    squash = WlzGeomTriangleAffineSolve(sE->trX, sE->trY, areaSn2,
    					sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
  }
  else
  {
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
* \brief	Creates a new value table for the destination object by
*		transforming those of the source object.
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
  WlzGreyP	olpBuf;
  int		tI[4];
  WlzGreyP	dGP;
  WlzGreyType	tGreyType;
  WlzPixelV	bkdV;
  WlzDVertex2	sPosD;
  WlzValues	newValues;
  WlzCMeshScanItv *mItv0,
  		*mItv1,
		*mItv2;
  WlzCMeshScanWSp *mSWSp = NULL;
  WlzCMeshScanElm2D *sElm;
  WlzGreyValueWSpace *srcGVWSp = NULL;
  WlzGreyWSpace dstGWSp;
  WlzIntervalWSpace dstIWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  olpBuf.inp = NULL;
  newValues.core = NULL;
  bkdV = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    tGreyType = WlzGreyTableTypeToGreyType(srcObj->values.v->type, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzValueConvertPixel(&bkdV, bkdV, tGreyType);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newValues.v = WlzNewValueTb(dstObj, srcObj->values.v->type,
    				bkdV, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj->values = WlzAssignValues(newValues, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bufWidth = dstObj->domain.i->lastkl - dstObj->domain.i->kol1 + 1;
    switch(tGreyType)
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
    }
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
    mItvIdx0 = 0;
    mSWSp = WlzCMeshScanWSpInit(mTr, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItv0 = mSWSp->itvs;
    errNum = WlzInitGreyScan(dstObj, &dstIWSp, &dstGWSp);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    srcGVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
  }
  while((errNum == WLZ_ERR_NONE) &&
	(WlzNextGreyInterval(&dstIWSp) == 0))
  {
    dGP = dstGWSp.u_grintptr;
    itvWidth = dstIWSp.rgtpos - dstIWSp.lftpos + 1;
    switch(tGreyType)
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
    }
    WlzValueSetInt(olpCnt, 0, itvWidth);
    /* Update the mesh interval pointer so that it points to the first
     * mesh interval on the which intersects the current grey interval. */
    while((mItv0->line < dstIWSp.linpos) && (mItvIdx0 < mSWSp->nItvs))
    {
      ++mItvIdx0;
      ++mItv0;
    }
    while((mItv0->line <= dstIWSp.linpos) &&
	  (mItv0->rgtI < dstIWSp.lftpos) &&
	  (mItvIdx0 < mSWSp->nItvs))
    {
      ++mItvIdx0;
      ++mItv0;
    }
    if((mItv0->line == dstIWSp.linpos) &&
       (dstIWSp.lftpos <= mItv0->rgtI) &&
       (dstIWSp.rgtpos >= mItv0->lftI))
    {
      /* Mesh interval mItv0 intersects the current grey interval find
       * the last mesh interval mItv1 which also intersects the current grey
       * interval. */
      mItv1 = mItv0;
      mItvIdx1 = mItvIdx0;
      while((mItv1->line == dstIWSp.linpos) &&
	    (mItv1->lftI <= dstIWSp.rgtpos) &&
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
	trXYC = (sElm->trX[1] * dstIWSp.linpos) + sElm->trX[2];
	trYX = sElm->trY[0];
	trYYC = (sElm->trY[1] * dstIWSp.linpos) + sElm->trY[2];
        /* Find length of intersection and set the grey pointer. */
	iLft = ALG_MAX(mItv1->lftI, dstIWSp.lftpos);
	iRgt = ALG_MIN(mItv1->rgtI, dstIWSp.rgtpos);
	idX = iLft;
	switch(interp)
	{
	  case WLZ_INTERPOLATION_NEAREST:
	    switch(dstGWSp.pixeltype)
	    {
	      case WLZ_GREY_INT:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
	            idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += srcGVWSp->gVal[0].inv;
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_SHORT:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += srcGVWSp->gVal[0].shv;
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_UBYTE:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += srcGVWSp->gVal[0].ubv;
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_FLOAT:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.dbp + idP) += srcGVWSp->gVal[0].flv;
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_DOUBLE:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.dbp + idP) += srcGVWSp->gVal[0].dbv;
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_RGBA:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGet(srcGVWSp, 0,
		  		  WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += WLZ_RGBA_RED_GET(
					   srcGVWSp->gVal[0].rgbv);
		    *(olpBuf.inp + bufWidth + idP) += WLZ_RGBA_GREEN_GET(
						  srcGVWSp->gVal[0].rgbv);
		    *(olpBuf.inp + (2 * bufWidth) + idP) += WLZ_RGBA_BLUE_GET(
						  srcGVWSp->gVal[0].rgbv);
		    *(olpBuf.inp + (3 * bufWidth) + idP) += WLZ_RGBA_ALPHA_GET(
						  srcGVWSp->gVal[0].rgbv);
		  }
		  ++idX;
		}
		break;
	    }
	    break;
	  case WLZ_INTERPOLATION_LINEAR:
	    switch(dstGWSp.pixeltype)
	    {
	      case WLZ_GREY_INT:
		while(idX <= iRgt)
		{
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD0 = ((srcGVWSp->gVal[0]).inv * tD2 * tD3) +
			  ((srcGVWSp->gVal[1]).inv * tD0 * tD3) +
			  ((srcGVWSp->gVal[2]).inv * tD2 * tD1) +
			  ((srcGVWSp->gVal[3]).inv * tD0 * tD1);
		    WLZ_CLAMP(tD0, 0.0, 255.0);
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += srcGVWSp->gVal[0].inv;
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
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD0 = ((srcGVWSp->gVal[0]).shv * tD2 * tD3) +
			  ((srcGVWSp->gVal[1]).shv * tD0 * tD3) +
			  ((srcGVWSp->gVal[2]).shv * tD2 * tD1) +
			  ((srcGVWSp->gVal[3]).shv * tD0 * tD1);
		    WLZ_CLAMP(tD0, 0.0, 255.0);
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += srcGVWSp->gVal[0].shv;
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
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD0 = ((srcGVWSp->gVal[0]).ubv * tD2 * tD3) +
			  ((srcGVWSp->gVal[1]).ubv * tD0 * tD3) +
			  ((srcGVWSp->gVal[2]).ubv * tD2 * tD1) +
			  ((srcGVWSp->gVal[3]).ubv * tD0 * tD1);
		    WLZ_CLAMP(tD0, 0.0, 255.0);
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += WLZ_NINT(tD0);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) += srcGVWSp->gVal[0].ubv;
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
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD0 = ((srcGVWSp->gVal[0]).flv * tD2 * tD3) +
			  ((srcGVWSp->gVal[1]).flv * tD0 * tD3) +
			  ((srcGVWSp->gVal[2]).flv * tD2 * tD1) +
			  ((srcGVWSp->gVal[3]).flv * tD0 * tD1);
		    WLZ_CLAMP(tD0, 0.0, 255.0);
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.dbp + idP) += WLZ_NINT(tD0);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += srcGVWSp->gVal[0].flv;
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
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD0 = ((srcGVWSp->gVal[0]).dbv * tD2 * tD3) +
			  ((srcGVWSp->gVal[1]).dbv * tD0 * tD3) +
			  ((srcGVWSp->gVal[2]).dbv * tD2 * tD1) +
			  ((srcGVWSp->gVal[3]).dbv * tD0 * tD1);
		    WLZ_CLAMP(tD0, 0.0, 255.0);
		    idP = idX - dstIWSp.lftpos;
		    ++*(olpCnt + idP);
		    *(olpBuf.dbp + idP) += WLZ_NINT(tD0);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.dbp + idP) += srcGVWSp->gVal[0].dbv;
		    }
		  }
		  ++idX;
		}
		break;
	      case WLZ_GREY_RGBA:
		while(idX <= iRgt)
		{
		  idP = idX - dstIWSp.lftpos;
		  sPosD.vtX = (trXX * idX) + trXYC;
		  sPosD.vtY = (trYX * idX) + trYYC;
		  WlzGreyValueGetCon(srcGVWSp, 0, sPosD.vtY, sPosD.vtX);
		  if(srcGVWSp->bkdFlag == 0)
		  {
		    tD0 = sPosD.vtX - floor(sPosD.vtX);
		    tD1 = sPosD.vtY - floor(sPosD.vtY);
		    tD2 = 1.0 - tD0;
		    tD3 = 1.0 - tD1;
		    tD4 = (WLZ_RGBA_RED_GET((srcGVWSp->gVal[0]).rgbv) *
			   tD2 * tD3) +
			  (WLZ_RGBA_RED_GET((srcGVWSp->gVal[1]).rgbv) *
			   tD0 * tD3) +
			  (WLZ_RGBA_RED_GET((srcGVWSp->gVal[2]).rgbv) *
			   tD2 * tD1) +
			  (WLZ_RGBA_RED_GET((srcGVWSp->gVal[3]).rgbv) *
			   tD0 * tD1);
		    WLZ_CLAMP(tD4, 0.0, 255.0);
		    ++*(olpCnt + idP);
		    *(olpBuf.inp + idP) += WLZ_NINT(tD4);
		    tD4 = (WLZ_RGBA_GREEN_GET((srcGVWSp->gVal[0]).rgbv) *
			   tD2 * tD3) +
			  (WLZ_RGBA_GREEN_GET((srcGVWSp->gVal[1]).rgbv) *
			   tD0 * tD3) +
			  (WLZ_RGBA_GREEN_GET((srcGVWSp->gVal[2]).rgbv) *
			   tD2 * tD1) +
			  (WLZ_RGBA_GREEN_GET((srcGVWSp->gVal[3]).rgbv) *
			   tD0 * tD1);
		    WLZ_CLAMP(tD4, 0.0, 255.0);
		    *(olpBuf.inp + bufWidth + idP) += WLZ_NINT(tD4);
		    tD4 = (WLZ_RGBA_BLUE_GET((srcGVWSp->gVal[0]).rgbv) *
			   tD2 * tD3) +
			  (WLZ_RGBA_BLUE_GET((srcGVWSp->gVal[1]).rgbv) *
			   tD0 * tD3) +
			  (WLZ_RGBA_BLUE_GET((srcGVWSp->gVal[2]).rgbv) *
			   tD2 * tD1) +
			  (WLZ_RGBA_BLUE_GET((srcGVWSp->gVal[3]).rgbv) *
			   tD0 * tD1);
		    WLZ_CLAMP(tD4, 0.0, 255.0);
		    *(olpBuf.inp + (2 * bufWidth) + idP) += WLZ_NINT(tD4);
		    tD4 = (WLZ_RGBA_ALPHA_GET((srcGVWSp->gVal[0]).rgbv) *
			   tD2 * tD3) +
			  (WLZ_RGBA_ALPHA_GET((srcGVWSp->gVal[1]).rgbv) *
			   tD0 * tD3) +
			  (WLZ_RGBA_ALPHA_GET((srcGVWSp->gVal[2]).rgbv) *
			   tD2 * tD1) +
			  (WLZ_RGBA_ALPHA_GET((srcGVWSp->gVal[3]).rgbv) *
			   tD0 * tD1);
		    WLZ_CLAMP(tD4, 0.0, 255.0);
		    *(olpBuf.inp + (3 * bufWidth) + idP) += WLZ_NINT(tD4);
		  }
		  else
		  {
		    WlzGreyValueGet(srcGVWSp, 0,
				    WLZ_NINT(sPosD.vtY), WLZ_NINT(sPosD.vtX));
		    if(srcGVWSp->bkdFlag == 0)
		    {
		      idP = idX - dstIWSp.lftpos;
		      ++*(olpCnt + idP);
		      *(olpBuf.inp + idP) +=
		      		WLZ_RGBA_RED_GET(srcGVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + bufWidth + idP) +=
		      		WLZ_RGBA_GREEN_GET(srcGVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + (2 * bufWidth) + idP) +=
		      		WLZ_RGBA_BLUE_GET(srcGVWSp->gVal[0].rgbv);
		      *(olpBuf.inp + (3 * bufWidth) + idP) +=
		      		WLZ_RGBA_ALPHA_GET(srcGVWSp->gVal[0].rgbv);
		    }
		  }
		  ++idX;
		}
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
      idX = dstIWSp.lftpos;
      switch(interp)
      {
	case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
	case WLZ_INTERPOLATION_LINEAR:
	  switch(dstGWSp.pixeltype)
	  {
	    case WLZ_GREY_INT:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
		*(dGP.inp + idP) = (*(olpCnt + idP) >= 1)?
		                   *(olpBuf.inp + idP) / *(olpCnt + idP):
				   bkdV.v.inv;
		++idX;
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
		*(dGP.shp + idP) = (*(olpCnt + idP) >= 1)?
				   *(olpBuf.inp + idP) / *(olpCnt + idP):
				   bkdV.v.shv;
		++idX;
	      }
	      break;
	    case WLZ_GREY_UBYTE:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
		*(dGP.ubp + idP) = (*(olpCnt + idP) >= 1)?
				   *(olpBuf.inp + idP) / *(olpCnt + idP):
				   bkdV.v.ubv;
		++idX;
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
		*(dGP.flp + idP) = (*(olpCnt + idP) >= 1)?
				   *(olpBuf.dbp + idP) / *(olpCnt + idP):
				   bkdV.v.flv;
		++idX;
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
		*(dGP.dbp + idP) = (*(olpCnt + idP) >= 1)?
				   *(olpBuf.dbp + idP) / *(olpCnt + idP):
				   bkdV.v.dbv;
		++idX;
	      }
	      break;
	    case WLZ_GREY_RGBA:
	      while(idX <= dstIWSp.rgtpos)
	      {
		idP = idX - dstIWSp.lftpos;
	        if(*(olpCnt + idP) >= 1)
		{
		  tI[0] = *(olpBuf.inp + idP) / *(olpCnt + idP);
		  tI[1] = *(olpBuf.inp + bufWidth + idP) / *(olpCnt + idP);
		  tI[2] = *(olpBuf.inp + (2 * bufWidth) + idP) /
		          *(olpCnt + idP);
		  tI[3] = *(olpBuf.inp + (3 * bufWidth) + idP) /
		          *(olpCnt + idP);
		  WLZ_RGBA_RGBA_SET(*(dGP.rgbp + idP),
		                    tI[0], tI[1], tI[2], tI[3]);
		}
		else
		{
		  *(dGP.rgbp + idP) = bkdV.v.rgbv;
		}
		++idX;
	      }
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
    }
  }
  if(errNum == WLZ_ERR_EOO)
  {
    errNum = WLZ_ERR_NONE;
  }
  AlcFree(olpBuf.inp);
  AlcFree(olpCnt);
  WlzCMeshScanWSpFree(mSWSp);
  WlzGreyValueFreeWSp(srcGVWSp);
  return(errNum);
}

/*!
* \return	New conforming mesh scan workspace.
* \ingroup	WlzTransform
* \brief	Allocate and initialise a conforming mesh scan workspace.
* \param	mTr			Conforming mesh transform.
* \param	dstErr			Destination error pointer.
*/
static WlzCMeshScanWSp *WlzCMeshScanWSpInit(WlzCMeshTransform *mTr,
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
  WlzCMeshScanWSp *mSWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmRes = &(mTr->mesh.m2->res.elm);
  if((mSWSp = (WlzCMeshScanWSp *)
  	      AlcCalloc(1, sizeof(WlzCMeshScanWSp))) == NULL)
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
	nod = elm->edg[0].nod;
	dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	eLnMin = eLnMax = nod->pos.vtY + dsp->vtY;
	nod = elm->edg[1].nod;
	dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	ndLn = nod->pos.vtY + dsp->vtY;
	if(ndLn < eLnMin)
	{
	  eLnMin = ndLn;
	}
	else if(ndLn > eLnMax)
	{
	  eLnMax = ndLn;
	}
	nod = elm->edg[2].nod;
	dsp = (WlzDVertex2 *)AlcVectorItemGet(mTr->dspVec, nod->idx);
	ndLn = nod->pos.vtY + dsp->vtY;
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
    if(((mSWSp->itvs = (WlzCMeshScanItv *)
    		       AlcMalloc(sizeof(WlzCMeshScanItv) *
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
        iIdx += WlzCMeshScanTriElm(mSWSp, elm, iIdx);
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
    qsort(mSWSp->itvs, mSWSp->nItvs, sizeof(WlzCMeshScanItv),
          WlzCMeshItvCmp);
#ifdef WLZ_CMESHTRANSFORM_DEBUG
    for(iIdx = 0; iIdx < mSWSp->nItvs; ++iIdx)
    {
      (void )fprintf(stderr,
      		     "%d %d %d %d\n",
		     (mSWSp->itvs + iIdx)->elmIdx,
		     (mSWSp->itvs + iIdx)->lftI, (mSWSp->itvs + iIdx)->rgtI,
		     (mSWSp->itvs + iIdx)->line);
    }
#endif
  }
  else
  {
    WlzCMeshScanWSpFree(mSWSp);
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
* \brief	Free's a conforming mesh scan workspace.
* \param	mSnWSp			Conforming mesh scan workspace.
*/
static void	WlzCMeshScanWSpFree(WlzCMeshScanWSp *mSWSp)
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
static int	WlzCMeshScanTriElm(WlzCMeshScanWSp *mSWSp, WlzCMeshElm2D *elm,
				   int iIdx)
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
  WlzCMeshScanItv *itv;

  /* Compute the integer displaced nodes of the element. */
  for(ndIdx0 = 0; ndIdx0 < 3; ++ndIdx0)
  {
    nod = elm->edg[ndIdx0].nod;
    dVx0 = nod->pos;
    dVx1 = *(WlzDVertex2 *)AlcVectorItemGet(mSWSp->mTr->dspVec, nod->idx);
    tD0 = dVx0.vtX + dVx1.vtX;
    tD1 = dVx0.vtY + dVx1.vtY;
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
* \return	Sorting value for qsort.
* \ingroup	WlzTransform
* \brief	Callback function for qsort(3) to sort conforming mesh
*		element intervals by line and then left column.
* \param	cmp0			Used to pass first mesh interval.
* \param	cmp1			Used to pass second mesh interval.
*/
static int	WlzCMeshItvCmp(const void *cmp0, const void *cmp1)
{
  int		rtn;
  WlzCMeshScanItv *itv0,
  		 *itv1;

  itv0 = (WlzCMeshScanItv *)cmp0;
  itv1 = (WlzCMeshScanItv *)cmp1;
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
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Applies a 2D conforming mesh transform to the given source
*		object which must be a 2D object.
* \param	srcObj			Object to be transformed.
* \param	mTr			Conforming mesh transform.
* \param	interp			Type of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshTransformObj2D(WlzObject *srcObj,
				     WlzCMeshTransform *mTr,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzDomain	dstDom;
  WlzValues	srcValues;
  WlzObject	*tObj0 = NULL,
  		*tObj1 = NULL,
		*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
 
  dstDom.core = NULL;
  srcValues.core = NULL;
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
	  tObj1 = WlzCMeshTransformObj2D(tObj0, mTr, interp, &errNum);
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
	errNum = WlzCMeshTransformValues2D(dstObj, srcObj, mTr, interp);
      }
      break;
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
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
	    edg[eId1++] = tbl[mElm->edg[0].nod->idx];
	    edg[eId1++] = tbl[mElm->edg[1].nod->idx];
	    edg[eId1++] = tbl[mElm->edg[2].nod->idx];
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

#ifdef WLZ_CMESH_DEBUG_MAIN
/*!
* \return       void
* \ingroup      WlzTransform
* \brief        Debuging function for 2D mesh output in VTK format.
* \param        fP                      Given file pointer.
* \param        mesh                    Given mesh.
* \param	dspFlg			Non zero for displaced mesh output.
*/
void            WlzCMeshTransDbgOutVTK2D(FILE *fP, WlzCMesh2D *mesh,
					 int dspFlg)
{
  int           idE,
                idN,
                bCnt,
                nElm,
                nVElm,
                nNod;
  WlzDVertex2	*dsp;
  WlzCMeshElm2D  *elm;
  WlzCMeshNod2D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D) &&
    ((nNod = mesh->res.nod.maxEnt) > 0) &&
    ((nElm = mesh->res.elm.maxEnt) > 0))
  {
    (void )fprintf(fP,
                   "# vtk DataFile Version 1.0\n"
                   "WlzCMesh2D 2D\n"
                   "ASCII\n"
                   "DATASET POLYDATA\n"
                   "POINTS %d float\n",
                   nNod);
    for(idN = 0; idN < nNod; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	if(dspFlg)
	{
	  dsp = (WlzDVertex2 *)(nod->prop);
	  (void )fprintf(fP, "%g %g 0\n",
			 nod->pos.vtX + dsp->vtX,
			 nod->pos.vtY + dsp->vtY);
	}
	else
	{
	  (void )fprintf(fP, "%g %g 0\n",
			 nod->pos.vtX, nod->pos.vtY);
        }
      }
      else
      {
        (void )fprintf(fP, "0 0 0\n");
      }
    }
    nVElm = 0;
    for(idE = 0; idE < nElm; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        ++nVElm;
      }
    }
    (void )fprintf(fP, "POLYGONS %d %d\n",
                   nVElm, nVElm * 4);
    for(idE = 0; idE < nElm; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        (void )fprintf(fP, "3 %d %d %d\n",
                       elm->edg[0].nod->idx, elm->edg[1].nod->idx,
                       elm->edg[2].nod->idx);
      }
    }
  }
}

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		dspFlg = 0,
  		ok = 1,
  		option,
  		usage = 0;
  double	minElmSz = 25.0,
  		maxElmSz = 100.0;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  WlzCMesh2D 	*mesh = NULL;
  WlzCMeshTransform *mTr = NULL;
  static char   optList[] = "dhm:M:o:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
        dspFlg = 1;
	break;
      case 'm':
        if(sscanf(optarg, "%lg", &minElmSz) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%lg", &maxElmSz) != 1)
	{
	  usage = 1;
	}
        break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    (void )WlzAssignObject(obj, NULL);
    mesh = WlzCMeshFromObj2D(obj, minElmSz, maxElmSz, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to create conforming mesh, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((mTr = WlzMakeCMeshTransform2D(mesh, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to make mesh transform from mesh, %s.\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    WlzCMeshTransDbgOutVTK2D(fP, mesh, dspFlg);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>] [-m#] [-M#] [<input object>]\n"
    	    "Computes a conforming mesh for the given input object.\n"
	    "Options are:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n"
	    "  -d  Output displaced mesh.\n"
	    "  -m  Minimum mesh element size.\n"
	    "  -M  Maximum mesh element size.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* WLZ_CMESH_DEBUG_MAIN */
