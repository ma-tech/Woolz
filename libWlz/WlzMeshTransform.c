#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMeshTransform.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz functions for computing mesh transforms.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

/*!
* \struct	_WlzMeshScanDElm
* \ingroup	WlzTransform
* \brief	Mesh scanning element.
*/
typedef struct _WlzMeshScanDElm
{
  int		valid;			/*! Non-zero if valid. */
  double	xTr[3];         	/*! Affine transform coefficients
  					    for columns. */
  double	yTr[3];           	/*! Affine transform coefficients
  					    for lines. */
} WlzMeshScanDElm;

/*!
* \struct	_WlzMeshScanItv
* \ingroup	WlzTransform
* \brief	Scan interval within an element.
*/
typedef struct _WlzMeshScanItv
{
  int		elmIdx;			/*! Element index. */
  int		line;			/*! Line of interval. */
  int		lftI;			/*! Start of interval. */
  int		rgtI;			/*! End of interval. */
} WlzMeshScanItv;

/*!
* \struct	_WlzMeshScanWSp
* \ingroup	WlzTransform
* \brief	Mesh scanning workspace.
*/
typedef struct _WlzMeshScanWSp
{
  WlzMeshTransform *mesh;		/*! The mesh transform. */
  int		nItvs; 			/*! Number of element intervals. */
  WlzMeshScanItv *itvs; 		/*! Element intervals sorted by line
  					    then left column. */
  WlzMeshScanDElm *dElm; 		/*! Destination mesh element data. */
} WlzMeshScanWSp;

/*!
* \struct
* \ingroup	WlzTransform
* \brief	Linked list based polygon data structure.
*/
typedef struct	_WlzMeshPolyVx
{
  struct _WlzMeshPolyVx *prev;		/*! Next vertex in polygon. */
  struct _WlzMeshPolyVx *next;		/*! Previous vertex in polygon. */
  int		id;			/*! Index of the element. */
  WlzDVertex2	vx;			/*! Vertex position. */
} WlzMeshPolyVx;


static int			WlzMeshItvCmp(
				  const void *,
				  const void *),
				WlzMeshScanTriElm(
				  WlzMeshScanWSp *,
				  int,
				  int);
static void			WlzMeshScanWSpFree(
				  WlzMeshScanWSp *);
static void			WlzMeshAfTrSolve(
				  double *, 
				  double *,
				  double,
				  WlzDVertex2 *,
				  WlzDVertex2 *);
static unsigned int 		WlzMeshTransFillBlockLnNod(
				  WlzMeshNode *,
				  WlzInterval *,
			       	  WlzIVertex2,
				  unsigned int);
static WlzErrorNum 		WlzMeshRemoveObjBox(
				  WlzMeshTransform *),
				WlzMeshBoundCvPolyAdd(
				  WlzMeshTransform *,
				  WlzObject *,
				  double,
				  double),
				WlzMeshBoundPolyFix(
				  WlzObject *,
				  double),
				WlzMeshScanDElmUpdate(
				  WlzMeshScanWSp *,
				  int),
				WlzMeshLineCvExtrema(
				  WlzInterval *,
				  WlzObject *,
				  unsigned int,
				  unsigned int),
				WlzMeshTransFillBlock(
				  WlzMeshTransform *mesh,
				  WlzInterval *,
				  WlzIVertex2,
				  unsigned int,
				  unsigned int),
				WlzMeshTransformVxVecI(
				  WlzMeshTransform *,
				  WlzIVertex2 *,
			  	  int),
				WlzMeshTransformVxVecF(
				  WlzMeshTransform *,
				  WlzFVertex2 *,
			  	  int),
				WlzMeshTransformVxVecD(
				  WlzMeshTransform *,
				  WlzDVertex2 *,
			  	  int),
				WlzMeshTransformValues2D(
				  WlzObject *,
				  WlzObject *,
				  WlzMeshTransform *,
				  WlzInterpolationType);
static WlzObject 		*WlzMeshTransformObjPrv(
				  WlzObject *,
				  WlzMeshTransform *,
				  WlzInterpolationType,
				  WlzErrorNum *);
static WlzMeshTransform 	*WlzMeshFromObjBlock(
				  WlzObject *,
				  unsigned int,
			          WlzErrorNum *),
        			*WlzMeshFromObjGrad(
				  WlzObject *,
				  unsigned int,
				  unsigned int,
				  WlzErrorNum *);
static WlzPolygonDomain		*WlzMeshTransformPoly(
				  WlzPolygonDomain *,
			      	  WlzMeshTransform *,
			      	  WlzErrorNum *);
static WlzBoundList 		*WlzMeshTransformBoundList(
				  WlzBoundList *,
				  WlzMeshTransform *,
				  WlzErrorNum *);
static WlzMeshScanWSp 		*WlzMeshScanWSpInit(
				  WlzMeshTransform *,
				  WlzErrorNum *);

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Free's the given mesh transform.
* \param	mesh			Given mesh transform.
*/
WlzErrorNum	WlzMeshFreeTransform(WlzMeshTransform *mesh)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh)
  {
    if((mesh->maxElem > 0) && mesh->elements)
    {
      AlcFree(mesh->elements);
    }
    if((mesh->maxNodes > 0) && mesh->nodes)
    {
      AlcFree(mesh->nodes);
    }
    AlcFree(mesh);
  }
  return(errNum);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a mesh transform data structure with the nodes,
*		elements and displacements allocated and initialized to zero.
* \param	nNode			Number of nodes (and displacements).
* \param	nElem			Number of elements.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform *WlzMeshTransformNew(unsigned int nNode,
				      unsigned int nElem,
				      WlzErrorNum *dstErr)
{
  WlzMeshTransform *mesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((mesh = (WlzMeshTransform *)
  	      AlcCalloc(1, sizeof(WlzMeshTransform))) == NULL) ||
     ((mesh->elements = (WlzMeshElem *)
     		        AlcCalloc(nElem, sizeof(WlzMeshElem))) == NULL) ||
     ((mesh->nodes = (WlzMeshNode *)
     		     AlcCalloc(nNode, sizeof(WlzMeshNode))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
    if(mesh)
    {
      if(mesh->elements)
      {
        AlcFree(mesh->elements);
      }
      if(mesh->nodes)
      {
        AlcFree(mesh->nodes);
      }
      AlcFree(mesh);
      mesh = NULL;
    }
  }
  else
  {
    mesh->type = WLZ_TRANSFORM_2D_MESH;
    mesh->maxNodes = nNode;
    mesh->maxElem = nElem;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Adapts the given mesh transform so that each of the
*		elements in the displaced mesh has an area greater
*		than the given minimum.
* \param	gMesh			Given mesh transform.
* \param	minArea			Minimum element area.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform *WlzMeshTransformAdapt(WlzMeshTransform *gMesh,
					double minArea, WlzErrorNum *dstErr)
{
  int		tI0,
  		eId0,
		eNodId,
		dNodId,
		dCnt,
  		nId0;
  double	tD0,
  		tD1,
		minArea2,
  		dArea2,
		sArea2;
  WlzMeshElem	*elm0;
  WlzMeshNode	*nod0;
  WlzMeshTransform *nMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  double	dSeg[3];
  WlzDVertex2	sVx[3],
  		dVx[3];
  const unsigned int maxSegTbl[8] = {0, 0, 1, 0, 2, 2, 1, 0},
  		nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				WLZ_MESH_ELEM_FLAGS_NBR_1,
				WLZ_MESH_ELEM_FLAGS_NBR_2};

  if(gMesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gMesh->nElem < 0) || (gMesh->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    if(minArea < 1.0)
    {
      minArea = 1.0;
    }
    nMesh = WlzMeshTransformCopy(gMesh, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    minArea2 = minArea * 2; 
    /* Find mesh element with area < given minimum area. */
    do
    {
      dCnt = 0;
      for(eId0 = 0; (errNum == WLZ_ERR_NONE) && (eId0 < nMesh->nElem); ++eId0)
      {
        elm0 = nMesh->elements + eId0;
	if((elm0->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
	{
	  /* Compute displaced vertices. */
	  for(nId0 = 0; nId0 < 3; ++nId0)
	  {
	    nod0 = nMesh->nodes + elm0->nodes[nId0];
	    sVx[nId0] = nod0->position;
	    dVx[nId0].vtX = nod0->position.vtX + nod0->displacement.vtX;
	    dVx[nId0].vtY = nod0->position.vtY + nod0->displacement.vtY;
	  }
	  sArea2 = WlzGeomTriangleSnArea2(sVx[0], sVx[1], sVx[2]);
	  dArea2 = WlzGeomTriangleSnArea2(dVx[0], dVx[1], dVx[2]);
	  if((dArea2 < minArea2) || (sArea2 < minArea2))
	  {
	    /* Compute lengths of the displacemed element segments. */
	    tD0 = dVx[2].vtX - dVx[1].vtX;
	    tD1 = dVx[2].vtY - dVx[1].vtY;
	    dSeg[0] = (tD0 * tD0) + (tD1 * tD1);
	    tD0 = dVx[0].vtX - dVx[2].vtX;
	    tD1 = dVx[0].vtY - dVx[2].vtY;
	    dSeg[1] = (tD0 * tD0) + (tD1 * tD1);
	    tD0 = dVx[1].vtX - dVx[0].vtX;
	    tD1 = dVx[1].vtY - dVx[0].vtY;
	    dSeg[2] = (tD0 * tD0) + (tD1 * tD1);
	    /* Find longest segment. */
	    tI0 = ((dSeg[2] > dSeg[0]) << 2) |
		  ((dSeg[1] > dSeg[2]) << 1) |
		  (dSeg[0] > dSeg[1]);
	    eNodId = maxSegTbl[tI0];
	    /* Choose a node to delete. */
	    if((elm0->flags & nbrFlgTbl[eNodId]) != 0)
	    {
	      /* Node is on mesh boundary. */
	      dNodId = elm0->nodes[eNodId];
	    }
	    else
	    {
	      /* Element is not on mesh boundary: Delete node opposite the
	       * longest element segment. */
	      dNodId = elm0->nodes[eNodId];
	    }
	    errNum = WlzMeshNodeDelIdx(nMesh, elm0->idx, &dNodId, 1);
	    ++dCnt;
	  }
	}
      }
    }
    while((errNum == WLZ_ERR_NONE) && (dCnt > 0));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzMeshSqueeze(nMesh);
  }
#ifdef WLZ_MESH_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzMeshTransformVerify(nMesh, 0, NULL, NULL);
  }
#endif /* WLZ_MESH_DEBUG */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nMesh);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Copies the given mesh transform. The copied mesh will have any
*		zombie elements squeezed out.
* \param	gMesh			Given mesh transform.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform *WlzMeshTransformCopy(WlzMeshTransform *gMesh,
				       WlzErrorNum *dstErr)
{
  WlzMeshTransform *nMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gMesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gMesh->nElem < 0) || (gMesh->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    nMesh = WlzMeshTransformNew(gMesh->nNodes, gMesh->nElem, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nMesh->nNodes = gMesh->nNodes;
    (void )memcpy(nMesh->nodes, gMesh->nodes,
    		  gMesh->nNodes * sizeof(WlzMeshNode));
    nMesh->nElem = gMesh->nElem;
    (void )memcpy(nMesh->elements, gMesh->elements,
    		  gMesh->nElem * sizeof(WlzMeshElem));
    errNum = WlzMeshSqueeze(nMesh);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nMesh);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a mesh transform for the given object with all mesh
* 		displacements zero.
* \param	srcObj			The given object.
* \param	method			Mesh generation method to use.
* \param	minDist			Minimum distance between mesh vertices.
* \param	maxDist			Maximum distance between mesh vertices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform *WlzMeshFromObj(WlzObject *srcObj, WlzMeshGenMethod method,
				 double minDist, double maxDist,
				 WlzErrorNum *dstErr)
{
  WlzMeshTransform *mesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
  else
  {
    switch(method)
    {
      case WLZ_MESH_GENMETHOD_BLOCK:
	mesh = WlzMeshFromObjBlock(srcObj,
				   (unsigned int )abs(WLZ_NINT(minDist)),
				   &errNum);
        break;
      case WLZ_MESH_GENMETHOD_GRADIENT:
        mesh = WlzMeshFromObjGrad(srcObj, 
				  (unsigned int )abs(WLZ_NINT(minDist)),
				  (unsigned int )abs(WLZ_NINT(maxDist)),
				   &errNum);
        break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup 	WlzTransform
* \brief	Transforms a woolz object using a the given mesh transform.
* \param	srcObj			Object to be transformed.
* \param	gMesh			Given mesh transform to apply.
* \param	interp			Level of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzMeshTransformObj(WlzObject *srcObj,
				     WlzMeshTransform *gMesh,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzMeshTransform *aMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double minElmArea = 1.0;

  if(srcObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gMesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gMesh->nElem < 0) || (gMesh->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    aMesh = WlzMeshTransformAdapt(gMesh, minElmArea, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dstObj = WlzMeshTransformObjPrv(srcObj, aMesh, interp, &errNum);
  }
  if(aMesh)
  {
    (void )WlzMeshFreeTransform(aMesh);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Computes a mesh transform for the given object and a
*		set of control points using both an affine and a basis
*		function transform to set the mesh displacements.
* \param	obj			Given object.
* \param	basisFnType		Required basis function type.
* \param	polyOrder		Order of polynomial, only used for
*					WLZ_FN_BASIS_2DPOLY.
* \param	nSPts			Number of source control points.
* \param	sPts			Source control points.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	meshGenMtd		Mesh generation method.
* \param	meshMinDist		Minimum mesh vertex distance.
* \param	meshMaxDist		Maximum mesh vertex distance.
* \param	dstErr			Destination error pointer, may be
*					NULL.
*/
WlzMeshTransform *WlzMeshTransformFromCPts(WlzObject *obj,
				WlzFnType basisFnType, int polyOrder,
				int nSPts, WlzDVertex2 *sPts,
				int nDPts, WlzDVertex2 *dPts,
				WlzMeshGenMethod meshGenMtd,
				double meshMinDist, double meshMaxDist,
				WlzErrorNum *dstErr)
{
  int		idx;
  WlzDVertex2	tDV0;
  WlzDVertex2	*dPtsT = NULL;
  WlzAffineTransform *aTr = NULL,
  		*aTrI = NULL;
  WlzBasisFnTransform *bTr = NULL;
  WlzMeshTransform *mTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((nDPts <= 0) || (nDPts != nSPts))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Compute least squares affine transform from the tie points. */
  if(errNum == WLZ_ERR_NONE)
  {
    aTr = WlzAffineTransformLSq2D(nSPts, sPts, nSPts, dPts,
				  WLZ_TRANSFORM_2D_AFFINE, &errNum);
  }
  /* Compute a mesh transform for the given object. */
  if(errNum == WLZ_ERR_NONE)
  {
    mTr = WlzMeshFromObj(obj, meshGenMtd, meshMinDist, meshMaxDist, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nSPts >= 4)
    {
      /* Create a new array of destination vertices which have the original
       * destination transformed by the inverse affine transform. */
      if((dPtsT = (WlzDVertex2 *)AlcMalloc(sizeof(WlzDVertex2) *
	      				   nSPts)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        aTrI = WlzAffineTransformInverse(aTr, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        for(idx = 0; idx < nDPts; ++idx)
	{
          dPtsT[idx] = WlzAffineTransformVertexD2(aTrI, dPts[idx], NULL);
	}
        bTr = WlzBasisFnTrFromCPts2D(basisFnType, polyOrder,
			  	   nSPts, sPts, nSPts, dPtsT, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        /* Set the mesh transform displacements and then apply the affine
	 * transform. */
        errNum = WlzBasisFnSetMesh(mTr, bTr);
      }
    }
    /* Apply the affine transform to the mesh transform. */
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzMeshAffineProduct(mTr, aTr);
    }
  }
  AlcFree(dPtsT);
  (void )WlzBasisFnFreeTransform(bTr);
  (void )WlzFreeAffineTransform(aTr);
  (void )WlzFreeAffineTransform(aTrI);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzMeshFreeTransform(mTr);
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
* \brief	Computes the product of the given affine and mesh transforms
*		in place, ie the mesh transform has it's displacements
* 		overwritten.
* \param	mTr		Given mesh transform.
* \param	aTr		Given affine transform.
*/
WlzErrorNum WlzMeshAffineProduct(WlzMeshTransform *mTr,
				 WlzAffineTransform *aTr)
{
  int           count;
  WlzDVertex2   tDV0;
  WlzMeshNode	*node;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((mTr == NULL) || (aTr == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mTr->nodes == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  /* Loop through nodes, resetting the displacement. */
  if(errNum == WLZ_ERR_NONE)
  {
    node = mTr->nodes;
    count = mTr->nNodes;
    while(count-- > 0)
    {
      WLZ_VTX_2_ADD(tDV0, node->position, node->displacement);
      tDV0 = WlzAffineTransformVertexD2(aTr, tDV0, &errNum);
      WLZ_VTX_2_SUB(node->displacement, tDV0, node->position);
      ++node;
    }
  }
  return(errNum);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Private version of WlzMeshTransformObj() which transforms
*		a woolz object using a the given mesh transform.
* \param	srcObj			Object to be transformed.
* \param	mesh			Given mesh transform to apply.
* \param	interp			Level of interpolation.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzMeshTransformObjPrv(WlzObject *srcObj,
				     WlzMeshTransform *mesh,
				     WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  WlzDomain	dstDom;
  WlzValues	srcValues;
  WlzObject	*tObj0,
  		*tObj1,
		*dstObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstDom.core = NULL;
  srcValues.core = NULL;
  switch(srcObj->type)
  {
    case WLZ_EMPTY_OBJ:
      dstObj = WlzMakeEmpty(&errNum);
      break;
    case WLZ_2D_POLYGON:
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
	    dstDom.poly = WlzMeshTransformPoly(srcObj->domain.poly,
					       mesh, &errNum);
	    break;
	  case WLZ_BOUNDLIST:
	    dstDom.b = WlzMeshTransformBoundList(srcObj->domain.b,
						 mesh, &errNum);
	    break;
	}
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
	  tObj1 = WlzMeshTransformObjPrv(tObj0, mesh, interp, &errNum);
	  WlzFreeObj(tObj0);
	  tObj0 = NULL;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzBoundToObj(tObj1->domain.b,
				 WLZ_EVEN_ODD_FILL, &errNum);
	  WlzFreeObj(tObj1);
	  tObj1 = NULL;
	}
	if((errNum == WLZ_ERR_NONE) &&
	   (srcObj->values.core))
	{
	  errNum = WlzMeshTransformValues2D(dstObj, srcObj, mesh, interp);
	}
	if(tObj0)
	{
	  WlzFreeObj(tObj0);
	}
	if(tObj1)
	{
	  WlzFreeObj(tObj1);
	}
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
* \return	Woolz error code.
* \ingroup 	WlzTransform
* \brief	Checks that the given mesh transform is valid.
* \param	mesh			Given mesh transform.
* \param	dispFlg			Verify displacements if non-zero.
* \param	badElm			Destination ptr for the index of the
* 					first bad mesh element.
* \param	dstErrMsk		Destination pointer to be set with a
* 					mesh error after verifying this
* 					element, may be NULL.
*/
WlzErrorNum	WlzMeshTransformVerify(WlzMeshTransform *mesh, int dispFlg,
				       int *badElm, WlzMeshError *dstErrMsk)
{
  int		eIdx = 0;
  WlzMeshElem	*elm;
  WlzMeshError  errMsk = WLZ_MESH_ERR_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh->nElem < 0) || (mesh->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    eIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (eIdx < mesh->nElem))
    {
      elm = mesh->elements + eIdx;
      if((elm->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
      {
	if((errNum = WlzMeshElemVerify(mesh, dispFlg,
				       elm, &errMsk)) == WLZ_ERR_NONE)
	{
	  ++eIdx;
	}
      }
      else
      {
        ++eIdx;
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(badElm)
    {
      *badElm = eIdx;
    }
    if(dstErrMsk)
    {
      *dstErrMsk = errMsk;
    }
  }
  return(errNum);
}

/*!
* \return	Transformed polygon domain or NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms the given polygon domain using the given mesh
* 		transform.
* \param	srcPoly			Given polygon domain.
* \param	mesh			Mesh transform to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzPolygonDomain	*WlzMeshTransformPoly(WlzPolygonDomain *srcPoly,
					      WlzMeshTransform *mesh,
					      WlzErrorNum *dstErr)
{
  WlzPolygonDomain *dstPoly = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
        errNum = WlzMeshTransformVxVecI(mesh, dstPoly->vtx,
					dstPoly->nvertices);
	break;
      case WLZ_POLYGON_FLOAT:
        errNum = WlzMeshTransformVxVecF(mesh, (WlzFVertex2 *)(dstPoly->vtx),
					dstPoly->nvertices);
	break;
      case WLZ_POLYGON_DOUBLE:
        errNum = WlzMeshTransformVxVecD(mesh, (WlzDVertex2 *)(dstPoly->vtx),
					dstPoly->nvertices);
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPoly);
}

/*!
* \return	Transformed boundary list or NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms the given boundary list using the given mesh
*		transform.
* \param	srcBound		Given boundary list.
* \param	mesh			Mesh transform to apply.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzBoundList *WlzMeshTransformBoundList(WlzBoundList *srcBound,
					       WlzMeshTransform *mesh,
					       WlzErrorNum *dstErr)
{
  WlzDomain     dumDom;
  WlzBoundList  *dstBound = NULL;
  WlzObject	*polyobj;
  WlzPolygonDomain	*pgdm1, *pgdm2;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((dstBound = (WlzBoundList *)AlcCalloc(sizeof(WlzBoundList), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* wrap set to 1 for closed lines by WlzPolyDecimate */
    dstBound->type = srcBound->type;
    dstBound->wrap = srcBound->wrap?1:0;

    /* Transform the polygon */
    if( polyobj = WlzPolyTo8Polygon(srcBound->poly, srcBound->wrap, &errNum) ){
      if( pgdm1 = WlzMeshTransformPoly(polyobj->domain.poly, mesh,&errNum) ){
	if( pgdm2 = WlzPolyDecimate(pgdm1, srcBound->wrap, 0.75, &errNum) ){
	  dstBound->poly = WlzAssignPolygonDomain(pgdm2, NULL);
	}
	else {
	  dstBound->poly = NULL;
	}
	WlzFreePolyDmn(pgdm1);
      }
      WlzFreeObj(polyobj);
    }

    /* transform the remaining boundlists */
    if(errNum == WLZ_ERR_NONE)
    {
      /* Transform next */
      if(srcBound->next)
      {
	if((dumDom.b = WlzMeshTransformBoundList(srcBound->next, mesh,
						 &errNum)) != NULL)
	{
	  (void )WlzAssignDomain(dumDom, &errNum);
	  dstBound->next = dumDom.b;
	}
      }
      /* Transform down */
      if(srcBound->down && (errNum == WLZ_ERR_NONE))
      {
	if((dumDom.b = WlzMeshTransformBoundList(srcBound->down, mesh,
						 &errNum)) != NULL)
	{
	  (void )WlzAssignDomain(dumDom, &errNum);
	  dstBound->down = dumDom.b;
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBound);
}

/*!
* \return	Pixel value.
* \ingroup	WlzTransform
* \brief	Interpolate pixel value using maximum probability.
* \param	g			Array of four values at integer
* 					coordinates.
* \param	p			Column offset from integer coordinate.
* \param	q			Line offset from integer coordinate.
*/
double 		WlzClassValCon4(double *g, double p, double q)
{
  double	classVal[4], classProb[4], maxProb;
  int		i, j, nClass;

  /* initialise probability matrix */
  for(i=0; i < 4; i++){classProb[i] = 0.0;}
  /* determine classes and probabilities */
  nClass = 0;
  for(i=0; i < 4; i++){
    for(j=0; j < nClass; j++){
      if( g[i] == classVal[j] ){
	break;
      }
    }
    if( j == nClass ){
      classVal[j] = g[i];
      nClass++;
    }
    switch( i ){
    case 0:
      classProb[j] += (1-p)*(1-q);
      break;
    case 1:
      classProb[j] += p*(1-q);
      break;
    case 2:
      classProb[j] += (1-p)*q;
      break;
    case 3:
      classProb[j] += p*q;
      break;
    }
  }
  /* find max probability */
  maxProb = 0.0;
  for(i=0; i < nClass; i++){
    if( classProb[i] > maxProb ){
      maxProb = classProb[i];
      j = i;
    }
  }
  return classVal[j];
}
  
/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Creates a new value table, fills in the values and adds it to
* 		the given new object.
* \param	dstObj			Partialy transformed object with a
* 					valid domain.
* \param	srcObj			2D domain object which is being
*					transformed.
* \param	mesh			Given mesh transform.
* \param	interp			Level of interpolation.
*/
static WlzErrorNum WlzMeshTransformValues2D(WlzObject *dstObj,
					    WlzObject *srcObj,
					    WlzMeshTransform *mesh,
					    WlzInterpolationType interp)
{
  int		mItvIdx, indx;
  double	tD0,
  		tD1,
		tD2,
		tD3,
		trXX,
		trXYC,
		trYX,
		trYYC;
  double	gTmp[4];
  WlzIVertex2	dPosI,
  		sPosI;
  WlzDVertex2	sPosD;
  WlzGreyP	dGP;
  WlzGreyType	newGreyType;
  WlzPixelV	bkdV;
  WlzValues	newValues;
#ifdef WLZ_MESH_DEBUG
  WlzMeshScanItv oldItvDebug;
#endif /* WLZ_MESH_DEBUG */
  WlzMeshScanItv *mItv;
  WlzMeshScanWSp *mSnWSp = NULL;
  WlzMeshScanDElm *dElm;
  WlzGreyValueWSpace *gVWSp = NULL;
  WlzGreyWSpace gWSp;
  WlzIntervalWSpace iWSp;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  newValues.core = NULL;
  bkdV = WlzGetBackground(srcObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    newGreyType = WlzGreyTableTypeToGreyType(srcObj->values.v->type,
    					     &errNum);
  }
  if((errNum == WLZ_ERR_NONE) && (bkdV.type != newGreyType))
  {
    errNum = WlzValueConvertPixel(&bkdV, bkdV, newGreyType);
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
    mItvIdx = 0;
    mSnWSp = WlzMeshScanWSpInit(mesh, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mItv = mSnWSp->itvs;
    errNum = WlzInitGreyScan(dstObj, &iWSp, &gWSp);
    if(errNum == WLZ_ERR_NONE)
    {
      gVWSp = WlzGreyValueMakeWSp(srcObj, &errNum);
    }
    while((errNum == WLZ_ERR_NONE) &&
          (WlzNextGreyInterval(&iWSp) == 0))
    {
      dPosI.vtX = iWSp.lftpos;
      dPosI.vtY = iWSp.linpos;
      dGP = gWSp.u_grintptr;
      while((errNum == WLZ_ERR_NONE) && (dPosI.vtX <= iWSp.rgtpos))
      {
	/* Find the appropriate mesh scan interval. */
#ifdef WLZ_MESH_DEBUG
	oldItvDebug = *mItv;
#endif /* WLZ_MESH_DEBUG */
        while((mItv->line < dPosI.vtY) && (mItvIdx < mSnWSp->nItvs))
	{
	  ++mItvIdx;
	  ++mItv;
	}
	while((mItv->line <= dPosI.vtY) &&
	      (mItv->rgtI < dPosI.vtX) && (mItvIdx < mSnWSp->nItvs))
	{
	  ++mItvIdx;
	  ++mItv;
	}
	if((mItv->line != dPosI.vtY) ||
	   (dPosI.vtX < mItv->lftI) || (dPosI.vtX > mItv->rgtI))
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
#ifdef WLZ_MESH_DEBUG
  (void )fprintf(stderr,
  "Debug oldItvDebug %d %d %d %d, mItv %d %d %d %d, dPos %d %d\n",
  oldItvDebug.elmIdx, oldItvDebug.line, oldItvDebug.lftI, oldItvDebug.rgtI,
  mItv->elmIdx, mItv->line, mItv->lftI, mItv->rgtI,
  dPosI.vtX, dPosI.vtY);
#endif /* WLZ_MESH_DEBUG */
	}
	else
	{
	  dElm = mSnWSp->dElm + mItv->elmIdx;
	  if(dElm->valid == 0)
	  {
	    errNum = WlzMeshScanDElmUpdate(mSnWSp, mItv->elmIdx);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    trXX = dElm->xTr[0];
	    trXYC = (dElm->xTr[1] * dPosI.vtY) + dElm->xTr[2];
	    trYX = dElm->yTr[0];
	    trYYC = (dElm->yTr[1] * dPosI.vtY) + dElm->yTr[2];
	  }
	}
	while((errNum == WLZ_ERR_NONE) && (dPosI.vtX <= mItv->rgtI) &&
	      (dPosI.vtX <= iWSp.rgtpos))
	{
	  sPosD.vtX = (trXX * dPosI.vtX) + trXYC;
	  sPosD.vtY = (trYX * dPosI.vtX) + trYYC;
	  switch(interp)
	  {
	    case WLZ_INTERPOLATION_NEAREST:
	      sPosI.vtX = WLZ_NINT(sPosD.vtX);
	      sPosI.vtY = WLZ_NINT(sPosD.vtY);
	      WlzGreyValueGet(gVWSp, 0, sPosI.vtY, sPosI.vtX);
	      switch(gWSp.pixeltype)
	      {
		case WLZ_GREY_INT:
		  *(dGP.inp)++ = (*(gVWSp->gVal)).inv;
		  break;
		case WLZ_GREY_SHORT:
		  *(dGP.shp)++ = (*(gVWSp->gVal)).shv;
		  break;
		case WLZ_GREY_UBYTE:
		  *(dGP.ubp)++ = (*(gVWSp->gVal)).ubv;
		  break;
		case WLZ_GREY_FLOAT:
		  *(dGP.flp)++ = (*(gVWSp->gVal)).flv;
		  break;
		case WLZ_GREY_DOUBLE:
		  *(dGP.dbp)++ = (*(gVWSp->gVal)).dbv;
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	      break;
	    case WLZ_INTERPOLATION_LINEAR:
	      WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
	      tD0 = sPosD.vtX - floor(sPosD.vtX);
	      tD1 = sPosD.vtY - floor(sPosD.vtY);
	      tD2 = 1.0 - tD0;
	      tD3 = 1.0 - tD1;
	      switch(gWSp.pixeltype)
	      {
		case WLZ_GREY_INT:
		  tD0 = ((gVWSp->gVal[0]).inv * tD2 * tD3) +
			((gVWSp->gVal[1]).inv * tD0 * tD3) +
			((gVWSp->gVal[2]).inv * tD2 * tD1) +
			((gVWSp->gVal[3]).inv * tD0 * tD1);
		  *(dGP.inp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_SHORT:
		  tD0 = ((gVWSp->gVal[0]).shv * tD2 * tD3) +
			((gVWSp->gVal[1]).shv * tD0 * tD3) +
			((gVWSp->gVal[2]).shv * tD2 * tD1) +
			((gVWSp->gVal[3]).shv * tD0 * tD1);
		  *(dGP.shp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_UBYTE:
		  tD0 = ((gVWSp->gVal[0]).ubv * tD2 * tD3) +
			((gVWSp->gVal[1]).ubv * tD0 * tD3) +
			((gVWSp->gVal[2]).ubv * tD2 * tD1) +
			((gVWSp->gVal[3]).ubv * tD0 * tD1);
		  WLZ_CLAMP(tD0, 0.0, 255.0);
		  *(dGP.ubp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_FLOAT:
		  tD0 = ((gVWSp->gVal[0]).flv * tD2 * tD3) +
			((gVWSp->gVal[1]).flv * tD0 * tD3) +
			((gVWSp->gVal[2]).flv * tD2 * tD1) +
			((gVWSp->gVal[3]).flv * tD0 * tD1);
		  *(dGP.flp)++ = tD0;
		  break;
		case WLZ_GREY_DOUBLE:
		  tD0 = ((gVWSp->gVal[0]).dbv * tD2 * tD3) +
			((gVWSp->gVal[1]).dbv * tD0 * tD3) +
			((gVWSp->gVal[2]).dbv * tD2 * tD1) +
			((gVWSp->gVal[3]).dbv * tD0 * tD1);
		  *(dGP.dbp)++ = tD0;
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	      break;
	    case WLZ_INTERPOLATION_CLASSIFY_1:
	      WlzGreyValueGetCon(gVWSp, 0, sPosD.vtY, sPosD.vtX);
	      tD0 = sPosD.vtX - floor(sPosD.vtX);
	      tD1 = sPosD.vtY - floor(sPosD.vtY);
	      switch(gWSp.pixeltype)
	      {
		case WLZ_GREY_INT:
		  for(indx=0; indx < 4; indx++){
		    gTmp[indx] = (gVWSp->gVal[indx]).inv;
		  }
		  tD0 = WlzClassValCon4(gTmp, tD0, tD1);
		  *(dGP.inp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_SHORT:
		  for(indx=0; indx < 4; indx++){
		    gTmp[indx] = (gVWSp->gVal[indx]).shv;
		  }
		  tD0 = WlzClassValCon4(gTmp, tD0, tD1);
		  *(dGP.shp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_UBYTE:
		  for(indx=0; indx < 4; indx++){
		    gTmp[indx] = (gVWSp->gVal[indx]).ubv;
		  }
		  tD0 = WlzClassValCon4(gTmp, tD0, tD1);
		  WLZ_CLAMP(tD0, 0.0, 255.0);
		  *(dGP.ubp)++ = WLZ_NINT(tD0);
		  break;
		case WLZ_GREY_FLOAT:
		  for(indx=0; indx < 4; indx++){
		    gTmp[indx] = (gVWSp->gVal[indx]).flv;
		  }
		  tD0 = WlzClassValCon4(gTmp, tD0, tD1);
		  *(dGP.flp)++ = tD0;
		  break;
		case WLZ_GREY_DOUBLE:
		  for(indx=0; indx < 4; indx++){
		    gTmp[indx] = (gVWSp->gVal[indx]).dbv;
		  }
		  tD0 = WlzClassValCon4(gTmp, tD0, tD1);
		  *(dGP.dbp)++ = tD0;
		  break;
		default:
		  errNum = WLZ_ERR_GREY_TYPE;
		  break;
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_INTERPOLATION_TYPE;
	      break;
	  }
	  ++(dPosI.vtX);
	}
      }
    }
    if(errNum == WLZ_ERR_EOO)           /* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  WlzMeshScanWSpFree(mSnWSp);
  WlzGreyValueFreeWSp(gVWSp);
  return(errNum);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a mesh transform for the given object where the
*               mesh has just two elements which enclose the objects
*               bounding box.  All mesh displacements are zero.
* \param	srcObj			The given object.
* \param	dstBox			Destination bounding box pointer,
*					may be NULL.
* \param	boxDilation		Dilation of box to be sure of enclosing
* 					any possible node.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzMeshTransform 	*WlzMeshFromObjBox(WlzObject *srcObj, WlzIBox2 *dstBox,
					   int boxDilation,
					   WlzErrorNum *dstErr)
{
  WlzMeshNode	*nod;
  WlzMeshElem	*elm;
  WlzMeshTransform *mesh = NULL;
  WlzIBox2	box;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	maxElems = 2,
  		maxNodes = 4;

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
  else
  {
    box = WlzBoundingBox2I(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    box.xMin -= boxDilation;
    box.yMin -= boxDilation;
    box.xMax += boxDilation;
    box.yMax += boxDilation;
    mesh = WlzMeshTransformNew(maxNodes, maxElems, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mesh->nElem = maxElems;
    mesh->nNodes = maxNodes;
    nod = mesh->nodes + 0;
    nod->flags = WLZ_MESH_NODE_FLAGS_BBOX;
    nod->position.vtX = box.xMin;
    nod->position.vtY = box.yMin;
    nod->displacement.vtX = 0.0;
    nod->displacement.vtY = 0.0;
    nod = mesh->nodes + 1;
    nod->flags = WLZ_MESH_NODE_FLAGS_BBOX;
    nod->position.vtX = box.xMax;
    nod->position.vtY = box.yMin;
    nod->displacement.vtX = 0.0;
    nod->displacement.vtY = 0.0;
    nod = mesh->nodes + 2;
    nod->flags = WLZ_MESH_NODE_FLAGS_BBOX;
    nod->position.vtX = box.xMin;
    nod->position.vtY = box.yMax;
    nod->displacement.vtX = 0.0; 
    nod->displacement.vtY = 0.0; 
    nod = mesh->nodes + 3;
    nod->flags = WLZ_MESH_NODE_FLAGS_BBOX;
    nod->position.vtX = box.xMax;
    nod->position.vtY = box.yMax;
    nod->displacement.vtX = 0.0; 
    nod->displacement.vtY = 0.0; 
    elm = mesh->elements + 0;
    elm->idx = 0;
    elm->type = WLZ_MESH_ELEM_TRILINEAR;
    elm->flags = WLZ_MESH_ELEM_FLAGS_NBR_0;
    elm->nodes[0] = 0; elm->nodes[1] = 1; elm->nodes[2] = 2;
    elm->neighbours[0] = 1;
    elm->strainU[0] = elm->strainU[1] = elm->strainU[2] = 0.0;
    elm->strainA[0] = elm->strainA[1] = elm->strainA[2] = 0.0;
    elm = mesh->elements + 1;
    elm->type = WLZ_MESH_ELEM_TRILINEAR;
    elm->idx = 1;
    elm->flags = WLZ_MESH_ELEM_FLAGS_NBR_0;
    elm->nodes[0] = 3; elm->nodes[1] = 2; elm->nodes[2] = 1;
    elm->neighbours[0] = 0;
    elm->strainU[0] = elm->strainU[1] = elm->strainU[2] = 0.0;
    elm->strainA[0] = elm->strainA[1] = elm->strainA[2] = 0.0;
  }
  if(dstBox)
  {
    *dstBox = box;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Removes the bounding box nodes from the mesh.
* \param	mesh			Given mesh transform.
*/
static WlzErrorNum	WlzMeshRemoveObjBox(WlzMeshTransform *mesh)
{
  int		nId,
  		nCnt;
  int		bBoxNodIds[4];
  WlzMeshNode	*nod;
  WlzErrorNum	 errNum = WLZ_ERR_NONE;

  /* Remove bounding box. */
  nId = 0;
  nCnt = 0;
  while((nCnt < 4) && (nCnt < mesh->nNodes))
  {
    nod = mesh->nodes + nId;
    if((nod->flags & WLZ_MESH_NODE_FLAGS_BBOX) != 0)
    {
      bBoxNodIds[nCnt++] = nId;
    }
    ++nId;
  }
#ifdef WLZ_MESH_DEBUG
  errNum = WlzMeshTransformVerify(mesh, 0, NULL, NULL);
#endif /* WLZ_MESH_DEBUG */
  errNum = WlzMeshNodeDelIdx(mesh, 0, bBoxNodIds, nCnt);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds nodes to the mesh for a bounding convex polygon.
*               The given boundary polygon is modified in place. All
*		mesh displacements are zero.
* \param	mesh			Given mesh transform.
* \param	polyObj			Convex boundary polygon.
* \param	minDist			Minimum distance between mesh nodes.
* \param	maxDist			Maximum distance between mesh nodes.
*/
static WlzErrorNum WlzMeshBoundCvPolyAdd(WlzMeshTransform *mesh,
				         WlzObject *polyObj,
				         double minDist, double maxDist)
{
  WlzDVertex2	scaleVx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((errNum = WlzMeshBoundPolyFix(polyObj, minDist)) == WLZ_ERR_NONE)
  {
    scaleVx.vtX = scaleVx.vtY = 1.0;
    errNum = WlzMeshDomainAdd(mesh, polyObj, minDist, scaleVx);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Fixes a boundary polygon so that it's minimum length
*		side is >= the given miniimum distance.
* \param	polyObj			Convex boundary polygon.
* \param	minDist			Minimum distance between mesh nodes.
*/
static WlzErrorNum WlzMeshBoundPolyFix(WlzObject *polyObj, double minDist)
{
  int		tstId,
		cntOut,
		distFlg,
		loopFlg,
		nVxCnt;
  double	tD0,
		minDistSq,
		segDist,
  		segDistSq;
  WlzDVertex2	newVx,
  		segVx;
  WlzIVertex2	*curVxP,
    		*tstVxP;
  WlzPolygonDomain *poly;
  WlzMeshPolyVx	*curVxElm,
		*cvxVxElm,
		*fstVxElm,
		*tstVxElm,
  		*nxtVxElm,
		*prvVxElm,
		*polyVxPool = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(polyObj->type != WLZ_2D_POLYGON)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((poly = polyObj->domain.poly) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(poly->type != WLZ_POLYGON_INT)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    /* Remove any wrap around and allocate a pool of linked list elements. */
    nVxCnt = poly->nvertices;
    curVxP = poly->vtx;
    tstVxP = poly->vtx + nVxCnt - 1;
    while((nVxCnt > 0) &&
          (curVxP->vtX == tstVxP->vtX) && (curVxP->vtY == tstVxP->vtY))
    {
      --nVxCnt;
      --tstVxP;
    }
    if(nVxCnt < 3)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if((polyVxPool = (WlzMeshPolyVx *)
    		          AlcMalloc(nVxCnt * sizeof(WlzMeshPolyVx))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    minDistSq = minDist * minDist;
    /* Copy the given polygon's vertices into a linked list data 
     * structure. */
    tstId = 0;
    prvVxElm = NULL;
    fstVxElm = polyVxPool;
    curVxElm = fstVxElm;
    curVxP = poly->vtx;
    while(tstId < nVxCnt)
    {
      curVxElm->id = tstId++;
      curVxElm->vx.vtX = curVxP->vtX;
      curVxElm->vx.vtY = curVxP->vtY;
      ++curVxP;
      curVxElm->prev = prvVxElm;
      prvVxElm = curVxElm;
      prvVxElm->next = ++curVxElm;
    }
    nxtVxElm = polyVxPool + nVxCnt - 1;
    fstVxElm->prev = nxtVxElm;
    nxtVxElm->next = fstVxElm;
    /* Do a modified Jarvis March around the polygon until all vertices
     * have been visited and the distance between the vertices is greater
     * than the given minimum distance. */
    loopFlg = 0;
    curVxElm = fstVxElm;
    do
    {
      /* Find next vertex on the convex hull of the polygon. */
      cvxVxElm = curVxElm->next;
      nxtVxElm = curVxElm->prev;
      tstVxElm = cvxVxElm->next;
      cntOut = nVxCnt;
      while((--cntOut > 0) && (tstVxElm->id != nxtVxElm->id))
      {
	/* Does the test vertex lie to the left or right of the line segment
	 * formed by the current and current convex hull vertex? */
        if((tstVxElm->id != curVxElm->id) && (tstVxElm->id != cvxVxElm->id))
	{
	  tD0 = WlzGeomTriangleSnArea2(curVxElm->vx, cvxVxElm->vx,
	  			       tstVxElm->vx);
	  if(tD0 < -(DBL_EPSILON))
	  {
	    cvxVxElm = tstVxElm;
	  }
	}
	tstVxElm = tstVxElm->next;
      }
      /* Find distance from current to new convex hull vertex. */
      segVx.vtX = cvxVxElm->vx.vtX - curVxElm->vx.vtX; 
      segVx.vtY = cvxVxElm->vx.vtY - curVxElm->vx.vtY; 
      segDistSq = (segVx.vtX * segVx.vtX) + (segVx.vtY * segVx.vtY);
      /* If distance is less that the minimum distance then move the
       * vertex along the segment from it to the current vertex until
       * it is minDist away. */
      distFlg = segDistSq < minDistSq;
      if(distFlg)
      {
        /* Compute new vertex position. */
	newVx = cvxVxElm->vx;
	if((segVx.vtX >= -(DBL_EPSILON)) && (segVx.vtX <= DBL_EPSILON))
	{
	  newVx.vtY = (segVx.vtY > 0.0)?
		      curVxElm->vx.vtY + minDist:
		      curVxElm->vx.vtY - minDist;
	}
	else if((segVx.vtY >= -(DBL_EPSILON)) && (segVx.vtY <= DBL_EPSILON))
	{
	  newVx.vtX = (segVx.vtX > 0.0)?
		      curVxElm->vx.vtX + minDist:
		      curVxElm->vx.vtX - minDist;
	}
	else
	{
	  segDist = sqrt(segDistSq);
	  tD0 = (minDist / segDist);
	  newVx.vtX = curVxElm->vx.vtX +
	  	      (tD0 * (newVx.vtX - curVxElm->vx.vtX));
	  newVx.vtY = curVxElm->vx.vtY +
	  	      (tD0 * (newVx.vtY - curVxElm->vx.vtY));
	}
	cvxVxElm->vx = newVx;
      }
      curVxElm->next = cvxVxElm;
      cvxVxElm->prev = curVxElm;
      /* Check for loop around polygon. */
      if(loopFlg || (cvxVxElm->id < curVxElm->id))
      {
        ++loopFlg;
      }
      curVxElm = cvxVxElm;
    } while(((cntOut > 0) && (loopFlg < 3)) || distFlg);
    if(cntOut > 0)
    {
      /* Copy vertices back into the polygon. */
      tstId = 0;
      curVxP = poly->vtx;
      fstVxElm = curVxElm;
      do
      {
	curVxP->vtX = curVxElm->vx.vtX;
	curVxP->vtY = curVxElm->vx.vtY;
	++curVxP;
	curVxElm = curVxElm->next;
	++tstId;
      }
      while(curVxElm != fstVxElm);
      poly->nvertices = tstId;
    }
    else
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(polyVxPool)
  {
    AlcFree(polyVxPool);
  }
  return(errNum);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a mesh transform for the given object where the
*               mesh has a constant node separation. All mesh displacements
*               are zero.
* \param	srcObj			The given object.
* \param	lDist			Mesh vertex inter-line and 
*					inter-column distance.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzMeshTransform *WlzMeshFromObjBlock(WlzObject *srcObj,
					     unsigned int lDist,
					     WlzErrorNum *dstErr)
{
  int		tI0;
  unsigned int	lIdx,
  		nLines,
		maxNodes,
		maxElems;
  WlzInterval	*itv0,
  		*lnItvs = NULL;
  WlzObject	*tObj0,
  		*cvPolyObj = NULL,
  		*cvObj = NULL,
		*dCvObj = NULL;
  WlzMeshTransform *mesh = NULL;
  WlzIVertex2	org;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
  else if(lDist < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    cvPolyObj = WlzObjToConvexPolygon(srcObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cvObj = WlzPolyToObj(cvPolyObj->domain.poly, WLZ_SIMPLE_FILL,
    			 &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj0 = WlzAssignObject(
    	    WlzMakeCircleObject(lDist, 0.0, 0.0, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      dCvObj= WlzAssignObject(
      	      WlzStructDilation(srcObj, tObj0, &errNum), NULL);
    }
    if(tObj0)
    {
      WlzFreeObj(tObj0);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Calculate interline spacing from vertex distance. */
    org.vtX = dCvObj->domain.i->kol1;
    org.vtY = dCvObj->domain.i->line1;
    nLines = ((dCvObj->domain.i->lastln - dCvObj->domain.i->line1 + lDist) /
    	     lDist) + 1;
    /* Allocate space for line start end points using intervals */
    if((lnItvs = (WlzInterval *)AlcMalloc(sizeof(WlzInterval) *
    					  nLines)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzMeshLineCvExtrema(lnItvs, dCvObj, lDist, nLines);
  }
  if(dCvObj)
  {
    WlzFreeObj(dCvObj);
  }
  if(cvPolyObj)
  {
    WlzFreeObj(cvPolyObj);
  }
  if(cvObj)
  {
    WlzFreeObj(cvObj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Calculate the number of vertex nodes and triangular elements */
    lIdx = 0;
    itv0 = lnItvs;
    maxNodes = 0;
    maxElems = 0;
    while(lIdx < nLines)
    {
      tI0 = (itv0->iright - itv0->ileft) / lDist;
      maxElems += tI0;
      maxNodes += tI0 + 1;
      ++itv0;
      ++lIdx;
    }
    maxElems *= 2; 		      /* Two triangles per rectangular block */
    mesh = WlzMeshTransformNew(maxNodes, maxElems, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzMeshTransFillBlock(mesh, lnItvs, org, lDist, nLines);
  }
  if(lnItvs)
  {
    AlcFree(lnItvs);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	New mesh transform.
* \ingroup	WlzTransform
* \brief	Creates a mesh transform for the given object where the
*               mesh vertices are placed with a higher density where the given
*               object's gradient is highest. All mesh displacements are zero.
* \param	srcObj			The given object.
* \param	minDist			Minimum distance between mesh vertices.
* \param	maxDist			Maximum distance between mesh vertices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzMeshTransform *WlzMeshFromObjGrad(WlzObject *srcObj,
				            unsigned int minDist,
				            unsigned int maxDist,
				     	    WlzErrorNum *dstErr)
{
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
  		*cvPolyObj = NULL,
  		*grdObj = NULL,
		*ssObjMin = NULL,
		*ssObjMax = NULL,
		*thrObjMin = NULL,
		*thrObjMax = NULL;
  WlzMeshTransform *mesh = NULL;
  WlzIVertex3	ssSzMin,
  		ssSzMax;
  WlzDVertex2	scale;
  WlzPixelV	minV,
  		maxV,
		thrV;
  WlzIBox2	box;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(minDist < 2)
  {
    minDist = 2;
  }
  if(maxDist < minDist)
  {
    maxDist = minDist;
  }
  /* Compute boundary. */
  if(((mesh = WlzMeshFromObjBox(srcObj, &box, maxDist * 4,
  				&errNum)) != NULL) &&
     (errNum == WLZ_ERR_NONE))
  {
    /* Compute dilated convex hull. */
    tObj0 = WlzAssignObject(
    	    WlzMakeCircleObject(minDist, 0.0, 0.0, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      tObj1= WlzAssignObject(
      	     WlzStructDilation(srcObj, tObj0, &errNum), NULL);
    }
    if(tObj0)
    {
      WlzFreeObj(tObj0);
      tObj0 = NULL;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      cvPolyObj = WlzAssignObject(
      		  WlzObjToConvexPolygon(tObj1, &errNum), NULL);
    }
    if(tObj1)
    {
      WlzFreeObj(tObj1);
      tObj1 = NULL;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Add nodes along the original object's convex hull. */
      if(cvPolyObj && (cvPolyObj->type == WLZ_2D_POLYGON))
      {
	errNum = WlzMeshBoundCvPolyAdd(mesh, cvPolyObj, minDist, maxDist);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzMeshRemoveObjBox(mesh);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzMeshSqueeze(mesh);
    }
    if(cvPolyObj)
    {
      WlzFreeObj(cvPolyObj);
    }
  }
  if((errNum == WLZ_ERR_NONE) && (srcObj->values.core != NULL))
  {
    /* Fill boundary triangulation with nodes depending on image. */
    /* Compute a gradient image. */
    grdObj = WlzAssignObject(
    	     WlzLaplacian(srcObj, 3, 1, 1, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      /* Get the grey range, set background and compute threshold value. */
      errNum = WlzGreyRange(grdObj, &minV, &maxV);
      if(errNum == WLZ_ERR_NONE)
      {
	(void )WlzValueConvertPixel(&minV, minV, WLZ_GREY_INT);
	(void )WlzValueConvertPixel(&maxV, maxV, WLZ_GREY_INT);
	thrV.type = WLZ_GREY_INT;
	thrV.v.inv = minV.v.inv + ((maxV.v.inv - minV.v.inv) / 4);
	errNum = WlzSetBackground(grdObj, minV);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      ssSzMin.vtX = ssSzMin.vtY = (minDist / 5) + 1;
      if(ssSzMin.vtX > 1.0)
      {
	ssObjMin = WlzAssignObject(
		   WlzSampleObj(grdObj, ssSzMin, WLZ_SAMPLEFN_MIN,
				&errNum), NULL);
      }
      else
      {
	ssObjMin = WlzAssignObject(grdObj, NULL);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      ssSzMax.vtX = ssSzMax.vtY = (maxDist / 5) + 1;
      if(ssSzMax.vtX > 1.0)
      {
	ssObjMax = WlzAssignObject(
		   WlzSampleObj(grdObj, ssSzMax, WLZ_SAMPLEFN_MIN,
				&errNum), NULL);
      }
      else
      {
	ssObjMax = WlzAssignObject(grdObj, NULL);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      thrObjMin = WlzAssignObject(
	      WlzThreshold(ssObjMin, thrV, WLZ_THRESH_LOW, &errNum), NULL);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      thrObjMax = WlzAssignObject(
		  WlzThreshold(ssObjMax, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
    }
    if(grdObj)
    {
      WlzFreeObj(grdObj);
    }
    if(ssObjMin)
    {
      WlzFreeObj(ssObjMin);
    }
    if(ssObjMax)
    {
      WlzFreeObj(ssObjMax);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Add nodes in regions of high gradient to the mesh. */
      if(thrObjMin && (thrObjMin->type == WLZ_2D_DOMAINOBJ))
      {
	scale.vtX = ssSzMin.vtX;
	scale.vtY = ssSzMin.vtY;
	errNum = WlzMeshDomainAdd(mesh, thrObjMin, minDist, scale);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Add nodes in regions of low gradient to the mesh. */
      if(thrObjMax && (thrObjMax->type == WLZ_2D_DOMAINOBJ))
      {
	scale.vtX = ssSzMax.vtX;
	scale.vtY = ssSzMax.vtY;
	errNum = WlzMeshDomainAdd(mesh, thrObjMax, maxDist, scale);
      }
    }
    if(thrObjMin)
    {
      WlzFreeObj(thrObjMin);
    }
    if(thrObjMax)
    {
      WlzFreeObj(thrObjMax);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mesh);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Calculates the block mesh line extrema for the given (convex)
* 		object.
* \param	lnItvs			Vector of intervals to be set to block
* 					line extrema.
* \param	cvObj			The given object.
* \param	lDist			Distance between block lines.
* \param	nLines			Number of block lines.
*/
static WlzErrorNum WlzMeshLineCvExtrema(WlzInterval *lnItvs,
					WlzObject *cvObj,
				        unsigned int lDist,
					unsigned int nLines)
{
  int		tI0;
  unsigned int	lIdx,
  		lLftIdx,
		lRgtIdx;
  WlzInterval	*itv0,
		*itv1,
  		*itv2;
  WlzIVertex2	org;
  WlzIntervalWSpace iWSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  org.vtX = cvObj->domain.i->kol1;
  org.vtY = cvObj->domain.i->line1;
  /* Preset impossible line start and end points for faster scanning */
  lIdx = 0;
  itv0 = lnItvs;
  lnItvs->ileft = cvObj->domain.i->lastkl;
  lnItvs->iright = cvObj->domain.i->kol1;
  while(++lIdx < nLines)
  {
    *++itv0 = *lnItvs;
  }
  /* Scan through the convex object building line start and end points */
  if((errNum = WlzInitRasterScan(cvObj, &iWSp,
				 WLZ_RASTERDIR_ILIC)) == WLZ_ERR_NONE)
  {
    while((errNum == WLZ_ERR_NONE) &&
	  ((errNum = WlzNextInterval(&iWSp)) == WLZ_ERR_NONE))
    {
      lIdx = (iWSp.linpos - org.vtY) / lDist;
      itv0 = (lnItvs + lIdx);
      if(itv0->ileft > iWSp.lftpos)
      {
	itv0->ileft = iWSp.lftpos;
      }
      if(itv0->iright < iWSp.rgtpos)
      {
	itv0->iright = iWSp.rgtpos;
      }
    }
    if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make all intervals relative to the domains first column and
       quantize using the inter-line distance.
       Find indicies of the lines with left/right extrema */
    lIdx = lLftIdx = lRgtIdx = 0;
    itv0 = itv1 = itv2 = lnItvs;
    while(lIdx < nLines)
    {
      tI0 = itv0->ileft - org.vtX;
      itv0->ileft = tI0 - (tI0 % lDist);
      tI0 = itv0->iright - org.vtX;
      itv0->iright = tI0 + (lDist - (tI0 % lDist));
      if(itv0->ileft < itv1->ileft)
      {
        itv1 = itv0;
	lLftIdx = lIdx;
      }
      if(itv0->iright > itv2->iright)
      {
        itv2 = itv0;
	lRgtIdx = lIdx;
      }
      ++lIdx;
      ++itv0;
    }
    /* Make last line of intervals same as last - 1 */
    *(lnItvs + nLines - 1) = *(lnItvs + nLines - 2);
    /* Propagate left/right extrema from their maximal lines */
    lIdx = nLines - 1;
    itv0 = lnItvs + nLines - 1;
    itv1 = lnItvs + nLines - 2;
    while(lIdx-- > lLftIdx)
    {
      itv0--->ileft = itv1--->ileft;
    }
    lIdx = nLines - 1;
    itv0 = lnItvs + nLines - 1;
    itv1 = lnItvs + nLines - 2;
    while(lIdx-- > lRgtIdx)
    {
      itv0--->iright = itv1--->iright;
    }
  }
/* #define WLZ_MESH_DEBUG */
#ifdef WLZ_MESH_DEBUG
  lIdx = 0;
  itv0 = lnItvs;
  while(lIdx < nLines)
  {
    fprintf(stderr, "%4d %4d %4d\n", lIdx, itv0->ileft, itv0->iright);
    ++lIdx;
    ++itv0;
  }
  printf("\n\n");
#endif /* WLZ_MESH_DEBUG */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Calculates the mesh transform (already allocated) using the
* 		given block mesh line interval extrema.
* \param	mesh			Mesh transform to fill.
* \param	lnItvs			Vector of block line intervals.
* \param	org			Origin used by intervals.
* \param	lDist			Distance between block lines.
* \param	nLines			Number of block lines.
*/
static WlzErrorNum WlzMeshTransFillBlock(WlzMeshTransform *mesh,
				         WlzInterval *lnItvs,
					 WlzIVertex2 org,
					 unsigned int lDist,
					 unsigned int nLines)
{
  int		tI0,
  		lIdx,
		nNod0,
		nNod1,
		elmIdx,
		elmIdxLo,
		elmIdxHi,
		elmIdxHiPrev,
  		nodIdx0,
		nodIdx1,
		nodIdxLo,
		nodIdxLoLft,
		nodIdxLoRgt,
		nodIdxHi,
		nodIdxHiLft,
		nodIdxHiRgt;
  double	tD0;
  WlzMeshElem	*elm,
  		*elmLo,
  		*elmHi,
		*elmHiPrev;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  lIdx = 1;
  elmIdx = 0;
  elmIdxHiPrev = -1;
  elm = mesh->elements;
  nodIdx1 = 0;
  nNod1 = WlzMeshTransFillBlockLnNod(mesh->nodes, lnItvs, org, lDist);
  mesh->nNodes = nNod1;
  while(lIdx < nLines)
  {
    nodIdx0 = nodIdx1;
    nNod0 = nNod1;
    nodIdx1 = nodIdx0 + nNod0;
    org.vtY += lDist;
    nNod1 = WlzMeshTransFillBlockLnNod(mesh->nodes + nodIdx1, lnItvs + lIdx,
    				       org, lDist);
    mesh->nNodes += nNod1;
    nodIdxLoLft = nodIdx0;
    nodIdxHiLft = nodIdx1;
    tD0 = (mesh->nodes + nodIdxHiLft)->position.vtX -
          (mesh->nodes + nodIdxLoLft)->position.vtX;
    if(tD0 > 0.5)
    {
      tI0 = (tD0 + 0.5) / lDist;
      nodIdxLoLft += tI0;
    }
    else if(tD0 < -0.5)
    {
      tI0 = (-tD0 + 0.5) / lDist;
      nodIdxHiLft += tI0;
    }
    nodIdxLoRgt = nodIdx0 + nNod0 - 1;
    nodIdxHiRgt = nodIdx1 + nNod1 - 1;
    tD0 = (mesh->nodes + nodIdxHiRgt)->position.vtX -
          (mesh->nodes + nodIdxLoRgt)->position.vtX;
    if(tD0 > 0.5)
    {
      tI0 = (tD0 + 0.5) / lDist;
      nodIdxHiRgt -= tI0;
    }
    else if(tD0 < -0.5)
    {
      tI0 = (tD0 - 0.5) / lDist;
      nodIdxLoRgt += tI0;
    }
    nodIdxLo = nodIdxLoLft;
    nodIdxHi = nodIdxHiLft;
    if(elmIdxHiPrev > 0)
    {
      elmHiPrev = mesh->elements + elmIdxHiPrev;
      if((tI0 = nodIdxLoLft - elmHiPrev->nodes[1]) > 0)
      {
	elmIdxHiPrev += 2 * tI0;
      }
    }
    while(nodIdxLo < nodIdxLoRgt)
    {
      elmLo = elm++;
      elmHi = elm++;
      elmIdxLo = elmIdx++;
      elmIdxHi = elmIdx++;
      elmLo->type = WLZ_MESH_ELEM_TRILINEAR;
      elmLo->idx = elmIdxLo;
      elmLo->nodes[0] = nodIdxLo;
      elmLo->nodes[1] = nodIdxLo + 1;
      elmLo->nodes[2] = nodIdxHi;
      elmLo->flags |= WLZ_MESH_ELEM_FLAGS_NBR_0;
      elmLo->neighbours[0] = elmIdxHi;
      if(nodIdxLo != nodIdxLoLft)
      {
        elmLo->flags |= WLZ_MESH_ELEM_FLAGS_NBR_1;
	elmLo->neighbours[1] = elmIdxLo - 1;
	(elmLo - 1)->flags |= WLZ_MESH_ELEM_FLAGS_NBR_1;
	(elmLo - 1)->neighbours[1] = elmIdxLo;
      }
      if(elmIdxHiPrev > 0)
      {
	elmHiPrev = mesh->elements + elmIdxHiPrev;
	if((elmHiPrev->nodes[1] == nodIdxLo) &&
	   (elmHiPrev->nodes[1] < nodIdxLoRgt))
	{
	  elmLo->flags |= WLZ_MESH_ELEM_FLAGS_NBR_2;
	  elmLo->neighbours[2] = elmIdxHiPrev;
	  elmHiPrev->flags |= WLZ_MESH_ELEM_FLAGS_NBR_2;
	  elmHiPrev->neighbours[2]  = elmIdxLo;
	  elmIdxHiPrev += 2;
	}
      }
      elmHi->type = WLZ_MESH_ELEM_TRILINEAR;
      elmHi->idx = elmIdxHi;
      elmHi->nodes[0] = nodIdxHi + 1;
      elmHi->nodes[1] = nodIdxHi;
      elmHi->nodes[2] = nodIdxLo + 1;
      elmHi->flags |= WLZ_MESH_ELEM_FLAGS_NBR_0;
      elmHi->neighbours[0] = elmIdxLo;
      ++nodIdxLo;
      ++nodIdxHi;
    }
    elmIdxHiPrev = elmIdx - (2 * (nodIdxLoRgt - nodIdxLoLft)) + 1;
    ++lIdx;
  }
  mesh->nElem = elmIdx;
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Calculates the mesh transform node vertices for a single
*		line between extrema of the interval.
* \param	nod			Ptr to first node available.
* \param	lnItv			The line interval.
* \param	org			Origin used by this interval.
* \param	lDist			Distance between blocks.
*/
static unsigned int WlzMeshTransFillBlockLnNod(WlzMeshNode *nod,
					       WlzInterval *lnItv,
					       WlzIVertex2 org,
					       unsigned int lDist)
{
  WlzIVertex2	vx0,
  		vx1;
  unsigned int	nNod = 0;

  vx0.vtX = org.vtX + lnItv->ileft;
  vx1.vtX = org.vtX + lnItv->iright;
  vx0.vtY = vx1.vtY = org.vtY;
  while(vx0.vtX <= vx1.vtX)
  {
    nod->flags = WLZ_MESH_NODE_FLAGS_BLOCK;
    nod->position.vtX = vx0.vtX;
    nod->position.vtY = vx0.vtY;
    nod->displacement.vtX = 0.0;
    nod->displacement.vtY = 0.0;
    vx0.vtX += lDist;
    ++nNod;
    ++nod;
  }
  return(nNod);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given integer vertex
*               vector in place and using the given mesh transform. It is an
*               error if any vertex is not in the mesh or if the signed area
*               of any mesh element is <= zero.
* \param	mesh			Given mesh transform.
* \param	vxVec			Given integer vertex vector.
* \param	vxCount			Number of vertices.
*/
static WlzErrorNum WlzMeshTransformVxVecI(WlzMeshTransform *mesh,
					  WlzIVertex2 *vxVec,
					  int vxCount)
{
  int		elmIdx,
  		neighbourId,
  		trValid;
  unsigned int	neighbourMask;
  double	elmArea2;
  WlzMeshElem	*elm;
  WlzDVertex2	vxD,
		vxP;
  double	pArea2[3],
  		xTr[3], /* Affine transform for destination triangle:
		           xd = xTr[0] * xs + xTr[1] * ys + xTr[2] */
		yTr[3]; /* yd = yTr[0] * xs + yTr[1] * ys + yTr[2] */
  WlzMeshNode	*nod[3];
  WlzDVertex2	dspVx[3],
  		elmVx[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmIdx = 0;
  trValid = 0;
  elm = mesh->elements;
  vxP.vtX = vxVec->vtX;
  vxP.vtY = vxVec->vtY;
  nod[0] = mesh->nodes + elm->nodes[0];
  nod[1] = mesh->nodes + elm->nodes[1];
  nod[2] = mesh->nodes + elm->nodes[2];
  elmVx[0] = nod[0]->position;
  elmVx[1] = nod[1]->position;
  elmVx[2] = nod[2]->position;
  if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
  					elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  while((vxCount > 0) && (errNum == WLZ_ERR_NONE))
  {
    /* Find neighbour which is in the direction of the vertex,
       if none then the vertex is contained within this element */
    neighbourMask = WLZ_MESH_ELEM_FLAGS_NONE;
    if((pArea2[0] = WlzGeomTriangleSnArea2(elmVx[1], elmVx[2],
					   vxP)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 0;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_0;
    }
    else if((pArea2[1] = WlzGeomTriangleSnArea2(elmVx[2], elmVx[0],
					       vxP)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 1;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_1;
    }
    else if((pArea2[2] = elmArea2 - pArea2[0] -
    			 pArea2[1]) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 2;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_2;
    }
    if(neighbourMask == WLZ_MESH_ELEM_FLAGS_NONE)
    {
      /* Compute new affine transform for interpolation from the source
         to the displaced element if required. */
      if(trValid == 0)
      {
	vxD = nod[0]->displacement;
        (dspVx + 0)->vtX = vxD.vtX + (elmVx + 0)->vtX;
        (dspVx + 0)->vtY = vxD.vtY + (elmVx + 0)->vtY;
	vxD = nod[1]->displacement;
        (dspVx + 1)->vtX = vxD.vtX + (elmVx + 1)->vtX;
        (dspVx + 1)->vtY = vxD.vtY + (elmVx + 1)->vtY;
	vxD = nod[2]->displacement;
        (dspVx + 2)->vtX = vxD.vtX + (elmVx + 2)->vtX;
        (dspVx + 2)->vtY = vxD.vtY + (elmVx + 2)->vtY;
	WlzMeshAfTrSolve(xTr, yTr, elmArea2, elmVx, dspVx);
        trValid = 1;
      }
      /* Vertex is in this element interpolate new displaced vertex using
      the affine transform. */
      vxD.vtX = (xTr[0] * vxP.vtX) + (xTr[1] * vxP.vtY) + xTr[2];
      vxD.vtY = (yTr[0] * vxP.vtX) + (yTr[1] * vxP.vtY) + yTr[2];
      vxVec->vtX = WLZ_NINT(vxD.vtX);
      vxVec->vtY = WLZ_NINT(vxD.vtY);
      if(--vxCount > 0)
      {
        ++vxVec;
	vxP.vtX = vxVec->vtX;
	vxP.vtY = vxVec->vtY;
      }
    }
    else
    {
      /* Vertex is NOT in this element so walk towards the element
	 which does contain it. */
      if(elm->flags & neighbourMask)
      {
	trValid = 0;
	elmIdx = elm->neighbours[neighbourId];
	elm = mesh->elements + elmIdx;
	nod[0] = mesh->nodes + elm->nodes[0];
	nod[1] = mesh->nodes + elm->nodes[1];
	nod[2] = mesh->nodes + elm->nodes[2];
	elmVx[0] = nod[0]->position;
	elmVx[1] = nod[1]->position;
	elmVx[2] = nod[2]->position;
	if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
					    elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
      else
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex vector, in
* 		place and using the given mesh transform.  It is an error if
* 		any vertex is not in the mesh or if the signed area of any
* 		mesh element is <= zero.
* \param	mesh			Given mesh transform.
* \param	vxVec			Given double vertex vector.
* \param	vxCount			Number of vertices.
*/
static WlzErrorNum WlzMeshTransformVxVecF(WlzMeshTransform *mesh,
					  WlzFVertex2 *vxVec,
					  int vxCount)
{
  int		elmIdx,
  		neighbourId,
  		trValid;
  unsigned int	neighbourMask;
  double	elmArea2;
  WlzMeshElem	*elm;
  WlzDVertex2	vxD;
  double	pArea2[3],
  		xTr[3], /* Affine transform for destination triangle:
		           xd = xTr[0] * xs + xTr[1] * ys + xTr[2] */
		yTr[3]; /* yd = yTr[0] * xs + yTr[1] * ys + yTr[2] */
  WlzMeshNode	*nod[3];
  WlzDVertex2	dspVx[3],
  		elmVx[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmIdx = 0;
  trValid = 0;
  elm = mesh->elements;
  nod[0] = mesh->nodes + elm->nodes[0];
  nod[1] = mesh->nodes + elm->nodes[1];
  nod[2] = mesh->nodes + elm->nodes[2];
  elmVx[0] = nod[0]->position;
  elmVx[1] = nod[1]->position;
  elmVx[2] = nod[2]->position;
  if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
  					elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  while((vxCount > 0) && (errNum == WLZ_ERR_NONE))
  {
    /* Find neighbour which is in the direction of the vertex,
       if none then the vertex is contained within this element */
    vxD.vtX = vxVec->vtX;
    vxD.vtY = vxVec->vtY;
    neighbourMask = WLZ_MESH_ELEM_FLAGS_NONE;
    if((pArea2[0] = WlzGeomTriangleSnArea2(elmVx[1], elmVx[2],
					   vxD)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 0;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_0;
    }
    else if((pArea2[1] = WlzGeomTriangleSnArea2(elmVx[2], elmVx[0],
					      vxD)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 1;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_1;
    }
    else if((pArea2[2] = elmArea2 - pArea2[0] -
    			 pArea2[1]) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 2;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_2;
    }
    if(neighbourMask == WLZ_MESH_ELEM_FLAGS_NONE)
    {
      /* Compute new affine transform for interpolation from the source
         to the displaced element if required. */
      if(trValid == 0)
      {
	vxD = nod[0]->displacement;
        (dspVx + 0)->vtX = vxD.vtX + (elmVx + 0)->vtX;
        (dspVx + 0)->vtY = vxD.vtY + (elmVx + 0)->vtY;
	vxD = nod[1]->displacement;
        (dspVx + 1)->vtX = vxD.vtX + (elmVx + 1)->vtX;
        (dspVx + 1)->vtY = vxD.vtY + (elmVx + 1)->vtY;
	vxD = nod[2]->displacement;
        (dspVx + 2)->vtX = vxD.vtX + (elmVx + 2)->vtX;
        (dspVx + 2)->vtY = vxD.vtY + (elmVx + 2)->vtY;
	WlzMeshAfTrSolve(xTr, yTr, elmArea2, elmVx, dspVx);
        trValid = 1;
      }
      /* Vertex is in this element interpolate new displaced vertex using
      the affine transform. */
      vxVec->vtX = (xTr[0] * vxVec->vtX) + (xTr[1] * vxVec->vtY) + xTr[2];
      vxVec->vtY = (yTr[0] * vxVec->vtX) + (yTr[1] * vxVec->vtY) + yTr[2];
      if(vxCount-- > 0)
      {
        ++vxVec;
      }
    }
    else
    {
      /* Vertex is NOT in this element so walk towards the element
	 which does contain it. */
      if(elm->flags & neighbourMask)
      {
	trValid = 0;
	elmIdx = elm->neighbours[neighbourId];
	elm = mesh->elements + elmIdx;
	nod[0] = mesh->nodes + elm->nodes[0];
	nod[1] = mesh->nodes + elm->nodes[1];
	nod[2] = mesh->nodes + elm->nodes[2];
	elmVx[0] = nod[0]->position;
	elmVx[1] = nod[1]->position;
	elmVx[2] = nod[2]->position;
	if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
					    elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
      else
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Ttransformed vertex.
* \ingroup	WlzTransform
* \brief	Transform the vertex using the given mesh transform.
* \param	vtx			Given double vertex.
* \param	mesh			Given mesh transform.
* \param	dstErr			Error return.
*/
WlzDVertex2 WlzMeshTransformVtx(WlzDVertex2 vtx, WlzMeshTransform *mesh,
  				WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDVertex2	rtnVtx;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh->nElem < 0) || (mesh->nNodes < 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else {
    rtnVtx = vtx;
    errNum = WlzMeshTransformVxVecD(mesh, &rtnVtx, 1);
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rtnVtx);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Transforms the vertices in the given double vertex vector,
*		in place and using the given mesh transform. It is an error if
*		any vertex is not in the mesh or if the signed area of any mesh
*		element is <= zero.
* \param	mesh			Given mesh transform.
* \param	vxVec			Given double vertex vector.
* \param	vxCount			Number of vertices.
*/
static WlzErrorNum WlzMeshTransformVxVecD(WlzMeshTransform *mesh,
					  WlzDVertex2 *vxVec,
					  int vxCount)
{
  int		elmIdx,
  		neighbourId,
  		trValid;
  unsigned int	neighbourMask;
  double	elmArea2;
  WlzMeshElem	*elm;
  WlzDVertex2	vxD;
  double	pArea2[3],
  		xTr[3], /* Affine transform for destination triangle:
		           xd = xTr[0] * xs + xTr[1] * ys + xTr[2] */
		yTr[3]; /* yd = yTr[0] * xs + yTr[1] * ys + yTr[2] */
  WlzMeshNode	*nod[3];
  WlzDVertex2	dspVx[3],
  		elmVx[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elmIdx = 0;
  trValid = 0;
  elm = mesh->elements;
  nod[0] = mesh->nodes + elm->nodes[0];
  nod[1] = mesh->nodes + elm->nodes[1];
  nod[2] = mesh->nodes + elm->nodes[2];
  elmVx[0] = nod[0]->position;
  elmVx[1] = nod[1]->position;
  elmVx[2] = nod[2]->position;
  if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
  					elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  while((vxCount > 0) && (errNum == WLZ_ERR_NONE))
  {
    /* Find neighbour which is in the direction of the vertex,
       if none then the vertex is contained within this element */
    neighbourMask = WLZ_MESH_ELEM_FLAGS_NONE;
    if((pArea2[0] = WlzGeomTriangleSnArea2(elmVx[1], elmVx[2],
					   *vxVec)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 0;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_0;
    }
    else if((pArea2[1] = WlzGeomTriangleSnArea2(elmVx[2], elmVx[0],
					     *vxVec)) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 1;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_1;
    }
    else if((pArea2[2] = elmArea2 - pArea2[0] -
    			 pArea2[1]) < -(WLZ_MESH_TOLERANCE_SQ))
    {
      neighbourId = 2;
      neighbourMask = WLZ_MESH_ELEM_FLAGS_NBR_2;
    }
    if(neighbourMask == WLZ_MESH_ELEM_FLAGS_NONE)
    {
      /* Compute new affine transform for interpolation from the source
         to the displaced element if required. */
      if(trValid == 0)
      {
	vxD = nod[0]->displacement;
        (dspVx + 0)->vtX = vxD.vtX + (elmVx + 0)->vtX;
        (dspVx + 0)->vtY = vxD.vtY + (elmVx + 0)->vtY;
	vxD = nod[1]->displacement;
        (dspVx + 1)->vtX = vxD.vtX + (elmVx + 1)->vtX;
        (dspVx + 1)->vtY = vxD.vtY + (elmVx + 1)->vtY;
	vxD = nod[2]->displacement;
        (dspVx + 2)->vtX = vxD.vtX + (elmVx + 2)->vtX;
        (dspVx + 2)->vtY = vxD.vtY + (elmVx + 2)->vtY;
	WlzMeshAfTrSolve(xTr, yTr, elmArea2, elmVx, dspVx);
        trValid = 1;
      }
      /* Vertex is in this element interpolate new displaced vertex using
      the affine transform. */
      vxVec->vtX = (xTr[0] * vxVec->vtX) + (xTr[1] * vxVec->vtY) + xTr[2];
      vxVec->vtY = (yTr[0] * vxVec->vtX) + (yTr[1] * vxVec->vtY) + yTr[2];
      if(vxCount-- > 0)
      {
        ++vxVec;
      }
    }
    else
    {
      /* Vertex is NOT in this element so walk towards the element
	 which does contain it. */
      if(elm->flags & neighbourMask)
      {
	trValid = 0;
	elmIdx = elm->neighbours[neighbourId];
	elm = mesh->elements + elmIdx;
	nod[0] = mesh->nodes + elm->nodes[0];
	nod[1] = mesh->nodes + elm->nodes[1];
	nod[2] = mesh->nodes + elm->nodes[2];
	elmVx[0] = nod[0]->position;
	elmVx[1] = nod[1]->position;
	elmVx[2] = nod[2]->position;
	if((elmArea2 = WlzGeomTriangleSnArea2(elmVx[0], elmVx[1],
					    elmVx[2])) < WLZ_MESH_TOLERANCE_SQ)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
      else
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Solve's a system of linear equations for the coefficients
*		of a 2D affine transform from the source triangle to the
*		destination triangle. Because we know that the area of the
*		triangle is NOT zero and that we have a small system, Cramer's
*		rule is used.
* \param	xTr			Transform coordinates for x.
* \param	yTr			Transform coordinates for y.
* \param	dd			Twice the area of the source triangle.
* \param	sVx			Source triangle vertices.
* \param	dVx			Destination triangle vertices.
*/
static void	WlzMeshAfTrSolve(double *xTr, double *yTr, double dd,
				 WlzDVertex2 *sVx, WlzDVertex2 *dVx)
{
  double	tD0,
  		tD1,
		tD2,
		dx0, dy0, sx0, sy0,
		dx1, dy1, sx1, sy1,
		dx2, dy2, sx2, sy2;

  dd = 1.0 / dd;
  dx0 = dVx->vtX; dy0 = dVx->vtY; ++dVx;
  dx1 = dVx->vtX; dy1 = dVx->vtY; ++dVx;
  dx2 = dVx->vtX; dy2 = dVx->vtY;
  sx0 = sVx->vtX; sy0 = sVx->vtY; ++sVx;
  sx1 = sVx->vtX; sy1 = sVx->vtY; ++sVx;
  sx2 = sVx->vtX; sy2 = sVx->vtY;
  tD0 = sy1 - sy2; tD1 = sy2 - sy0; tD2 = sy0 - sy1;
  *(xTr + 0) = ((dx0 * tD0) + (dx1 * tD1) + (dx2 * tD2)) * dd;
  *(yTr + 0) = ((dy0 * tD0) + (dy1 * tD1) + (dy2 * tD2)) * dd;
  tD0 = sx2 - sx1; tD1 = sx0 - sx2; tD2 = sx1 - sx0;
  *(xTr + 1) = ((dx0 * tD0) + (dx1 * tD1) + (dx2 * tD2)) * dd;
  *(yTr + 1) = ((dy0 * tD0) + (dy1 * tD1) + (dy2 * tD2)) * dd;
  tD0 = (sx1 * sy2) - (sx2 * sy1);
  tD1 = (sx2 * sy0) - (sx0 * sy2);
  tD2 = (sx0 * sy1) - (sx1 * sy0);
  *(xTr + 2) = ((dx0 * tD0) + (dx1 * tD1) + (dx2 * tD2)) * dd;
  *(yTr + 2) = ((dy0 * tD0) + (dy1 * tD1) + (dy2 * tD2)) * dd;
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Updates the destination mesh element data for the given
* 		element index.
* \param	mSnWSp			Mesh scan workspace.
* \param	eIdx			Element index.
*/
static WlzErrorNum WlzMeshScanDElmUpdate(WlzMeshScanWSp *mSnWSp, int eIdx)
{
  int		nodCnt;
  double	areaSn2;
  int		*nodIdxP;
  WlzDVertex2	*dstVxP,
  		*srcVxP;
  WlzMeshElem	*meshElems;
  WlzMeshNode	*nodP,
  		*meshNodes;
  WlzMeshScanDElm *dElm;
  WlzDVertex2	dstVx[3],
		srcVx[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dElm = mSnWSp->dElm + eIdx;
  meshElems = mSnWSp->mesh->elements;
  meshNodes = mSnWSp->mesh->nodes;
  nodIdxP = (meshElems + eIdx)->nodes;
  nodCnt = 3;
  dstVxP = dstVx;
  srcVxP = srcVx;
  while(nodCnt-- > 0)
  {
    nodP = meshNodes + *nodIdxP++;
    *srcVxP++ = *dstVxP = nodP->position;
    dstVxP->vtX += nodP->displacement.vtX;
    dstVxP++->vtY += nodP->displacement.vtY;
  }
  if((areaSn2 = WlzGeomTriangleSnArea2(dstVx[0], dstVx[1],
  				       dstVx[2])) < WLZ_MESH_TOLERANCE_SQ)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    WlzMeshAfTrSolve(dElm->xTr, dElm->yTr, areaSn2, dstVx, srcVx);
    dElm->valid = 1;
  }
  return(errNum);
}

/*!
* \return	New mesh scan workspace.
* \ingroup	WlzTransform
* \brief	Allocate and initialise a mesh scan workspace.
* \param	mesh			Mesh transform.
* \param	dstErr			Destination error pointer.
*/
static WlzMeshScanWSp *WlzMeshScanWSpInit(WlzMeshTransform *mesh,
				    	  WlzErrorNum *dstErr)
{
  int		iIdx,
  		eIdx,
		ndIdx;
  double	ndLn,
  		eLnMin,
		eLnMax;
  WlzMeshElem	*elm;
  WlzMeshScanWSp *meshSnWSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((meshSnWSp = (WlzMeshScanWSp *)
  		  AlcCalloc(1, sizeof(WlzMeshScanWSp))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    meshSnWSp->mesh = mesh;
    /* Compute the total number of intervals in the displaced mesh. */
    eIdx = 0;
    while((errNum == WLZ_ERR_NONE) && (eIdx < mesh->nElem))
    {
      elm = mesh->elements + eIdx;
      ndIdx = elm->nodes[0];
      eLnMin = eLnMax = ndLn = (mesh->nodes + ndIdx)->position.vtY +
			       (mesh->nodes + ndIdx)->displacement.vtY;
      ndIdx = elm->nodes[1];
      ndLn = (mesh->nodes + ndIdx)->position.vtY +
      	     (mesh->nodes + ndIdx)->displacement.vtY;
      if(ndLn < eLnMin)
      {
	eLnMin = ndLn;
      }
      else if(ndLn > eLnMax)
      {
	eLnMax = ndLn;
      }
      ndIdx = elm->nodes[2];
      ndLn = (mesh->nodes + ndIdx)->position.vtY +
             (mesh->nodes + ndIdx)->displacement.vtY;
      if(ndLn < eLnMin)
      {
	eLnMin = ndLn;
      }
      else if(ndLn > eLnMax)
      {
	eLnMax = ndLn;
      }
      meshSnWSp->nItvs += WLZ_NINT(eLnMax) - WLZ_NINT(eLnMin) + 1;
      ++eIdx;
    }
    if((meshSnWSp->itvs = (WlzMeshScanItv *)AlcMalloc(sizeof(WlzMeshScanItv) *
						    meshSnWSp->nItvs)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else if((meshSnWSp->dElm = (WlzMeshScanDElm *)
			       AlcCalloc(mesh->nElem,
					 sizeof(WlzMeshScanDElm))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill in the mesh scan intervals */
    eIdx = 0;
    iIdx = 0;
    while(eIdx < mesh->nElem)
    {
      iIdx += WlzMeshScanTriElm(meshSnWSp, eIdx, iIdx);
      ++eIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Sort the mesh scan intervals by line and then left column */
    qsort(meshSnWSp->itvs, meshSnWSp->nItvs, sizeof(WlzMeshScanItv),
          WlzMeshItvCmp);
  }
  else
  {
    WlzMeshScanWSpFree(meshSnWSp);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(meshSnWSp);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Free's a mesh scan workspace.
* \param	mSnWSp			Mesh scan workspace.
*/
static void	WlzMeshScanWSpFree(WlzMeshScanWSp *mSnWSp)
{
  if(mSnWSp)
  {
    if(mSnWSp->itvs)
    {
      AlcFree(mSnWSp->itvs);
    }
    if(mSnWSp->dElm)
    {
      AlcFree(mSnWSp->dElm);
    }
    AlcFree(mSnWSp);
  }
}

/*!
* \return	Number of intervals added from the given mesh element.
* \ingroup	WlzTransform
* \brief	Scans a single triangular mesh element into mesh intervals.
* \param	mSnWSp			Mesh scan workspace.
* \param	eIdx			Element index.
* \param	iIdx			Mesh element interval index.
*/
static int	WlzMeshScanTriElm(WlzMeshScanWSp *mSnWSp, int eIdx, int iIdx)
{
  int		count,
		kolI,
  		lineI,
		ndIdx0,
  		ndIdx1,
		ndIdx2,
		iCnt = 0;
  double	tD0,
  		tD1,
		kolD,
		inc;
  WlzIVertex2	dNd[3],
  		sNd[3];
  WlzMeshNode	*nod;
  WlzMeshElem	*elm;
  WlzDVertex2	dVx0,
  		dVx1;
  WlzMeshScanItv *itv;

  elm = mSnWSp->mesh->elements + eIdx;
  /* Compute the integer displaced nodes of the element. */
  for(ndIdx0 = 0; ndIdx0 < 3; ++ndIdx0)
  {
    ndIdx1 = elm->nodes[ndIdx0];
    nod = mSnWSp->mesh->nodes + ndIdx1;
    dVx0 = nod->position;
    dVx1 = nod->displacement;
    tD0 = dVx0.vtX + dVx1.vtX;
    tD1 = dVx0.vtY + dVx1.vtY;
    dNd[ndIdx0].vtX = WLZ_NINT(tD0);
    dNd[ndIdx0].vtY = WLZ_NINT(tD1);
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
  dNd[0].vtX = sNd[0].vtX - sNd[1].vtX;
  dNd[0].vtY = sNd[0].vtY - sNd[1].vtY;
  dNd[1].vtX = sNd[1].vtX - sNd[2].vtX;
  dNd[1].vtY = sNd[1].vtY - sNd[2].vtY;
  dNd[2].vtX = sNd[2].vtX - sNd[0].vtX;
  dNd[2].vtY = sNd[2].vtY - sNd[0].vtY;
  /* If the element's nodes are not coincident scan convert it. */
  if(dNd[2].vtY && (dNd[0].vtX || dNd[1].vtX))
  {
    /* Nodes: min -> max */
    itv = mSnWSp->itvs + iIdx;
    kolD = sNd[0].vtX;
    lineI = sNd[0].vtY;
    inc = (double )(dNd[2].vtX) / (double )(dNd[2].vtY);
    count = iCnt = sNd[2].vtY - sNd[0].vtY + 1;
    while(count-- > 0)
    {
      itv->elmIdx = eIdx;
      itv->line = lineI++;
      itv->lftI = itv->rgtI = WLZ_NINT(kolD);
      kolD += inc;
      ++itv;
    }
    if(dNd[0].vtY)
    {
      /* Nodes: mid -> min */
      itv = mSnWSp->itvs + iIdx;
      kolD = sNd[0].vtX;
      inc = (double )(dNd[0].vtX) / (double )(dNd[0].vtY);
      count = sNd[1].vtY - sNd[0].vtY + 1;
      while(count-- > 0)
      {
	kolI = WLZ_NINT(kolD);
        if(kolI > itv->lftI)
	{
	  itv->rgtI = kolI;
	}
	else
	{
	  itv->lftI = kolI;
	}
        kolD += inc;
	++itv;
      }
    }
    if(dNd[1].vtY)
    {
      /* Nodes: max -> mid */
      itv = mSnWSp->itvs + iIdx + iCnt - 1;
      kolD = sNd[2].vtX;
      inc = (double )(dNd[1].vtX) / (double )(dNd[1].vtY);
      count = sNd[2].vtY - sNd[1].vtY + 1;
      while(count-- > 0)
      {
	kolI = WLZ_NINT(kolD);
        if(kolI > itv->lftI)
	{
	  itv->rgtI = kolI;
	}
	else
	{
	  itv->lftI = kolI;
	}
        kolD -= inc;
	--itv;
      }
    }
  }
  return(iCnt);
}

/*!
* \return	Sorting value for qsort.
* \ingroup	WlzTransform
* \brief	Callback function for qsort(3) to sort mesh element
*		intervals by line and then left left column.
* \param	cmp0			Used to pass first mesh interval.
* \param	cmp1			Used to pass second mesh interval.
*/
static int	WlzMeshItvCmp(const void *cmp0, const void *cmp1)
{
  int		rtn;
  WlzMeshScanItv *itv0,
  		 *itv1;

  itv0 = (WlzMeshScanItv *)cmp0;
  itv1 = (WlzMeshScanItv *)cmp1;
  if((rtn = (itv0->line - itv1->line)) == 0)
  {
    rtn = itv0->lftI - itv1->lftI;
  }
  return(rtn);
}
