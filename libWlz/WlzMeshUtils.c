#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMeshUtils.c
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
* \brief	Utility functions for manipulating Woolz mesh transforms.
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

typedef struct _WlzMeshIntVec
{
  int		*vector;
  int		size;
  int		count;
  int		minSize;
  int		mulSize;
} WlzMeshIntVec;

typedef struct _WlzMeshEar
{
  int		flags;
  double	power;
  double	snArea2;
  int		nodes[3];
  int		neighbours[3];
  int		selfNeighbour[3];
  struct _WlzMeshEar *next;
  struct _WlzMeshEar *prev;
} WlzMeshEar;

typedef struct _WlzMeshEarList
{
  int		inEars;
  int		maxEars;
  WlzMeshEar	*topIn;
  WlzMeshEar	*pool;
} WlzMeshEarList;

static int	WlzMeshElemConflict(WlzMeshTransform *, int, WlzDVertex2);
static void	WlzMeshElemFindVxForce(WlzMeshTransform *, WlzDVertex2,
				int *, int *, int *),
		WlzMeshElemReplace1With1(WlzMeshTransform *, int, int),
		WlzMeshElemReplace1With2(WlzMeshTransform *, int, int,
				WlzDVertex2, unsigned int),
		WlzMeshElemReplace1With3(WlzMeshTransform *, int,
				WlzDVertex2, unsigned int),
		WlzMeshElemUnlink(WlzMeshTransform *, int),
		WlzMeshEarPowerSet(WlzMeshTransform *, WlzMeshEar *, int),
		WlzMeshEarMatchElm(WlzMeshTransform *, WlzMeshEar *,
				WlzMeshIntVec *, int *),
		WlzMeshNodeDelInit(WlzMeshIntVec *, WlzMeshIntVec *,
			   	WlzMeshEarList *),
		WlzMeshNodeDelFree(WlzMeshIntVec *, WlzMeshIntVec *, 
				WlzMeshEarList *);
static WlzErrorNum WlzMeshAddToIntVec(WlzMeshIntVec *, int),
		WlzMeshNodeDel(WlzMeshTransform *, 
				WlzMeshIntVec *, WlzMeshIntVec *,
				WlzMeshEarList *, int, WlzDVertex2),
		WlzMeshIDomAdd(WlzMeshTransform *, WlzObject *,
				double, WlzDVertex2),
		WlzMeshPolyDomAdd(WlzMeshTransform *, WlzObject *,
				double, WlzDVertex2),
		WlzMeshElemFindVxWalk(WlzMeshTransform *, WlzDVertex2,
				int *, int *, int *),
		WlzMeshQueConflictElem(WlzMeshIntVec *,
				WlzMeshTransform *, int, int, WlzDVertex2),
		WlzMeshElemReplaceN(WlzMeshTransform *, int *, int,
				WlzDVertex2, unsigned int),
		WlzMeshElemReplace1(WlzMeshTransform *, int, WlzDVertex2,
				unsigned int),
		WlzMeshElemReplaceNWithN(WlzMeshTransform *, int *,
				int, WlzDVertex2, unsigned int),
		WlzMeshNodeDelVecBuild(WlzMeshIntVec *, WlzMeshIntVec *,
				WlzMeshTransform *, int, int),
		WlzMeshEarsCreate(WlzMeshEarList *,
				WlzMeshIntVec *, WlzMeshIntVec *,
				WlzMeshTransform *, int),
		WlzMeshEarListRealloc(WlzMeshEarList *, int);
static WlzMeshEar *WlzMeshEarGetMinPower(WlzMeshEarList *);

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds mesh nodes within the domain of the given 2D domain
* 		object.
* \param	mesh			Given mesh transform.
* \param	obj			Given object with domain.
* \param	minDist			Minimum distance between mesh nodes.
* \param	scaleVx			Object scale factor.
*/
WlzErrorNum	WlzMeshDomainAdd(WlzMeshTransform *mesh, WlzObject *obj,
				 double minDist, WlzDVertex2 scaleVx)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((obj->domain.core == NULL) || (mesh == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->nElem < 1)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	errNum = WlzMeshIDomAdd(mesh, obj, minDist * minDist, scaleVx);
	break;
      case WLZ_2D_POLYGON:
	errNum = WlzMeshPolyDomAdd(mesh, obj, minDist, scaleVx);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds mesh nodes at the vertices of the given (2D) vertex
*		vector.
* \param	mesh			Given mesh transform.
* \param	vxVec			Given vector of vertices.
* \param	nVx			Number of vertices.
* \param	minDistSq		Square of the minimum distance
*					between mesh nodes.
* \param	nodeFlags		Node flags to set (eg source).
*/
WlzErrorNum	WlzMeshVxVecAdd(WlzMeshTransform *mesh, WlzDVertex2 *vxVec,
				int nVx, double minDistSq,
				unsigned int nodeFlags)
{
  int		vxCnt,
		eId0,
		nId,
		qId,
		insFlg,
		tryVxFlg,
		extFlg,
		firstVxFlg = 1,
		startElm = 0;
  double	tD0,
		tD1;
  WlzMeshNode	*nod;
  WlzMeshElem	*elm;
  WlzDVertex2	*vxVecP0;
  WlzDVertex2	ndVx0,
		lastVx;
  WlzMeshIntVec	eCnfQVec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vxCnt = nVx;
  vxVecP0 = vxVec;
  eCnfQVec.vector = NULL;
  eCnfQVec.size = 0;
  eCnfQVec.minSize = 8;
  eCnfQVec.mulSize = 2;
  lastVx.vtX = lastVx.vtY = 0; 			 /* Just to keep lint happy. */
  while((errNum == WLZ_ERR_NONE) && (vxCnt-- > 0))
  {
    if(firstVxFlg)
    {
      tryVxFlg = 1;
    }
    else
    {
      tD0 = lastVx.vtX - vxVecP0->vtX;
      tD1 = lastVx.vtY - vxVecP0->vtY;
      tryVxFlg = (tD0 * tD0) + (tD1 * tD1) >= minDistSq;
    }
    if(tryVxFlg)
    {
      /* Find new vertex in the mesh. */
      eId0 = WlzMeshElemFindVx(mesh, *vxVecP0, startElm, &extFlg, &errNum);
      if((errNum == WLZ_ERR_NONE) && (eId0 >= 0))
      {
	elm = mesh->elements + eId0;
	if(extFlg)
	{
	  /* Node already exists find which node of element and modify
	   * it's flags. */
	  if((nId = WlzMeshElemNodeIdxFromVx(mesh, elm, *vxVecP0)) >= 0)
	  {
	    nod = mesh->nodes + *(elm->nodes + nId);
	    nod->flags = nodeFlags;
	  }
	  else
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
	else
	{
	  tD0 = WlzGeomTriangleSnArea2((mesh->nodes + elm->nodes[0])->position,
				       (mesh->nodes + elm->nodes[1])->position,
				       (mesh->nodes + elm->nodes[2])->position);
	  if(tD0 >= minDistSq);
	  {
	    /* Find elements in conflict with the new node vertex. */
	    eCnfQVec.count = 0;
	    errNum = WlzMeshQueConflictElem(&eCnfQVec, mesh, eId0, eId0,
					    *vxVecP0);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      /* Check new vertex is outside of minimum distance from node. */
	      qId = 0;
	      insFlg = 1;
	      while(insFlg && (qId < eCnfQVec.count))
	      {
		nId = 0;
		elm = mesh->elements + *(eCnfQVec.vector + qId);
		while(insFlg && (nId < 3))
		{
		  nod = mesh->nodes + *(elm->nodes + nId);
		  if((nod->flags & WLZ_MESH_NODE_FLAGS_BBOX) == 0)
		  {
		    ndVx0 = nod->position;
		    tD0 = ndVx0.vtX - vxVecP0->vtX;
		    tD1 = ndVx0.vtY - vxVecP0->vtY;
		    insFlg = (tD0 * tD0) + (tD1 * tD1) >= minDistSq;
		  }
		  ++nId;
		}
		++qId;
	      }
	      if(insFlg)
	      {
		/* Replace conflicting elements and add new elements. */
		errNum = WlzMeshElemReplaceN(mesh, eCnfQVec.vector,
					     eCnfQVec.count, *vxVecP0,
					     nodeFlags);
		firstVxFlg = 0;
		lastVx = *vxVecP0;
	      }
	      else
	      {
		/* Clear zombie flags if elements not replaced. */
		qId = 0;
		while(qId < eCnfQVec.count)
		{
		  elm = mesh->elements + *(eCnfQVec.vector + qId);
		  elm->flags = elm->flags & ~(WLZ_MESH_ELEM_FLAGS_ZOMBIE);
		  ++qId;
		}
	      }
	      startElm = eId0;
	    }
	  }
        }
      }
    }
    ++vxVecP0;
  }
  if(eCnfQVec.vector)
  {
    AlcFree(eCnfQVec.vector);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds mesh nodes within the interval domain of the given
*		2D domain object.
* \param	mesh			Given mesh transform.
* \param	obj			Given object with domain.
* \param	minDistSq		Square of the minimum distance
*					between mesh nodes.
* \param	scaleVx			Object scale factor.
*/
WlzErrorNum	WlzMeshIDomAdd(WlzMeshTransform *mesh, WlzObject *obj,
			       double minDistSq, WlzDVertex2 scaleVx)
{
  int		eId0,
		nId,
		qId,
		insFlg,
		tryVxFlg,
		extFlg,
		firstVxFlg = 1,
		startElm = 0;
  double	tD0,
		tD1;
  WlzMeshElem	*elm;
  WlzMeshNode	*nod;
  WlzDVertex2	ndVx,
		ivVx,
		newVx,
		lastVx;
  WlzMeshIntVec	eCnfQVec;
  WlzIntervalWSpace iWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  eCnfQVec.vector = NULL;
  eCnfQVec.size = 0;
  eCnfQVec.minSize = 8;
  eCnfQVec.mulSize = 2;
  errNum = WlzInitRasterScan(obj, &iWsp, WLZ_RASTERDIR_ILIC);
  while((errNum == WLZ_ERR_NONE) &&
	((errNum = WlzNextInterval(&iWsp)) == WLZ_ERR_NONE))
  {
    ivVx.vtX = iWsp.lftpos;
    ivVx.vtY = iWsp.linpos;
    newVx.vtX = scaleVx.vtX * ivVx.vtX;
    newVx.vtY = scaleVx.vtY * ivVx.vtY;
    while((errNum == WLZ_ERR_NONE) && (ivVx.vtX <= iWsp.rgtpos))
    {
      if(firstVxFlg)
      {
	tryVxFlg = 1;
      }
      else
      {
	tD0 = lastVx.vtX - newVx.vtX;
	tD1 = lastVx.vtY - newVx.vtY;
	tryVxFlg = (tD0 * tD0) + (tD1 * tD1) >= minDistSq;
      }
      if(tryVxFlg)
      {
	/* Find new vertex in the mesh. */
	eId0 = WlzMeshElemFindVx(mesh, newVx, startElm, &extFlg, &errNum);
	if((errNum == WLZ_ERR_NONE) && (eId0 >= 0))
	{
	  elm = mesh->elements + eId0;
	  if(extFlg)
	  {
	    /* Node already exists find which node of element and modify
	     * it's flags. */
	    if((nId = WlzMeshElemNodeIdxFromVx(mesh, elm, newVx)) > 0)
	    {
	      nod = mesh->nodes + *(elm->nodes + nId);
	      nod->flags = WLZ_MESH_NODE_FLAGS_IDOM;
	    }
	    else
	    {
	      errNum = WLZ_ERR_DOMAIN_DATA;
	    }
	  }
	  else
	  {
	    tD0 = WlzGeomTriangleSnArea2((mesh->nodes +
	    				  elm->nodes[0])->position,
					 (mesh->nodes +
					  elm->nodes[1])->position,
					 (mesh->nodes +
					  elm->nodes[2])->position);
	    if(tD0 >= minDistSq);
	    {
	      /* Find elements in conflict with the new node vertex. */
	      eCnfQVec.count = 0;
	      errNum = WlzMeshQueConflictElem(&eCnfQVec, mesh, eId0, eId0,
	      				      newVx);
	      if(errNum == WLZ_ERR_NONE)
	      {
		/* Check newVx is outside of minimum distance from node. */
		qId = 0;
		insFlg = 1;
		while(insFlg && (qId < eCnfQVec.count))
		{
		  nId = 0;
		  elm = mesh->elements + *(eCnfQVec.vector + qId);
		  while(insFlg && (nId < 3))
		  {
		    nod = mesh->nodes + *(elm->nodes + nId);
		    if((nod->flags & WLZ_MESH_NODE_FLAGS_BBOX) == 0)
		    {
		      ndVx = nod->position;
		      tD0 = ndVx.vtX - newVx.vtX;
		      tD1 = ndVx.vtY - newVx.vtY;
		      insFlg = (tD0 * tD0) + (tD1 * tD1) >= minDistSq;
		    }
		    ++nId;
		  }
		  ++qId;
		}
		if(insFlg)
		{
		  /* Replace conflicting elements and add new elements. */
		  errNum = WlzMeshElemReplaceN(mesh, eCnfQVec.vector,
					       eCnfQVec.count, newVx,
					       WLZ_MESH_NODE_FLAGS_IDOM);
		  firstVxFlg = 0;
		  lastVx = newVx;
		}
		else
		{
		  /* Clear zombie flags if elements not replaced. */
		  qId = 0;
		  while(qId < eCnfQVec.count)
		  {
		    elm = mesh->elements + *(eCnfQVec.vector + qId);
		    elm->flags = elm->flags & ~(WLZ_MESH_ELEM_FLAGS_ZOMBIE);
		    ++qId;
		  }
		}
		startElm = eId0;
	      }
	    }
	  }
	}
      }
      newVx.vtX += scaleVx.vtX;
      ++(ivVx.vtX);
    }
  }
  if(errNum == WLZ_ERR_EOO)		/* Reset error from end of intervals */ 
  {
    errNum = WLZ_ERR_NONE;
  }
  if(eCnfQVec.vector)
  {
    AlcFree(eCnfQVec.vector);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds mesh nodes along the polygon domain of the given
*		2D domain object.
* \param	mesh			Given mesh transform.
* \param	obj			Given object with domain.
* \param	minDist			Minimum distance between mesh nodes.
* \param	scaleVx			Object scale factor.
*/
WlzErrorNum	WlzMeshPolyDomAdd(WlzMeshTransform *mesh, WlzObject *obj,
				  double minDist, WlzDVertex2 scaleVx)
{
  int		vxCnt0,
		vxCnt2,
		vxId0,
		dVxCnt,
		sVxCnt;
  double	tD0,
		tD1,
		tD2,
		tD3,
		tD4;
  WlzDVertex2	*tVxP0,
		*tVxP1,
		*tVxP2,
		*dVxVec = NULL,
		*sVxVec = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((sVxCnt = obj->domain.poly->nvertices) > 0)
  {
    if((sVxVec = (WlzDVertex2 *)
		 AlcMalloc(sizeof(WlzDVertex2) * sVxCnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      /* Convert to double vertices. */
      switch(obj->domain.poly->type)
      {
	case WLZ_POLYGON_INT:
	  WlzValueCopyIVertexToDVertex(sVxVec,
				       obj->domain.poly->vtx,
				       sVxCnt);
	  break;
	case WLZ_POLYGON_FLOAT:
	  WlzValueCopyFVertexToDVertex(sVxVec, 
				       (WlzFVertex2 *)(obj->domain.poly->vtx),
				       sVxCnt);
	  break;
	case WLZ_POLYGON_DOUBLE:
	  WlzValueCopyDVertexToDVertex(sVxVec, 
				       (WlzDVertex2 *)(obj->domain.poly->vtx),
				       sVxCnt);
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      /* Scale vertices. */
      if((fabs(scaleVx.vtX - 1.0) > DBL_EPSILON) ||
	 (fabs(scaleVx.vtY - 1.0) > DBL_EPSILON))
      {
	vxCnt0 = sVxCnt;
	tVxP0 = sVxVec;
	while(vxCnt0-- > 0)
	{
	  tVxP0->vtX *= scaleVx.vtX;
	  tVxP0->vtY *= scaleVx.vtY;
	  ++tVxP0;
	}
      }
      if(sVxCnt > 1)
      {
	/* Make a new vector of vertices in-between the given vertices. */
	/* Compute in-between vector size. */
	dVxCnt = 0;
	vxCnt0 = sVxCnt;
	tVxP0 = sVxVec + sVxCnt - 1;
	tVxP1 = sVxVec;
	tD2 = 4 * minDist * minDist;
	while(vxCnt0-- > 0)
	{
	  tD0 = tVxP1->vtX - tVxP0->vtX;
	  tD1 = tVxP1->vtY - tVxP0->vtY;
	  dVxCnt += floor(sqrt(((tD0 * tD0) + (tD1 * tD1)) / tD2));
	  tVxP0 = tVxP1++;
	}

      }
      if(dVxCnt > 0)
      {
	if((dVxVec = (WlzDVertex2 *)
		     AlcMalloc(sizeof(WlzDVertex2) * dVxCnt)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  vxCnt0 = sVxCnt;
	  tVxP0 = sVxVec;
	  tVxP1 = sVxVec + sVxCnt - 1;
	  tVxP2 = dVxVec;
	  tD2 = 4 * minDist * minDist;
	  while(vxCnt0-- > 0)
	  {
	    tD0 = tVxP1->vtX - tVxP0->vtX;
	    tD1 = tVxP1->vtY - tVxP0->vtY;
	    vxCnt2 = floor(sqrt(((tD0 * tD0) + (tD1 * tD1)) / tD2));
	    tD3 = vxCnt2 + 1;
	    for(vxId0 = 1; vxId0 <= vxCnt2; ++vxId0)
	    {
	      tD4 = vxId0 / tD3;
	      tVxP2->vtX = tVxP0->vtX + (tD0 * tD4);
	      tVxP2->vtY = tVxP0->vtY + (tD1 * tD4);
	      ++tVxP2;
	    }
	    tVxP1 = tVxP0++;
	  }
	}
      }
    }
    if((errNum == WLZ_ERR_NONE) && (sVxCnt > 0))
    {
      /* Add given vertices to the mesh. */
      errNum = WlzMeshVxVecAdd(mesh, sVxVec, sVxCnt, minDist * minDist,
      			       WLZ_MESH_NODE_FLAGS_POLY);
    }
    if((errNum == WLZ_ERR_NONE) && (dVxCnt > 0))
    {
      /* Add in-between vertices to the mesh. */
      errNum = WlzMeshVxVecAdd(mesh, dVxVec, dVxCnt, minDist * minDist,
      			       WLZ_MESH_NODE_FLAGS_POLY);
    }
    if(sVxVec)
    {
      AlcFree(sVxVec);
    }
    if(dVxVec)
    {
      AlcFree(dVxVec);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Checks that the given mesh transform is valid.
* \param	mesh			Given mesh transform.
* \param	dispFlg			Verify displacements if non-zero.
* \param	elm			Element to verify.
* \param	dstErrMsk		Destination pointer to be set with
*					a mesh error after verifying this
*					element, may be NULL.
*/
WlzErrorNum	WlzMeshElemVerify(WlzMeshTransform *mesh, int dispFlg,
			          WlzMeshElem *elm, WlzMeshError *dstErrMsk)
{
  int		nEId,
		nId0,
		nId1,
		ndId0,
		ndId1,
		nnId;
  WlzDVertex2	sNd[3],
		dNd[3];
  WlzMeshNode	*nod;
  WlzMeshElem	*nElm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzMeshError	errMsk = WLZ_MESH_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(elm == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(elm->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE)
  {
    /* This element is an unreclaimed zombie. */
    errNum = WLZ_ERR_DOMAIN_DATA;
    errMsk = WLZ_MESH_ERR_ELEM_ZOMBIE;
  }
  else
  {
    /* Check nodes are valid. */
    nId0 = 0;
    while((errNum == WLZ_ERR_NONE) && (nId0 < 3))
    {
      ndId0 = elm->nodes[nId0];
      if((ndId0 < 0) || (ndId0 > mesh->nNodes))
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
	errMsk = WLZ_MESH_ERR_ELEM_NODE;
      }
      else
      {
	nod = mesh->nodes + ndId0;
	sNd[nId0].vtX = nod->position.vtX;
	sNd[nId0].vtY = nod->position.vtY;
	if(dispFlg)
	{
	  dNd[nId0].vtX = sNd[nId0].vtX + nod->displacement.vtX;
	  dNd[nId0].vtY = sNd[nId0].vtY + nod->displacement.vtY;
	}
      }
      ++nId0;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check for CCW node and displaced node order. */
    if(dispFlg && (WlzGeomTriangleSnArea2(dNd[0], dNd[1],
			      		  dNd[2]) < WLZ_MESH_TOLERANCE_SQ))
    {
      /* Element nodes are not CCW. */
      errNum = WLZ_ERR_DOMAIN_DATA;
      errMsk = WLZ_MESH_ERR_DELEM_CW;
    }
    else if (WlzGeomTriangleSnArea2(sNd[0], sNd[1],
				    sNd[2]) < WLZ_MESH_TOLERANCE_SQ)
    {
      /* Displaced element nodes are not CCW. */
      errNum = WLZ_ERR_DOMAIN_DATA;
      errMsk = WLZ_MESH_ERR_ELEM_CW;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check connectivity with each of the elements neighbours. */
    nId0 = 0;
    while((nId0 < 3) && (errNum == WLZ_ERR_NONE))
    {
      if(elm->flags & nbrFlgTbl[nId0])      /* Does nId0'th neighbour exist. */
      {
	nEId = elm->neighbours[nId0];
	if((nEId < 0) || (nEId > mesh->nElem))
	{
	  /* Invalid neighbouring element index. */
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  errMsk = WLZ_MESH_ERR_NELEM_INDEX;
	}
	else
	{
	  nElm = mesh->elements + elm->neighbours[nId0];
	  if(nElm->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE)
	  {
	    /* Neighbour is an unreclaimed zombie. */
	    errNum = WLZ_ERR_DOMAIN_DATA;
	    errMsk = WLZ_MESH_ERR_NELEM_ZOMBIE;
	  }
	  else
	  {
	    ndId0 = elm->nodes[(nId0 + 1) % 3];
	    ndId1 = elm->nodes[(nId0 + 2) % 3];
	    /* Find which of the neighbours neighbours shares these nodes. */
	    nnId = 0;
	    while((nnId < 3) && (ndId0 != nElm->nodes[nnId]))
	    {
	      ++nnId;
	    }
	    if(nnId > 2)
	    {
	      /* 1st node which should be shared with neighbour isn't. */
	      errNum = WLZ_ERR_DOMAIN_DATA;
	      errMsk = WLZ_MESH_ERR_NELEM_NODE;
	    }
	    else
	    {
	      if(ndId1 == nElm->nodes[(nnId + 1) % 3])
	      {
		nId1 = (nnId + 2) % 3;
	      }
	      else if(ndId1 == nElm->nodes[(nnId + 2) % 3])
	      {
		nId1 = (nnId + 1) % 3;
	      }
	      else
	      {
		/* 2nd node which should be shared with neighbour isn't. */
		errNum = WLZ_ERR_DOMAIN_DATA;
		errMsk = WLZ_MESH_ERR_NELEM_NODE;
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		if(((nElm->flags & nbrFlgTbl[nId1]) == 0) ||
		   (nElm->neighbours[nId1] != elm->idx))
		{
		  /* Neighbouring element which should share these nodes
		   * with this element shares them with another element.  */
		  errNum = WLZ_ERR_DOMAIN_DATA;
		  errMsk = WLZ_MESH_ERR_NELEM_NOTNBR;
		}
	      }
	    }
	  }
	}
      }
      ++nId0;
    }
  }
  if(dstErrMsk)
  {
    *dstErrMsk = errMsk;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Expands a mesh to make sure that there are enough mesh
*		elements and nodes available.
* \param	mesh			Given mesh transform.
* \param	nElem			Minimum number of mesh elements
* 					required.
* \param	nNodes			Minimum number of mesh nodes required.
*/
WlzErrorNum	WlzMeshExpand(WlzMeshTransform *mesh, int nElem, int nNodes)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minElem = 64,
		minNodes = 64;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_TRANSFORM_2D_MESH)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(mesh->maxElem < nElem)
    {
      mesh->maxElem = ((mesh->maxElem < minElem)? minElem: mesh->maxElem) * 2;
      if((mesh->elements = (WlzMeshElem *)
			   AlcRealloc(mesh->elements,
				      sizeof(WlzMeshElem) *
				      mesh->maxElem)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(mesh->maxNodes < nNodes)
      {
	mesh->maxNodes = ((mesh->maxNodes < minNodes)? minNodes:
						       mesh->maxNodes) * 2;
	if((mesh->nodes = (WlzMeshNode *)
			  AlcRealloc(mesh->nodes,
				     sizeof(WlzMeshNode) *
				     mesh->maxNodes)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Squeeze out any zombie nodes and/or elements so that they
*		are available for reuse.
* \param	mesh			Given mesh transform.
*/
WlzErrorNum	WlzMeshSqueeze(WlzMeshTransform *mesh)
{
  int		eId0,
  		eId1,
		nId0,
  		nId1,
		sqElem,
		sqNodes,
		tblSz;
  int		*tbl = NULL,
  		*nIdP;
  WlzMeshNode	*nod0,
  		*nod1;
  WlzMeshElem	*elm0,
  		*elm1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  /* Make a table for squeezing node and element id's. */
  tblSz = (mesh->nElem > mesh->nNodes)? mesh->nElem: mesh->nNodes;
  if((tbl = AlcMalloc(tblSz * sizeof(int))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Squeeze out zombie nodes and build node table. */
    nId0 = 0;
    nod0 = mesh->nodes;
    while((nId0 < mesh->nNodes) &&
	  ((nod0->flags & WLZ_MESH_NODE_FLAGS_ZOMBIE) == 0))
    {
      *(tbl + nId0) = nId0;
      ++nId0;
      ++nod0;
    }
    nId1 = nId0 + 1;
    nod1 = nod0 + 1;
    while(nId1 < mesh->nNodes)
    {
      *(tbl + nId1) = nId0;
      if((nod1->flags & WLZ_MESH_NODE_FLAGS_ZOMBIE) == 0)
      {
	*nod0++ = *nod1;
	++nId0;
      }
      ++nId1;
      ++nod1;
    }
    sqNodes = nId0;
    /* Set element nodes using table. */
    eId0 = 0;
    elm0 = mesh->elements;
    while((eId0 < mesh->nElem) && (errNum == WLZ_ERR_NONE))
    {
      if((elm0->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
      {
        nId0 = 0;
	nIdP = elm0->nodes;
	while(nId0 < 3)
	{
	  if((*nIdP < 0) || (*nIdP >= mesh->nNodes))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    *nIdP = tbl[*nIdP];
	  }
	  ++nIdP;
	  ++nId0;
	}
      }
      ++elm0;
      ++eId0;
    }
    /* Set new value for number of nodes in the mesh. */
    mesh->nNodes = sqNodes;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Squeeze out zombie elements and build element table. */
    elm0 = mesh->elements;
    eId0 = 0;
    while((eId0 < mesh->nElem) &&
	  ((elm0->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0))
    {
      *(tbl + eId0) = eId0;
      ++eId0;
      ++elm0;
    }
    eId1 = eId0 + 1;
    elm1 = elm0 + 1;
    while(eId1 < mesh->nElem)
    {
      *(tbl + eId1) = eId0;
      if((elm1->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
      {
	*elm0++ = *elm1;
	++eId0;
      }
      ++eId1;
      ++elm1;
    }
    sqElem = eId0;
    /* Set element indicies and element neighbours using table. */
    eId0 = 0;
    elm0 = mesh->elements;
    while((eId0 < sqElem) && (errNum == WLZ_ERR_NONE))
    {
      elm0->idx = eId0;
      nId0 = 0;
      while(nId0 < 3)
      {
	if(elm0->flags & nbrFlgTbl[nId0])
	{
          nIdP = &(elm0->neighbours[nId0]);
	  if((*nIdP < 0) || (*nIdP >= mesh->nElem))
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    *nIdP = tbl[*nIdP];
	  }
	}
	++nId0;
      }
      ++elm0;
      ++eId0;
    }
    /* Set new value for number of elements in the mesh. */
    mesh->nElem = sqElem;
  }
  if(tbl)
  {
    AlcFree(tbl);
  }
  return(errNum);
}

/*!
* \return	Index to the node in the given element or -1 if the vertex
*		is not a node of the given element or two nodes are coincident
*		at the given vertex.
* \ingroup	WlzTransform
* \brief	Finds which of the given elements nodes are coincident with
* 		the given vertex.
* \param	mesh			Mesh transform.
* \param	elm			Mesh element.
* \param	gVx			Given vertex.
*/
int		WlzMeshElemNodeIdxFromVx(WlzMeshTransform *mesh,
					 WlzMeshElem *elm, WlzDVertex2 gVx)
{
  int		idx,
		wchCnt = 0;
  WlzDVertex2	vx;

  if(mesh && elm)
  {
    vx = (mesh->nodes + elm->nodes[0])->position;
    if((fabs(vx.vtX - gVx.vtX) < DBL_EPSILON) &&
       (fabs(vx.vtY - gVx.vtY) < DBL_EPSILON))
    {
      ++wchCnt;
      idx = 0;
    }
    vx = (mesh->nodes + elm->nodes[1])->position;
    if((fabs(vx.vtX - gVx.vtX) < DBL_EPSILON) &&
       (fabs(vx.vtY - gVx.vtY) < DBL_EPSILON))
    {
      ++wchCnt;
      idx = 1;
    }
    vx = (mesh->nodes + elm->nodes[2])->position;
    if((fabs(vx.vtX - gVx.vtX) < DBL_EPSILON) &&
       (fabs(vx.vtY - gVx.vtY) < DBL_EPSILON))
    {
      ++wchCnt;
      idx = 2;
    }
  }
  if(wchCnt != 1)
  {
    idx = -1;
  }
  return(idx);
}

/*!
* \return	Index to the node in the given element or  -1 if the vertex
*		is not a node of the given element or two nodes are coincident
*		at the given vertex.
* \ingroup	WlzTransform
* \brief	Finds which of the given ear/element nodes has the given
*		mesh node index.
* \param	nodes			Mesh ear/element nodes.
* \param	mNodId			Given mesh node index.
*/
int		WlzMeshElemNodeIdxFromNodeIdx(int *nodes, int mNodId)
{
  int		idx,
		wchCnt = 0;

  if(nodes)
  {
    if(*(nodes + 0) == mNodId)
    {
      ++wchCnt;
      idx = 0;
    }
    if(*(nodes + 1) == mNodId)
    {
      ++wchCnt;
      idx = 1;
    }
    if(*(nodes + 2) == mNodId)
    {
      ++wchCnt;
      idx = 2;
    }
  }
  if(wchCnt != 1)
  {
    idx = -1;
  }
  return(idx);
}

/*!
* \return	Index to the neighbour of the given element which would (if
* 		exists) share given nodes, or -1 if both nodes not shared with
* 		the element.
* \ingroup	WlzTransform
* \brief	Finds which neighbour of the given element whould share
*		the nodes has the given pair of nodes.
* \param	elm			Mesh element.
* \param	nodId0			First node.
* \param	nodId1			Second node.
*/
int		WlzMeshElemNbrIdxFromNodes(WlzMeshElem *elm,
					   int nodId0, int nodId1)
{
  int		nbrId;

  if(((elm->nodes[0] == nodId0) && (elm->nodes[1] == nodId1)) ||
     ((elm->nodes[0] == nodId1) && (elm->nodes[1] == nodId0)))
  {
    nbrId = 2;
  }
  else if(((elm->nodes[1] == nodId0) && (elm->nodes[2] == nodId1)) ||
          ((elm->nodes[1] == nodId1) && (elm->nodes[2] == nodId0)))
  {
    nbrId = 0;
  }
  else if(((elm->nodes[2] == nodId0) && (elm->nodes[0] == nodId1)) ||
          ((elm->nodes[2] == nodId1) && (elm->nodes[0] == nodId0)))
  {
    nbrId = 1;
  }
  else
  {
    nbrId = -1;
  }
  return(nbrId);
}

/*!
* \return	Element index, negative if not found.
* \ingroup	WlzTransform
* \brief	Searches the mesh for the element which contains the given
* 		vertex. It is NOT an error if the vertex is not within the
* 		mesh.
* \param	mesh			Given mesh transform.
* \param	gvnVx			Given vertex.
* \param	startElm		If >= 0, the index of the element from
* 					which to start the search.
* \param	existsFlg		Destination ptr for vertex already
* 					exists flag.
* \param	dstErr			Destination pointer for error number,
* 					may be NULL.
*/
int		WlzMeshElemFindVx(WlzMeshTransform *mesh, WlzDVertex2 gvnVx,
				  int startElm, int *existsFlg,
				  WlzErrorNum *dstErr)
{
  int		elmId,
		extFlg = 0,
		fndFlg = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_TRANSFORM_2D_MESH)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(startElm >= mesh->nElem)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(mesh->nElem > 0)
  {
    elmId = ((startElm < 0) || (startElm >= mesh->nElem))? 0: startElm;
    errNum = WlzMeshElemFindVxWalk(mesh, gvnVx, &elmId, &fndFlg, &extFlg);
    if((fndFlg == 0) && (errNum == WLZ_ERR_NONE))
    {
      WlzMeshElemFindVxForce(mesh, gvnVx, &elmId, &fndFlg, &extFlg);
    }
  }
  if(fndFlg == 0)
  {
    elmId = -1;
  }
  if(existsFlg)
  {
    *existsFlg = extFlg;
  }
  if(*dstErr)
  {
    *dstErr = errNum;
  }
  return(elmId);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Splits the given mesh element by placing a new node at
*		it's circumcentre.
* \param	mesh			Given mesh transform.
* \param	sElmId			Index of element to split.
*/
WlzErrorNum	WlzMeshElemSplit(WlzMeshTransform *mesh, int sElmId)
{
  int		bndFlg = 0,
		fndFlg = 0,
		countOut;
  double	sElmArea,
		sElmArea0,
		sElmArea1;
  WlzMeshElem	*sElm;
  WlzMeshNode	*nodes;
  WlzDVertex2	nVx,
		sVx0,
		sVx1,
		sVx2;
  WlzMeshIntVec	eCnfQVec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  eCnfQVec.vector = NULL;
  eCnfQVec.size = 0;
  eCnfQVec.minSize = 8;
  eCnfQVec.mulSize = 2;
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_TRANSFORM_2D_MESH)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(sElmId >= mesh->nElem)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(mesh->nElem > 0)
  {
    nodes = mesh->nodes;
    sElm = mesh->elements + sElmId;
    sVx0 = (nodes + sElm->nodes[0])->position;
    sVx1 = (nodes + sElm->nodes[1])->position;
    sVx2 = (nodes + sElm->nodes[2])->position;
    if(WlzGeomTriangleCircumcentre(&nVx, sVx0, sVx1, sVx2) == 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find which element the circumcentre lies in and if it lies across
     * the boundary of the mesh. If the circumcentre is across the mesh
     * boundary then use the midpoint of the boundary segment as the new
     * vertex. */
    countOut = mesh->nNodes;
    do
    {
      if((sElmArea = WlzGeomTriangleSnArea2(sVx0, sVx1,
					 sVx2)) < WLZ_MESH_TOLERANCE_SQ)
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
      else
      {
	if((sElmArea0 = WlzGeomTriangleSnArea2(sVx1, sVx2, nVx)) < 0.0)
	{
	  if(sElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0)
	  {
	    sElmId = sElm->neighbours[0];
	  }
	  else
	  {
	    nVx.vtX = (sVx1.vtX + sVx2.vtX) / 2.0;
	    nVx.vtY = (sVx1.vtY + sVx2.vtY) / 2.0;
	    bndFlg = 1;
	  }
	}
	else if((sElmArea1 = WlzGeomTriangleSnArea2(sVx2, sVx0, nVx)) < 0.0)
	{
	  if(sElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1)
	  {
	    sElmId = sElm->neighbours[1];
	  }
	  else
	  {
	    nVx.vtX = (sVx2.vtX + sVx0.vtX) / 2.0;
	    nVx.vtY = (sVx2.vtY + sVx0.vtY) / 2.0;
	    bndFlg = 1;
	  }
	}
	else if(sElmArea - sElmArea0 - sElmArea1 < 0.0)
	{
	  if(sElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2)
	  {
	    sElmId = sElm->neighbours[2];
	  }
	  else
	  {
	    nVx.vtX = (sVx0.vtX + sVx1.vtX) / 2.0;
	    nVx.vtY = (sVx0.vtY + sVx1.vtY) / 2.0;
	    bndFlg = 1;
	  }
	}
	else
	{
	  fndFlg = 1;
	}
	if((fndFlg == 0) && (bndFlg == 0))
	{
	  sElm = mesh->elements + sElmId;
	  sVx0 = (nodes + sElm->nodes[0])->position;
	  sVx1 = (nodes + sElm->nodes[1])->position;
	  sVx2 = (nodes + sElm->nodes[2])->position;
	}
      }
    } while((fndFlg == 0) && (bndFlg == 0) && (errNum == WLZ_ERR_NONE) &&
            (countOut-- > 0));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(countOut <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find elements in conflict with the new node vertex. */
    eCnfQVec.count = 0;
    errNum = WlzMeshQueConflictElem(&eCnfQVec, mesh, sElmId, sElmId, nVx);
    if(errNum == WLZ_ERR_NONE)
    {
      /* Replace conflicting and add new elements. */
      errNum = WlzMeshElemReplaceN(mesh, eCnfQVec.vector,
      				   eCnfQVec.count, nVx,
				   WLZ_MESH_NODE_FLAGS_NONE);
    }
  }
  if(eCnfQVec.vector)
  {
    AlcFree(eCnfQVec.vector);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds a node to the given mesh, the new node must lie inside
*		the existing mesh.
* \param	mesh			Given mesh transform.
* \param	startElm		If >= 0, the index of the element
*					from which to start the search.
* \param	newVx			New node vertex.
* \param	nodeFlags		Node flags to set (eg source).
*/
WlzErrorNum	WlzMeshNodeAdd(WlzMeshTransform *mesh, int startElm,
			       WlzDVertex2 newVx, unsigned int nodeFlags)
{
  int		eId0,
		eId1,
		nId0,
		extFlg;
  WlzMeshElem	*elm;
  WlzMeshNode	*nod;
  WlzMeshIntVec	eCnfQVec;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  eCnfQVec.vector = NULL;
  eCnfQVec.size = 0;
  eCnfQVec.minSize = 8;
  eCnfQVec.mulSize = 2;
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_TRANSFORM_2D_MESH)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(startElm >= mesh->nElem)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(mesh->nElem > 0)
  {
    /* Find vertex in the mesh. */
    eId0 = WlzMeshElemFindVx(mesh, newVx, startElm, &extFlg, &errNum);
    if(eId0 < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(extFlg)
      {
	/* Node already exists find which node of element and modify
	 * it's flags. */
	elm = mesh->elements + eId0;
	if((nId0 = WlzMeshElemNodeIdxFromVx(mesh, elm, newVx)) > 0)
	{
	  nod = mesh->nodes + *(elm->nodes + nId0);
	  nod->flags = nodeFlags;
	}
	else
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
      else
      {
	eId1 = eId0;
	/* Find elements in conflict with the new node vertex. */
	eCnfQVec.count = 0;
	errNum = WlzMeshQueConflictElem(&eCnfQVec, mesh, eId0, eId1, newVx);
	if(errNum == WLZ_ERR_NONE)
	{
	  /* Replace conflicting and add new elements. */
	  errNum = WlzMeshElemReplaceN(mesh, eCnfQVec.vector,
	  			       eCnfQVec.count, newVx,
				       WLZ_MESH_NODE_FLAGS_NONE);
	}
      }
    }
  }
  if(eCnfQVec.vector)
  {
    AlcFree(eCnfQVec.vector);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Deletes nodes, specified by their index, from the given mesh.
* \param	mesh			Given mesh transform.
* \param	startElm		If >= 0, the index of the element
*					from which to start the search an
*					element using the node.
* \param	nodIdP			Vector of node indicies.
* \param	nNod			Number of indicies in vector.
*/
WlzErrorNum	WlzMeshNodeDelIdx(WlzMeshTransform *mesh, int startElm,
			          int *nodIdP, int nNod)
{
  int		nId0,
  		nId1;
  WlzDVertex2	nodVx;
  WlzMeshIntVec	elmVec,
  		nodVec;
  WlzMeshEarList earList;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_TRANSFORM_2D_MESH)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    nId0 = 0;
    WlzMeshNodeDelInit(&elmVec, &nodVec, &earList);
    while((nId0 < nNod) && (errNum == WLZ_ERR_NONE))
    {
      nId1 = *(nodIdP + nId0);
      if((nId1 < 0) || (nId1 >= mesh->nNodes))
      {
	errNum = WLZ_ERR_PARAM_DATA;
      }
      else
      {
	if((startElm < 0) || (startElm >= mesh->nElem))
	{
	  startElm = 0;
	}
	nodVx = (mesh->nodes + nId1)->position;
	errNum = WlzMeshNodeDel(mesh, &elmVec, &nodVec, &earList, startElm,
				nodVx);
      }
      ++nId0;
    }
    WlzMeshNodeDelFree(&elmVec, &nodVec, &earList);
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Initialize the element vector, node vector and ear list
*		ready for node deletion.
* \param	elmVec			Element vector.
* \param	nodVec			Node vector.
* \param	earList			Ear list.
*/
static void	WlzMeshNodeDelInit(WlzMeshIntVec *elmVec,
				   WlzMeshIntVec *nodVec,
				   WlzMeshEarList *earList)
{
  elmVec->count = 0;
  elmVec->size = 0;
  elmVec->vector = NULL;
  elmVec->minSize = 256;
  elmVec->mulSize = 2;
  nodVec->count = 0;
  nodVec->size = 0;
  nodVec->vector = NULL;
  nodVec->minSize = 256;
  nodVec->mulSize = 2;
  earList->maxEars = 0;
  earList->pool = NULL;
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Free storage allocated for the element vector, node vector
*		and ear list after node deletion.
* \param	elmVec			Element vector.
* \param	nodVec			Node vector.
* \param	earList			Ear list.
*/
static void	WlzMeshNodeDelFree(WlzMeshIntVec *elmVec, 
				   WlzMeshIntVec *nodVec, 
				   WlzMeshEarList *earList)
{
  if(elmVec->size && elmVec->vector)
  {
    AlcFree(elmVec->vector);
  }
  if(nodVec->size && nodVec->vector)
  {
    AlcFree(nodVec->vector);
  }
  if(earList->maxEars && earList->pool)
  {
    AlcFree(earList->pool);
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Deletes a node, specified by it's position, from the
*		given mesh.
* \param	mesh			Given mesh transform.
* \param	elmVec			Element vector.
* \param	nodVec			Node vector.
* \param	earList			Ear list.
* \param	startElm		Index of the element from which to
*					start the search for an element using
*					the node.
* \param	nodVx			Node position.
*/
static WlzErrorNum WlzMeshNodeDel(WlzMeshTransform *mesh,
				  WlzMeshIntVec *elmVec,
				  WlzMeshIntVec *nodVec, 
				  WlzMeshEarList *earList,
				  int startElm, WlzDVertex2 nodVx)
{
  int		eId0,
		eId1,
		eId2,
		eNodId0,
		delNodId,
  		extFlg;
  WlzMeshEar	*ear0,
  		*ear1,
  		*ear2;
  WlzMeshElem	*elm0,
  		*elm1;
  WlzMeshNode	*dNod;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  eId0 = WlzMeshElemFindVx(mesh, nodVx, startElm, &extFlg, &errNum);
  if((eId0 < 0) || (extFlg == 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    /* Find index of node within in the element. */
    if((eNodId0 = WlzMeshElemNodeIdxFromVx(mesh, mesh->elements + eId0,
    					   nodVx)) < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Build a vector of the id's of all the mesh elements which include
     * the node to be deleted. Where the element and node vectors are in
     * ccw order around the node to be deleted. All the elements are marked
     * as zombies. */
    errNum = WlzMeshNodeDelVecBuild(elmVec, nodVec, mesh, eId0, eNodId0);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Build a circular list of possible ears and assign each ear a
       power value. */
    delNodId = *((mesh->elements + eId0)->nodes + eNodId0);
    errNum = WlzMeshEarsCreate(earList, elmVec, nodVec, mesh, delNodId);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Make sure that each of the elements in the element vector is
     * unlinked from its neighbours and marked as a zombie. */
    for(eId0 = 0; eId0 < elmVec->count; ++eId0)
    {
      elm0 = mesh->elements + *(elmVec->vector + eId0);
      for(eId1 = 0; eId1 < 3; ++eId1)
      {
        if(elm0->flags & nbrFlgTbl[eId1])
	{
	  elm1 = mesh->elements + elm0->neighbours[eId1];
	  for(eId2 = 0; eId2 < 3; ++eId2)
	  {
	    if((elm1->flags & nbrFlgTbl[eId2]) &&
	       (elm1->neighbours[eId2] == elm0->idx))
	    {
	      elm1->flags &=  ~(nbrFlgTbl[eId2]);
	    }
	  }
	}
	elm0->flags &= ~(nbrFlgTbl[eId1]);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* While the in-list has more than three entries: Pull off the minimum
     * power entry from it and add it to the out stack, modifying the
     * entry's previous and next to reflect its removal. */
    while((earList->inEars > 3) &&
          ((ear0 = WlzMeshEarGetMinPower(earList))->snArea2 > 
	   WLZ_MESH_ELEM_AREA_TOLERANCE * 2))
    {
      --(earList->inEars);
      ear1 = ear0->prev;
      ear2 = ear0->next;
      ear1->nodes[2] = ear0->nodes[2];
      ear1->next = ear2;
      ear2->nodes[0] = ear0->nodes[0];
      ear2->prev = ear1;
      WlzMeshEarPowerSet(mesh, ear1, delNodId);
      WlzMeshEarPowerSet(mesh, ear2, delNodId);
      /* Add this ear to the mesh recycling the last of the mesh elements
       * in the element vector. */
      ear0->prev = ear0->next = NULL;
      /* Make sure that the element for recycling is available as the
       * element vector may have duplicate entries! */
      do
      {
	eId0 = *(elmVec->vector + --(elmVec->count));
	elm0 = mesh->elements + eId0;
      } while(((elm0->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0) &&
              (elmVec->count > 0));
      elm0->flags = ear0->flags;
      for(eId1 = 0; eId1 < 3; ++eId1)
      {
	elm0->nodes[eId1] = ear0->nodes[eId1];
	if(ear0->flags & nbrFlgTbl[eId1])
	{
	  elm0->neighbours[eId1] = ear0->neighbours[eId1];
	  elm1 = mesh->elements + elm0->neighbours[eId1];
	  elm1->neighbours[ear0->selfNeighbour[eId1]] = elm0->idx;
	  elm1->flags |= nbrFlgTbl[ear0->selfNeighbour[eId1]];
	}
      }
      ear1->flags |=  WLZ_MESH_ELEM_FLAGS_NBR_0;
      ear1->neighbours[0] = elm0->idx;
      ear1->selfNeighbour[0] = 1;
      ear2->flags |=  WLZ_MESH_ELEM_FLAGS_NBR_2;
      ear2->neighbours[2] = elm0->idx;
      ear2->selfNeighbour[2] = 1;
    }
    if((ear0 = WlzMeshEarGetMinPower(earList))->snArea2 > 
       WLZ_MESH_ELEM_AREA_TOLERANCE * 2)
    {
      /* The three ears left are rotations around the same three nodes,
       * but have incomplete neighbour flags/ids. Add thisear to the mesh,
       * recycling the last of the mesh elements in the element vector. */
      ear0->prev = ear0->next = NULL;
      /* Add this ear to the mesh recycling the last of the mesh elements
       * in the element vector. */
      eId0 = *(elmVec->vector + --(elmVec->count));
      elm0 = mesh->elements + eId0;
      elm0->flags = ear0->flags;
      for(eId1 = 0; eId1 < 3; ++eId1)
      {
        if((ear0->flags & nbrFlgTbl[eId1]) != 0)
	{
	  elm0->nodes[eId1] = ear0->nodes[eId1];
	  if((ear0->flags & nbrFlgTbl[eId1]) != 0)
	  {
	    elm0->flags |= nbrFlgTbl[eId1];
	    elm0->neighbours[eId1] = ear0->neighbours[eId1];
	  }
	  elm1 = mesh->elements + elm0->neighbours[eId1];
	  elm1->neighbours[ear0->selfNeighbour[eId1]] = elm0->idx;
	  elm1->flags |= nbrFlgTbl[ear0->selfNeighbour[eId1]];
	}
	else
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
    }
    dNod = mesh->nodes + delNodId;
    dNod->flags = WLZ_MESH_NODE_FLAGS_ZOMBIE;
  }
#ifdef WLZ_MESH_DEBUG
  errNum = WlzMeshTransformVerify(mesh, 0, NULL, NULL);
#endif /* WLZ_MESH_DEBUG */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Searches the mesh for the element which contains the given
* 		vertex. Walks from the given element in the direction of the
* 		given vertex until the element which encloses the given vertex
* 		is found. If the walk oscillates back and forth between two
* 		elements then either may be returned as the enclosing element.
* \param	mesh			Given mesh transform.
* \param	gVx			Given vertex.
* \param	elmId			Source and destination element index
*					pointer.
* \param	foundFlg		Destination pointer for element found
*					flag.
* \param	existsFlg		Destination pointer for node exists
*					flag.
*/
static WlzErrorNum WlzMeshElemFindVxWalk(WlzMeshTransform *mesh,
					 WlzDVertex2 gVx,
					 int *elmId, int *foundFlg,
					 int *existsFlg)
{
  int		eId,
		bndFlg = 0,
		extFlg = 0,
		fndFlg = 0,
		countOut;
  double	tD0,
		eArea,
		pArea0,
		pArea1;
  WlzMeshElem	*elm = NULL,
		*lastElm = NULL;
  WlzDVertex2	eVx0,
		eVx1,
		eVx2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  eId = *elmId;
  while((eId < mesh->nElem) &&
        ((elm = mesh->elements + eId)->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE))
  {
    ++eId;
  }
  if(eId >= mesh->nElem)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  countOut = mesh->nNodes;
  while((fndFlg == 0) && (bndFlg == 0) && (errNum == WLZ_ERR_NONE) &&
        (countOut-- > 0))
  {
    elm = mesh->elements + eId;
    eVx0 = (mesh->nodes + elm->nodes[0])->position;
    eVx1 = (mesh->nodes + elm->nodes[1])->position;
    eVx2 = (mesh->nodes + elm->nodes[2])->position;
    if((eArea = WlzGeomTriangleSnArea2(eVx0, eVx1,
				       eVx2)) < WLZ_MESH_TOLERANCE_SQ)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      if(((((tD0 = (eVx0.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	  (((tD0 = (eVx0.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)) ||
	 ((((tD0 = (eVx1.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	  (((tD0 = (eVx1.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)) ||
	 ((((tD0 = (eVx2.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	  (((tD0 = (eVx2.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)))
      {
	fndFlg = 1;
	extFlg = 1;
      }
      else if((pArea0 = WlzGeomTriangleSnArea2(eVx1, eVx2, gVx)) < 0.0)
      {
	if(elm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0)
	{
	  eId = elm->neighbours[0];
	  if(elm && lastElm && (eId == lastElm->idx))
	  {
	    fndFlg = 1;
	  }
	}
	else
	{
	  bndFlg = 1;
	}
      }
      else if((pArea1 = WlzGeomTriangleSnArea2(eVx2, eVx0, gVx)) < 0.0)
      {
	if(elm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1)
	{
	  eId = elm->neighbours[1];
	  if(elm && lastElm && (eId == lastElm->idx))
	  {
	    fndFlg = 1;
	  }
	}
	else
	{
	  bndFlg = 1;
	}
      }
      else if(eArea - pArea0 - pArea1 < 0.0)
      {
	if(elm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2)
	{
	  eId = elm->neighbours[2];
	  if(elm && lastElm && (eId == lastElm->idx))
	  {
	    fndFlg = 1;
	  }
	}
	else
	{
	  bndFlg = 1;
	}
      }
      else
      {
	fndFlg = 1;
      }
    }
    lastElm = elm;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(countOut <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  *elmId = eId,
  *foundFlg = fndFlg;
  *existsFlg = extFlg;
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Finds the mesh element which contains the given vertex by
*		brute force.
* \param	mesh			Given mesh transform.
* \param	gVx			Given vertex.
* \param	elmId			Source and destination element index
*					pointer.
* \param	foundFlg		Destination pointer for element
*					found flag.
* \param	existsFlg		Destination pointer for node exists
*					flag.
*/
static void	 WlzMeshElemFindVxForce(WlzMeshTransform *mesh,
					  WlzDVertex2 gVx,
					  int *elmId, int *foundFlg,
					  int *existsFlg)
{
  int		tstId = 0,
		extFlg = 0,
		fndId = -1;
  double	tD0;
  int		*eNdP;
  WlzMeshElem	*elm;
  WlzDVertex2	eVx0,
		eVx1,
		eVx2;


  elm = mesh->elements;
  while((tstId < mesh->nElem) && (fndId < 0))
  {
    if((elm->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
    {
      eNdP = elm->nodes;
      eVx0 = (mesh->nodes + *(eNdP + 0))->position;
      eVx1 = (mesh->nodes + *(eNdP + 1))->position;
      eVx2 = (mesh->nodes + *(eNdP + 2))->position;
      if(WlzGeomVxInTriangle(eVx0, eVx1, eVx2, gVx))
      {
	fndId = tstId;
	if(((((tD0 = (eVx0.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	    (((tD0 = (eVx0.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)) ||
	   ((((tD0 = (eVx1.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	    (((tD0 = (eVx1.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)) ||
	   ((((tD0 = (eVx2.vtX - gVx.vtX)) * tD0) < WLZ_MESH_TOLERANCE_SQ) &&
	    (((tD0 = (eVx2.vtY - gVx.vtY)) * tD0) < WLZ_MESH_TOLERANCE_SQ)))
	{
	  extFlg = 1;
	}
      }
    }
    ++elm;
    ++tstId;
  }
  *elmId = fndId;
  *existsFlg = extFlg;
  *foundFlg = fndId != -1;
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Recursively builds a que of elements in conflict with the
*		given new node vertex, all elements in conflict being marked
*		using a flag. An element is in conflict with a new node vertex
*		if the new node lies within the circumcircle of the element.
* \param	eCnfQVec		The mesh element que to build.
* \param	mesh			Given mesh transform.
* \param	srcElmId		Index of element which contains the
*					vertex.
* \param	elmId			Index of element to test and possibly
*					add to the queue.
* \param	newVx			New node vertex.
*/
static WlzErrorNum WlzMeshQueConflictElem(WlzMeshIntVec *eCnfQVec,
					  WlzMeshTransform *mesh,
					  int srcElmId, int elmId,
					  WlzDVertex2 newVx)
{
  int		nbrId,
		nbrElmId;
  WlzMeshElem	*elm0,
		*elm1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  if((srcElmId == elmId) || WlzMeshElemConflict(mesh, elmId, newVx))
  {
    if(errNum == WLZ_ERR_NONE)
    {
      /* Put this element in the queue and mark it a zombie. */
      elm0 = mesh->elements + elmId;
      if((elm0->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
      {
	elm0->flags |= WLZ_MESH_ELEM_FLAGS_ZOMBIE;
	errNum = WlzMeshAddToIntVec(eCnfQVec, elmId);
	nbrId = 0;
	/* Test and possibly enqueue each neighbour. */
	while((nbrId < 3) && (errNum == WLZ_ERR_NONE))
	{
	  if(elm0->flags & nbrFlgTbl[nbrId])
	  {
	    nbrElmId = elm0->neighbours[nbrId];
	    elm1 = mesh->elements + nbrElmId;
	    if((elm1->flags & WLZ_MESH_ELEM_FLAGS_ZOMBIE) == 0)
	    {
	      errNum = WlzMeshQueConflictElem(eCnfQVec, mesh,
	      				      srcElmId, nbrElmId, newVx);
	    }
	  }
	  ++nbrId;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Non-zero if element is in conflict with the new node vertex.
* \ingroup	WlzTransform
* \brief	Tests the given element for an in-circumcircle conflict with
* 		the given new node vertex.
* \param	mesh			Given mesh transform.
* \param	elmId			Index of element in mesh.
* \param	newVx			New node vertex.
*/
static int	WlzMeshElemConflict(WlzMeshTransform *mesh, int elmId,
				    WlzDVertex2 newVx)
{
  int		conflict;
  int		*eNIP;
  WlzMeshNode	*nod;

  nod = mesh->nodes;
  eNIP = (mesh->elements + elmId)->nodes;
  conflict = WlzGeomInTriangleCircumcircle((nod + *(eNIP + 0))->position,
					   (nod + *(eNIP + 1))->position,
					   (nod + *(eNIP + 2))->position,
					   newVx);
  return(conflict);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Replaces those mesh elements that are queued and flaged as
* 		zombies with new mesh elements which include the given new
* 		node vertex.
* \param	mesh			Given mesh transform.
* \param	eCnfQ			The mesh element que to build.
* \param	qCnt			Number of elements in queue.
* \param	newVx			New node vertex.
* \param	nodeFlags		Node flags to set (eg source).
*/
static WlzErrorNum WlzMeshElemReplaceN(WlzMeshTransform *mesh, int *eCnfQ,
				       int qCnt, WlzDVertex2 newVx,
				       unsigned int nodeFlags)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(qCnt == 1)
  {
    errNum = WlzMeshElemReplace1(mesh, *eCnfQ, newVx, nodeFlags);
  }
  else
  {
    errNum = WlzMeshElemReplaceNWithN(mesh, eCnfQ, qCnt, newVx, nodeFlags);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Replaces a single mesh element with 1 element (no action),
*		2 elements (1 new and 1 recycled) or 3 elements (2 new and 1
*		recycled) using the new node to split the element.
* \param	mesh			Given mesh transform.
* \param	eId			Index of the element containing the
*					new node vertex.
* \param	newVx			New node vertex.
* \param 	nodeFlags		Node flags to set (eg source).
*/
static WlzErrorNum WlzMeshElemReplace1(WlzMeshTransform *mesh,
				       int eId, WlzDVertex2 newVx,
				       unsigned int nodeFlags)
{
  int		invId,
		withCnt = 0;
  WlzMeshElem	*elm;
  WlzDVertex2	vx0,
		vx1,
		vx2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  elm = mesh->elements + eId;
  vx0 = (mesh->nodes + elm->nodes[0])->position;
  vx1 = (mesh->nodes + elm->nodes[1])->position;
  vx2 = (mesh->nodes + elm->nodes[2])->position;
  if(WlzGeomTriangleSnArea2(vx1, vx2, newVx) < WLZ_MESH_TOLERANCE_SQ)
  {
    invId = 0;
  }
  else
  {
    ++withCnt;
  }
  if(WlzGeomTriangleSnArea2(vx2, vx0, newVx) < WLZ_MESH_TOLERANCE_SQ)
  {
    invId = 1;
  }
  else
  {
    ++withCnt;
  }
  if(WlzGeomTriangleSnArea2(vx0, vx1, newVx) < WLZ_MESH_TOLERANCE_SQ)
  {
    invId = 2;
  }
  else
  {
    ++withCnt;
  }
  errNum = WlzMeshExpand(mesh, mesh->nElem + withCnt, mesh->nNodes + 1);
  if(errNum == WLZ_ERR_NONE)
  {
    if(withCnt == 3)
    {
      WlzMeshElemReplace1With3(mesh, eId, newVx, nodeFlags);
    }
    else if(withCnt == 2)
    {
      WlzMeshElemReplace1With2(mesh, eId, invId, newVx, nodeFlags);
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Replaces a single mesh element with another mesh element.
* \param	mesh			Given mesh transform.
* \param	eId0			Index of the element to be replaced in
*					the mesh.
* \param	eId1			Index of the element to replace
*					element with index eId0.
*/
static void	WlzMeshElemReplace1With1(WlzMeshTransform *mesh,
					 int eId0, int eId1)
{
  int		nId;
  WlzMeshElem	*elm0,
		*elm1,
		*nElm;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  elm0 = mesh->elements + eId0;
  elm1 = mesh->elements + eId1;
  *elm1 = *elm0;
  elm1->idx = eId1;
  for(nId = 0; nId < 3; ++nId)
  {
    if((elm1->flags & nbrFlgTbl[nId]) != 0)
    {
      nElm = mesh->elements + elm1->neighbours[nId];
      if(nElm->neighbours[0] == eId0)
      {
	nElm->neighbours[0] = eId1;
      }
      else if(nElm->neighbours[1] == eId0)
      {
	nElm->neighbours[1] = eId1;
      }
      else if(nElm->neighbours[2] == eId0)
      {
	nElm->neighbours[2] = eId1;
      }
    }
  }
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Replaces a single mesh element with 2 mesh elements,
*		(1 new and 1 recycled) which include the given new node
*		vertex.
* \param	mesh			Given mesh transform.
* \param	eId			Index of the element containing
*					the new node vertex.
* \param	nod0			Which of the existing elements nodes is
* 					shared between the two new elements.
* \param	newVx			New node vertex.
* \param	nodeFlags		Node flags to set (eg source).
*/
static void	WlzMeshElemReplace1With2(WlzMeshTransform *mesh,
					 int eId, int nod0,
					 WlzDVertex2 newVx,
					 unsigned int nodeFlags)
{
  int		nod1,
		nod2,
		newNodId,
		nElem,
		nNodes;
  WlzMeshNode	*nNod;
  WlzMeshElem	*nElm;
  WlzMeshElem	*elm[2];
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  newNodId = mesh->nNodes;
  nElem = mesh->nElem + 1;
  nNodes = mesh->nNodes + 1;
  nod1 = (nod0 + 1) % 3;
  nod2 = (nod0 + 2) % 3;
  /* Mesh has been expanded already for nElem elements and nNodes nodes. */
  /* Set element indicies. */
  elm[0] = mesh->elements + eId;
  nNod = mesh->nodes + newNodId;
  nNod->flags = nodeFlags;
  nNod->position = newVx;
  nNod->displacement.vtX = 0.0;
  nNod->displacement.vtY = 0.0;
  elm[1] = mesh->elements + mesh->nElem;
  elm[1]->type = WLZ_MESH_ELEM_TRILINEAR;
  elm[1]->idx = mesh->nElem;
  elm[1]->flags = nbrFlgTbl[nod2];
  elm[1]->nodes[nod0] = elm[0]->nodes[nod0];
  elm[1]->nodes[nod1] = newNodId;
  elm[1]->nodes[nod2] = elm[0]->nodes[nod2];
  elm[1]->neighbours[nod2] = elm[0]->idx;
  if((elm[0]->flags & nbrFlgTbl[nod1]) != 0)
  {
    elm[1]->flags |= nbrFlgTbl[nod1];
    elm[1]->neighbours[nod1] = elm[0]->neighbours[nod1];
    nElm = mesh->elements + elm[0]->neighbours[nod1];
    if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0) &&
       (nElm->neighbours[0] == elm[0]->idx))
    {
      nElm->neighbours[0] = elm[1]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0) &&
	    (nElm->neighbours[1] == elm[0]->idx))
    {
      nElm->neighbours[1] = elm[1]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2) != 0) &&
	    (nElm->neighbours[2] == elm[0]->idx))
    {
      nElm->neighbours[2] = elm[1]->idx;
    }
  }
  elm[0]->flags &= nbrFlgTbl[nod2];
  elm[0]->flags |= nbrFlgTbl[nod1];
  elm[0]->nodes[nod2] = newNodId;
  elm[0]->neighbours[nod1] = elm[1]->idx;
  elm[0]->strainU[0] = elm[0]->strainU[1] = elm[0]->strainU[2] = 0.0;
  elm[0]->strainA[0] = elm[0]->strainA[1] = elm[0]->strainA[2] = 0.0;
  elm[1]->strainU[0] = elm[1]->strainU[1] = elm[1]->strainU[2] = 0.0;
  elm[1]->strainA[0] = elm[1]->strainA[1] = elm[1]->strainA[2] = 0.0;
  mesh->nNodes = nNodes;
  mesh->nElem = nElem;
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Replaces a single mesh element with 3 mesh elements,
*		(2 new and 1 recycled) which include the given new node
*		vertex.
* \param	mesh			Given mesh transform.
* \param	eId			Index of the element containing the
*					new node vertex.
* \param	newVx			New node vertex.
* \param	nodeFlags		Node flags to set (eg source).
*/
static void	 WlzMeshElemReplace1With3(WlzMeshTransform *mesh,
					  int eId, WlzDVertex2 newVx,
					  unsigned int nodeFlags)
{
  int		nId,
		newNodId,
		nElem,
		nNodes;
  WlzMeshElem	*nElm;
  WlzMeshNode	*nNod;
  WlzMeshElem	*elm[3];

  newNodId = mesh->nNodes;
  nElem = mesh->nElem + 2;
  nNodes = mesh->nNodes + 1;
  /* Mesh has been expanded already for nElem elements and nNodes nodes. */
  /* Set element indicies. */
  elm[0] = mesh->elements + eId;
  nNod = mesh->nodes + newNodId; 
  nNod->flags = nodeFlags;
  nNod->position = newVx;
  nNod->displacement.vtX = 0.0;
  nNod->displacement.vtY = 0.0; 
  elm[1] = mesh->elements + mesh->nElem;
  elm[1]->type = WLZ_MESH_ELEM_TRILINEAR;
  elm[1]->idx = mesh->nElem;
  elm[2] = mesh->elements + mesh->nElem + 1;
  elm[2]->type = WLZ_MESH_ELEM_TRILINEAR;
  elm[2]->idx = mesh->nElem + 1;
  /* Create 1st new element. */
  elm[1]->flags = WLZ_MESH_ELEM_FLAGS_NBR_1 |
		  WLZ_MESH_ELEM_FLAGS_NBR_2;
  elm[1]->nodes[0] = newNodId;
  elm[1]->nodes[1] = elm[0]->nodes[1];
  elm[1]->nodes[2] = elm[0]->nodes[2];
  if((elm[0]->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0)
  {
    nId = elm[0]->neighbours[0];
    elm[1]->flags |= WLZ_MESH_ELEM_FLAGS_NBR_0;
    elm[1]->neighbours[0] = nId;
    nElm = mesh->elements + nId;
    if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0) &&
       (nElm->neighbours[0] == elm[0]->idx))
    {
      nElm->neighbours[0] = elm[1]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0) &&
	    (nElm->neighbours[1] == elm[0]->idx))
    {
      nElm->neighbours[1] = elm[1]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2) != 0) &&
	    (nElm->neighbours[2] == elm[0]->idx))
    {
      nElm->neighbours[2] = elm[1]->idx;
    }
  }
  elm[1]->neighbours[1] = elm[2]->idx;
  elm[1]->neighbours[2] = elm[0]->idx;
  elm[1]->strainU[0] = elm[1]->strainU[1] = elm[1]->strainU[2] = 0.0;
  elm[1]->strainA[0] = elm[1]->strainA[1] = elm[1]->strainA[2] = 0.0;
  /* Create 2nd new element. */
  elm[2]->flags = WLZ_MESH_ELEM_FLAGS_NBR_0 |
		  WLZ_MESH_ELEM_FLAGS_NBR_2;
  elm[2]->nodes[0] = elm[0]->nodes[0];
  elm[2]->nodes[1] = newNodId;
  elm[2]->nodes[2] = elm[0]->nodes[2];
  if((elm[0]->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0)
  {
    nId = elm[0]->neighbours[1];
    elm[2]->flags |= WLZ_MESH_ELEM_FLAGS_NBR_1;
    elm[2]->neighbours[1] = nId;
    nElm = mesh->elements + nId;
    if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0) &&
       (nElm->neighbours[0] == elm[0]->idx))
    {
      nElm->neighbours[0] = elm[2]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0) &&
	    (nElm->neighbours[1] == elm[0]->idx))
    {
      nElm->neighbours[1] = elm[2]->idx;
    }
    else if(((nElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2) != 0) &&
	    (nElm->neighbours[2] == elm[0]->idx))
    {
      nElm->neighbours[2] = elm[2]->idx;
    }
  }
  elm[2]->neighbours[0] = elm[1]->idx;
  elm[2]->neighbours[2] = elm[0]->idx;
  elm[2]->strainU[0] = elm[2]->strainU[1] = elm[2]->strainU[2] = 0.0;
  elm[2]->strainA[0] = elm[2]->strainA[1] = elm[2]->strainA[2] = 0.0;
  elm[0]->flags = WLZ_MESH_ELEM_FLAGS_NBR_0 |
		  WLZ_MESH_ELEM_FLAGS_NBR_1 |
		  (elm[0]->flags & WLZ_MESH_ELEM_FLAGS_NBR_2);
  /* elm[0]->nodes[0] = elm[0]->nodes[0]; */
  /* elm[0]->nodes[1] = elm[0]->nodes[1]; */
  elm[0]->nodes[2] = newNodId;
  elm[0]->neighbours[0] = elm[1]->idx;
  elm[0]->neighbours[1] = elm[2]->idx;
  /* elm[0]->neighbours[2] = elm[2]->neighbours[2]; */
  elm[0]->strainU[0] = elm[0]->strainU[1] = elm[0]->strainU[2] = 0.0;
  elm[0]->strainA[0] = elm[0]->strainA[1] = elm[0]->strainA[2] = 0.0;
  mesh->nNodes = nNodes;
  mesh->nElem = nElem;
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Replaces all queued mesh elements (which always have a
*		convex hull) with elements that use the new node vertex.
*		This function is only ever called with more than one mesh
*		element enqueued. Computing the number of nodes in the enqueued
*		elements is easy because the region is always convex with all
*		it's nodes on it's boundary. Every triangulation of a polygon
*		with N nodes has N - 2 elements.
* \param	mesh			Given mesh transform.
* \param	zElmIdVec		The vector (queue) of zombie mesh
*					elements to replace.
* \param	zElmCnt			Number of zombie elements.
* \param	newVx			New node vertex.
* \param	nodeFlags		Node flags to set (eg source).
*/
static WlzErrorNum WlzMeshElemReplaceNWithN(WlzMeshTransform *mesh,
					    int *zElmIdVec,
					    int zElmCnt, WlzDVertex2 newVx,
					    unsigned int nodeFlags)
{
  int		bId,
		nId,
		rId,
		wId,
		zId,
		zNId,
		fndFlg,
		wElmCnt,
		mElmCnt,
		mNodCnt,
		newNodId;
  WlzMeshElem	*nElm,
		*rElm,
		*wElm,
		*zElm,
		*zNElm;
  WlzMeshNode	*nNod;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  /* Make sure that there's room in the mesh for some working elements
   * and two new elements. Compute the number of working elements required
   (wElmCnt), the total number of elements required in the mesh (mElmCnt)
   and the total number of nodes required in the mesh. */
  wElmCnt = zElmCnt + 2;
  mElmCnt = mesh->nElem + wElmCnt;
  mNodCnt = mesh->nNodes + 1;
  if((mesh->maxElem < mElmCnt) || (mesh->maxNodes < mNodCnt))
  {
    errNum = WlzMeshExpand(mesh, mElmCnt, mNodCnt);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Add the new node vertex. */
    newNodId = mesh->nNodes;
    nNod = mesh->nodes + newNodId; 
    nNod->flags = nodeFlags;
    nNod->position = newVx;
    nNod->displacement.vtX = 0.0;
    nNod->displacement.vtY = 0.0; 
    /* Find a zombie element which has an edge not shared with another
     * zombie. At the end of the search the zombie element zElm has
     * index zId and either a non-zombie neighbour or no neighbour for
     * for neighbour zNId. */
    zId = 0;
    fndFlg = 0;
    while((fndFlg == 0) && (zId < zElmCnt))
    {
      zElm = mesh->elements + *(zElmIdVec + zId);
      zNId = 0;
      while((zNId < 3) &&
	    ((zElm->flags & nbrFlgTbl[zNId]) != 0) &&
	    (((zNElm = mesh->elements + zElm->neighbours[zNId])->flags &
	      WLZ_MESH_ELEM_FLAGS_ZOMBIE) != 0))
      {
	++zNId;
      }
      if(zNId < 3)
      {
	fndFlg = 1;
      }
      else
      {
	++zId;
      }
    }
    /* Walk CCW around the zombie elements filling in the work elements. */
    wId = 0;
    while(wId < wElmCnt)
    {
      /* Set up the wId'th working element. */
      wElm = mesh->elements + mesh->nElem + wId;
      wElm->type = WLZ_MESH_ELEM_TRILINEAR;
      wElm->idx = mesh->nElem + wId;
      wElm->flags = WLZ_MESH_ELEM_FLAGS_NBR_1 | WLZ_MESH_ELEM_FLAGS_NBR_2;
      if((zElm->flags & nbrFlgTbl[zNId]) != 0)
      {
	nId = zElm->neighbours[zNId];
	wElm->flags |= WLZ_MESH_ELEM_FLAGS_NBR_0;
	wElm->neighbours[0] = nId;
	zNElm = mesh->elements + nId;
	if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0) &&
	   (zNElm->neighbours[0] == zElm->idx))
	{
	  zNElm->neighbours[0] = wElm->idx;
	}
	else if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0) &&
	        (zNElm->neighbours[1] == zElm->idx))
	{
	  zNElm->neighbours[1] = wElm->idx;
	}
	else if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2) != 0) &&
	        (zNElm->neighbours[2] == zElm->idx))
	{
	  zNElm->neighbours[2] = wElm->idx;
	}
      }
      wElm->nodes[0] = newNodId;
      wElm->nodes[1] = zElm->nodes[(zNId + 1) % 3];
      wElm->nodes[2] = zElm->nodes[(zNId + 2) % 3];
      if(++wId < wElmCnt)
      {
	zNId = (zNId + 1) % 3;      /* Try this zombies next CCW neighbour */
	while(((zElm->flags & nbrFlgTbl[zNId]) != 0) &&
	      (((zNElm = mesh->elements + zElm->neighbours[zNId])->flags &
		WLZ_MESH_ELEM_FLAGS_ZOMBIE) != 0))
	{
	  /* This zombie element has a neighbour on zNId and it's a zombie,
	   * so walk to the next neighbour which has a non-zombie neighbour
	   * or no neighbour. */
	  nId = zElm->nodes[(zNId + 1) % 3];
	  if(zNElm->nodes[0] == nId)
	  {
	    zNId = 2;
	  }
	  else if(zNElm->nodes[1] == nId)
	  {
	    zNId = 0;
	  }
	  else /* zNElm->nodes[2] == nId */
	  {
	    zNId = 1;
	  }
	  zElm = zNElm;
	}
      }
    }
    /* If the new node lies on an element edge segment then one of the working
     * elements may have a zero area: Check for this. */
    bId = -1;
    wId = 0;
    while((wId < wElmCnt) && (bId < 0))
    {
      /* Check the wId'th working element. */
      wElm = mesh->elements + mesh->nElem + wId;
      if(WlzGeomTriangleSnArea2((mesh->nodes + wElm->nodes[0])->position,
				(mesh->nodes + wElm->nodes[1])->position,
				(mesh->nodes + wElm->nodes[2])->position) <
	 WLZ_MESH_TOLERANCE_SQ)
      {
	bId = wId;
      }
      ++wId;
    }
    /* Rebuild the mesh reusing the zombies and the first one (if new node
     * lies on and edge segment) or two of the working elements. */
    wElm = mesh->elements + mesh->nElem /* + 0 */ ;	      /* 1st element */
    wElm->strainU[0] = wElm->strainU[1] = wElm->strainU[2] = 0.0;
    wElm->strainA[0] = wElm->strainA[1] = wElm->strainA[2] = 0.0;
    wElm->neighbours[1] = mesh->nElem + 1;
    wElm->neighbours[2] = *(zElmIdVec + zElmCnt - 1);
    wElm = mesh->elements + mesh->nElem + 1;		      /* 2nd element */
    wElm->strainU[0] = wElm->strainU[1] = wElm->strainU[2] = 0.0;
    wElm->strainA[0] = wElm->strainA[1] = wElm->strainA[2] = 0.0;
    wElm->neighbours[1] = *(zElmIdVec /* + 0 */ );
    wElm->neighbours[2] = mesh->nElem /* + 0 */ ;
    wId = 2;
    zId = 0;
    while(zId < zElmCnt)			 /* The rest of the elements */
    {
      wElm = mesh->elements + mesh->nElem + wId;
      zElm = mesh->elements + *(zElmIdVec + zId);
      zElm->flags = wElm->flags;
      zElm->nodes[0] = newNodId;
      zElm->nodes[1] = wElm->nodes[1];
      zElm->nodes[2] = wElm->nodes[2];
      if(wElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0)
      {
	nId = wElm->neighbours[0];
	zElm->neighbours[0] = nId;
	zNElm = mesh->elements + nId;
	if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_0) != 0) &&
	   (zNElm->neighbours[0] == wElm->idx))
	{
	  zNElm->neighbours[0] = zElm->idx;
	}
	else if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_1) != 0) &&
	        (zNElm->neighbours[1] == wElm->idx))
	{
	  zNElm->neighbours[1] = zElm->idx;
	}
	else if(((zNElm->flags & WLZ_MESH_ELEM_FLAGS_NBR_2) != 0) &&
	        (zNElm->neighbours[2] == wElm->idx))
	{
	  zNElm->neighbours[2] = zElm->idx;
	}
      }
      if(zId < (zElmCnt - 1))
      {
	zElm->neighbours[1] = *(zElmIdVec + zId + 1);
      }
      else
      {
	zElm->neighbours[1] = mesh->nElem;
      }
      if(zId > 0)
      {
	zElm->neighbours[2] = *(zElmIdVec + zId - 1);
      }
      else
      {
	zElm->neighbours[2] = mesh->nElem + 1;
      }
      zElm->strainU[0] = zElm->strainU[1] = zElm->strainU[2] = 0.0;
      zElm->strainA[0] = zElm->strainA[1] = zElm->strainA[2] = 0.0;
      ++wId;
      ++zId;
    }
    if(bId < 0)
    {
      mElmCnt = mesh->nElem + 2;
    }
    else
    {
      mElmCnt = mesh->nElem + 1;
      /* Remove the bad element due to the new node being on the bad elements
       * edge segment from the mesh. */
      rId = (bId < 2)? (mesh->nElem + bId): *(zElmIdVec + bId - 2);
      rElm = mesh->elements + rId;
      nElm = mesh->elements + rElm->neighbours[1];
      nElm->flags = nElm->flags & ~(WLZ_MESH_ELEM_FLAGS_NBR_2);
      nElm = mesh->elements + rElm->neighbours[2];
      nElm->flags = nElm->flags & ~(WLZ_MESH_ELEM_FLAGS_NBR_1);
      if(bId != 1)
      {
	/* Recycle the removed element by using it in-place of the second
	 * working element. */
	wId = mesh->nElem + 1;
	WlzMeshElemReplace1With1(mesh, wId, rId);
      }
    }
    /* Finished: Update the mesh node and element counts. */
    mesh->nNodes = mNodCnt;
    mesh->nElem = mElmCnt;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Builds a vector of the id's of all the mesh elements
*		which include given the node (to be deleted), and a vector to
*		the id's of all the mesh nodes which lie around the perimeter
*		of the polygon formed by the vector of elements. The element
*		and node vectors are in ccw order around the given node. The
*		vector of nodes does not include the given node. All the
*		elements in the vector are marked as zombies ready for
*		recycling.
* \param	elmVec			Mesh element vector to build.
* \param	nodVec			Mesh node vector to build.
* \param	mesh			Given mesh transform.
* \param	dElmId			Index of an element known to use
*					the given node.
* \param	dElmNodId		Index of node to be deleted within
* 					the given element.
*/
static WlzErrorNum WlzMeshNodeDelVecBuild(WlzMeshIntVec *elmVec,
					  WlzMeshIntVec *nodVec,
					  WlzMeshTransform *mesh,
					  int dElmId,
					  int dElmNodId)
{
  int		eNodId0,
		eNodId1,
  		eNodId2,
		eNodId3,
		dNodId,
  		sNodId,
		nodId0,
		nodId1,
		nodId2,
		eId0,
		countOut;
  WlzMeshElem	*dElm,
  		*sElm,
		*elm0,
		*elm1,
		*elm2;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  dElm = mesh->elements + dElmId;
  dNodId = dElm->nodes[dElmNodId];
  /* First walk from the given node CW around it's enclosing polygon
   * until the walk returns either to the given node or to the first
   * after the given node. */
  elm0 = dElm;
  eNodId0 = (dElmNodId + 2) % 3;      /* First node known to be on polygon edge
  				       * and CW from the given node. */
  nodId0 = elm0->nodes[eNodId0];
  elm2 = elm0;
  eNodId2 = eNodId0;
  nodId2 = nodId0;
  countOut = mesh->nElem;
  do
  {
    elm1 = elm2;
    eNodId1 = eNodId2;
    nodId1 = nodId2;
    eNodId2 = (eNodId1 + 2) % 3;
    nodId2 = elm1->nodes[eNodId2];
    if(nodId2 != dNodId)
    {
      elm2 = elm1;
    }
    else
    {
      eNodId3 = (eNodId2 + 2) % 3;
      if((elm1->flags & nbrFlgTbl[eNodId3]) != 0)
      {
	eId0 = elm1->neighbours[eNodId3];
        elm2 = mesh->elements + eId0;
	if((eNodId3 = WlzMeshElemNodeIdxFromNodeIdx(elm2->nodes, nodId1)) < 0)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  eNodId2 = (eNodId3 + 2) % 3;
	  nodId2 = elm2->nodes[eNodId2];
	}
      }
    }
  } while((errNum == WLZ_ERR_NONE) &&
          (nodId2 != dNodId) && (nodId2 != nodId0) &&
	  (countOut-- > 0));
  if(errNum == WLZ_ERR_NONE)
  {
    if(countOut <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Add initial element and node to vectors. */
  if(errNum == WLZ_ERR_NONE)
  {
    sElm = elm1;
    sNodId = nodId1;
    nodVec->count = 0;
    errNum = WlzMeshAddToIntVec(nodVec, sNodId);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    elmVec->count = 0;
    sElm->flags |= WLZ_MESH_ELEM_FLAGS_ZOMBIE;
    errNum = WlzMeshAddToIntVec(elmVec, sElm->idx);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Now walk from the node found CCW around the given node until the walk
     * returns either to the given node of the first node of this walk. */
    elm2 = sElm;
    nodId2 = sNodId;
    if((eNodId2 = WlzMeshElemNodeIdxFromNodeIdx(elm2->nodes, nodId2)) < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      countOut = mesh->nElem;
      do
      {
        elm1 = elm2;
	eNodId1 = eNodId2;
	nodId1 = nodId2;
	eNodId2 = (eNodId1 + 1) % 3;
	nodId2 = elm1->nodes[eNodId2];
	if((elm1->flags & nbrFlgTbl[eNodId1]) != 0)
	{
	  eId0 = elm1->neighbours[eNodId1];
	  elm2 = mesh->elements + eId0;
	  if((eNodId2 = WlzMeshElemNodeIdxFromNodeIdx(elm2->nodes, nodId2)) < 0)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	}
	if((errNum == WLZ_ERR_NONE) && (elm1->idx != elm2->idx) &&
	   (elm2->idx != sElm->idx))
	{
          elm2->flags |= WLZ_MESH_ELEM_FLAGS_ZOMBIE;
	  errNum = WlzMeshAddToIntVec(elmVec, elm2->idx);
	}
	if((errNum == WLZ_ERR_NONE) && (nodId1 != nodId2) &&
	   (nodId2 != dNodId) &&
	   (nodId2 != *(nodVec->vector + 0)))
	{
	  errNum = WlzMeshAddToIntVec(nodVec, nodId2);
	}
      } while((errNum == WLZ_ERR_NONE) && (nodId2 != sNodId) &&
              (nodId2 != dNodId) &&
	      (countOut-- > 0));
      if(errNum == WLZ_ERR_NONE)
      {
	if(countOut <= 0)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
      }
    }
  }
/* #define WLZ_MESH_DEBUG */
#ifdef WLZ_MESH_DEBUG
  if(errNum == WLZ_ERR_NONE)
  {
    WlzMeshNode	*nod0;

    fprintf(stderr, "elements\n");
    for(eId0 = 0; eId0 < elmVec->count; ++eId0)
    {
      elm0 = mesh->elements + *(elmVec->vector + eId0);
      fprintf(stderr, "%d 0x%x (%d, %d, %d)\n",
              elm0->idx, elm0->flags,
	      elm0->nodes[0], elm0->nodes[1], elm0->nodes[2]);
    }
    fprintf(stderr, "nodes\n");
    for(nodId0 = 0; nodId0 < nodVec->count; ++nodId0)
    {
      nod0 = mesh->nodes + *(nodVec->vector + nodId0);
      fprintf(stderr, "%d (%g %g)\n",
      	      *(nodVec->vector + nodId0),
	      nod0->position.vtX, nod0->position.vtY);
    }
      fprintf(stderr, "\n");
  }
#endif /* WLZ_MESH_DEBUG */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Adds an into to a vector of int's expanding the vector as
*		required.
* \param	vec			Given int vector.
* \param	val			Value to add to vector.
*/
static WlzErrorNum WlzMeshAddToIntVec(WlzMeshIntVec *vec, int val)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(vec->count >= vec->size)
  {
    vec->size = (vec->size < vec->minSize)? vec->minSize:
    					    vec->size * vec->mulSize;
    if((vec->vector = (int *)AlcRealloc(vec->vector,
    					sizeof(int) * vec->size)) == NULL) 
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *(vec->vector + (vec->count)++) = val;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Realloc's a mesh ear list pool.
* \param	earList			Given ear list for recycling.
* \param	maxEars			Maximum number of ears that can be
*					used in the stack.
*/
static WlzErrorNum WlzMeshEarListRealloc(WlzMeshEarList *earList,
					 int maxEars)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  earList->inEars = 0;
  earList->topIn = NULL;
  if(earList->maxEars < maxEars)
  {
    if(earList->pool)
    {
      AlcFree(earList->pool);
    }
    if((earList->pool = (WlzMeshEar *)
			AlcMalloc(maxEars * sizeof(WlzMeshEar))) == NULL)
    {
      earList->maxEars = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      earList->maxEars = maxEars;
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Computes a power for the given ear, where the power is
*		given by:
*		\f[
		  \frac
		  {
		    \left|
		    \begin{array}{cccc}
		      x_o & x_1 & x_2 & x_3 \\
		      y_0 & y_1 & y_2 & y_3 \\
		      x_0^2 + y_0^2 & x_1^2 + y_1^2 & x_2^2 + y_2^2 &
		      x_p^2 + y_p^2 \\
		      1   & 1   & 1   & 1
		    \end{array}
		    \right|
		  }
		  {
		    \left|
		    \begin{array}{ccc}
		       x_0 & x_1 & x_2 \\
		       y_0 & y_1 & y_2 \\
		       1   & 1   & 1
		    \end{array}
		    \right|
		  }
		\f]
*		Where the given (CCW order) nodes are (x0,y0), (x1,y1)
*		and (x2,y2). The node to be deleted is (xp, yp). In practice
*		this function first checks that the nodesare not co-linear,
*		then that the nodes are CCW and if all ok so far then computes
*		the power.
* \param	mesh			Given mesh transform.
* \param	ear			Given mesh ear.
* \param	delNodId		Id of node to be deleted.
*/
static void	WlzMeshEarPowerSet(WlzMeshTransform *mesh, WlzMeshEar *ear,
				   int delNodId)
{
  double	x0y1,
  		x0y2,
		x0yp,
		x1y0,
		x1y2,
		x1yp,
		x2y0,
		x2y1,
		x2yp,
		xpy0,
		xpy1,
		xpy2;
  WlzDVertex2	nodVx0,
  		nodVx1,
		nodVx2,
		nodVxP;

  nodVx0 = (mesh->nodes + ear->nodes[0])->position;
  nodVx1 = (mesh->nodes + ear->nodes[1])->position;
  nodVx2 = (mesh->nodes + ear->nodes[2])->position;
  if(fabs(((nodVx0.vtY - nodVx2.vtY) *
  	   (nodVx1.vtX - nodVx0.vtX)) -
          ((nodVx0.vtX - nodVx2.vtX) *
	   (nodVx1.vtY - nodVx0.vtY))) < DBL_EPSILON)
  {
    ear->snArea2 = 0.0;
    ear->power = DBL_MAX;
  }
  else
  {
    x0y1 = nodVx0.vtX * nodVx1.vtY;
    x0y2 = nodVx0.vtX * nodVx2.vtY;
    x1y0 = nodVx1.vtX * nodVx0.vtY;
    x1y2 = nodVx1.vtX * nodVx2.vtY;
    x2y0 = nodVx2.vtX * nodVx0.vtY;
    x2y1 = nodVx2.vtX * nodVx1.vtY;
    if((ear->snArea2 = x0y1 - x0y2 - x1y0 + x1y2 + x2y0 - x2y1) < 
       WLZ_MESH_ELEM_AREA_TOLERANCE * 2)
    {
      ear->power = DBL_MAX;
    }
    else
    {
      nodVxP = (mesh->nodes + delNodId)->position;
      x0yp = nodVx0.vtX * nodVxP.vtY;
      x1yp = nodVx1.vtX * nodVxP.vtY;
      x2yp = nodVx2.vtX * nodVxP.vtY;
      xpy0 = nodVxP.vtX * nodVx0.vtY;
      xpy1 = nodVxP.vtX * nodVx1.vtY;
      xpy2 = nodVxP.vtX * nodVx2.vtY;
      ear->power = ((((nodVx0.vtX * nodVx0.vtX) +
			 (nodVx0.vtY * nodVx0.vtY))  *
			(x1y2 - x1yp - x2y1 + x2yp + xpy1 - xpy2)) -
		       (((nodVx1.vtX * nodVx1.vtX) +
			 (nodVx1.vtY * nodVx1.vtY))  *
			(x0y2 - x0yp - x2y0 + x2yp + xpy0 - xpy2)) +
		       (((nodVx2.vtX * nodVx2.vtX) +
			 (nodVx2.vtY * nodVx2.vtY))  *
			(x0y1 - x0yp - x1y0 + x1yp + xpy0 - xpy1)) -
		       (((nodVxP.vtX * nodVxP.vtX) +
			 (nodVxP.vtY * nodVxP.vtY))  * ear->snArea2)) /
		      ear->snArea2;
    }
  }
}

/*!
* \return	Ear with minimum power.
* \ingroup	WlzTransform
* \brief	Searches the input ear list for the ear with the minimum
*		power.
* \param	earList			Given ear list.
*/
static WlzMeshEar *WlzMeshEarGetMinPower(WlzMeshEarList *earList)
{
  int		eId0,
  		eId1,
		count;
  WlzMeshEar	*ear0,
		*ear1,
  		*minEar = NULL;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  if(earList && earList->topIn)
  {
    count = earList->inEars;
    minEar = earList->topIn;
    if(count > 3)
    {
      ear0 = minEar;
      while((--count > 0) && ear0)
      {
	ear0 = ear0->next;
	if(ear0 && (ear0->power < minEar->power))
	{
	  minEar = ear0;
	}
      }
    }
    else if(count == 3)
    {
      /* Special case when only three ears are left: Each of the three
       * remaining ears is a possibly incomplete rotation of the same
       * ear. */
      eId0 = 0;
      minEar = earList->topIn;
      ear0 = minEar->next;
      ear1 = ear0->next;
      while(minEar && (eId0 < 3))
      {
        if((minEar->flags & nbrFlgTbl[eId0]) == 0)
	{
	  /* Need to get this neighbour from one of the other two ears. */
	  eId1 = WlzMeshElemNodeIdxFromNodeIdx(ear0->nodes,
	  				       *(minEar->nodes + eId0));
	  if((ear0->flags & nbrFlgTbl[eId1]) != 0)
	  {
	    minEar->flags |= nbrFlgTbl[eId0];
	    minEar->neighbours[eId0] = ear0->neighbours[eId1];
	    minEar->selfNeighbour[eId0] = ear0->selfNeighbour[eId1];
	  }
	  else
	  {
	    eId1 = WlzMeshElemNodeIdxFromNodeIdx(ear1->nodes,
	    					 *(minEar->nodes + eId0));
	    if((ear1->flags & nbrFlgTbl[eId1]) != 0)
	    {
	      minEar->flags |= nbrFlgTbl[eId0];
	      minEar->neighbours[eId0] = ear1->neighbours[eId1];
	      minEar->selfNeighbour[eId0] = ear1->selfNeighbour[eId1];
	    }
	    else
	    {
	      minEar->power = DBL_MAX;
	      minEar->snArea2 = -1.0;
	    }
	  }
	}
        ++eId0;
      }
    }
    else if(count < 3)
    {
      minEar->power = DBL_MAX;
      minEar->snArea2 = -1.0;
    }
    if(earList->topIn == minEar)
    {
      earList->topIn = minEar->next;
    }
  }
  return(minEar);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Creates a list  of mesh ears from a vector of nodes and
*		the node to be deleted. The list is in the order of given
*		nodes, ie CCW around the node to be deleted.
* \param	earList			Given ear list for recycling.
* \param	elmVec			Vector of mesh elements which runs
*					CCW around the node to be deleted.
* \param	nodVec			Vector of mesh nodes which runs
*					CCW around the enclosing polygon of
*					the node to be deleted.
* \param	mesh			Given mesh transform.
* \param	delNodId		Index of an node to be deleted.
*/
static WlzErrorNum WlzMeshEarsCreate(WlzMeshEarList *earList,
				     WlzMeshIntVec *elmVec,
				     WlzMeshIntVec *nodVec,
				     WlzMeshTransform *mesh,
				     int delNodId)
{
  int		earId0,
  		elmVecId0 = 0;
  WlzMeshEar	*ear0,
  		*ear1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzMeshEarListRealloc(earList, nodVec->size);
  if(errNum == WLZ_ERR_NONE)
  {
    ear1 = NULL;
    for(earId0 = 0; earId0 < nodVec->count; ++earId0)
    {
      ear0 = earList->pool + (earList->inEars)++;
      /* Set ear nodes. */
      ear0->nodes[0] = *(nodVec->vector + 
      		        (earId0 + nodVec->count - 1) % nodVec->count);
      ear0->nodes[1] = *(nodVec->vector + earId0);
      ear0->nodes[2] = *(nodVec->vector + 
      		        (earId0 + nodVec->count + 1) % nodVec->count);
      /* Set ear power. */
      WlzMeshEarPowerSet(mesh, ear0, delNodId);
      /* Set ear flags and neighbours. */
      ear0->flags = WLZ_MESH_ELEM_FLAGS_NONE;
      WlzMeshEarMatchElm(mesh, ear0, elmVec, &elmVecId0);
      /* Add ear to the list. */
      if(ear1)
      {
        ear0->prev = ear1;
	ear1->next = ear0;
      }
      else
      {
        earList->topIn = ear0;
      }
      ear1 = ear0;
    }
    if(ear1)
    {
      /* Make the in ears list a circular list. */
      ear1->next = earList->topIn;
      earList->topIn->prev = ear1;
    }
  }
  return(errNum);
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Finds the indicies into the given element vector for the
*		elements which connect to the given ear and sets the ear's
*		neighbours and flags.
* \param	mesh			Given mesh transform.
* \param	ear			Ear (with nodes already set).
* \param	elmVec			Element vector.
* \param	elmVecIdP		Used to remember element vector index
*					between calls.
*/
static void	WlzMeshEarMatchElm(WlzMeshTransform *mesh,
				   WlzMeshEar *ear, WlzMeshIntVec *elmVec,
				   int *elmVecIdP)
{
  int		elmCnt,
  		elmVecId0,
		elmNbrId,
		earNbrId,
		earNodId0,
		earNodId1,
  		foundCnt;
  WlzMeshElem	*elm0,
  		*elm1;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  foundCnt = 0;
  elmCnt = elmVec->count;
  elmVecId0 = *elmVecIdP % elmVec->count;
  while((foundCnt < 2) && (elmCnt-- > 0))
  {
    elm0 = mesh->elements + *(elmVec->vector + elmVecId0);
    for(earNodId0 = 0; earNodId0 < 3; ++earNodId0)
    {
      earNodId1 = (earNodId0 + 1) % 3;
      elmNbrId = WlzMeshElemNbrIdxFromNodes(elm0, ear->nodes[earNodId0],
      				 	    ear->nodes[earNodId1]);
      if(elmNbrId >= 0)
      {
	++foundCnt;
        earNbrId = (earNodId0 + 2) % 3;
	if((elm0->flags & nbrFlgTbl[elmNbrId]) != 0)
	{
	  ear->flags |= nbrFlgTbl[earNbrId];
	  ear->neighbours[earNbrId] = elm0->neighbours[elmNbrId];
	  elm1 = mesh->elements + elm0->neighbours[elmNbrId];
	  ear->selfNeighbour[earNbrId] = WlzMeshElemNbrIdxFromNodes(elm1,
	  					ear->nodes[earNodId0],
						ear->nodes[earNodId1]);
	}
	else
	{
	  ear->flags &= ~(nbrFlgTbl[earNbrId]);
	}
      }
    }
    elmVecId0 = (elmVecId0 + 1) % elmVec->count;
  }
  if(elmVec->count > 1)
  {
    *elmVecIdP = (elmVecId0 + elmVec->count - 2) % elmVec->count;
  }
  else
  {
    *elmVecIdP = elmVecId0;
  }
}

/*!
* \return	<void>
* \ingroup	WlzTransform
* \brief	Unlinks the element with the given index from it's neighbours
*		in the mesh.
* \param	mesh			Given mesh transform.
* \param	eId			Index of element to be unlinked.
*/
static void	WlzMeshElemUnlink(WlzMeshTransform *mesh, int eId)
{
  int		eNId,
  		eNNId;
  WlzMeshElem	*elm,
  		*nElm;
  const unsigned int nbrFlgTbl[3] = {WLZ_MESH_ELEM_FLAGS_NBR_0,
				     WLZ_MESH_ELEM_FLAGS_NBR_1,
				     WLZ_MESH_ELEM_FLAGS_NBR_2};

  elm = mesh->elements + eId;
  for(eNId = 0; eNId < 3; ++eNId)	  /* Check each neighbour of element */
  {
    if((elm->flags & nbrFlgTbl[eNId]) != 0)
    {
      eNNId = 0;
      nElm = mesh->elements + elm->neighbours[eNId];
      for(eNNId = 0; eNNId < 3; ++eNNId)
      {
        if(((nElm->flags & nbrFlgTbl[eNNId]) != 0) &&
	   (nElm->neighbours[eNNId] = eId))
        {
	  nElm->flags = nElm->flags & ~(nbrFlgTbl[eNNId]);
	  eNNId = 3; 	      /* Can only be a neighbour once, so break loop */
	}
      }
    }
  }
}
