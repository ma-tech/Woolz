#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshFMar_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshFMar.c
* \author       Bill Hill
* \date         February 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Fast marching methods within conforming meshes.
* \ingroup	WlzMesh
* \todo         -
* \bug          None known.
*/
#include <limits.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

/* #define WLZ_CMESH_FMAR_DEBUG */
/* #define WLZ_CMESH_DEBUG_LOCATION */

/*!
* \struct	_WlzCMeshFMarQEnt
* \ingroup	WlzMesh
* \brief	An entry of a AlcHeap based queue.
* 		Typedef: ::WlzCMeshFMarQEnt.
*/
typedef struct _WlzCMeshFMarQEnt
{
  double		priority;	/*!< Entry priority highest priority
  					     at the head of the queue. */
  void			*entity;	/*!< Pointer to mesh entity. */
} WlzCMeshFMarQEnt;

/*!
* \struct	_WlzCMeshFMarElmQEnt
* \ingroup	WlzMesh
* \brief	An entry of an element queue.
* 		Typedef: ::WlzCMeshFMarElmQEnt.
*/
typedef struct _WlzCMeshFMarElmQEnt
{
  double		priority;	/*!< Priority of queue entry: The
  					     priority is the simple sum of
					     the priority of the nodes (2
					     for an upwind node or the current
					     node, 1 for any other known node
					     and zero for an unknown (down
					     wind) node). */
  WlzCMeshElmP		elm;		/*!< Element pointer. */
} WlzCMeshFMarElmQEnt;

static int			WlzCMeshFMarElmQCalcPriority2D(
				  WlzCMeshElm2D *elm,
				  WlzCMeshNod2D *cNod);
static int			WlzCMeshFMarElmQCalcPriority3D(
				  WlzCMeshElm3D *elm,
				  WlzCMeshNod3D *cNod);
static double			WlzCMeshFMarQSElmPriority2D(
				  WlzCMeshElm2D *elm,
				  double *dst,
				  WlzDVertex2 org);
static double			WlzCMeshFMarQSElmPriority3D(
				  WlzCMeshElm3D *elm,
				  double *dst,
				  WlzDVertex3 org);
static double 			WlzCMeshFMarSolve2D2(
				  WlzDVertex2 p0,
				  WlzDVertex2 p1,
				  WlzDVertex2 p2,
				  double d0,
				  double d1);
static int		 	WlzCMeshFMarCompute2D(
				  WlzCMeshNod2D *nod0,
				  WlzCMeshNod2D *nod1,
				  WlzCMeshNod2D *nod2,
				  double *distances,
				  WlzCMeshElm2D *elm);
static int			WlzCMeshFMarCompute3D(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances,
				  WlzCMeshElm3D *elm);
static int			WlzCMeshFMarCompute3D1(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static int			WlzCMeshFMarCompute3D2(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static int			WlzCMeshFMarCompute3D3(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static WlzErrorNum 		WlzCMeshFMarAddSeeds2D(
				  AlcHeap *queue,
				  WlzCMesh2D *mesh, 
				  int qMin,
				  double *distances,
				  int nSeeds,
				  WlzDVertex2 *seeds);
static WlzErrorNum 		WlzCMeshFMarAddSeeds3D(
				  AlcHeap *queue,
				  WlzCMesh3D *mesh,
				  int qMin,
				  double *distances,
				  int nSeeds,
				  WlzDVertex3 *seeds);
static WlzErrorNum 		WlzCMeshFMarAddSeed2D(
				  AlcHeap  *edgQ,
                                  WlzCMesh2D *mesh,
				  double *distances,
				  double *sDists,
				  WlzDVertex2 seed,
				  WlzUByte *eFlgs);
static WlzErrorNum 		WlzCMeshFMarAddSeed3D(
				  AlcHeap  *fceQ,
                                  WlzCMesh3D *mesh,
				  double *distances,
				  double *sDists,
				  WlzDVertex3 seed,
				  WlzUByte *eFlgs);
static WlzErrorNum 		WlzCMeshFMarQInsertNod2D(
				  AlcHeap *queue,
				  WlzCMeshNod2D *nod,
				  double dist);
static WlzErrorNum 		WlzCMeshFMarQInsertNod3D(
				  AlcHeap *queue,
				  WlzCMeshNod3D *nod,
				  double dist);
static WlzErrorNum 		WlzCMeshFMarSElmQInsert2D(
				  AlcHeap *edgQ,
				  WlzCMeshElm2D *elm,
				  double *dst,
				  WlzDVertex2 org);
static WlzErrorNum 		WlzCMeshFMarSElmQInsert3D(
				  AlcHeap *fceQ,
				  WlzCMeshElm3D *elm,
				  double *dst,
				  WlzDVertex3 org);
static WlzErrorNum 		WlzCMeshFMarElmQInit2D(
				  AlcHeap *queue,
				  WlzCMeshNod2D *nod);
static WlzErrorNum 		WlzCMeshFMarElmQInit3D(
				  AlcHeap *queue,
				  WlzCMeshNod3D *nod);
#ifdef WLZ_CMESH_DEBUG_LOCATION
extern	void  			WlzCMeshDebugResetNElmQuery(
				  void);
extern void			WlzCMeshDebugReportNElmQuery(
				  void);
#endif /* WLZ_CMESH_DEBUG_LOCATION */

/*!
* \return	A 2D domain object, an empty object if the mesh has
* 		no elements or NULL on error.
* \ingroup	WlzMesh
* \brief	Computes a new 2D domain object with values that are the
* 		distance from the given seeds within the given mesh.
* \param	mesh			Given mesh.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seeds			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
*					are not within the mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDistance2D(WlzCMesh2D *mesh,
				int nSeeds, WlzDVertex2 *seeds,
				WlzErrorNum *dstErr)
{
  int		idE,
  		idK,
		idN;
  double	d;
  double	*distances = NULL,
  		*dst;
  WlzObject	*obj0 = NULL,
  		*obj1 = NULL;
  WlzCMeshElm2D	*elm;
  WlzCMeshNod2D	*nod[3];
  WlzCMeshTransform *mTr = NULL;
  WlzObjectType	vTT;
  WlzValues	val;
  WlzPixelV	bgdV;
  WlzDVertex2	pos;
  WlzGreyWSpace gWsp;
  WlzIntervalWSpace iWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_DOUBLE;
  bgdV.v.dbv = DBL_MAX;
  vTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(mesh->res.elm.numEnt == 0)
  {
    obj1 = WlzMakeEmpty(&errNum);
  }
  else
  {
    if((distances = AlcMalloc(sizeof(double) * mesh->res.nod.maxEnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzCMeshFMarNodes2D(mesh, distances, nSeeds, seeds);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_2D_CMESH, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      mTr->mesh.m2 = mesh;
      obj0 = WlzCMeshToDomObj(mTr, 0, &errNum);
    }
    if(mTr)
    {
      mTr->mesh.m2 = NULL;
      (void )WlzFreeCMeshTransform(mTr);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      val.v = WlzNewValueTb(obj0, vTT, bgdV, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
			 obj0->domain, val, NULL, NULL, &errNum);
    }
    (void )WlzFreeObj(obj0);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzInitGreyScan(obj1, &iWsp, &gWsp);
      while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
      {
        dst = gWsp.u_grintptr.dbp;
	pos.vtY = iWsp.linpos;
	idE = -1;
	for(idK = iWsp.lftpos; idK <= iWsp.rgtpos; ++idK)
	{
	  pos.vtX = idK;
	  if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE,
	                                      pos.vtX, pos.vtY,
					      0, &idN)) >= 0)
	  {
	    elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	    nod[0] = elm->edu[0].nod;
	    nod[1] = elm->edu[1].nod;
	    nod[2] = elm->edu[2].nod;
	    d = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos, nod[2]->pos,
	        		        distances[nod[0]->idx],
					distances[nod[1]->idx],
					distances[nod[2]->idx],
					pos);
	  }
	  else if((idN >= 0) && (idN < mesh->res.nod.maxEnt))
	  {
	    d = distances[idN];
	  }
	  else
	  {
	    d = DBL_MAX;
	  }
	  *dst++ = d;
	}
      }
      if(errNum == WLZ_ERR_EOO)
      {
        errNum = WLZ_ERR_NONE;
      }
    }
  }
  AlcFree(distances);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj1);
    obj1 = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(obj1);
}


/*!
* \return	A 3D domain object, an empty object if the mesh has
* 		no elements or NULL on error.
* \ingroup	WlzMesh
* \brief	Computes a new 3D domain object with values that are the
* 		distance from the given seeds within the given mesh.
* \param	mesh			Given mesh.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seeds			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
*					are not within the mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDistance3D(WlzCMesh3D *mesh,
				int nSeeds, WlzDVertex3 *seeds,
				WlzErrorNum *dstErr)
{
  int		idE,
  		idK,
		idP,
		idN,
		pCnt;
  double	d;
  double	*distances = NULL,
  		*dst;
  WlzObject	*obj0 = NULL,
  		*obj1 = NULL,
		*obj2 = NULL;
  WlzCMeshElm3D	*elm;
  WlzCMeshNod3D	*nod[3];
  WlzCMeshTransform *mTr = NULL;
  WlzObjectType	vTT;
  WlzValues	val;
  WlzPixelV	bgdV;
  WlzDVertex3	pos;
  WlzGreyWSpace gWsp;
  WlzIntervalWSpace iWsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bgdV.type = WLZ_GREY_DOUBLE;
  bgdV.v.dbv = DBL_MAX;
  vTT = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(mesh->res.elm.numEnt == 0)
  {
    obj1 = WlzMakeEmpty(&errNum);
  }
  else
  {
    if((distances = AlcMalloc(sizeof(double) * mesh->res.nod.maxEnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzCMeshFMarNodes3D(mesh, distances, nSeeds, seeds);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      mTr = WlzMakeCMeshTransform(WLZ_TRANSFORM_3D_CMESH, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      mTr->mesh.m3 = mesh;
      obj0 = WlzCMeshToDomObj(mTr, 0, &errNum);
    }
    if(mTr)
    {
      mTr->mesh.m3 = NULL;
      (void )WlzFreeCMeshTransform(mTr);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      val.vox = WlzNewValuesVox(obj0, vTT, bgdV, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ,
			 obj0->domain, val, NULL, NULL, &errNum);
    }
    (void )WlzFreeObj(obj0);
    if(errNum == WLZ_ERR_NONE)
    {
      pCnt = obj1->domain.p->lastpl - obj1->domain.p->plane1 + 1;
      for(idP = 0; idP < pCnt; ++idP)
      {
        obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ,
	                   *(obj1->domain.p->domains + idP),
			   *(obj1->values.vox->values + idP),
			   NULL, NULL, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  pos.vtZ = obj1->domain.p->plane1 + idP;
          errNum = WlzInitGreyScan(obj2, &iWsp, &gWsp);
	}
	while((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE)
	{
	  dst = gWsp.u_grintptr.dbp;
	  pos.vtY = iWsp.linpos;
	  idE = -1;
	  for(idK = iWsp.lftpos; idK <= iWsp.rgtpos; ++idK)
	  {
	    pos.vtX = idK;
	    if((idE = WlzCMeshElmEnclosingPos3D(mesh, idE,
						pos.vtX, pos.vtY, pos.vtZ,
						0, &idN)) >= 0)
	    {
	      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	      WlzCMeshElmGetNodes3D(elm, nod + 0, nod + 1, nod + 2, nod + 3);
	      d = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
	                                  nod[2]->pos, nod[3]->pos,
					  distances[nod[0]->idx],
					  distances[nod[1]->idx],
					  distances[nod[2]->idx],
					  distances[nod[3]->idx],
					  pos);
	    }
	    else if((idN >= 0) && (idN < mesh->res.nod.maxEnt))
	    {
	      d = distances[idN];
	    }
	    else
	    {
	      d = DBL_MAX;
	    }
	    *dst++ = d;
	  }
	}
        (void )WlzFreeObj(obj2);
	if(errNum == WLZ_ERR_EOO)
	{
	  errNum = WLZ_ERR_NONE;
	}
	else if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
  }
  AlcFree(distances);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(obj1);
    obj1 = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(obj1);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes constrained distances within a mesh by
* 		propagating wavefronts within a 2D conforming mesh. The
* 		wavefronts are propagated from either the mesh boundary
* 		or a number of seed positions within the mesh.
* 		The given mesh will have modified node and element flags
* 		on return. 
* \param	mesh			Given mesh.
* \param	distances		Array for computed distances.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seeds			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
* 					any seeds are not within the
* 					mesh.
*/
WlzErrorNum	WlzCMeshFMarNodes2D(WlzCMesh2D *mesh, double *distances,
				int nSeeds, WlzDVertex2 *seeds)
{
  int		idM,
  		idN,
  		idS,
		cnt;
  WlzCMeshNod2D	*nod0,
                *nod1;
  WlzCMeshNod2D	*nodes[3];
  WlzCMeshElm2D	*elm;
  AlcHeap *nodQ = NULL;
  AlcHeap *elmQ = NULL;
  WlzCMeshFMarQEnt *nodQEntP;
  WlzCMeshFMarElmQEnt *elmQEntP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((distances == NULL) || ((nSeeds > 0) && (seeds == NULL)))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set distance for all mesh nodes to maximum value. */
    WlzValueSetDouble(distances, DBL_MAX, mesh->res.nod.maxEnt);
    /* Clear mesh node flags, set boundary node flags and count number of
     * boundary nodes. */
    WlzCMeshClearNodFlags2D(mesh, WLZ_CMESH_NOD_FLAG_ALL);
    if((cnt = WlzCMeshSetBoundNodFlags2D(mesh)) <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Create and then initialise the active node queue using the given seed
   * or boundary nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodQ = AlcHeapNew(sizeof(WlzCMeshFMarQEnt), cnt, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nSeeds > 0)
    {
      errNum = WlzCMeshFMarAddSeeds2D(nodQ, mesh, cnt + 1,
                                      distances, nSeeds, seeds);
    }
    else
    {
      nSeeds = cnt;
      if((seeds = (WlzDVertex2 *)
                  AlcMalloc(nSeeds * sizeof(WlzDVertex2))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	idS = 0;
        for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
	{
	  nod0 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if((nod0->idx >= 0) &&
	     ((nod0->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0))
	  {
	    seeds[idS] = nod0->pos;
	    ++idS;
	  }
	}
	errNum = WlzCMeshFMarAddSeeds2D(nodQ, mesh, cnt + 1,
				        distances, nSeeds, seeds);
	AlcFree(seeds);
      }
    }
  }
  /* Create element queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((elmQ = AlcHeapNew(sizeof(WlzCMeshFMarElmQEnt), 1024, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Until the queue is empty: Pop the node with lowest priority (ie
   * lowest distance) from the queue, and process it. */
  while((errNum == WLZ_ERR_NONE) &&
	((nodQEntP = (WlzCMeshFMarQEnt *)AlcHeapTop(nodQ)) != NULL))
  {
    /* Find all neighbouring nodes that are neither active nor upwind.
     * For each of these neighbouring nodes, compute their distance, set
     * them to active and insert them into the queue.*/
    nod0 = (WlzCMeshNod2D *)(nodQEntP->entity);
    AlcHeapEntFree(nodQ);
    if((nod0->flags & WLZ_CMESH_NOD_FLAG_UPWIND) == 0)
    {
      errNum = WlzCMeshFMarElmQInit2D(elmQ, nod0);
      if(errNum == WLZ_ERR_NONE)
      {
	/* While element list is not empty, remove element and compute all
	 * node distances for it. */
	while((elmQEntP = (WlzCMeshFMarElmQEnt *)AlcHeapTop(elmQ)) != NULL)
	{
	  elm = elmQEntP->elm.e2;
	  AlcHeapEntFree(elmQ);
	  /* Compute distances. */
	  for(idN = 0; idN < 3; ++idN)
	  {
	    if(elm->edu[idN].nod == nod0)
	    {
	      break;
	    }
	  }
	  cnt = 0;
	  for(idM = 0; idM < 3; ++idM)
	  {
	    nodes[idM] = elm->edu[(idN + idM) % 3].nod;
	  }
	  if((nodes[2]->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0)
	  {
	    nod1 = nodes[2];
	    nodes[2] = nodes[1];
	    nodes[1] = nod1;
	  }
	  if(WlzCMeshFMarCompute2D(nodes[0], nodes[1], nodes[2],
				   distances, elm) != 0)
	  {
	    errNum = WlzCMeshFMarQInsertNod2D(nodQ, nodes[2],
					      distances[nodes[2]->idx]);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
	/* Set the current node to be upwind. */
	if(errNum == WLZ_ERR_NONE)
	{
	  nod0->flags = (nod0->flags & ~(WLZ_CMESH_NOD_FLAG_ACTIVE)) |
					 WLZ_CMESH_NOD_FLAG_KNOWN |
		                         WLZ_CMESH_NOD_FLAG_UPWIND;
	}
      }
      AlcHeapAllEntFree(elmQ, 0);
    }
  }
  /* Clear up. */
  AlcHeapFree(elmQ);
  AlcHeapFree(nodQ);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes constrained distances within a mesh by
* 		propagating wavefronts within a 3D conforming mesh. The
* 		wavefronts are propagated from either the mesh boundary
* 		or a number of seed positions within the mesh.
* 		The given mesh will have modified node and element flags
* 		on return. 
* \param	mesh			Given mesh.
* \param	distances		Array for computed distances.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seeds			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
* 					any seeds are not within the
* 					mesh.
*/
WlzErrorNum	WlzCMeshFMarNodes3D(WlzCMesh3D *mesh, double *distances,
				int nSeeds, WlzDVertex3 *seeds)
{
  int		idM,
  		idN,
		idP,
  		idS,
		cnt;
  WlzCMeshNod3D	*nod0,
                *nod1;
  WlzCMeshNod3D	*nodes[4];
  WlzCMeshElm3D	*elm;
  AlcHeap *nodQ = NULL;
  AlcHeap *elmQ = NULL;
  WlzCMeshFMarQEnt *nodQEntP;
  WlzCMeshFMarElmQEnt *elmQEntP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((distances == NULL) || ((nSeeds > 0) && (seeds == NULL)))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set distance for all mesh nodes to maximum value. */
    WlzValueSetDouble(distances, DBL_MAX, mesh->res.nod.maxEnt);
    /* Clear mesh node flags, set boundary node flags and count number of
     * boundary nodes. */
    WlzCMeshClearNodFlags3D(mesh, WLZ_CMESH_NOD_FLAG_ALL);
    if((cnt = WlzCMeshSetBoundNodFlags3D(mesh)) <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Create and then initialise the active node queue using the given seed
   * or boundary nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodQ = AlcHeapNew(sizeof(WlzCMeshFMarQEnt), cnt, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_CMESH_DEBUG_LOCATION
    WlzCMeshDebugResetNElmQuery();
#endif /* WLZ_CMESH_DEBUG_LOCATION */
    if(nSeeds > 0)
    {
      errNum = WlzCMeshFMarAddSeeds3D(nodQ, mesh, cnt + 1,
                                      distances, nSeeds, seeds);
    }
    else
    {
      nSeeds = cnt;
      if((seeds = (WlzDVertex3 *)
                  AlcMalloc(nSeeds * sizeof(WlzDVertex3))) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	idS = 0;
        for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
	{
	  nod0 = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if((nod0->idx >= 0) &&
	     ((nod0->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0))
	  {
	    seeds[idS] = nod0->pos;
	    ++idS;
	  }
	}
	errNum = WlzCMeshFMarAddSeeds3D(nodQ, mesh, cnt + 1,
				        distances, nSeeds, seeds);
	AlcFree(seeds);
      }
    }
  }
  /* Create element queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((elmQ = AlcHeapNew(sizeof(WlzCMeshFMarElmQEnt), 1024, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Until the queue is empty: Pop the node with lowest priority (ie
   * lowest distance) from the queue, and process it. */
  while((errNum == WLZ_ERR_NONE) &&
	((nodQEntP = (WlzCMeshFMarQEnt *)AlcHeapTop(nodQ)) != NULL))
  {
    /* Find all neighbouring nodes that are neither active nor upwind.
     * For each of these neighbouring nodes, compute their distance, set
     * them to active and insert them into the queue.*/
    nod0 = (WlzCMeshNod3D *)(nodQEntP->entity);
    AlcHeapEntFree(nodQ);
    if((nod0->flags & WLZ_CMESH_NOD_FLAG_UPWIND) == 0)
    {
      errNum = WlzCMeshFMarElmQInit3D(elmQ, nod0);
      if(errNum == WLZ_ERR_NONE)
      {
	/* While element list is not empty, remove element and compute all
	 * node distances for it. */
	while((elmQEntP = (WlzCMeshFMarElmQEnt *)AlcHeapTop(elmQ)) != NULL)
	{
	  elm = elmQEntP->elm.e3;
	  AlcHeapEntFree(elmQ);
	  /* Compute distances: Find the elements nodes, sort them by
	   * distance and then compute the unknown node distances. */
	  WlzCMeshElmGetNodes3D(elm, nodes + 0, nodes + 1,
	                             nodes + 2, nodes + 3);
	  for(idM = 0; idM < 3; ++idM)
	  {
	    idP = idM;
	    for(idN = idM + 1; idN < 4; ++idN)
	    {
	      if(distances[nodes[idN]->idx] < distances[nodes[idP]->idx])
	      {
	        idP = idN;
	      }
	    }
	    nod1 = nodes[idM]; nodes[idM] = nodes[idP]; nodes[idP] = nod1;
	  }
	  if(WlzCMeshFMarCompute3D(nodes[0], nodes[1], nodes[2], nodes[3],
				   distances, elm) != 0)
	  {
	    errNum = WlzCMeshFMarQInsertNod3D(nodQ, nodes[3],
					      distances[nodes[3]->idx]);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
	/* Set the current node to be upwind. */
	if(errNum == WLZ_ERR_NONE)
	{
	  nod0->flags = (nod0->flags & ~(WLZ_CMESH_NOD_FLAG_ACTIVE)) |
					 WLZ_CMESH_NOD_FLAG_KNOWN |
		                         WLZ_CMESH_NOD_FLAG_UPWIND;
	}
      }
      AlcHeapAllEntFree(elmQ, 0);
    }
  }
#ifdef WLZ_CMESH_DEBUG_LOCATION
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshDebugReportNElmQuery();
  }
#endif /* WLZ_CMESH_DEBUG_LOCATION */
  /* Clear up. */
  AlcHeapFree(elmQ);
  AlcHeapFree(nodQ);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts the given node from a 2D constraied mesh into
* 		the given priority queue. The given distances must be valid
* 		for the given node and any upwind nodes.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	nod			Given node to insert into the queue.
* \param	dist			Node distance.
*/
static WlzErrorNum WlzCMeshFMarQInsertNod2D(AlcHeap *queue,
					    WlzCMeshNod2D *nod, double dist)
{
  WlzCMeshFMarQEnt ent;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  nod->flags |= WLZ_CMESH_NOD_FLAG_KNOWN | WLZ_CMESH_NOD_FLAG_ACTIVE;
  ent.entity = nod;
  ent.priority = dist;
  if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts the given node from a 3D constraied mesh into
* 		the given priority queue. The given distances must be valid
* 		for the given node and any upwind nodes.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	nod			Given node to insert into the queue.
* \param	dist			Node distance.
*/
static WlzErrorNum WlzCMeshFMarQInsertNod3D(AlcHeap *queue,
				WlzCMeshNod3D *nod, double dist)
{
  WlzCMeshFMarQEnt ent;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  nod->flags |= WLZ_CMESH_NOD_FLAG_KNOWN | WLZ_CMESH_NOD_FLAG_ACTIVE;
  ent.entity = nod;
  ent.priority = dist;
  if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  return(errNum);
}

/*!
* \return	Edge priority.
* \ingroup	WlzMesh
* \brief	Priority for 2D seed element in queue. The priority increases
* 		with the distance of the unknown node from the seed.
* \param	elm			Element use to compute priority for.
* \param	dst			Known node distances from the seed.
* \param	org			Origin from which to compute minimum
* 					distance of node on edge.
*/
static double	WlzCMeshFMarQSElmPriority2D(WlzCMeshElm2D *elm,
					   double *dst,
					   WlzDVertex2 org)
{
  int		idE;
  double	d = 0.0;

  for(idE = 0; idE < 3; ++idE)
  {
    if(dst[elm->edu[idE].nod->idx] > DBL_MAX / 2.0)
    {
      break;
    }
  }
  if(idE < 3)
  {
    d = WlzGeomDist2D(org, elm->edu[idE].nod->pos);
  }
  return(d);
}

/*!
* \return	Edge priority.
* \ingroup	WlzMesh
* \brief	Priority for 3D seed element in queue. The priority increases
* 		with the distance of the unknown node from the seed.
* \param	elm			Element use to compute priority for.
* \param	dst			Known node distances from the seed.
* \param	org			Origin from which to compute minimum
* 					distance of node on edge.
*/
static double	WlzCMeshFMarQSElmPriority3D(WlzCMeshElm3D *elm,
					   double *dst,
					   WlzDVertex3 org)
{
  int		idN;
  double	d = 0.0;
  WlzCMeshNod3D	*nodes[4];

  WlzCMeshElmGetNodes3D(elm, nodes + 0, nodes + 1, nodes + 2, nodes + 3 );
  for(idN = 0; idN < 4; ++idN)
  {
    if(dst[nodes[idN]->idx] > DBL_MAX / 2.0)
    {
      break;
    }
  }
  if(idN < 4)
  {
    d = WlzGeomDist3D(org, nodes[idN]->pos);
  }
  return(d);
}

/*!
* \return	Non zero if distance computed and less than current distance.
* \ingroup	WlzMesh
* \brief	Computes wavefront distance for the given unknown node.
* 		Given a pair of nodes that are connected by a single edge
* 		in a 2D conforming mesh, the first of which has an unknown
* 		and the second a known wavefront propagation distance, this
* 		function computes the unknown distance.
* 		A triangle specified by it's nodes (nod0, nod1 and nod2)
* 		is given and the distance of the unknown node (nod2)
* 		is computed. Internal angles at the nodes of phi0, phi1 and
* 		phi2. Edge lengths opposite to the similarly numbered nodes
* 		of len0, len1 and len2. The solution is similar to that in
* 		"Fast Sweeping Methods For Eikonal equations On triangular
* 		meshes", Jianliang Qian, etal, SIAM journal on Mumerical
* 		Analysis, Vol 45, pp 83-107, 2007. But uses simple
* 		interpolation in the case of phi2 being obtuse.
* \param	nod0			A known node directly connected to
* 					the unknown node.
* \param	nod1			A second node which shares an element
* 					with nod2 and nod0. The distance for
* 					nod1 is >= the distance for  nod0.
* \param	nod2			Unknown node.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
* 		elm			Element using these nodes.
*/
static int	WlzCMeshFMarCompute2D(WlzCMeshNod2D *nod0,
				WlzCMeshNod2D *nod1, WlzCMeshNod2D *nod2, 
				double *distances, WlzCMeshElm2D *elm)
/* #define NEW_CODE_HACK */
#ifdef NEW_CODE_HACK
{
  int		idN,
		rtn = 0;
  double	d0,
		d2;
  WlzCMeshNod2D *nod3;
  WlzCMeshEdgU2D *edu;
  const double  maxCosAng = 0.996195;		        /* Cosine 5 degrees. */

  /* For each element which is a neighbour of the given one find a
   * node in the element which is not a node of the current element.
   * If this node has a smaller distance than nod1 use it as nod1. */
  for(idN = 0; idN < 3; ++idN)
  {
    edu = elm->edu + idN;
    if((edu->opp != NULL) && (edu->opp != edu))
    {
      nod3 = edu->opp->next->next->nod;
      if(distances[nod3->idx] < distances[nod1->idx])
      {
	d0 = WlzGeomCos3V(nod0->pos, nod3->pos, nod2->pos);
	if(d0 < maxCosAng)
	{
	  nod1 = nod3;
	}
      }
    }
  }
  d2 = WlzCMeshFMarSolve2D2(nod0->pos, nod1->pos, nod2->pos,
			    *(distances + nod0->idx),
			    *(distances + nod1->idx));
  if(d2 < *(distances + nod2->idx))
  {
    *(distances + nod2->idx) = d2;
    rtn = 1;
  }
  return(rtn);
}
#else /* NEW_CODE_HACK */
{
  int		idN,
  		flg,
		rtn = 0;
  double	d0,
		d1,
		theta;
  double	len[3],
		lenSq[3],
  		phi[3],
		dist[3];
  WlzDVertex2	del,
  		pos;
  WlzCMeshNod2D *nod3;
  WlzCMeshEdgU2D *edu;
  const double  maxCosAng = 0.996195;		        /* Cosine 5 degrees. */

  /* For each element which is a neighbour of the given one find a
   * node in the element which is not a node of the current element.
   * If this node has a smaller distance than nod1 use it as nod1. */
  for(idN = 0; idN < 3; ++idN)
  {
    edu = elm->edu + idN;
    if((edu->opp != NULL) && (edu->opp != edu))
    {
      nod3 = edu->opp->next->next->nod;
      if(distances[nod3->idx] < distances[nod1->idx])
      {
	d0 = WlzGeomCos3V(nod0->pos, nod3->pos, nod2->pos);
	if(d0 < maxCosAng)
	{
	  nod1 = nod3;
	}
      }
    }
  }
  dist[0] = distances[nod0->idx];
  dist[1] = distances[nod1->idx];
  if(dist[1] < dist[0])
  {
    nod3 = nod0;
    nod0 = nod1;
    nod1 = nod3;
    dist[0] = dist[1];
    dist[1] = distances[nod1->idx];
  }
  if((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) == 0)
  {
    /* Only have one node with known distance in this element. Compute the
     * distance using just the edge length. */
    WLZ_VTX_2_SUB(del, nod2->pos, nod0->pos);
    lenSq[1] = WLZ_VTX_2_SQRLEN(del);
    len[1] = sqrt(lenSq[1]);
    dist[2] = *(distances + nod0->idx) + len[1];
  }
  else
  {
    /* Have two known nodes (nod0 and nod1) and one unknown node (nod2). */
    dist[0] = *(distances + nod0->idx);
    dist[1] = *(distances + nod1->idx);
    WLZ_VTX_2_SUB(del, nod0->pos, nod1->pos);
    lenSq[2] = WLZ_VTX_2_SQRLEN(del);
    len[2] = sqrt(lenSq[2]);
    if(len[2] < DBL_EPSILON)
    {
      /* Element is degenerate so just use edge length to compute distance. */
      dist[2] = *(distances + nod0->idx) + len[1];
    }
    else
    {
      /* Nod0 and nod1 both have known distances and len2 is greater than
       * zero. */
      WLZ_VTX_2_SUB(del, nod1->pos, nod2->pos);
      lenSq[0] = WLZ_VTX_2_SQRLEN(del);
      len[0] = sqrt(lenSq[0]);
      WLZ_VTX_2_SUB(del, nod2->pos, nod0->pos);
      lenSq[1] = WLZ_VTX_2_SQRLEN(del);
      len[1] = sqrt(lenSq[1]);
      phi[2] = acos((lenSq[0] + lenSq[1] - lenSq[2]) /
                    (2.0 * len[0] * len[1]));
      if(ALG_M_PI_2 < phi[2])
      {
        /* Angle at unknown node is obtuse so interpolate a virtual node
	 * half way between nod0 and nod1 which will bisect phi2  with
	 * dist1 = (dist0 + dist1)/2, pos1 = (pos1 + pos0)/2. */
        WLZ_VTX_2_ADD(pos, nod0->pos, nod1->pos);
	WLZ_VTX_2_SCALE(pos, pos, 0.5);
	dist[1] = 0.5 * (dist[1] + dist[0]);
	WLZ_VTX_2_SUB(del, nod0->pos, pos);
	lenSq[2] = WLZ_VTX_2_SQRLEN(del);
	len[2] = sqrt(lenSq[2]);
	WLZ_VTX_2_SUB(del, pos, nod2->pos);
	lenSq[0] = WLZ_VTX_2_SQRLEN(del);
	len[0] = sqrt(lenSq[0]);
      }
      flg = 0;
      if(len[2] > (dist[1] - dist[0]))
      {
        theta = asin((dist[1] - dist[0]) / len[2]);
	phi[0] = acos((lenSq[1] + lenSq[2] - lenSq[0]) /
	              (2.0 * len[1] * len[2]));
	phi[1] = acos((lenSq[2] + lenSq[0] - lenSq[1]) /
	              (2.0 * len[2] * len[0]));
	if((theta > DBL_EPSILON) &&
	   (theta > (phi[1] - ALG_M_PI_2)) && 
	   ((ALG_M_PI_2 - phi[0]) > theta))
	{
	  flg = 1;
	  d0 = len[0] * sin(phi[1] - theta); /* h0 */
	  d1 = len[1] * sin(phi[0] + theta); /* h1 */
	  dist[2] = 0.5 * ((d0 + dist[0]) + (d1 + dist[1]));
	}
      }
      if(flg == 0)
      {
	d0 = dist[0] + len[1];
	d1 = dist[1] + len[0];
	dist[2] = ALG_MIN(d0, d1);
      }
    }
  }
  if(dist[2] < *(distances + nod2->idx))
  {
    *(distances + nod2->idx) = dist[2];
    rtn = 1;
  }
  return(rtn);
}
#endif /* NEW_CODE_HACK */

/*!
* \return	Computed distance of the third vertex.
* \ingroup	WlzMesh
* \brief	Computes the distances of the third vertex in a triangle
* 		of three vertices, with the first two having known distances
* 		which are given.
* \param	p0			Position of first vertex.
* \param	p1			Position of second vertex.
* \param	p2			Position of third vertex.
* \param	d0			Known distance of first vertex.
* \param	d1			Known distance of second vertex.
*/
static double 	WlzCMeshFMarSolve2D2(WlzDVertex2 p0, WlzDVertex2 p1,
			             WlzDVertex2 p2, double d0, double d1)
{
  int		flg;
  double	d2,
  		d20,
		d21,
		theta;
  double	len[3],
		lenSq[3],
  		phi[3];
  WlzDVertex2	del,
  		pos;

  if(d1 < d0)
  {
    pos = p0; p0 = p1; p1 = pos;
    d20 = d0; d0 = d1; d1 = d20;
  }
  WLZ_VTX_2_SUB(del, p0, p1);
  lenSq[2] = WLZ_VTX_2_SQRLEN(del);
  len[2] = sqrt(lenSq[2]);
  if(len[2] < DBL_EPSILON)
  {
    /* Element is degenerate so just use edge length to compute distance. */
    d2 = d0 + len[1];
  }
  else
  {
    WLZ_VTX_2_SUB(del, p1, p2);
    lenSq[0] = WLZ_VTX_2_SQRLEN(del);
    len[0] = sqrt(lenSq[0]);
    WLZ_VTX_2_SUB(del, p2, p0);
    lenSq[1] = WLZ_VTX_2_SQRLEN(del);
    len[1] = sqrt(lenSq[1]);
    phi[2] = acos((lenSq[0] + lenSq[1] - lenSq[2]) / (2.0 * len[0] * len[1]));
    if(ALG_M_PI_2 < phi[2])
    {
      /* Angle at unknown node is obtuse so interpolate a virtual node
       * half way between nod0 and nod1 which will bisect phi2  with
       * d1 = (d0 + d1)/2, pos1 = (pos1 + pos0)/2. */
      WLZ_VTX_2_ADD(pos, p0, p1);
      WLZ_VTX_2_SCALE(pos, pos, 0.5);
      d1 = 0.5 * (d1 + d0);
      WLZ_VTX_2_SUB(del, p0, pos);
      lenSq[2] = WLZ_VTX_2_SQRLEN(del);
      len[2] = sqrt(lenSq[2]);
      WLZ_VTX_2_SUB(del, pos, p2);
      lenSq[0] = WLZ_VTX_2_SQRLEN(del);
      len[0] = sqrt(lenSq[0]);
    }
    flg = 0;
    if(len[2] > (d1 - d0))
    {
      theta = asin((d1 - d0) / len[2]);
      phi[0] = acos((lenSq[1] + lenSq[2] - lenSq[0])/(2.0 * len[1] * len[2]));
      phi[1] = acos((lenSq[2] + lenSq[0] - lenSq[1])/(2.0 * len[2] * len[0]));
      if((theta > DBL_EPSILON) && (theta > (phi[1] - ALG_M_PI_2)) && 
	 ((ALG_M_PI_2 - phi[0]) > theta))
      {
	flg = 1;
	d20 = len[0] * sin(phi[1] - theta); /* h0 */
	d21 = len[1] * sin(phi[0] + theta); /* h1 */
	d2 = 0.5 * ((d20 + d0) + (d21 + d1));
      }
    }
    if(flg == 0)
    {
      d20 = d0 + len[1];
      d21 = d1 + len[0];
      d2 = ALG_MIN(d20, d21);
    }
  }
  return(d2);
}

/*!
* \ingroup	WlzMesh
* \return	Non zero if distance computed and less than current distance.
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element.
* 		This function just classifies the problem according to the
* 		number of known nodes and passes then calls the appropriate
* 		function.
* \param	nod0			Current node (first).
* \param	nod1			Second node of the element.
* \param	nod2			Third node of the element.
* \param	nod3			Fourth node of the element.
* \param	distances		Array of distanes indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
* 		elm			Element using these nodes.
*/
static int	WlzCMeshFMarCompute3D(WlzCMeshNod3D *nod0,
				      WlzCMeshNod3D *nod1,
				      WlzCMeshNod3D *nod2,
				      WlzCMeshNod3D *nod3,
				      double *distances,
				      WlzCMeshElm3D *elm)
{
  int		rtn = 0,
  		kwnMsk = 0;

  kwnMsk = (((nod0->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 0) |
           (((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 1) |
           (((nod2->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 2) |
           (((nod3->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 3);
  switch(kwnMsk)
  {
    case  1: /* 0001 */
      rtn = WlzCMeshFMarCompute3D1(nod0, nod1, nod2, nod3, distances);
      break;
    case  2: /* 0010 */
      rtn = WlzCMeshFMarCompute3D1(nod1, nod2, nod3, nod0, distances);
      break;
    case  3: /* 0011 */
      rtn = WlzCMeshFMarCompute3D2(nod0, nod1, nod2, nod3, distances);
      break;
    case  4: /* 0100 */
      rtn = WlzCMeshFMarCompute3D1(nod2, nod3, nod0, nod1, distances);
      break;
    case  5: /* 0101 */
      rtn = WlzCMeshFMarCompute3D2(nod0, nod2, nod1, nod3, distances);
      break;
    case  6: /* 0110 */
      rtn = WlzCMeshFMarCompute3D2(nod1, nod2, nod0, nod3, distances);
      break;
    case  7: /* 0111 */
      rtn = WlzCMeshFMarCompute3D3(nod0, nod1, nod2, nod3, distances);
      break;
    case  8: /* 1000 */
      rtn = WlzCMeshFMarCompute3D1(nod3, nod0, nod1, nod2, distances);
      break;
    case  9: /* 1001 */
      rtn = WlzCMeshFMarCompute3D2(nod0, nod3, nod1, nod2, distances);
      break;
    case 10: /* 1010 */
      rtn = WlzCMeshFMarCompute3D2(nod1, nod3, nod0, nod2, distances);
      break;
    case 11: /* 1011 */
      rtn = WlzCMeshFMarCompute3D3(nod0, nod1, nod3, nod2, distances);
      break;
    case 12: /* 1100 */
      rtn = WlzCMeshFMarCompute3D2(nod2, nod3, nod0, nod1, distances);
      break;
    case 13: /* 1101 */
      rtn = WlzCMeshFMarCompute3D3(nod0, nod2, nod3, nod1, distances);
      break;
    case 14: /* 1110 */
      rtn = WlzCMeshFMarCompute3D3(nod1, nod2, nod3, nod0, distances);
      break;
    default:
      break;
  }
  return(rtn);
}

/*!
* \return	Non zero if distance computed and less than current distance.
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element.
* 		This function is given just one known node and computes the
* 		times for the other nodes by assuming that the propagation
* 		is along the edges of the element.
* \param	nod0			Known (current) node.
* \param	nod1			First unknown node.
* \param	nod2			Second unknown node.
* \param	nod3			Third unknown node.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static int	WlzCMeshFMarCompute3D1(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  int		idx,
  		rtn = 0;
  double	d;
  double	*dP;
  WlzDVertex3	del;
  WlzCMeshNod3D *uNodes[3];

  uNodes[0] = nod1;
  uNodes[1] = nod2;
  uNodes[2] = nod3;
  for(idx = 0; idx < 3; ++idx)
  {
    WLZ_VTX_3_SUB(del, uNodes[idx]->pos, nod0->pos);
    d = *(distances + nod0->idx) + WLZ_VTX_3_LENGTH(del);
    dP = distances + uNodes[idx]->idx;
    if(*dP > d)
    {
      *dP = d;
      rtn = 1;
    }
  }
  return(rtn);
}

/*!
* \return	Non zero if distance computed and less than current distance.
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element along the faces of the element.
* \param	nod0			First known (current) node.
* \param	nod1			Second known node.
* \param	nod2			First unknown node.
* \param	nod3			Second unknown node.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static int	WlzCMeshFMarCompute3D2(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  int		rtn = 0;
  double	d;
  WlzDVertex2	q0,
  		q1,
		q2;

  WLZ_VTX_2_SET(q0, 0.0, 0.0);
  WlzGeomMap3DTriangleTo2D(nod0->pos, nod1->pos, nod2->pos, &q1, &q2);
  d = WlzCMeshFMarSolve2D2(q0, q1, q2,
                          distances[nod0->idx], distances[nod1->idx]);
  if(d < distances[nod2->idx])
  {
    distances[nod2->idx] = d;
    ++rtn;
  }
  WlzGeomMap3DTriangleTo2D(nod0->pos, nod1->pos, nod3->pos, &q1, &q2);
  d = WlzCMeshFMarSolve2D2(q0, q1, q2,
                          distances[nod0->idx], distances[nod1->idx]);
  if(d < distances[nod3->idx])
  {
    distances[nod3->idx] = d;
    ++rtn;
  }
  return(rtn);
}

/*!
* \return	Non zero if distance computed and less than current distance.
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element.
* 		This function is given three known nodes and computes the
* 		times for the remaining node using the following method.
* 		The solution is similar to that in "Fast Sweeping Methods
* 		For Eikonal equations On triangular meshes", Jianliang Qian,
* 		etal, SIAM journal on Mumerical Analysis, Vol 45, pp 83-107,
* 		2007. But is given in greater detail here than in the paper.
* 		Given a tetrahedron with four nodes \f$(n_0, n_1, n_2, n_3)\f$
* 		and known front arival times at the first three of these
* 		nodes \f$(t_0, t_1, t_2)\f$ such that \f$t0 \leq t_1, t_2\f$.
* 		This method is only used provided causality constraints are
* 		satisfied: ie \f$n_1 - n_0 > (t_1 - t_0)s_3\f$,
* 		\f$n_2 - n_0 > (t_2 - t_0)s_3\f$ and the normal (see below)
* 		passes through the triangle fromed by \f$n_0, n_1, n_2\f$.
*
* 		The normal vector of the propagating front through \f$n_3\f$
* 		is given by solving:
* 		\f[
		\overline{n_0 n_1} . \mathbf{n} = \frac{t_1 - t_0}{s_3},
		\overline{n_0 n_2} . \mathbf{n} = \frac{t_2 - t_0}{s_3},
		|\mathbf{n}| = 1
                \f]
*		for \f$n\f$, the normal vector of the front which has speed
*		\f$s_3\f$ at \f$n_3\f$.
*		Using:
*               \f[
*               \mathbf{l_1} = \overline{n_0 n_1},
		\mathbf{l_2} = \overline{n_0 n_2},
		d_1 = \frac{t_1 - t_0}{s_3},
		d_2 = \frac{t_2 - t_0}{s_3}
		a   =   l_{1y} l_{2x} - l_{1x} l_{2y}
		b   =   l_{1z} l_{2y} - l_{1y} l_{2z}
		c   =   l_{1z} l_{2x} - l_{1x} l_{2z}
		d   = - d_1 l_{2y}    + d_2 l_{1y}
		e   = - d_1 l_{2x}    + d_2 l_{1x}
                \f]
*		and solving the quadratic gives
*		\f[
		n_z = \frac{-(d b + e c) \pm
		            \sqrt{(d b + e c)^2 - 
			          (d^2 + e^2 - a^2)(a^2 + b^2 + c^2)}}
		           {d^2 + e^2 - a^2}
		n_x =   \frac{d + b n_z}{a}
		n_y = - \frac{e + c n_z}{a}
		\f]
*		With the normal of the advancing front at \f$n_3\f$
*		\f$\mathbf{n}\f$ found. The intersection of this vector
*		with the triangle formed by vertices
*		\f$(\mathbf{P_0}, \mathbf{P_1}, \mathbf{P_2})\f$ corresponding
*		to the nodes\f$(n_0, n_1, n_2)\f$ at vertex \f$\mathbf{Q}\f$
*		determines the causality of the solution. For causality
*		\f$\mathbf{Q}\f$ must be within the triangle.  If this is
*		not satisfied the time value is the minimum for the path
*		along the other three faces.
* \param	nod0			First (current) known node.
* \param	nod1			Second known node.
* \param	nod2			Third known node.
* \param	nod3			Unknown node.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static int	WlzCMeshFMarCompute3D3(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  int		id0,
  		id1,
		hit = 0,
		par = 0,
		rtn = 0;
  double	a,
		a2,
  		b,
		b2,
		c,
		c2,
		d,
		d1,
  		d2,
		e,
		e2,
		f,
		f2,
		g,
		h;
  WlzDVertex2	q0,
  		q1,
		q2;
  WlzDVertex3	n0,
		n1,
		t0,
		t1;
  WlzDVertex3	l[4];
  WlzCMeshNod3D	*tNod;
  WlzCMeshNod3D	*nod[8];

  /* Sort nodes 0 - 2, by time st distances[nod[0]->idx] <=
   * distances[nod[1]->idx] <= distances[nod[2]->idx]. */
  nod[0] = nod0;
  nod[1] = nod1;
  nod[2] = nod2;
  nod[3] = nod3;
  for(id0 = 0; id0 < 3; ++id0)
  {
    for(id1 = id0 + 1; id1 < 3; ++id1)
    {
      if(*(distances + nod[id1]->idx) < *(distances + nod[id0]->idx))
      {
        tNod = nod[id0]; nod[id0] = nod[id1]; nod[id1] = tNod;
      }
    }
  }
  /* Compute vectors and distances relative to nod[0]. */
  WLZ_VTX_3_SUB(l[1], nod[1]->pos, nod[0]->pos);
  WLZ_VTX_3_SUB(l[2], nod[2]->pos, nod[0]->pos);
  d1 = *(distances + nod[1]->idx) - *(distances + nod[0]->idx);
  d2 = *(distances + nod[2]->idx) - *(distances + nod[0]->idx);
  a = WLZ_VTX_3_LENGTH(l[1]);
  b = WLZ_VTX_3_LENGTH(l[2]);
  if((a < d1) && (b < d2))
  {
    /* Compute the unit vector which is normal to the propagation front
     * by solving:
     *   l[1] . n = d1
     *   l[2] . n = d2
     *   n . n  = 1
     */
    a  = l[1].vtY * l[2].vtX - l[1].vtX * l[2].vtY;
    b  = l[1].vtZ * l[2].vtY - l[1].vtY * l[2].vtZ;
    c  = l[1].vtZ * l[2].vtX - l[1].vtX * l[2].vtZ;
    d  = d2 * l[1].vtY - d1 * l[2].vtY;
    e  = d2 * l[1].vtX - d1 * l[2].vtX;
    f  = -(d * b + e * c);
    a2 = a * a;
    b2 = b * b;
    c2 = c * c;
    d2 = d * d;
    e2 = e * e;
    f2 = f * f;
    g  = sqrt(f2 - (d2 + e2 - a2) * (a2 + b2 + c2));
    h  = d2 + e2 - a2;
    n0.vtZ = (f + g) / h;
    n1.vtZ = (f - g) / h;
    n0.vtY = -(e + c * n0.vtZ) / a;
    n1.vtY = -(e + c * n1.vtZ) / a;
    n0.vtX =  (d + b * n0.vtZ) / a;
    n1.vtX =  (d + b * n1.vtZ) / a;
    /* Have two solutions for the normal: n0 and n1, choose the one that runs
     * from the centre of the triangle formed by nodes 0, 1 and 2 to node 3. */
    WLZ_VTX_3_ADD3(t0, nod[0]->pos, nod[1]->pos, nod[2]->pos);
    WLZ_VTX_3_SCALE(t0, t0, 1.0 / 3.0);
    WLZ_VTX_3_SUB(t1, nod[3]->pos, t0);
    a = WLZ_VTX_3_DOT(n0, t1);
    if(a < 0)
    {
      n0 = n1;
    }
    hit = WlzGeomLineTriangleIntersect3D(nod[3]->pos, n0,
					nod[0]->pos, nod[1]->pos, nod[2]->pos,
					&par, NULL, NULL, NULL);
    if(par != 0)
    {
      hit = 0;
    }
    else if(hit != 0)
    {
      /* Normal is through the triangle (nodes 0, 1 and 2), so compute the
       * distance at node 3: t_3 = t_0 + n . (n_3 - n_0). */
      WLZ_VTX_3_SUB(l[3], nod[3]->pos, nod[0]->pos);
      d = WLZ_VTX_3_DOT(n0, l[3]);
      if(d > 0.0)
      {
	d = *(distances + nod[0]->idx) + d;
	if(*(distances + nod[3]->idx) > d)
	{
	  *(distances + nod[3]->idx) = d;
	  rtn = 1;
	}
      }
    }
  }
  else
  {
    /* Find minimum distance using path along the faces of the tetrahedron. */
    WLZ_VTX_2_SET(q0, 0.0, 0.0);
    for(id0 = 0; id0 < 3; ++id0)
    {
      id1 = (id0 + 1) % 3;
      WlzGeomMap3DTriangleTo2D(nod[id0]->pos, nod[id1]->pos,
                               nod[3]->pos, &q1, &q2);
      d = WlzCMeshFMarSolve2D2(q0, q1, q2,
                          distances[nod[id0]->idx],
			  distances[nod[id1]->idx]);
      if(d < distances[nod[3]->idx])
      {
	distances[nod[3]->idx] = d;
	hit = 1;
	rtn = 1;
      }
    }
  }
  if(hit == 0)
  {
    /* TIf all else fails use the minimum distance along the edges. */
    for(id0 = 0; id0 < 3; ++id0)
    {
      WLZ_VTX_3_SUB(t0, nod[3]->pos, nod[id0]->pos);
      d = WLZ_VTX_3_LENGTH(t0) + *(distances + nod[id0]->idx);
      if(*(distances + nod[3]->idx) > d)
      {
	rtn = 1;
	*(distances + nod[3]->idx) = d;
      }
    }
  }
  return(rtn);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given seeds with zero distance, setting known
* 		distances and initialising the node queue.
* \param	nodQ			The node queue.
* \param	mesh			The constrained mesh.
* \param	edgQMin			Initial number of entries for the
* 					edge queue. Appropriate value can
* 					avoid reallocation.
* \param	distances		Array of distances.
* \param	nSeed			Number of seeds.
* \param	seeds			Array of seeds.
*/
static WlzErrorNum WlzCMeshFMarAddSeeds2D(AlcHeap *nodQ,
				WlzCMesh2D *mesh, int edgQMin,
				double *distances,
				int nSeeds, WlzDVertex2 *seeds)
{
  int		idN,
  		idS,
		hit;
  WlzCMeshNod2D	*nod0,
  		*nod1;
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  AlcHeap	*sElmQ = NULL;
  WlzUByte	*eFlgs = NULL;
  double	*sDists = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idS = 0;
  if(((sDists = (double *)
                AlcMalloc(sizeof(double) * mesh->res.nod.maxEnt)) == NULL) ||
     ((eFlgs = (WlzUByte *)
                AlcMalloc(sizeof(WlzUByte) * mesh->res.elm.maxEnt)) == NULL) ||
     ((sElmQ = AlcHeapNew(sizeof(WlzCMeshFMarQEnt), edgQMin, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* For each seed: Update the Euclidean distances. */
    while((errNum == WLZ_ERR_NONE) && (idS < nSeeds))
    {
      errNum = WlzCMeshFMarAddSeed2D(sElmQ, mesh,  distances, sDists,
                                     seeds[idS], eFlgs);
      AlcHeapAllEntFree(sElmQ, 0);
      ++idS;
    }
  }
  AlcFree(sDists);
  AlcFree(eFlgs);
  AlcHeapFree(sElmQ);
  /* Add nodes which have known distance but are not surrounded by nodes
   * with known distance to the node queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      if(distances[idN] < DBL_MAX / 2.0)
      {
        nod0 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	nod0->flags |= WLZ_CMESH_NOD_FLAG_KNOWN;
	edu1 = edu0 = nod0->edu;
	hit = 0;
	do
	{
	  nod1 = edu1->next->nod;
	  if(distances[nod1->idx] > DBL_MAX / 2.0)
	  {
	    hit = 1;
	    break;
	  }
	  edu1 = edu1->nnxt;
	} while(edu1 != edu0);
	if(hit == 0)
	{
	  nod0->flags |= WLZ_CMESH_NOD_FLAG_UPWIND;
	}
	else
	{
          if((errNum = WlzCMeshFMarQInsertNod2D(nodQ, nod0,
					distances[nod0->idx])) != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given seeds with zero distance, setting known
* 		distances and initialising the node queue.
* \param	nodQ			The node queue.
* \param	mesh			The constrained mesh.
* \param	fceQMin			Initial number of entries for the
* 					face queue. Appropriate value can
* 					avoid reallocation.
* \param	distances		Array of distances.
* \param	nSeed			Number of seeds.
* \param	seeds			Array of seeds.
*/
static WlzErrorNum WlzCMeshFMarAddSeeds3D(AlcHeap *nodQ,
				WlzCMesh3D *mesh, int fceQMin,
				double *distances,
				int nSeeds, WlzDVertex3 *seeds)
{
  int		idN,
  		idS,
		hit;
  WlzCMeshNod3D	*nod0,
  		*nod1;
  WlzCMeshEdgU3D *edu0,
  		*edu1;
  AlcHeap	*sElmQ = NULL;
  WlzUByte	*eFlgs = NULL;
  double	*sDists = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idS = 0;
  if(((sDists = (double *)
                AlcMalloc(sizeof(double) * mesh->res.nod.maxEnt)) == NULL) ||
     ((eFlgs = (WlzUByte *)
                AlcMalloc(sizeof(WlzUByte) * mesh->res.elm.maxEnt)) == NULL) ||
     ((sElmQ = AlcHeapNew(sizeof(WlzCMeshFMarQEnt), fceQMin, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* For each seed: Update the Euclidean distances. */
    while((errNum == WLZ_ERR_NONE) && (idS < nSeeds))
    {
      errNum = WlzCMeshFMarAddSeed3D(sElmQ, mesh,  distances, sDists,
                                     seeds[idS], eFlgs);
      AlcHeapAllEntFree(sElmQ, 0);
      ++idS;
    }
  }
  AlcFree(sDists);
  AlcFree(eFlgs);
  AlcHeapFree(sElmQ);
  /* Add nodes which have known distance but are not surrounded by nodes
   * with known distance to the node queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      if(distances[idN] < DBL_MAX / 2.0)
      {
        nod0 = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	if(nod0->idx >= 0)
	{
	  nod0->flags |= WLZ_CMESH_NOD_FLAG_KNOWN;
	  edu1 = edu0 = nod0->edu;
	  hit = 0;
	  do
	  {
	    nod1 = edu1->next->nod;
	    if(distances[nod1->idx] > DBL_MAX / 2.0)
	    {
	      hit = 1;
	      break;
	    }
	    edu1 = edu1->nnxt;
	  } while(edu1 != edu0);
	  if(hit == 0)
	  {
	    nod0->flags |= WLZ_CMESH_NOD_FLAG_KNOWN |
	                   WLZ_CMESH_NOD_FLAG_UPWIND;
	  }
	  else
	  {
	    if((errNum = WlzCMeshFMarQInsertNod3D(nodQ, nod0,
					distances[nod0->idx])) != WLZ_ERR_NONE)
	    {
	      break;
	    }
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds a single seed to the element queue.
* \param	sElmQ			The queue to use for elements while
* 					propagating out the region within
* 					which all node to seed straight
* 					line paths are within the mesh.
* \param	mesh			The mesh.
* \param	distances		Minimum distances from all seeds.
* \param	sDst			Distances from this seed.
* \param	seed			A single seed position.
* \param	eFlgs			Array of element flags which are
* 					non zero when element has been
* 					visited.
*/
static WlzErrorNum WlzCMeshFMarAddSeed2D(AlcHeap  *sElmQ,
                                         WlzCMesh2D *mesh,
					 double *distances, double *sDst,
					 WlzDVertex2 seed,
					 WlzUByte *eFlgs)
{
  int 		idE,
  		idF,
		idM,
		ilos;
  double	d;
  WlzDVertex2	del;
  WlzCMeshNod2D *nod0;
  WlzCMeshEdgU2D *edu0,
  		*edu1,
		*edu2;
  WlzCMeshElm2D *elm0,
  		*elm2;
  WlzCMeshFMarQEnt *sElmQEnt;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find any element which encloses the seed (there may be more than one
   * if the seed is on an element's edge or at a node. */
  idM = WlzCMeshElmEnclosingPos2D(mesh, -1, seed.vtX, seed.vtY, 0, NULL);
  if(idM < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    elm0 = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idM);
    if(elm0->idx < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueSetDouble(sDst, DBL_MAX, mesh->res.nod.maxEnt);
    WlzValueSetUByte(eFlgs, 0, mesh->res.elm.maxEnt);
    /* Compute the distances of this element's nodes from the seed and
     * update the minimum distances, then initialize the edge queue
     * using this elements edges. */
    for(idE = 0; idE < 3; ++idE)
    {
      edu0 = elm0->edu + idE;
      nod0 = edu0->nod;
      WLZ_VTX_2_SUB(del, seed, nod0->pos);
      d = WLZ_VTX_2_LENGTH(del);
      sDst[nod0->idx] = d;
      if(distances[nod0->idx] > d)
      {
	distances[nod0->idx] = d;
      }
      if((edu0->opp != NULL) && (edu0->opp != edu0))
      {
	if((errNum = WlzCMeshFMarSElmQInsert2D(sElmQ, edu0->opp->elm,
					       sDst, seed)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    eFlgs[elm0->idx] = 1;
  }
  /* Pop element that has the min(maximum node distance) from the edge
   * queue. */
  while((errNum == WLZ_ERR_NONE) &&
        ((sElmQEnt = (WlzCMeshFMarQEnt *)AlcHeapTop(sElmQ)) != NULL))
  {
    elm0 = (WlzCMeshElm2D *)(sElmQEnt->entity);
    AlcHeapEntFree(sElmQ);
    /* Find node with unknown distance for the element. */
    for(idE = 0; idE < 3; ++idE)
    {
      edu0 = elm0->edu + idE;
      nod0 = edu0->nod;
      if(sDst[nod0->idx] > DBL_MAX / 2.0)
      {
        break;
      }
    }
    if(idE >= 3)
    {
      ilos = 1;
    }
    else
    {
      /* If the line segment from the node with unknown distance to the
       * seed passes through this element, look at the element which
       * shares the edge opposide to this node to see if it is in line
       * of sight to the seed.
       * If the line segment does not pass through this element, then
       * check all elements which use the node to see if they are in
       * line of sight. */
      ilos = 0;
      if(WlzGeomLineSegmentsIntersect(seed, nod0->pos,
                                      edu0->next->nod->pos,
				      edu0->next->next->nod->pos, NULL) > 0)
      {
	/* Line segment passes though this element. */
        edu2 = edu0->next;
	if((edu2->opp != NULL) && (edu2->opp != edu2))
	{
	  elm2 = edu2->opp->elm;
	  if((sDst[elm2->edu[0].nod->idx] < DBL_MAX / 2.0) &&
	     (sDst[elm2->edu[1].nod->idx] < DBL_MAX / 2.0) &&
	     (sDst[elm2->edu[2].nod->idx] < DBL_MAX / 2.0))
	  {
	    ilos = 1;
	  }
	}
      }
      else
      {
        edu1 = edu0->nnxt;
	while((ilos == 0) && (edu1 != edu0))
	{
	  if(WlzGeomLineSegmentsIntersect(seed, nod0->pos,
	                                  edu1->next->nod->pos,
					  edu1->next->next->nod->pos, NULL) > 0)
	  {
	    /* Line segment passes though element with this edge use. */
	    edu2 = edu1->next;
	    if((edu2->opp != NULL) && (edu2->opp != edu2))
	    {
	      elm2 = edu2->opp->elm;
	      if((sDst[elm2->edu[0].nod->idx] < DBL_MAX / 2.0) &&
	         (sDst[elm2->edu[1].nod->idx] < DBL_MAX / 2.0) &&
		 (sDst[elm2->edu[2].nod->idx] < DBL_MAX / 2.0))
	      {
	        ilos = 1;
	      }
	    }
	  }
	  edu1 = edu1->nnxt;
	}
      }
    }
    if(ilos != 0)
    {
      eFlgs[elm0->idx] = 1;
      WLZ_VTX_2_SUB(del, seed, nod0->pos);
      d = WLZ_VTX_2_LENGTH(del);
      sDst[nod0->idx] = d;
      if(distances[nod0->idx] > d)
      {
	distances[nod0->idx] = d;
      }
      for(idF = 0; idF < 3; ++idF)
      {
	edu1 = elm0->edu + idF;
	if((edu1->opp != NULL) && (edu1->opp != edu1) &&
	   (eFlgs[edu1->opp->elm->idx] == 0))
	{
	  if((errNum = WlzCMeshFMarSElmQInsert2D(sElmQ, edu1->opp->elm,
					     sDst, seed)) != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds a single seed to the element queue.
* \param	sElmQ			The queue to use for elements while
* 					propagating out the region within
* 					which all node to seed straight
* 					line paths are within the mesh.
* \param	mesh			The mesh.
* \param	distances		Minimum distances from all seeds.
* \param	sDst			Distances from this seed.
* \param	seed			A single seed position.
* \param	eFlgs			Array of element flags which are
* 					non zero when element has been
* 					visited.
*/
static WlzErrorNum WlzCMeshFMarAddSeed3D(AlcHeap  *sElmQ,
                                         WlzCMesh3D *mesh,
					 double *distances, double *sDst,
					 WlzDVertex3 seed,
					 WlzUByte *eFlgs)
{
  int 		idE,
  		idF,
		idM,
		idN,
		ilos;
  WlzDVertex3	dir;
  WlzCMeshNod3D *nod0;
  WlzCMeshFace	*fce0;
  WlzCMeshElm3D *elm0,
  		*elm2;
  WlzCMeshFMarQEnt *sElmQEnt;
  WlzCMeshNod3D *nodes[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find any element which encloses the seed (there may be more than one
   * if the seed is on an element's edge or at a node. */
  idM = WlzCMeshElmEnclosingPos3D(mesh, -1, seed.vtX, seed.vtY, seed.vtZ,
  			          0, NULL);
  if(idM < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    elm0 = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idM);
    if(elm0->idx < 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueSetDouble(sDst, DBL_MAX, mesh->res.nod.maxEnt);
    WlzValueSetUByte(eFlgs, 0, mesh->res.elm.maxEnt);
    /* Compute the distances of this element's nodes from the seed and
     * update the minimum distances, then initialize the face queue
     * using this elements faces. */
    WlzCMeshElmGetNodes3D(elm0, nodes + 0, nodes + 1, nodes + 2, nodes + 3);
    for(idN = 0; idN < 4; ++idN)
    {
      nod0 = nodes[idN];
      sDst[nod0->idx] = WlzGeomDist3D(seed, nod0->pos);
      if(distances[nod0->idx] > sDst[nod0->idx])
      {
	distances[nod0->idx] = sDst[nod0->idx];
      }
    }
    for(idN = 0; idN < 4; ++idN)
    {
      fce0 = elm0->face + 3 - idN;
      if((fce0->opp != NULL) && (fce0->opp != fce0))
      {
	if((errNum = WlzCMeshFMarSElmQInsert3D(sElmQ, fce0->opp->elm,
					       sDst, seed)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    eFlgs[elm0->idx] = 1;
  }
  /* Pop element that has the min(maximum node distance) from the edge
   * queue. */
  while((errNum == WLZ_ERR_NONE) &&
        ((sElmQEnt = (WlzCMeshFMarQEnt *)AlcHeapTop(sElmQ)) != NULL))
  {
    ilos = 0;
    elm0 = (WlzCMeshElm3D *)(sElmQEnt->entity);
    AlcHeapEntFree(sElmQ);
    WlzCMeshElmGetNodes3D(elm0, nodes + 0, nodes + 1, nodes + 2, nodes + 3);
    /* Find node with unknown distance for the element. */
    for(idN = 0; idN < 4; ++idN)
    {
      if(sDst[nodes[idN]->idx] > DBL_MAX / 2.0)
      {
        break;
      }
    }
    if(idN >= 4)
    {
      ilos = 4;
    }
    else
    {
      /* If the line segment from the node with unknown distance to the
       * seed passes through this element, look at the element which
       * shares the edge opposide to this node to see if it is in line
       * of sight to the seed.
       * If the line segment does not pass through this element, then
       * check all elements which use the node to see if they are in
       * line of sight. */
      ilos = 0;
      idF = 3 - idN;
      nod0 = nodes[idN];
      fce0 = elm0->face + idF;
      WLZ_VTX_3_SUB(dir, seed, nod0->pos);
      if(WlzGeomLineTriangleIntersect3D(seed, dir,
                                       fce0->edu[0].nod->pos,
                                       fce0->edu[1].nod->pos,
                                       fce0->edu[2].nod->pos,
                                       NULL, NULL, NULL, NULL) > 0)
      {
	/* Line segment passes through this element. */
	if((fce0->opp != NULL) && (fce0->opp != fce0))
	{
	  elm2 = fce0->opp->elm;
	  if((sDst[elm2->face[0].edu[0].nod->idx] < DBL_MAX / 2.0) &&
	     (sDst[elm2->face[0].edu[1].nod->idx] < DBL_MAX / 2.0) &&
	     (sDst[elm2->face[0].edu[2].nod->idx] < DBL_MAX / 2.0) &&
	     (sDst[elm2->face[1].edu[1].nod->idx] < DBL_MAX / 2.0))
	  {
	    ilos = 1;
	  }
	}
      }
      else
      {
	idE = idF;
	while((idE = (idE + 1) % 4) != idF)
	{
	  fce0 = elm0->face + idF;
	  if(WlzGeomLineTriangleIntersect3D(seed, dir,
	                                   fce0->edu[0].nod->pos,
					   fce0->edu[1].nod->pos,
					   fce0->edu[2].nod->pos,
					   NULL, NULL, NULL, NULL) > 0)
	  {
	    elm2 = fce0->opp->elm;
	    if((sDst[elm2->face[0].edu[0].nod->idx] < DBL_MAX / 2.0) &&
	       (sDst[elm2->face[0].edu[1].nod->idx] < DBL_MAX / 2.0) &&
	       (sDst[elm2->face[0].edu[2].nod->idx] < DBL_MAX / 2.0) &&
	       (sDst[elm2->face[1].edu[1].nod->idx] < DBL_MAX / 2.0))
	    {
	      ilos = 1;
	      break;
	    }
	  }
	}
      }
    }
    if(ilos != 0)
    {
      eFlgs[elm0->idx] = 1;
      if(ilos < 4)
      {
	sDst[nod0->idx] = WlzGeomDist3D(seed, nod0->pos);
	if(distances[nod0->idx] > sDst[nod0->idx])
	{
	  distances[nod0->idx] = sDst[nod0->idx];
	}
      }
      for(idF = 0; idF < 4; ++idF)
      {
	fce0 = elm0->face + idF;
	if((fce0->opp != NULL) && (fce0->opp != fce0) &&
	   (eFlgs[fce0->opp->elm->idx] == 0))
	{
	  if((errNum = WlzCMeshFMarSElmQInsert3D(sElmQ, fce0->opp->elm,
					         sDst, seed)) != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts an element into the 2D seed element queue.
* \param	queue			The seed element queue.
* \param	elm			Element use to insert into the queue.
* \param	dst			Distances for the seed.
* \param	org			Seed for Euclidean distances.
*/
static WlzErrorNum WlzCMeshFMarSElmQInsert2D(AlcHeap *queue,
					     WlzCMeshElm2D *elm,
					     double *dst,
					     WlzDVertex2 org)
{
  WlzCMeshFMarQEnt ent;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(elm != NULL)
  {
    /* Insert entry into the queue and add hash table entry. */
    ent.priority = WlzCMeshFMarQSElmPriority2D(elm, dst, org);
    ent.entity = elm;
    if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts an element into the 3D seed element queue.
* \param	queue			The seed element queue.
* \param	elm			Element use to insert into the queue.
* \param	dst			Distances for the seed.
* \param	org			Seed for Euclidean distances.
*/
static WlzErrorNum WlzCMeshFMarSElmQInsert3D(AlcHeap *queue,
					     WlzCMeshElm3D *elm,
					     double *dst,
					     WlzDVertex3 org)
{
  WlzCMeshFMarQEnt ent;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  if(elm != NULL)
  {
    /* Insert entry into the queue and add hash table entry. */
    ent.priority = WlzCMeshFMarQSElmPriority3D(elm, dst, org);
    ent.entity = elm;
    if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Clears the existing queue entries and then, allocating
* 		entries as required, populates the queue. The entries
* 		are formed from all elements which use the given node
* 		have at least 1 known nodes and at least one node which
* 		is not upwind.
* \param	queue			Element queue to use/reuse.
* \param	nod			Node with which to populate the
* 					queue.
*/
static WlzErrorNum WlzCMeshFMarElmQInit2D(AlcHeap *queue, WlzCMeshNod2D *nod)
{
  int		priority;
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  WlzCMeshFMarElmQEnt ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  queue->nEnt = 0;
  queue->data = nod;
  edu0 = edu1 = nod->edu;
  do
  {
    ent.elm.e2 = edu1->elm;
    priority = WlzCMeshFMarElmQCalcPriority2D(edu1->elm, nod);
    if(priority > 0)
    {
      ent.priority = priority;
      if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
    }
    edu1 = edu1->nnxt;
  } while(edu1 != edu0);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Clears the existing queue entries and then, allocating
* 		entries as required, populates the queue. The entries
* 		are formed from all elements which use the given node
* 		have at least 1 known nodes and at least one node which
* 		is not upwind.
* \param	queue			Element queue to use/reuse.
* \param	nod			Node with which to populate the
* 					queue.
*/
static WlzErrorNum WlzCMeshFMarElmQInit3D(AlcHeap *queue, WlzCMeshNod3D *nod)
{
  int		priority;
  WlzCMeshEdgU3D *edu0,
  		*edu1;
  WlzCMeshFMarElmQEnt ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  queue->nEnt = 0;
  queue->data = nod;
  edu0 = edu1 = nod->edu;
  do
  {
    ent.elm.e3 = edu1->face->elm;
    priority = WlzCMeshFMarElmQCalcPriority3D(edu1->face->elm, nod);
    if(priority > 0)
    {
      ent.priority = priority;
      if(AlcHeapInsertEnt(queue, &ent) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
	break;
      }
    }
    edu1 = edu1->nnxt;
  } while(edu1 != edu0);
  return(errNum);
}

/*!
* \return	Element priority value.
* \return	Element priority.
* \ingroup	W;zMesh
* \brief	Computes the priority value of an element queue entry
* 		given the element and it's current node about which
* 		the elements are clustered.
*		The priority given by:
*		\f$p = \sum_{i=0}^{2}{p_i}\f$
*		where
*		\f$p_i = 2\f$ if the node is upwind or the current node,
*		\f$p_i = 1\f$ if the node is an other known node
*		\f$p_i = 0\f$ if the node is unknown (downwind).
* \param	elm			Given element.
* \param	cNod			The current node around wich the
* 					elements in the queue are clustered.
*/
static int	WlzCMeshFMarElmQCalcPriority2D(WlzCMeshElm2D *elm,
				               WlzCMeshNod2D *cNod)
{
  int		idx,
  		priority = 6;
  WlzCMeshNod2D *nod;

  for(idx = 0; idx < 3; ++idx)
  {
    nod = elm->edu[idx].nod;
    priority -= ((nod->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) +
                ((nod->idx == cNod->idx) ||
		 ((nod->flags & WLZ_CMESH_NOD_FLAG_UPWIND) != 0));
  }
  return(priority);
}

/*!
* \return	Element priority value.
* \ingroup	WlzMesh
* \brief	Computes the priority value of an element queue entry
* 		given the element and it's current node about which
* 		the elements are clustered.
*		The priority given by:
*		\f$p = 8 - \sum_{i=0}^{2}{p_i}\f$
*		where
*		\f$p_i = 2\f$ if the node is upwind or the current node,
*		\f$p_i = 1\f$ if the node is an other known node
*		\f$p_i = 0\f$ if the node is unknown (downwind).
* \param	elm			The element.
* \param	cNod			The current node around wich the
* 					elements in the queue are clustered.
*/
static int	WlzCMeshFMarElmQCalcPriority3D(WlzCMeshElm3D *elm,
				WlzCMeshNod3D *cNod)
{
  int		idx,
  		priority = 8;
  WlzCMeshNod3D *nod;

  for(idx = 0; idx < 3; ++idx)
  {
    nod = elm->face[0].edu[idx].nod;
    priority -= ((nod->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) +
                ((nod->idx == cNod->idx) ||
		 ((nod->flags & WLZ_CMESH_NOD_FLAG_UPWIND) != 0));
  }
  nod = elm->face[1].edu[1].nod;
  priority -= ((nod->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) +
              ((nod->idx == cNod->idx) ||
	       ((nod->flags & WLZ_CMESH_NOD_FLAG_UPWIND) != 0));
  return(priority);
}
