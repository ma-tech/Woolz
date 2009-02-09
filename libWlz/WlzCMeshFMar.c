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

/*!
* \struct	_WlzCMeshFMarQEnt
* \ingroup	WlzMesh
* \brief	An entry of a WlzCMeshFMarQ queue.
* 		Typedef: ::WlzCMeshFMarQEnt.
*/
typedef struct _WlzCMeshFMarQEnt
{
#ifdef DEBUG
  int			dbgIdx;		/*!< For debug only. */
#endif
  int			next;		/*!< Index of next towards tail,
  					     -ve iff at tail. */
  int			prev;		/*!< Index of next towards head,
                                             -ve iff at head. */
  int			hashNxt;	/*!< Index of next in hash bucket. */
  void			*entity;	/*!< Pointer to mesh entity. */
  double		priority;	/*!< Entry priority highest priority
  					     at the head of the queue. */
} WlzCMeshFMarQEnt;

/*!
* \struct	_WlzCMeshFMarQ
* \ingroup	WlzMesh
* \brief	A priority queue for the nodes of an advancing front in
* 		WlzCMesh based fast marching algorithms.
* 		Typedef: ::WlzCMeshFMarQ.
*
* 		The queue is encoded as a linked list of entries with
* 		fast access to the highest and lowest priority  entries
* 		via the head/tail fields. There is also fast random
* 		access to the entries for a node via hash based buckets.
* 		A key assumption, that is true for it's intended
* 		purpose (a queue of front nodes in fast marchinga),
* 		is that a change to a single node can only change the
* 		priority of the entry holding that node. Because of
* 		this it is possible to keep the queue sorted through
* 		incremental actions.
*/
typedef struct _WlzCMeshFMarQ
{
  int			head;           /*!< Index of the queue head, with
                                             the highest priority entry at
					     the head. */
  int			tail;           /*!< Index of the queue tail, with
  					     the lowest priority entry at
					     the tail. */
  int			last;		/*!< Index of last entry inserted. */
  int			cnt;		/*!< Number of entries in use:
                                             incremented when entry inserted,
					     no change if unlinked,
					     decremented when entry freed. */
  int			max;		/*!< Number of entries allocated. */
  int 			free;		/*!< Index of first free entry, rest
                                             via next index. */
  WlzCMeshFMarQEnt	*entries;	/*!< Array of allocated entries. */
  int			*buckets;	/*!< Indices for hash buckets, used
  					     fast for access to entry by node
					     index. */
  int			(*hashFn)(struct _WlzCMeshFMarQ *, void *);
  					/*!< Hash function which when passed
					     an entry entity pointer returns
					     a hash table index in the range
					     [0 - (queue->max - 1)]. */
} WlzCMeshFMarQ;

/*!
* \struct	_WlzCMeshFMarElmQ
* \ingroup	WlzMesh
* \brief	A queue for the elements surrounding a node of an
* 		advancing front in WlzCMesh based fast marching algorithms.
* 		Typedef: ::WlzCMeshFMarElmQ.
*
* 		The queue is maintained sorted by the number of known
* 		nodes in the element. The highest priority element is
* 		the last in the list and this is easily removed by
* 		decrementing the number of entries counter.
*/
typedef struct _WlzCMeshFMarElmQ
{
  int			nEnt;		/*!< Number of entries. */
  int			maxEnt;		/*!< Number of entries space is
  					     allocated for. */
  WlzCMeshNodP		nod;		/*!< Current node around which the
  					     queue is formed. */
  struct _WlzCMeshFMarElmQEnt *entries;	/*!< The queue entries. */
} WlzCMeshFMarElmQ;

/*!
* \struct	_WlzCMeshFMarElmQEnt
* \ingroup	WlzMesh
* \brief	An entry of an element queue.
* 		Typedef: ::WlzCMeshFMarElmQEnt.
*/
typedef struct _WlzCMeshFMarElmQEnt
{
  WlzCMeshElmP		elm;		/*!< Element pointer. */
  WlzCMeshNodP		nod[4];		/*!< Pointers to nodes of element. */
  int			priority;	/*!< Priority of queue entry: The
  					     priority is the simple sum of
					     the priority of the nodes (2
					     for an upwind node or the current
					     node, 1 for any other known node
					     and zero for an unknown (down
					     wind) node). */
} WlzCMeshFMarElmQEnt;

static int			WlzCMeshFMarQHashFn(
				  int value,
				  int maxVal);
static int			WlzCMeshFMarHashFnElm2D(
				  WlzCMeshFMarQ *queue,
				  void *entity);
static int			WlzCMeshFMarHashFnNod2D(
				  WlzCMeshFMarQ *queue,
				  void *entity);
static int			WlzCMeshFMarHashFnNod3D(
				  WlzCMeshFMarQ *queue,
				  void *entity);
static int			WlzCMeshFMarElmQElmIdxCmp2D(
				  const void *p0,
				  const void *p1);
static int			WlzCMeshFMarElmQElmIdxCmp3D(
				  const void *p0,
				  const void *p1);
static int			WlzCMeshFMarElmQPriorityCmp(
				  const void *p0,
				  const void *p1);
static double		   	WlzCMeshFMarQNodPriority2D(
				  WlzCMeshNod2D *nod,
				  double dist);
static double			WlzCMeshFMarQSElmPriority2D(
				  WlzCMeshElm2D *elm,
				  double *dst,
				  WlzDVertex2 org);
static double		   	WlzCMeshFMarQNodPriority3D(
				  WlzCMeshNod3D *nod,
				  double *dist);
static void			WlzCMeshFMarQFree(
				  WlzCMeshFMarQ *queue);
static void			WlzCMeshFMarElmQFree(
				  WlzCMeshFMarElmQ *queue);
static void			WlzCMeshFMarQEntFree(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEnt *qEnt);
static void			WlzCMeshFMarQEntFreeAll(
				  WlzCMeshFMarQ *queue);
static void			WlzCMeshFMarQUnlinkEntFromList(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEnt *ent);
static void			WlzCMeshFMarQUnlinkEntFromHash(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEnt *ent);
static void			WlzCMeshFMarQRehash(
				  WlzCMeshFMarQ *queue);
static void			WlzCMeshFMarQInsertEnt(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEnt *gEnt);
static void		 	WlzCMeshFMarCompute2D(
				  WlzCMeshNod2D *nod2,
				  WlzCMeshNod2D *nod0,
				  WlzCMeshNod2D *nod1,
				  double *distances,
				  WlzCMeshElm2D *elm);
static void			WlzCMeshFMarCompute3D(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static void			WlzCMeshFMarCompute3D1(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static void			WlzCMeshFMarCompute3D2(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static void			WlzCMeshFMarCompute3D3(
                                  WlzCMeshNod3D *nod0,
				  WlzCMeshNod3D *nod1,
				  WlzCMeshNod3D *nod2,
				  WlzCMeshNod3D *nod3,
				  double *distances);
static void			WlzCMeshFMarElmQSqueeze2D(
				  WlzCMeshFMarElmQ *queue);
static void			WlzCMeshFMarElmQSqueeze3D(
				  WlzCMeshFMarElmQ *queue);
static void			WlzCMeshFMarElmQSort2D(
				  WlzCMeshFMarElmQ* queue);
static void			WlzCMeshFMarElmQSort3D(
				  WlzCMeshFMarElmQ* queue);
static void			WlzCMeshFMarElmQCalcPriority2D(
				  WlzCMeshFMarElmQEnt *ent,
				  WlzCMeshNod2D *cNod);
static void			WlzCMeshFMarElmQCalcPriority3D(
				  WlzCMeshFMarElmQEnt *ent,
				  WlzCMeshNod3D *cNod);
static WlzErrorNum 		WlzCMeshFMarQRealloc(
				  WlzCMeshFMarQ *queue,
				  int minEnt);
static WlzErrorNum 		WlzCMeshFMarAddSeeds2D(
				  WlzCMeshFMarQ *queue,
				  WlzCMesh2D *mesh, 
				  int edgQMin,
				  double *distances,
				  int nSeeds,
				  WlzDVertex2 *seeds);
static WlzErrorNum 		WlzCMeshFMarAddSeeds3D(
				  WlzCMeshFMarQ *queue,
				  WlzCMesh3D *mesh,
				  double *distances,
				  int nSeeds,
				  WlzDVertex3 *seeds);
static WlzErrorNum 		WlzCMeshFMarAddSeed2D(
				  WlzCMeshFMarQ  *edgQ,
                                  WlzCMesh2D *mesh,
				  double *distances,
				  double *sDists,
				  WlzDVertex2 seed,
				  WlzUByte *eFlgs);
static WlzErrorNum 		WlzCMeshFMarInsertSeed3D(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshNod3D *nod,
				  double *distances,
				  double *speeds,
				  WlzDVertex3 seedPos);
static WlzErrorNum 		WlzCMeshFMarQInsertNod2D(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshNod2D *nod,
				  double dist);
static WlzErrorNum 		WlzCMeshFMarQInsertNod3D(
				  WlzCMeshFMarQ *queue,
				  double *times,
				  WlzCMeshNod3D *nod);
static WlzErrorNum 		WlzCMeshFMarSElmQInsert2D(
				  WlzCMeshFMarQ *edgQ,
				  WlzCMeshElm2D *elm,
				  double *dst,
				  WlzDVertex2 org);
static WlzErrorNum 		WlzCMeshFMarElmQInit2D(
				  WlzCMeshFMarElmQ *queue,
				  WlzCMeshNod2D *nod);
static WlzErrorNum 		WlzCMeshFMarElmQInit3D(
				  WlzCMeshFMarElmQ *queue,
				  WlzCMeshNod3D *nod);
static WlzCMeshFMarQ 		*WlzCMeshFMarQNew(
				  int nEnt,
				  int (*hashFn)(WlzCMeshFMarQ *, void *));
static WlzCMeshFMarQEnt 	*WlzCMeshFMarQPopHead(
				  WlzCMeshFMarQ *queue);
static WlzCMeshFMarQEnt 	*WlzCMeshFMarQPopTail(
				  WlzCMeshFMarQ *queue);
static WlzCMeshFMarQEnt 	*WlzCMeshFMarQNewEnt(
				  WlzCMeshFMarQ *queue,
				  WlzErrorNum *dstErr);
static WlzCMeshFMarElmQ 	*WlzCMeshFMarElmQNew(
				  void);
static WlzCMeshFMarElmQEnt	*WlzCMeshFMarElmQPopTail2D(
				  WlzCMeshFMarElmQ *queue);
static WlzCMeshFMarElmQEnt	*WlzCMeshFMarElmQPopTail3D(
				  WlzCMeshFMarElmQ *queue);

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
  else if(mesh->res.elm.numEnt < 0)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
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
      obj0 = WlzCMeshToDomObj2D(mesh, &errNum);
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
	                                      pos.vtX, pos.vtY, &idN)) >= 0)
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
  int		
  		idN,
  		idS,
		nBnd;
  WlzCMeshNod2D	*nod0,
                *nod1;
  WlzCMeshFMarQ *nodQ = NULL;
  WlzCMeshFMarElmQ *elmQ = NULL;
  WlzCMeshFMarQEnt *nodQEnt;
  WlzCMeshFMarElmQEnt *elmQEnt = NULL;
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
    if((nBnd = WlzCMeshSetBoundNodFlags2D(mesh)) <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Create and then initialise the active node queue using the given seed
   * or boundary nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodQ = WlzCMeshFMarQNew(nBnd, WlzCMeshFMarHashFnNod2D)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nSeeds > 0)
    {
      errNum = WlzCMeshFMarAddSeeds2D(nodQ, mesh, nBnd + 1,
                                      distances, nSeeds, seeds);
    }
    else
    {
      nSeeds = nBnd;
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
	errNum = WlzCMeshFMarAddSeeds2D(nodQ, mesh, nBnd + 1,
				        distances, nSeeds, seeds);
	AlcFree(seeds);
      }
    }
  }
  /* Create element queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((elmQ = WlzCMeshFMarElmQNew()) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Until the queue is empty: Pop the node with lowest priority (ie
   * lowest distance) from the queue, and process it. */
  while((errNum == WLZ_ERR_NONE) &&
	((nodQEnt = WlzCMeshFMarQPopTail(nodQ)) != NULL))
  {
    /* Find all neighbouring nodes that are neither active nor upwind.
     * For each of these neighbouring nodes, compute their distance, set
     * them to active and insert them into the queue.*/
    nod0 = (WlzCMeshNod2D *)(nodQEnt->entity);
    errNum = WlzCMeshFMarElmQInit2D(elmQ, nod0);
    if(errNum == WLZ_ERR_NONE)
    {
      /* While element list is not empty, remove element and compute all
       * node distances for it. */
      while((elmQEnt = WlzCMeshFMarElmQPopTail2D(elmQ)) != NULL)
      {
	/* Entry nodes must be ordered 
	 *   0 - current node
	 *   1 - node with least distance of remaining nodes
	 *   2 - node with greatest distance of remaining nodes
	 */
	if(distances[elmQEnt->nod[2].n2->idx] <
	   distances[elmQEnt->nod[1].n2->idx])
	{
	  nod1 = elmQEnt->nod[2].n2;
	  elmQEnt->nod[2].n2 = elmQEnt->nod[1].n2;
	  elmQEnt->nod[1].n2 = nod1;
	}
	/* Compute distances. */
	WlzCMeshFMarCompute2D(elmQEnt->nod[2].n2, elmQEnt->nod[0].n2,
			      elmQEnt->nod[1].n2, distances, elmQEnt->elm.e2);
	elmQEnt->nod[2].n2->flags = WLZ_CMESH_NOD_FLAG_KNOWN |
	                            WLZ_CMESH_NOD_FLAG_ACTIVE;
	if((errNum = WlzCMeshFMarQInsertNod2D(nodQ, elmQEnt->nod[2].n2,
			  distances[elmQEnt->nod[2].n2->idx])) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
      /* Set the current node to be upwind. */
      if(errNum == WLZ_ERR_NONE)
      {
	nod0->flags = (nod0->flags & ~(WLZ_CMESH_NOD_FLAG_ACTIVE)) |
		     WLZ_CMESH_NOD_FLAG_UPWIND;
      }
    }
    WlzCMeshFMarQEntFree(nodQ, nodQEnt);
  }
  /* Clear up. */
  WlzCMeshFMarElmQFree(elmQ);
  WlzCMeshFMarQFree(nodQ);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Propagates wavefronts within a 3D conforming mesh. The
* 		wavefronts are propagated from either the mesh boundary
* 		or a number of seed nodes within the mesh.
* 		The given mesh will have modified node and element flags
* 		on return. 
* \param	mesh			Given mesh.
* \param	distances		Array for propagated distances.
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
  int		idN,
  		idS,
		nBnd;
  WlzCMeshNod3D	*nod0;
  WlzCMeshFMarQ *nodQ = NULL;
  WlzCMeshFMarQEnt *nodQEnt;
  WlzCMeshFMarElmQ *elmQ = NULL;
  WlzCMeshFMarElmQEnt *elmQEnt;
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
    /* Clear mesh element/node flags, set boundary element/node flags and
     * count number of  boundary nodes. */
    WlzCMeshClearElmFlags3D(mesh, WLZ_CMESH_NOD_FLAG_ALL);
    WlzCMeshClearNodFlags3D(mesh, WLZ_CMESH_NOD_FLAG_ALL);
    (void )WlzCMeshSetBoundElmFlags3D(mesh);
    nBnd = WlzCMeshSetBoundNodFlags3D(mesh);
    if(nBnd <= 0)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Create and then initialise the active node queue using the given seed
   * or boundary nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodQ = WlzCMeshFMarQNew(nBnd, WlzCMeshFMarHashFnNod3D)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nSeeds > 0)
    {
      idS = 0;
      while((errNum == WLZ_ERR_NONE) && (idS < nSeeds))
      {
	errNum = WlzCMeshFMarAddSeeds3D(nodQ, mesh, distances, nSeeds, seeds);
        ++idS;
      }
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
      errNum = WlzCMeshFMarAddSeeds3D(nodQ, mesh, distances, nSeeds, seeds);
      AlcFree(seeds);
    }
  }
  /* Create element queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((elmQ = WlzCMeshFMarElmQNew()) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Until the queue is empty: Pop the node with lowest priority (ie
   * lowest distance) from the queue, and process it. */
  while((errNum == WLZ_ERR_NONE) &&
	((nodQEnt = WlzCMeshFMarQPopTail(nodQ)) != NULL))
  {
    /* Find all neighbouring nodes that are neither active nor upwind.
     * For each of these neighbouring nodes, compute their distance, set
     * them to active and insert them into the queue.
     * This is done by forming a queue of all the elements that use
     * the current node from the node queue. */

    /* Get current node from node queue. */
    nod0 = (WlzCMeshNod3D *)(nodQEnt->entity);
    errNum = WlzCMeshFMarElmQInit3D(elmQ, nod0);
    if(errNum == WLZ_ERR_NONE)
    {
      /* While element list is not empty, remove element and compute all
       * node distances for it. */
      while((elmQEnt = WlzCMeshFMarElmQPopTail3D(elmQ)) != NULL)
      {
	/* Compute distances. */
	WlzCMeshFMarCompute3D(elmQEnt->nod[0].n3, elmQEnt->nod[1].n3,
			      elmQEnt->nod[2].n3, elmQEnt->nod[3].n3,
			      distances);
      }
      /* Set the current node to be upwind. */
      if(errNum == WLZ_ERR_NONE)
      {
	nod0->flags = (nod0->flags & ~(WLZ_CMESH_NOD_FLAG_ACTIVE)) |
		     WLZ_CMESH_NOD_FLAG_UPWIND;
      }
    }
  }
  /* Clear up. */
  WlzCMeshFMarQFree(nodQ);
  WlzCMeshFMarElmQFree(elmQ);
  return(errNum);
}

/*!
* \return	New constrained mesh node priority queue, NULL on error.
* \ingroup	WlzMesh
* \brief	Constructs a new constrained mesh node priority queue
* 		with room allocated for at least the given number of
* 		entries.
* \param	nEnt			Initial number of entries allocated
* 					for the queue.
*/
static WlzCMeshFMarQ *WlzCMeshFMarQNew(int nEnt,
                                       int (*hashFn)(WlzCMeshFMarQ *, void *))
{
  int		idE;
  WlzCMeshFMarQ	*queue;
  WlzCMeshFMarQEnt *ent;
  const size_t	minEntries = 1024; /* Just to avoid costly reallocation for
  				    * small queues. */

  queue = (WlzCMeshFMarQ *)AlcMalloc(sizeof(WlzCMeshFMarQ));
  if(queue)
  {
    queue->hashFn = hashFn;
    queue->max = (nEnt < minEntries)? minEntries: nEnt;
    if(((queue->entries = (WlzCMeshFMarQEnt *)
                          AlcMalloc(sizeof(WlzCMeshFMarQEnt) *
		                    queue->max)) == NULL) ||
       ((queue->buckets = (int *)
                          AlcMalloc(sizeof(int) * queue->max)) == NULL))
    {
      AlcFree(queue->entries);
      AlcFree(queue);
      queue = NULL;
    }
    else
    {
      queue->head = queue->tail = queue->last = -1;
      queue->cnt = 0;
      queue->free = 0;
      for(idE = 0; idE < queue->max; ++idE)
      {
	ent = queue->entries + idE;
#ifdef DEBUG
        ent->dbgIdx = idE;	 	/* For debug only. */
#endif
	ent->next = idE + 1;
	ent->prev = idE - 1;
        ent->hashNxt = -1;
	ent->entity = NULL;
	ent->priority = 0.0;
        *(queue->buckets + idE) = -1;
      }
      ent = queue->entries + queue->max - 1;
      ent->next = -1;
    }
  }
  return(queue);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Reallocates the queue entries and hash buckets so that
* 		there is room for at least the epecified minimum number.
* 		After the reallocatiion the hash table is recomputed.
* \param	queue			The queue to reallocate.
* \param	minEnt			Minimum number of entries.
*/
static WlzErrorNum WlzCMeshFMarQRealloc(WlzCMeshFMarQ *queue, int minEnt)
{
  int		idE;
  WlzCMeshFMarQEnt *ent;
  size_t	newMax;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const size_t	minEntryInc = 1024;

  /* Avoid frequent costly reallocation by having a minimum reallocation
   * size.  */
  newMax = queue->max + minEntryInc;
  if(minEnt > newMax)
  {
    newMax = minEnt;
  }
  /* Realllocate the data structures. */
  if(((queue->entries = (WlzCMeshFMarQEnt *)
                        AlcRealloc(queue->entries,
				   sizeof(WlzCMeshFMarQEnt) *
				   newMax)) == NULL) ||
     ((queue->buckets = (int *)
			AlcRealloc(queue->buckets,
				   sizeof(int) * newMax)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* Setup the new entries in the free list and rebuild the node
     * index hash table. */
    for(idE = queue->max; idE < newMax; ++idE)
    {
      ent = queue->entries + idE;
      ent->next = idE + 1;
      ent->prev = idE - 1;
      ent->hashNxt = -1;
      ent->entity = NULL;
      ent->priority = 0.0;
    }
    ent = queue->entries + newMax - 1;
    ent->next = queue->free;
    queue->free = queue->max;
    queue->max = newMax;
    WlzCMeshFMarQRehash(queue);
  }
  return(errNum);
}

/*!
* \ingroup	WlzMesh
* \brief	Frees the given mesh fast marching queue.
* \param	queue			The queue to free.
*/
static void	WlzCMeshFMarQFree(WlzCMeshFMarQ *queue)
{
  if(queue != NULL)
  {
    AlcFree(queue->entries);
    AlcFree(queue->buckets);
    AlcFree(queue);
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Frees the given mesh fast marching queue entry.
* 		The entry must already have been unlinked from
* 		both the queue's list and hash table.
* \param	queue			The queue.
* \param	ent			The queue entry to free.
*/
static void	WlzCMeshFMarQEntFree(WlzCMeshFMarQ *queue,
				     WlzCMeshFMarQEnt *qEnt)
{
  qEnt->next = queue->free;
  queue->free = qEnt - queue->entries;
  --(queue->cnt);
}

/*!
* \ingroup	WlzMesh
* \brief	Unlinks and frees all the mesh fast marching queue
* 		entries. This resets the queue for reuse.
* \param	queue			The queue.
*/
static void	WlzCMeshFMarQEntFreeAll(WlzCMeshFMarQ *queue)
{
  int		idE;
  WlzCMeshFMarQEnt *ent;

  queue->head = queue->tail = queue->last = -1;
  queue->cnt = 0;
  queue->free = 0;
  for(idE = 0; idE < queue->max; ++idE)
  {
    ent = queue->entries + idE;
#ifdef DEBUG
    ent->dbgIdx = idE;	 	/* For debug only. */
#endif
    ent->next = idE + 1;
    ent->prev = idE - 1;
    ent->hashNxt = -1;
    ent->entity = NULL;
    ent->priority = 0.0;
    *(queue->buckets + idE) = -1;
  }
  ent = queue->entries + queue->max - 1;
  ent->next = -1;
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
static WlzErrorNum WlzCMeshFMarQInsertNod2D(WlzCMeshFMarQ *queue,
					    WlzCMeshNod2D *nod, double dist)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEnt *ent,
  		*prevEnt;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  /* Check for node already in queue. If the entry already exists remove it
   * from the queue and the hash table, but don't put it on the free stack. */
  ent = NULL;
  idH = queue->hashFn(queue, nod);
  idE = queue->buckets[idH];
  if(idE >= 0)
  {
    /* Search through the hash table's list for an entry matching the node
     * index. If found unlink it from the queue's list and the queue's node
     * index hash table. */
    prevEnt = NULL;
    ent = queue->entries + idE;
    while(((WlzCMeshNod2D *)(ent->entity) != nod) && (ent->hashNxt > 0))
    {
      prevEnt = ent;
      ent = queue->entries + ent->hashNxt;
    }
    /* If entry found matching the node index unlink it and remove the
     * hash table entry. */
    if((WlzCMeshNod2D *)(ent->entity) == nod)
    {
      WlzCMeshFMarQUnlinkEntFromList(queue, ent);
      WlzCMeshFMarQUnlinkEntFromHash(queue, ent);
      if(prevEnt != NULL)
      {
        prevEnt->hashNxt = ent->hashNxt;
      }
      else
      {
        queue->buckets[idH] = ent->hashNxt;
      }
    }
    else
    {
      ent = NULL;
    }
  }
  if(ent == NULL)
  {
    /* Node entry not found so create a new entry. */
    ent = WlzCMeshFMarQNewEnt(queue, &errNum);
  }
  /* Insert entry into the queue and add hash table entry. */
  if(errNum == WLZ_ERR_NONE)
  {
    idE = ent - queue->entries;
    ent->priority = WlzCMeshFMarQNodPriority2D(nod, dist);
    WlzCMeshFMarQInsertEnt(queue, ent);
    ent->entity = nod;
    ent->hashNxt = queue->buckets[idH];
    queue->buckets[idH] = idE; 
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
* \param	times			Node times.
* \param	nod			Given node to insert into the queue.
*/
static WlzErrorNum WlzCMeshFMarQInsertNod3D(WlzCMeshFMarQ *queue,
				double *times, WlzCMeshNod3D *nod)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEnt *ent,
  		*prevEnt;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  /* Check for node already in queue. If the entry already exists remove it
   * from the queue and the hash table, but don't put it on the free stack. */
  ent = NULL;
  idH = queue->hashFn(queue, nod);
  idE = queue->buckets[idH];
  if(idE >= 0)
  {
    /* Search through the hash table's list for an entry matching the node
     * index. If found unlink it from the queue's list and the queue's node
     * index hash table. */
    prevEnt = NULL;
    ent = queue->entries + idE;
    while(((WlzCMeshNod3D *)(ent->entity) != nod) && (ent->hashNxt > 0))
    {
      prevEnt = ent;
      ent = queue->entries + ent->hashNxt;
    }
    /* If entry found matching the node index unlink it and remove the
     * hash table entry. */
    if((WlzCMeshNod3D *)(ent->entity) == nod)
    {
      WlzCMeshFMarQUnlinkEntFromList(queue, ent);
      WlzCMeshFMarQUnlinkEntFromHash(queue, ent);
      if(prevEnt != NULL)
      {
        prevEnt->hashNxt = ent->hashNxt;
      }
      else
      {
        queue->buckets[idH] = ent->hashNxt;
      }
    }
    else
    {
      ent = NULL;
    }
  }
  if(ent == NULL)
  {
    /* Node entry not found so create a new entry. */
    ent = WlzCMeshFMarQNewEnt(queue, &errNum);
  }
  /* Insert entry into the queue and add hash table entry. */
  if(errNum == WLZ_ERR_NONE)
  {
    idE = ent - queue->entries;
    ent->priority = WlzCMeshFMarQNodPriority3D(nod, times);
    WlzCMeshFMarQInsertEnt(queue, ent);
    ent->entity = nod;
    ent->hashNxt = queue->buckets[idH];
    queue->buckets[idH] = idE; 
  }
  return(errNum);
}

/*!
* \return	New (or recycled) entry.
* \ingroup	WlzMesh
* \brief	Gets a new priority queue which is neither in the queue
* 		nor hash table buckets.
* \param	queue			The priority queue.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshFMarQEnt *WlzCMeshFMarQNewEnt(WlzCMeshFMarQ *queue,
				WlzErrorNum *dstErr)
{
  int		idE;
  WlzCMeshFMarQEnt *ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(queue->cnt >= queue->max)
  {
    /* Reallocate entries to get more and update the free list. */
    errNum = WlzCMeshFMarQRealloc(queue, queue->cnt);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Pop entry from free list. */
    idE = queue->free;
    ent = queue->entries + idE;
    queue->free = ent->next;
    ent->next = ent->prev = ent->hashNxt = -1;
    ent->entity = NULL;
    ent->priority = 0.0;
    ++(queue->cnt);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(ent);
}

/*!
* \ingroup	WlzMesh
* \brief	Unlinks the given entry from the priority queue's list.
* 		If the head, last or tail index entry is unlinked the
* 		corresponding index is changed to the next or previous
* 		entry in the queue.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	ent			Entry to unlink.
*/
static void	WlzCMeshFMarQUnlinkEntFromList(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEnt *ent)
{
  /* If unlinked entry is the queue last entry set queue last entry to
   * next or previous entry. */
  if(queue->last == ent - queue->entries)
  {
    if(ent->next >= 0)
    {
      queue->last = ent->next;
    }
    else
    {
      queue->last = ent->prev;
    }
  }
  /* Break prev link. */
  if(ent->prev >= 0)
  {
    (queue->entries + ent->prev)->next = ent->next;
  }
  else
  {
    queue->head = ent->next;
  }
  /* Break next link. */
  if(ent->next >= 0)
  {
    (queue->entries + ent->next)->prev = ent->prev;
  }
  else
  {
    queue->tail = ent->prev;
  }
  ent->prev = ent->next = -1;
}

/*!
* \ingroup	WlzMesh
* \brief	Unlinks the given entry from the priority queue's node
* 		index hash table.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	ent			Entry to unlink.
*/
static void	WlzCMeshFMarQUnlinkEntFromHash(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEnt *gEnt)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEnt *ent,
  		*prevEnt;

  if((gEnt != NULL) && (gEnt->entity != NULL))
  {
    idH = queue->hashFn(queue, gEnt->entity);
    idE = queue->buckets[idH];
    if(idE >= 0)
    {
      prevEnt = NULL;
      ent = queue->entries + idE;
      while((ent != gEnt) && (ent->hashNxt > 0))
      {
	prevEnt = ent;
	ent = queue->entries + ent->hashNxt;
      }
      /* If entry found for node remove it. */
      if(ent == gEnt)
      {
	if(prevEnt != NULL)
	{
	  prevEnt->hashNxt = ent->hashNxt;
	}
	else
	{
	  queue->buckets[idH] = ent->hashNxt;
	}
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Recomputes the hash table of the priority queue.
* \param	queue			Given constrained mesh node priority
* 					queue.
*/
static void	WlzCMeshFMarQRehash(WlzCMeshFMarQ *queue)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEnt *ent;

  /* Clear hash table. */
  for(idE = 0; idE < queue->max; ++idE)
  {
    queue->buckets[idE] = -1;
  }
  /* Add all entries in the queue's list to the hash table. */
  idE = queue->head;
  while(idE >= 0)
  {
    ent = queue->entries + idE;
    idH = queue->hashFn(queue, ent->entity);
    ent->hashNxt = queue->buckets[idH];
    queue->buckets[idH] = idE;
    idE = ent->next;
  }
}

/*!
* \return	Node priority.
* \ingroup	WlzMesh
* \brief	Priority for node in queue. The priority is increases with
* 		distance.
* \param	nod			Node to compute priority for.
* \param	dist			Distances of the node.
*/
static double	WlzCMeshFMarQNodPriority2D(WlzCMeshNod2D *nod, double dist)
{
  return(dist);
}

/*!
* \return	Edge priority.
* \ingroup	WlzMesh
* \brief	Priority for seed element in queue. The priority increases
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
  WlzDVertex2	del;

  for(idE = 0; idE < 3; ++idE)
  {
    if(dst[elm->edu[idE].nod->idx] > DBL_MAX / 2.0)
    {
      break;
    }
  }
  if(idE < 3)
  {
    WLZ_VTX_2_SUB(del, org, elm->edu[idE].nod->pos);
    d = WLZ_VTX_2_LENGTH(del);
  }
  return(d);
}

/*!
* \return	Node priority.
* \ingroup	WlzMesh
* \brief	Priority for node in queue. The priority is increases with
* 		distance.
* \param	nod			Node to compute priority for.
* \param	dist			Array of distances indexed by the
* 					mesh node indices.
*/
static double	WlzCMeshFMarQNodPriority3D(WlzCMeshNod3D *nod, double *dist)
{
  return(*(dist + nod->idx));
}

/*!
* \ingroup	WlzMesh
* \brief	Inserts the given entry into the given queue, knowing
* 		that the entry is valid and that it's not already in
* 		the queue. The queue is kept sorted with the highest
* 		priority entry at the head and the lowest at the tail.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	ent0			Entry to insert.
*/
static void	WlzCMeshFMarQInsertEnt(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEnt *gEnt)
{
  int		idG;
  WlzCMeshFMarQEnt *ent0,
  		*ent1;

  idG = gEnt - queue->entries;
  if(queue->head < 0)
  {
    /* Queue empty. */
    gEnt->next = gEnt->prev = -1;
    queue->head = queue->tail = queue->last = idG;
  }
  else
  {
    ent1 = ent0 = queue->entries + queue->last;
    if(gEnt->priority > ent0->priority)
    {
      /* Insert entry above last, towards the head of the queue. */
      while((gEnt->priority > ent1->priority) && (ent1->prev >= 0))
      {
	ent0 = ent1;
	ent1 = queue->entries + ent0->prev;
      }
      if(ent1->prev < 0)
      {
	gEnt->next = queue->head;
	gEnt->prev = -1;
        queue->head = idG;
	ent1->prev = idG;
      }
      else
      {
	gEnt->prev = ent0->prev;
	gEnt->next = ent0 - queue->entries;
	ent0->prev = idG;
	(queue->entries + gEnt->prev)->next = idG;
      }
    }
    else
    {
      /* Insert entry below last, towards the tail of the queue. */
      while((gEnt->priority < ent1->priority) && (ent1->next >= 0))
      {
	ent0 = ent1;
	ent1 = queue->entries + ent0->next;
      }
      if(ent1->next < 0)
      {
        gEnt->next = -1;
	gEnt->prev = queue->tail;
	queue->tail = idG;
	ent1->next = idG;
      }
      else
      {
	gEnt->prev = ent0 - queue->entries;
	gEnt->next = ent0->next;
	ent0->next = idG;
        (queue->entries + gEnt->next)->prev = idG;
      }
    }
  }
  queue->last = idG;
#ifdef WLZ_CMESH_FMAR_QUEUE_DEBUG
  (void )fprintf(stderr, "queue - h % 8d t % 8d l % 8d c % 8d m % 8d f % 8d\n",
                 queue->head,
		 queue->tail,
		 queue->last,
		 queue->cnt,
		 queue->max,
		 queue->free);
  idG = queue->head;
  while(idG >= 0)
  {
    ent0 = queue->entries + idG;
    (void )fprintf(stderr, "  entry - i % 8d n % 8d p % 8d % 8d 0x%08lx %g\n",
    		   ent0->dbgIdx,
		   ent0->next,
		   ent0->prev,
		   ent0->hashNxt,
		   (unsigned long )(ent0->entity),
		   ent0->priority);
    idG = ent0->next;
  }
#endif
}

/*!
* \return	Unlinked entry.
* \ingroup	WlzMesh
* \brief	Unlinks the entry with the lowest priority from the
* 		given queue but does not free it.
* \param	queue			Given priority queue.
*/
static WlzCMeshFMarQEnt *WlzCMeshFMarQPopHead(WlzCMeshFMarQ *queue)
{
  WlzCMeshFMarQEnt *ent = NULL;

  if(queue->head >= 0)
  {
    ent = queue->entries + queue->head;
    WlzCMeshFMarQUnlinkEntFromList(queue, ent);
    WlzCMeshFMarQUnlinkEntFromHash(queue, ent);
  }
  return(ent);
}

/*!
* \return	Unlinked entry.
* \ingroup	WlzMesh
* \brief	Unlinks the entry with the highest priority from the
* 		given queue but does not free it.
* \param	queue			Given priority queue.
*/
static WlzCMeshFMarQEnt *WlzCMeshFMarQPopTail(WlzCMeshFMarQ *queue)
{
  WlzCMeshFMarQEnt *ent = NULL;

  if(queue->tail >= 0)
  {
    ent = queue->entries + queue->tail;
    WlzCMeshFMarQUnlinkEntFromList(queue, ent);
    WlzCMeshFMarQUnlinkEntFromHash(queue, ent);
  }
  return(ent);
}

/*!
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
* \param	nod2			Unknown node.
* \param	nod0			A known node directly connected to
* 					the unknown node.
* \param	nod1			A second node which shares an element
* 					with nod2 and nod0. The distance for
* 					nod1 is >= the distance for  nod0.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
* 		elm			Element using these nodes.
*/
static void	WlzCMeshFMarCompute2D(WlzCMeshNod2D *nod2, WlzCMeshNod2D *nod0,
				WlzCMeshNod2D *nod1, double *distances,
				WlzCMeshElm2D *elm)
{
  int		idN,
  		flg;
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
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element.
* 		This function just classifies the problem according to the
* 		number of known nodes and passes then calls the appropriate
* 		function.
* \param	nod0			Current node (first).
* \param	nod1			Second node of the element.
* \param	nod2			Third node of the element.
* \param	nod3			Fourth node of the element.
* \param	nod			Array of the four nodes used by
* 					this element.
* \param	distances		Array of distanes indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static void	WlzCMeshFMarCompute3D(WlzCMeshNod3D *nod0,
				      WlzCMeshNod3D *nod1,
				      WlzCMeshNod3D *nod2,
				      WlzCMeshNod3D *nod3,
				      double *distances)
{
  int		kwnMsk = 0;

  kwnMsk = (((nod0->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 0) |
           (((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 1) |
           (((nod2->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 2) |
           (((nod3->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) << 3);
  switch(kwnMsk)
  {
    case  1: /* 0001 */
      WlzCMeshFMarCompute3D1(nod0, nod1, nod2, nod3, distances);
      break;
    case  2: /* 0010 */
      WlzCMeshFMarCompute3D1(nod1, nod2, nod3, nod0, distances);
      break;
    case  3: /* 0011 */
      WlzCMeshFMarCompute3D2(nod0, nod1, nod2, nod3, distances);
      break;
    case  4: /* 0100 */
      WlzCMeshFMarCompute3D1(nod2, nod3, nod0, nod1, distances);
      break;
    case  5: /* 0101 */
      WlzCMeshFMarCompute3D2(nod0, nod2, nod1, nod3, distances);
      break;
    case  6: /* 0110 */
      WlzCMeshFMarCompute3D2(nod1, nod2, nod0, nod3, distances);
      break;
    case  7: /* 0111 */
      WlzCMeshFMarCompute3D3(nod0, nod1, nod2, nod3, distances);
      break;
    case  8: /* 1000 */
      WlzCMeshFMarCompute3D1(nod3, nod0, nod1, nod2, distances);
      break;
    case  9: /* 1001 */
      WlzCMeshFMarCompute3D2(nod0, nod3, nod1, nod2, distances);
      break;
    case 10: /* 1010 */
      WlzCMeshFMarCompute3D2(nod1, nod3, nod0, nod2, distances);
      break;
    case 11: /* 1011 */
      WlzCMeshFMarCompute3D3(nod0, nod1, nod3, nod2, distances);
      break;
    case 12: /* 1100 */
      WlzCMeshFMarCompute3D2(nod2, nod3, nod0, nod1, distances);
      break;
    case 13: /* 1101 */
      WlzCMeshFMarCompute3D3(nod0, nod2, nod3, nod1, distances);
      break;
    case 14: /* 1110 */
      WlzCMeshFMarCompute3D3(nod1, nod2, nod3, nod0, distances);
      break;
    default:
      break;
  }
}

/*!
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
static void	WlzCMeshFMarCompute3D1(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  int		idx;
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
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the unknown nodes of
* 		the given element.
* 		This function is given just two known nodes and computes the
* 		times for the other nodes by assuming that the propagation
* 		is along the edges of the element.
* \param	nod0			First known (current) node.
* \param	nod1			Second known node.
* \param	nod2			First unknown node.
* \param	nod3			Second unknown node.
* \param	distances		Array of distances indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static void	WlzCMeshFMarCompute3D2(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  /* TODO restrict to face of element. */
  /* HACK WlzCMeshFMarCompute3D1(nod0, nod1, nod2, nod3, s0, times); */
}

/*!
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
*		With the normal of the advancing front at \f%n_3\f$
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
static void	WlzCMeshFMarCompute3D3(WlzCMeshNod3D *nod0,
                                       WlzCMeshNod3D *nod1,
                                       WlzCMeshNod3D *nod2,
                                       WlzCMeshNod3D *nod3,
				       double *distances)
{
  int		id0,
  		id1,
  		onEdge = 1,
		hit = 0;
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
     * from node 3 towards the face formed by nodes 0, 1 and 2. */
    WLZ_VTX_3_ADD3(t0, nod[0]->pos, nod[1]->pos, nod[2]->pos);
    WLZ_VTX_3_SCALE(t0, t0, 1.0 / 3.0);
    WLZ_VTX_3_SUB(t1, t0, nod[3]->pos);
    a = WLZ_VTX_3_DOT(n0, t1);
    if(a < 0)
    {
      n0 = n1;
    }
    hit = WlzGeomRayTriangleIntersect3D(n0,
					nod[0]->pos, nod[1]->pos, nod[2]->pos,
					nod[3]->pos, NULL, NULL);
    if(hit)
    {
      /* Normal from node 3 passes through the triangle formed by the three
       * known nodes: 0, 1 and 2. Calculate the distance from the plane of
       * propagation (which passes through node 0) to node 3.
       * This is just: \frac{1}{s_3}(p_0 - p_3) . n + p3. */
      WLZ_VTX_3_SUB(l[3], nod[3]->pos, nod[0]->pos);
      WLZ_VTX_3_SUB(t0, n0, nod[3]->pos);
      d = WLZ_VTX_3_LENGTH(t0);
      if(d > 0.0)
      {
	d = *(distances + nod[0]->idx) + d;
	if(*(distances + nod[3]->idx) > d)
	{
	  *(distances + nod[3]->idx) = d;
	  onEdge = 0;
	}
      }
    }
  }
  if(onEdge != 0)
  {
    /* TODO Do something more sophisticated here if needed. 
    for(id0 = 0; id0 < 3; ++id0)
    {
      WLZ_VTX_3_SUB(t0, nod[id0]->pos, nod[3]->pos);
      d = WLZ_VTX_3_LENGTH(t0);
      if(d > 0.0)
      {
	d += *(distances + nod[0]->idx);
	if(*(distances + nod[3]->idx) > d)
	{
	  *(distances + nod[3]->idx) = d;
	}
      }
    }
    */
  }
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
static WlzErrorNum WlzCMeshFMarAddSeeds2D(WlzCMeshFMarQ *nodQ,
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
  WlzCMeshFMarQ	*sElmQ = NULL;
  WlzUByte	*eFlgs = NULL;
  double	*sDists = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idS = 0;
  if(((sDists = (double *)
                AlcMalloc(sizeof(double) * mesh->res.nod.maxEnt)) == NULL) ||
     ((eFlgs = (WlzUByte *)
                AlcMalloc(sizeof(WlzUByte) * mesh->res.elm.maxEnt)) == NULL) ||
     ((sElmQ = WlzCMeshFMarQNew(edgQMin, WlzCMeshFMarHashFnElm2D)) == NULL))
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
      WlzCMeshFMarQEntFreeAll(sElmQ);
      ++idS;
    }
  }
  AlcFree(sDists);
  AlcFree(eFlgs);
  WlzCMeshFMarQFree(sElmQ);
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
	  nod0->flags |= WLZ_CMESH_NOD_FLAG_ACTIVE;
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
* \brief	Adds a single seed to the 
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
static WlzErrorNum WlzCMeshFMarAddSeed2D(WlzCMeshFMarQ  *sElmQ,
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
  idM = WlzCMeshElmEnclosingPos2D(mesh, -1, seed.vtX, seed.vtY, NULL);
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
        ((sElmQEnt = WlzCMeshFMarQPopTail(sElmQ)) != NULL))
  {
    elm0 = (WlzCMeshElm2D *)(sElmQEnt->entity);
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
    WlzCMeshFMarQEntFree(sElmQ, sElmQEnt);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts an element into the seed element queue.
* \param	queue			The seed element queue.
* \param	elm			Element use to insert into the queue.
* \param	dst			Distances for the seed.
* \param	org			Seed for Euclidean distances.
*/
static WlzErrorNum WlzCMeshFMarSElmQInsert2D(WlzCMeshFMarQ *queue,
					     WlzCMeshElm2D *elm,
					     double *dst,
					     WlzDVertex2 org)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEnt *ent;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  if(elm != NULL)
  {
    idH = queue->hashFn(queue, elm);
    ent = WlzCMeshFMarQNewEnt(queue, &errNum);
    /* Insert entry into the queue and add hash table entry. */
    if(errNum == WLZ_ERR_NONE)
    {
      idE = ent - queue->entries;
      ent->entity = elm;
      ent->priority = WlzCMeshFMarQSElmPriority2D(elm, dst, org);
      WlzCMeshFMarQInsertEnt(queue, ent);
      ent->hashNxt = queue->buckets[idH];
      queue->buckets[idH] = idE; 
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the queue as part of the initial
* 		seeding.
* \param	queue			The queue.
* \param	mesh			The constrained mesh.
* \param	distances		Array of distances.
* \param	nSeeds			Number of seeds.
* \param	seeds			Seed positions.
*/
static WlzErrorNum WlzCMeshFMarAddSeeds3D(WlzCMeshFMarQ *queue,
				WlzCMesh3D *mesh, double *distances,
				int nSeeds, WlzDVertex3 *seeds)
{
#ifdef HACK_UNDEF
  int		cls = 0,
  		idE,
		idF,
		idM;
  WlzCMeshEdgU3D *edu0,
  		*edu1;
  WlzCMeshNod3D	*nod0,
  		*nod1;
  WlzCMeshFace	*fce0;
  WlzCMeshElm3D *elm0,
  		*elm1;
  WlzCMeshNod3D	*nodes0[4],
  		*nodes1[4];
#endif
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef HACK_UNDEF
  idM = WlzCMeshElmEnclosingPos3D(mesh, -1,
                                  seedPos.vtX, seedPos.vtY, seedPos.vtZ,
				  NULL);
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
    /* Is seed coincident with a mesh node? If so set cls == 1 and nod0 to
     * node. */
    WlzCMeshElmGetNodes3D(elm0, nodes0 + 0, nodes0 + 1, nodes0 + 2,
                          nodes0 + 3);
    for(idE = 0; idE < 3; ++idE)
    {
      nod0 = nodes0[idE];
      if(WlzGeomVtxEqual3D(nod0->pos, seedPos, WLZ_MESH_TOLERANCE_SQ))
      {
	cls = 1;
        break;
      }
    }
    if(cls == 0)
    {
      /* Does seed lie on an edge of the element? If so set cls == 2 and nod0,
       * nod1 to nodes at either end of the edge. */
      for(idE = 0; idE < 3; ++idE)
      {
	idF = (idE + 1) % 3;
	nod1 = nodes0[idF];
	if(WlzGeomVtxOnLineSegment3D(seedPos,
	                              (nod0 = nodes0[idE])->pos, nod1->pos,
			              WLZ_MESH_TOLERANCE))
	{
	  cls = 2;
	  break;
	}
	else if(WlzGeomVtxOnLineSegment3D(seedPos,
	                              (nod0 = nodes0[3])->pos, nod1->pos,
			              WLZ_MESH_TOLERANCE))
	{
	  cls = 2;
	  break;
	}
      }
    }
    if(cls == 0)
    {
      /* Does seed lie on a face of the element? If so set cls == 3 and fce0
       * to the face. */
      for(idF = 0; idF < 4; ++idF)
      {
	fce0 = &(elm0->face[idF]);
        if(WlzGeomVxInTriangle3D(seedPos, fce0->edu[0].nod->pos,
	                                  fce0->edu[1].nod->pos,
					  fce0->edu[2].nod->pos) >= 0)
	{
	  cls = 3;
	  break;
	}
      }
    }
    switch(cls)
    {
      case 0:                           /* Seed is contained within element. */
        /* Add each of the elements nodes to the queue. */
	idE = 0;
	do
        {
	  errNum = WlzCMeshFMarInsertSeed3D(queue, nodes0[idE],
					    times, speeds, seedPos);
	} while((errNum == WLZ_ERR_NONE) && (++idE < 4));
        break;
      case 1:                          /* Seed is coincident with mesh node. */
        /* Add all nodes of all elements that use this node. */
	edu0 = edu1 = nod0->edu;
	do
	{
	  edu1 = edu1->nnxt;
	  elm1 = edu1->face->elm;
	  WlzCMeshElmGetNodes3D(elm1, nodes1 + 0, nodes1 + 1, nodes1 + 2,
				nodes1 + 3);
          idE = 0;
	  do
	  {
	    errNum = WlzCMeshFMarInsertSeed3D(queue,
					      nodes1[idE],
					      times, speeds, seedPos);
	  } while((errNum == WLZ_ERR_NONE) && (++idE < 4));
	} while((errNum == WLZ_ERR_NONE) && (edu1 != edu0));
        break;
      case 2:                          /* Seed is on an edge of the element. */
        /* Add all nodes of all elements that use this edge. */
	edu1 = edu0;
	nod0 = edu0->next->nod;
	do
	{
	  edu1 = edu1->nnxt;
	  if(edu1->next->nod == nod0)
	  {
	    elm1 = edu1->face->elm;
	    WlzCMeshElmGetNodes3D(elm1, nodes1 + 0, nodes1 + 1, nodes1 + 2,
				  nodes1 + 3);
	    idE = 0;
	    do
	    {
	      errNum = WlzCMeshFMarInsertSeed3D(queue,
						nodes1[idE],
						times, speeds, seedPos);
	    } while((errNum == WLZ_ERR_NONE) && (++idE < 4));
	  }
	} while((errNum == WLZ_ERR_NONE) && (edu1 != edu0));
        break;
      case 3:                           /* Seed is on a face of the element. */
        /* Add all nodes of those elements that share this face. */
	idE = 0;
	do
        {
	  errNum = WlzCMeshFMarInsertSeed3D(queue, nodes0[idE],
					    times, speeds, seedPos);
	} while((errNum == WLZ_ERR_NONE) && (++idE < 4));
	if((errNum == WLZ_ERR_NONE) &&
	   (fce0->opp != NULL) && (fce0->opp != fce0))
	{
	  elm1 = fce0->opp->elm;
	  idE = 0;
	  do
	  {
	    errNum = WlzCMeshFMarInsertSeed3D(queue, nodes0[idE],
					      times, speeds, seedPos);
	  } while((errNum == WLZ_ERR_NONE) && (++idE < 4));
	}
        break;

      default:
        break;
    }
  }
#endif /* HACK_TODO */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the queue as part of the initial
* 		seeding.
* \param	queue			The queue.
* \param	nod			Node to be added to the queue.
* \param	times			Array of wavefront arrival times
* 					indexed by the mesh node indices,
* 					which will be set on return.
* \param	speeds			Array of wavefront propagation
* 					speeds indexed by the mesh node
* 					indices, may be NULL in which
* 					case a constant speed of 1.0 is used.
* 					Speeds must all be > zero.
* \param	seedPos			Seed position.
*/
static WlzErrorNum WlzCMeshFMarInsertSeed3D(WlzCMeshFMarQ *queue,
				WlzCMeshNod3D *nod,
				double *times, double *speeds,
				WlzDVertex3 seedPos)
{
  double	newT;
  WlzDVertex3	dsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  WLZ_VTX_3_SUB(dsp, seedPos, nod->pos);
  newT = WLZ_VTX_3_LENGTH(dsp);
  if(speeds != NULL)
  {
    newT /= *(speeds + nod->idx);
  }
  if(newT < *(times + nod->idx))
  {
    *(times + nod->idx) = newT;
    nod->flags |= WLZ_CMESH_NOD_FLAG_ACTIVE | WLZ_CMESH_NOD_FLAG_KNOWN;
    errNum = WlzCMeshFMarQInsertNod3D(queue, times, nod);
  }
  return(errNum);
}

/*!
* \return	New element queue with no elements allocated.
* \ingroup	WlzMesh
* \brief	Allocates a new element queue, WlzCMeshFMarElmQFree()
* 		should be used to free the queue. No queue entries
* 		are allocated by this function.
*/
static WlzCMeshFMarElmQ *WlzCMeshFMarElmQNew(void)
{
  WlzCMeshFMarElmQ *queue;

  queue = (WlzCMeshFMarElmQ *)AlcCalloc(1, sizeof(WlzCMeshFMarElmQ));
  return(queue);
}

/*!
* \ingroup	WlzMesh
* \brief	Frees an element queue created by WlzCMeshFMarElmQNew().
* \param	queue			Queue to free.
*/
static void		WlzCMeshFMarElmQFree(WlzCMeshFMarElmQ *queue)
{
  if(queue)
  {
    AlcFree(queue->entries);
    AlcFree(queue);
  }
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
static WlzErrorNum WlzCMeshFMarElmQInit2D(WlzCMeshFMarElmQ *queue,
					  WlzCMeshNod2D *nod)
{
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  WlzCMeshFMarElmQEnt *ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	queueEntInc = 1024;

  queue->nEnt = 0;
  queue->nod.n2 = nod;
  edu0 = edu1 = nod->edu;
  do
  {
    /* Add entry to the queue for each edge directed from the node. */
    if(queue->nEnt >= queue->maxEnt)
    {
      queue->maxEnt += queueEntInc;
      if((queue->entries = (WlzCMeshFMarElmQEnt *)
                           AlcRealloc(queue->entries,
				      sizeof(WlzCMeshFMarElmQEnt) *
				      queue->maxEnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      ent = queue->entries + queue->nEnt++;
      ent->elm.e2 = edu1->elm;
      /* Node pointers are set in WlzCMeshFMarElmQSqueeze2D() after
       * removing redundant entries. */
      edu1 = edu1->nnxt;
    }
  } while((errNum == WLZ_ERR_NONE) && (edu1 != edu0));
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshFMarElmQSqueeze2D(queue);
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
static WlzErrorNum WlzCMeshFMarElmQInit3D(WlzCMeshFMarElmQ *queue,
					  WlzCMeshNod3D *nod)
{
  WlzCMeshEdgU3D *edu0,
  		*edu1;
  WlzCMeshFMarElmQEnt *ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	queueEntInc = 1024;

  queue->nEnt = 0;
  queue->nod.n3 = nod;
  edu0 = edu1 = nod->edu;
  do
  {
    /* Add entry to the queue for each edge directed from the node. */
    if(queue->nEnt >= queue->maxEnt)
    {
      queue->maxEnt += queueEntInc;
      if((queue->entries = (WlzCMeshFMarElmQEnt *)
                           AlcRealloc(queue->entries,
				      sizeof(WlzCMeshFMarElmQEnt) *
				      queue->maxEnt)) == NULL)
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      ent = queue->entries + queue->nEnt++;
      ent->elm.e3 = edu1->face->elm;
      /* Node pointers are set in WlzCMeshFMarElmQSqueeze3D() after
       * removing redundant entries. */
      edu1 = edu1->nnxt;
    }
  } while((errNum == WLZ_ERR_NONE) && (edu1 != edu0));
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshFMarElmQSqueeze3D(queue);
  }
  return(errNum);
}

/*!
* \ingroup	WlzMesh
* \brief	Squeezes out queue entries which have duplicate elements
* 		or have all four nodes known.
* \param	queue			Element queue.
*/
static void	WlzCMeshFMarElmQSqueeze2D(WlzCMeshFMarElmQ *queue)
{
  int		idx0,
  		idx1,
		idx2;
  WlzCMeshNod2D *nod[3];
  WlzCMeshFMarElmQEnt *ent0,
  		*ent1;

  /* Squeeze out the duplicate element index entries. */
  qsort(queue->entries, queue->nEnt, sizeof(WlzCMeshFMarElmQEnt),
        WlzCMeshFMarElmQElmIdxCmp2D);
  idx0 = 0;
  idx1 = 1;
  ent0 = queue->entries + 0;
  ent1 = queue->entries + 1;
  while(idx1 < queue->nEnt)
  {
    if(ent0->elm.e2->idx != ent1->elm.e2->idx)
    {
      ++idx0;
      ++ent0;
      *ent0 = *ent1;
    }
    ++idx1;
    ++ent1;
  }
  queue->nEnt = idx0 + 1;
  /* Set node pointers and compute priority for the elements. */
  for(idx0 = 0; idx0 < queue->nEnt; ++idx0)
  {
    ent0 = queue->entries + idx0;
    WlzCMeshElmGetNodes2D(ent0->elm.e2, nod + 0, nod + 1, nod + 2);
    /* Make sure that the current node is the first one. */
    idx1 = 0;
    if(nod[1] == queue->nod.n2)
    {
      idx1 = 1;
    }
    else if(nod[2] == queue->nod.n2)
    {
      idx1 = 2;
    }
    for(idx2 = 0; idx2 < 3; ++idx2)
    {
       ent0->nod[idx2].n2 = nod[(idx1 + idx2) % 3];
    }
  }
  WlzCMeshFMarElmQSort2D(queue);
}

/*!
* \ingroup	WlzMesh
* \brief	Squeezes out queue entries which have duplicate elements
* 		or have all four nodes known.
* \param	queue			Element queue.
*/
static void	WlzCMeshFMarElmQSqueeze3D(WlzCMeshFMarElmQ *queue)
{
  int		idx0,
  		idx1,
		idx2;
  WlzCMeshNod3D *nod[4];
  WlzCMeshFMarElmQEnt *ent0,
  		*ent1;

  /* Squeeze out the duplicate element index entries. */
  qsort(queue->entries, queue->nEnt, sizeof(WlzCMeshFMarElmQEnt),
        WlzCMeshFMarElmQElmIdxCmp3D);
  idx0 = 0;
  idx1 = 1;
  ent0 = queue->entries + 0;
  ent1 = queue->entries + 1;
  while(idx1 < queue->nEnt)
  {
    if(ent0->elm.e3->idx != ent1->elm.e3->idx)
    {
      ++idx0;
      ++ent0;
      *ent0 = *ent1;
    }
    ++idx1;
    ++ent1;
  }
  queue->nEnt = idx0 + 1;
  /* Set node pointers and compute priority for the elements. */
  for(idx0 = 0; idx0 < queue->nEnt; ++idx0)
  {
    ent0 = queue->entries + idx0;
    WlzCMeshElmGetNodes3D(ent0->elm.e3, nod + 0, nod + 1, nod + 2, nod + 3);
    /* Make sure that the current node is the first one. */
    idx1 = 0;
    if(nod[1] == queue->nod.n3)
    {
      idx1 = 1;
    }
    else if(nod[2] == queue->nod.n3)
    {
      idx1 = 2;
    }
    else if(nod[3] == queue->nod.n3)
    {
      idx1 = 3;
    }
    for(idx2 = 0; idx2 < 4; ++idx2)
    {
       ent0->nod[idx2].n3 = nod[(idx1 + idx2) % 4];
    }
  }
  WlzCMeshFMarElmQSort3D(queue);
  /* Squeeze out unwanted entries, ie those with a priority > 7. */
  idx0 =  queue->nEnt - 1;
  while(idx0 >= 0)
  {
    ent0 = queue->entries + idx0;
    if(ent0->priority < 7)
    {
      break;
    }
    --idx0;
  }
  queue->nEnt = idx0 + 1;
}

/*!
* \return	List element index entry which will be posative for all
* 		valid elements), or a negative value if the list is empty.
* \ingroup	WlzMesh
* \brief	Removes the entry at the tail of the list and retuurns it.
* \param	queue			Element queue to get entry from.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshFMarElmQEnt *WlzCMeshFMarElmQPopTail2D(WlzCMeshFMarElmQ *queue)
{
  WlzCMeshFMarElmQEnt *ent = NULL;

  WlzCMeshFMarElmQSort2D(queue);
  if(queue->nEnt > 0)
  {
    ent = queue->entries + --(queue->nEnt);
  }
  return(ent);
}

/*!
* \return	List element index entry which will be posative for all
* 		valid elements), or a negative value if the list is empty.
* \ingroup	WlzMesh
* \brief	Removes the entry at the tail of the list and retuurns it.
* \param	queue			Element queue to get entry from.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshFMarElmQEnt *WlzCMeshFMarElmQPopTail3D(WlzCMeshFMarElmQ *queue)
{
  WlzCMeshFMarElmQEnt *ent = NULL;

  WlzCMeshFMarElmQSort3D(queue);
  if(queue->nEnt > 0)
  {
    ent = queue->entries + --(queue->nEnt);
  }
  return(ent);
}

/*!
* \ingroup	WlzMesh
* \brief	Sorts the mesh element queue by priority so that highest
* 		priority elements are last ani  the list.
* \param	queue			The mesh element queue.
*/
static void	WlzCMeshFMarElmQSort2D(WlzCMeshFMarElmQ *queue)
{
  int		idx;
  WlzCMeshFMarElmQEnt *ent = NULL;

  for(idx = 0; idx < queue->nEnt; ++idx)
  {
    ent = queue->entries + idx;
    WlzCMeshFMarElmQCalcPriority2D(ent, queue->nod.n2);
  }
  qsort(queue->entries, queue->nEnt, sizeof(WlzCMeshFMarElmQEnt),
        WlzCMeshFMarElmQPriorityCmp);
  idx =  queue->nEnt - 1;
  while(idx >= 0)
  {
    ent = queue->entries + idx;
    if(ent->priority < 5)
    {
      break;
    }
    --idx;
  }
  queue->nEnt = idx + 1;
}

/*!
* \ingroup	WlzMesh
* \brief	Sorts the mesh element queue by priority so that highest
* 		priority elements are last ani  the list.
* \param	queue			The mesh element queue.
*/
static void	WlzCMeshFMarElmQSort3D(WlzCMeshFMarElmQ *queue)
{
  int		idx;
  WlzCMeshFMarElmQEnt *ent = NULL;

  for(idx = 0; idx < queue->nEnt; ++idx)
  {
    ent = queue->entries + idx;
    WlzCMeshFMarElmQCalcPriority3D(ent, queue->nod.n3);
  }
  qsort(queue->entries, queue->nEnt, sizeof(WlzCMeshFMarElmQEnt),
        WlzCMeshFMarElmQPriorityCmp);
  idx =  queue->nEnt - 1;
  while(idx >= 0)
  {
    ent = queue->entries + idx;
    if(ent->priority < 7)
    {
      break;
    }
    --idx;
  }
  queue->nEnt = idx + 1;
}
  
/*!
* \ingroup	W;zMesh
* \brief	Computes the priority value of an element queue entry.
*		The priority is the simple sum of the priority of the
*		nodes (2 for an upwind node or the current node,
*		1 for any other known node and zero for an unknown
*		(downwind) node.
* \param	ent			The element queue entry.
* \param	cNod			The current node around wich the
* 					elements in the queue are clustered.
*/
static void	WlzCMeshFMarElmQCalcPriority2D(WlzCMeshFMarElmQEnt *ent,
				WlzCMeshNod2D *cNod)
{
  int		idx;
  WlzCMeshNod2D *nod;

  ent->priority = 0;
  for(idx = 0; idx < 3; ++idx)
  {
    nod = ent->nod[idx].n2;
    ent->priority += ((nod->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) +
                     ((nod->idx == cNod->idx) ||
		      ((nod->flags & WLZ_CMESH_NOD_FLAG_UPWIND) != 0));
  }
}

/*!
* \ingroup	W;zMesh
* \brief	Computes the priority value of an element queue entry.
*		The priority is the simple sum of the priority of the
*		nodes (2 for an upwind node or the current node,
*		1 for any other known node and zero for an unknown
*		(downwind) node.
* \param	ent			The element queue entry.
* \param	cNod			The current node around wich the
* 					elements in the queue are clustered.
*/
static void	WlzCMeshFMarElmQCalcPriority3D(WlzCMeshFMarElmQEnt *ent,
				WlzCMeshNod3D *cNod)
{
  int		idx;
  WlzCMeshNod3D *nod;

  ent->priority = 0;
  for(idx = 0; idx < 4; ++idx)
  {
    nod = ent->nod[idx].n3;
    ent->priority += ((nod->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0) +
                     ((nod->idx == cNod->idx) ||
		      ((nod->flags & WLZ_CMESH_NOD_FLAG_UPWIND) != 0));
  }
}

/*!
* \return	Comparison value for qsort().
* \ingroup	WlzMesh
* \brief	Compares the priority the two given element queue entries.
* \param	p0			Pointer for first entry.
* \param	p1			Pointer for second entry.
*/
static int	WlzCMeshFMarElmQPriorityCmp(const void *p0, const void *p1)
{
  int		cmp;
  WlzCMeshFMarElmQEnt *ent0,
  		*ent1;

  ent0 = (WlzCMeshFMarElmQEnt *)p0;
  ent1 = (WlzCMeshFMarElmQEnt *)p1;
  cmp = ent0->priority - ent1->priority;
  return(cmp);
}

/*!
* \return	Comparison value for qsort().
* \ingroup	WlzMesh
* \brief	Compares the element index the two given element queue entries.
* \param	p0			Pointer for first entry.
* \param	p1			Pointer for second entry.
*/
static int	WlzCMeshFMarElmQElmIdxCmp2D(const void *p0, const void *p1)
{
  int		cmp;
  WlzCMeshFMarElmQEnt *ent0,
  		*ent1;

  ent0 = (WlzCMeshFMarElmQEnt *)p0;
  ent1 = (WlzCMeshFMarElmQEnt *)p1;
  cmp = ent0->elm.e2->idx - ent1->elm.e2->idx;
  return(cmp);
}
/*!
* \return	Comparison value for qsort().
* \ingroup	WlzMesh
* \brief	Compares the element index the two given element queue entries.
* \param	p0			Pointer for first entry.
* \param	p1			Pointer for second entry.
*/
static int	WlzCMeshFMarElmQElmIdxCmp3D(const void *p0, const void *p1)
{
  int		cmp;
  WlzCMeshFMarElmQEnt *ent0,
  		*ent1;

  ent0 = (WlzCMeshFMarElmQEnt *)p0;
  ent1 = (WlzCMeshFMarElmQEnt *)p1;
  cmp = ent0->elm.e3->idx - ent1->elm.e3->idx;
  return(cmp);
}

/*!
* \return	Hash key.
* \ingroup	WlzMesh
* \brief	Simple hash function for integer values. The hash value
* 		returned is in the range [0-maxVal], but maxVal must be
* 		less than the largest of the generator primes (99999989).
* 		All this function has to do is map the given integers
* 		fairly uniformly over the integer range 0 - maxVal.
* \todo		TODO Check the distribution is fairly uniform, without
* 		to many spikes.
* \param	value			Given integer value.
* \param	maxVal			Maximum hash value.
*/
static int	WlzCMeshFMarQHashFn(int value, int maxVal)
{
  int		key;
  const long long p0 = 15486157, 	   /* 4 different increasing primes. */
  		p1 = 32453039,
		p2 = 46048241,
		p3 = 99999989;

  key = (((((long long)value * p0) ^ p1) + p2) % p3) % maxVal;
  return(key);
}

/*!
* \return	Hash table index.
* \ingroup	WlzMesh
* \brief	Computes a hash table index for a 2D CMesh node.
* 		The hash index is computed using the nodes index.
* \param	queue			Given queue.
* \param	ent			Entity which is a 2D node.
*/
static int	WlzCMeshFMarHashFnNod2D(WlzCMeshFMarQ *queue, void *ent)
{
  int		idx;
  WlzCMeshNod2D	*nod;

  nod = (WlzCMeshNod2D *)ent;
  idx = WlzCMeshFMarQHashFn(nod->idx, queue->max);
  return(idx);
}

/*!
* \return	Hash table index.
* \ingroup	WlzMesh
* \brief	Computes a hash table index for a 3D CMesh node.
* 		The hash index is computed using the nodes index.
* \param	queue			Given queue.
* \param	ent			Entity which is a 3D node.
*/
static int	WlzCMeshFMarHashFnNod3D(WlzCMeshFMarQ *queue, void *ent)
{
  int		idx;
  WlzCMeshNod3D	*nod;

  nod = (WlzCMeshNod3D *)ent;
  idx = WlzCMeshFMarQHashFn(nod->idx, queue->max);
  return(idx);
}

/*!
* \return	Hash table index.
* \ingroup	WlzMesh
* \brief	Computes a hash table index for a 2D CMesh  use.
* 		The hash table index is computed using the index of the
* 		element.
* \param	queue			Given queue.
* \param	ent			Entity which is a 2D element.
*/
static int	WlzCMeshFMarHashFnElm2D(WlzCMeshFMarQ *queue, void *ent)
{
  int		idx;
  WlzCMeshElm2D *elm;

  elm = (WlzCMeshElm2D *)ent;
  idx = WlzCMeshFMarQHashFn(elm->idx, queue->max);
  return(idx);
}
