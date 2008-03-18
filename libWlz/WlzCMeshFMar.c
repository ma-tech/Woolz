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
#include <Wlz.h>

typedef struct _WlzCMeshFMarQEntry
{
  int			next;		/* Index of next towards tail,
  					   -ve iff at tail. */
  int			prev;		/* Index of next towards head,
                                           -ve iff at head. */
  int			hashNxt;	/* Index of next in hash bucket. */
  int			nodIdx;		/* Index of mesh node. */
  double		priority;	/* Entry priority highest priority
  					   at the head of the queue. */
} WlzCMeshFMarQEntry;

typedef struct _WlzCMeshFMarQ
{
  int			head;           /* Index of the queue head, with
                                           the highest priority entry at
					   the head. */
  int			tail;           /* Index of the queue tail, with
  					   the lowest priority entry at
					   the tail. */
  int			last;		/* Index of last entry inserted. */
  int			cnt;		/* Number of entries in queue. */
  int			max;		/* Number of entries allocated. */
  int 			free;		/* Index of first free entry, rest
                                           via next index. */
  WlzCMeshFMarQEntry    *entries;	/* Array of allocated entries. */
  int			*buckets;	/* Indices for hash buckets, used
  					   fast for access to entry by node
					   index. */
} WlzCMeshFMarQ;

static WlzCMeshFMarQ 		*WlzCMeshFMarQNew(
				  int nEntries);
static WlzErrorNum 		WlzCMeshFMarQRealloc(
				  WlzCMeshFMarQ *queue,
				  int minEnt);
static void			WlzCMeshFMarQFree(
				  WlzCMeshFMarQ *queue);
static void			WlzCMeshFMarQEntryFree(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEntry *qEnt);
static WlzErrorNum 		WlzCMeshFMarQInsert2D(
				  WlzCMeshFMarQ *queue,
				  double *times,
				  WlzCMeshNod2D *nod);
static WlzCMeshFMarQEntry 	*WlzCMeshFMarQNewEntry(
				  WlzCMeshFMarQ *queue,
				  WlzErrorNum *dstErr);
static void			WlzCMeshFMarQEntryUnlink(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEntry *ent);
static void			WlzCMeshFMarQRehash(
				  WlzCMeshFMarQ *queue);
static int			WlzCMeshFMarQHashFn(
				  int value,
				  int maxVal);
static void			WlzCMeshFMarQInsertEntry(
				  WlzCMeshFMarQ *queue,
				  WlzCMeshFMarQEntry *gEnt);
static WlzCMeshFMarQEntry 	*WlzCMeshFMarQPopTail(
				  WlzCMeshFMarQ *queue);
static void		 	WlzCMeshFMarCompute2D(
				  WlzCMeshNod2D *nod0,
				  WlzCMeshNod2D *nod1,
				  double s2,
				  double *times);
static double		   	WlzCMeshFMarQPriority(
				  WlzCMeshNod2D *nod,
				  double *times);

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Propagates wavefronts within a 2D conforming mesh. The
* 		wavefronts are propagated from either the mesh boundary
* 		or a number of seed positions within the mesh.
* 		The given mesh will have modified node and element flags
* 		on return. 
* \param	mesh			Given mesh.
* \param	times			Array of wavefront arrival times
* 					indexed by the mesh node indices,
* 					which will be set on return.
* \param	speeds			Array of wavefront propagation
* 					speeds indexed by the mesh node
* 					indices, may be NULL in which
* 					case a constant speed of 1.0 is used.
* 					Speeds must all be > zero.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seedPos			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
* 					any seeds are not within the
* 					mesh.
*/
WlzErrorNum	WlzCMeshFMarNodes2D(WlzCMesh2D *mesh,
				double *times, double *speeds,
				int nSeeds, WlzDVertex2 *seedPos)
{
  int		idE,
  		idN,
  		idS,
		nBnd;
  WlzDVertex2	dsp;
  WlzCMeshElm2D	*elm;
  WlzCMeshNod2D	*nod0,
  		*nod1;
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  WlzCMeshNod2D	*nodes[3];
  WlzCMeshFMarQ *queue = NULL;
  WlzCMeshFMarQEntry *qEnt;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((times == NULL) || ((nSeeds > 0) && (seedPos == NULL)))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Set time for all mesh nodes to maximum value. */
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      *(times + idN) = DBL_MAX;
    }
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
    if((queue = WlzCMeshFMarQNew(nBnd)) == NULL)
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
        idE = WlzCMeshElmEnclosingPos2D(mesh, -1,
	                                (seedPos + idS)->vtX,
					(seedPos + idS)->vtY);
	if(idE < 0)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
	else
	{
	  elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	  if(elm->idx > 0)
	  {
	    WlzCMeshElmGetNodes2D(elm, nodes + 0, nodes + 1, nodes + 2);
	    for(idN = 0; idN < 3; ++idN)
	    {
	      nod0 = nodes[idN];
	      if((nod0->flags & WLZ_CMESH_NOD_FLAG_KNOWN) == 0)
	      {
		WLZ_VTX_2_SUB(dsp, *(seedPos + idS), nod0->pos);
		*(times + nod0->idx) = WLZ_VTX_2_LENGTH(dsp);
		if(speeds != NULL)
		{
		  *(times + nod0->idx) /= *(speeds + nod0->idx);
		}
		nod0->flags |= WLZ_CMESH_NOD_FLAG_ACTIVE |
			       WLZ_CMESH_NOD_FLAG_KNOWN;
		errNum = WlzCMeshFMarQInsert2D(queue, times, nod0);
	      }
	    }
	  }
	}
        ++idS;
      }
    }
    else
    {
      idN = 0;
      while((errNum == WLZ_ERR_NONE) && (idN < mesh->res.nod.maxEnt))
      {
	nod0 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	if((nod0->idx >= 0) &&
	   ((nod0->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0))
	{
          *(times + idN) = 0.0;
	  nod0->flags |= WLZ_CMESH_NOD_FLAG_ACTIVE |
	                 WLZ_CMESH_NOD_FLAG_KNOWN;
	  errNum = WlzCMeshFMarQInsert2D(queue, times, nod0);
	}
        ++idN;
      }
    }
  }
  /* Until the queue is empty: Pop the node with lowest priority (ie
   * lowest time) from the queue, and process it. */
  while((errNum == WLZ_ERR_NONE) &&
	((qEnt = WlzCMeshFMarQPopTail(queue)) != NULL))
  {
    /* Find all neighbouring nodes that are neither active nor upwind.
     * For each of these neighbouring nodes, compute their time, set
     * them to active and insert them into the queue.*/
    nod0 = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, qEnt->nodIdx);
    edu0 = edu1 = nod0->edu;
    WlzCMeshFMarQEntryFree(queue, qEnt);
    do
    {
      nod1 = edu1->next->nod;
      if((nod1->flags & WLZ_CMESH_NOD_FLAG_UPWIND) == 0)
      {
	/* Node may have known time if active or unknown time if downwind. */
        WlzCMeshFMarCompute2D(nod1, nod0,
	                      (speeds == NULL)? 1.0: *(speeds + nod1->idx),
			      times);
        nod1->flags |= WLZ_CMESH_NOD_FLAG_ACTIVE;
        errNum = WlzCMeshFMarQInsert2D(queue, times, nod1);
      }
      edu1 = edu1->nnxt;
    } while(edu1 != edu0);
    /* Set the current node to upwind. */
    if(errNum == WLZ_ERR_NONE)
    {
      nod0->flags = (nod0->flags & ~(WLZ_CMESH_NOD_FLAG_ACTIVE)) |
		   WLZ_CMESH_NOD_FLAG_UPWIND;
    }
  }
  /* Clear up. */
  WlzCMeshFMarQFree(queue);
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
* \param	times			Array of wavefront arrival times
* 					indexed by the mesh node indices,
* 					which will be set on return.
* \param	speeds			Array of wavefront propagation
* 					speeds indexed by the mesh node
* 					indices, may be NULL in which case
* 					a constant speed of 1.0 is used.
* 					Speeds must all be > 0.
* \param	nSeeds			Number of seed nodes, if \f$<\f$ 1
* 					then all boundary nodes of the
* 					given mesh are used as seed nodes.
* \param	seedPos			Array of seed positions, may be
* 					NULL iff the number of seed nodes
* 					is \f$<\f$ 1. It is an error if
* 					any seeds are not within the
* 					mesh.
*/
WlzErrorNum	WlzCMeshFMarNodes3D(WlzCMesh3D *mesh,
				double *times, double *speeds,
				int nSeeds, WlzDVertex3 *seedPos)
{
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  return(errNum);
}


/*!
* \return	New constrained mesh node priority queue.
* \ingroup	WlzMesh
* \brief	Constructs a new constrained mesh node priority queue
* 		with room allocated for at least the given number of
* 		entries.
* \param	nEntries		Minimum number of entries for
* 					the queue.
*/
static WlzCMeshFMarQ *WlzCMeshFMarQNew(int nEntries)
{
  int		idE;
  WlzCMeshFMarQ	*queue;
  WlzCMeshFMarQEntry *ent;
  const size_t	minEntries = 1024;

  queue = (WlzCMeshFMarQ *)AlcMalloc(sizeof(WlzCMeshFMarQ));
  if(queue)
  {
    queue->max = (nEntries < minEntries)? minEntries: nEntries;
    if(((queue->entries = (WlzCMeshFMarQEntry *)
                          AlcMalloc(sizeof(WlzCMeshFMarQEntry) *
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
	ent->next = idE + 1;
	ent->prev = idE - 1;
        ent->hashNxt = -1;
	ent->nodIdx = -1;
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
  WlzCMeshFMarQEntry *ent;
  size_t	newMax;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const size_t	minEntryInc = 1024;

  newMax = queue->max + minEntryInc;
  if(minEnt > newMax)
  {
    newMax = minEnt;
  }
  if(((queue->entries = (WlzCMeshFMarQEntry *)
                        AlcRealloc(queue->entries,
				   sizeof(WlzCMeshFMarQEntry) *
				   newMax)) == NULL) ||
     ((queue->buckets = (int *)
			AlcRealloc(queue->buckets,
				   sizeof(int) * newMax)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    for(idE = queue->max; idE < newMax; ++idE)
    {
      ent = queue->entries + idE;
      ent->next = idE + 1;
      ent->prev = idE - 1;
      ent->hashNxt = -1;
      ent->nodIdx = -1;
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
* 		The entry must already have been unlinked.
* \param	queue			The queue.
* \param	ent			The queue entry to free.
*/
static void	WlzCMeshFMarQEntryFree(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEntry *qEnt)
{
  qEnt->next = queue->free;
  queue->free = qEnt - queue->entries;
  --(queue->cnt);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Inserts the given  node from a 2D constraied mesh into
* 		the given priority queue. The given times must be valid
* 		for the given node and any upwind nodes.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	times			Node times.
* \param	nod			Given node to insert into the queue.
*/
static WlzErrorNum WlzCMeshFMarQInsert2D(WlzCMeshFMarQ *queue,
				double *times, WlzCMeshNod2D *nod)
{
  int		idE,
  		idH;
  WlzCMeshFMarQEntry *ent,
  		*prevEnt;
  WlzErrorNum errNum = WLZ_ERR_NONE;

  /* Check for node already in queue. If the entry already exists remove it
   * from the queue and the hash table, but don't put it on the free stack. */
  ent = NULL;
  idH = WlzCMeshFMarQHashFn(nod->idx, queue->max);
  idE = queue->buckets[idH];
  if(idE >= 0)
  {
    prevEnt = NULL;
    ent = queue->entries + idE;
    while(ent->nodIdx != nod->idx)
    {
      prevEnt = ent;
      ent = queue->entries + ent->hashNxt;
    }
    /* If entry fouund for node unlink it and remove the hash table entry. */
    if(ent->nodIdx == nod->idx)
    {
      WlzCMeshFMarQEntryUnlink(queue, ent);
      if(prevEnt != NULL)
      {
        prevEnt->hashNxt = ent->hashNxt;
      }
      else
      {
        queue->buckets[idH] = -1;
      }
    }
    else
    {
      ent = NULL;
    }
  }
  if(ent == NULL)
  {
    /* Node entry not found: Create a new entry. */
    ent = WlzCMeshFMarQNewEntry(queue, &errNum);
  }
  /* Insert entry into the queue and add hash table entry. */
  if(errNum == WLZ_ERR_NONE)
  {
    idE = ent - queue->entries;
    ent->priority = WlzCMeshFMarQPriority(nod, times);
    WlzCMeshFMarQInsertEntry(queue, ent);
    ent->hashNxt = queue->buckets[idH];
    ent->nodIdx = nod->idx;
    queue->buckets[idH] = idE; 
    queue->last = idE;
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
static WlzCMeshFMarQEntry *WlzCMeshFMarQNewEntry(WlzCMeshFMarQ *queue,
				WlzErrorNum *dstErr)
{
  int		idE;
  WlzCMeshFMarQEntry *ent;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(queue->cnt >= queue->max)
  {
    /* Reallocate entries to get more and update the free list. */
    errNum = WlzCMeshFMarQRealloc(queue, queue->cnt);
  }
  else
  {
    /* Pop entry from free list. */
    idE = queue->free;
    ent = queue->entries + idE;
    queue->free = ent->next;
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
* \brief	Unlinks the given entry from the priority queue.
* 		If the head, last or tail index entry is unlinked the
* 		corresponding index is changed to the next or previous
* 		entry in the queue.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	ent			Entry to unlink.
*/
static void	WlzCMeshFMarQEntryUnlink(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEntry *ent)
{
  if(ent->prev >= 0)
  {
    (queue->entries + ent->prev)->next = ent->next;
  }
  else
  {
    queue->head = ent->next;
  }
  if(ent->next >= 0)
  {
    (queue->entries + ent->next)->prev = ent->prev;
  }
  else
  {
    queue->tail = ent->prev;
  }
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
  WlzCMeshFMarQEntry *ent;

  for(idE = 0; idE < queue->max; ++idE)
  {
    queue->buckets[idE] = -1;
  }
  idE = queue->head;
  while(idE >= 0)
  {
    ent = queue->entries + idE;
    idH = WlzCMeshFMarQHashFn(ent->nodIdx, queue->max);
    ent->hashNxt = queue->buckets[idH];
    queue->buckets[idH] = idE;
    idE = ent->next;
  }
}

/*!
* \return	Hash key.
* \ingroup	WlzMesh
* \brief	Simple hash function for integer values. The hash value
* 		returned is in the range [0-maxVal], but maxVal must be
* 		less than the largest of the generator primes (99999989).
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
* \return	Node priority.
* \ingroup	WlzMesh
* \brief	Priority for node in queue. The priority is increases with
* 		time.
* \param	nod			Node to compute priority for.
* \param	times			Array of times indexed by the
* 					mesh node indices.
*/
static double	WlzCMeshFMarQPriority(WlzCMeshNod2D *nod, double *times)
{
  return(*(times + nod->idx));
}

/*!
* \ingroup	WlzMesh
* \brief	Inserts the given entry into the given queue, knowing
* 		that the entry is valid and is not already in the queue.
* \param	queue			Given constrained mesh node priority
* 					queue.
* \param	ent0			Entry to insert.
*/
static void	WlzCMeshFMarQInsertEntry(WlzCMeshFMarQ *queue,
				WlzCMeshFMarQEntry *gEnt)
{
  int		idE,
  		idG;
  WlzCMeshFMarQEntry *ent;

  idG = gEnt - queue->entries;
  if(queue->last < 0)
  {
    /* Queue empty. */
    gEnt->next = gEnt->prev = -1;
    queue->head = queue->tail = queue->last = idG;
  }
  else
  {
    ent = queue->entries + queue->last;
    if(gEnt->priority > ent->priority)
    {
      /* Insert entry above last towards the head of the queue. */
      while((gEnt->priority > ent->priority) && (ent->prev >= 0))
      {
	ent = queue->entries + ent->prev;
      }
      idE = ent - queue->entries;
      if(ent->prev >= 0)
      {
	/* gEnt->priority <= ent->priority --GE--*/
	gEnt->prev = ent->prev;
	(queue->entries + gEnt->prev)->next = idG;
	gEnt->next = ent->next;
	ent->prev = idG;
      }
      else
      {
	/* gEnt->priority >= ent->priority */
	gEnt->prev = -1;
	if((gEnt->next = queue->head) >= 0)
	{
	  (queue->entries + gEnt->next)->prev = idG;
	}
	queue->head = idG;
      }
    }
    else
    {
      /* Insert entry at or below last (towards the tail). */
      while((gEnt->priority < ent->priority) && (ent->next >= 0))
      {
	ent = queue->entries + ent->next;
      }
      idE = ent - queue->entries;
      if(ent->next >= 0)
      {
	/* gEnt->priority <= ent->priority: --EG-- */
        gEnt->next = ent->next;
	(queue->entries + gEnt->next)->prev = idG;
	gEnt->prev = idE;
	ent->next = idG;
      }
      else
      {
	gEnt->next = -1;
	if((gEnt->prev = queue->tail) >= 0)
	{
	  (queue->entries + gEnt->prev)->next = idG;
	}
	queue->tail = idG;
      }
    }
  }
}

/*!
* \return	Unlinked entry.
* \ingroup	WlzMesh
* \brief	Unlinks the entry with the highest priority from the
* 		given queue but does not free it.
* \param	queue			Given priority queue.
*/
static WlzCMeshFMarQEntry *WlzCMeshFMarQPopTail(WlzCMeshFMarQ *queue)
{
  WlzCMeshFMarQEntry *ent = NULL;

  if(queue->tail >= 0)
  {
    ent = queue->entries + queue->tail;
    WlzCMeshFMarQEntryUnlink(queue, ent);
  }
  return(ent);
}

/*!
* \ingroup	WlzMesh
* \brief	Computes wavefront propagation time for the given unknown
* 		node.
* 		Given a pair of nodes that are connected by a single edge
* 		in a 2D conforming mesh, the first of which has an unknown
* 		and the second a known wavefront propagation time, this
* 		function computes the unknown time.
* 		A mesh element is found which includes the two given nodes
* 		and if possible has a third node for which the time is known.
* 		This element has nodes nod0, nod1 and nod2. Internal angles
* 		at the nodes of phi0, phi1 and ph2. Edge lengths opposite
* 		to the similarly numbered nodes of len0, len1 and len2.
* 		The solution is similar to that in "Fast Sweeping Methods
* 		For Eikonal equations On triangular meshes", Jianliang Quian,
* 		etal, SIAM journal on Mumerical Analysis, Vol 45, pp 83-107,
* 		2007. But uses simple interpolation in the case of phi2
* 		being obtuse.
* \param	nod2			Unknown node.
* \param	nod0			A known node directly connected to
* 					the unknown node.
* \param	s2			Wavefront propagation speed at
* 					the unknown node, must be greater
* 					than zero.
* \param	times			Array of times indexed by the
* 					mesh node indices, which will be
* 					set for the unknown node on return.
*/
static void	WlzCMeshFMarCompute2D(WlzCMeshNod2D *nod2, WlzCMeshNod2D *nod0,
				double s2, double *times)
{
  int		idM,
  		idN;
  double	d0,
		d1,
		is2,
		theta;
  double	len[3],
		lenSq[3],
  		phi[3],
		time[3];
  WlzDVertex2	del,
  		pos;
  WlzCMeshNod2D	*nod1,
		*nod3;
  WlzCMeshEdgU2D *edu;
  WlzCMeshElm2D	*elm;

  /* Find an element which is common to both the given nodes. */
  edu = nod2->edu;
  is2 = 1.0 / s2;
  do
  {
    edu = edu->nnxt;
  } while((edu != nod2->edu) && (edu->next->nod != nod0));
  if(edu->next->nod != nod0)
  {
    /* Nodes are boundary and directed edge runs from nod0 to nod2
     * so missed it when looking for directed edge from nod2 to nod0
     * so look from nod0 to nod2 instead. */
    edu = nod0->edu;
    do
    {
      edu = edu->nnxt;
    }
    while((edu != nod0->edu) && (edu->next->nod != nod2));
  }
  elm = edu->elm;
  /* Get nod1 from element. */
  for(idN = 0; idN < 3; ++idN)
  {
    nod1 = elm->edu[idN].nod;
    if((nod1 != nod2) && (nod1 != nod0))
    {
      break;
    }
  }
  if(((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) == 0) &&
     (edu->opp != NULL) && (edu->opp != edu))
  {
    /* If nod1 has an unknown time then check for elm on opposite edge. */
    for(idN = 0; idN < 3; ++idN)
    {
      edu = elm->edu + idN;
      if(((edu->nod == nod0) && (edu->next->nod == nod2)) ||
         ((edu->nod == nod2) && (edu->next->nod == nod0)))
      {
	if((edu->opp != NULL) && (edu->opp != edu))
	{
          elm = edu->opp->elm;
	  /* Get nod1 from element. */
	  for(idM = 0; idM < 3; ++idM)
	  {
	    nod3 = elm->edu[idM].nod;
	    if((nod3 != nod2) && (nod3 != nod0))
	    {
	      if((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) != 0)
	      {
		nod1 = nod3;
	      }
	      break;
	    }
	  }
	}
        break;
      }
    }
  }
  if((nod1->flags & WLZ_CMESH_NOD_FLAG_KNOWN) == 0)
  {
    /* Only have one node with known time in this element. Compute the time
     * using just the edge length. */
    WLZ_VTX_2_SUB(del, nod2->pos, nod0->pos);
    lenSq[1] = WLZ_VTX_2_SQRLEN(del);
    len[1] = sqrt(lenSq[1]);
    time[2] = *(times + nod0->idx) + is2 * len[1];
  }
  else
  {
    /* Have two known nodes (nod0 and nod1) and one unknown node (nod2).
     * Make sure the time for nod0 is <= time for nod1. */
    time[0] = *(times + nod0->idx);
    time[1] = *(times + nod1->idx);
    if(time[0] > time[1])
    {
      nod3 = nod0;
      nod0 = nod1;
      nod1 = nod3;
      time[0] = time[1];
      time[1] = *(times + nod1->idx);
    }
    WLZ_VTX_2_SUB(del, nod0->pos, nod1->pos);
    lenSq[2] = WLZ_VTX_2_SQRLEN(del);
    len[2] = sqrt(lenSq[2]);
    if(len[2] < DBL_EPSILON)
    {
      /* Element is degenerate so just use edge length to compute time. */
      time[2] = *(times + nod0->idx) + is2 * len[1];
    }
    else
    {
      /* nod0 and nod1 both have known times and len2 is greater than
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
	 * time1 = (time0 + time1)/2, pos1 = (pos1 + pos0)/2. */
        WLZ_VTX_2_ADD(pos, nod0->pos, nod1->pos);
	WLZ_VTX_2_SCALE(pos, pos, 0.5);
	time[1] = 0.5 * (time[1] + time[0]);
	WLZ_VTX_2_SUB(del, nod0->pos, pos);
	lenSq[2] = WLZ_VTX_2_SQRLEN(del);
	len[2] = sqrt(lenSq[2]);
	WLZ_VTX_2_SUB(del, pos, nod2->pos);
	lenSq[0] = WLZ_VTX_2_SQRLEN(del);
	len[0] = sqrt(lenSq[0]);
      }
      if((len[2] * is2) > (time[1] - time[0]))
      {
        theta = asin(is2 * (time[1] -
	                    time[0]) / len[2]);
	phi[0] = acos((lenSq[1] + lenSq[2] - lenSq[0]) /
	              (2.0 * len[1] * len[2]));
	phi[1] = acos((lenSq[2] + lenSq[0] - lenSq[1]) /
	              (2.0 * len[2] * len[0]));
	if(((theta > DBL_EPSILON) &&
	    (theta > (phi[1] - ALG_M_PI_2)) && 
	    ((ALG_M_PI_2 - phi[0]) > theta)) ||
	   ((DBL_EPSILON > theta) &&
	    (theta > phi[1] - ALG_M_PI_2) &&
	    ((ALG_M_PI_2 - phi[0]) > theta)))
	{
	  d0 = len[0] * sin(phi[1] - theta); /* h0 */
	  d1 = len[1] * sin(phi[0] + theta); /* h1 */
	  time[2] = 0.5 * ((d0 * is2 + time[1]) + (d1 * is2 + time[0]));
	}
	else
	{
	  d0 = time[0] + len[1] * is2;
	  d1 = time[1] + len[0] * is2;
	  time[2] = ALG_MIN(d0, d1);
	}
      }
      else
      {
	d0 = time[0] + len[1] * is2;
	d1 = time[1] + len[0] * is2;
	time[2] = ALG_MIN(d0, d1);
      }
    }
  }
  if(time[2] < *(times + nod2->idx))
  {
    *(times + nod2->idx) = time[2];
  }
}
