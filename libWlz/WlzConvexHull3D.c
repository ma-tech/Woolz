#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvexHull3D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzConvexHull3D.c
* \author       Bill Hill
* \date         April 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Functions for computing 3D convex hull domains.
* \ingroup	WlzConvexHull
*/

#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <Wlz.h>


/*!
* \def	WLZ_CONVHULL_EPS	(1.0e-06)
* \brief	Tollerance value for lengths.
*/
#define WLZ_CONVHULL_EPS	(1.0e-06)

/*!
* \struct	_WlzConvHullArc
* \brief	A conflict arc connecting a face and a vertex these form
* 		circular doubly linked lists for the faces and vertices.
* 		Typedef: ::WlzConvHullArc
*/
typedef struct _WlzConvHullArc
{
  struct _WlzConvHullFce *fce;		/*!< Face of conflict. */
  struct _WlzConvHullVtx *vtx;		/*!< Vertex in conflict. */
  struct _WlzConvHullArc *prvVtx; 	/*!< Previous vertex arc. */
  struct _WlzConvHullArc *nxtVtx; 	/*!< Next vertex arc, also used for
  					     the stack of free arcs in the
					     arc pool. */
  struct _WlzConvHullArc *prvFce;	/*!< Previous face arc. */
  struct _WlzConvHullArc *nxtFce;	/*!< Next face arc. */
} WlzConvHullArc;

/*!
* \struct	_WlzConvHullVtx
* \brief	A vertex of the convex hull for use in a convex hull
* 		workspace.
* 		Typedef: ::WlzConvHullVtx
*/
typedef struct _WlzConvHullVtx
{
  int			idx;		/*!< Index of the vertex. */
  int			cvx;		/*!< Vertex is on the convex hull. */
  struct _WlzConvHullArc *arc;		/*!< Conflict list. */
  struct _WlzConvHullVtx *nxt;		/*!< Next vertex in list. */
  struct _WlzConvHullVtx *prv;		/*!< Previous vertex in list. */
} WlzConvHullVtx;

/*!
* \struct	_WlzConvHullFce
* \brief	A face of the convex hull for use in a convex hull
* 		workspace.
* 		Typedef: ::WlzConvHullFce
*/
typedef struct _WlzConvHullFce
{
  int			vtx[3];		/*!< Indices of the vertices. */
  struct _WlzConvHullFce *opp[3]; 	/*!< Opposite faces, opp[i] is the
  					     face opposite on the edge in
					     this face directed from vtx[i]
					     to vtx[(i + 1)%3]. */
  WlzDVertex3		nrm;		/*!< Outward directed normal for
  					     this face. */
  struct _WlzConvHullArc *arc;		/*!< Conflict list. */
  struct _WlzConvHullFce *nxt;		/*!< Next face in list. */
  struct _WlzConvHullFce *prv;		/*!< Previous face in list. */
} WlzConvHullFce;

/*!
* \struct	_WlzConvHullArcPool
* \brief	A pool of conflict arcs.
* 		Typedef: ::WlzConvHullArcPool
*/
typedef struct _WlzConvHullArcPool
{
  int			freeCnt;	/*!< Number of arcs in free pool. */
  WlzConvHullArc 	*free;		/*!< Next arc available. */
  AlcBlockStack		*blkStk;	/*!< Block stack for allocation. */
} WlzConvHullArcPool;

/*!
* \struct	_WlzConvHullFcePool
* \brief	A pool of faces.
* 		Typedef: ::WlzConvHullFcePool
*/
typedef struct _WlzConvHullFcePool
{
  int			freeCnt;	/*!< Number of faces in free pool. */
  WlzConvHullFce 	*free;		/*!< Next face available. */
  AlcBlockStack		*blkStk;	/*!< Block stack for allocation. */
} WlzConvHullFcePool;

/*!
* \struct	_WlzConvHullVtxPool
* \brief	A pool of vertices.
* 		Typedef: ::WlzConvHullVtxPool
*/
typedef struct _WlzConvHullVtxPool
{
  int			freeCnt;	/*!< Number of arcs in free pool. */
  WlzConvHullVtx 	*free;		/*!< Next vertex available. */
  AlcBlockStack		*blkStk;	/*!< Block stack for allocation. */
} WlzConvHullVtxPool;

/*!
* \struct	_WlzConvHullHorEdg
* \brief	A horizon edge as seen from an external vertex.
* 		Typedef: ::WlzConvHullHorEdg
*/
typedef struct _WlzConvHullHorEdg
{
  WlzConvHullFce	*fce;		/*!< Face on horizon. */
  int			edg;		/*!< Edge on horizon face. */
} WlzConvHullHorEdg;

/*!
* \struct	_WlzConvHullWSp3
* \brief	A workspace for computing the 3D convex hull of vertices
* 		with integral coordinates.
* 		Typedef: ::WlzConvHullWSp3
*/
typedef struct _WlzConvHullWSp3
{
  int			nHorizon;	/*!< Number of horizon edges. */
  int			maxHorizon;	/*!< Maximum horizon edge space. */
  int			nFceBuf;	/*!< Number of faces in face buffer. */
  int			maxFceBuf;	/*!< Maximum face buffer space. */
  WlzVertexType		vtxType;	/*!< Vertex type, either WLZ_VERTEX_I3
  					     or WLZ_VERTEX_D3. */
  WlzVertexP		vtxPos;		/*!< Pointer to vertex positions. */
  int			*vtxPrm;	/*!< Vertex permutation table. */
  WlzConvHullFce	*fceLst;	/*!< List of faces in the convex
                                             hull so far. */
  WlzConvHullVtx	*vtxLst;	/*!< List of vertices in the convex
  					     hull so far. */
  WlzConvHullVtx	*vtxQue;	/*!< Queue of vertices waiting to
  					     be processed. */
  WlzConvHullFce	**fceBuf;	/*!< Buffer for faces to be added
  					     or deleted. */
  WlzConvHullHorEdg	*horizon;	/*!< Buffer for horizon edges. */
  WlzConvHullArcPool 	arcPool;	/*!< Pool for arc allocation. */
  WlzConvHullFcePool 	fcePool;	/*!< Pool for face allocation. */
  WlzConvHullVtxPool 	vtxPool;	/*!< Pool for vertex allocation. */
} WlzConvHullWSp3;

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Expand the pool of conflict arcs to at least the given
* 		minimum number of free arcs.
* \param	pool			The arc pool.
* \param	minElm			Given minimum number of elements.
*/
static WlzErrorNum		WlzConvHullArcPoolExpand(
				  WlzConvHullArcPool *pool,
				  int minElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	incStep = 4096;

  if(pool->freeCnt < minElm)
  {
    int		nElm;
    AlcBlockStack *newBlk;

    nElm = ((minElm / incStep) + 1) * incStep;
    newBlk = AlcBlockStackNew(nElm, sizeof(WlzConvHullArc),
                              pool->blkStk, NULL);
    if(newBlk == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int		i;
      WlzConvHullArc *arc;

      pool->freeCnt += nElm;
      arc = (WlzConvHullArc *)(newBlk->elements);
      pool->blkStk = newBlk;
      for(i = 0; i < nElm; ++i)
      {
	arc->nxtVtx = pool->free;
	pool->free = arc;
	++arc;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Expand the pool of faces to at least the given
* 		minimum number of free faces.
* \param	pool			The face pool.
* \param	minElm			Given minimum number of elements.
*/
static WlzErrorNum		WlzConvHullFcePoolExpand(
				  WlzConvHullFcePool *pool,
				  int minElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	incStep = 4096;

  if(pool->freeCnt < minElm)
  {
    int		nElm;
    AlcBlockStack *newBlk;

    nElm = ((minElm / incStep) + 1) * incStep;
    newBlk = AlcBlockStackNew(nElm, sizeof(WlzConvHullFce),
    		              pool->blkStk, NULL);
    if(newBlk == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int		i;
      WlzConvHullFce *fce;

      pool->freeCnt += nElm;
      fce = (WlzConvHullFce *)(newBlk->elements);
      pool->blkStk = newBlk;
      for(i = 0; i < nElm; ++i)
      {
	fce->nxt = pool->free;
	pool->free = fce;
	++fce;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Expand the pool of vertices to at least the given
* 		minimum number of free vertices.
* \param	pool			The vertex pool.
* \param	minElm			Given minimum number of elements.
*/
static WlzErrorNum		WlzConvHullVtxPoolExpand(
				  WlzConvHullVtxPool *pool,
				  int minElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	incStep = 4096;

  if(pool->freeCnt < minElm)
  {
    int		nElm;
    AlcBlockStack *newBlk;

    nElm = ((minElm / incStep) + 1) * incStep;
    newBlk = AlcBlockStackNew(nElm, sizeof(WlzConvHullVtx),
    			      pool->blkStk, NULL);
    if(newBlk == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int		i;
      WlzConvHullVtx *vtx;

      pool->freeCnt += nElm;
      vtx = (WlzConvHullVtx *)(newBlk->elements);
      pool->blkStk = newBlk;
      for(i = 0; i < nElm; ++i)
      {
	vtx->nxt = pool->free;
	pool->free = vtx;
	++vtx;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Expand the face buffer to at least the given
* 		minimum number of face pointers.
* \param	pool			Convex hull workspace.
* \param	minElm			Given minimum number of elements.
*/
static WlzErrorNum		WlzConvHullFceBufExpand(
				  WlzConvHullWSp3 *wSp,
				  int nMaxElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	bufInc = 4096;

  if(wSp->maxFceBuf < nMaxElm)
  {
    wSp->maxFceBuf += bufInc;
    if((wSp->fceBuf = (WlzConvHullFce **)
		      AlcRealloc(wSp->fceBuf,
			  sizeof(WlzConvHullFce *) *
			  wSp->maxFceBuf)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Expand the horizon edge buffer to at least the given
* 		minimum number of horizon edges.
* \param	pool			Convex hull workspace.
* \param	minElm			Given minimum number of elements.
*/
static WlzErrorNum		WlzConvHullHorizonExpand(
				  WlzConvHullWSp3 *wSp,
				  int nMaxElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	bufInc = 4096;

  if(wSp->maxHorizon < nMaxElm)
  {
    wSp->maxHorizon += bufInc;
    if((wSp->horizon = (WlzConvHullHorEdg *)
		      AlcRealloc(wSp->horizon,
			  sizeof(WlzConvHullHorEdg) *
			  wSp->maxHorizon)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  return(errNum);
}

/*!
* \return	New conflict arc.
* \ingroup	WlzConvexHull
* \brief	Get a conflict arc from the conflict arc pool. All fields of
* 		the arc are cleared to zero.
* \param	wSp			Convex hull workspace.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullArc		*WlzConvHullNewArc(
				  WlzConvHullWSp3 *wSp,
				  WlzErrorNum *dstErr)
{
  WlzConvHullArc *arc = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->arcPool.free == NULL)
  {
    errNum = WlzConvHullArcPoolExpand(&(wSp->arcPool), 1);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    --(wSp->arcPool.freeCnt);
    arc = wSp->arcPool.free;
    wSp->arcPool.free = wSp->arcPool.free->nxtVtx;
    if(arc == wSp->arcPool.free)
    {
      wSp->arcPool.free = NULL;
    }
    (void )memset(arc, 0, sizeof(WlzConvHullArc));
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(arc);
}

/*!
* \return	New face.
* \ingroup	WlzConvexHull
* \brief	Get a face from the face pool. All fields of the face are
* 		cleared to zero.
* \param	wSp			Convex hull workspace.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullFce		*WlzConvHullNewFce(
				  WlzConvHullWSp3 *wSp,
				  WlzErrorNum *dstErr)
{
  WlzConvHullFce *fce = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->fcePool.free == NULL)
  {
    errNum = WlzConvHullFcePoolExpand(&(wSp->fcePool), 1);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    --(wSp->fcePool.freeCnt);
    fce = wSp->fcePool.free;
    wSp->fcePool.free = wSp->fcePool.free->nxt;
    if(fce == wSp->fcePool.free)
    {
      wSp->fcePool.free = NULL;
    }
    (void )memset(fce, 0, sizeof(WlzConvHullFce));
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(fce);
}

/*!
* \return	New vertex.
* \ingroup	WlzConvexHull
* \brief	Get a vertex from the vertex pool. All fields of the vertex
* 		are cleared to zero.
* \param	wSp			Convex hull workspace.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullVtx		*WlzConvHullNewVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzErrorNum *dstErr)
{
  WlzConvHullVtx *vtx = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->vtxPool.free == NULL)
  {
    errNum = WlzConvHullVtxPoolExpand(&(wSp->vtxPool), 1);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    --(wSp->vtxPool.freeCnt);
    vtx = wSp->vtxPool.free;
    wSp->vtxPool.free = wSp->vtxPool.free->nxt;
    if(vtx == wSp->vtxPool.free)
    {
      wSp->vtxPool.free = NULL;
    }
    (void )memset(vtx, 0, sizeof(WlzConvHullVtx));
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vtx);
}

/*!
* \return	New convex hull workspace.
* \ingroup	WlzConvexHull
* \brief	Creates a new 3D convex hull workspace.
* \param	vtxType			Type of vertex, must be either
* 					WLZ_VERTEX_I3 or WLZ_VERTEX_D3.
* \param	nVtx			Number of vertices.
* \param	vtx			The vertex positions.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullWSp3		*WlzConvHullMakeWSp3(
				  WlzVertexType vtxType,
				  int nVtx,
				  WlzVertexP vtx,
				  WlzErrorNum *dstErr)
{
  WlzConvHullWSp3 *wSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((vtxType != WLZ_VERTEX_I3) && (vtxType != WLZ_VERTEX_D3))
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else if(((wSp = AlcCalloc(1, sizeof(WlzConvHullWSp3))) == NULL) ||
          ((wSp->vtxPrm = (int *)AlcMalloc(nVtx * sizeof(int))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    wSp->vtxType = vtxType;
    wSp->vtxPos = vtx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(wSp);
}

/*!
* \ingroup	WlzConvexHull
* \brief	Frees the given workspace.
* \param	wSp			Given convex hull workspace.
*/
static void			WlzConvHullFreeWSp3(
				  WlzConvHullWSp3 *wSp)
{
  if(wSp)
  {
    AlcFree(wSp->vtxPrm);
    AlcFree(wSp->fceBuf);
    AlcFree(wSp->horizon);
    AlcBlockStackFree(wSp->arcPool.blkStk);
    AlcBlockStackFree(wSp->fcePool.blkStk);
    AlcBlockStackFree(wSp->vtxPool.blkStk);
    AlcFree(wSp);
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Frees an unused conflict arc by returning it to it's pool.
* \param	wSp			Convex hull workspace.
* \param	arc			Unused conflict arc to be freed.
*/
static void			WlzConvHullFreeArc(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullArc *arc)
{
  if(arc)
  {
    WlzConvHullArcPool *pool;

    pool = &(wSp->arcPool);
    (void )memset(arc, 0, sizeof(WlzConvHullArc));
    arc->nxtVtx = pool->free;
    pool->free = arc;
    ++(pool->freeCnt);
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Frees an unused face by returning it to it's pool.
* \param	wSp			Convex hull workspace.
* \param	arc			Unused face to be freed.
*/
static void			WlzConvHullFreeFce(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce)
{
  if(fce)
  {
    WlzConvHullFcePool *pool;

    pool = &(wSp->fcePool);
    (void )memset(fce, 0, sizeof(WlzConvHullFce));
    fce->nxt = pool->free;
    pool->free = fce;
    ++(pool->freeCnt);
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Frees an unused vertex by returning it to it's pool.
* \param	wSp			Convex hull workspace.
* \param	arc			Unused vertex to be freed.
*/
static void			WlzConvHullFreeVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx *vtx)
{
  if(vtx)
  {
    WlzConvHullVtxPool *pool;

    pool = &(wSp->vtxPool);
    (void )memset(vtx, 0, sizeof(WlzConvHullVtx));
    vtx->nxt = pool->free;
    pool->free = vtx;
    ++(pool->freeCnt);
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Push the given vertex onto the pending vertex queue.
* 		Simple randomisation seems better than a queue ordered
* 		by maximum distance from a face.
* \param	wSp			Convex hull workspace.
* \param	vtx			Vertex to push.
*/
static void			WlzConHullPushPendingVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx *vtx)
{
  if(wSp->vtxQue == NULL)
  {
    wSp->vtxQue = vtx->prv = vtx->nxt = vtx;
  }
  else
  {
    vtx->nxt = wSp->vtxQue;
    vtx->prv = wSp->vtxQue->prv;
    vtx->prv->nxt = vtx;
    vtx->nxt->prv = vtx;
    wSp->vtxQue = vtx;
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Adds the given face to the list of faces currently in the
* 		convex hull.
* \param	wSp			Convex hull workspace.
* \param	fce			Face to add.
*/
static void			WlzConHullAddConvexFce(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce)
{
  if(wSp->fceLst == NULL)
  {
    wSp->fceLst = fce->prv = fce->nxt = fce;
  }
  else
  {
    fce->nxt = wSp->fceLst;
    fce->prv = wSp->fceLst->prv;
    fce->nxt->prv = fce->prv->nxt = fce;
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Adds the given vertex to the list of vertices currently in the
* 		convex hull.
* \param	wSp			Convex hull workspace.
* \param	fce			Face to add.
*/
static void			WlzConHullAddConvexVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx *vtx)
{
  if(wSp->vtxLst == NULL)
  {
    wSp->vtxLst = vtx->prv = vtx->nxt = vtx;
  }
  else
  {
    vtx->cvx = 1;
    vtx->nxt = wSp->vtxLst;
    vtx->prv = wSp->vtxLst->prv;
    vtx->nxt->prv = vtx->prv->nxt = vtx;
  }
}

/*!
* \return	Non-zero if the vertex is behind the face, ie the vertex
* 		lies on the face or on the convex hull side of the face.
* \ingroup	WlzConvexHull
* \brief	Determines whether the vertex is behind the face.
* 		If a vertex is behind all faces of a convex hull then it
* 		is within (or possibly on) te convex hull.
* 		This function can consume most of the CPU cycles when
* 		computing a convex hull.
* \param	wSp			Convex hull workspace.
* \param	fce			Given face.
* \param	vtx			Given vertex index.
*/
static int 			WlzConvHullFceBehind(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce,
				  int vtx)
{
  int		b;
  double	p;
  WlzDVertex3	u;

  if(wSp->vtxType == WLZ_VERTEX_I3)
  {
    WlzIVertex3   v0,
		  v1;

    v0 = wSp->vtxPos.i3[fce->vtx[0]];
    v1 = wSp->vtxPos.i3[vtx];
    WLZ_VTX_3_SUB(u, v1, v0);
  }
  else /* wSp->vtxType == WLZ_VERTEX_D3 */
  {
    WlzDVertex3   v0,
		  v1;

    v0 = wSp->vtxPos.d3[fce->vtx[0]];
    v1 = wSp->vtxPos.d3[vtx];
    WLZ_VTX_3_SUB(u, v1, v0);
  }
  p = WLZ_VTX_3_DOT(fce->nrm, u);
  b = (p < WLZ_CONVHULL_EPS);
  return(b);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Adds a conflict arc between the given face and vertex.
* 		Arcs are maintained as joint circular doubly linked lists
* 		through the faces and vertices.
* \param	wSp			Convex hull workspace.
* \param	fce			Given face.
* \param	vtx			Given vertex.
*/
static WlzErrorNum		WlzConvHullConfAdd(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce,
				  WlzConvHullVtx *vtx)
{
  WlzConvHullArc *newArc;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  newArc = WlzConvHullNewArc(wSp, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    newArc->fce = fce;
    newArc->vtx = vtx;
    /* Insert arc into face list. */
    if(fce->arc == NULL)
    {
      fce->arc = newArc;
      newArc->nxtFce = newArc->prvFce = newArc;
    }
    else
    {
      newArc->nxtFce = fce->arc;
      newArc->prvFce = fce->arc->prvFce;
      fce->arc->prvFce->nxtFce = newArc;
      fce->arc->prvFce = newArc;
    }
    /* Insert arc into vertex list. */
    if(vtx->arc == NULL)
    {
      vtx->arc = newArc;
      newArc->nxtVtx = newArc->prvVtx = newArc;
    }
    else
    {
      newArc->nxtVtx = vtx->arc;
      newArc->prvVtx = vtx->arc->prvVtx;
      vtx->arc->prvVtx->nxtVtx = newArc;
      vtx->arc->prvVtx = newArc;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvHull
* \brief	Deletes (and recycles to the arc pool) a conflict arc
* 		from the doubly linked lists through the faces and vertices.
* 		If this is the only arc in a conflict list then the conflict
* 		list will be set to NULL.
* \param	wSp			Convex hull workspace.
* \param	arc			Arc to be deleted.
*/
static void			WlzConvHullDelArc(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullArc *arc)
{
  WlzConvHullArc *nxtArc;
  WlzConvHullFce *fce;
  WlzConvHullVtx *vtx;

  fce = arc->fce;
  vtx = arc->vtx;
  /* Unlink the arc from the confict list. */
  if(arc == fce->arc)
  {
    if(arc == arc->nxtFce)
    {
      fce->arc = NULL;
    }
    else
    {
      nxtArc = arc->nxtFce;
      fce->arc = nxtArc;
      arc->prvFce->nxtFce = nxtArc;
      nxtArc->prvFce = arc->prvFce;
    }
  }
  else
  {
    nxtArc = arc->nxtFce;
    arc->prvFce->nxtFce = nxtArc;
    nxtArc->prvFce = arc->prvFce;
  }
  if(arc == vtx->arc)
  {
    if(arc == arc->nxtVtx)
    {
      vtx->arc = NULL;
    }
    else
    {
      nxtArc = arc->nxtVtx;
      vtx->arc = nxtArc;
      arc->prvVtx->nxtVtx = nxtArc;
      nxtArc->prvVtx = arc->prvVtx;
    }
  }
  else
  {
    nxtArc = arc->nxtVtx;
    arc->prvVtx->nxtVtx = nxtArc;
    nxtArc->prvVtx = arc->prvVtx;
  }
  /* Put the arc back into the free pool. */
  WlzConvHullFreeArc(wSp, arc);
}

/*!
* \ingroup	WlzConvexHull
* \brief	Deletes a face and any conflicts it may have recycling
* 		the face and conflict arcs back to their pools.
* \param	wSp			Convex hull workspace.
* \param	head			Head of the list of faces, may
* 					be modified and may be NULL if the
* 					list becomes empty.
* \param	fce			Face to be deleted.
*/
static void			WlzConvHullDelFce(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce **head,
				  WlzConvHullFce *fce)
{
  WlzConvHullFce *nxtFce;

  /* Free all conflict arcs of this face. */
  while(fce->arc)
  {
    WlzConvHullDelArc(wSp, fce->arc->nxtFce);
  }
  /* Unlink the face from it's queue or list. */
  if(fce == *head)
  {
    if(fce == fce->nxt)
    {
      *head = NULL;
    }
    else
    {
      nxtFce = fce->nxt;
      *head = nxtFce;
      fce->prv->nxt = nxtFce;
      nxtFce->prv = fce->prv;
    }
  }
  else
  {
    nxtFce = fce->nxt;
    fce->prv->nxt = nxtFce;
    nxtFce->prv = fce->prv;
  }
  /* Free the face back to it's pool. */
  WlzConvHullFreeFce(wSp, fce);
}

/*!
* \ingroup	WlzConvexHull
* \brief	Unlink a vertex from a list but don't free it.
* \param	wSp			Convex hull workspace.
* \param	head			Head of the list of vertices, may
* 					be modified.
* \param	vtx			Vertex to be unlinked.
*/
static void 			WlzConvHullUnkVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx **head,
				  WlzConvHullVtx *vtx)
{
  WlzConvHullVtx *nxtVtx,
  		 *prvVtx;

  if(vtx == *head)
  {
    if(vtx == vtx->nxt)
    {
      *head = NULL;
    }
    else
    {
      nxtVtx = vtx->nxt;
      prvVtx = vtx->prv;
      prvVtx->nxt = nxtVtx;
      nxtVtx->prv = prvVtx;
      *head = nxtVtx;
    }
  }
  else
  {
    nxtVtx = vtx->nxt;
    prvVtx = vtx->prv;
    prvVtx->nxt = nxtVtx;
    nxtVtx->prv = prvVtx;
  }
}

/*!
* \ingroup	WlzConvexHull
* \brief	Deletes a vertex and any conflicts it may have recycling
* 		the vertex and conflict arcs back to their pools.
* \param	wSp			Convex hull workspace.
* \param	head			Head of the list of vertices, may
* 					be modified and may be NULL if the
* 					list becomes empty. The head may
* 					also be supplied as a NULL pointer
* 					which implies that the vertex is not
* 					in any list.
* \param	vtx			Vertex to be deleted.
*/
static void			WlzConvHullDelVtx(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx **head,
				  WlzConvHullVtx *vtx)
{
  /* Free all conflict arcs of this vertex. */
  while(vtx->arc)
  {
    WlzConvHullDelArc(wSp, vtx->arc->nxtVtx);
  }
  /* Unlink the vertex from it's queue or list. */
  if(head)
  {
    WlzConvHullUnkVtx(wSp, head, vtx);
  }
  /* Free the vertex back to it's pool. */
  WlzConvHullFreeVtx(wSp, vtx);
}

/*!
* \return	Top pending vertex unlinked from queue.
* \ingroup	WlzConvexHull
* \brief	Pops the top vertex from the pending vertex queue.
* \param	wSp			Convex hull workspace.
*/
static WlzConvHullVtx		*WlzConHullPopPendingVtx(
				  WlzConvHullWSp3 *wSp)
{
  WlzConvHullVtx	*vtx;

  if((vtx = wSp->vtxQue) != NULL)
  {
    if(vtx == vtx->nxt)
    {
      wSp->vtxQue = NULL;
    }
    else
    {
      wSp->vtxQue = vtx->nxt;
      wSp->vtxQue->prv = vtx->prv;
      wSp->vtxQue->prv->nxt = wSp->vtxQue;
    }
  }
  return(vtx);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Computes the face normal in place using the vertices of the
* 		face.
* \param	wSp			Convex hull workspace.
* \param	fce			Given face.
*/
static WlzErrorNum		WlzConvHullFceNrm(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce)
{
  double	len;
  WlzDVertex3 	nrm;
  WlzDVertex3 	u[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->vtxType == WLZ_VERTEX_I3)
  {
    WlzIVertex3 v[3];

    v[0] = wSp->vtxPos.i3[fce->vtx[0]];
    v[1] = wSp->vtxPos.i3[fce->vtx[1]];
    v[2] = wSp->vtxPos.i3[fce->vtx[2]];
    WLZ_VTX_3_SUB(u[0], v[1], v[0]);
    WLZ_VTX_3_SUB(u[1], v[2], v[1]);
  }
  else /* wSp->vtxType == WLZ_VERTEX_D3 */
  {
    WlzDVertex3 v[3];

    v[0] = wSp->vtxPos.d3[fce->vtx[0]];
    v[1] = wSp->vtxPos.d3[fce->vtx[1]];
    v[2] = wSp->vtxPos.d3[fce->vtx[2]];
    WLZ_VTX_3_SUB(u[0], v[1], v[0]);
    WLZ_VTX_3_SUB(u[1], v[2], v[1]);
  }
  WLZ_VTX_3_CROSS(nrm, u[0], u[1]);
  len = WLZ_VTX_3_LENGTH(nrm);
  if(len < WLZ_CONVHULL_EPS)
  {
    WLZ_VTX_3_ZERO(fce->nrm);
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    len = 1.0 / len;
    WLZ_VTX_3_SCALE(fce->nrm, nrm, len);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Builds the horizon edge array using the conflicting faces
* 		of the vertex which are put into the face buffer array.
* 		The horizon edges run in sequence and with a direction
* 		correct in the non conflicting faces (those which surround
* 		the conflicting ones).
* 		These arrays are expanded as required.
* \param	wSp			Convex hull workspace.
* \param	vtx			Given vertex.
*/
static WlzErrorNum		WlzConvHullBuildHorizon(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullVtx *vtx)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  wSp->nFceBuf = 0;
  wSp->nHorizon = 0;
  /* Make sure there's room for at least 3 horizon edges, the maximum
   * possible from a single face. This simplifies buffer expansion. */
  errNum = WlzConvHullHorizonExpand(wSp, 3);
  /* Build array of faces that conflict with the given vertex. */
  if((errNum == WLZ_ERR_NONE) && (vtx->arc != NULL))
  {
    WlzConvHullArc *arc;

    arc = vtx->arc;
    do
    {
      int 	f,
      		fceIn = 0;
      WlzConvHullFce *fce;

      fce = arc->fce;
      /* Check for face already in the face buffer. */
      for(f = 0; f < wSp->nFceBuf; ++f)
      {
        if(fce == wSp->fceBuf[f])
	{
	  fceIn = 1;
	  break;
	}
      }
      /* If not then add it. */
      if(fceIn == 0)
      {
	if(wSp->nFceBuf + 1 >= wSp->maxFceBuf)
	{
	  errNum = WlzConvHullFceBufExpand(wSp, wSp->nFceBuf + 1);
	  if(errNum != WLZ_ERR_NONE)
	  {
	      break;
	  }
	}
	wSp->fceBuf[wSp->nFceBuf] = fce;
	++(wSp->nFceBuf);
      }
      arc = arc->nxtVtx; 			/* Next arc for this vertex. */
    } while(arc != vtx->arc);
  }
  /* Build horizon from the faces that are in the face buffer. */
  if((errNum == WLZ_ERR_NONE) && (wSp->nFceBuf > 0))
  {
    int		f;

    /* For each face of the face buffer. */
    for(f = 0; (errNum == WLZ_ERR_NONE) && (f < wSp->nFceBuf); ++f)
    {
      int	o,
      		oo,
		prvNHorizon;
      WlzConvHullFce *fce;

      fce = wSp->fceBuf[f];
      prvNHorizon = wSp->nHorizon;
      /* For each face connected by an edge to the current face. */
      for(o = 0; o < 3; ++o)
      {
	WlzConvHullArc *oArc;
	WlzConvHullFce *oFce; 

        oFce = fce->opp[o];
	/* Is the opposite face in the face buffer, if not then edge on
	 * the opposite face must form part of the horizon. */
	if((oArc = oFce->arc) == NULL)
	{
	  /* There are no vertices in conflict with this face so it's edge
	   * must form part of the horizon. */
          for(oo = 0; oo < 3; ++oo)
	  {
	    if(oFce->opp[oo] == fce)
	    {
	      wSp->horizon[wSp->nHorizon].fce = oFce;
	      wSp->horizon[wSp->nHorizon].edg = oo;
	      ++(wSp->nHorizon);
	    }
	  }
	}
	else
	{
	  int 	b,
	  	inBuf = 0;

	  /* Look to see if the opposite face is in the face buffer,
	   * if not must if form part of the horizon.. */

	  for(b = 0; b < wSp->nFceBuf; ++b)
	  {
	    if(oFce == wSp->fceBuf[b])
	    {
	      inBuf = 1;
	    }
	  }
	  if(inBuf == 0)
	  {
	    for(oo = 0; oo < 3; ++oo)
	    {
	      if(oFce->opp[oo] == fce)
	      {
		wSp->horizon[wSp->nHorizon].fce = oFce;
		wSp->horizon[wSp->nHorizon].edg = oo;
		++(wSp->nHorizon);
	      }
	    }
	  }
	}
	if(wSp->nHorizon > prvNHorizon)
	{
	  /* Stay 3 horizon edges ahead to keep allocation simple. */
	  errNum = WlzConvHullHorizonExpand(wSp, wSp->nHorizon + 3);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  /* Squeeze out any repeated horizon edges and make them into a chain
   * so that each follows the previous. */
  if((errNum == WLZ_ERR_NONE) && (wSp->nHorizon > 1))
  {
    int		i,
    		vF, 			/* Vertex of first horizon edge, used
					 * to check for closed loop of vertex
					 * edges. */
		cnt,
		loopOK = 0;

    cnt = 1;
    vF = wSp->horizon[0].fce->vtx[wSp->horizon[0].edg];
    for(i = 0; i < wSp->nHorizon; ++i)
    {
      int	j,
      		vC;			/* Vertex of the current edge. */
      WlzConvHullHorEdg	*edgC, 		/* Current horizon edge. */
      			*edgN; 		/* Next horizon edge. */
      
      edgC = wSp->horizon + i;
      vC = edgC->fce->vtx[(edgC->edg + 1) % 3];
      if(vC == vF)
      {
	/* Have a complete loop of vertex edges. */
	loopOK = 1;
        break;
      }
      /* Find next edge and put it in after the current one. */
      if(i < (wSp->nHorizon - 1))
      {
	edgN = edgC + 1;
	for(j = i + 1; j < wSp->nHorizon; ++j)
	{
	  int	vT;
	  WlzConvHullHorEdg *edgT;

	  edgT = wSp->horizon + j;
	  vT = edgT->fce->vtx[edgT->edg];
	  if(vC == vT)
	  {
	    WlzConvHullHorEdg t;

	    ++cnt;
	    /* Found next edge swap it into position in the buffer. */
	    t = *edgN;
	    *edgN = *edgT;
	    *edgT = t;
	    break;
	  }
	}
      }
    }
    if(loopOK && (cnt >= 3))
    {
      wSp->nHorizon = cnt;
    }
    else
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  return(errNum);
}

/*!
* \return	New 3D convex hull domain.
* \ingroup	WlzConvexHull
* \brief	Creates a new 3D convex hull domain from the given
* 		3D convex hull workspace using it's vertex and faces lists
* 		together with it's vertex position pointer.
* \param	wSp			Convex hull workspace.
* \param	dstErr			Desination error pointer, may be NULL.
*/
static WlzConvHullDomain3	*WlzConvexHullFromWSp3(
				  WlzConvHullWSp3 *wSp,
				  WlzErrorNum *dstErr)
{
  int		nFce = 0,
  		nVtx = 0;
  WlzConvHullFce *fce;
  WlzConvHullVtx *vtx;
  WlzConvHullDomain3 *cvh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find number of vertices and faces, building a lookup table (using
   * the workspace permutation table) which maps given vertex array index
   * to the convex hull domain vertex array index. */
  if((vtx = wSp->vtxLst) != NULL)
  {
    int		i = 0;

    do
    {
      wSp->vtxPrm[vtx->idx] = i++;
      vtx = vtx->nxt;
    } while(vtx != wSp->vtxLst);
    nVtx = i;
  }
  if((nVtx > 0) && ((fce = wSp->fceLst) != NULL))
  {
    int		i = 0;

    do
    {
      ++i;
      fce = fce->nxt;
    } while(fce != wSp->fceLst);
    nFce = i;
  }
  if((nVtx <= 0) || (nFce <= 0))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  /* Make convex hull domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    cvh = WlzMakeConvexHullDomain3(nVtx, nFce, wSp->vtxType, &errNum);
  }
  /* Copy the vertices and translate the face indices into the new convex
   * hull domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    int		*face;
    WlzDVertex3 cen;

    i = 0;
    vtx = wSp->vtxLst;
    WLZ_VTX_3_ZERO(cen);
    if(wSp->vtxType == WLZ_VERTEX_I3)
    {
      do
      {
	WlzIVertex3 v;

	v = wSp->vtxPos.i3[vtx->idx];
	WLZ_VTX_3_ADD(cen, cen, v);
	cvh->vertices.i3[i] = v;
	++i;
	vtx = vtx->nxt;
      } while(vtx != wSp->vtxLst);
    }
    else /* wSp->vtxType == WLZ_VERTEX_D3 */
    {
      do
      {
	WlzDVertex3 v;

	v = wSp->vtxPos.d3[vtx->idx];
	WLZ_VTX_3_ADD(cen, cen, v);
	cvh->vertices.d3[i] = v;
	++i;
	vtx = vtx->nxt;
      } while(vtx != wSp->vtxLst);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      fce = wSp->fceLst;
      i = 0;
      do
      {
	int	j;

	face = cvh->faces + i;
	for(j = 0; j < 3; ++j)
	{
	  face[j] = wSp->vtxPrm[fce->vtx[j]];
	}
	i += 3;
      fce = fce->nxt;
      } while(fce != wSp->fceLst);
      WLZ_VTX_3_SCALE(cen, cen, (1.0 / nVtx));
      if(wSp->vtxType == WLZ_VERTEX_I3)
      {
	WLZ_VTX_3_NINT(cvh->centroid.i3, cen);
      }
      else /* wSp->vtxType == WLZ_VERTEX_D3 */
      {
	cvh->centroid.d3 = cen;
      }
      cvh->vtxType = wSp->vtxType;
      cvh->nVertices = nVtx;
      cvh->nFaces = nFce;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}

#ifdef WLZ_CONVHULL_DEBUG_FCE
/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Checks the topology  of the given face.
* 		This function is only intended for use in debugging
* 		it is very compute hungry and quite inefficient.
* \param	wSp			Given workspace.
* \param	fce			Given face.
*/
static WlzErrorNum		WlzConvHullChkFce(
				  WlzConvHullWSp3 *wSp,
				  WlzConvHullFce *fce)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fce == NULL) ||
     (fce->opp[0] == NULL) || (fce->opp[1] == NULL) || (fce->opp[2] == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Check the topology. */
  if(errNum == WLZ_ERR_NONE) 
  {
    int		i;

    for(i = 0; i < 3; ++i)
    {
      int	j,
      		v0,
      		v1;
      WlzConvHullFce *oFce;

      v0 = fce->vtx[i];
      v1 = fce->vtx[(i + 1) % 3];
      oFce = fce->opp[i];
      if((oFce->opp[0] == NULL) ||
         (oFce->opp[1] == NULL) ||
	 (oFce->opp[2] == NULL))
      {
        errNum = WLZ_ERR_PARAM_NULL;
	break;
      }
      for(j = 0; j < 3; ++j)
      {
        if(oFce->vtx[j] == v1)
	{
	  if((oFce->vtx[(j + 1) % 3] != v0) || (oFce->opp[j] != fce))
	  {
	    errNum = WLZ_ERR_PARAM_DATA;
	  }
	  break;
	}
      }
      if(j == 4)
      {
        errNum = WLZ_ERR_PARAM_DATA;
	break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Checks the topology  of all the faces in the face
* 		list of the workspace.
* 		This function is only intended for use in debugging
* 		it is very compute hungry and quite inefficient.
* \param	wSp			Given workspace.
*/
static WlzErrorNum		WlzConvHullChkAll(
				  WlzConvHullWSp3 *wSp)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp && wSp->fceLst && wSp->vtxLst)
  {
    int		cnt = 0;
    const int	cntMax = 1000000;
    WlzConvHullFce *fce;

    fce = wSp->fceLst;
    do
    {
      if((errNum = WlzConvHullChkFce(wSp, fce)) != WLZ_ERR_NONE)
      {
	break;
      }
      else if(cnt > cntMax)
      {
        errNum = WLZ_ERR_PARAM_DATA;
	break;
      }
      fce = fce->nxt;
    } while(fce != wSp->fceLst);
  }
  return(errNum);
}
#endif /* WLZ_CONVHULL_DEBUG_FCE */

/*!
* \return	Woolz error code. If the volume of the tetrahedron is
* 		zero WLZ_ERR_DEGENERATE will be returned.
* \ingroup	WlzConvexHull
* \brief	Sets the vertex permutation indices so that the first
* 		four are those of a tetrahedron with z_max, z_min
* 		then the vertex giving the max area triangle followed
* 		by the vertex giving the max volume tetrahedron.
* \param	nPnt			Number of given vertex positions.
* \param	wSp			Given workspace.
* \param	dstZVol			Destination pointer used to set
* 					for zero volume tetrahedron.
*/
static WlzErrorNum		WlzConvHullInitTet(
				  int nPnt,
				  WlzConvHullWSp3 *wSp)
{
  int		i,
	      	m,
	        n;
  double	a,
  		aMax;
  WlzDVertex3	u0,
  		pos;
  int	      	tetIdx[4], /* Buffer of 4 indices to zMax, zMin, yMin, xMin. */
		oldIdx[4],
		addIdx[4],
		delIdx[4];
  WlzDVertex3   tetPos[4]; /* Buffer of 4 positions, zMax, zMin, yMin, xMin. */
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const double	eps = WLZ_CONVHULL_EPS;

  tetIdx[0] = tetIdx[1] = tetIdx[2] = tetIdx[3] = 0;
  /* Find vertices with max (tetPos[0]) and min (tetPos[1]) z */
  if(wSp->vtxType == WLZ_VERTEX_I3)
  {
    WlzIVertex3 v;

    v = wSp->vtxPos.i3[wSp->vtxPrm[0]];
    WLZ_VTX_3_SET(tetPos[0], v.vtX, v.vtY, v.vtZ);
  }
  else
  {
    tetPos[0] = wSp->vtxPos.d3[wSp->vtxPrm[0]];
  }
  tetPos[1] = tetPos[2] = tetPos[3] = tetPos[0];
  for(i = 1; i < nPnt; ++i)
  {
    if(wSp->vtxType == WLZ_VERTEX_I3)
    {
      WlzIVertex3 v;

      v = wSp->vtxPos.i3[wSp->vtxPrm[i]];
      WLZ_VTX_3_SET(pos, v.vtX, v.vtY, v.vtZ);
    }
    else
    {
      pos = wSp->vtxPos.d3[wSp->vtxPrm[i]];
    }
    if(pos.vtZ > tetPos[0].vtZ)
    {
      tetIdx[0] = i;
      tetPos[0] = pos;
    }
    else if(pos.vtZ < tetPos[1].vtZ)
    {
      tetIdx[1] = i;
      tetPos[1] = pos;
    }
  }
  WLZ_VTX_3_SUB(u0, tetPos[1], tetPos[0]);
  a = WLZ_VTX_3_SQRLEN(u0);
  if(a < eps)
  {
    errNum = WLZ_ERR_DEGENERATE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find vertex creates max area triangle with this line segment. */
    aMax = 0.0;
    for(i = 0; i < nPnt; ++i)
    {
      if((i != tetIdx[0]) && (i != tetIdx[1]))
      {
	WlzDVertex3 u1,
		    u2;

	if(wSp->vtxType == WLZ_VERTEX_I3)
	{
	  WlzIVertex3 v;

	  v = wSp->vtxPos.i3[wSp->vtxPrm[i]];
	  WLZ_VTX_3_SET(pos, v.vtX, v.vtY, v.vtZ);
	}
	else
	{
	  pos = wSp->vtxPos.d3[wSp->vtxPrm[i]];
	}
	WLZ_VTX_3_SUB(u1, pos, tetPos[0]);
	WLZ_VTX_3_CROSS(u2, u0, u1);
	a = WLZ_VTX_3_SQRLEN(u2);
	if(a > aMax)
	{
	  tetIdx[2] = i;
	  tetPos[2] = pos;
	  aMax = a;
	}
      }
    }
    if(aMax < eps)
    {
      errNum = WLZ_ERR_DEGENERATE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Find vertex creates max volume tetrahedron with this triangle. */
    aMax = 0.0;
    for(i = 0; i < nPnt; ++i)
    {
      if((i != tetIdx[0]) && (i != tetIdx[1]) && (i != tetIdx[2]))
      {
	if(wSp->vtxType == WLZ_VERTEX_I3)
	{
	  WlzIVertex3 v;

	  v = wSp->vtxPos.i3[wSp->vtxPrm[i]];
	  WLZ_VTX_3_SET(pos, v.vtX, v.vtY, v.vtZ);
	}
	else
	{
	  pos = wSp->vtxPos.d3[wSp->vtxPrm[i]];
	}
	a = WlzGeomTetraSnVolume6(tetPos[0], tetPos[1], tetPos[2], pos);
	a *= a;
        if(a > aMax)	
	{
	  tetIdx[3] = i;
	  tetPos[3] = pos;
	  aMax = a;
	}
      }
    }
    if(aMax < eps)
    {
      errNum = WLZ_ERR_DEGENERATE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    a = WlzGeomTetraSnVolume6(tetPos[0], tetPos[1], tetPos[2], tetPos[3]);
    if(a < 0)
    {
      m = tetIdx[3]; tetIdx[3] = tetIdx[2]; tetIdx[2] = m;
      pos = tetPos[3]; tetPos[3] = tetPos[2]; tetPos[2] = pos;
    }
    /* Make the first four vertex indices of the permutation buffer the
     * indices of the tetrahedron. */
    for(i = 0; i < 4; ++i)
    {
      oldIdx[i] = wSp->vtxPrm[i];
    }
    for(i = 0; i < 4; ++i)
    {
      int	t;

      t = tetIdx[i];
      if(t < 4)
      {
	wSp->vtxPrm[i] = oldIdx[t];
      }
      else
      {
	wSp->vtxPrm[i] = wSp->vtxPrm[t];
      }
    }
    m = n = 0;
    for(i = 0; i < 4; ++i)
    {
      int	t;

      t = wSp->vtxPrm[i];
      if((t != oldIdx[0]) && (t != oldIdx[1]) &&
	 (t != oldIdx[2]) && (t != oldIdx[3]))
      {
	addIdx[n++] = t;
      }
      t = oldIdx[i];
      if((t != wSp->vtxPrm[0]) && (t != wSp->vtxPrm[1]) &&
	 (t != wSp->vtxPrm[2]) && (t != wSp->vtxPrm[3]))
      {
	delIdx[m++] = t;
      }
    }
    for(i = 4; (n > 0) && (i < nPnt); ++i)
    {
      int	j,
		t;

      t = wSp->vtxPrm[i];
      for(j = n - 1; j >= 0; --j)
      {
	if(t == addIdx[j])
	{
	  wSp->vtxPrm[i] = delIdx[--m];
	  break;
	}
      }
      n = m;
    }
    /* Randomise the indices above the first four of the tetrahderon. */
    n = nPnt - 4;
    srand(0);
    for(i = 0; i < n; ++i)
    {
      int		r,
		  t;

      r = (rand() % n) + 4;
      t = wSp->vtxPrm[r];
      wSp->vtxPrm[r] = wSp->vtxPrm[i + 4];
      wSp->vtxPrm[i + 4] = t;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Build the initial tetrahedron faces, check for vertices inside
* 		the initial tetrahedron, adding vertices to the vertex list if
* 		they're not in the tetrahedron and initialising the conflict
* 		graph with all visible pairs of remaining vertices and faces.
* \param	nPnt			Number of given vertex positions.
* \param	wSp			Given workspace.
*/
static WlzErrorNum		WlzConvHullBuildTet(
				  int nPnt,
				  WlzConvHullWSp3 *wSp)
{
  int		i;
  WlzConvHullVtx *vtx[4];
  WlzConvHullFce **fce;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  fce = wSp->fceBuf;
  for(i = 0; i < 4; ++i)
  {
    fce[i] = WlzConvHullNewFce(wSp, NULL);
    vtx[i] = WlzConvHullNewVtx(wSp, NULL);
    vtx[i]->idx = wSp->vtxPrm[i];
  }
  /* 0 */
  fce[0]->vtx[0] = vtx[0]->idx;
  fce[0]->vtx[1] = vtx[2]->idx;
  fce[0]->vtx[2] = vtx[1]->idx;
  fce[0]->opp[0] = fce[2];
  fce[0]->opp[1] = fce[3];
  fce[0]->opp[2] = fce[1];
  /* 1 */
  fce[1]->vtx[0] = vtx[1]->idx;
  fce[1]->vtx[1] = vtx[3]->idx;
  fce[1]->vtx[2] = vtx[0]->idx;
  fce[1]->opp[0] = fce[3];
  fce[1]->opp[1] = fce[2];
  fce[1]->opp[2] = fce[0];
  /* 2 */
  fce[2]->vtx[0] = vtx[2]->idx;
  fce[2]->vtx[1] = vtx[0]->idx;
  fce[2]->vtx[2] = vtx[3]->idx;
  fce[2]->opp[0] = fce[0];
  fce[2]->opp[1] = fce[1];
  fce[2]->opp[2] = fce[3];
  /* 3 */
  fce[3]->vtx[0] = vtx[3]->idx;
  fce[3]->vtx[1] = vtx[1]->idx;
  fce[3]->vtx[2] = vtx[2]->idx;
  fce[3]->opp[0] = fce[1];
  fce[3]->opp[1] = fce[0];
  fce[3]->opp[2] = fce[2];
  for(i = 0; i < 4; ++i)
  {
    if((errNum = WlzConvHullFceNrm(wSp, fce[i])) != WLZ_ERR_NONE)
    {
      break;
    }
    WlzConHullAddConvexFce(wSp, fce[i]);
    WlzConHullAddConvexVtx(wSp, vtx[i]);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 4; i < nPnt; ++i)
    {
      int	behind = 1;
      WlzConvHullVtx *newVtx;

      newVtx = WlzConvHullNewVtx(wSp, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int		f;

	newVtx->idx = wSp->vtxPrm[i];
	for(f = 0; f < 4; ++f)
	{
	  WlzConvHullFce *qFce;

	  qFce = wSp->fceBuf[f];
	  if(WlzConvHullFceBehind(wSp, wSp->fceBuf[f], newVtx->idx) == 0)
	  {
	    behind = 0;
	  }
	  else
	  {
	    errNum = WlzConvHullConfAdd(wSp, qFce, newVtx);
	    if(errNum != WLZ_ERR_NONE)
	    {
	      i = nPnt; 	    	    /* Break from outer vertex loop. */
	      break;			      /* Break from inner face loop. */
	    }
	  }
	}
	if(behind)
	{
	  /* Remove vertex as behind all faces and conflicts of the current
	   * convex hull. */
	  WlzConvHullDelVtx(wSp, NULL, newVtx);
	}
	else
	{
	  /* Push vertex onto the pending vertex queue. */
	  WlzConHullPushPendingVtx(wSp, newVtx);
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	New 3D convex hull domain.
* \ingroup	WlzConvexHull
* \brief	Creates a new 3D convex hull domain which encloses the
* 		given known degenerate (all vertices on a single plane,
* 		all on a single line or all coincident).
* \param	pType			Type of vertex given, must be either
* 					WLZ_VERTEX_I3 or WLZ_VERTEX_D3.
* \param	nPnt			Number of given vertices.
* \param	pnt			The given vertices.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzConvHullDomain3	*WlzConvexHullDegenerate3(
				  WlzVertexType pType,
				  int nPnt,
				  WlzVertexP pnt,
				  WlzErrorNum *dstErr)
{
  int		i;
  double	t,
		l,
  		l0,
  	        l1,
		l2;
  WlzDVertex3	c,
		u,
  		v0,
		v1,
		v2;
  WlzConvHullDomain3 *cvh = NULL;
  const double	eps = 1.0e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
				  
  /* Find centroid of vertices \f$c\f$. */
  WLZ_VTX_3_ZERO(c);
  for(i = 0; i < nPnt; ++i)
  {
    if(pType == WLZ_VERTEX_I3)
    {
      c.vtX += pnt.i3[i].vtX;
      c.vtY += pnt.i3[i].vtY;
      c.vtZ += pnt.i3[i].vtZ;
    }
    else
    {
      WLZ_VTX_3_ADD(c, c, pnt.d3[i]);
    }
  }
  WLZ_VTX_3_SCALE(c, c, (1.0 / nPnt));
  /* Find vertices closest and furthest from \f$c\f$ the centroid,
   * \f$v_0\f$ and \f$v_1\f$. */
  for(i = 0; i < nPnt; ++i)
  {
    WlzDVertex3 d,
    		v;

    if(pType == WLZ_VERTEX_I3)
    {
      WLZ_VTX_3_SET(v, pnt.i3[i].vtX, pnt.i3[i].vtY, pnt.i3[i].vtZ);
    }
    else
    {
      v = pnt.d3[i];
    }
    WLZ_VTX_3_SUB(d, c, v);
    l = WLZ_VTX_3_SQRLEN(d);
    if(i == 0)
    {
      l0 = l1 = l;
      v0 = v1 = v;
    }
    else if(l < l0)
    {
      l0 = l;
      v0 = v;
    }
    else if(l > l1)
    {
      l1 = l;
      v1 = v;
    }
  }
  /* If length \f$|v_0 - v_1 | < \epsilon\f$ then have all vertices
   * coincident. */
  WLZ_VTX_3_SUB(u, v0, v1);
  l = WLZ_VTX_3_SQRLEN(u);
  if(l < eps)
  {
    /* Set degenerate convex hull with a single vertex and a single
     * degenerate face. */
    cvh = WlzMakeConvexHullDomain3(1, 1, pType, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      cvh->nVertices = 1;
      cvh->nFaces = 1;
      if(pType == WLZ_VERTEX_I3)
      { 
        cvh->vertices.i3[0] = pnt.i3[0];
      }
      else
      {
        cvh->vertices.d3[0] = pnt.d3[0];
      }
      cvh->faces[0] = cvh->faces[1] = cvh->faces[2] = 0;
    }
  }
  else
  {
    double	uu;

    /* Now set \f$v_0\f$ to be the vertex furthest from \f$v_1\f$. */
    l0 = 0.0;
    v0 = v1;
    for(i = 0; i < nPnt; ++i)
    {
      double	l;
      WlzDVertex3 d,
		  v;

      if(pType == WLZ_VERTEX_I3)
      {
	WLZ_VTX_3_SET(v, pnt.i3[i].vtX, pnt.i3[i].vtY, pnt.i3[i].vtZ);
      }
      else
      {
	v = pnt.d3[i];
      }
      WLZ_VTX_3_SUB(d, v1, v);
      l = WLZ_VTX_3_SQRLEN(d);
      if(l > l0)
      {
	l0 = l;
	v0 = v;
      }
    }
    /* Find vertex \f$v_2\f$ furthest from the line segment \f$v_0\f$,
     * \f$v_1\f$.
     * To do this consider the line segment to be parameterised
     * \f$n(t) = v_0 + t u\f$,
     * where \f$u = v_1 - v_0\f$ and \f$t \in \mathbb{R}\f$.
     * At the closest point on the line segment to \f$v_2\f$ we have
     * \f$t = (u . (v_2 - v_0))/(u . u)\f$
     * then the distance of \f$v_2\f$ from the line segment is
     * 
     * \f[
       l = \left{
	   \begin{array}{ll}
	   |v_2 - v_0|         & t < 0 \\
	   |v_2 - (v_0 + t u)| & 0 < t < 1 \\
	   |v_2 - v_1|         & t > 1
	   \end{array}
	   \right.
       \f]
     */ 
    l2 = 0; 		/* These two assignments keep the compiller happy */
    v2 = v0;		/* allowing for (impossible ?) coincident vertices. */
    for(i = 0; i < nPnt; ++i)
    {
      l2 = 0;
      WLZ_VTX_3_SUB(u, v1, v0);
      uu = WLZ_VTX_3_SQRLEN(u);
      for(i = 0; i < nPnt; ++i)
      {
	double      l;
	WlzDVertex3 d,
		    v;

	if(pType == WLZ_VERTEX_I3)
	{
	  WLZ_VTX_3_SET(v, pnt.i3[i].vtX, pnt.i3[i].vtY, pnt.i3[i].vtZ);
	}
	else
	{
	  v = pnt.d3[i];
	}
	WLZ_VTX_3_SUB(d, v, v0);
	t = WLZ_VTX_3_SQRLEN(d);
	if((t > -eps) && (t < (uu + eps)))
	{
	  t = t / uu;
	  WLZ_VTX_3_SCALE_ADD(d, u, t, v0);
	  WLZ_VTX_3_SUB(d, v, d);
	  l = WLZ_VTX_3_SQRLEN(d);
	  if(l > l2)
	  {
	    l2 = l;
	    v2 = v;
	  }
	}
      }
      /* If length \f$v_2\f$ to line segment \f$v_0\f$, \f$v_1\f$ is less than
       * \f$\epsilon\f$ have all vertices on the line segment \f$v_0\f$,
       * \f$v_1\f$. */
      if(l2 < eps)
      {
	/* Set degenerate convex hull with two vertices and a single degenerate
	 * face. */
	cvh = WlzMakeConvexHullDomain3(2, 1, pType, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  cvh->nVertices = 1;
	  cvh->nFaces = 1;
	  if(pType == WLZ_VERTEX_I3)
	  { 
	    WLZ_VTX_3_NINT(cvh->vertices.i3[0], v0);
	    WLZ_VTX_3_NINT(cvh->vertices.i3[1], v1);
	  }
	  else
	  {
	    cvh->vertices.d3[0] = v0;
	    cvh->vertices.d3[1] = v1;
	  }
	  cvh->faces[0] = 0,
	  cvh->faces[1] = 1;
	  cvh->faces[2] = 0;
	}
      }
      else
      {
        WlzVertexP pnt2;

	if((pnt2.v = AlcMalloc(nPnt * sizeof(WlzDVertex2))) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzDVertex3 b0,
		      b1,
		      t0,
		      t1;
	  WlzConvHullDomain2 *cvh2 = NULL;

	  /* The convex hull is planar.
	   * Create a new set of unit basis vectors for the plane, \f$b_0\f$
	   * and \f$b_1\f$
	   * \f[
	     b_0 = (v_1 - v_0) / |v_1 - v_0|
	     b_1 = ((v_1 - v_0) \times ((v_1 - v_0) \times (v_2 - v_0))) /
		   |(v_1 - v_0) \times ((v_1 - v_0) \times (v_2 - v_0))|
	     \f]
	   */
	  WLZ_VTX_3_SUB(b0, v1, v0);
	  l = WLZ_VTX_3_LENGTH(b0);
	  WLZ_VTX_3_SCALE(b0, b0, (1.0/l));
	  WLZ_VTX_3_SUB(t0, v2, v0);
	  WLZ_VTX_3_CROSS(t1, b0, t0);
	  WLZ_VTX_3_CROSS(b1, b0, t1);
	  l = WLZ_VTX_3_LENGTH(b1);
	  WLZ_VTX_3_SCALE(b1, b1, (1.0/l));
	  /* Map the 3D vertices to 2D by taking scalar products with \f$b_0\f$
	   * and \f$b_1\f$, ie for \f$p \in \mathbb{R}^3\f$ and
	   * \f$q \in \mathbb{R}^2\f$ then
	   * \f$q = (b_0 . (p - v_0), b_1 . (p - v_0))\f$,
	   * where
	   * \f$q = (q0, q1)\f$. */
	  for(i = 0; i < nPnt; ++i)
	  {
	    WlzDVertex3 v;

	    if(pType == WLZ_VERTEX_I3)
	    {
	      WLZ_VTX_3_SET(v, pnt.i3[i].vtX, pnt.i3[i].vtY, pnt.i3[i].vtZ);
	    }
	    else
	    {
	      v = pnt.d3[i];
	    }
	    WLZ_VTX_3_SUB(v, v, v0);
	    pnt2.d2[i].vtX = WLZ_VTX_3_DOT(b0, v);
	    pnt2.d2[i].vtY = WLZ_VTX_3_DOT(b1, v);
	  }
	  /* Compute convex hull in 2D. */
	  cvh2 = WlzConvexHullFromVtx2(WLZ_VERTEX_D2, nPnt, pnt2, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    /* Create a 3D convex hull with the given vertex type, the same
	     * number of vertices and faces as the number of vertices in 2D.
	     * Make faces from the 2D vertices and implied edges by using the
	     * centroid as the first vertex and all faces then using it,
	     * cf spokes on a wheel. */
	    cvh = WlzMakeConvexHullDomain3(cvh2->nVertices + 1,
	                                   cvh2->nVertices,
					   pType, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    WlzDVertex2 *p2;
	    WlzDVertex3 c;

	    cvh->nVertices = cvh2->nVertices + 1;
	    cvh->nFaces = cvh2->nVertices;
	    /* Now map the vertices of the 2D convex hull back into the 3D
	     * space using
	     * \f$p = b0 q0 + b1 q1\f$.
	     */ 
	    if(pType == WLZ_VERTEX_I3)
	    {
	      WlzIVertex3 *p3;

	      WLZ_VTX_3_SCALE(t0, b0, cvh2->centroid.d2.vtX);
	      WLZ_VTX_3_SCALE(t1, b1, cvh2->centroid.d2.vtY);
	      WLZ_VTX_3_ADD3(c, t0, t1, v0);
	      WLZ_VTX_3_NINT(cvh->centroid.i3, c);
	      cvh->vertices.i3[0] = cvh->centroid.i3;
	      for(i = 0; i < cvh2->nVertices; ++i)
	      {
		p2 = cvh2->vertices.d2 + i;
		p3 = cvh->vertices.i3 + i + 1;

		WLZ_VTX_3_SCALE(t0, b0, p2->vtX);
		WLZ_VTX_3_SCALE(t1, b1, p2->vtY);
		WLZ_VTX_3_ADD3(t0, t0, t1, v0);
		*p3 = WlzConvexHullVtxD3ToI3(t0, c);
	      }
	    }
	    else
	    {
	      WlzDVertex3 *p3;

	      WLZ_VTX_3_SCALE(t0, b0, cvh2->centroid.d2.vtX);
	      WLZ_VTX_3_SCALE(t1, b1, cvh2->centroid.d2.vtY);
	      WLZ_VTX_3_ADD3(c, t0, t1, v0);
	      cvh->centroid.d3 = c;
	      cvh->vertices.d3[0] = c;
	      for(i = 0; i < cvh2->nVertices; ++i)
	      {
		p2 = cvh2->vertices.d2 + i;
		p3 = cvh->vertices.d3 + i + 1;
		WLZ_VTX_3_SCALE(t0, b0, p2->vtX);
		WLZ_VTX_3_SCALE(t1, b1, p2->vtY);
		WLZ_VTX_3_ADD3(*p3, t0, t1, v0);
	      }
	    }
	    /* Now set the faces. */
	    for(i = 0; i < cvh2->nVertices; ++i)
	    {
	      int *f;

	      f = cvh->faces + (i * 3);
	      f[0] = 0;
	      f[1] = i + 1;
	      f[2] = ((i + 1) % cvh2->nVertices) + 1;
	    }
	  }
	  (void )WlzFreeConvexHullDomain2(cvh2);
	  AlcFree(pnt2.v);
        }
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}

/*!
* \return	New 3D convex hull domain.
* \ingroup	WlzConvexHull
* \brief	Creates a new 3D convex hull domain which encloses the
* 		given vertices using a randomized incremental algorithm
* 		based on Clarkson and Shor's algorithm:
* 		Kenneth L. Clarkson and PeterW. Shor. "Applications of Random
* 		Sampling in Computational Geometry, II", Discrete &
* 		Computational Geometry 4(1) 1989.
* 		This algorithm is also described in a way more suited
* 		for implementation in:
* 		Mark de Berg, Otfried Cheong, Marc van Kreveld and
* 		Mark Overmars. "Computational Geometry: Algorithms and
* 		Applications", Springer.
* 		The algorithm makes use of randomized insertion,
* 		lists, queues and a conflict graph of conflicts between
* 		vertices and faces. A vertex and face are said to be in
* 		conflict when a vertex does not lie behind the face,
* 		where behind means in the half-plane defined by the
* 		face and the body of the convex hull. In practice the
* 		vertex-behind-face test accounts for the almost all
* 		the CPU time. The expected run time for this algorithm
* 		is O(n log n). On a 3GHz Intel i7 the computation time
* 		is a ~10ms for 10^4 vertices randomly distributed over
* 		a cube but ~1s for 10^5 vertices.
* 		When given a degenerate set of vertices (all on a single plane,
* 		all on a single line or all coincident) this function will
* 		still compute a 3D convex hull domain but will set the error
* 		code to WLZ_ERR_DEGENERATE to indicate this.
* \param	pType			Type of vertex given, must be either
* 					WLZ_VERTEX_I3 or WLZ_VERTEX_D3.
* \param	nPnt			Number of given vertices.
* \param	pnt			The given vertices.
* \param	dstErr			Destination error pointer, may be NULL.
* 					If the volume of the tetrahedron
* 					with maximum/minimum z coordinate
* 					and minimum x and y coordiantes is
* 					zero the erro code will be
* 					WLZ_ERR_DEGENERATE.
*/
WlzConvHullDomain3		*WlzConvexHullFromVtx3(
				  WlzVertexType pType,
				  int nPnt,
				  WlzVertexP pnt,
				  WlzErrorNum *dstErr)
{
  WlzConvHullWSp3 *wSp = NULL;
  WlzConvHullDomain3 *cvh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = WLZ_CONVHULL_EPS;

  if(nPnt < 4)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(pnt.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((pType != WLZ_VERTEX_I3) && (pType != WLZ_VERTEX_D3))
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate and initialise datastructures. */
    wSp = WlzConvHullMakeWSp3(pType, nPnt, pnt, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	i,
    		j;

    /* Sort the given vertices using the permutation buffer to remove
     * duplicates. */
    for(i = 0; i < nPnt; ++i)
    {
      wSp->vtxPrm[i] = i;
    }
    if(pType == WLZ_VERTEX_I3)
    {
      (void )AlgHeapSortIdx(pnt.i3, wSp->vtxPrm, nPnt,
                            WlzVertexHeapSortIdxFnI3);
    }
    else /* vtxType == WLZ_VERTEX_D3 */
    {
      (void )AlgHeapSortIdx(pnt.d3, wSp->vtxPrm, nPnt,
                            WlzVertexHeapSortIdxFnD3);
    }
    i = 0;
    j = 0;
    if(pType == WLZ_VERTEX_I3)
    {
      while(j < nPnt)
      {
	WlzIVertex3 v0,
		    v1;

	v0 = pnt.i3[wSp->vtxPrm[i]];
	v1 = pnt.i3[wSp->vtxPrm[j]];
	if((v0.vtX != v1.vtX) ||
	   (v0.vtY != v1.vtY) ||
	   (v0.vtZ != v1.vtZ))
	{
	  ++i;
	}
	wSp->vtxPrm[i] = wSp->vtxPrm[j];
	++j;
      }
    }
    else /* pType == WLZ_VERTEX_D3 */
    {
      while(j < nPnt)
      {
	WlzDVertex3 v0,
		    v1;

	v0 = pnt.d3[wSp->vtxPrm[i]];
	v1 = pnt.d3[wSp->vtxPrm[j]];
	if((fabs(v0.vtX - v1.vtX) > eps) ||
	   (fabs(v0.vtY - v1.vtY) > eps) ||
	   (fabs(v0.vtZ - v1.vtZ) > eps))
	{
	  ++i;
	}
	wSp->vtxPrm[i] = wSp->vtxPrm[j];
	++j;
      }
    }
    nPnt = i + 1;
    if(nPnt < 4)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Make sure there's room for at least the 4 initial vertices and faces
   * in the respective pools. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzConvHullFcePoolExpand(&(wSp->fcePool), 4);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzConvHullVtxPoolExpand(&(wSp->vtxPool), 4);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzConvHullFceBufExpand(wSp, 4);
  }
  /* Find four vertices that form a tetrahedron with the max/min z and min
   * y and x coordinates. Make these the first for in the permutation
   * index table. Randomise the remaining indices. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzConvHullInitTet(nPnt, wSp);
    if(errNum != WLZ_ERR_DEGENERATE)
    {
      /* Build the initial faces, check for vertices inside the initial
       * tetrahedron, adding vertices to the vertex list if they're not
       * in the tetrahedron and initialise the conflict graph with all
       * visible pairs of remaining vertices and faces. */
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzConvHullBuildTet(nPnt, wSp);
      }
      if(errNum == WLZ_ERR_NONE)
      {
#ifdef WLZ_CONVHULL_DEBUG_ITR
	int	itr = 0,
		maxItr = INT_MAX;
#endif /* WLZ_CONVHULL_DEBUG_ITR */
	WlzConvHullVtx *vtx;

	/* For each vertex in the pending vertex queue add it to the convex
	 * hull. */ 
	while((errNum == WLZ_ERR_NONE) &&
#ifdef WLZ_CONVHULL_DEBUG_ITR
	      (itr++ < maxItr) &&
#endif /* WLZ_CONVHULL_DEBUG_ITR */
	      ((vtx = WlzConHullPopPendingVtx(wSp)) != NULL))
	{
	  /* Are thre conflicts?. */
	  if(vtx->arc)
	  {
	    int	f;

	    /* Determine the horizon as an array with those faces/edges that
	     * are to remain and while doing so build a list of faces to be
	     * deleted. */
	    errNum = WlzConvHullBuildHorizon(wSp, vtx);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      /* Delete all faces which conflict with the vertex. */
	      for(f = 0; f < wSp->nFceBuf; ++f)
	      {
		WlzConvHullDelFce(wSp, &(wSp->fceLst), wSp->fceBuf[f]);
	      }
	      wSp->nFceBuf = 0;
	      /* Make sure there's space allocated for a face per horizon edge
	       * (more than enough). */
	      errNum = WlzConvHullFcePoolExpand(&(wSp->fcePool),
	                                        wSp->nHorizon);
	    }
	    /* Make sure there's enough space in the face buffer for the new
	     * faces. */
	    if(errNum == WLZ_ERR_NONE)
	    {
	      wSp->nFceBuf = wSp->nHorizon;
	      errNum = WlzConvHullFceBufExpand(wSp, wSp->nHorizon);
	    }
	    if(errNum == WLZ_ERR_NONE)
	    {
	      /* Create new faces. */
	      for(f = 0; f < wSp->nFceBuf; ++f)
	      {
		wSp->fceBuf[f] = WlzConvHullNewFce(wSp, NULL);
	      }
	      /* Link up new faces. The direction of the horizon edges is
	       * correct for the remaining faces but the opposite of direction
	       * in the new faces. */
	      for(f = 0; f < wSp->nFceBuf; ++f)
	      {
		WlzConvHullFce *fce;
		WlzConvHullHorEdg *edg;

		fce = wSp->fceBuf[f];
		edg = wSp->horizon + f;
		fce->vtx[0] = edg->fce->vtx[(edg->edg + 1) % 3];
		fce->vtx[1] = edg->fce->vtx[edg->edg];
		fce->vtx[2] = vtx->idx;
		fce->opp[0] = edg->fce;
		fce->opp[1] = wSp->fceBuf[(f + wSp->nFceBuf - 1) %
		                          wSp->nFceBuf];
		fce->opp[2] = wSp->fceBuf[(f + 1) % wSp->nFceBuf];
		edg->fce->opp[edg->edg] = fce;
		if((errNum = WlzConvHullFceNrm(wSp, fce)) != WLZ_ERR_NONE)
		{
		  break;
		}
	      }
	    }
	    /* For each of the pending vertices, check for conflicts with
	     * each of the new faces, delete those inside the new faces and
	     * add new conflicts. */
	    if(errNum == WLZ_ERR_NONE)
	    {
	      WlzConvHullVtx *pVtx;

	      if((pVtx = wSp->vtxQue) != NULL)
	      {
		do
		{
		  WlzConvHullVtx *nxtPVtx;

		  nxtPVtx = pVtx->nxt;
		  for(f = 0; f < wSp->nFceBuf; ++f)
		  {
		    WlzConvHullFce *fce;

		    fce = wSp->fceBuf[f];
		    if(WlzConvHullFceBehind(wSp, fce, pVtx->idx) != 0)
		    {
		      errNum = WlzConvHullConfAdd(wSp, fce, pVtx);
		    }
		  }
		  if((pVtx->arc == NULL) && (pVtx->cvx == 0) &&
		     (pVtx != wSp->vtxQue))
		  {
		    WlzConvHullDelVtx(wSp, &(wSp->vtxQue), pVtx);
		  }
		  pVtx = nxtPVtx;
		} while((errNum == WLZ_ERR_NONE) && (pVtx != wSp->vtxQue));
	      }
	    }
#ifdef WLZ_CONVHULL_DEBUG_FCE
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzConvHullChkAll(wSp);
	    }
#endif
	    /* Add new vertex and faces to the convex hull. */
	    if(errNum == WLZ_ERR_NONE)
	    {
	      WlzConHullAddConvexVtx(wSp, vtx);
	      for(f = 0; f < wSp->nFceBuf; ++f)
	      {
		WlzConHullAddConvexFce(wSp, wSp->fceBuf[f]);
	      }
	    }
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	/* Make convex hull domain from the convex hull workspace. */
	cvh = WlzConvexHullFromWSp3(wSp, &errNum);
      }
    }
  }
  /* Free the workspace. */
  WlzConvHullFreeWSp3(wSp);
  /* If the vertices were degenerate try mapping to a plane. */
  if(errNum == WLZ_ERR_DEGENERATE)
  {
    cvh = WlzConvexHullDegenerate3(pType, nPnt, pnt, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WLZ_ERR_DEGENERATE;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}
