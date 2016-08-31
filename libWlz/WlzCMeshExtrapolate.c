#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshExtrapolate_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCMeshExtrapolate.c
* \author       Bill Hill
* \date         July 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Functions to extrapolate values within conforming meshes.
* \ingroup	WlzMesh
*/

#include <limits.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

/*!                                             
* \enum		_WlzCMeshExpFlag
* \ingroup	WlzMesh
* \brief	Bit flags for extrapolation via an expanding front.
* 		Typedef: WlzCMeshExpFlag
*/
typedef enum _WlzCMeshExpFlag
{
  WLZ_CMESHEXP_FLAG_NONE    = (0),
  WLZ_CMESHEXP_FLAG_UNKNOWN = (1),   	/*!< Has unknown value to be
  				             extrapolated. */
  WLZ_CMESHEXP_FLAG_ACTIVE  = (1<<1),   /*!< Active, being used for
                                             extrapolation. If unknown then
					     it's in the queue, if known then
					     it's part of the propagating
					     front and being used for
					     extrapolation. */
  WLZ_CMESHEXP_FLAG_UPDATED = (1<<2)    /*!< Value has been updated through
  					     extrapolation. */
} WlzCMeshExpFlag;

/*!
* \struct	_WlzCMeshExpEnt
* \ingroup	WlzMesh
* \brief	A FIFO queue active node/element entity for mesh value
* 		extrapolation workspace.
* 		Typedef: ::WlzCMeshExpEnt
*/
typedef struct _WlzCMeshExpEnt
{
  struct _WlzCMeshExpEnt *nxt;
  struct _WlzCMeshExpEnt *prv;
  WlzCMeshEntP	ent;

} WlzCMeshExpEnt;

/*!
* \struct	_WlzCMeshExpWSp
* \ingroup	WlzMesh
* \brief	A mesh value extrapolation workspace with a FIFO queue for
* 		active nodes (or elements) during mesh value extrapolation
* 		along with matrices and vectors for SVD.
*		Typedef: ::WlzCMeshExpWSp
*/
typedef struct _WlzCMeshExpWSp
{
  int		nEnt;		  /*!< Number of entities (nodes or elements)
  				       in the mesh. */
  WlzCMeshP	mesh;		  /*!< The conforming mesh. */
  WlzIndexedValues *ixv;	  /*!< Indexed values attached to the mesh. */
  WlzUByte	*flags;		  /*!< Per entity flags. */
  double	*pDst;		  /*!< Approximate propagation distance for
  				       extrapolated entities. */
  int		nKNbr;            /*!< Number of known neighbours. */
  int		nUNbr;            /*!< Number of unknown neighbours. */
  int		maxNbr;		  /*!< Space allocated for neighbour arrays. */
  WlzCMeshEntPP	kNbr;		  /*!< Buffer for known neighbours. */
  WlzCMeshEntPP	uNbr;		  /*!< Buffer for unknown neighbours. */
  double	*nVec;		  /*!< Used for computing normal vector using
  				       SVD and the normal vector itself. */
  AlgMatrix	aMat;		  /*!< Used for computing normal vector using
  				       SVD. */
  AlgMatrix	vMat;		  /*!< Used for computing normal vector using
  				       SVD. */
  WlzCMeshExpEnt *head;           /*!< Head of queue of active entities. */
  WlzCMeshExpEnt *tail;           /*!< Tail of queue of active entities. */
  WlzCMeshExpEnt *pool;		  /*!< Entities available for reuse. */
  AlcBlockStack  *blocks;  	  /*!< Block stack for allocation. */
} WlzCMeshExpWSp;

#define WLZ_CMESH_EXP_MEMINC	(1024)

static WlzCMeshExpWSp		*WlzCMeshExpWSpInit2D(
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  WlzUByte *unk,
				  WlzErrorNum *dstErr);
static WlzCMeshExpWSp		*WlzCMeshExpWSpInit3D(
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  WlzUByte *unk,
				  WlzErrorNum *dstErr);
static WlzCMeshExpEnt 		*WlzCMeshExpWSpGetEnt(
				  WlzCMeshExpWSp *wSp);
static WlzErrorNum		WlzCMeshExpValues2D(
				  WlzObject *gObj,
				  WlzUByte *unk);
static WlzErrorNum		WlzCMeshExpValues3D(
				  WlzObject *gObj,
				  WlzUByte *unk);
static WlzErrorNum		WlzCMeshExpGntVector2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static WlzErrorNum		WlzCMeshExpGntVector3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static WlzErrorNum		WlzCMeshExpUpdate2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static WlzErrorNum		WlzCMeshExpUpdate3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static WlzErrorNum		WlzCMeshExpWSpAddNod(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNodP nod);
static WlzErrorNum		WlzCMeshExpWSpAddNod2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod2D *nod);
static WlzErrorNum		WlzCMeshExpWSpAddNod3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod3D *nod);
static WlzErrorNum		WlzCMeshExpReallocate(
				  WlzCMeshExpWSp *wSp,
				  int dim,
				  int nNbr);
static WlzErrorNum		WlzCMeshExpWSpPoolExpand(
				  WlzCMeshExpWSp *wSp);
static void 			WlzCMeshExpRmDuplicates(
				  int *n,
				  void **ary);
static void			WlzCMeshExpWSpRecycleEnt(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static void			WlzCMeshExpWSpFree(
				  WlzCMeshExpWSp *wSp);
static void			WlzCMeshExpDebugPrintQueue(
				  WlzCMeshExpWSp *wSp);

/*!
* \return	New Woolz object sharing the domain of the given object
* 		but with new values which are all known or NULL on error.
* \ingroup	WlzMesh
* \brief	Given a conforming mesh object with attached values
* 		and an array of know node flags this function extrapolates
* 		the value of the unknown values.
* \param	gObj			Given conforming mesh object.
* \param	unk			Array of unknown value flags which
* 					will be modified and must have space
* 					for at least the maximum index of the
* 					conforming mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzCMeshExpValues(
				  WlzObject *gObj,
				  WlzUByte *unk,
				  WlzErrorNum *dstErr)
{
  WlzObject 	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gObj->type != WLZ_CMESH_2D) && (gObj->type != WLZ_CMESH_3D))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(gObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(gObj->values.core->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if(unk == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    WlzValues	newVal;

    /* Create a new object for return using the domain of the given object,
     * but a new value table copying the values from the given object. */
    newVal = WlzCopyValues(gObj->type, gObj->values, gObj->domain, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      rObj = WlzMakeMain(gObj->type, gObj->domain, newVal, NULL, NULL,
                         &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(gObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WlzCMeshExpValues2D(rObj, unk);
	break;
      case WLZ_CMESH_3D:
	errNum = WlzCMeshExpValues3D(rObj, unk);
	break;
      default:
	break;
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
* \ingroup	WlzMesh
* \brief	Extrapolates values in a conforming 2D mesh object in
* 		place. See WlzCMeshExpValues().
* \param	gObj			Given 2D conforming mesh object.
* \param	unk			Array of unknown value flags.
*/
static WlzErrorNum		WlzCMeshExpValues2D(
				  WlzObject *gObj,
				  WlzUByte *unk)
{
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshExpWSp *wSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if((mesh = gObj->domain.cm2)->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv = gObj->values.x)->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    wSp = WlzCMeshExpWSpInit2D(mesh, ixv, unk, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshExpEnt *ent;

    /* Pull an entity from the active known queue. */
    while((errNum == WLZ_ERR_NONE) &&
          ((ent = WlzCMeshExpWSpGetEnt(wSp)) != NULL))
    {
      /* Compute normal vector at this known entity using it's known
       * neighbours. */
      errNum = WlzCMeshExpGntVector2D(wSp, ent);
      if(errNum == WLZ_ERR_NONE)
      {
        /* Update all unknown neighbouring entities. */
        errNum = WlzCMeshExpUpdate2D(wSp, ent);
      }
      /* Recycle the known entity's queue entry. */
      WlzCMeshExpWSpRecycleEnt(wSp, ent);
    }
  }
  WlzCMeshExpWSpFree(wSp);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Extrapolates values in a conforming 3D mesh object in
* 		place. See WlzCMeshExpValues().
* \param	gObj			Given 3D conforming mesh object.
* \param	unk			Array of unknown value flags.
*/
static WlzErrorNum		WlzCMeshExpValues3D(
				  WlzObject *gObj,
				  WlzUByte *unk)
{
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzCMeshExpWSp *wSp = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if((mesh = gObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv = gObj->values.x)->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    wSp = WlzCMeshExpWSpInit3D(mesh, ixv, unk, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshExpEnt *ent;

    /* Pull an entity from the active known queue. */
    while((errNum == WLZ_ERR_NONE) &&
          ((ent = WlzCMeshExpWSpGetEnt(wSp)) != NULL))
    {
      /* Compute normal vector at this known entity using it's known
       * neighbours. */
      errNum = WlzCMeshExpGntVector3D(wSp, ent);
      if(errNum == WLZ_ERR_NONE)
      {
        /* Update all unknown neighbouring entities. */
        errNum = WlzCMeshExpUpdate3D(wSp, ent);
      }
      /* Recycle the known entity's queue entry. */
      WlzCMeshExpWSpRecycleEnt(wSp, ent);
    }
  }
  WlzCMeshExpWSpFree(wSp);
  return(errNum);
}

/*!
* \return	New mesh extrapolation workspace or NULL on error.
* \ingroup	WlzMesh
* \brief	Allocates and initialises a mesh extrapolation workspace.
* \param	mesh			Given 2D conforming mesh.
* \param	ixv			Values attached to the mesh.
* \param	unk			Flags for unknowns.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshExpWSp		*WlzCMeshExpWSpInit2D(
				  WlzCMesh2D *mesh,
				  WlzIndexedValues *ixv,
				  WlzUByte *unk,
				  WlzErrorNum *dstErr)
{
  WlzCMeshExpWSp *wSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	nVal = 3,
  		initMaxNbr = WLZ_CMESH_EXP_MEMINC;

  if((wSp = (WlzCMeshExpWSp *)AlcCalloc(1, sizeof(WlzCMeshExpWSp))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    wSp->ixv = ixv;
    wSp->mesh.m2 = mesh;
    switch(ixv->vType)
    {
      case WLZ_GREY_INT:    /* FALLTHROUGH */
      case WLZ_GREY_SHORT:  /* FALLTHROUGH */
      case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	break;
      default:
	errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(ixv->attach)
    {
      case WLZ_VALUE_ATTACH_NOD:
	if(ixv->rank > 0)
	{
	  /* Extrapolation of vector and tensor values is not implemented
	   * yet. */
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	}
	else
	{
	  wSp->nEnt = mesh->res.nod.maxEnt;
	}
	break;
      case WLZ_VALUE_ATTACH_ELM:
	/* Extrapolation of element attached values is not implemented yet. */
	errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
	errNum = WLZ_ERR_VALUES_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((wSp->nVec = (double *)AlcMalloc(nVal * sizeof(double))) == NULL) ||
       ((wSp->kNbr.v = AlcMalloc(initMaxNbr *
				 sizeof(WlzCMeshNod2D *))) == NULL) ||
       ((wSp->uNbr.v = AlcMalloc(initMaxNbr *
				 sizeof(WlzCMeshNod2D *))) == NULL) ||
       ((wSp->aMat.rect = AlgMatrixRectNew(initMaxNbr,
       				           nVal, NULL)) == NULL) ||
       ((wSp->vMat.rect = AlgMatrixRectNew(nVal, nVal, NULL)) == NULL) ||
       ((wSp->pDst = (double *)AlcCalloc(wSp->nEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Add all entities which are known but have an unknown neighbour to
   * the active entity queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    wSp->flags = unk;
    wSp->maxNbr = initMaxNbr;
    /* First set the flags using the bit mask WLZ_CMESHEXP_FLAG_UNKNOWN
     * only. */
    for(i = 0; i < wSp->nEnt; ++i)
    {
      unk[i] = (unk[i])? WLZ_CMESHEXP_FLAG_UNKNOWN: 0;
    }
    /* For each unknown node check if it's neighbour is known and not in
     * the queue yet. If so add the neighbour to the active node/element
     * queue. */
    for(i = 0; (errNum == WLZ_ERR_NONE) && (i < wSp->nEnt); ++i)
    {
      if((unk[i] & WLZ_CMESHEXP_FLAG_UNKNOWN) != 0)
      {
	WlzCMeshNod2D *un;
	WlzCMeshEdgU2D *eu0,
		       *eu1;

	un = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, i);
	eu0 = eu1 = un->edu;
	do
	{
	  WlzCMeshNod2D *tn;

	  tn = eu1->next->nod;
	  if((unk[tn->idx] &
	      (WLZ_CMESHEXP_FLAG_UNKNOWN | WLZ_CMESHEXP_FLAG_ACTIVE)) == 0)
	  {
	    errNum = WlzCMeshExpWSpAddNod2D(wSp, tn);
	    unk[tn->idx] |= WLZ_CMESHEXP_FLAG_ACTIVE;
	  }
	  eu1 = eu1->nnxt;
	} while((errNum == WLZ_ERR_NONE) && (eu1 != eu0));
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzCMeshExpWSpFree(wSp);
    wSp = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(wSp);
}

/*!
* \return	New mesh extrapolation workspace or NULL on error.
* \ingroup	WlzMesh
* \brief	Allocates and initialises a mesh extrapolation workspace.
* \param	mesh			Given 3D conforming mesh.
* \param	ixv			Values attached to the mesh.
* \param	unk			Flags for unknowns.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzCMeshExpWSp		*WlzCMeshExpWSpInit3D(
				  WlzCMesh3D *mesh,
				  WlzIndexedValues *ixv,
				  WlzUByte *unk,
				  WlzErrorNum *dstErr)
{
  WlzCMeshExpWSp *wSp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	nVal = 4,
  		initMaxNbr = WLZ_CMESH_EXP_MEMINC;

  if((wSp = (WlzCMeshExpWSp *)AlcCalloc(1, sizeof(WlzCMeshExpWSp))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    wSp->ixv = ixv;
    wSp->mesh.m3 = mesh;
    switch(ixv->vType)
    {
      case WLZ_GREY_INT:    /* FALLTHROUGH */
      case WLZ_GREY_SHORT:  /* FALLTHROUGH */
      case WLZ_GREY_UBYTE:  /* FALLTHROUGH */
      case WLZ_GREY_FLOAT:  /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	break;
      default:
	errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(ixv->attach)
    {
      case WLZ_VALUE_ATTACH_NOD:
	if(ixv->rank > 0)
	{
	  /* Extrapolation of vector and tensor values is not implemented
	   * yet. */
	  errNum = WLZ_ERR_UNIMPLEMENTED;
	}
	else
	{
	  wSp->nEnt = mesh->res.nod.maxEnt;
	}
	break;
      case WLZ_VALUE_ATTACH_ELM:
	/* Extrapolation of element attached values is not implemented yet. */
	errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      default:
	errNum = WLZ_ERR_VALUES_DATA;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((wSp->nVec = (double *)AlcMalloc(nVal * sizeof(double))) == NULL) ||
       ((wSp->kNbr.v = AlcMalloc(initMaxNbr *
				 sizeof(WlzCMeshNod3D *))) == NULL) ||
       ((wSp->uNbr.v = AlcMalloc(initMaxNbr *
				 sizeof(WlzCMeshNod3D *))) == NULL) ||
       ((wSp->aMat.rect = AlgMatrixRectNew(initMaxNbr,
       				           nVal, NULL)) == NULL) ||
       ((wSp->vMat.rect = AlgMatrixRectNew(nVal, nVal, NULL)) == NULL) ||
       ((wSp->pDst = (double *)AlcCalloc(wSp->nEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Add all entities which are known but have an unknown neighbour to
   * the active entity queue. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    wSp->flags = unk;
    wSp->maxNbr = initMaxNbr;
    /* First set the flags using the bit mask WLZ_CMESHEXP_FLAG_UNKNOWN
     * only. */
    for(i = 0; i < wSp->nEnt; ++i)
    {
      unk[i] = (unk[i])? WLZ_CMESHEXP_FLAG_UNKNOWN: 0;
    }
    /* For each unknown node check if it's neighbour is known and not in
     * the queue yet. If so add the neighbour to the active node/element
     * queue. */
    for(i = 0; (errNum == WLZ_ERR_NONE) && (i < wSp->nEnt); ++i)
    {
      if((unk[i] & WLZ_CMESHEXP_FLAG_UNKNOWN) != 0)
      {
	WlzCMeshNod3D *un;
	WlzCMeshEdgU3D *eu0,
		       *eu1;

	un = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, i);
	eu0 = eu1 = un->edu;
	do
	{
	  WlzCMeshNod3D *tn;

	  tn = eu1->next->nod;
	  if((unk[tn->idx] &
	      (WLZ_CMESHEXP_FLAG_UNKNOWN | WLZ_CMESHEXP_FLAG_ACTIVE)) == 0)
	  {
	    errNum = WlzCMeshExpWSpAddNod3D(wSp, tn);
	    unk[tn->idx] |= WLZ_CMESHEXP_FLAG_ACTIVE;
	  }
	  eu1 = eu1->nnxt;
	} while((errNum == WLZ_ERR_NONE) && (eu1 != eu0));
      }
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    WlzCMeshExpWSpFree(wSp);
    wSp = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(wSp);
}

/*!
* \return	Next extrapolation queue entity or NULL either when the
* 		queue is empty or on error.
* \ingroup	WlzMesh
* \brief	Gets the next extrapolation queue entity from the queue.
* \param	wSp			Given extrapolation workspace with
* 					queue.
*/
static WlzCMeshExpEnt 		*WlzCMeshExpWSpGetEnt(
				  WlzCMeshExpWSp *wSp)
{
  WlzCMeshExpEnt *ent = NULL;

  if((ent = wSp->tail) != NULL)
  {
    if((wSp->tail = ent->prv) != NULL)
    {
      wSp->tail->nxt = NULL;
    }
  }
  if(ent == wSp->head)
  {
    wSp->head = NULL;
    wSp->tail = NULL;
  }
  return(ent);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes the mesh extrapolation normal vector at the
* 		given extrapolation queue entity with known value
* 		using it's immediate neighbours.
* \param	wSp			Given extrapolation workspace.
* \param	ent			Given extrapolation queue entity
* 					with known value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum		WlzCMeshExpGntVector2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent)
{
  WlzCMeshNod2D *nod;
  WlzCMeshEdgU2D *eu0,
  		 *eu1;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const int	nVal = 3;

  /* Collect this node and it's neighbours in known and unknown arrays. */
  wSp->nKNbr = 1;
  wSp->nUNbr = 0;
  nod = ent->ent.n2;
  wSp->kNbr.n2[0] = nod;
  eu0 = eu1 = nod->edu;
  do
  {
    int		maxNNbr;

    maxNNbr = ALG_MAX(wSp->nKNbr, wSp->nUNbr);
    if(maxNNbr >= wSp->maxNbr)
    {
      errNum = WlzCMeshExpReallocate(wSp, 2, maxNNbr);
      break;
    }
    nod = eu1->next->nod;
    if((wSp->flags[nod->idx] & WLZ_CMESHEXP_FLAG_UNKNOWN) == 0)
    {
      /* This neighbouring node is known. */
      wSp->kNbr.n2[wSp->nKNbr] = nod;
      ++(wSp->nKNbr);
    }
    else
    {
      /* This neighbouring node is unknown. */
      wSp->uNbr.n2[wSp->nUNbr] = nod;
      ++(wSp->nUNbr);
    }
    eu1 = eu1->nnxt;
  } while(eu1 != eu0);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Squeeze out duplicates in the known and unknown arrays. */
    WlzCMeshExpRmDuplicates(&(wSp->nKNbr), wSp->kNbr.v);
    WlzCMeshExpRmDuplicates(&(wSp->nUNbr), wSp->uNbr.v);
    /* Compute gradient vector. */
    wSp->nVec[0] = wSp->nVec[1] = wSp->nVec[2] = 0.0;
    switch(wSp->nKNbr)
    {
      case 1:
	/* Only this node so assume constant. */
        break;
      case 2:
	/* Only two nodes known so assume simple linear gradient across
	 * the nodes. */
	{
	  double    len;
          double    buf[4];
	  WlzGreyP  gP0,
		    gP1,
		    gP2;
	  const double eps = 1.0e-06;

	  buf[0] = wSp->kNbr.n2[1]->pos.vtX - wSp->kNbr.n2[0]->pos.vtX;
	  buf[1] = wSp->kNbr.n2[1]->pos.vtY - wSp->kNbr.n2[0]->pos.vtY;
	  gP0.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n2[0]->idx);
	  gP1.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n2[1]->idx);
	  gP2.dbp = &(buf[2]);
	  WlzValueCopyGreyToGrey(gP2, 0, WLZ_GREY_DOUBLE,
	                         gP0, 0, wSp->ixv->vType, 1);
	  WlzValueCopyGreyToGrey(gP2, 1, WLZ_GREY_DOUBLE,
	                         gP1, 0, wSp->ixv->vType, 1);
	  buf[2] = buf[3] - buf[2];
	  len = sqrt((buf[0] * buf[0]) +
	             (buf[1] * buf[1]));
	  if(len > eps)
	  {
	    wSp->nVec[0] = buf[0] / len;
	    wSp->nVec[1] = buf[1] / len;
	    wSp->nVec[2] = buf[2] / len;
	  }
	}
	break;
      default:
	/* Three or more nodes known so compute the gradient vector
	 * using SVD by fitting a plane. */
	if(errNum == WLZ_ERR_NONE)
	{
	  int	 i,
		 k;
	  double *row;

	  wSp->aMat.rect->nR = wSp->nKNbr;
	  /* Collect positions and values. */
	  for(k = 0; k < wSp->nKNbr; ++k)
	  {
	    WlzGreyP	gP0,
	    		gP1;
	    WlzCMeshNod2D *nod;

	    nod = wSp->kNbr.n2[k];
	    row = wSp->aMat.rect->array[k];
	    row[0] = nod->pos.vtX;
	    row[1] = nod->pos.vtY;
	    gP1.dbp = &(row[2]);
	    gP0.v = WlzIndexedValueGet(wSp->ixv, nod->idx);
	    WlzValueCopyGreyToGrey(gP1, 0, WLZ_GREY_DOUBLE,
	                           gP0, 0, wSp->ixv->vType, 1);
	  }
	  /* Remove centroid. */
	  for(i = 0; i < nVal; ++i)
	  {
	    double	cen = 0.0;

	    for(k = 0; k < wSp->nKNbr; ++k)
	    {
	      cen += wSp->aMat.rect->array[k][i];
	    }
	    cen /= wSp->nKNbr;
	    for(k = 0; k < wSp->nKNbr; ++k)
	    {
	      wSp->aMat.rect->array[k][i] -= cen;
	    }
	  }
	  /* Compute SVD. */
	  errNum = WlzErrorFromAlg(
		   AlgMatrixSVDecomp(wSp->aMat, wSp->nVec, wSp->vMat));
	}
	/* Find minimum singular value, this corresponds to the in (hyper)plane
	 * normal vector. */
	if(errNum == WLZ_ERR_NONE)
	{
	  int	i,
	        m;

	  m = 0;
	  for(i = 1; i < nVal; ++i)
	  {
	    if(wSp->nVec[i] < wSp->nVec[m])
	    {
	      m = i;
	    }
	  }
	  for(i = 0; i < nVal; ++i)
	  {
	    wSp->nVec[i] = wSp->vMat.rect->array[i][m];
	  }
	}
	break;
    }
  }
  return(errNum);
}


/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes the mesh extrapolation normal vector at the
* 		given extrapolation queue entity with known value
* 		using it's immediate neighbours.
* \param	wSp			Given extrapolation workspace.
* \param	ent			Given extrapolation queue entity
* 					with known value.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzErrorNum		WlzCMeshExpGntVector3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent)
{
  WlzCMeshNod3D *nod;
  WlzCMeshEdgU3D *eu0,
  		 *eu1;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const int	nVal = 4;

  /* Collect this node and it's neighbours in known and unknown arrays. */
  wSp->nKNbr = 1;
  wSp->nUNbr = 0;
  nod = ent->ent.n3;
  wSp->kNbr.n3[0] = nod;
  eu0 = eu1 = nod->edu;
  do
  {
    int		maxNNbr;

    maxNNbr = ALG_MAX(wSp->nKNbr, wSp->nUNbr);
    if(maxNNbr >= wSp->maxNbr)
    {
      errNum = WlzCMeshExpReallocate(wSp, 3, maxNNbr);
      break;
    }
    nod = eu1->next->nod;
    if((wSp->flags[nod->idx] & WLZ_CMESHEXP_FLAG_UNKNOWN) == 0)
    {
      /* This neighbouring node is known. */
      wSp->kNbr.n3[wSp->nKNbr] = nod;
      ++(wSp->nKNbr);
    }
    else
    {
      /* This neighbouring node is unknown. */
      wSp->uNbr.n3[wSp->nUNbr] = nod;
      ++(wSp->nUNbr);
    }
    eu1 = eu1->nnxt;
  } while(eu1 != eu0);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Squeeze out duplicates in the known and unknown arrays. */
    WlzCMeshExpRmDuplicates(&(wSp->nKNbr), wSp->kNbr.v);
    WlzCMeshExpRmDuplicates(&(wSp->nUNbr), wSp->uNbr.v);
    /* Compute gradient vector. */
    switch(wSp->nKNbr)
    {
      case 1:
	/* Only this node so assume constant. */
	wSp->nVec[0] = wSp->nVec[1] = wSp->nVec[2] = wSp->nVec[3] = 0.0;
	break;
      case 2:
	/* Only two nodes known so assume simple linear gradient across
	 * the nodes. */
	{
	  double    len;
          double    buf[5];
	  WlzGreyP  gP0,
		    gP1,
		    gP2;
	  const double eps = 1.0e-06;

	  buf[0] = wSp->kNbr.n3[1]->pos.vtX - wSp->kNbr.n3[0]->pos.vtX;
	  buf[1] = wSp->kNbr.n3[1]->pos.vtY - wSp->kNbr.n3[0]->pos.vtY;
	  buf[2] = wSp->kNbr.n3[1]->pos.vtY - wSp->kNbr.n3[0]->pos.vtY;
	  gP0.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n3[0]->idx);
	  gP1.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n3[1]->idx);
	  gP2.dbp = &(buf[3]);
	  WlzValueCopyGreyToGrey(gP2, 0, WLZ_GREY_DOUBLE,
	                         gP0, 0, wSp->ixv->vType, 1);
	  WlzValueCopyGreyToGrey(gP2, 1, WLZ_GREY_DOUBLE,
	                         gP1, 0, wSp->ixv->vType, 1);
	  buf[3] = buf[4] - buf[3];
	  len = sqrt((buf[0] * buf[0]) +
	             (buf[1] * buf[1]) +
		     (buf[2] * buf[2]));
	  if(len > eps)
	  {
	    wSp->nVec[0] = buf[0] / len;
	    wSp->nVec[1] = buf[1] / len;
	    wSp->nVec[2] = buf[2] / len;
	    wSp->nVec[3] = buf[3] / len;
	  }
	}
	break;
      case 3:
        /* Find the two vectors directed from the current node and
	 * compute their mean. */
	{
	  int	    tst[2];
	  double    len[2];
          double    buf[9];
	  WlzGreyP  gP0,
		    gP1,
		    gP2,
		    gP3;
	  const double eps = 1.0e-06;

	  buf[0] = wSp->kNbr.n3[1]->pos.vtX - wSp->kNbr.n3[0]->pos.vtX;
	  buf[1] = wSp->kNbr.n3[1]->pos.vtY - wSp->kNbr.n3[0]->pos.vtY;
	  buf[2] = wSp->kNbr.n3[1]->pos.vtZ - wSp->kNbr.n3[0]->pos.vtZ;
	  buf[3] = wSp->kNbr.n3[2]->pos.vtX - wSp->kNbr.n3[0]->pos.vtX;
	  buf[4] = wSp->kNbr.n3[2]->pos.vtY - wSp->kNbr.n3[0]->pos.vtY;
	  buf[5] = wSp->kNbr.n3[2]->pos.vtZ - wSp->kNbr.n3[0]->pos.vtZ;
	  gP0.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n3[0]->idx);
	  gP1.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n3[1]->idx);
	  gP2.v = WlzIndexedValueGet(wSp->ixv, wSp->kNbr.n3[2]->idx);
	  gP3.dbp = &(buf[6]);
	  WlzValueCopyGreyToGrey(gP3, 0, WLZ_GREY_DOUBLE,
	                         gP0, 0, wSp->ixv->vType, 1);
	  WlzValueCopyGreyToGrey(gP3, 1, WLZ_GREY_DOUBLE,
	                         gP1, 0, wSp->ixv->vType, 1);
	  WlzValueCopyGreyToGrey(gP3, 2, WLZ_GREY_DOUBLE,
	                         gP2, 0, wSp->ixv->vType, 1);
	  buf[7] = buf[7] - buf[6];
	  buf[8] = buf[8] - buf[6];
	  len[0] = sqrt((buf[0] * buf[0]) +
	                (buf[1] * buf[1]) +
		        (buf[2] * buf[2]));
	  len[1] = sqrt((buf[3] * buf[3]) +
	                (buf[4] * buf[4]) +
		        (buf[5] * buf[5]));
	  tst[0] = len[0] > eps;
	  tst[1] = len[1] > eps;
	  if(tst[0] && tst[1])
	  {
	    len[0] *= 2.0;
	    len[1] *= 2.0;
	    wSp->nVec[0] = (buf[0] / len[0]) + (buf[3] / len[1]);
	    wSp->nVec[1] = (buf[1] / len[0]) + (buf[4] / len[1]);
	    wSp->nVec[2] = (buf[2] / len[0]) + (buf[5] / len[1]);
	    wSp->nVec[3] = (buf[7] / len[0]) + (buf[8] / len[1]);
	  }
	  else if(tst[0])
	  {
	    wSp->nVec[0] = buf[0] / len[0];
	    wSp->nVec[1] = buf[1] / len[0];
	    wSp->nVec[2] = buf[2] / len[0];
	    wSp->nVec[3] = buf[7] / len[0];
	  }
	  else if(tst[1])
	  {
	    wSp->nVec[0] = buf[3] / len[1];
	    wSp->nVec[1] = buf[4] / len[1];
	    wSp->nVec[2] = buf[5] / len[1];
	    wSp->nVec[3] = buf[8] / len[1];
	  }
	}
	break;
      default:
	/* Four or more nodes known so compute the gradient vector
	 * using SVD by fitting a hyperplane. */
	if(errNum == WLZ_ERR_NONE)
	{
	  int	    i,
		    k;
	  double    *row;

	  wSp->aMat.rect->nR = wSp->nKNbr;
	  /* Collect positions and values. */
	  for(k = 0; k < wSp->nKNbr; ++k)
	  {
	    WlzGreyP	gP0,
	    		gP1;
	    WlzCMeshNod3D *nod;

	    nod = wSp->kNbr.n3[k];
	    row = wSp->aMat.rect->array[k];
	    row[0] = nod->pos.vtX;
	    row[1] = nod->pos.vtY;
	    row[2] = nod->pos.vtZ;
	    gP1.dbp = &(row[3]);
	    gP0.v = WlzIndexedValueGet(wSp->ixv, nod->idx);
	    WlzValueCopyGreyToGrey(gP1, 0, WLZ_GREY_DOUBLE,
	                           gP0, 0, wSp->ixv->vType, 1);
	  }
	  /* Remove centroid. */
	  for(i = 0; i < nVal; ++i)
	  {
	    double cen = 0.0;

	    for(k = 0; k < wSp->nKNbr; ++k)
	    {
	      cen += wSp->aMat.rect->array[k][i];
	    }
	    cen /= wSp->nKNbr;
	    for(k = 0; k < wSp->nKNbr; ++k)
	    {
	      wSp->aMat.rect->array[k][i] -= cen;
	    }
	  }
	  /* Compute SVD. */
	  errNum = WlzErrorFromAlg(
		   AlgMatrixSVDecomp(wSp->aMat, wSp->nVec, wSp->vMat));
	}
	/* Find minimum singular value, this corresponds to the in hyperplane
	 * normal vector. */
	if(errNum == WLZ_ERR_NONE)
	{
	  int	i,
	 	m;

	  m = 0;
	  for(i = 1; i < nVal; ++i)
	  {
	    if(wSp->nVec[i] < wSp->nVec[m])
	    {
	      m = i;
	    }
	  }
	  for(i = 0; i < nVal; ++i)
	  {
	    wSp->nVec[i] = wSp->vMat.rect->array[i][m];
	  }
	}
        break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Updates the unknown entities (which are first order neighbours
* 		of the given mesh entity) using the normal vector and if
* 		these are now known adds them to the queue.
* \param	wSp			Given mesh extrapolation workspace.
* \param	ent			Given active extrapolation queue entity
* 					with known value.
* \param	gV			Gradient vector.
*/
static WlzErrorNum		WlzCMeshExpUpdate2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent)
{
  int		idU;
  WlzCMeshNod2D *gnod;
  WlzGreyP	gval;
  const double	eps = 1.0e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gnod = ent->ent.n2;                        /* Given node with known value. */
  gval.v = WlzIndexedValueGet(wSp->ixv, gnod->idx);
  /* For each of the unknown neighbours of the given known entity. */
  for(idU = 0; idU < wSp->nUNbr; ++ idU)
  {
    int		   fstKn,
    		   lstKn;
    double	   gv,
		   uv,
		   pDst;
    WlzDVertex2    del;
    WlzCMeshEdgU2D *eu0,
		   *eu1;
    WlzCMeshNod2D  *uNod;
    WlzGreyP	   tval,
		   uval;

    uNod = wSp->uNbr.n2[idU];
    /* Compute new value and if appropriate set value. 
     * There normal vector has been computed using AlgMatrixSVDecomp()
     * and so is known to be a unit length vector. */
    WLZ_VTX_2_SUB(del, uNod->pos, gnod->pos);
    uval.v = WlzIndexedValueGet(wSp->ixv, uNod->idx);
    tval.dbp = &gv;
    WlzValueCopyGreyToGrey(tval, 0, WLZ_GREY_DOUBLE,
			   gval, 0, wSp->ixv->vType, 1);
    tval.dbp = &uv;
    WlzValueCopyGreyToGrey(tval, 0, WLZ_GREY_DOUBLE,
			   uval, 0, wSp->ixv->vType, 1);
    pDst = wSp->pDst[gnod->idx] + WLZ_VTX_2_LENGTH(del);
    if(fabs(wSp->nVec[2]) < eps)
    {
      if((wSp->flags[uNod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0)
      {
	uv = gv;
	wSp->pDst[uNod->idx] = pDst;
      }
    }
    else
    {
      double v;

      if(wSp->nKNbr == 2)
      {
	v = gv + wSp->nVec[2] * ((del.vtX * wSp->nVec[0]) +
				 (del.vtY * wSp->nVec[1]));
      }
      else /* wSp->nKNbr >= 3 */
      {
	v = gv - 
	    (((wSp->nVec[0] * del.vtX) +
	      (wSp->nVec[1] * del.vtY)) * wSp->nVec[2]);
      }
      if(((wSp->flags[uNod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0) ||
	 (pDst < wSp->pDst[uNod->idx]))
      {
	uv = v;
	wSp->pDst[uNod->idx] = pDst;
	/* Only set the updated flag if the extrapolation is good. */
	if(wSp->nKNbr >= 3)
	{
	  wSp->flags[uNod->idx] |= WLZ_CMESHEXP_FLAG_UPDATED;
	}
      }
    }
    tval.dbp = &uv;
    WlzValueClampGreyIntoGrey(uval, 0, wSp->ixv->vType,
			      tval, 0, WLZ_GREY_DOUBLE, 1);
    /* If the unknown node has just one active known neighbour (ie this
     * given node) then make the unknown node known and active. */
    fstKn = -1;
    lstKn = 1;
    eu0 = eu1 = uNod->edu;
    do
    {
      WlzCMeshNod2D *tnod;

      tnod = eu1->next->nod;
      if((wSp->flags[tnod->idx] & 
	  (WLZ_CMESHEXP_FLAG_UNKNOWN | WLZ_CMESHEXP_FLAG_ACTIVE)) == 
	  WLZ_CMESHEXP_FLAG_ACTIVE)
      {
	if(fstKn == -1)
	{
	  fstKn = tnod->idx;
        }
	else if(fstKn != tnod->idx)
	{
	  lstKn = 0;
	  break;
	}
      }
      eu1 = eu1->nnxt;
    } while(eu1 != eu0);
    if(lstKn)
    {
      wSp->flags[uNod->idx] = WLZ_CMESHEXP_FLAG_ACTIVE;
      errNum = WlzCMeshExpWSpAddNod2D(wSp, uNod);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  wSp->flags[gnod->idx] = 0;      /* The given node is now no longer active. */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Updates the unknown entities (which are first order neighbours
* 		of the given mesh entity) using the normal vector and if
* 		these are now known adds them to the queue.
* \param	wSp			Given mesh extrapolation workspace.
* \param	ent			Given active extrapolation queue entity
* 					with known value.
* \param	gV			Gradient vector.
*/
static WlzErrorNum		WlzCMeshExpUpdate3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent)
{
  int		idU;
  WlzCMeshNod3D *gnod;
  WlzGreyP	gval;
  const double	eps = 1.0e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gnod = ent->ent.n3;                        /* Given node with known value. */
  gval.v = WlzIndexedValueGet(wSp->ixv, gnod->idx);
  /* For each of the unknown neighbours of the given known entity. */
  for(idU = 0; idU < wSp->nUNbr; ++ idU)
  {
    int		   fstKn,
    		   lstKn;
    double	   gv,
		   uv,
		   pDst;
    WlzDVertex3    del;
    WlzCMeshEdgU3D *eu0,
		   *eu1;
    WlzCMeshNod3D  *uNod;
    WlzGreyP	   tval,
		   uval;

    uNod = wSp->uNbr.n3[idU];
    /* Compute new value and if appropriate set value. 
     * There normal vector has been computed using AlgMatrixSVDecomp()
     * and so is known to be a unit length vector. */
    WLZ_VTX_3_SUB(del, uNod->pos, gnod->pos);
    uval.v = WlzIndexedValueGet(wSp->ixv, uNod->idx);
    tval.dbp = &gv;
    WlzValueCopyGreyToGrey(tval, 0, WLZ_GREY_DOUBLE,
			   gval, 0, wSp->ixv->vType, 1);
    tval.dbp = &uv;
    WlzValueCopyGreyToGrey(tval, 0, WLZ_GREY_DOUBLE,
			   uval, 0, wSp->ixv->vType, 1);
    pDst = wSp->pDst[gnod->idx] + WLZ_VTX_2_LENGTH(del);
    if(fabs(wSp->nVec[2]) < eps)
    {
      if((wSp->flags[uNod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0)
      {
	uv = gv;
	wSp->pDst[uNod->idx] = pDst;
      }
    }
    else
    {
      double v;

      if((wSp->nKNbr == 2) || (wSp->nKNbr == 3))
      {
	v = gv + wSp->nVec[2] * ((del.vtX * wSp->nVec[0]) +
				 (del.vtY * wSp->nVec[1]) +
				 (del.vtZ * wSp->nVec[2]));
      }
      else /* wSp->nKNbr >= 4 */
      {
	v = gv - 
	    (((wSp->nVec[0] * del.vtX) +
	      (wSp->nVec[1] * del.vtY) +
	      (wSp->nVec[2] * del.vtZ)) * wSp->nVec[3]);
      }
      if(((wSp->flags[uNod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0) ||
	 (pDst < wSp->pDst[uNod->idx]))
      {
	uv = v;
	wSp->pDst[uNod->idx] = pDst;
	/* Only set the updated flag if the extrapolation is good. */
	if(wSp->nKNbr >= 4)
	{
	  wSp->flags[uNod->idx] |= WLZ_CMESHEXP_FLAG_UPDATED;
	}
      }
    }
    tval.dbp = &uv;
    WlzValueClampGreyIntoGrey(uval, 0, wSp->ixv->vType,
			      tval, 0, WLZ_GREY_DOUBLE, 1);
    /* If the unknown node has just one active known neighbour (ie this
     * given node) then make the unknown node known and active. */
    fstKn = -1;
    lstKn = 1;
    eu0 = eu1 = uNod->edu;
    do
    {
      WlzCMeshNod3D *tnod;

      tnod = eu1->next->nod;
      if((wSp->flags[tnod->idx] & 
	  (WLZ_CMESHEXP_FLAG_UNKNOWN | WLZ_CMESHEXP_FLAG_ACTIVE)) == 
	  WLZ_CMESHEXP_FLAG_ACTIVE)
      {
	if(fstKn == -1)
	{
	  fstKn = tnod->idx;
        }
	else if(fstKn != tnod->idx)
	{
	  lstKn = 0;
	  break;
	}
      }
      eu1 = eu1->nnxt;
    } while(eu1 != eu0);
    if(lstKn)
    {
      wSp->flags[uNod->idx] = WLZ_CMESHEXP_FLAG_ACTIVE;
      errNum = WlzCMeshExpWSpAddNod3D(wSp, uNod);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      break;
    }
  }
  wSp->flags[gnod->idx] = 0;      /* The given node is now no longer active. */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the extrapolation queue. The given
* 		node is known not to be already in the quque.
* \param	wSp			Given mesh extrapolation workspace.
* \param	nod			Given 2 or 3D node.
*/
static WlzErrorNum		WlzCMeshExpWSpAddNod2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod2D *nod)
{
  WlzCMeshNodP	nodP;
  WlzErrorNum	errNum;

  nodP.n2 = nod;
  errNum = WlzCMeshExpWSpAddNod(wSp, nodP);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the extrapolation queue. The given
* 		node is known not to be already in the quque.
* \param	wSp			Given mesh extrapolation workspace.
* \param	nod			Given 3D node.
*/
static WlzErrorNum		WlzCMeshExpWSpAddNod3D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod3D *nod)
{
  WlzCMeshNodP	nodP;
  WlzErrorNum	errNum;

  nodP.n3 = nod;
  errNum = WlzCMeshExpWSpAddNod(wSp, nodP);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the extrapolation queue. The given
* 		node is known not to be already in the quque.
* \param	wSp			Given mesh extrapolation workspace.
* \param	nod			Given 2 or 3D node.
*/
static WlzErrorNum		WlzCMeshExpWSpAddNod(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNodP nod)
{
  WlzCMeshExpEnt *nEnt;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(wSp->pool == NULL)
  {
    errNum = WlzCMeshExpWSpPoolExpand(wSp);
  }
  /* Take an entry from the pool and set it up as the new entry. */
  if(errNum == WLZ_ERR_NONE)
  {
    nEnt = wSp->pool;
    wSp->pool = nEnt->nxt;
    switch(wSp->mesh.core->type)
    {
      case WLZ_CMESH_2D:
        nEnt->ent.n2 = nod.n2;
	break;
      case WLZ_CMESH_3D:
        nEnt->ent.n3 = nod.n3;
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  /* Put the new entry into the queue so that the previous queue entry is
   * either NULL or has a propagation distance <= the new entry. */
  if(errNum == WLZ_ERR_NONE);
  {
    if(wSp->head == NULL)
    {
      wSp->head = nEnt;
      wSp->tail = nEnt;
    }
    else
    {
      double	nDst;
      WlzCMeshExpEnt *pEnt;

      pEnt = wSp->tail;
      switch(wSp->mesh.core->type)
      {
        case WLZ_CMESH_2D:
	  nDst = wSp->pDst[nEnt->ent.n2->idx];
	  while(pEnt && (wSp->pDst[pEnt->ent.n2->idx] < nDst))
	  {
	    pEnt = pEnt->prv;
	  }
	  break;
        case WLZ_CMESH_3D:
	  nDst = wSp->pDst[nEnt->ent.n3->idx];
	  while(pEnt && (wSp->pDst[pEnt->ent.n3->idx] < nDst))
	  {
	    pEnt = pEnt->prv;
	  }
	  break;
        default:
	  break;
      }
      if(pEnt == NULL)
      {
	nEnt->nxt = wSp->head;
	nEnt->nxt->prv = nEnt;
	wSp->head = nEnt;
      }
      else
      {
	nEnt->nxt = pEnt->nxt;
	nEnt->prv = pEnt;
	if(nEnt->nxt)
	{
	  nEnt->nxt->prv = nEnt;
	}
	else
	{
	  wSp->tail = nEnt;
	}
	pEnt->nxt = nEnt;
	if(nEnt->prv == NULL)
	{
	  wSp->head = nEnt;
	}
      }
    }
    wSp->head->prv = NULL;
    wSp->tail->nxt = NULL;
  } 
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Reallocates the 2D neighbour buffer and SVD matrix.
* \param	wSp			Given mesh extrapolation workspace.
* \param	dim			Dimension which must be either 2 or 3.
* \param	nNbr			New number of neighbours.
*/
static WlzErrorNum		WlzCMeshExpReallocate(
				  WlzCMeshExpWSp *wSp,
				  int dim, int nNbr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nNbr >= wSp->maxNbr)
  {
    if((dim != 2) && (dim != 3))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else
    {
      int	n,
		nVal,
		maxNbr;
      const int inc = WLZ_CMESH_EXP_MEMINC;

      nVal = dim + 1;
      n = (nNbr + inc - 1) / inc;
      maxNbr = n * inc;
      if(((wSp->kNbr.v = (void *)AlcRealloc(wSp->kNbr.v,
			     maxNbr * (sizeof(WlzCMeshEntP)))) == NULL) ||
         ((wSp->uNbr.v = (void *)AlcRealloc(wSp->uNbr.v,
			     maxNbr * (sizeof(WlzCMeshEntP)))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
	AlgMatrixRectFree(wSp->aMat.rect);
	if((wSp->aMat.rect = AlgMatrixRectNew(wSp->maxNbr,
		nVal, NULL)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  wSp->maxNbr = maxNbr;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Expands the extrapolation queue pool.
* \param	wSp			Given mesh extrapolation workspace.
*/
static WlzErrorNum		WlzCMeshExpWSpPoolExpand(
				  WlzCMeshExpWSp *wSp)
{
  AlcBlockStack *newBlk;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	inc = WLZ_CMESH_EXP_MEMINC;

  if((newBlk = AlcBlockStackNew(inc, sizeof(WlzCMeshExpEnt),
				wSp->blocks, NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    int		i;
    WlzCMeshExpEnt *ent;

    ent = (WlzCMeshExpEnt *)(newBlk->elements);
    wSp->blocks = newBlk;
    for(i = 0; i < inc; ++i)
    {
      ent->nxt = wSp->pool;
      wSp->pool = ent;
      ++ent;
    }
  }
  return(errNum);
}

/*!
* \ingroup	WlzMesh
* \brief	Recycles the given mesh extrapolation queue entity for reuse.
* \param	wSp			Given mesh extrapolation workspace.
* \param	ent			Given extrapolation queue entity.
*/
static void			WlzCMeshExpWSpRecycleEnt(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent)
{
  ent->nxt = wSp->pool;
  wSp->pool = ent;
}

/*!
* \ingroup	WlzMesh
* \brief	Frees the given mesh extrapolation workspace.
* \param	wSp			Given mesh extrapolation workspace.
*/
static void			WlzCMeshExpWSpFree(
				  WlzCMeshExpWSp *wSp)
{
  if(wSp)
  {
    AlcFree(wSp->nVec);
    AlcFree(wSp->pDst);
    AlcFree(wSp->kNbr.v);
    AlcFree(wSp->uNbr.v);
    AlcBlockStackFree(wSp->blocks);
    AlgMatrixRectFree(wSp->aMat.rect);
    AlgMatrixRectFree(wSp->vMat.rect);
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Squeezes out duplicates from the given array of pointers.
* \param	nAry			Number of elements in the array on
* 					entry and return.
* \param	ary			The array of pointers.
*/
static void 			WlzCMeshExpRmDuplicates(
				  int *nAry,
				  void **ary)
{
  int		n,
		m = 0;

  n = *nAry;
  if(n > 0)
  {
    int		i;

    m = 1;
    for(i = 1; i < n; ++i)
    {
      int	j,
      		dup = 0;

      for(j = 0; j < m; ++j)
      {
	if(ary[i] == ary[j])
	{
	  dup = 1;
	  break;
	}
      }
      if(!dup)
      {
	ary[m++] = ary[i];
      }
    }
  }
  *nAry = m;
}

static void			WlzCMeshExpDebugPrintQueue(
				  WlzCMeshExpWSp *wSp)
{
  if(wSp)
  {
    (void )fprintf(stderr, "WlzCMeshExpDebugPrintQueue()\n");
    if(wSp->head == NULL)
    {
      (void )fprintf(stderr, "wSp->head == NULL\n");
    }
    else if(wSp->tail == NULL)
    {
      (void )fprintf(stderr, "wSp->tail == NULL\n");
    }
    else
    {
      int	cnt,
		idx;
      WlzCMeshExpEnt	*ent;

      cnt = 0;
      ent = wSp->head;
      while(ent)
      {
	switch(wSp->mesh.core->type)
	{
	  case WLZ_CMESH_2D:
	    idx = ent->ent.n2->idx;
	    break;
	  case WLZ_CMESH_3D:
	    idx = ent->ent.n2->idx;
	    break;
	  default:
	    break;
	}
	(void )fprintf(stderr, "  % 6d % 6d % 8g\n",
	                   cnt, idx, wSp->pDst[idx]);
	++cnt;
        ent = ent->nxt;
      }
    }
  }
}
