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
  int		nQ;		  /*!< Number of queue entries in use. */
  int		maxQ;		  /*!< Number of queue entries allocated. */
  int		nEnt;		  /*!< Number of entities (nodes or elements)
  				       in the mesh. */
  WlzCMeshP	mesh;		  /*!< The conforming mesh. */
  WlzIndexedValues *ixv;	  /*!< Indexed values attached to the mesh. */
  WlzUByte	*flags;		  /*!< Per entity flags. */
  int		nNbr;             /*!< Number of neighbours. */
  int		maxNbr;		  /*!< Space allocated for neighbours. */
  WlzCMeshEntPP	nbr;		  /*!< Buffer for neighbours. */
  double	*nVec;		  /*!< Used for computing normal vector using
  				       SVD and the normal vector itself. */
  AlgMatrix	aMat;		  /*!< Used for computing normal vector using
  				       SVD. */
  AlgMatrix	vMat;		  /*!< Used for computing normal vector using
  				       SVD. */
  WlzCMeshExpEnt *head;           /*!< Head of queue of active entities. */
  WlzCMeshExpEnt *tail;           /*!< Tail of queue of active entities. */
  WlzCMeshExpEnt *pool;		  /*!< Entities available for reuse. */
  AlcBlockStack *blocks;	  /*!< Block stack for allocation. */
} WlzCMeshExpWSp;

#define WLZ_CMESH_EXP_MEMINC	(1024)

static WlzCMeshExpWSp		*WlzCMeshExpWSpInit2D(
				  WlzCMesh2D *mesh,
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
static WlzErrorNum		WlzCMeshExpUpdate2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static WlzErrorNum		WlzCMeshExpWSpAddNod2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod2D *nod);
static WlzErrorNum		WlzCMeshExpReallocate2D(
				  WlzCMeshExpWSp *wSp,
				  int nNbr);
static WlzErrorNum		WlzCMeshExpWSpPoolExpand(
				  WlzCMeshExpWSp *wSp);
static void			WlzCMeshExpWSpRecycleEnt(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshExpEnt *ent);
static void			WlzCMeshExpWSpFree(
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
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  /* TODO */
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
       ((wSp->nbr.v = AlcMalloc(initMaxNbr *
				sizeof(WlzCMeshNod2D *))) == NULL) ||
       ((wSp->aMat.rect = AlgMatrixRectNew(initMaxNbr,
       				           nVal, NULL)) == NULL) ||
       ((wSp->vMat.rect = AlgMatrixRectNew(nVal, nVal, NULL)) == NULL))
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
    /* For each neighbour check if it's known and not in the queue yet.
     * If so add it to the active node/element queue. */
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

  if((ent = wSp->head) != NULL)
  {
    if((wSp->head = ent->nxt) != NULL)
    {
      wSp->head->prv = NULL;
    }
  }
  if(ent == wSp->tail)
  {
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
  int		nKnN;
  WlzCMeshNod2D *nod;
  WlzCMeshEdgU2D *eu0,
  		 *eu1;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const int	nVal = 3;

  /* Collect this node and it's known neighbours. */
  nKnN = 1;
  nod = ent->ent.n2;
  wSp->nbr.n2[0] = nod;
  eu0 = eu1 = nod->edu;
  do
  {
    nod = eu1->next->nod;
    /* If this node is known. */
    if((wSp->flags[nod->idx] & WLZ_CMESHEXP_FLAG_UNKNOWN) == 0)
    {
      /* Are the neighbour buffers and matricies big enough? */
      if(nKnN >= wSp->maxNbr)
      {
        errNum = WlzCMeshExpReallocate2D(wSp, nKnN);
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
      wSp->nbr.n2[nKnN] = nod;
      ++nKnN;
    }
    eu1 = eu1->nnxt;
  } while(eu1 != eu0);
  /* Compute plane using SVD. */
  if(errNum == WLZ_ERR_NONE)
  {
    wSp->nNbr = nKnN;
    if(nKnN == 1)
    {
      wSp->nVec[0] = wSp->nVec[1] = wSp->nVec[2] = 0.0;
    }
    else if(nKnN == 2)
    {
      double	len;
      WlzGreyP	gP0,
      		gP1;
      const double eps = 1.0e-06;

      wSp->nVec[0] = wSp->nbr.n2[1]->pos.vtX - wSp->nbr.n2[0]->pos.vtX;
      wSp->nVec[1] = wSp->nbr.n2[1]->pos.vtY - wSp->nbr.n2[0]->pos.vtY;
      gP0.v = WlzIndexedValueGet(wSp->ixv, wSp->nbr.n2[0]->idx);
      gP1.v = WlzIndexedValueGet(wSp->ixv, wSp->nbr.n2[1]->idx);
      switch(wSp->ixv->vType)
      {
	case WLZ_GREY_INT:
	  wSp->nVec[2] = gP1.inp[0] - gP0.inp[0];
	  break;
	case WLZ_GREY_SHORT:
	  wSp->nVec[2] = gP1.shp[0] - gP0.shp[0];
	  break;
	case WLZ_GREY_UBYTE:
	  wSp->nVec[2] = gP1.ubp[0] - gP0.ubp[0];
	  break;
	case WLZ_GREY_FLOAT:
	  wSp->nVec[2] = gP1.flp[0] - gP0.flp[0];
	  break;
	case WLZ_GREY_DOUBLE:
	  wSp->nVec[2] = gP1.dbp[0] - gP0.dbp[0];
	  break;
	default:
	  break;
      }
      len = sqrt((wSp->nVec[0] * wSp->nVec[0]) +
                 (wSp->nVec[1] * wSp->nVec[1]));
      if(len > eps)
      {
        wSp->nVec[0] /= len;
	wSp->nVec[1] /= len;
	wSp->nVec[2] /= len;
      }
      else
      {
        wSp->nVec[0] = wSp->nVec[1] = wSp->nVec[2] = 0.0;
      }
    }
    else
    {
      /* Find normal vector using SVD. */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i,
		k;
	double	*row;

	wSp->aMat.rect->nR = nKnN;
	/* Collect positions and values. */
	for(k = 0; k < nKnN; ++k)
	{
	  WlzGreyP	gP;
	  WlzCMeshNod2D *nod;

	  nod = wSp->nbr.n2[k];
	  row = wSp->aMat.rect->array[k];
	  row[0] = nod->pos.vtX;
	  row[1] = nod->pos.vtY;
	  gP.v = WlzIndexedValueGet(wSp->ixv, nod->idx);
	  switch(wSp->ixv->vType)
	  {
	    case WLZ_GREY_INT:
	      row[2] = gP.inp[0];
	      break;
	    case WLZ_GREY_SHORT:
	      row[2] = gP.shp[0];
	      break;
	    case WLZ_GREY_UBYTE:
	      row[2] = gP.ubp[0];
	      break;
	    case WLZ_GREY_FLOAT:
	      row[2] = gP.flp[0];
	      break;
	    case WLZ_GREY_DOUBLE:
	      row[2] = gP.dbp[0];
	      break;
	    default:
	      break;
	  }
	}
	/* Remove centroid. */
	for(i = 0; i < nVal; ++i)
	{
	  double	cen = 0.0;

	  for(k = 0; k < nKnN; ++k)
	  {
	    cen += wSp->aMat.rect->array[k][i];
	  }
	  cen /= nKnN;
	  for(k = 0; k < nKnN; ++k)
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
  WlzCMeshNod2D *gnod;
  WlzCMeshEdgU2D *eu0,
  		 *eu1;
  WlzGreyP	gval;
  const double	eps = 1.0e-06;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  gnod = ent->ent.n2;                        /* Given node with known value. */
  gval.v = WlzIndexedValueGet(wSp->ixv, gnod->idx);
  /* For each of the neighbours of the given known entity which is unknown. */
  eu0 = eu1 = gnod->edu;
  do
  {
    int		  nKN;
    WlzCMeshEdgU2D *eu2,
		   *eu3;
    WlzCMeshNod2D *unod;

    unod = eu1->next->nod;
    /* If this node is unknown. */
    if((wSp->flags[unod->idx] & WLZ_CMESHEXP_FLAG_UNKNOWN) != 0)
    {
      double	  gv,
      		  uv;
      WlzDVertex2 del;
      WlzGreyP	  uval;

      /* Compute new value and if appropriate set value. 
       * There normal vector has been computed using AlgMatrixSVDecomp()
       * and so is known to be a unit length vector. */
      WLZ_VTX_2_SUB(del, unod->pos, gnod->pos);
      uval.v = WlzIndexedValueGet(wSp->ixv, unod->idx);
      switch(wSp->ixv->vType)
      {
	case WLZ_GREY_INT:
	  gv = gval.inp[0];
	  uv = uval.inp[0];
	  break;
	case WLZ_GREY_SHORT:
	  gv = gval.shp[0];
	  uv = uval.shp[0];
	  break;
	case WLZ_GREY_UBYTE:
	  gv = gval.ubp[0];
	  uv = uval.ubp[0];
	  break;
	case WLZ_GREY_FLOAT:
	  gv = gval.flp[0];
	  uv = uval.flp[0];
	  break;
	case WLZ_GREY_DOUBLE:
	  gv = gval.dbp[0];
	  uv = uval.dbp[0];
	  break;
        default:
	  break;
      }
      if(fabs(wSp->nVec[2]) < eps)
      {
	if((wSp->flags[unod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0)
	{
	  uv = gv;
	}
      }
      else
      {
	double v;

	if(wSp->nNbr == 2)
	{
	  v = gv + wSp->nVec[2] * ((del.vtX * wSp->nVec[0]) +
				   (del.vtY * wSp->nVec[1]));
	}
	else /* wSp->nNbr >= 3 */
	{
	  v = gv - 
	      (((wSp->nVec[0] * del.vtX) +
		(wSp->nVec[1] * del.vtY)) / wSp->nVec[2]);
	}
	if(((wSp->flags[unod->idx] & WLZ_CMESHEXP_FLAG_UPDATED) == 0) ||
	   ((v < gv) && (v > uv)) ||
	   ((v > gv) && (v < uv)))
	{
	  uv = v;
	  /* Only set the updated flag if the extrapolation is good. */
	  if(wSp->nNbr >= 3)
	  {
	    wSp->flags[unod->idx] |= WLZ_CMESHEXP_FLAG_UPDATED;
	  }
	}
      }
      switch(wSp->ixv->vType)
      {
	case WLZ_GREY_INT:
	  uval.inp[0] = (int )WLZ_CLAMP(uv, INT_MIN, INT_MAX);
	  break;
	case WLZ_GREY_SHORT:
	  uval.shp[0] = (short )WLZ_CLAMP(uv, SHRT_MIN, SHRT_MAX);
	  break;
	case WLZ_GREY_UBYTE:
	  uval.ubp[0] = (WlzUByte )WLZ_CLAMP(uv, 0, 255);
	  break;
	case WLZ_GREY_FLOAT:
	  uval.flp[0] = (float )WLZ_CLAMP(uv, -(FLT_MAX), FLT_MAX);
	  break;
	case WLZ_GREY_DOUBLE:
	  uval.dbp[0] = uv;
	  break;
        default:
	  break;
      }
      /* Find number of neighbouring known active entities of this
       * unkown neighbour. */
      nKN = 0;
      eu2 = eu3 = unod->edu;
      do
      {
	WlzCMeshNod2D *tnod;

	tnod = eu3->next->nod;
        if((wSp->flags[tnod->idx] & 
	    (WLZ_CMESHEXP_FLAG_UNKNOWN | WLZ_CMESHEXP_FLAG_ACTIVE)) ==
	    WLZ_CMESHEXP_FLAG_ACTIVE)
	{
	  ++nKN;
	}
        eu3 = eu3->nnxt;
      } while(eu3 != eu2);
      /* If the given entity is the unknown neighbour's only active neighbour
       * then make the unknown neighbour known and add it to the active entity
       * queue. */
      if(nKN == 1)
      {
        wSp->flags[unod->idx] = WLZ_CMESHEXP_FLAG_ACTIVE;
	errNum = WlzCMeshExpWSpAddNod2D(wSp, unod);
      }
    }
    eu1 = eu1->nnxt;
  } while((errNum == WLZ_ERR_NONE) && (eu1 != eu0));
  wSp->flags[gnod->idx] = 0;      /* The given node is now no longer active. */
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Adds the given node to the extrapolation queue. The given
* 		node is known not to be already in the quque.
* \param	wSp			Given mesh extrapolation workspace.
* \param	nod			Given 2D node.
*/
static WlzErrorNum		WlzCMeshExpWSpAddNod2D(
				  WlzCMeshExpWSp *wSp,
				  WlzCMeshNod2D *nod)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((wSp->maxQ - wSp->nQ) < 1)
  {
    errNum = WlzCMeshExpWSpPoolExpand(wSp);
  }
  /* Take an entry from the pool, set it up and put it at the queue's tail. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshExpEnt *ent;

    ++(wSp->nQ);
    ent = wSp->pool;
    wSp->pool = ent->nxt;
    ent->ent.n2 = nod;
    ent->nxt = NULL;
    if((ent->prv = wSp->tail) != NULL)
    {
      ent->prv->nxt = ent;
    }
    wSp->tail = ent;
    if(wSp->head == NULL)
    {
      wSp->head = ent;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Reallocates the 2D neighbour buffer and SVD matrix.
* \param	wSp			Given mesh extrapolation workspace.
* \param	nNbr			New number of neighbours.
*/
static WlzErrorNum		WlzCMeshExpReallocate2D(
				  WlzCMeshExpWSp *wSp,
				  int nNbr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(nNbr >= wSp->maxNbr)
  {
    int		n,
    		maxNbr;
    const int   nVal = 3,
    	        inc = WLZ_CMESH_EXP_MEMINC;

    n = (nNbr + inc - 1) / inc;
    maxNbr = n * inc;
    if((wSp->nbr.v = AlcRealloc(wSp->nbr.v,
	    maxNbr * sizeof(WlzCMeshNod2D *))) == NULL)
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
    wSp->maxQ += inc;
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
  --(wSp->nQ);
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
    AlcFree(wSp->nbr.v);
    AlcBlockStackFree(wSp->blocks);
    AlgMatrixRectFree(wSp->aMat.rect);
    AlgMatrixRectFree(wSp->vMat.rect);
  }
}
