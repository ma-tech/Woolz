#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshSurfMap_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCMeshSurfMap.c
* \author       Bill Hill
* \date         May 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
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
* \brief	Functions for computing surface mappings that are based
* 		on conformal transformations and other mappings between
* 		surfaces and planes.
* \ingroup	WlzTransform
*/

#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

static int			WlzCMeshSurfMapIdxCmpFn(
				  const void *p0,
				  const void *p1);
static WlzGMModel 		*WlzCMeshToGMModel2D(
				  WlzObject *mObj,
				  double disp,
				  WlzErrorNum *dstErr);
static WlzGMModel 		*WlzCMeshToGMModel2D5(
				  WlzObject *mObj,
				  double disp,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshCompSurfMapFT(
				  WlzCMesh2D5 *mesh,
				  int nP,
				  WlzDVertex3 *pV,
				  int *pIdx,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzCMeshCompSurfMapLevy(
				  WlzCMesh2D5 *mesh,
				  int nP,
				  WlzDVertex3 *dPV,
				  int *pIdx,
				  WlzErrorNum *dstErr);

/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a isomophic mapping of the given surface patch
* 		(which must be isomorphic to a disk) to the given spatial
* 		domain object (which must have a single 2D convex boundary
* 		polygon.
* \param	mshObj			Input conforming mesh object which
* 					must be of type WLZ_CMESH_2D5.
* \param	domObj			Given planar spatial domain, which
* 					must be a single piece.
* \param	nDV			Number of spatial domain vertices
* 					which must be the same as nMV.
* \param	dV			Spatial domain vertices which must
* 					be on or near the domain boundary.
* 					These are 3D vertices with the Z
* 					component zero.
* \param	nMV			Number of in mesh vertices.
* \param	mV			Mesh vertices.
* \param	dstErr			Woolz error code, may be NULL.
*/
WlzObject			*WlzCMeshCompSurfMapToDomain(
				  WlzObject *mshObj,
				  WlzObject *domObj,
				  int nDV, WlzDVertex3 *dV,
				  int nMV, WlzDVertex3 *mV,
				  WlzErrorNum *dstErr)
{
  int		nPin = 0,	/* Number of pin vertex, node index pairs. */
		nPlyBV = 0,	/* Number of boundary polygon vertices. */
  		nMshBN = 0,     /* Number of boundary mesh nodes. */
  		mshBNInc = 1;   /* +/-1 direction of mesh boundary nodes. */
  int	 	*mshNI = NULL,  /* Indices to the boundary mesh node indices
                                 * which are closest to the given mesh
				 * vertices, ie
				 * min((mv[i] - nod->pos)^2), where nod->idx
				 * is in mshBNI[]. */
		*prmI  = NULL,  /* Permutation to make given points cyclic. */
		*pinNI = NULL,  /* Pin node indices. */
  		*plyVI = NULL,  /* Given point indices in boundary polygon. */
		*mshBNI = NULL; /* Mesh boundary node indices. */
  AlcVector	*mshNV;         /* Mesh node vector. */
  WlzIVertex2	*plyV = NULL;
  WlzDVertex2	*mshV2 = NULL,
  		*plyV2 = NULL;
  WlzDVertex3	*pinV = NULL;
  WlzBasisFnTransform *bFn = NULL;
  WlzPolygonDomain *ply = NULL;
  WlzObject	*prmObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;      
  WlzCMesh2D5	*mesh = NULL;

  /* Check object. */
  if((mshObj == NULL) || (domObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((mshObj->type != WLZ_CMESH_2D5) ||
          (domObj->type != WLZ_2D_DOMAINOBJ))
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(((mesh = mshObj->domain.cm2d5) == NULL) ||
          (domObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((nMV <= 3) || (nMV != nDV))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((mV == NULL) || (dV == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    mshNV = mesh->res.nod.vec;
    /* Allocate polygon and mesh node index arrays for the closest
     * vertices and nodes to the given vertices. */
    if(((prmI = (int *)AlcCalloc(nMV, sizeof(int))) == NULL) ||
       ((mshNI = (int *)AlcCalloc(nMV, sizeof(int))) == NULL) ||
       ((plyVI = (int *)AlcCalloc(nMV, sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject	*bndObj;

    /* Get the boundary polygon of the given 2D spatial domain object,
     * making sure that it is a simple chain of equi-spaced vertices
     * (but including the given vertices). */
    bndObj = WlzObjToBoundary(domObj, 0, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      if((bndObj == NULL) || (bndObj->type != WLZ_BOUNDLIST) ||
         (bndObj->domain.b == NULL) ||
	 (bndObj->domain.b->up != NULL) ||
	 (bndObj->domain.b->next != NULL) ||
	 (bndObj->domain.b->down != NULL) ||
	 (bndObj->domain.b->poly == NULL))
      {
        errNum = WLZ_ERR_DOMAIN_DATA;
      }
      else
      {
        ply = WlzPolyEquispace(bndObj->domain.b->poly,
			       0, 1.0, 1, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  plyV = ply->vtx;
          nPlyBV = ply->nvertices;
	}
      }
    }
    (void )WlzFreeObj(bndObj);
  }
  /* Find the indices of the boundary nodes of the mesh and get them as
   * a cyclic ordered array. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshGetBoundNodes2D5(mesh, &nMshBN, &mshBNI, 1);
  }
  /* Find the indices of the closest 2D domain boundary polygon vertices to
   * each of the given 2D vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nDV; ++i)
    {
      int	j;
      double	d;
      WlzDVertex2 gV,
      		  fV;

      gV.vtX = dV[i].vtX;
      gV.vtY = dV[i].vtY;
      WLZ_VTX_2_SUB(fV, gV, plyV[0]);
      d = WLZ_VTX_2_SQRLEN(fV);
      for(j = 1; j < nPlyBV; ++j)
      {
	double	    d1;
        WlzDVertex2 bV,
	            tV;

	bV.vtX = plyV[j].vtX;
	bV.vtY = plyV[j].vtY;
	WLZ_VTX_2_SUB(tV, bV, gV);
	d1 = WLZ_VTX_2_SQRLEN(tV);
	if(d1 < d)
	{
	  d = d1;
	  plyVI[i] = j;
	}
      }
    }
  }
  /* Find the indices of the closest mesh boundary nodes to each of the
   * given 3D vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nMV; ++i)
    {
      int	j;
      double	d = DBL_MAX;
      WlzDVertex3 gV;

      gV = mV[i];
      for(j = 0; j < nMshBN; ++j)
      {
	double    d1;
	WlzDVertex3 tV;
        WlzCMeshNod2D5 *nod;

	/* All nodes will be valid on boundary. */
	nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[j]);
	WLZ_VTX_3_SUB(tV, nod->pos, gV);
	d1 = WLZ_VTX_3_SQRLEN(tV);
	if(d1 < d)
	{
	  d = d1;
	  mshNI[i] = j;
	}
      }
    }
    /* Determine boundary node direction increment. +/-1. */
    for(i = 1; i < nMshBN; ++i)
    {
      int	j;

      j = (mshNI[0] + i) % nMshBN;
      if(mshBNI[j] == mshBNI[mshNI[1]])
      {
        mshBNInc = 1;
	break;
      }
      else if(mshBNI[j] == mshBNI[mshNI[nMV - 1]])
      {
        mshBNInc = -1;
	break;
      }
    }
#ifdef WLZ_CMESH_SURFMAP_DEBUG
    (void )fprintf(stderr, "WLZ_CMESH_SURFMAP_DEBUG\n");
    for(i = 0; i < nMV; ++i)
    {
      WlzCMeshNod2D5 *nod;

      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[mshNI[i]]);
      (void )fprintf(stderr, "%d % 6d    %g %g %g\n",
                     i, nod->idx,
		     nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ);
    }
#endif
  }
  /* Compute the permutation which makes the given 2D vertices a cycle on
   * the given domain's boundary. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nMV; ++i)
    {
      prmI[i] = i;
    }
    AlgHeapSortIdx(plyVI, prmI, nMV, AlgHeapSortCmpIdxIFn);
  }
  /* Allocate vertex and mesh node pin arrays such that they will always have
   * less elements that the number allocated. */
  if(errNum == WLZ_ERR_NONE)
  {

    nPin = ALG_MAX(nPlyBV, nMshBN);
    if(((pinNI = (int *)AlcMalloc(sizeof(int) * nPin)) == NULL) ||
       ((pinV = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * nPin)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Using the permutation index array to ensure that the spatial domain
   * boundary vertices are in a cyclic ordering (around the domain):
   * For each inter mesh boundary node segment interpolate using the
   * inter-node and inter-polygon distances to find the spatial domain
   * boundary polygon vertex. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    nPin = 0;
    for(i = 0; i < nMV; ++i)
    {
      int	iN,
		iV0,
      		vI0,    /* Cyclic ordered index of 1st vertex/node. */
      		vI1,   	/* Cyclic ordered index of 2nd vertex/node. */
		nI0,	/* 1st boundary node index, ie index into mshBNI. */
		nI1;	/* 2nd boundary node index, ie index into mshBNI. */
      double	d,
		dpl,	/* Distance between current and next matched
		         * boundary polygon vertices. */
      		dmn;    /* Distance between current and next matched
                         * boundary mesh nodes. */
      WlzCMeshNod2D5 *nod[2];

      dmn = 0.0;
      vI0 = prmI[i % nMV];
      vI1 = prmI[(i + 1) % nMV];
      nI0 = mshNI[vI0];
      nI1 = mshNI[vI1];
      /* Compute distance between current and next vertex (easy because
       * the poly equi spaced. */
      dpl = (nPlyBV + plyVI[vI1] - plyVI[vI0]) % nPlyBV;
      /* Compute distance between current and next matched node. */
      iN = nI0;
      nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[nI0]);
      while(iN != nI1)
      {
	WlzDVertex3 t;

        nod[0] = nod[1];
	iN = (iN + mshBNInc + nMshBN) % nMshBN;
	nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[iN]);
	WLZ_VTX_3_SUB(t, nod[0]->pos, nod[1]->pos);
	dmn += WLZ_VTX_3_LENGTH(t);
      }
      /* Now interpolate the polgon vertices corresponding to the intermediate
       * boundary mesh nodes if there are any. */
      d = 0;              
      iV0 = -1;
      iN = nI0;           
      nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[nI0]);
      while(iN != nI1)    
      {
	int	iV1;
	double	s;
	WlzIVertex2 *v2;
        WlzDVertex3 t3;

	nod[0] = nod[1];
	s = d * dpl / dmn;
	iV1 = ((plyVI[vI0] + (int )((floor)(s))) % nPlyBV);
	if((iV1 - iV0) > 0)
	{
	  iV0 = iV1;
	  v2 = plyV + ((plyVI[vI0] + (int )((floor)(s))) % nPlyBV);
	  pinNI[nPin] = nod[0]->idx;
	  pinV[nPin].vtX = v2->vtX;
	  pinV[nPin].vtY = v2->vtY;
	  pinV[nPin].vtZ = 0.0;
	  ++nPin;
	}
	iN = (iN + mshBNInc + nMshBN) % nMshBN;
	nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[iN]);
	WLZ_VTX_3_SUB(t3, nod[0]->pos, nod[1]->pos);
	d += WLZ_VTX_3_LENGTH(t3);
      }
    }
    /* Avoid possible duplicated start node at end. */
    if((nPin > 1) && (pinNI[nPin - 1] == pinNI[0]))
    {
      --nPin;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nPin < 3)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  AlcFree(prmI);
  AlcFree(mshNI);
  AlcFree(plyVI);
  AlcFree(mshBNI);
  (void )WlzFreePolyDmn(ply);
  /* Compute the parameterised mesh using the pinned nodes and their
   * parameterised positions. */
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_CMESH_SURFMAP_DEBUG
    int		i;

    (void )fprintf(stderr, "WLZ_CMESH_SURFMAP_DEBUG\n");
    for(i = 0; i < nPin; ++i)
    {
      WlzCMeshNod2D5 *nod;

      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, pinNI[i]);
      (void )fprintf(stderr, "%d    %g %g    % 6d %g %g %g\n",
                     i, pinV[i].vtX, pinV[i].vtY,
		     nod->idx, nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ);
    }
#endif
    prmObj = WlzCMeshCompSurfMapFT(mesh, nPin, pinV, pinNI, &errNum);
  }
  AlcFree(pinV);
  AlcFree(pinNI);
  AlcFree(mshV2);
  AlcFree(plyV2);
  (void )WlzBasisFnFreeTransform(bFn);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(prmObj);
    prmObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(prmObj);
}

/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a least squares conformal transformation which
* 		maps the source surface to a destination plane with z
* 		coordinate zero. When matching the mesh vertices from
* 		their position, the closest vertices to the given positions
* 		are used. See WlzCMeshCompSurfMapLevy().
* \param	inObj			Input conforming mesh object which
* 					must be of type WLZ_CMESH_2D5.
* \param	nDV			Number of destination vertices.
* \param	dV			Destination vertices.
* \param	nSV			Number of source vertices
* 					which must be the same as nDV.
* \param	sV			Source vertices.
* \param	dstErr			Woolz error code, may be NULL.
*/
WlzObject	*WlzCMeshCompSurfMap(WlzObject *inObj,
				int nDV, WlzDVertex3 *dV,
				int nSV, WlzDVertex3 *sV,
				WlzErrorNum *dstErr)
{
  int		*nodTb = NULL;
  WlzCMesh2D5	*mesh = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((nDV < 1) || (nDV != nSV) || (dV == NULL) || (sV == NULL))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    WlzDomain	dom;

    switch(inObj->type)
    {
      case WLZ_CMESH_2D5:
        mesh = inObj->domain.cm2d5;
        break;
      case WLZ_CONTOUR:
        mesh = WlzCMeshFromGM(inObj->domain.ctr->model, &errNum);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    dom.cm2d5 = mesh;
    (void )WlzAssignDomain(dom, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((mesh->res.nod.numEnt < 3) || (mesh->res.elm.numEnt < 1))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if((nodTb = (int *)AlcMalloc(sizeof(int) * nSV)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int 	i;

      for(i = 0; i < nSV; ++i)
      {
        if((nodTb[i] = WlzCMeshClosestNod2D5(mesh, sV[i])) < 0)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnObj = WlzCMeshCompSurfMapLevy(mesh, nDV, dV, nodTb, &errNum);
  }
  AlcFree(nodTb);
  if(mesh)
  {
    (void )WlzCMeshFree2D5(mesh);
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
}


/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a parameterisation of the mesh using the
* 		barycenric method of Tutte and Floater, ie:
* 		W. T. Tutte "How to Draw a Graph", Proceedings of the
* 		London Mathematical Society 13: 743–767 1963
*		and
*		Michael S. Floater "Parametrization and Smooth
*		Approximation of Surface Triangulations", Computer Aided
*		Geometric Design 14(3): 231-250 1997.
*
* 		Given a WlzCMesh2D5 with disk topology and a set of paired
* 		pin boundary vertex, planar pin locations; this function
* 		computes displacements for the mesh nodes so that they
* 		are arranged in the plane at the barycentres of enclosing
* 		polygons formed by the edge connected neighbours of the
* 		nodes. The given pins establish a mapping of the mesh
* 		boundary nodes to a planar (z = 0) convex polygon.
*
* 		For \f$j \in \{0, ..., n_i - 1\}\f$ a set of coefficients
* 		are defined with:
* 		\f$-a_{jj} = \sum_{k \neq j}{a_{jk}\vec{u}}\f$
* 		and
* 		\f$a_{jk} = 1\f$ (\f$j \neqk\f$),
* 		\f$a_{jj} = -|N(j)|\f$
* 		where \f$N(j)\f$ is the set of nodes which ar edge connected
* 		neighbours of node \f$j\f$ and \f$|N(j)|\f$ is the cardinality
* 		of the set.
*
* 		To compute the solution two matrices are computed
* 		\f$A\f$ and \f$B\f$ with:
* 		\f$A = [a_{jk}]\f$ with \f$j,k\f$ covering the set of
* 		internal nodes and
* 		\f$B = [b_{jk}]\f$ here with \f$j,k\f$ covering the set of
* 		pinned boundary nodes.
* 		The positions of the boundary nodes \f$\vec{u}_b\f$ are
* 		known so the positions of the interior nodes
* 		\f$\vec{u}_i\f$ are found by solving
* 		\f$A\vec{u}_i = B\vec{u}_b\f$ which is of the form
* 		\f$A\vec{x} = \vec{b\}\f$
* 		This equation is solved for both components of \f$\vec{u}\f$,
* 		ie \f$U\f$ and \f$v\f$.
* \param	mesh			Input conforming mesh which must be
* 					of type WLZ_CMESH_2D5.
* \param	nP			Number of pinned nodes (these are the
* 					boundary nodes of the mesh).
* \param	pV			Coordinates of the pinned nodes.
* 					All z components should be equal and
* 					usualy set to zero. The coordinates
* 					must correspond to the indices of
* 					pIdx.
* \param	pIdx			Indices of the pinned nodes (ie all
* 					the mesh boundary nodes. The indices
* 					must correspond to the pin locations.
* \param	dstErr			Woolz error code, may be NULL.
*/
WlzObject	*WlzCMeshCompSurfMapFT(WlzCMesh2D5 *mesh,
				int nP, WlzDVertex3 *pV, int *pIdx,
				 WlzErrorNum *dstErr)
{
  int		nI,			/* Number of interior nodes. */
                nN,			/* Total number of nodes. */
		nM;			/* Maximum node index. */
  int		*nIdx = NULL,		/* Partitioned array of nP pinned
                                           nodes first and then nI internal
					   nodes. Each of the two partitions
					   are sorted such that for the
					   pinned nodes:
					   nIdx[i] = pIdx[sIdx[i]]
					   pIdx[i] = nIdx[rIdx[i]]. */
  		*rIdx = NULL,		/* Array of reverse node indices
		                           from node index to offset into
					   nIdx array for each node. */
                *sIdx = NULL;		/* Indices such that pinned nodes
					   are sorted when accessed via this
					   index, ie:
					   pIdx[sPidx[i+1]] > pIdx[sPidx[i]]. */
  double	*bV = NULL,
		*lV = NULL,		/* The lambda values. */
		*xV = NULL;
  AlgMatrix	aM,			/* Matrix A of nI x nI weights. */
                wM; 			/* Rectangular matrix 4 x nI for CG. */
  WlzUByte	**qTab = NULL;		/* Table with qTab[i][j] set to 1
                                           iff i'th internal node is a
					   neighbour of the j'th pinned
					   node. */
  AlcVector	*nV = NULL;		/* Mesh node vector. */
  WlzObject	*mapObj = NULL;
  WlzIndexedValues *ixv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	itr = 1000;
  const double	tol = 0.000001;

  aM.core = NULL;
  wM.core = NULL;
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(nP < 3)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((pV == NULL) || (pIdx == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  /* Create return object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	dom;
    WlzValues	val;

    dom.cm2d5 = mesh;
    val.core = NULL;
    mapObj = WlzMakeMain(WLZ_CMESH_2D5, dom, val, NULL, NULL, &errNum);
  }
  /* Allocate indexed values. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 3;
    WlzValues	val;

    val.x = WlzMakeIndexedValues(mapObj, 1, &dim, WLZ_GREY_DOUBLE,
	                         WLZ_VALUE_ATTACH_NOD, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      ixv = val.x;
      nV = mesh->res.nod.vec;
      mapObj->values = WlzAssignValues(val, NULL);
    }
  }
  /* Set offsets for the pinned vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    for(i = 0; i < nP; ++i)
    {
      double    *dsp;
      WlzDVertex3 *p,
      		  *q;
      WlzCMeshNod2D5 *nod;

      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nV, pIdx[i]);
      dsp = (double *)WlzIndexedValueGet(ixv, pIdx[i]);
      p = &(pV[i]);
      q = &(nod->pos);
      dsp[0] = p->vtX - q->vtX;
      dsp[1] = p->vtY - q->vtY;
      dsp[2] = -(q->vtZ);
    }
  }
  /* Compute offsets for the interior vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    nM = mesh->res.nod.maxEnt;
    nN = mesh->res.nod.numEnt;
    nI = nN - nP;
    if(nI < 0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
    else if(nI > 0)
    {
      /* Allocate index arrays and the matrices. */
      if(((nIdx    = (int *)AlcMalloc(sizeof(int) * (nN + nP))) == NULL) ||
         ((rIdx    = (int *)AlcMalloc(sizeof(int) * nM)) == NULL) ||
         ((sIdx    = (int *)AlcMalloc(sizeof(int) * nP)) == NULL) ||
         ((bV      = (double *)AlcMalloc(sizeof(double) * nI)) == NULL) ||
         ((lV      = (double *)AlcMalloc(sizeof(double) * nI)) == NULL) ||
         ((xV      = (double *)AlcMalloc(sizeof(double) * nI)) == NULL) ||
         ((aM.llr = AlgMatrixLLRNew(nI, nI, 10 * nI, tol, NULL)) == NULL) ||
         ((wM.rect = AlgMatrixRectNew(4, nI, NULL)) == NULL) ||
         (AlcUnchar2Calloc(&qTab, nI, nP) != ALC_ER_NONE))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        int 	i,
		j,
		k;

	/* Compute a sorted pinned node index permutation table in the last
	 * nP of the node index table. */
	for(i = 0, j = 0; i < nM; ++i)
	{
	  WlzCMeshNod2D5 *nod;

	  nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nV, i);
	  if(nod->idx >= 0)
	  {
	    nIdx[j++]= nod->idx;
	  }
	}
	if(j != nN)
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else
	{
	  AlgSort((void *)nIdx, nN, sizeof(int), WlzCMeshSurfMapIdxCmpFn);
	  for(i = 0; i < nP; ++i)
	  {
	    sIdx[i] = i;
	  }
	  AlgHeapSortIdx(pIdx, sIdx, nP, AlgHeapSortCmpIdxIFn);
	  for(i = 0; i < nP; ++i)
	  {
	    nIdx[nN + i] = pIdx[sIdx[i]];
	  }
	  /* Add the internal nodes. */
	  for(i = 0, j = 0, k = 0; i < nN; ++i)
	  {
	    int	t;

	    t = nIdx[i] - nIdx[nN + j];
	    if(t == 0)
	    {
	      if(++j >= nP)
	      {
		j = nP - 1;
	      }
	    }
	    else
	    {
	      nIdx[k++] = nIdx[i];
	    }
	  }
	  if(k != nI)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  else
	  {
	    for(i = 0; i < nP; ++i)
	    {
	      nIdx[nI + i] = nIdx[nN + i];
	    }
	    for(i = 0; i < nM; ++i)
	    {
	      rIdx[nIdx[i]] = i;
	    }
	  }
	}
      }
      /* Set values for matrix A (aM) and table (qTab) used to
       * compute b. */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i;

	AlgMatrixSetAll(aM, 0.0);
        for(i = 0; i < nI; ++i)
	{
	  WlzCMeshNod2D5 *iNod;
	  WlzCMeshEdgU2D5 *edu;

	  iNod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nV, nIdx[i]);
	  if(iNod && (iNod->idx >= 0) && (iNod->edu != NULL))
	  {
	    int	   nCnt;

	    /* First count number of edge connected neighbours. */
	    nCnt = 0;
	    edu = iNod->edu;
	    do
	    {
	      ++nCnt;
	      edu = edu->nnxt;
	    }
	    while(edu != iNod->edu);
	    /* Now set elements of matrix A. */
	    if(nCnt > 0)
	    {
	      lV[i] = 1.0 / nCnt;
	      edu = iNod->edu;
	      do
	      {
		int	j;
		WlzCMeshEdgU2D5 *nnxt;
		nnxt = edu->nnxt;
                if(nnxt->opp)
		{
		  WlzCMeshNod2D5 *jNod;

		  jNod = nnxt->opp->nod;
		  if((j = rIdx[jNod->idx]) >= nI)
		  {
		    int	k;

		    /* Neighbour is pinned to boundary. */
		    k = j - nI;
		    *(*(qTab + i) + k) = 1;
		  }
		  else
		  {
		    /* Neighbour is internal. */
		    (void )AlgMatrixSet(aM, i, j, -lV[i]);
		  }
		}
		edu = nnxt;
	      }
	      while(edu != iNod->edu);
	      (void )AlgMatrixSet(aM, i, i, 1.0);
	    }
	    else
	    {
	      errNum = WLZ_ERR_DOMAIN_DATA;
	      break;
	    }
	  }
	}
      }
      /* Set bV and xV using the pinned node X coordinates and the
       * table (qTab) then solve for xV */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i;

	for(i = 0; i < nI; ++i)
	{
	  int	j;
	  WlzUByte *q;
	  double b = 0.0;

	  q = qTab[i];
	  for(j = 0; j < nP; ++j)
	  {
	    if(q[j])
	    {
	      b += lV[i] * pV[sIdx[j]].vtX;
	    }
	  }
	  bV[i] = b;
	  xV[i] = 0.0;
	}
	errNum = WlzErrorFromAlg(
	    AlgMatrixCGSolve(aM, xV, bV, wM, NULL, NULL,
	                     tol, itr, NULL, NULL));
      }
      /* Set X coordinates of the interior nodes in the new map object. */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i;

        for(i = 0; i < nI; ++i)
	{
	  double    *dsp;
	  WlzCMeshNod2D5 *nod;

	  nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nV, nIdx[i]);
	  dsp = (double *)WlzIndexedValueGet(ixv, nIdx[i]);
	  dsp[0] = xV[i] - nod->pos.vtX;
	  dsp[2] = -(nod->pos.vtZ);
	}
      }
      /* Set bV and xV using the pinned node Y coordinates and the
       * table (qTab) then solve for xV */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i;

	for(i = 0; i < nI; ++i)
	{
	  int	j;
	  WlzUByte *q;
	  double b = 0.0;

	  q = qTab[i];
	  for(j = 0; j < nP; ++j)
	  {
	    if(q[j])
	    {
	      b += lV[i] * pV[sIdx[j]].vtY;
	    }
	  }
	  bV[i] = b;
	  xV[i] = 0.0;
	}
	errNum = WlzErrorFromAlg(
	    AlgMatrixCGSolve(aM, xV, bV, wM, NULL, NULL,
	                     tol, itr, NULL, NULL));
      }
      /* Set Y coordinates of the interior nodes in the new map object. */
      if(errNum == WLZ_ERR_NONE)
      {
	int	i;

        for(i = 0; i < nI; ++i)
	{
	  double    *dsp;
	  WlzCMeshNod2D5 *nod;

	  nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nV, nIdx[i]);
	  dsp = (double *)WlzIndexedValueGet(ixv, nIdx[i]);
	  dsp[1] = xV[i] - nod->pos.vtY;
	}
      }
    }
  }
  AlcFree(bV);
  AlcFree(lV);
  AlcFree(xV);
  AlcFree(nIdx);
  AlcFree(rIdx);
  AlcFree(sIdx);
  AlgMatrixFree(aM);
  AlgMatrixFree(wM);
  (void )Alc2Free((void **)qTab);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(mapObj);
    mapObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mapObj);
}

/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a least squares conformal transformation which
* 		maps the source surface to a destination plane with z
* 		coordinate zero.
* 		The algorithm used here is based on the paper: Bruno L'evy,
*               etal "Least Squares Conformal Maps for Automatic Texture
*               Atlas Generation" SIGGRAPH 2002.
* \param	mesh			Input conforming mesh which must be
* 					of type WLZ_CMESH_2D5.
* \param	nP			Number of pinned nodes.
* \param	dPV			Destination coordinates of the
* 					pinned nodes. All z components
* 					should be equal and usualy set to
* 					zero. The coordinates must correspond
* 					to the indices of pIdx.
* \param	pIdx			Indices of the pinned nodes
* 					which must all be valid. The indices
* 					must be sorted soo that the indices
* 					increase monotonically.
* \param	dstErr			Woolz error code, may be NULL.
*/
WlzObject	*WlzCMeshCompSurfMapLevy(WlzCMesh2D5 *mesh,
				int nP, WlzDVertex3 *dPV, int *pIdx,
				 WlzErrorNum *dstErr)
{
  int		mE,
		nE2,
		nP2,
		mN,
		nN2,
		nE = 0,	      /* Just to avoid invalid uninitialized warning. */
		nN = 0;       /* Just to avoid invalid uninitialized warning. */
  int		*pIdxSortTb = NULL,
     		*pIdxIdxTb = NULL,
  		*eIdxTb = NULL,
  		*nIdxTb = NULL,
		*pIdxSorted = NULL;
  double	*bV = NULL,
  		*bUV = NULL,
		*xV = NULL;
  AlgMatrix	aM,
  		bPM;
  WlzObject	*mapObj = NULL;
  WlzIndexedValues *ixv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 0.000001;

  aM.core = NULL;
  bPM.core = NULL;
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((nP < 2) || (mesh->res.nod.numEnt < nP) ||
          (mesh->res.elm.numEnt < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  /* Allocate matrices, element and node index look up tables. */
  if(errNum == WLZ_ERR_NONE)
  {
    mE = mesh->res.elm.maxEnt;
    mN = mesh->res.nod.maxEnt;
    nN = mesh->res.nod.numEnt;
    nE = mesh->res.elm.numEnt;
    nE2 = 2 * nE;
    nN2 = 2 * nN;
    nP2 = 2 * nP;
    if(((aM.llr = AlgMatrixLLRNew(nE2, nN2, 6 * nE, tol, NULL)) == NULL) ||
       ((bPM.llr = AlgMatrixLLRNew(nE2, nP2, 6 * nE, tol, NULL)) == NULL) ||
       ((bUV = (double *)AlcMalloc(sizeof(double) * nP2)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nE2)) == NULL) ||
       ((xV = (double *)AlcCalloc(nE2, sizeof(double))) == NULL) ||
       ((pIdxSortTb = (int *)AlcMalloc(sizeof(int) * nP)) == NULL) ||
       ((pIdxIdxTb = (int *)AlcMalloc(sizeof(int) * mN)) == NULL) ||
       ((pIdxSorted = (int *)AlcMalloc(sizeof(int) * nP)) == NULL) ||
       ((eIdxTb = (int *)AlcMalloc(sizeof(int) * mE)) == NULL) ||
       ((nIdxTb = (int *)AlcMalloc(sizeof(int) * mN)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Create return object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	dom;
    WlzValues	val;

    dom.cm2d5 = mesh;
    val.core = NULL;
    mapObj = WlzMakeMain(WLZ_CMESH_2D5, dom, val, NULL, NULL, &errNum);
  }
  /* Allocate indexed values. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 3;
    WlzValues	val;

    val.x = WlzMakeIndexedValues(mapObj, 1, &dim, WLZ_GREY_DOUBLE,
	                         WLZ_VALUE_ATTACH_NOD, &errNum);
    ixv = val.x;
    mapObj->values = WlzAssignValues(val, NULL);
  }
  /* Create the mapping. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE,
		idN;

    for(idN = 0; idN < nP; ++idN)
    {
      pIdxIdxTb[pIdx[idN]] = idN;    /* LUT from node index to given arrays. */
      pIdxSortTb[idN] = idN;  /* LUT for pinned nodes sorted by index value. */
    }
    (void )AlgHeapSortIdx(pIdx, pIdxSortTb, nP, AlgHeapSortCmpIdxIFn);
    for(idN = 0; idN < nP; ++idN)
    {
      pIdxSorted[idN] = pIdx[pIdxSortTb[idN]];
    }
    /* Fill in the element and node index tables. */
    (void )WlzCMeshSetNodIdxTbl2D5(mesh, nIdxTb);
    (void )WlzCMeshSetElmIdxTbl2D5(mesh, eIdxTb);
    /* Compute matrix values. */
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm2D5 *elm;

      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	int  	idN;
	double	a2;
	WlzDVertex3 v[3],
		    p[3];
	WlzCMeshNod2D5 *nod[3];

	/* Get element vertices. */
	nod[0] = WLZ_CMESH_ELM2D5_GET_NODE_0(elm); p[0] = nod[0]->pos;
	nod[1] = WLZ_CMESH_ELM2D5_GET_NODE_1(elm); p[1] = nod[1]->pos;
	nod[2] = WLZ_CMESH_ELM2D5_GET_NODE_2(elm); p[2] = nod[2]->pos;
	/* Compute useful vectors and lengths. */
	WLZ_VTX_3_SUB(v[0], p[1], p[0]);
	WLZ_VTX_3_SUB(v[1], p[2], p[0]);
	WLZ_VTX_3_CROSS(v[2], v[0], v[1]);
	a2 = WLZ_VTX_3_LENGTH(v[2]);
	if(a2 < WLZ_MESH_TOLERANCE)
	{
	  /* This should never occur! */
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  break;
	}
	else
	{
	  int    idT;
	  double d,
		 l0,
		 l2;
	  double wR[3],
		 wI[3];
	  WlzDVertex2 q2;
	  WlzDVertex3 u[3]; /* Basis vectors within the plane of the current
			       triangular element. */

	  idT = eIdxTb[elm->idx];
	  d = 1.0 / sqrt(a2);
	  l0 = WLZ_VTX_3_LENGTH(v[0]);
	  l2 = a2; /* WLZ_VTX_3_LENGTH(v[2]) */
	  /* Compute the orthonormal basis vectors for this element. */
	  WLZ_VTX_3_SCALE(u[0], v[0], 1.0 / l0);
	  WLZ_VTX_3_SCALE(u[2], v[2], 1.0 / l2);
	  WLZ_VTX_3_CROSS(u[1], u[2], u[0]);
	  q2.vtX = WLZ_VTX_3_DOT(v[1], u[0]);
	  q2.vtY = WLZ_VTX_3_DOT(v[1], u[1]);
	  wR[0] = d * (q2.vtX - l0);
	  wI[0] = d * q2.vtY; 
	  wR[1] = d * -q2.vtX;
	  wI[1] = d * -q2.vtY;
	  wR[2] = d * l0;
	  wI[2] = 0.0;
	  for(idN = 0; idN < 3; ++idN)
	  {
	    int idV;
	    int *idPP;

	    idV = nIdxTb[nod[idN]->idx];
	    if((idPP = bsearch(&(nod[idN]->idx), pIdxSorted, nP, sizeof(int),
		    WlzCMeshSurfMapIdxCmpFn)) == NULL)
	    {
	      /* Node is free. */
	      (void )AlgMatrixSet(aM, idT,      idV,       wR[idN]);
	      (void )AlgMatrixSet(aM, idT + nE, idV,      -wI[idN]);
	      (void )AlgMatrixSet(aM, idT,      idV + nN,  wI[idN]);
	      (void )AlgMatrixSet(aM, idT + nE, idV + nN,  wR[idN]);
	    }
	    else
	    {
	      int idQ,
	      idP;

	      /* Node is pinned. */
	      idQ = idPP - pIdxSorted;    /* Index into table pinned node. */
	      idP = pIdxIdxTb[nod[idN]->idx];
	      (void )AlgMatrixSet(bPM, idT,      idQ,       wR[idN]);
	      (void )AlgMatrixSet(bPM, idT + nE, idQ,      -wI[idN]);
	      (void )AlgMatrixSet(bPM, idT,      idQ + nP,  wI[idN]);
	      (void )AlgMatrixSet(bPM, idT + nE, idQ + nP,  wR[idN]);
	      bUV[idQ     ] = dPV[idP].vtX;
	      bUV[idQ + nP] = dPV[idP].vtY;
	    }
	  }
	}
      }
    }
    /* Compute bV and solve for mapped vertices. */
    if(errNum == WLZ_ERR_NONE)
    {
      AlgMatrixVectorMul(bV, bPM, bUV);
#ifdef WLZ_CMESH_SM_DEBUG
      (void )fprintf(stderr, "WlzCMeshCompSurfMapLevy() aM =\n[\n");
      AlgMatrixWriteAscii(aM, stderr);
      (void )fprintf(stderr, "]\n");
      (void )fprintf(stderr, "WlzCMeshCompSurfMapLevy() bV =\n[\n");
      (void )AlcDouble1WriteAsci(stderr, bV, nE * 2);
      (void )fprintf(stderr, "]\n");
#endif 
      errNum = WlzErrorFromAlg(
	  AlgMatrixSolveLSQR(aM, bV, xV, 1.0e-3, 1.0e-9, 1.0e-9, 10000, 0,
	    NULL, NULL, NULL, NULL, NULL, NULL, NULL));
    }
    if(errNum == WLZ_ERR_NONE)
    {
      int		idN;
      double 	*dsp;

      /* Set displacements for the indexed values. */
      for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
      {
	WlzCMeshNod2D5 *nod;

	nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	if(nod->idx >= 0)
	{
	  int	idV;

	  dsp = (double *)WlzIndexedValueGet(ixv, nod->idx);
	  idV = nIdxTb[nod->idx];
	  if(bsearch(&(nod->idx), pIdxSorted, nP, sizeof(int),
		WlzCMeshSurfMapIdxCmpFn) == NULL)
	  {
	    /* Node is free. */
	    dsp[0] = -(xV[idV     ] + nod->pos.vtX);
	    dsp[1] = -(xV[idV + nN] + nod->pos.vtY);
	  }
	  else
	  {
	    /* Node is pinned. */
	    WlzDVertex3 *pV;

	    pV = dPV + pIdxIdxTb[idV];
	    dsp[0] = pV->vtX - nod->pos.vtX;
	    dsp[1] = pV->vtY - nod->pos.vtY;
	  }
	  dsp[2] = -(nod->pos.vtZ);
	}
      }
    }
  }
  AlcFree(xV);
  AlcFree(bV);
  AlcFree(bUV);
  AlcFree(nIdxTb);
  AlcFree(eIdxTb);
  AlgMatrixFree(aM);
  AlgMatrixFree(bPM);
  AlcFree(pIdxIdxTb);
  AlcFree(pIdxSortTb);
  AlcFree(pIdxSorted);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mapObj);
}

/*!
* \return	Comparison value for AlgSort().
* \ingroup	WlzMesh
* \brief	Called by AlgSort() to sort integers into acending order.
* \param	p0			Pointer to first int.
* \param	p1			Pointer to second int.
*/
static int	WlzCMeshSurfMapIdxCmpFn(const void *p0, const void *p1)
{
  int		*i0,
  		*i1;

  i0 = (int *)p0;
  i1 = (int *)p1;
  return(*i0 - *i1);
}

/*!
* \return       New contour object corresponding to the given mesh.
* \ingroup      WlzMesh
* \brief        Creates a contour corresponding to the given conforming mesh
*               which must be a 2D5 mesh, ie a surface.
* \param        mObj                    Given conforming mesh.
* \param        disp                    Scale factor for the displacements, 0.0
* 					implies no displacements, 1.0 implies
* 					full displacement.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToContour(WlzObject *mObj, double disp,
                                   WlzErrorNum *dstErr)
{
  WlzDomain dom;
  WlzValues val;
  WlzObject     *cObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    dom.ctr = WlzMakeContour(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr->model = WlzAssignGMModel(WlzCMeshToGMModel(mObj, disp, &errNum),
    			              NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cObj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeContour(dom.ctr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cObj);
}

/*!
* \return       New geometric model corresponding to the given mesh.
* \ingroup      WlzMesh
* \brief        Creates a geometric model corresponding to the given
* 		conforming mesh which must be either a 2D or 2D5 mesh,
* 		ie a surface. The resulting model will have either
* 		have type WLZ_GMMOD_3D (from 2D5) or WLZ_GMMOD_2D (from
* 		2D).
* \param        mObj                    Given conforming mesh.
* \param        disp                    Scale factor for the displacements, 0.0
* 					implies no displacements, 1.0 implies
* 					full displacement.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzGMModel	*WlzCMeshToGMModel(WlzObject *mObj, double disp,
				   WlzErrorNum *dstErr)
{
  WlzGMModel	*model = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(disp != 0)
  {
    if(mObj->values.core == NULL)
    {
      errNum = WLZ_ERR_VALUES_NULL;
    }
    else if(mObj->values.core->type != (WlzObjectType )WLZ_INDEXED_VALUES)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(mObj->type)
    {
      case WLZ_CMESH_2D:
	model = WlzCMeshToGMModel2D(mObj, disp, &errNum);
        break;
      case WLZ_CMESH_2D5:
	model = WlzCMeshToGMModel2D5(mObj, disp, &errNum);
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
  return(model);
}

/*!
* \return	New 2D conforming mesh object.
* \ingroup	WlzMesh
* \brief	Creates a new 2D conforming mesh object by flattening the given
* 		2D5 conforming mesh object. This is done by applying the
* 		2D5 object's indexed values which are assumed to be valid
* 		displacements to a plane. See WlzCMeshCompSurfMap().
* \param	gObj			Given 2D5 conforming mesh object
* 					with valid displacements.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshFlatten2D5(WlzObject *gObj, WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*gMesh = NULL;
  WlzCMesh2D	*rMesh = NULL;
  WlzObject	*rObj = NULL;
  WlzIndexedValues *gIxv;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((gMesh = gObj->domain.cm2d5) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((gMesh->res.nod.numEnt < 3) || (gMesh->res.elm.numEnt < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if((gIxv = gObj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((gIxv->rank != 1) || (gIxv->dim[0] < 2))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    rMesh = WlzCMeshNew2D(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((AlcVectorExtendAndGet(rMesh->res.nod.vec,
                              gMesh->res.nod.maxEnt) == NULL) ||
       (AlcVectorExtendAndGet(rMesh->res.elm.vec,
                              gMesh->res.elm.maxEnt) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)                
  {
    int	idN;

    for(idN = 0; idN < gMesh->res.nod.maxEnt; ++idN)
    {
      WlzCMeshNod2D  *rNod;
      WlzCMeshNod2D5 *gNod;

      gNod = (WlzCMeshNod2D5 *)AlcVectorItemGet(gMesh->res.nod.vec, idN);
      if(gNod->idx >= 0)
      {
	double	*dsp;
	WlzDVertex2 pos;

	dsp = (double *)WlzIndexedValueGet(gIxv, gNod->idx);
	pos.vtX = gNod->pos.vtX + dsp[0];
	pos.vtY = gNod->pos.vtY + dsp[1];
        rNod = WlzCMeshNewNod2D(rMesh, pos, NULL);
      }
      else
      {
	/* Insert matching invalid mesh node. */
        rNod = (WlzCMeshNod2D *)AlcVectorItemGet(rMesh->res.nod.vec, idN);
	rNod->idx = -1;
	++(rMesh->res.nod.nextIdx);
      }
    }
    WlzCMeshUpdateBBox2D(rMesh);
    errNum = WlzCMeshReassignGridCells2D(rMesh, rMesh->res.nod.numEnt);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int	idE;

    for(idE = 0; idE < gMesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm2D *rElm;
      WlzCMeshElm2D5 *gElm;

      gElm = (WlzCMeshElm2D5 *)AlcVectorItemGet(gMesh->res.elm.vec, idE);
      if(gElm->idx >= 0)
      {
	int	idN;
        WlzCMeshNod2D *nod[3];

	for(idN = 0; idN < 3; ++idN)
	{
	  nod[idN] = (WlzCMeshNod2D *)AlcVectorItemGet(rMesh->res.nod.vec,
						       gElm->edu[idN].nod->idx);
	}
	(void )WlzCMeshNewElm2D(rMesh, nod[0], nod[1], nod[2], 1, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
      else
      {
	/* Insert matching invalid mesh element. */
        rElm = (WlzCMeshElm2D *)AlcVectorItemGet(rMesh->res.elm.vec, idE);
	rElm->idx = -1;
	++(rMesh->res.elm.nextIdx);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshReassignGridCells2D(rMesh, rMesh->res.nod.numEnt);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain dom;
    WlzValues val;

    dom.cm2 = rMesh;
    val.core = NULL;
    WlzCMeshUpdateMaxSqEdgLen2D(rMesh);
    rObj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(rObj != NULL)
    {
      (void )WlzFreeObj(rObj);
      rObj = NULL;
    }
    else if(rMesh != NULL)
    {
      (void )WlzCMeshFree2D(rMesh);
      rMesh = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return       New geometric model corresponding to the given mesh.
* \ingroup      WlzMesh
* \brief        Creates a geometric model corresponding to the given
* 		conforming mesh which is assumed to be a 2D mesh.
* 		The resulting model will be a WLZ_GMMOD_2D model.
* \param        mObj                    Given conforming mesh.
* \param        disp                    Scale factor for the displacements, 0.0
* 					implies no displacements, 1.0 implies
* 					full displacement.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzGMModel *WlzCMeshToGMModel2D(WlzObject *mObj, double disp,
				       WlzErrorNum *dstErr)
{
  int		nBkSz,
  		nHTSz,
		useDisp,
		nNod = 0;
  WlzCMesh2D	*mesh;
  WlzGMModel	*model = NULL;
  WlzIndexedValues *ixv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minBkSz = 1024,
  		minHTSz = 1024;

  mesh = mObj->domain.cm2;
  useDisp = (fabs(disp) > WLZ_MESH_TOLERANCE)? 1: 0;
  if(((nNod = mesh->res.nod.numEnt) < 3) || (mesh->res.elm.numEnt < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if(useDisp != 0)
  {
    if((ixv = mObj->values.x) == NULL)
    {
      errNum = WLZ_ERR_VALUES_NULL;
    }
    else if((ixv->rank != 1) || (ixv->dim[0] < 2))
    {
      errNum = WLZ_ERR_VALUES_DATA;
    }
  }
  /* Create a new geometric model. */
  if(errNum == WLZ_ERR_NONE)        
  {
    if((nBkSz = nNod / 16) < minBkSz)
    {
      nBkSz = minBkSz;
    }
    if((nHTSz = nNod / 4) < minHTSz)
    {
      nHTSz = minHTSz;
    }
    model = WlzGMModelNew(WLZ_GMMOD_2D, nBkSz, nHTSz, &errNum);
  }
  /* Add the simplices (edges) to the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm2D *elm;
      
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	int	idN;
	WlzDVertex2 pos[4];

	for(idN = 0; idN < 3; ++idN)
	{
	  pos[idN] = elm->edu[idN].nod->pos;
	  if(useDisp)
	  {
	    double *dsp;
	    
	    dsp = (double *)WlzIndexedValueGet(ixv, elm->edu[idN].nod->idx);
	    pos[idN].vtX += disp * dsp[0];
	    pos[idN].vtY += disp * dsp[1];
	  }
	}
	pos[3] = pos[0];
	for(idN = 0; idN < 3; ++idN)  
	{
	  errNum = WlzGMModelConstructSimplex2D(model, pos + idN);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelRehashVHT(model, 0);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFree(model);
    model = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
}

/*!
* \return       New geometric model corresponding to the given mesh.
* \ingroup      WlzMesh
* \brief        Creates a geometric model corresponding to the given
* 		conforming mesh which is assumed to be a 2D5 mesh.
* 		The resulting model will be a WLZ_GMMOD_3D model.
* \param        mObj                    Given conforming mesh.
* \param        disp                    Scale factor for the displacements, 0.0
* 					implies no displacements, 1.0 implies
* 					full displacement.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
static WlzGMModel *WlzCMeshToGMModel2D5(WlzObject *mObj, double disp,
				        WlzErrorNum *dstErr)
{
  int		nBkSz,
  		nHTSz,
		useDisp,
		nNod = 0;
  WlzCMesh2D5	*mesh;
  WlzGMModel	*model = NULL;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	minBkSz = 1024,
  		minHTSz = 1024;

  ixv = mObj->values.x;
  mesh = mObj->domain.cm2d5;
  useDisp = (fabs(disp) > WLZ_MESH_TOLERANCE)? 1: 0;
  if(((nNod = mesh->res.nod.numEnt) < 3) || (mesh->res.elm.numEnt < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else if(ixv == NULL)
  {
    useDisp = 0;
  }
  else if((ixv->rank != 1) || (ixv->dim[0] < 3))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  /* Create a new geometric model. */
  if(errNum == WLZ_ERR_NONE)        
  {
    if((nBkSz = nNod / 16) < minBkSz)
    {
      nBkSz = minBkSz;
    }
    if((nHTSz = nNod / 4) < minHTSz)
    {
      nHTSz = minHTSz;
    }
    model = WlzGMModelNew(WLZ_GMMOD_3D, nBkSz, nHTSz, &errNum);
  }
  /* Add the simplices (faces) to the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm2D5 *elm;
      
      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	int	idN;
	WlzDVertex3 pos[3];

	for(idN = 0; idN < 3; ++idN)
	{
	  pos[idN] = elm->edu[idN].nod->pos;
	}
	if(useDisp)
	{
	  for(idN = 0; idN < 3; ++idN)
	  {
	    double *dsp;
	    
	    dsp = (double *)WlzIndexedValueGet(ixv, elm->edu[idN].nod->idx);
	    pos[idN].vtX += disp * dsp[0];
	    pos[idN].vtY += disp * dsp[1];
	    pos[idN].vtZ += disp * dsp[2];
	  }
	}
        errNum = WlzGMModelConstructSimplex3D(model, pos);
	if(errNum != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelRehashVHT(model, 0);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzGMModelFree(model);
    model = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
}
