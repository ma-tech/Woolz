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
* 		on conformal transformations.
* \ingroup	WlzTransform
*/

#include<Wlz.h>
#include <stdlib.h>
#include <float.h>

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

/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a least squares conformal transformation which
* 		maps the mesh surface to the the given planar spatial domain.
* 		The boundaries of the mesh and spatial domain are
* 		considered to be equivalent and the given vertices must
* 		lie on these boundaries. This function can give invalid
* 		meshes with flipped/crossing elements when the given
* 		domain or mesh boundary are non-convex.
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
  		nMshBN = 0;     /* Number of boundary mesh nodes. */
  int	 	*mshNI = NULL,  /* Given point indices of nodes in mesh. */
		*prmI  = NULL,  /* Permutation to make given points cyclic. */
		*pinNI = NULL,  /* Pin node indices. */
  		*plyVI = NULL,  /* Given point indices in boundary polygon. */
		*mshBNI = NULL; /* Mesh boundary node indices. */
  AlcVector	*mshNV;         /* Mesh node vector. */
  WlzIVertex2	*plyV;
  WlzDVertex2	*mshV2 = NULL,
  		*plyV2 = NULL;
  WlzDVertex3	*pinV = NULL;
  WlzBasisFnTransform *bFn = NULL;
  WlzPolygonDomain *ply = NULL;
  WlzObject	*prmObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;      
  WlzCMesh2D5	*mesh;

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

    for(i = 0; i < nMV; ++i)
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
	WlzDVertex2 tV;
        WlzCMeshNod2D5 *nod;

	/* All nodes will be valid on boundary. */
	nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[j]);
	WLZ_VTX_2_SUB(tV, nod->pos, gV);
	d1 = WLZ_VTX_2_SQRLEN(tV);
	if(d1 < d)
	{
	  d = d1;
	  mshNI[i] = j;
	}
      }
    }
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
  /* Allocate vertex and mesh node pin arrays such that theye will always have
   * less elements that the number allocated. */
  if(errNum == WLZ_ERR_NONE)
  {

    nPin = ALG_MAX(ply->nvertices, nMshBN);
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
      int	j,
      		j0,     /* Cyclic ordered index of 1st vertex/node. */
      		j1,    	/* Cyclic ordered index of 2nd vertex/node. */
		k,
		nSPV, 	/* Number of polygon vertices in segment. */
		nSMN; 	/* Number of mesh nodes in segment. */
      double	d,
      		d3;     /* Distance from 1st to 2nd node of mesh. */
      WlzCMeshNod2D5 *nod[2];

      j0 = prmI[i];
      j1 = prmI[(i + 1) % nMV];
      nSPV = ((plyVI[j1] - plyVI[j0] + ply->nvertices) % ply->nvertices) + 1;
      nSMN = ((mshNI[j1] - mshNI[j0] + nMshBN) % nMshBN) + 1;
      /* Don't need to find the polygon segment distance from 1st to 2nd
       * given vertex of the polygon because they're equispaced. *
       * Find distance from 1st to 2nd given node of mesh boundary. */
      d3 = 0.0;
      k = mshNI[j0];
      nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[k]);
      for(j = 1; j < nSMN; ++j)
      {
	WlzDVertex3 t3;

	k = (k + 1) % nMshBN;
	nod[0] = nod[1];
        nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[k]);
        WLZ_VTX_3_SUB(t3, nod[0]->pos, nod[1]->pos);
	d3 += WLZ_VTX_3_LENGTH(t3);
      }
      /* Set the pins using the mesh boundary nodes and the corresponding
       * spatial domain boundary polygon vertices. */
      d = 0.0;
      k = mshNI[j0];                                 
      pinV[nPin].vtX = plyV[plyVI[j0]].vtX;
      pinV[nPin].vtY = plyV[plyVI[j0]].vtY;
      pinV[nPin].vtZ = 0.0;
      pinNI[nPin] = mshBNI[k];
      ++nPin;
      nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[k]);
      for(j = 1; j < nSMN; ++j)
      {
	int	    p;
	WlzDVertex3 t3;

	k = (k + 1) % nMshBN;
	nod[0] = nod[1];
        nod[1] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mshNV, mshBNI[k]);
        WLZ_VTX_3_SUB(t3, nod[0]->pos, nod[1]->pos);
	d += WLZ_VTX_3_LENGTH(t3);
	p = (plyVI[j0] + ALG_NINT(nSPV * (d / d3))) % nPlyBV;
	pinV[nPin].vtX = plyV[p].vtX;
	pinV[nPin].vtY = plyV[p].vtY;
	pinV[nPin].vtZ = 0.0;
        pinNI[nPin] = mshBNI[k];
        ++nPin;
      }
    }
    if(pinNI[nPin - 1] == pinNI[0])   /* Avoid duplicated start node at end. */
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
  AlcFree(mshNI);
  AlcFree(plyVI);
  (void )WlzFreePolyDmn(ply);
  /* Compute the parameterised mesh using the pinned nodes and their
   * parameterised positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    prmObj = WlzCMeshCompSurfMapIdx(mesh, nPin, pinV, pinNI, &errNum);
  }
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
* 		are used. See WlzCMeshCompSurfMapIdx().
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
  WlzCMesh2D5	*mesh;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(inObj->domain.core->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((nDV < 1) || (nDV != nSV) || (dV == NULL) || (sV == NULL))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    mesh = inObj->domain.cm2d5;
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
    rtnObj = WlzCMeshCompSurfMapIdx(mesh, nDV, dV, nodTb, &errNum);
  }
  AlcFree(nodTb);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
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
WlzObject	*WlzCMeshCompSurfMapIdx(WlzCMesh2D5 *mesh,
				int nP, WlzDVertex3 *dPV, int *pIdx,
				 WlzErrorNum *dstErr)
{
  int		mE,
  		nE,
		nE2,
		nP2,
		mN,
		nN,
		nN2;
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
      (void )fprintf(stderr, "WlzCMeshCompSurfMapIdx() aM =\n[\n");
      AlgMatrixWriteAscii(aM, stderr);
      (void )fprintf(stderr, "]\n");
      (void )fprintf(stderr, "WlzCMeshCompSurfMapIdx() bV =\n[\n");
      (void )AlcDouble1WriteAsci(stderr, bV, nE * 2);
      (void )fprintf(stderr, "]\n");
#endif 
      errNum = WlzErrorFromAlg(
	  AlgMatrixSolveLSQR(aM, bV, xV, 1.0e-10, 1.0e-10, 1.0e-10, 10000, 0,
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
* \return	Comparison value for qsort().
* \ingroup	WlzMesh
* \brief	Called by qsort() to sort integers into acending order.
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
    else if(mObj->values.core->type != WLZ_INDEXED_VALUES)
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
  if((gMesh->res.nod.numEnt < 3) || (gMesh->res.elm.numEnt < 1))
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
