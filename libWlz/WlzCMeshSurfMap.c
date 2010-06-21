#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshSurfMap_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshSurfMap.c
* \author       Bill Hill
* \date         May 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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

static int			WlzCMeshSurfMapIdxCmpFn(
				  const void *p0,
				  const void *p1);

/*!
* \return	New mesh object with displacements set or NULL on error.
* \ingroup	WlzTransform
* \brief	Computes a least squares conformal transformation which
* 		maps the source surface to a destination plane with z
* 		coordinate zero. See WlzCMeshCompSurfMapConformalIdx().
* \param	inObj			Input conforming mesh object which
* 					must be of type WLZ_CMESH_2D5.
* \param	nDV			Number of destination vertices.
* \param	dV			Destination vertices.
* \param	nSV			Number of destination vertices
* 					which must be the same as nDV.
* \param	sV			Source vertices.
* \param	dstErr			Woolz error code, may be NULL.
*/
WlzObject	*WlzCMeshCompSurfMapConformal(WlzObject *inObj,
				int nDV, WlzDVertex3 *dV,
				int nSV, WlzDVertex3 *sV,
				WlzErrorNum *dstErr)
{
  int		*nodTb = NULL;
  WlzCMesh2D5	*mesh;
  WlzValues	val;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
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
      if(WlzCMeshMatchNNodIdx2D5(mesh, nSV, sV, nodTb) != nDV)
      {
        errNum = WLZ_ERR_DOMAIN_DATA;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rtnObj = WlzCMeshCompSurfMapConformalIdx(mesh, nDV, dV, nodTb, &errNum);
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
* 		TODO
* \param	mesh			Input conforming mesh which must be
* 					of type WLZ_CMESH_2D5.
* \param	nPN			Number of pinned nodes.
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
WlzObject	*WlzCMeshCompSurfMapConformalIdx(WlzCMesh2D5 *mesh,
				int nP, WlzDVertex3 *dPV, int *pIdx,
				WlzErrorNum *dstErr)
{
  int		nE,
		nE2,
		nN,
		nN2,
  		nF,
		nF2,
		nP2;
  int		*pIdxSortTb = NULL,
  		*eIdxTb = NULL,
  		*nIdxTb = NULL,
		*pIdxSorted = NULL;
  double	*bM = NULL,
  		*bUM = NULL,
		*xM = NULL;
  double	**aM = NULL,
  		**bPM = NULL;
  WlzObject	*rtnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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
    nN = mesh->res.nod.numEnt;
    nE = mesh->res.elm.numEnt;
    nF = nN - nP;
    nE2 = nE * 2;
    nN2 = nN * 2;
    nF2 = nF * 2;
    nP2 = nP * 2;
    if((AlcDouble2Calloc(&aM, nE2, nN2) != ALC_ER_NONE) ||
       (AlcDouble2Calloc(&bPM, nE2, nP2) != ALC_ER_NONE) ||
       ((bUM = AlcMalloc(sizeof(double) * nP2)) == NULL) ||
       ((bM = AlcMalloc(sizeof(double) * nE2)) == NULL) ||
       ((xM = AlcMalloc(sizeof(double) * nE2)) == NULL) ||
       ((pIdxSortTb = AlcMalloc(sizeof(int) * nP)) == NULL) ||
       ((pIdxSorted = AlcMalloc(sizeof(int) * nP)) == NULL) ||
       ((eIdxTb = AlcMalloc(sizeof(int) * mesh->res.elm.maxEnt)) == NULL) ||
       ((nIdxTb = AlcMalloc(sizeof(int) * mesh->res.nod.maxEnt)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;

    /* Sort the pinned nodes by index value. */
    for(idE = 0; idE < nP; ++idE)
    {
      pIdxSortTb[idE] = idE;
    }
    (void )AlgHeapSortIdx(pIdx, pIdxSortTb, nP, AlgHeapSortCmpIdxIFn);
    for(idE = 0; idE < nP; ++idE)
    {
      pIdxSorted[idE] = pIdx[pIdxSortTb[idE]];
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
	WlzDVertex2 v2;
	WlzDVertex3 v[3],
		    p[3];
	WlzCMeshNod2D5 *nod[3];

	/* Get element vertives. */
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
	  int    idT,
	  	 idT2;
	  double l0,
	  	 l2;
	  double wR[3],
		 wI[3];
	  WlzDVertex3 u[3]; /* Basis vectors within the plane of the current
	                       triangular element. */

	  a2 = 1.0 / sqrt(a2);
	  l0 = WLZ_VTX_3_LENGTH(v[0]);
	  l2 = WLZ_VTX_3_LENGTH(v[2]);
	  idT = eIdxTb[elm->idx];
	  idT2 = 2 * idT;
	  /* Compute the orthonormal basis vectors for this element. */
	  WLZ_VTX_3_SCALE(u[0], v[0], 1.0 / l0);
	  WLZ_VTX_3_SCALE(u[2], v[2], 1.0 / l2);
	  WLZ_VTX_3_CROSS(u[1], u[2], u[0]);
	  v2.vtX = a2 * WLZ_VTX_3_DOT(v[1], u[0]);
	  v2.vtY = a2 * WLZ_VTX_3_DOT(v[1], u[1]);
	  l0 *= a2;
	  wR[0] = v2.vtX - l0;
	  wI[0] = v2.vtY; 
	  wR[1] = -v2.vtX;
	  wI[1] = -v2.vtY;
	  wR[2] = l0;
	  wI[2] = 0.0;
	  for(idN = 0; idN < 3; ++idN)
	  {
	    int idV,
	    	idV2;
	    int *idPP;

	    idV = nIdxTb[nod[idN]->idx];
	    idV2 = idV * 2;
	    if((idPP = bsearch(&(nod[idN]->idx), pIdxSorted, nP, sizeof(int),
	                       WlzCMeshSurfMapIdxCmpFn)) != NULL)
	    {
	      int idP,
	          idP2;

	      /* Node is pinned. */
	      idP = pIdxSortTb[idPP - pIdxSorted];
	      idP2 = 2 * idP;
	      bPM[idT2    ][idP2    ] =  wR[idN];
	      bPM[idT2    ][idP2 + 1] = -wI[idN];
	      bPM[idT2 + 1][idP2    ] =  wI[idN];
	      bPM[idT2 + 1][idP2 + 1] =  wR[idN];
	      bUM[idP2    ] = dPV[idP].vtX;
	      bUM[idP2 + 1] = dPV[idP].vtY;
	    }
	    else
	    {
	      /* Node is free. */
	      aM[idT2    ][idV2    ] =  wR[idN];
	      aM[idT2    ][idV2 + 1] = -wI[idN];
	      aM[idT2 + 1][idV2    ] =  wI[idN];
	      aM[idT2 + 1][idV2 + 1] =  wR[idN];
	    }
	  }
	}
      }
    }
  }
  AlcFree(pIdxSortTb);
  /* Compute bM and solve for mapped vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    AlgMatrixVectorMul(bM, ALG_MATRIX_RECT, bPM, bUM, nE2, nP2);
    /* AlgMatrixScale(&bM, &bM, -1.0, 1, nE2); */
#ifdef WLZ_CMESH_SM_DEBUG
    (void )fprintf(stderr, "WlzCMeshCompSurfMapConformalIdx() aM =\n[\n");
    (void )AlcDouble2WriteAsci(stderr, aM, nE2, nN2);
    (void )fprintf(stderr, "]\n");
    (void )fprintf(stderr, "WlzCMeshCompSurfMapConformalIdx() bM =\n[\n");
    (void )AlcDouble1WriteAsci(stderr, bM, nE2);
    (void )fprintf(stderr, "]\n");
#endif 
    errNum = WlzErrorFromAlg(
	     AlgMatrixSolveLSQR(ALG_MATRIX_RECT, aM, nE2, nN2, bM, xM,
	                        0.0, 1.0e-9, 1.0e-9, 1000, 0,
				NULL, NULL, NULL, NULL, NULL, NULL, NULL));
  }
  AlcFree(bUM);
  AlcFree(eIdxTb);
  Alc2Free((void **)bPM);
  Alc2Free((void **)aM);
  /* Create return object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	dom;
    WlzValues	val;

#ifdef WLZ_CMESH_SM_DEBUG
    (void )fprintf(stderr, "WlzCMeshCompSurfMapConformalIdx() xM =\n[\n");
    (void )AlcDouble1WriteAsci(stderr, xM, nN2);
    (void )fprintf(stderr, "]\n");
#endif 
    dom.cm2d5 = mesh;
    val.core = NULL;
    rtnObj = WlzMakeMain(WLZ_CMESH_2D5, dom, val, NULL, NULL, &errNum);
  }
  /* Allocate indexed values. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim = 3;
    WlzValues	val;

    val.x = WlzMakeIndexedValues(rtnObj, 1, &dim, WLZ_GREY_DOUBLE,
                                 WLZ_VALUE_ATTACH_NOD, &errNum);
    rtnObj->values = WlzAssignValues(val, NULL);
  }
  /* Set displacements for the indexed values. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;
    double 	*dsp;
    WlzIndexedValues *ixv;

    ixv = rtnObj->values.x;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      WlzCMeshNod2D5 *nod;

      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	int	idV,
	        idV2;

	dsp = (double *)WlzIndexedValueGet(ixv, nod->idx);
	idV = nIdxTb[nod->idx];
	if(bsearch(&(nod->idx), pIdxSorted, nP, sizeof(int),
		   WlzCMeshSurfMapIdxCmpFn) == NULL)
	{
	  dsp[0] = 0.0;
	  dsp[1] = 0.0;
	  dsp[2] = 0.0;
	}
	else
	{
	  idV2 = 2 * idV;
	  dsp[0] = xM[idV2    ] - nod->pos.vtX;
	  dsp[1] = xM[idV2 + 1] - nod->pos.vtY;
	  dsp[3] = 0.0;
	}
      }
    }
  }
  AlcFree(xM);
  AlcFree(bM);
  AlcFree(nIdxTb);
  AlcFree(pIdxSorted);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rtnObj);
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
* \param        disp                    Non zero if the mesh displacements
*                                       are to be applied.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject       *WlzCMeshToContour(WlzObject *mObj, int disp,
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
  else if(mObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
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
* \return       New contour object corresponding to the given mesh.
* \ingroup      WlzMesh
* \brief        Creates a contour corresponding to the given conforming mesh
*               which must be a 2D5 mesh, ie a surface.
* \param        mObj                    Given conforming mesh.
* \param        disp                    Non zero if the mesh displacements
*                                       are to be applied.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzGMModel	*WlzCMeshToGMModel(WlzObject *mObj, int disp,
				   WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*mesh;
  WlzGMModel	*model = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(mObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if((mesh = mObj->domain.cm2d5) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh->res.elm.numEnt < 1) || (mesh->res.nod.numEnt < 3))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    /* Create a new geometric model. */
    /* Create the vertices of the model. */
    /* Create the loops in the model. */

    /* TODO */
    errNum = WLZ_ERR_UNIMPLEMENTED;
    
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
}

/*!
* \return	New 2D domain object with values.
* \ingroup	WlzMesh
* \brief	Creates a 2D domain object with values corresponding to
* 		the Gaussian curvature of the given mesh object with
* 		displacements to a plane.
* \param	mObj			Given conforming mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzCMeshCurvToImage(WlzObject *mObj, WlzErrorNum *dstErr)
{
  WlzObject	*iObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* TODO */
  errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iObj);
}
