#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTensor_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzTensor.c
* \author       Bill Hill
* \date         January 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Functions which derive and manipulate tensor quantities.
* \ingroup	WlzFeatures
*/
#include <stdio.h>
#include <Wlz.h>

static WlzObject		*WlzCMeshDGTensor3D(
				  WlzObject *cObj,
				  int invert,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzCMeshDGTensorAtPts3D(
				  WlzObject *cObj,
				  int invert,
				  WlzDVertex3 sd,
				  int dither,
				  WlzErrorNum *dstErr);
static void	 		WlzCMeshDGToStrainTensor3D(
				  WlzObject *tObj);
static void 			WlzPtsDGToStrainTensor3D(
				  WlzObject *tObj);
static void			WlzCMeshElmSetDGTensor3D(
				  WlzCMeshElm3D *elm,
				  int invert,
				  WlzIndexedValues *cIxv,
				  double *ten);

/*!
* \return	Conforming mesh object with tensor values attached to
* 		elements or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor for each of it's valid
* 		elements.
* 		Given displacement \f$\vec{u}(\vec{r})\f$ with position
* 		vector \f$\vec{r}\f$ which maps a point from a space
* 		\f$\vec{x}\f$ to a space \f$\vec{u}\f$ the displacement
* 		gradient tensor is defined in 3D as
  		\f[
		{
		\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}
		u_{i,j} =
		\left [
		\begin{array}{ccc}
		  \pd{u_0}{x_0} & \pd{u_0}{x_1} & \pd{u_0}{x_2} \\
		  \pd{u_1}{x_0} & \pd{u_1}{x_1} & \pd{u_1}{x_2} \\
		  \pd{u_2}{x_0} & \pd{u_2}{x_1} & \pd{u_2}{x_2}
		\end{array}
		\right]
		}
  		\f]
*		with
  		\f[
		\Delta r_i = u_{ij} r_j
  		\f]
*		where \f$\vec{u} = \left[u_0, u_1, u_2\right]^T\f$
*		and \f$\Delta \vec{r} = \left[x_0, x_1, x_2\right]^T\f$.
*		The displacement gradient tensor matrix is just the
*		rotation and independent scaling part of the  affine
*		transform that displaces the element.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDGTensor(WlzObject *cObj, int invert,
				  WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	tObj = WlzCMeshDGTensor3D(cObj, invert, &errNum);
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
  return(tObj);
}

/*!
* \return	Points object with tensor values at the point locations
* 		or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		The tensor values at the sample points are computed at
* 		each point by computing an iverse distance weighted
* 		least squares general affine transform for the ring of
* 		nodes surrounding the closest node.
* 		See WlzCMeshDGTensor() for the description of the tensor.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshDGTensorAtPts(WlzObject *cObj, int invert,
				       WlzDVertex3 sd, int dither,
				       WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if((sd.vtX < WLZ_MESH_TOLERANCE) ||
          (sd.vtY < WLZ_MESH_TOLERANCE) ||
          (sd.vtZ < WLZ_MESH_TOLERANCE))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	tObj = WlzCMeshDGTensorAtPts3D(cObj, invert, sd, dither, &errNum);
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
  return(tObj);
}

/*!
* \return	Conforming mesh object with tensor values attached to
* 		elements or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the strain tensor for each of it's valid elements.
* 		This function uses WlzCMeshDGTensor() to compute the
* 		displacement gradient tensor and the derives the strain
* 		tensor \f$e_{ij}\f$ from this using:
		\f[
		e_{ij} = \frac{1}{2} (u_{ij} + u_{ji})
		\f]
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshStrainTensor(WlzObject *cObj, int invert,
 				      WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum  = WLZ_ERR_NONE;

  tObj = WlzCMeshDGTensor(cObj, invert, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	WlzCMeshDGToStrainTensor3D(tObj);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Points object with tensor values or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		The tensor values at the sample points are computed using
* 		WlzCMeshDGTensorAtPts(). The strain tensor is then
* 		computed from the displacement gradient tensor as for
* 		WlzCMeshStrainTensor().
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshStrainTensorAtPts(WlzObject *cObj, int invert,
					   WlzDVertex3 sd,
					   int dither,
					   WlzErrorNum *dstErr)
{
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum  = WLZ_ERR_NONE;

  tObj = WlzCMeshDGTensorAtPts(cObj, invert, sd, dither, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    switch(cObj->type)
    {
      case WLZ_CMESH_2D:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      case WLZ_CMESH_3D:
	WlzPtsDGToStrainTensor3D(tObj);
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	New Woolz point object or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform this function computes
* 		the displacement gradient tensor for each of it's valid
* 		elements. See WlzCMeshDGTensor().
* \param	cObj			Given 3D conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshDGTensor3D(WlzObject *cObj, int invert,
				     WlzErrorNum *dstErr)
{
  WlzObject	*iObj = NULL,
  		*tObj = NULL;
  WlzCMesh3D	*mesh = NULL;
  WlzIndexedValues *cIxv = NULL,
  		   *tIxv = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mesh = cObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((cIxv = cObj->values.x)->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((cIxv->rank < 1) || (cIxv->dim[0] < 3) ||
          (cIxv->vType != WLZ_GREY_DOUBLE) ||
	  (cIxv->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    if(invert)
    {
      iObj = WlzCMeshTransformInvert(cObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	cObj = iObj;
        mesh = cObj->domain.cm3;
	cIxv = cObj->values.x;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		dim[2];

    dim[0] = dim[1] = 3;
    tIxv = WlzMakeIndexedValues(cObj, 2, dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_ELM, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      WlzValues	tVal;

      tVal.x = tIxv;
      tObj = WlzMakeMain(WLZ_CMESH_3D, cObj->domain, tVal, NULL, NULL,
                         &errNum);
      if(errNum != WLZ_ERR_NONE)
      {
        (void )WlzFreeIndexedValues(tIxv);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idE,
    		maxElm;
    AlcVector	*elmVec;

    elmVec = mesh->res.elm.vec;
    maxElm = mesh->res.elm.maxEnt;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idE = 0; idE < maxElm; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
      if(elm->idx >= 0)
      {
	double	*ten;

	ten = (double *)WlzIndexedValueGet(tIxv, elm->idx);
	WlzCMeshElmSetDGTensor3D(elm, 0, cIxv, ten);
      }
    }
  }
  if(invert)
  {
    WlzFreeObj(iObj);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(tObj);
    tObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \return	Points object with tensor values at the point locations
* 		or NULL on error.
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform this function computes
* 		the displacement gradient tensor at regular cartesian
* 		grid sample points throughout the mesh.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGTensorAtPts() with all checks having been done.
* \param	cObj			Given conforming mesh object.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	sd			Sample distance.
* \param	dither			Dither the sample locations.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzCMeshDGTensorAtPts3D(WlzObject *cObj, int invert,
					  WlzDVertex3 sd, int dither,
					  WlzErrorNum *dstErr)
{
  WlzDomain	dom;
  WlzValues	val;
  WlzDVertex3	sOrg;
  WlzIVertex3	sNum;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzObject	*tObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if((mesh = cObj->domain.cm3)->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((ixv = cObj->values.x)->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv->rank < 1) || (ixv->dim[0] < 3) ||
          (ixv->vType != WLZ_GREY_DOUBLE) ||
	  (ixv->attach != WLZ_VALUE_ATTACH_NOD))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    WlzCMeshUpdateBBox3D(mesh);
    sOrg.vtX = ceil(mesh->bBox.xMin / sd.vtX) * sd.vtX;
    sOrg.vtY = ceil(mesh->bBox.yMin / sd.vtY) * sd.vtY;
    sOrg.vtZ = ceil(mesh->bBox.zMin / sd.vtZ) * sd.vtZ;
    sNum.vtX = floor((mesh->bBox.xMax - sOrg.vtX) / sd.vtX) + 1;
    sNum.vtY = floor((mesh->bBox.yMax - sOrg.vtY) / sd.vtY) + 1;
    sNum.vtZ = floor((mesh->bBox.zMax - sOrg.vtZ) / sd.vtZ) + 1;
    if((sNum.vtX <= 0) || (sNum.vtY <= 0) || (sNum.vtZ <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      int	maxPts;
      WlzVertexP dum;

      dum.v = NULL;
      maxPts = sNum.vtX * sNum.vtY * sNum.vtZ;
      dom.pts = WlzMakePoints(WLZ_POINTS_3D, 0, dum, maxPts, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
	int	dim[2];

	dim[0] = dim[1] = 3;
        val.pts = WlzMakePointValues(maxPts, 2, dim, WLZ_GREY_DOUBLE,
	                             &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj = WlzMakeMain(WLZ_POINTS, dom, val, NULL, NULL, &errNum);
      }
      else
      {
        (void )WlzFreeDomain(dom);
        (void )WlzFreePointValues(val.pts);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		pIdx;
    WlzIVertex3 sIdx;
    WlzDVertex3 sPos;

    pIdx = 0;
    AlgRandSeed(0);
    for(sIdx.vtZ = 0; sIdx.vtZ <sNum.vtZ; ++(sIdx.vtZ))
    {
      sPos.vtZ = sOrg.vtZ + (sIdx.vtZ * sd.vtZ);
      for(sIdx.vtY = 0; sIdx.vtY <sNum.vtY; ++(sIdx.vtY))
      {
        sPos.vtY = sOrg.vtY + (sIdx.vtY * sd.vtY);
	for(sIdx.vtX = 0; sIdx.vtX <sNum.vtX; ++(sIdx.vtX))
	{
	  int	eIdx;
	  WlzDVertex3 dPos;

          sPos.vtX = sOrg.vtX + (sIdx.vtX * sd.vtX);
	  dPos = sPos;
	  if(dither)
	  {
	    dPos.vtX = AlgRandZigNormal(sPos.vtX, sd.vtX * 0.2);
	    dPos.vtY = AlgRandZigNormal(sPos.vtY, sd.vtY * 0.2);
	    dPos.vtZ = AlgRandZigNormal(sPos.vtZ, sd.vtZ * 0.2);
	  }
	  /* If dPos is in an element of the mesh compute it's displacement
	   * gradient tensor. */
          eIdx = WlzCMeshElmEnclosingPos3D(mesh, -1,
	                                   dPos.vtX, dPos.vtY, dPos.vtZ,
					   0, NULL);
	  if(eIdx >= 0)
	  {
	    double	*ten;
	    WlzCMeshElm3D *elm;

	    elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, eIdx);
	    if(elm->idx >= 0)
	    {
	      dom.pts->points.d3[pIdx] = dPos;
	      ten = WlzPointValueGet(val.pts, pIdx);
	      WlzCMeshElmSetDGTensor3D(elm, invert, ixv, ten);
	      ++pIdx;
	    }
	  }
	}
      }
    }
    dom.pts->nPoints = pIdx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(tObj);
}

/*!
* \ingroup 	WlzFeatures
* \brief	Given a 3D conforming mesh transform with displacement gradient
* 		tensor values this function converts these values to strain
* 		tensor values.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGToStrainTensor() with all checks having been done.
* \param	cObj			Given 3D conforming mesh object.
*/
static void 	WlzCMeshDGToStrainTensor3D(WlzObject *tObj)
{
  int 		idE,
    		maxElm;
  AlcVector	*elmVec;
  WlzCMesh3D	*mesh = NULL;
  WlzIndexedValues *tIxv = NULL;

  mesh = tObj->domain.cm3;
  tIxv = tObj->values.x;
  elmVec = mesh->res.elm.vec;
  maxElm = mesh->res.elm.maxEnt;
  for(idE = 0; idE < maxElm; ++idE) /* TODO use OpenMP */
  {
    WlzCMeshElm3D *elm;

    elm = (WlzCMeshElm3D *)AlcVectorItemGet(elmVec, idE);
    if(elm->idx >= 0)
    {
      double	*ten;

      ten = (double *)WlzIndexedValueGet(tIxv, elm->idx);
      ten[1] = ten[3] = 0.5 * (ten[1] + ten[3]);
      ten[2] = ten[6] = 0.5 * (ten[2] + ten[6]);
      ten[5] = ten[7] = 0.5 * (ten[5] + ten[7]);
    }
  }
}

/*!
* \ingroup 	WlzFeatures
* \brief	Given a points object with 3D vertices and displacement
* 		gradient tensor values this function converts these values
* 		to strain tensor values.
* 		This is a static function which was written to be called by
* 		WlzCMeshDGToStrainTensor() with all checks having been done.
* \param	cObj			Given 3D conforming mesh object.
*/
static void 	WlzPtsDGToStrainTensor3D(WlzObject *tObj)
{
  int		idx,
  		nPts;
  WlzPoints	*pd;
  WlzPointValues *pv;

  pd = tObj->domain.pts;
  pv = tObj->values.pts;
  nPts = pd->nPoints;
  for(idx = 0; idx < nPts; ++idx) /* TODO use OpenMP */
  {
    double	*ten;

    ten = (double *)WlzPointValueGet(pv, idx);
    ten[1] = ten[3] = 0.5 * (ten[1] + ten[3]);
    ten[2] = ten[6] = 0.5 * (ten[2] + ten[6]);
    ten[5] = ten[7] = 0.5 * (ten[5] + ten[7]);
  }
}

/*!
* \ingroup	WlzFeatures
* \brief	Given a 3D conforming mesh node, an indexed values with
* 		displacements for the mesh and a tensor pointer; this
* 		function computes the displacement gradient tensor
* 		value.
* \param	elm			Given mesh element.
* \param	invert			Invert if non-zero, by default the
* 					tensors are computed for the
* 					transform from the mesh to the
* 					displaced mesh.
* \param	cIxv			Mesh displacment indexed values.
* \param	ten			Tensor pointer.
*/
static void	WlzCMeshElmSetDGTensor3D(WlzCMeshElm3D *elm,
					 int invert,
					 WlzIndexedValues *cIxv,
					 double *ten)
{
  int		idN;
  double	tr[16];
  WlzCMeshNod3D *nod[4];
  WlzDVertex3 sVx[4],
	      dVx[4];

  nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
  nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
  nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
  nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
  for(idN = 0; idN < 4; ++idN)
  {
    double	*dsp;

    sVx[idN] = nod[idN]->pos;
    dsp = (double *)WlzIndexedValueGet(cIxv, nod[idN]->idx);
    dVx[idN].vtX = sVx[idN].vtX + dsp[0];
    dVx[idN].vtY = sVx[idN].vtY + dsp[1];
    dVx[idN].vtZ = sVx[idN].vtZ + dsp[2];
  }
  if(invert)
  {
    WlzGeomTetraAffineSolve(tr, dVx, sVx, WLZ_MESH_TOLERANCE_SQ);
  }
  else
  {
    WlzGeomTetraAffineSolve(tr, sVx, dVx, WLZ_MESH_TOLERANCE_SQ);
  }
  ten[0] = tr[0]; ten[1] = tr[1]; ten[2] = tr[2]; 
  ten[3] = tr[4]; ten[4] = tr[5]; ten[5] = tr[6]; 
  ten[6] = tr[8]; ten[7] = tr[9]; ten[8] = tr[10]; 
}
