#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshCurvature_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCMeshCurvature.c
* \author       Bill Hill
* \date         September 2010
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
* \brief	Functions to compute curvatures on conforming simplical
* 		meshes.
* \ingroup	WlzMesh
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
* \return	Interpolated curvature.
* \ingroup	WlzMesh
* \brief	Interpolates the curvature at some position within the
* 		given mesh element using a distance weighting of the
* 		curvatures at the nodes of the element and those at
* 		thier (single) edge connected neighbours.
* 		Squares of distances are used rather than the distances
* 		themselves to decrease the relative contribution of
* 		distant nodes.
* \param	elm			Given element.
* \param	ixv			Node indexed values, used for
* 					curvatures.
* \param	elmNod			The nodes of the element.
* \param	meanCrv			Non-zero if mean curvature is required
* 					(otherwise Gaussian curvature).
* \param	dPos			Position within the element to be
* 					interpolated.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static double	WlzCMeshElmRingInterpolateNodCrvValue(WlzCMeshElm2D *elm,
					WlzIndexedValues *ixv,
					WlzCMeshNod2D **elmNod,
					int meanCrv, WlzDVertex2 dPos,
					WlzErrorNum *dstErr)
{
  int		i,
  		j,
		k,
		nCnt = 0;
  double	crv = 0.0,
  		sWgt = 0.0;
  /* int	nNI[3]; */
  double	*wgt = NULL;
  WlzCMeshNod2D **nod = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Count number of edge connected nodes that the given element has including
   * the nodes of the element itself (this is a safe over estimate). */
  nCnt = 3;
  for(i = 0; i < 3; ++i)
  {
    WlzCMeshEdgU2D *e0,
    		   *e1;

    ++nCnt;
    e0 = e1 = elmNod[i]->edu;
    do
    {
      ++nCnt;
      e1 = e1->nnxt;
    } while(e1 != e0);
  }
  if(((nod = (WlzCMeshNod2D **)
	     AlcMalloc(nCnt * sizeof(WlzCMeshNod2D *))) == NULL) ||
     ((wgt = (double *)AlcMalloc(nCnt * sizeof(double))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    /* Collect nodes and possibly reduce the node count by checking for
     * nodes counted more than once in the initial estimate. */
    for(k = 0; k < 3; ++k)
    {
      nod[k] = elmNod[k];
    }
    for(i = 0; i < 3; ++i)
    {
      WlzCMeshEdgU2D *e0,
		     *e1;

      /* nNI[i] = k; */
      e0 = e1 = nod[i]->edu;
      do
      {
        int	skip = 0;
        WlzCMeshNod2D *n;

	n = e1->next->nod;
	for(j = 0; j < i; ++j)
	{
	  if(n == nod[j])
	  {
	    skip = 1;
	    break;
	  }
	}
	if(skip == 0)
	{
	  nod[k] = n;
	  ++k;
	}
        e1 = e1->nnxt;
      } while(e1 != e0);
    }
    nCnt = k;
    /* Compute approximate geodesic distance for each node and the total
     * approximate geodesic distance. This is the distance to a node of
     * the enclosing element plus a edge length for the neighbouring
     * nodes. */
    for(i = 0; i < nCnt; ++i)
    {
      WlzDVertex2 d;

      WLZ_VTX_2_SUB(d, nod[i]->pos, dPos);
      wgt[i] = WLZ_VTX_2_LENGTH(d);
      sWgt += wgt[i];
    }
    /*
    for(i = 0; i < 3; ++i)
    {
      WlzDVertex2 d;

      WLZ_VTX_2_SUB(d, nod[i]->pos, dPos);
      wgt[i] = WLZ_VTX_2_LENGTH(d);
    }
    for(i = 3; i < nNI[1]; ++i)
    {
      WlzDVertex2 d;

      WLZ_VTX_2_SUB(d, nod[i]->pos, nod[0]->pos);
      wgt[i] = wgt[0] + WLZ_VTX_2_LENGTH(d);
      sWgt += wgt[i];
    }
    for(i = nNI[1]; i < nNI[2]; ++i)
    {
      WlzDVertex2 d;

      WLZ_VTX_2_SUB(d, nod[i]->pos, nod[1]->pos);
      wgt[i] = wgt[1] + WLZ_VTX_2_LENGTH(d);
      sWgt += wgt[i];
    }
    for(i = nNI[2]; i < nCnt; ++i)
    {
      WlzDVertex2 d;

      WLZ_VTX_2_SUB(d, nod[i]->pos, nod[2]->pos);
      wgt[i] = wgt[2] + WLZ_VTX_2_LENGTH(d);
      sWgt += wgt[i];
    }
    sWgt += wgt[0] + wgt[1] + wgt[2];
    */
  }
  /* Compute distances weights and the distance weighted curvature. */
  if(sWgt < WLZ_MESH_TOLERANCE)
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    sWgt = 1.0 / sWgt;
    for(i = 0; i < nCnt; ++i)
    {
      double	*c;

      c = (double *)WlzIndexedValueGet(ixv, nod[i]->idx);
      crv += ((meanCrv != 0)? 0.5 * (c[0] + c[1]): (c[0] * c[1])) *
             (1.0 - (wgt[i] / sWgt));
    }
  }
  AlcFree(nod);
  AlcFree(wgt);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(crv);
}

/*!
* \return	New 2D domain object with double values corresponding to the
* 		interpolated mesh curvature or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new 2D domain object with double values that covers
* 		the given WLZ_CMESH_2D5 after it has been flattened by applying
* 		it's displacements. The displacements of the mesh must be
* 		valid for flattening the mesh as computed by
* 		WlzCMeshCompSurfMapConformal(). The curvature values are not
* 		normalised and are either the Gaussian or mean mesh curvatures.
* \param	inObj			Input object which must be a
* 					WLZ_CMESH_2D5 object.
* \param	scale			Scale factor to use from cmesh to
* 					2D spatial domain.
* \param 	meanCrv			If non zero the curvatures are the
* 					mean rather than the Gaussian
* 					curvatures.
* \param	interp			Type of interpolation to use, valid
* 					types are: WLZ_INTERPOLATION_NEAREST
* 					(value at nearest node of enclosing
* 					element), WLZ_INTERPOLATION_LINEAR
* 					(barycentric weighted value within
* 					the element) or
* 					WLZ_INTERPOLATION_ORDER_2 (values of
* 					in or on 2nd ring weighted by
* 					approximate geodesic distance).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzCMeshCurvToImage(WlzObject *inObj, double scale,
				     int meanCrv, WlzInterpolationType interp,
				     WlzErrorNum *dstErr)
{
  double	iScale = 1.0;
  WlzObject	*crvObj = NULL,
		*domObj = NULL,
		*fltObj = NULL,
  		*outObj = NULL;
  WlzValues 	val;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const double	eps = 0.000001;

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
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(fabs(scale) < eps)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(interp)
    {
      case WLZ_INTERPOLATION_NEAREST: /* FALLTHROUGH */
      case WLZ_INTERPOLATION_LINEAR:  /* FALLTHROUGH */
      case WLZ_INTERPOLATION_ORDER_2:
        break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    fltObj = WlzCMeshExtract2D(inObj, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    crvObj = WlzCMeshComputeCurvatures(inObj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    domObj = WlzCMeshToDomObj(fltObj, 0, scale, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzPixelV bgd;
    WlzObjectType gTabType;

    iScale = 1.0 / scale;
    bgd.type = WLZ_GREY_DOUBLE;
    bgd.v.dbv = 0.0;
    gTabType = WlzGreyTableType(WLZ_GREY_TAB_RAGR, WLZ_GREY_DOUBLE, NULL);
    val.v = WlzNewValueTb(domObj, gTabType, bgd, &errNum);
  }
  /* Make 2D the domain object for return with WLZ_GREY_DOUBLE values that
   * covers the flattened mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    outObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domObj->domain, val, NULL, NULL,
    		         &errNum);
  }
  /* Set either the gaussian or mean curvatures within the 2D object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMesh2D	*mesh;
    WlzIndexedValues *ixv = NULL;
    WlzGreyWSpace gWsp;
    WlzIntervalWSpace iWsp;

    ixv = crvObj->values.x;
    mesh = fltObj->domain.cm2;
    errNum = WlzInitGreyScan(outObj, &iWsp, &gWsp);
    if(errNum == WLZ_ERR_NONE)
    {
      while((errNum == WLZ_ERR_NONE) &&
	    ((errNum = WlzNextGreyInterval(&iWsp)) == WLZ_ERR_NONE))
      {
	int    idE,
	       idK,
	       idN;
	double d = 0;
	double *c,
	       *dst;
	WlzDVertex2 dPos;
	WlzDVertex3 del;

	idE = -1;
	dst = gWsp.u_grintptr.dbp;
	dPos.vtY = iWsp.linpos * iScale;
	idK = iWsp.lftpos;
	while((errNum == WLZ_ERR_NONE) && (idK <= iWsp.rgtpos))
	{
	  dPos.vtX = idK * iScale;
	  if((idE = WlzCMeshElmEnclosingPos2D(mesh, idE, dPos.vtX, dPos.vtY,
					      0, &idN)) >= 0)
	  {
	    int		idC;
	    WlzCMeshElm2D *elm;
	    double	crv[3],
			  dst[3];
	    WlzCMeshNod2D *nod[3];

	    elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	    nod[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm);
	    nod[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm);
	    nod[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm);
	    for(idC = 0; idC < 3; ++idC)
	    {
	      c = (double *)WlzIndexedValueGet(ixv, nod[idC]->idx);
	      crv[idC] = (meanCrv != 0)? 0.5 * (c[0] + c[1]): c[0] * c[1];
	    }
	    switch(interp)
	    {
	      case WLZ_INTERPOLATION_NEAREST:
		WLZ_VTX_2_SUB(del, nod[0]->pos, dPos);
		dst[0] = WLZ_VTX_2_SQRLEN(del);
		WLZ_VTX_2_SUB(del, nod[1]->pos, dPos);
		dst[1] = WLZ_VTX_2_SQRLEN(del);
		WLZ_VTX_2_SUB(del, nod[2]->pos, dPos);
		dst[2] = WLZ_VTX_2_SQRLEN(del);
		if(dst[0] < dst[1])
		{
		  d = (dst[0] < dst[2])? crv[0]: crv[2];
		}
		else
		{
		  d = (dst[1] < dst[2])? crv[1]: crv[2];
		}
		break;
	      case WLZ_INTERPOLATION_LINEAR:
		d = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos,
					    nod[2]->pos, crv[0], crv[1], crv[2],
					    dPos);
		break;
	      case WLZ_INTERPOLATION_ORDER_2:
		d = WlzCMeshElmRingInterpolateNodCrvValue(elm, ixv, nod,
					  meanCrv, dPos, &errNum);
		break;
	      default:
		break;
	    }
	  }
	  else
	  {
	    c = (double *)WlzIndexedValueGet(ixv, idN);
	    if(c == NULL)
	    {
	      d = 0.0;
	    }
	    else
	    {
	      d = (meanCrv != 0)? 0.5 * (c[0] + c[1]): c[0] * c[1];
	    }
	  }
	  *dst++ = d;
	  ++idK;
	}
      }
      (void )WlzEndGreyScan(&gWsp);
      if(errNum == WLZ_ERR_EOO)
      {
	errNum = WLZ_ERR_NONE;
      }
    }
  }
  (void )WlzFreeObj(crvObj);
  (void )WlzFreeObj(domObj);
  (void )WlzFreeObj(fltObj);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(outObj);
    outObj = NULL;
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(outObj);
}

/*!
* \return	New object with the domain of the given object and curvatures
* 		in values or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the domain of the given object and
* 		a new indexed value table which has the curvatures. The
* 		curvatures are organized as a pair of double values per
* 		node, total curvature first and then the mean curvature.
* \param	inObj			Input object which must have a valid
* 					WLZ_CMESH_2D5 domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject     	*WlzCMeshComputeCurvatures(WlzObject *inObj,
					   WlzErrorNum *dstErr)
{
  WlzObject	*eNObj = NULL,
		*nNObj = NULL,
  		*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Let WlzCMeshComputeNormalsElm() do the checking of the given object. */
  eNObj = WlzCMeshComputeNormalsElm(inObj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    nNObj = WlzCMeshComputeNormalsNod(eNObj, &errNum);
  }
  (void )WlzFreeObj(eNObj);
  if(errNum == WLZ_ERR_NONE)
  {
    outObj = WlzCMeshComputeCurvaturesFromNodNorm(nNObj, &errNum);
  }
  (void )WlzFreeObj(nNObj);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(outObj);
}

/*!
* \return	New object with the domain of the given object and element
* 		normals in values or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the domain of the given object and
* 		a new indexed value table which has the element normals. The
* 		normals are organized as a 3 double values (x, y, z order)
* 		per element.
* \param	inObj			Input object which must have a valid
* 					WLZ_CMESH_2D5 domain.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject     	*WlzCMeshComputeNormalsElm(WlzObject *inObj,
					   WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*mesh;
  WlzIndexedValues *ixv = NULL;
  WlzObject	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if (inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh = inObj->domain.cm2d5)->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    int	dim = 3;

    ixv = WlzMakeIndexedValues(inObj, 1, &dim, WLZ_GREY_DOUBLE,
                               WLZ_VALUE_ATTACH_ELM, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm2D5 *elm;

      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	double	*n;
	WlzDVertex3 nrm;
        WlzCMeshNod2D5 *nod[3];

        WlzCMeshElmGetNodes2D5(elm, nod + 0, nod + 1, nod + 2);
	nrm = WlzGeomTriangleNormal(nod[0]->pos, nod[1]->pos, nod[2]->pos);
	n = (double *)WlzIndexedValueGet(ixv, elm->idx);
	n[0] = nrm.vtX; n[1] = nrm.vtY; n[2] = nrm.vtZ;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	val;

    val.x = ixv;
    outObj = WlzMakeMain(WLZ_CMESH_2D5, inObj->domain, val, NULL, NULL,
    			 &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(outObj != NULL)
    {
      (void )WlzFreeObj(outObj);
      outObj = NULL;
    }
    else if(ixv != NULL)
    {
      (void )WlzFreeIndexedValues(ixv);
    }
  }
  return(outObj);
}

/*!
* \return	New object with the domain of the given object and node
* 		normals in values or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the domain of the given object and
* 		a new indexed value table which has the node normals. The
* 		normals are organized as a 3 double values (x, y, z order)
* 		per node.
* \param	inObj			Input object which must have a valid
* 					WLZ_CMESH_2D5 domain and have an
* 					indexed value table with a single
* 					normal (3 doubles) attached to each
* 					element.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject     	*WlzCMeshComputeNormalsNod(WlzObject *inObj,
					   WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*mesh;
  WlzIndexedValues *eIxv,
  		*nIxv = NULL;
  WlzObject	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if (inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh = inObj->domain.cm2d5)->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    eIxv = inObj->values.x;

    if((eIxv->attach != WLZ_VALUE_ATTACH_ELM) ||
       (eIxv->rank < 1) || (eIxv->dim[0] < 3) ||
       (eIxv->vType != WLZ_GREY_DOUBLE))
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int	dim = 3;

    nIxv = WlzMakeIndexedValues(inObj, 1, &dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;
    WlzCMeshNod2D5 *nod;

    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	double l;
	double *n;
	WlzDVertex3 nrm;
	WlzCMeshEdgU2D5 *edu0,
			*edu1;

	WLZ_VTX_3_ZERO(nrm);
        edu1 = edu0 = nod->edu;
	do
	{
	  n = (double *)WlzIndexedValueGet(eIxv, edu1->elm->idx);
	  nrm.vtX += n[0]; nrm.vtY += n[1]; nrm.vtZ += n[2];
	  edu1 = edu1->nnxt;
	} while(edu1 != edu0);
	l = WLZ_VTX_3_LENGTH(nrm);
	if(l > DBL_EPSILON)
	{
	  l = 1.0 / l;
	  WLZ_VTX_3_SCALE(nrm, nrm, l);
	}
	n = (double *)WlzIndexedValueGet(nIxv, nod->idx);
	n[0] = nrm.vtX; n[1] = nrm.vtY; n[2] = nrm.vtZ;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	val;

    val.x = nIxv;
    outObj = WlzMakeMain(WLZ_CMESH_2D5, inObj->domain, val, NULL, NULL,
    			 &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(outObj != NULL)
    {
      (void )WlzFreeObj(outObj);
      outObj = NULL;
    }
    else if(nIxv != NULL)
    {
      (void )WlzFreeIndexedValues(nIxv);
    }
  }
  return(outObj);
}

/*!
* \return	New object with the domain of the given object and curvatures
* 		in values or NULL on error.
* \ingroup	WlzMesh
* \brief	Creates a new object with the domain of the given object and
* 		a new indexed value table which has the curvatures. The
* 		curvatures are organized as a pair of double values per
* 		node, total curvature first and then the mean curvature.
* \param	inObj			Input object which must have a valid
* 					WLZ_CMESH_2D5 domain and have an
* 					indexed value table with a single
* 					normal (3 doubles) attached to each
* 					node.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject     	*WlzCMeshComputeCurvaturesFromNodNorm(WlzObject *inObj,
				        WlzErrorNum *dstErr)
{
  WlzCMesh2D5	*mesh;
  WlzIndexedValues *nIxv = NULL,
  		   *cIxv = NULL;
  WlzObject	*outObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(inObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(inObj->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if (inObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((mesh = inObj->domain.cm2d5)->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(inObj->values.core == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    nIxv = inObj->values.x;

    if((nIxv->attach != WLZ_VALUE_ATTACH_NOD) ||
       (nIxv->rank < 1) || (nIxv->dim[0] < 3) ||
       (nIxv->vType != WLZ_GREY_DOUBLE))
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int	dim = 2;

    cIxv = WlzMakeIndexedValues(inObj, 1, &dim, WLZ_GREY_DOUBLE,
                                WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
    		idxBufMax = 0,
		posBufMax = 0;
    int		*idxBuf = NULL;
    WlzDVertex3	*posBuf = NULL;
    AlcVector	*nodVec;
    WlzCMeshNod2D5 *nod;

    nodVec = mesh->res.nod.vec;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(nodVec, idN);
      if(nod->idx >= 0)
      {
	double *crv;

	crv = (double *)WlzIndexedValueGet(cIxv, idN);
	if(WlzCMeshNodIsBoundary2D5(nod))
	{
	  crv[0] = crv[1] = 0.0;
	}
	else
	{
	  int	nN;
	  double *n;
	  WlzDVertex3 nrm;

	  /* Get edge connected neighbours of the node. */
	  n = (double *)WlzIndexedValueGet(nIxv, nod->idx);
	  nrm.vtX = n[0]; nrm.vtY = n[1]; nrm.vtZ = n[2];
	  nN = WlzCMeshNodRingNodIndices2D5(nod, &idxBufMax, &idxBuf, &errNum);
	  if(posBufMax < idxBufMax)
	  {
	    posBufMax = idxBufMax;
	    if((posBuf = (WlzDVertex3 *)
			 AlcRealloc(posBuf,
				    sizeof(WlzDVertex3) * posBufMax)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int	 idN;

	    /* Get the positions of the node and it's edge connected
	     * neighbours. */
	    for(idN = 0; idN < nN; ++idN)
	    {
	      WlzCMeshNod2D5 *nod1;

	      nod1 = (WlzCMeshNod2D5 *)AlcVectorItemGet(nodVec, idxBuf[idN]);
	      posBuf[idN] = nod1->pos;
	    }
	    /* Compute the curvature values for the node. */
	    errNum = WlzGeomCurvature(2, crv, nrm, nN, posBuf);
	  }
	}
      }
    }
    AlcFree(idxBuf);
    AlcFree(posBuf);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	val;

    val.x = cIxv;
    outObj = WlzMakeMain(WLZ_CMESH_2D5, inObj->domain, val, NULL, NULL,
    			 &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(outObj != NULL)
    {
      (void )WlzFreeObj(outObj);
      outObj = NULL;
    }
    else if(nIxv != NULL)
    {
      (void )WlzFreeIndexedValues(nIxv);
    }
  }
  return(outObj);
}
