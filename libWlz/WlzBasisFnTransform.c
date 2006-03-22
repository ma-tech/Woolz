#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzBasisFnTransform.c
* \author       Bill Hill, Jianguo Rao
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Functions for computing and applying basis function
* 		transforms.
* \ingroup	WlzTransform
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>


/*!
* \return	New Basis function transform.
* \ingroup	WlzTransform
* \brief	Makes a new basis function transform.
*		The transform will be returned with a NULL basis function
*		pointer which needs to be set along with all the other
*		fields of the transform structure.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzMakeBasisFnTransform(WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if((basisTr = (WlzBasisFnTransform *)
                AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisTr);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzTransform
* \brief	Free's the given basis function transform.
* \param	basisTr			Given basis function transform,
					may be NULL.
*/
WlzErrorNum	WlzBasisFnFreeTransform(WlzBasisFnTransform *basisTr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basisTr)
  {
    errNum = WlzBasisFnFree(basisTr->basisFn);
    AlcFree(basisTr);
  }
  return(errNum);
}

/*!
* \return	New basis function transform.
* \ingroup	WlzTransform
* \brief	Creates a new basis function transform of the given
*		type, which will transform an object with the given
*		source verticies into an object with the given destination
*		verticies.
*		If a constraining object is given all distances will be
*		computed within the given object for those basis functions
*		which support constrained evaluation (Gauss, multi-quadric
*		and thin-plate spline).
* \param	type			Required basis function type.
* \param	order			Order of polynomial, only used for
* 					WLZ_FN_BASIS_2DPOLY.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts).
* \param	sPts			Source control points.
* \param	cObj			Constraining object, within which all
*					distances are constrained. If NULL
*					Euclidean distances are used in place
*					of constrained distances.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts2D(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex2 *dPts,
					  int nSPts,
					  WlzDVertex2 *sPts,
					  WlzObject *cObj,
					  WlzErrorNum *dstErr)
{
  int		idx;
  WlzBasisFnTransform *basisTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	deltaMQ = 0.001;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(cObj)
  {
    if(cObj->type != WLZ_2D_DOMAINOBJ)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(cObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((basisTr = WlzMakeBasisFnTransform(NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    basisTr->type = WLZ_TRANSFORM_2D_BASISFN;
    switch(type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
	basisTr->basisFn = WlzBasisFnGauss2DFromCPts(nDPts,
					dPts, sPts, 0.9, cObj,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DPOLY:
	basisTr->basisFn = WlzBasisFnPoly2DFromCPts(nDPts, order,
					dPts, sPts,
				         &errNum);
	break;
      case WLZ_FN_BASIS_2DMQ:
	basisTr->basisFn = WlzBasisFnMQ2DFromCPts(nDPts,
					dPts, sPts, deltaMQ, cObj,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DTPS:
	basisTr->basisFn = WlzBasisFnTPS2DFromCPts(nDPts,
					dPts, sPts, cObj,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DCONF_POLY:
	basisTr->basisFn = WlzBasisFnConf2DFromCPts(nDPts, order,
					dPts, sPts,
				        &errNum);
	break;
      default:
	 errNum = WLZ_ERR_TRANSFORM_TYPE;
	 break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzBasisFnFreeTransform(basisTr);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisTr);
}


/*!
* \return	New basis function transform.
* \ingroup	WlzTransform
* \brief	Creates a new basis function transform of the given
*		type, which will transform an object with the given
*		source verticies into an object with the given 	
*		destination verticies.			
* \param	type	        	Required basis function type.
* \param	order			Order of polynomial, only 
*					used for WLZ_BASISFN_POLY.
* \param	nDPts			Number of destination control
*					points.			
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts)
* \param	sPts			Source control points.	
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts3(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex3 *dPts,
					  int nSPts,
					  WlzDVertex3 *sPts,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if((basisTr = WlzMakeBasisFnTransform(NULL)) == NULL)
  {
   errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    switch(type)
    { 
      case WLZ_FN_BASIS_3DMQ:
        {
	  basisTr->basisFn = WlzBasisFnMQ3DFromCPts(nDPts, dPts, sPts, 0.2,
	  					    &errNum);
	  basisTr->linkcount = 0;
	  basisTr->freeptr   = NULL;
	  
	}
	break;
      default:
	 errNum = WLZ_ERR_TRANSFORM_TYPE;
	 break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisTr);
}

/*!
* \return	Error number.
* \ingroup	WlzTransform
* \brief	Sets the displacements of the given mesh transform according
* 		to the basis function transform.
* \param	mesh			Given mesh transform.
* \param	basisTr			Given basis functiontransform.
*/
WlzErrorNum    	WlzBasisFnSetMesh(WlzMeshTransform *mesh,
				  WlzBasisFnTransform *basisTr)
{
  int		nodCnt;
  WlzMeshNode	*nod;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mesh == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((mesh->type != WLZ_TRANSFORM_2D_MESH) ||
	  (basisTr->type != WLZ_TRANSFORM_2D_BASISFN))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    nod = mesh->nodes;
    nodCnt = mesh->nNodes;
    switch(basisTr->basisFn->type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValueGauss2D(basisTr->basisFn,
	  				nod->position);
	  ++nod;
	}
	break;
      case WLZ_FN_BASIS_2DPOLY:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValuePoly2D(basisTr->basisFn,
	  				nod->position);
	  ++nod;
	}
	break;
      case WLZ_FN_BASIS_2DMQ:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValueMQ2D(basisTr->basisFn,
	  				nod->position);
	  ++nod;
	}
	break;
      case WLZ_FN_BASIS_2DTPS:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValueTPS2D(basisTr->basisFn,
	  				nod->position);
	  ++nod;
	}
	break;
    case WLZ_FN_BASIS_2DCONF_POLY:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValueConf2D(basisTr->basisFn,
	  				nod->position);
	  ++nod;
	}
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzTransform
* \brief	Sets the displacements of the given conforming mesh
*		transform according to the basis function transform.
* \param	meshTr			Given conforming mesh transform.
* \param	basisTr			Given basis functiontransform.
*/
WlzErrorNum    	WlzBasisFnSetCMesh(WlzCMeshTransform *meshTr,
				   WlzBasisFnTransform *basisTr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((meshTr == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(meshTr->type)
    {
      case WLZ_TRANSFORM_2D_CMESH:
        errNum = WlzBasisFnSetCMesh2D(meshTr, basisTr);
	break;
      case WLZ_TRANSFORM_3D_CMESH:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzTransform
* \brief	Sets the displacements of the given 2D conforming mesh
*		transform according to the basis function transform.
* \param	meshTr			Given mesh transform.
* \param	basisTr			Given basis functiontransform.
*/
WlzErrorNum    	WlzBasisFnSetCMesh2D(WlzCMeshTransform *meshTr,
				     WlzBasisFnTransform *basisTr)
{
  int		idN;
  WlzDVertex2	*dspP;
  WlzCMeshNod2D	*nod;
  WlzCMesh2D	*mesh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((meshTr == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((meshTr->type != WLZ_TRANSFORM_2D_CMESH) ||
	  (basisTr->type != WLZ_TRANSFORM_2D_BASISFN))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((mesh = meshTr->mesh.m2) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
        dspP = (WlzDVertex2 *)AlcVectorItemGet(meshTr->dspVec, idN);
	switch(basisTr->basisFn->type)
	{
	  case WLZ_FN_BASIS_2DGAUSS:
	    *dspP = WlzBasisFnValueGauss2D(basisTr->basisFn, nod->pos);
	    break;
	  case WLZ_FN_BASIS_2DMQ:
	    *dspP = WlzBasisFnValueMQ2D(basisTr->basisFn, nod->pos);
	    break;
	  case WLZ_FN_BASIS_2DTPS:
	    *dspP = WlzBasisFnValueTPS2D(basisTr->basisFn, nod->pos);
	    break;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Transformed object, NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms a woolz object using a the given basis function
*		transform. This function has been written as an example of how
*		to transform an object using a basis function and mesh. In most
*		cases WlzMeshFromObj(), WlzBasisFnSetMesh() and
*		WlzMeshTransformObj() would be called allowing a mesh to be
*		reused.
* \param	srcObj			Object to be transformed.
* \param	basisTr			Basis function transform to apply.
* \param	interp			Level of interpolation to use.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzBasisFnTransformObj(WlzObject *srcObj,
					WlzBasisFnTransform *basisTr,
					WlzInterpolationType interp,
					WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzDomain	dstDom;
  WlzValues	dumVal;
  WlzMeshTransform *mesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dumVal.core = NULL;
  if((srcObj == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(srcObj->type)
    {
      case WLZ_EMPTY_OBJ:
        if((dstObj = WlzMakeEmpty(&errNum)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	break;
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
	/* TODO: need better mesh generation */
	mesh = WlzMeshFromObj(srcObj, WLZ_MESH_GENMETHOD_BLOCK, 100.0, 100.0,
	                      &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzBasisFnSetMesh(mesh, basisTr);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dstObj = WlzMeshTransformObj(srcObj, mesh, interp, &errNum);
	}
	if(mesh)
	{
	  WlzMeshFreeTransform(mesh);
	}
        break;
      case WLZ_2D_POLYGON: /* FALLTHROUGH */
      case WLZ_BOUNDLIST: /* FALLTHROUGH */
      case WLZ_CONTOUR:
	if(srcObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  switch(srcObj->type)
	  {
	    case WLZ_2D_POLYGON:
	      dstDom.poly = WlzBasisFnTransformPoly2(srcObj->domain.poly,
	      					     basisTr, 1, &errNum);
	      break;
	    case WLZ_BOUNDLIST:
	      dstDom.b = WlzBasisFnTransformBoundList(srcObj->domain.b,
	      					      basisTr, 1, &errNum);
	      break;
	    case WLZ_CONTOUR:
	      dstDom.ctr = WlzBasisFnTransformContour(srcObj->domain.ctr,
	      					      basisTr, 1, &errNum);
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  if((dstObj = WlzMakeMain(srcObj->type, dstDom, dumVal,
	  			   NULL, NULL, NULL)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    (void )WlzFreeDomain(dstDom);
	  }
	}
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
  return(dstObj);
}

/*!
* \return	Transformed 2D polygon domain, NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms a 2D polygon domain using a the given basis function
*		transform.
* \param	srcPoly			Polygon domain to be transformed.
* \param	basisTr			Basis function transform to apply.
* \param	newPoly			Makes a new polygon domain if non-zero
*					otherwise the given polygon domain
*					will be transformed in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPolygonDomain *WlzBasisFnTransformPoly2(WlzPolygonDomain *srcPoly,
					   WlzBasisFnTransform *basisTr,
					   int newPoly,
					   WlzErrorNum *dstErr)
{
  int		idN;
  WlzVertexP	sVP,
  		dVP;
  WlzDVertex2	tD;
  WlzPolygonDomain *dstPoly = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcPoly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((srcPoly->type != WLZ_POLYGON_INT) &&
          (srcPoly->type != WLZ_POLYGON_FLOAT) &&
          (srcPoly->type != WLZ_POLYGON_DOUBLE))
  {
    errNum = WLZ_ERR_POLYGON_TYPE;
  }
  else if(newPoly &&
          ((dstPoly = WlzMakePolygonDomain(srcPoly->type, 0, NULL,
  					   srcPoly->nvertices, 1,
					   &errNum)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    sVP.v = srcPoly->vtx;
    if(newPoly)
    {
      dVP.v = dstPoly->vtx;
      dstPoly->nvertices = srcPoly->nvertices;
    }
    else
    {
      dVP.v = srcPoly->vtx;
    }
    switch(srcPoly->type)
    {
      case WLZ_POLYGON_INT:
	tD.vtX = sVP.i2->vtX;
	tD.vtY = sVP.i2->vtY;
	tD = WlzBasisFnTransformVertexD(basisTr, tD, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  dVP.i2->vtX = tD.vtX;
	  dVP.i2->vtY = tD.vtY;
	  for(idN = 1; idN < srcPoly->nvertices; ++idN)
	  {
	    ++(sVP.i2);
	    ++(dVP.i2);
	    tD.vtX = sVP.i2->vtX;
	    tD.vtY = sVP.i2->vtY;
	    tD = WlzBasisFnTransformVertexD(basisTr, tD, NULL);
	    dVP.i2->vtX = tD.vtX;
	    dVP.i2->vtY = tD.vtY;
	  }
	}
        break;
      case WLZ_POLYGON_FLOAT:
	tD.vtX = sVP.f2->vtX;
	tD.vtY = sVP.f2->vtY;
	tD = WlzBasisFnTransformVertexD(basisTr, tD, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  dVP.f2->vtX = tD.vtX;
	  dVP.f2->vtY = tD.vtY;
	  for(idN = 1; idN < srcPoly->nvertices; ++idN)
	  {
	    ++(sVP.f2);
	    ++(dVP.f2);
	    tD.vtX = sVP.f2->vtX;
	    tD.vtY = sVP.f2->vtY;
	    tD = WlzBasisFnTransformVertexD(basisTr, tD, NULL);
	    dVP.f2->vtX = tD.vtX;
	    dVP.f2->vtY = tD.vtY;
	  }
	}
        break;
      case WLZ_POLYGON_DOUBLE:
	*(dVP.d2) = WlzBasisFnTransformVertexD(basisTr, *(sVP.d2), &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  for(idN = 1; idN < srcPoly->nvertices; ++idN)
	  {
	    ++(sVP.d2);
	    ++(dVP.d2);
	    *(dVP.d2) = WlzBasisFnTransformVertexD(basisTr, *(sVP.d2), NULL);
	  }
	}
        break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreePolyDmn(dstPoly);
      dstPoly = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstPoly);
}

/*!
* \return	Transformed 2D boundary list, NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms a 2D boundary list using a the given basis function
*		transform.
* \param	srcBnd			Boundary list to be transformed.
* \param	basisTr			Basis function transform to apply.
* \param	newBnd			Makes a new bound list if non-zero
*					otherwise the given bound list
*					will be transformed in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBoundList *WlzBasisFnTransformBoundList(WlzBoundList *srcBnd,
					   WlzBasisFnTransform *basisTr,
					   int newBnd,
					   WlzErrorNum *dstErr)
{
  WlzDomain	tDom;
  WlzBoundList	*dstBnd = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcBnd == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(srcBnd->type != WLZ_BOUNDLIST_PIECE)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    if(newBnd)
    {
      dstBnd->type = srcBnd->type;
      dstBnd->wrap = srcBnd->wrap;
      if((dstBnd = (WlzBoundList *)AlcCalloc(sizeof(WlzBoundList), 1)) == NULL)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    else
    {
      dstBnd = srcBnd;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Transform the polygon */
    if((dstBnd->poly = WlzBasisFnTransformPoly2(srcBnd->poly, basisTr,
                                                newBnd, &errNum)) != NULL)
    {
      /* Transform next */
      if(srcBnd->next)
      {
        if((tDom.b = WlzBasisFnTransformBoundList(srcBnd->next, basisTr,
					          newBnd, &errNum)) != NULL)
        {
	  if(newBnd)
	  {
            (void )WlzAssignDomain(tDom, NULL);
	  }
          dstBnd->next = tDom.b;
        }
      }
      /* Transform down */
      if(srcBnd->down && (errNum == WLZ_ERR_NONE))
      {
        if((tDom.b = WlzBasisFnTransformBoundList(srcBnd->down, basisTr,
                                                  newBnd, &errNum)) != NULL)
        {
	  if(newBnd)
	  {
            (void )WlzAssignDomain(tDom, NULL);
	  }
          dstBnd->down = tDom.b;
        }
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      if(newBnd)
      {
        (void )WlzFreeBoundList(dstBnd);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstBnd);
}

/*!
* \return	Transformed contour, NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms a contour using a the given basis function
*		transform. See WlzBasisFnTransformGMModel().
* \param	srcCtr			Contour to be transformed.
* \param	basisTr			Basis function transform to apply.
* \param	newCtr			Makes a new contour if non-zero
*					otherwise the given contour will be
*					transformed in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzContour	*WlzBasisFnTransformContour(WlzContour *srcCtr,
					    WlzBasisFnTransform *basisTr,
					    int newCtr,
					    WlzErrorNum *dstErr)
{
  WlzContour	*dstCtr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstCtr = (newCtr)? WlzMakeContour(&errNum): srcCtr;
  if(errNum == WLZ_ERR_NONE)
  {
    dstCtr->model = WlzBasisFnTransformGMModel(srcCtr->model, basisTr,
    					       newCtr, &errNum);
    if(newCtr)
    {
      (void )WlzAssignGMModel(dstCtr->model, NULL);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(newCtr && dstCtr)
    {
      (void )WlzFreeContour(dstCtr);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstCtr);
}

/*!
* \return	Transformed model, NULL on error.
* \ingroup	WlzTransform
* \brief	Transforms a Woolz GMModel using a the given basis function
*		transform. This function assumes that the transformation
*		does not change the topology of the model, watch out
*		that this willnot always be true!
* \param	srcM			Model to be transformed.
* \param	basisTr			Basis function transform to apply.
* \param	newModel		Makes a new model if non-zero
*					otherwise the given model will be
*					transformed in place.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzGMModel	*WlzBasisFnTransformGMModel(WlzGMModel *srcM,
					    WlzBasisFnTransform *basisTr,
					    int newModel, WlzErrorNum *dstErr)
{
  int		idx,
  		cnt;
  WlzDVertex2	tD;
  AlcVector	*vec;
  WlzGMElemP	elmP;
  WlzGMModel	*dstM = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dstM = (newModel)? WlzGMModelCopy(srcM, &errNum): srcM;
  /* Transform vertex geometries. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    vec = dstM->res.vertexG.vec;
    cnt = dstM->res.vertexG.numIdx;
    while((idx < cnt) && (errNum == WLZ_ERR_NONE))
    {
      elmP.core = (WlzGMCore *)AlcVectorItemGet(vec, idx);
      if(elmP.core && (elmP.core->idx >= 0))
      {
        switch(dstM->type)
        {
          case WLZ_GMMOD_2I:
	    tD.vtX = elmP.vertexG2I->vtx.vtX;
	    tD.vtY = elmP.vertexG2I->vtx.vtY;
	    tD = WlzBasisFnTransformVertexD(basisTr, tD, &errNum);
	    elmP.vertexG2I->vtx.vtX = WLZ_NINT(tD.vtX);
	    elmP.vertexG2I->vtx.vtY = WLZ_NINT(tD.vtY);
            break;
          case WLZ_GMMOD_2D:
	    elmP.vertexG2D->vtx = WlzBasisFnTransformVertexD(basisTr,
	    				elmP.vertexG2D->vtx, &errNum);
            break;
          case WLZ_GMMOD_2N:
	    elmP.vertexG2N->nrm = WlzBasisFnTransformNormalD(basisTr,
					elmP.vertexG2N->vtx,
					elmP.vertexG2N->nrm,
					&(elmP.vertexG2N->vtx),
					&errNum);
            break;
          default:
            errNum = WLZ_ERR_DOMAIN_TYPE;
            break;
        }
      }
      ++idx;
    }
  }
  /* Compute shell geometries. */
  if(errNum == WLZ_ERR_NONE)
  {
    idx = 0;
    vec = dstM->res.shell.vec;
    cnt = dstM->res.shell.numIdx;
    while((idx < cnt) && (errNum == WLZ_ERR_NONE))
    {
      elmP.core = (WlzGMCore *)AlcVectorItemGet(vec, idx);
      if(elmP.core && (elmP.core->idx >= 0))
      {
        errNum = WlzGMShellComputeGBB(elmP.shell);
      }
      ++idx;
    }
  }
  /* Rehash the vertex location table. */
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzGMModelRehashVHT(dstM, 0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstM);
}

/*!
* \return	Transformed vertex.
* \ingroup	WlzTransform
* \brief	Transforms the given vertex and it's normal which is
*		unit length and directed from the given vertex.
* \param	basisTr			Basis function transform to apply.
* \param	srcVx			Given vertex.
* \param	srcNr			Given normal.
* \param	dstVx			Destination pointer for the
*					transformed vertex. Must not
*					be NULL.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
WlzDVertex2	WlzBasisFnTransformNormalD(WlzBasisFnTransform *basisTr,
					   WlzDVertex2 srcVx,
					   WlzDVertex2 srcNr,
					   WlzDVertex2 *dstVx,
					   WlzErrorNum *dstErr)
{
  double	len;
  WlzDVertex2	tVx,
  		dVx,
  		dNr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basisTr == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(basisTr->type != WLZ_TRANSFORM_2D_BASISFN)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    dVx = WlzBasisFnTransformVertexD(basisTr, srcVx, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WLZ_VTX_2_ADD(tVx, srcVx, srcNr);
    tVx = WlzBasisFnTransformVertexD(basisTr, tVx, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WLZ_VTX_2_SUB(tVx, tVx, dVx);
    len = WLZ_VTX_2_SQRLEN(tVx);
    if(len > DBL_EPSILON)
    {
      len = 1.0 / sqrt(len);
      WLZ_VTX_2_SCALE(tVx, tVx, len);
    }
    else
    {
      WLZ_VTX_2_ZERO(tVx);
      errNum = WLZ_ERR_TRANSFORM_DATA;
    }
    dNr = tVx;
    *dstVx = dVx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dNr);
}

/*!
* \return	Transformed vertex.
* \ingroup	WlzTransform
* \brief	Transforms the given WlzDVertex2.
* \param	basisTr			Basis function transform to apply.
* \param	srcVx			Vertex to be transformed.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
WlzDVertex2	WlzBasisFnTransformVertexD(WlzBasisFnTransform *basisTr,
					   WlzDVertex2 srcVx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dstVx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basisTr == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(basisTr->type != WLZ_TRANSFORM_2D_BASISFN)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    switch(basisTr->basisFn->type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
	dstVx = WlzBasisFnValueGauss2D(basisTr->basisFn, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_FN_BASIS_2DPOLY:
	dstVx = WlzBasisFnValuePoly2D(basisTr->basisFn, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_FN_BASIS_2DMQ:
	dstVx = WlzBasisFnValueMQ2D(basisTr->basisFn, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_FN_BASIS_2DTPS:
	dstVx = WlzBasisFnValueTPS2D(basisTr->basisFn, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_FN_BASIS_2DCONF_POLY:
	dstVx = WlzBasisFnValueConf2D(basisTr->basisFn, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstVx);
}

/*!
* \return	Transformed vertex.
* \ingroup	WlzTransform
* \brief	Transforms the given WlzFVertex2.
* \param	basisTr			Basis function transform to apply.
* \param	srcVxF			Vertex to be transformed.
* \param	dstErr			Destination pointer for error, may be
* 					NULL.
*/
WlzFVertex2	WlzBasisFnTransformVertexF(WlzBasisFnTransform *basisTr,
					   WlzFVertex2 srcVxF,
					   WlzErrorNum *dstErr)
{
  WlzFVertex2	dstVxF;
  WlzDVertex2	dstVxD,
		srcVxD;

  srcVxD.vtX = srcVxF.vtX;
  srcVxD.vtY = srcVxF.vtY;
  dstVxD = WlzBasisFnTransformVertexD(basisTr, srcVxD, dstErr);
  dstVxF.vtX = dstVxD.vtX;
  dstVxF.vtY = dstVxD.vtY;
  return(dstVxF);
}

/*!
* \return	Transformed vertex.
* \ingroup	WlzTransform
* \brief	Transforms the given WlzIVertex2.
* \param	basisTr			Basis function transform to apply.
* \param	srcVxI			Vertex to be transformed.
* \param	dstErr			Destination pointer for error, may be
*					NULL.
*/
WlzIVertex2	WlzBasisFnTransformVertexI(WlzBasisFnTransform *basisTr,
					   WlzIVertex2 srcVxI,
					   WlzErrorNum *dstErr)
{
  WlzIVertex2	dstVxI;
  WlzDVertex2	dstVxD,
		srcVxD;

  srcVxD.vtX = srcVxI.vtX;
  srcVxD.vtY = srcVxI.vtY;
  dstVxD = WlzBasisFnTransformVertexD(basisTr, srcVxD, dstErr);
  dstVxI.vtX = WLZ_NINT(dstVxD.vtX);
  dstVxI.vtY = WLZ_NINT(dstVxD.vtY);
  return(dstVxI);
}
