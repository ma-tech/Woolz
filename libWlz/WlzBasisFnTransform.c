#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBasisFnTransform.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz functions for computing and applying basis function
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
*		type, which will transform an object with the given source
*		verticies into an object with the given destination verticies.
* \param	type			Required basis function type.
* \param	order			Order of polynomial, only used for
* 					WLZ_FN_BASIS_2DPOLY.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts).
* \param	sPts			Source control points.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts2D(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex2 *dPts,
					  int nSPts,
					  WlzDVertex2 *sPts,
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
    basisTr->type = WLZ_TRANSFORM_2D_BASISFN;
    switch(type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
	basisTr->basisFn = WlzBasisFnGauss2DFromCPts(nDPts,
					dPts, sPts, 0.9,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DPOLY:
	basisTr->basisFn = WlzBasisFnPoly2DFromCPts(nDPts, order,
					dPts, sPts,
				         &errNum);
	break;
      case WLZ_FN_BASIS_2DMQ:
	basisTr->basisFn = WlzBasisFnMQ2DFromCPts(nDPts,
					dPts, sPts, 0.1,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DTPS:
	basisTr->basisFn = WlzBasisFnTPS2DFromCPts(nDPts,
					dPts, sPts,
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
  WlzMeshTransform *mesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcObj == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    /* TODO: need better mesh generation */
    mesh = WlzMeshFromObj(srcObj, WLZ_MESH_GENMETHOD_BLOCK, 100.0, 100.0,
			  &errNum);
  }
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
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dstObj);
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
