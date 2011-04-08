#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBasisFnTransform_c[] = "MRC HGU $Id$";
#endif
#endif
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
*		verticies. See WlzBasisFnTrFromCPts2DParam().
* \param	type			Required basis function type.
* \param	order			Order of polynomial, only used for
* 					WLZ_FN_BASIS_2DPOLY.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts).
* \param	sPts			Source control points.
* \param        mesh                    Mesh which is used to compute
*                                       constrained distances. If non NULL
*                                       and the mesh type is
*                                       WLZ_CMESH_2D then  constrained
*                                       distances are used and these are
*                                       computed using the mesh.
*                                       If NULL or the transform is
*                                       some other type then Euclidean
*                                       distances are used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts2D(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex2 *dPts,
					  int nSPts,
					  WlzDVertex2 *sPts,
					  WlzCMesh2D *mesh,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr = NULL;

  basisTr = WlzBasisFnTrFromCPts2DParam(type, order, nDPts, dPts,
                                        nSPts, sPts, mesh, 0, NULL,
					dstErr);
  return(basisTr);
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
*		Additional basis functions parameters may be supplied via
*		the nParam and param parameters. Currently this is only used to
*		supply the multi-quadric delta or gauss parameter scaling.
*		The default values of multi-quadric delta = 0.001 and
*		gauss param = 0.9 are used if nParam <= 0 or param == NULL.
* \param	type			Required basis function type.
* \param	order			Order of polynomial, only used for
* 					WLZ_FN_BASIS_2DPOLY.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts).
* \param	sPts			Source control points.
* \param        mesh                    Mesh which is used to compute
*                                       constrained distances. If non NULL
*                                       and the mesh type is
*                                       WLZ_CMESH_2D then  constrained
*                                       distances are used and these are
*                                       computed using the mesh.
*                                       If NULL or the transform is
*                                       some other type then Euclidean
*                                       distances are used.
* \param	nParam			Number of additional parameters.
* \param	param			Array of additional parameters.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts2DParam(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex2 *dPts,
					  int nSPts,
					  WlzDVertex2 *sPts,
					  WlzCMesh2D *mesh,
					  int nParam,
					  double *param,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	deltaMQ = 0.001,
		deltaIMQ = 0.300,
  		paramGauss = 0.9;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(mesh)
  {
    if(mesh->type != WLZ_CMESH_2D)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
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
					dPts, sPts,
					((nParam > 0) && (param != NULL))?
					*param: paramGauss,
					NULL, mesh, &errNum);
	break;
      case WLZ_FN_BASIS_2DPOLY:
	basisTr->basisFn = WlzBasisFnPoly2DFromCPts(nDPts,
					 order, dPts, sPts,
				         &errNum);
	break;
      case WLZ_FN_BASIS_2DIMQ:
	basisTr->basisFn = WlzBasisFnIMQ2DFromCPts(nDPts,
					dPts, sPts,
					((nParam > 0) && (param != NULL))?
					*param: deltaIMQ,
					NULL, mesh, &errNum);
	break;
      case WLZ_FN_BASIS_2DMQ:
	basisTr->basisFn = WlzBasisFnMQ2DFromCPts(nDPts,
					dPts, sPts,
					((nParam > 0) && (param != NULL))?
					*param: deltaMQ,
					NULL, mesh, &errNum);
	break;
      case WLZ_FN_BASIS_2DTPS:
	basisTr->basisFn = WlzBasisFnTPS2DFromCPts(nDPts,
					dPts, sPts, NULL, mesh,
					&errNum);
	break;
      case WLZ_FN_BASIS_2DCONF_POLY:
	basisTr->basisFn = WlzBasisFnConf2DFromCPts(nDPts,
					order, dPts, sPts,
				        &errNum);
	break;
      default:
	 errNum = WLZ_ERR_TRANSFORM_TYPE;
	 break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzBasisFnFreeTransform(basisTr);
      basisTr = NULL;
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
*		source verticies into an object with the given destination
*		verticies.
*		If a constraining object is given all distances will be
*		computed within the given object for those basis functions
*		which support constrained evaluation (Gauss, multi-quadric
*		and thin-plate spline).
*		Additional basis functions parameters may be supplied via
*		the nParam and param parameters. Currently this is only used to
*		supply the multi-quadric delta or gauss parameter scaling.
*		The default values of multi-quadric delta = 0.001 and
*		gauss param = 0.9 are used if nParam <= 0 or param == NULL.
* \param	type			Required basis function type.
* \param	order			Order of polynomial, only used for
* 					WLZ_FN_BASIS_3DPOLY.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control points
*					(must be same as nDPts).
* \param	sPts			Source control points.
* \param        mesh                    Mesh which is used to compute
*                                       constrained distances. If non NULL
*                                       and the mesh type is
*                                       WLZ_CMESH_3D then  constrained
*                                       distances are used and these are
*                                       computed using the mesh.
*                                       If NULL or the transform is
*                                       some other type then Euclidean
*                                       distances are used.
* \param	nParam			Number of additional parameters.
* \param	param			Array of additional parameters.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts3DParam(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex3 *dPts,
					  int nSPts,
					  WlzDVertex3 *sPts,
					  WlzCMesh3D *mesh,
					  int nParam,
					  double *param,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	deltaMQ = 0.001,
		deltaIMQ = 0.100;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(mesh)
  {
    if(mesh->type != WLZ_CMESH_3D)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
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
    basisTr->type = WLZ_TRANSFORM_3D_BASISFN;
    switch(type)
    {
      case WLZ_FN_BASIS_3DIMQ:
	basisTr->basisFn = WlzBasisFnIMQ3DFromCPts(nDPts,
					dPts, sPts,
					((nParam > 0) && (param != NULL))?
					*param: deltaIMQ,
					NULL, mesh, &errNum);
	break;
      case WLZ_FN_BASIS_3DMQ:
	basisTr->basisFn = WlzBasisFnMQ3DFromCPts(nDPts,
					dPts, sPts,
					((nParam > 0) && (param != NULL))?
					*param: deltaMQ,
					NULL, mesh, &errNum);
	break;
      default:
	 errNum = WLZ_ERR_TRANSFORM_TYPE;
	 break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzBasisFnFreeTransform(basisTr);
      basisTr = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisTr);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Changes control points in an existing basis function transform.
*		Using this function to add, move or delete control points
*		avoids recomputing distance transforms when using basis
*		functions which use constrained distances. Because distances
*		transforms are very expensive to compute calling this function
*		can be far more efficient, but when non-constrained (Euclidean)
*		distances are used then there is no benefit in using this
*		function as opposed to WlzBasisFnTPS2DFromCPts().
*		The full list of control points must be given.
* \param	basisTr			Existing basis function transform.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control
*					points (must be same as nDPts).
* \param	sPts			Source control points.
* \param	cObj			Constraining object, within which all
*					distances are constrained. If NULL
*					Euclidean distances are used in place
*					of constrained distances.
*/
WlzErrorNum	WlzBasisFnTPS2DChangeCPts(WlzBasisFnTransform *basisTr,
				int nDPts, WlzDVertex2 *dPts,
				int nSPts, WlzDVertex2 *sPts,
				WlzObject *cObj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzBasisFnTPS2DChangeCPtsParam(basisTr, nDPts, dPts, nSPts, sPts,
                                          cObj, 0, NULL);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Changes control points in an existing basis function transform.
*		Using this function to add, move or delete control points
*		avoids recomputing distance transforms when using basis
*		functions which use constrained distances. Because distances
*		transforms are expensive to compute calling this function
*		can be more efficient, but when non-constrained (Euclidean)
*		distances are used then there is no benefit in using this
*		function as opposed to WlzBasisFnTPS2DFromCPts().
*		The full list of control points must be given.
*		Additional basis functions parameters may be supplied via
*		the nParam and param parameters. Currently this is only used to
*		supply the multi-quadric delta or gauss parameter scaling.
*		The default values of multi-quadric delta = 0.001 and
*		gauss param = 0.9 are used if nParam <= 0 or param == NULL.
* \param	basisTr			Existing basis function transform.
* \param	nDPts			Number of destination control points.
* \param	dPts			Destination control points.
* \param	nSPts			Number of source control
*					points (must be same as nDPts).
* \param	sPts			Source control points.
* \param	cObj			Constraining object, within which all
*					distances are constrained. If NULL
*					Euclidean distances are used in place
*					of constrained distances.
* \param	nParam			Number of additional parameters.
* \param	param			Array of additional parameters.
*/
WlzErrorNum	WlzBasisFnTPS2DChangeCPtsParam(WlzBasisFnTransform *basisTr,
				int nDPts, WlzDVertex2 *dPts,
				int nSPts, WlzDVertex2 *sPts,
				WlzObject *cObj, int nParam,
				double *param)
{
  WlzBasisFn    *newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	deltaMQ = 0.001,
		deltaIMQ = 0.300,
  		paramGauss = 0.900;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(basisTr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(basisTr->type != WLZ_TRANSFORM_2D_BASISFN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(basisTr->basisFn->mesh.v != NULL)
  {
    if(basisTr->basisFn->mesh.m2->type != WLZ_CMESH_2D)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(basisTr->basisFn->type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
	newBasisFn = WlzBasisFnGauss2DFromCPts(nDPts, dPts, sPts,
			  ((nParam > 0) && (param != NULL))?
			  *param: paramGauss, basisTr->basisFn,
			  basisTr->basisFn->mesh.m2, &errNum);
	break;
      case WLZ_FN_BASIS_2DIMQ:
	newBasisFn = WlzBasisFnIMQ2DFromCPts(nDPts, dPts, sPts,
			  ((nParam > 0) && (param != NULL))?
			  *param: deltaIMQ, basisTr->basisFn,
			  basisTr->basisFn->mesh.m2, &errNum);
	break;
      case WLZ_FN_BASIS_2DMQ:
	newBasisFn = WlzBasisFnMQ2DFromCPts(nDPts, dPts, sPts,
			  ((nParam > 0) && (param != NULL))?
			  *param: deltaMQ, basisTr->basisFn,
			  basisTr->basisFn->mesh.m2, &errNum);
	break;
      case WLZ_FN_BASIS_2DTPS:
	newBasisFn = WlzBasisFnTPS2DFromCPts(nDPts, dPts, sPts,
					basisTr->basisFn,
					basisTr->basisFn->mesh.m2, &errNum);
	break;
      case WLZ_FN_BASIS_2DPOLY:
	newBasisFn = WlzBasisFnPoly2DFromCPts(nDPts, basisTr->basisFn->nPoly,
					 dPts, sPts,
				         &errNum);
	break;
      case WLZ_FN_BASIS_2DCONF_POLY:
	newBasisFn = WlzBasisFnConf2DFromCPts(nDPts, basisTr->basisFn->nPoly,
					dPts, sPts,
				        &errNum);
	break;
      default:
	 errNum = WLZ_ERR_TRANSFORM_TYPE;
	 break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )WlzBasisFnFree(basisTr->basisFn);
    basisTr->basisFn = newBasisFn;
  }
  return(errNum);
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
* \param        mesh                    Mesh which is used to compute
*                                       constrained distances. If non NULL
*                                       and the mesh type is
*                                       WLZ_CMESH_3D then  constrained
*                                       distances are used and these are
*                                       computed using the mesh.
*                                       If NULL or the transform is
*                                       some other type then Euclidean
*                                       distances are used.
* \param	dstErr			Destination error pointer,
*					may be NULL.
*/
WlzBasisFnTransform *WlzBasisFnTrFromCPts3D(WlzFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex3 *dPts,
					  int nSPts,
					  WlzDVertex3 *sPts,
					  WlzCMesh3D *mesh,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basisTr = NULL;

  basisTr = WlzBasisFnTrFromCPts3DParam(type, order, nDPts, dPts,
                                        nSPts, sPts, mesh, 0, NULL,
					dstErr);
  return(basisTr);
}

/*!
* \return	Error number.
* \ingroup	WlzTransform
* \brief	Sets the displacements of the given mesh transform according
* 		to the basis function transform.
* \param	mesh			Given mesh transform.
* \param	basisTr			Given basis function transform.
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
      case WLZ_FN_BASIS_2DIMQ:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnValueIMQ2D(basisTr->basisFn,
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
*		This function just calls either WlzBasisFnSetCMesh2D()
*		or WlzBasisFnSetCMesh3D() see these functions for their
*		conforming mesh transform object requirements.
* \param	mObj			Given conforming mesh transform object.
* \param	basisTr			Given basis function transform.
*/
WlzErrorNum    	WlzBasisFnSetCMesh(WlzObject *mObj,
				   WlzBasisFnTransform *basisTr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((mObj->domain.core == NULL) || (basisTr == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mObj->type)
    {
      case WLZ_CMESH_2D:
        errNum = WlzBasisFnSetCMesh2D(mObj, basisTr);
	break;
      case WLZ_CMESH_3D:
        errNum = WlzBasisFnSetCMesh3D(mObj, basisTr);
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
* \brief	Sets the displacements of the given 2D conforming mesh
*		transform object according to the basis function transform.
*		The conforming mesh object must have a valid 2D conforming
*		mesh and indexed values. The indexed values will be expanded
*		to cover the nodes of the mesh if required, but they must
*		have a rank of 1, a dimension of >= 2 and be of type double.
* \param	mObj			Given mesh transform object.
* \param	basisTr			Given basis function transform.
*/
WlzErrorNum    	WlzBasisFnSetCMesh2D(WlzObject *mObj,
				     WlzBasisFnTransform *basisTr)
{
  int		idN,
  		maxNodIdx;
  double	*dsp;
  WlzDVertex2	dspV;
  WlzCMeshNod2D	*nod;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((basisTr == NULL) || ((mesh = mObj->domain.cm2) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((basisTr->type != WLZ_TRANSFORM_2D_BASISFN) ||
          (mesh->type != WLZ_CMESH_2D))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((ixv = mObj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(ixv->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv->rank != 1) || (ixv->dim[0] < 2) ||
          (ixv->vType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    maxNodIdx = mesh->res.nod.maxEnt;
    ixv->attach = WLZ_VALUE_ATTACH_NOD;
    if(WlzIndexedValueExtGet(ixv, maxNodIdx) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(basisTr->basisFn->type)
    {
      case WLZ_FN_BASIS_2DGAUSS:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
	for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueGauss2D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	  }
	}
        break;
      case WLZ_FN_BASIS_2DIMQ:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
	for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueIMQ2D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	  }
	}
        break;
      case WLZ_FN_BASIS_2DMQ:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
	for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueMQ2D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	  }
	}
        break;
      case WLZ_FN_BASIS_2DTPS:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
	for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueTPS2D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	  }
	}
        break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return	Error number.
* \ingroup	WlzTransform
* \brief	Sets the displacements of the given 3D conforming mesh
*		transform object according to the basis function transform.
*		The conforming mesh object must have a valid 3D conforming
*		mesh and indexed values. The indexed values will be expanded
*		to cover the nodes of the mesh if required, but they must
*		have a rank of 1, a dimension of >= 3 and be of type double.
* \param	mObj			Given mesh transform object.
* \param	basisTr			Given basis function transform.
*/
WlzErrorNum    	WlzBasisFnSetCMesh3D(WlzObject *mObj,
				     WlzBasisFnTransform *basisTr)
{
  int		idN,
  		maxNodIdx;
  double	*dsp;
  WlzDVertex3	dspV;
  WlzCMeshNod3D	*nod;
  WlzCMesh3D	*mesh;
  WlzIndexedValues *ixv;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((basisTr == NULL) || ((mesh = mObj->domain.cm3) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((basisTr->type != WLZ_TRANSFORM_3D_BASISFN) ||
          (mesh->type != WLZ_CMESH_3D))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else if((ixv = mObj->values.x) == NULL)
  {
    errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(ixv->type != WLZ_INDEXED_VALUES)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else if((ixv->rank != 1) || (ixv->dim[0] < 3) ||
          (ixv->vType != WLZ_GREY_DOUBLE))
  {
    errNum = WLZ_ERR_VALUES_DATA;
  }
  else
  {
    maxNodIdx = mesh->res.nod.maxEnt;
    ixv->attach = WLZ_VALUE_ATTACH_NOD;
    if(WlzIndexedValueExtGet(ixv, maxNodIdx) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(basisTr->basisFn->type)
    {
      case WLZ_FN_BASIS_3DIMQ:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
        for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueIMQ3D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	    dsp[2] = dspV.vtZ;
	  }
	}
	break;
      case WLZ_FN_BASIS_3DMQ:
#ifdef _OPENMP
#pragma omp parallel for private(dsp, dspV, nod)
#endif
        for(idN = 0; idN < maxNodIdx; ++idN)
	{
	  nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	  if(nod->idx >= 0)
	  {
	    dsp = (double *)WlzIndexedValueGet(ixv, idN);
	    dspV = WlzBasisFnValueMQ3D(basisTr->basisFn, nod->pos);
	    dsp[0] = dspV.vtX;
	    dsp[1] = dspV.vtY;
	    dsp[2] = dspV.vtZ;
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_VALUES_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Uses the given source mesh to create a mesh transform
* 		object which transforms the source to target. See the
* 		functions WlzBasisFnMakeCMeshTr2D() and
* 		WlzBasisFnMakeCMeshTr3D() for details of the
* 		returned objects.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnMakeCMeshTr(WlzBasisFnTransform *basisTr,
					WlzCMeshP mesh,
				      	WlzErrorNum *dstErr)
{
  WlzObject	*mObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((basisTr == NULL) || (mesh.v == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(basisTr->type)
    {
      case WLZ_TRANSFORM_2D_BASISFN:
	if(mesh.m2->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          mObj = WlzBasisFnMakeCMeshTr2D(basisTr, mesh.m2, &errNum);
	}
	break;
      case WLZ_TRANSFORM_3D_BASISFN:
	if(mesh.m3->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          mObj = WlzBasisFnMakeCMeshTr3D(basisTr, mesh.m3, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Uses the given 2D target mesh to create a 2D mesh transform
* 		which transforms the source to target. Unlike
* 		WlzBasisFnMakeCMeshTr() this function does not check
* 		it's given parameters. The new object's indexed values
* 		have; rank = 1, dim = 2, attach = WLZ_VALUE_ATTACH_NOD
* 		and double values.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnMakeCMeshTr2D(WlzBasisFnTransform *basisTr,
				        WlzCMesh2D *mesh,
				        WlzErrorNum *dstErr)
{
  int		dim = 2;
  WlzObject	*mObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  if(mesh->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    dom.cm2 = mesh;
    mObj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.x = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
				    WLZ_VALUE_ATTACH_NOD, &errNum);
    mObj->values = WlzAssignValues(val, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzBasisFnSetCMesh2D(mObj, basisTr);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(mObj != NULL)
    {
      (void )WlzFreeObj(mObj);
      mObj = NULL;
    }
    else
    {
      (void )WlzFreeIndexedValues(val.x);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Uses the given 3D target mesh to create a 3D mesh transform
* 		which transforms the source to target. Unlike
* 		WlzBasisFnMakeCMeshTr() this function does not check
* 		it's given parameters. The new object's indexed values
* 		have; rank = 1, dim = 3, attach = WLZ_VALUE_ATTACH_NOD
* 		and double values.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnMakeCMeshTr3D(WlzBasisFnTransform *basisTr,
				        WlzCMesh3D *mesh,
				        WlzErrorNum *dstErr)
{
  int		dim = 3;
  WlzObject	*mObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  if(mesh->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    dom.cm3 = mesh;
    mObj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.x = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
				    WLZ_VALUE_ATTACH_NOD, &errNum);
    mObj->values = WlzAssignValues(val, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzBasisFnSetCMesh3D(mObj, basisTr);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(mObj != NULL)
    {
      (void )WlzFreeObj(mObj);
      mObj = NULL;
    }
    else
    {
      (void )WlzFreeIndexedValues(val.x);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Copies the given target mesh and uses it to create a
* 		mesh transform object which transforms the source to
* 		target. See the functions WlzBasisFnInvertMakeCMeshTr2D()
* 		and WlzBasisFnInvertMakeCMeshTr3D() for details of the
* 		returned objects.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnInvertMakeCMeshTr(WlzBasisFnTransform *basisTr,
					WlzCMeshP mesh,
				      	WlzErrorNum *dstErr)
{
  WlzObject	*mObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((basisTr == NULL) || (mesh.v == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(basisTr->type)
    {
      case WLZ_TRANSFORM_2D_BASISFN:
	if(mesh.m2->type != WLZ_CMESH_2D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          mObj = WlzBasisFnInvertMakeCMeshTr2D(basisTr, mesh.m2, &errNum);
	}
	break;
      case WLZ_TRANSFORM_3D_BASISFN:
	if(mesh.m3->type != WLZ_CMESH_3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
          mObj = WlzBasisFnInvertMakeCMeshTr3D(basisTr, mesh.m3, &errNum);
	}
	break;
      default:
        errNum = WLZ_ERR_TRANSFORM_TYPE;
	break;
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Copies the given 2D target mesh and uses it to create a
* 		2D mesh transform which transforms the source to target.
* 		Unlike WlzBasisFnInvertMakeCMeshTr() this function does
* 		not check it's given parameters. The new object's indexed
* 		values have; rank = 1, dim = 2, attach = WLZ_VALUE_ATTACH_NOD
* 		and double values.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnInvertMakeCMeshTr2D(
					WlzBasisFnTransform *basisTr,
				        WlzCMesh2D *mesh,
				        WlzErrorNum *dstErr)
{
  int		dim = 2;
  WlzObject	*mObj = NULL,
  		*rObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  dom.cm2 = WlzCMeshCopy2D(mesh, 1, 0, NULL, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.x = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
				    WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj->values = WlzAssignValues(val, NULL);
    errNum = WlzBasisFnSetCMesh2D(mObj, basisTr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzCMeshTransformInvert(mObj, &errNum);
  }
  (void )WlzFreeObj(mObj);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New constrained mesh transform object.
* \ingroup	WlzTransform
* \brief	Copies the given 3D target mesh and uses to create a
* 		3D mesh transform which transforms the source to target.
* 		Unlike WlzBasisFnInvertMakeCMeshTr() this function does
* 		not check it's given parameters. The new object's indexed
* 		values have; rank = 1, dim = 3, attach = WLZ_VALUE_ATTACH_NOD
* 		and double values.
* \param	basisTr			Given basis function transform
* 					which transforms target to source.
* \param	mesh			Given conforming mesh for target.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 	*WlzBasisFnInvertMakeCMeshTr3D(
					WlzBasisFnTransform *basisTr,
				        WlzCMesh3D *mesh,
				        WlzErrorNum *dstErr)
{
  int		dim = 3;
  WlzObject	*mObj = NULL,
  		*rObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  val.core = NULL;
  dom.cm3 = WlzCMeshCopy3D(mesh, 1, 0, NULL, NULL, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    val.x = WlzMakeIndexedValues(mObj, 1, &dim, WLZ_GREY_DOUBLE,
				    WLZ_VALUE_ATTACH_NOD, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj->values = WlzAssignValues(val, NULL);
    errNum = WlzBasisFnSetCMesh3D(mObj, basisTr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzCMeshTransformInvert(mObj, &errNum);
  }
  (void )WlzFreeObj(mObj);
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(rObj);
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
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
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
	  dVP.i2->vtX = (int )(tD.vtX);
	  dVP.i2->vtY = (int )(tD.vtY);
	  for(idN = 1; idN < srcPoly->nvertices; ++idN)
	  {
	    ++(sVP.i2);
	    ++(dVP.i2);
	    tD.vtX = sVP.i2->vtX;
	    tD.vtY = sVP.i2->vtY;
	    tD = WlzBasisFnTransformVertexD(basisTr, tD, NULL);
	    dVP.i2->vtX = (int )(tD.vtX);
	    dVP.i2->vtY = (int )(tD.vtY);
	  }
	}
        break;
      case WLZ_POLYGON_FLOAT:
	tD.vtX = sVP.f2->vtX;
	tD.vtY = sVP.f2->vtY;
	tD = WlzBasisFnTransformVertexD(basisTr, tD, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  dVP.f2->vtX = (float )(tD.vtX);
	  dVP.f2->vtY = (float )(tD.vtY);
	  for(idN = 1; idN < srcPoly->nvertices; ++idN)
	  {
	    ++(sVP.f2);
	    ++(dVP.f2);
	    tD.vtX = sVP.f2->vtX;
	    tD.vtY = sVP.f2->vtY;
	    tD = WlzBasisFnTransformVertexD(basisTr, tD, NULL);
	    dVP.f2->vtX = (float )(tD.vtX);
	    dVP.f2->vtY = (float )(tD.vtY);
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
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreePolyDmn(dstPoly);
      dstPoly = NULL;
    }
  }
  if(newPoly == 0)
  {
    dstPoly = srcPoly;
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

  dNr.vtX = 0.0;
  dNr.vtY = 0.0;
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

  dstVx.vtX = 0.0;
  dstVx.vtY = 0.0;
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
      case WLZ_FN_BASIS_2DIMQ:
	dstVx = WlzBasisFnValueIMQ2D(basisTr->basisFn, srcVx);
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
  dstVxF.vtX = (float )(dstVxD.vtX);
  dstVxF.vtY = (float )(dstVxD.vtY);
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
