#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzBasisFnTransform.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Woolz functions for computing and applying basis
*		function transforms.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

static WlzBasisFnTransform *WlzBasisFnGaussFromCPts(int, WlzDVertex2 *,
				WlzDVertex2 *, double, WlzErrorNum *),
		*WlzBasisFnPolyFromCPts(int, int, WlzDVertex2 *, 
				WlzDVertex2 *, WlzErrorNum *),
		*WlzBasisFnMQFromCPts(int, WlzDVertex2 *, 
				WlzDVertex2 *, double, WlzErrorNum *),
		*WlzBasisFnTPSFromCPts(int, WlzDVertex2 *, 
				WlzDVertex2 *, WlzErrorNum *);
static WlzDVertex2 WlzBasisFnDisplacementGauss(WlzBasisFnTransform *,
				WlzDVertex2),
		WlzBasisFnDisplacementPoly(WlzBasisFnTransform *,
				WlzDVertex2),
		WlzBasisFnDisplacementMQ(WlzBasisFnTransform *,
				WlzDVertex2),
		WlzBasisFnDisplacementTPS(WlzBasisFnTransform *,
				WlzDVertex2),
		WlzBasisFnDispRedPoly(WlzDVertex2 *, WlzDVertex2);
static void	WlzBasisFnVxExtent(WlzDBox2 *, WlzDVertex2 *, WlzDVertex2 *,
				int),
		WlzBasisFnGaussCoeff(WlzBasisFnTransform *, double *, int),
		WlzBasisFnMQCoeff(WlzBasisFnTransform *, double *,
				WlzDBox2 *, double, int),
		WlzBasisFnTPSCoeff(WlzBasisFnTransform *, double *,
				WlzDBox2 *, double, int);
static WlzBasisFnTransform *WlzBasisFnConfFromCPts(int nPts, int order,
						   WlzDVertex2 *dPts,
						   WlzDVertex2 *sPts,
						   WlzErrorNum *dstErr);
static WlzDVertex2 WlzBasisFnDisplacementConf(WlzBasisFnTransform *basis,
					     WlzDVertex2 srcVx);
				   

/************************************************************************
* Function:	WlzBasisFnFreeTransform					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Free's the given basis function transform.		*
* Global refs:	-							*
* Parameters:	WlzBasisFnTransform *basis: Given basis function	*
*					transform.			*
************************************************************************/
WlzErrorNum	WlzBasisFnFreeTransform(WlzBasisFnTransform *basis)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basis)
  {
    if(basis->poly)
    {
      AlcFree(basis->poly);
    }
    if(basis->basis)
    {
      AlcFree(basis->basis);
    }
    if(basis->verticies)
    {
      AlcFree(basis->verticies);
    }
    AlcFree(basis);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzBasisFnTrFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new basis function transform of the given	*
*		type, which will transform an object with the given	*
*		source verticies into an object with the given 		*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	WlzBasisFnType type:	Required basis function type.	*
*		int order:		Order of polynomial, only 	*
*					used for WLZ_BASISFN_POLY.	*
*		int nDPts:		Number of destination control	*
*					points.				*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		int nSPts:		Number of source control points	*
*					(must be same as nDPts).	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
					may be NULL.			*
************************************************************************/
WlzBasisFnTransform *WlzBasisFnTrFromCPts(WlzBasisFnType type,
					  int order,
					  int nDPts,
					  WlzDVertex2 *dPts,
					  int nSPts,
					  WlzDVertex2 *sPts,
					  WlzErrorNum *dstErr)
{
  WlzBasisFnTransform *basis = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nDPts != nSPts) || (nDPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((dPts == NULL) || (sPts == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(type)
    {
      case WLZ_BASISFN_GAUSS:
	basis = WlzBasisFnGaussFromCPts(nDPts, dPts, sPts, 0.9, &errNum);
	break;
      case WLZ_BASISFN_POLY:
	basis = WlzBasisFnPolyFromCPts(nDPts, order, dPts, sPts, &errNum);
	break;
      case WLZ_BASISFN_MQ:
	basis = WlzBasisFnMQFromCPts(nDPts, dPts, sPts, 0.1, &errNum);
	break;
      case WLZ_BASISFN_TPS:
	basis = WlzBasisFnTPSFromCPts(nDPts, dPts, sPts, &errNum);
	break;
      case WLZ_BASISFN_CONF_POLY:
	basis = WlzBasisFnConfFromCPts(nDPts, order, dPts, sPts, &errNum);
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
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnSetMesh					*
* Returns:	WlzErrorNum:		Error number.			*
* Purpose:	Sets the displacements of the given mesh transform	*
*		according to the basis function transform.		*
* Global refs:	-							*
* Parameters:	WlzMeshTransform *mesh:	Given mesh transform.		*
*		WlzBasisFnTransform *basis: Given basis function	*
*					transform.			*
************************************************************************/
WlzErrorNum    	WlzBasisFnSetMesh(WlzMeshTransform *mesh,
				  WlzBasisFnTransform *basis)
{
  int		nodCnt;
  WlzMeshNode	*nod;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((mesh == NULL) || (basis == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((mesh->type != WLZ_TRANSFORM_2D_MESH) ||
	  (basis->type != WLZ_TRANSFORM_2D_BASISFN))
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    nod = mesh->nodes;
    nodCnt = mesh->nNodes;
    switch(basis->basisFn)
    {
      case WLZ_BASISFN_GAUSS:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnDisplacementGauss(basis,
	  						  nod->position);
	  ++nod;
	}
	break;
      case WLZ_BASISFN_POLY:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnDisplacementPoly(basis,
						         nod->position);
	  ++nod;
	}
	break;
      case WLZ_BASISFN_MQ:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnDisplacementMQ(basis,
						       nod->position);
	  ++nod;
	}
	break;
      case WLZ_BASISFN_TPS:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnDisplacementTPS(basis,
							nod->position);
	  ++nod;
	}
	break;
    case WLZ_BASISFN_CONF_POLY:
	while(nodCnt-- > 0)
	{
	  nod->displacement = WlzBasisFnDisplacementConf(basis,
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

/************************************************************************
* Function:	WlzBasisFnTransformObj					*
* Returns:	WlzObject *:		Transformed object, NULL on	*
*					error.				*
* Purpose:	Transforms a woolz object using a the given basis	*
*		function transform.					*
*		This function has been written as an example of how	*
*		to transform an object using a basis function and mesh.	*
*		In most cases WlzMeshFromObj(), WlzBasisFnSetMesh()	*
*		and WlzMeshTransformObj() would be called allowing a	*
*		mesh to be reused.					*
* Global refs:	-							*
* Parameters:	WlzObject *srcObj:	Object to be transformed.	*
*		WlzBasisFnTransform *mesh: Basis function transform	*
*					to apply.			*
*		WlzInterpolationType interp: Level of interpolation to	*
*					use.				*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzBasisFnTransformObj(WlzObject *srcObj,
					WlzBasisFnTransform *basis,
					WlzInterpolationType interp,
					WlzErrorNum *dstErr)
{
  WlzObject	*dstObj = NULL;
  WlzMeshTransform *mesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((srcObj == NULL) || (basis == NULL))
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
    errNum = WlzBasisFnSetMesh(mesh, basis);
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

/************************************************************************
* Function:     WlzBasisFnTransformVertexD                              *
* Returns:      WlzDVertex2:             Transformed vertex.             *
* Purpose:      Transforms the given WlzDVertex2.                        *
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform to	*
*					apply.				*
*               WlzDVertex2 srcVx:	Vertex to be transformed.       *
*               WlzErrorNum *dstErr:    Destination pointer for error   *
*                                       number, may be NULL.            *
************************************************************************/
WlzDVertex2	WlzBasisFnTransformVertexD(WlzBasisFnTransform *basis,
					   WlzDVertex2 srcVx,
					   WlzErrorNum *dstErr)
{
  WlzDVertex2	dstVx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basis == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(basis->type != WLZ_TRANSFORM_2D_BASISFN)
  {
    errNum = WLZ_ERR_TRANSFORM_TYPE;
  }
  else
  {
    switch(basis->basisFn)
    {
      case WLZ_BASISFN_GAUSS:
	dstVx = WlzBasisFnDisplacementGauss(basis, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_BASISFN_POLY:
	dstVx = WlzBasisFnDisplacementPoly(basis, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_BASISFN_MQ:
	dstVx = WlzBasisFnDisplacementMQ(basis, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_BASISFN_TPS:
	dstVx = WlzBasisFnDisplacementTPS(basis, srcVx);
	dstVx.vtX += srcVx.vtX;
	dstVx.vtY += srcVx.vtY;
	break;
      case WLZ_BASISFN_CONF_POLY:
	dstVx = WlzBasisFnDisplacementConf(basis, srcVx);
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

/************************************************************************
* Function:     WlzBasisFnTransformVertexF                              *
* Returns:      WlzFVertex2:             Transformed vertex.             *
* Purpose:      Transforms the given WlzFVertex2.                        *
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform to	*
*					apply.				*
*               WlzFVertex2 srcVxF:	Vertex to be transformed.       *
*               WlzErrorNum *dstErr:    Destination pointer for error   *
*                                       number, may be NULL.            *
************************************************************************/
WlzFVertex2	WlzBasisFnTransformVertexF(WlzBasisFnTransform *basis,
					   WlzFVertex2 srcVxF,
					   WlzErrorNum *dstErr)
{
  WlzFVertex2	dstVxF;
  WlzDVertex2	dstVxD,
		srcVxD;

  srcVxD.vtX = srcVxF.vtX;
  srcVxD.vtY = srcVxF.vtY;
  dstVxD = WlzBasisFnTransformVertexD(basis, srcVxD, dstErr);
  dstVxF.vtX = dstVxD.vtX;
  dstVxF.vtY = dstVxD.vtY;
  return(dstVxF);
}

/************************************************************************
* Function:     WlzBasisFnTransformVertexI                              *
* Returns:      WlzIVertex2:             Transformed vertex.             *
* Purpose:      Transforms the given WlzIVertex2.                        *
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform to	*
*					apply.				*
*               WlzIVertex2 srcVxI:	Vertex to be transformed.       *
*               WlzErrorNum *dstErr:    Destination pointer for error   *
*                                       number, may be NULL.            *
************************************************************************/
WlzIVertex2	WlzBasisFnTransformVertexI(WlzBasisFnTransform *basis,
					   WlzIVertex2 srcVxI,
					   WlzErrorNum *dstErr)
{
  WlzIVertex2	dstVxI;
  WlzDVertex2	dstVxD,
		srcVxD;

  srcVxD.vtX = srcVxI.vtX;
  srcVxD.vtY = srcVxI.vtY;
  dstVxD = WlzBasisFnTransformVertexD(basis, srcVxD, dstErr);
  dstVxI.vtX = WLZ_NINT(dstVxD.vtX);
  dstVxI.vtY = WLZ_NINT(dstVxD.vtY);
  return(dstVxI);
}

/************************************************************************
* Function:     WlzBasisFnDisplacementPoly                              *
* Returns:      WlzDVertex2:             Displacement vertex.            *
* Purpose:      Calculates the displacement for the given vertex using	*
*		a polynomial basis function.				*
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform.	*
*               WlzDVertex2 srcVx:	Source vertex.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDisplacementPoly(WlzBasisFnTransform *basis,
					     WlzDVertex2 srcVx)
{
  int		idX,
		idY;
  double	tD0;
  WlzDVertex2	*polyP;
  WlzDVertex2	powVx,
  		dspVx;


  powVx.vtY = 1.0;
  polyP = basis->poly;
  dspVx.vtY = dspVx.vtX = 0.0;
  for(idY = 0; idY <= basis->nPoly; ++idY)
  {
    powVx.vtX = 1.0;
    for(idX = 0; idX <= basis->nPoly; ++idX)
    {
      tD0 = powVx.vtX * powVx.vtY;
      dspVx.vtX += polyP->vtX * tD0;
      dspVx.vtY += polyP->vtY * tD0;
      powVx.vtX *= srcVx.vtX;
      ++polyP;
    }
    powVx.vtY *= srcVx.vtY;
  }
  return(dspVx);
}

/************************************************************************
* Function:     WlzBasisFnDisplacementGauss				*
* Returns:      WlzDVertex2:             Displacement vertex.            *
* Purpose:      Calculates the displacement for the given vertex using	*
*		a Gaussian basis function.				*
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform.	*
*               WlzDVertex2 srcVx:	Source vertex.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDisplacementGauss(WlzBasisFnTransform *basis,
					      WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		dspVx;

  dspVx.vtX = 0.0;
  dspVx.vtY = 0.0;
  count = basis->nVtx;
  cPts = basis->verticies;
  basisCo = basis->basis;
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    tD1 = (tD0 > DBL_EPSILON)? exp(tD0 * basis->delta): 1.0;
    dspVx.vtX += basisCo->vtX * tD1;
    dspVx.vtY += basisCo->vtY * tD1;
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnDispRedPoly(basis->poly, srcVx);
  dspVx.vtX = dspVx.vtX + polyVx.vtX;
  dspVx.vtY = dspVx.vtY + polyVx.vtY;
  return(dspVx);
}

/************************************************************************
* Function:     WlzBasisFnDisplacementMQ				*
* Returns:      WlzDVertex2:             Displacement vertex.            *
* Purpose:      Calculates the displacement for the given vertex using	*
*		a multiquadric basis function.				*
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform.	*
*               WlzDVertex2 srcVx:	Source vertex.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDisplacementMQ(WlzBasisFnTransform *basis,
					   WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		dspVx;

  dspVx.vtX = 0.0;
  dspVx.vtY = 0.0;
  count = basis->nVtx;
  cPts = basis->verticies;
  basisCo = basis->basis;
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    if(tD0 > DBL_EPSILON)
    {
      tD0 = sqrt(tD0 + basis->delta);
      dspVx.vtX += basisCo->vtX * tD0;
      dspVx.vtY += basisCo->vtY * tD0;
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnDispRedPoly(basis->poly, srcVx);
  dspVx.vtX = dspVx.vtX + polyVx.vtX;
  dspVx.vtY = dspVx.vtY + polyVx.vtY;
  return(dspVx);
}

/************************************************************************
* Function:     WlzBasisFnDisplacementTPS                               *
* Returns:      WlzDVertex2:             Displacement vertex.            *
* Purpose:      Calculates the displacement for the given vertex using	*
*		a thin plate spline basis function.			*
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform.	*
*               WlzDVertex2 srcVx:	Source vertex.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDisplacementTPS(WlzBasisFnTransform *basis,
					    WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		dspVx;

  dspVx.vtX = 0.0;
  dspVx.vtY = 0.0;
  count = basis->nVtx;
  cPts = basis->verticies;
  basisCo = basis->basis;
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    if(tD0 > DBL_EPSILON)
    {
      tD0 *= log(tD0);
      dspVx.vtX += basisCo->vtX * tD0;
      dspVx.vtY += basisCo->vtY * tD0;
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnDispRedPoly(basis->poly, srcVx);
  dspVx.vtX = (dspVx.vtX / 2) + polyVx.vtX;
  dspVx.vtY = (dspVx.vtY / 2) + polyVx.vtY;
  return(dspVx);
}

/************************************************************************
* Function:     WlzBasisFnDisplacementConf                              *
* Returns:      WlzDVertex2:             Displacement vertex.           *
* Purpose:      Calculates the displacement for the given vertex using	*
*		a conformal polynomial basis function.			*
* Global refs:  -                                                       *
* Parameters:   WlzBasisFnTransform *basis: Basis function transform.	*
*               WlzDVertex2 srcVx:	Source vertex.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDisplacementConf(WlzBasisFnTransform *basis,
					     WlzDVertex2 srcVx)
{
  int		i;
  ComplexD	z, w, powW, a, b;
  WlzDVertex2	dspVx;
  WlzDVertex2	*polyP;

  polyP = basis->poly;
  w.re = srcVx.vtX;
  w.im = srcVx.vtY;
  z.re = polyP[0].vtX;
  z.im = polyP[0].vtY;
  powW.re = 1.0;
  powW.im = 0.0;

  for(i=1; i <= basis->nPoly; i++){
    powW = AlgCMult(powW, w);
    a.re = polyP[i].vtX;
    a.im = polyP[i].vtY;
    b = AlgCMult(a, powW);
    z.re += b.re;
    z.im += b.im;
  }
  dspVx.vtX = z.re;
  dspVx.vtY = z.im;

  return(dspVx);
}

/************************************************************************
* Function:	WlzBasisFnDispRedPoly					*
* Returns:	WlzDVertex2:		Displacement due to poly.	*
* Purpose:	Computes the displacement due to the reduced polynomial	*
*		used by the TPS, MQ and Gauss basis functions.		*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
					may be NULL.			*
************************************************************************/
static WlzDVertex2 WlzBasisFnDispRedPoly(WlzDVertex2 *poly, WlzDVertex2 srcVx)
{
  WlzDVertex2	dspVx;

  dspVx.vtX = poly->vtX;
  dspVx.vtY = poly->vtY;
  ++poly;
  dspVx.vtX += poly->vtX * srcVx.vtX;
  dspVx.vtY += poly->vtY * srcVx.vtX;
  ++poly;
  dspVx.vtX += poly->vtX * srcVx.vtY;
  dspVx.vtY += poly->vtY * srcVx.vtY;
  return(dspVx);
}

/************************************************************************
* Function:	WlzBasisFnGaussFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new Gaussian basis function transform,	*
*		which will transform an object with the given		*
*		source verticies into an object with the given 		*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		double delta:		Normalized delta value in range	*
*					[> 0.0 , < 1.0 ].		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzBasisFnTransform *WlzBasisFnGaussFromCPts(int nPts,
						    WlzDVertex2 *dPts,
						    WlzDVertex2 *sPts,
						    double delta,
						    WlzErrorNum *dstErr)
{
  int		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys;
  double	tD0,
		tD1,
		deltaSq,
		range,
		thresh,
		wMax;
  double	*bMx = NULL,
  		*wMx = NULL;
  double	**aMx = NULL,
  		**vMx = NULL;
  WlzDVertex2	tDVx0;
  WlzDBox2	extentDB;
  WlzBasisFnTransform *basis = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  deltaSq = delta * delta;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basis = (WlzBasisFnTransform *)
	       AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL) ||
     ((basis->poly = (WlzDVertex2 *)
		     AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basis->basis = (WlzDVertex2 *)
		      AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basis->verticies = (WlzDVertex2 *)
			  AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent(&extentDB, dPts, sPts, nPts);
    tD0 = extentDB.xMax - extentDB.xMin;
    tD1 = extentDB.yMax - extentDB.yMin;
    range = (tD0 > tD1)? tD0: tD1;
    if(range <= 1.0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    basis->type = WLZ_TRANSFORM_2D_BASISFN;
    basis->linkcount = 0;
    basis->freeptr = NULL;
    basis->basisFn = WLZ_BASISFN_GAUSS;
    basis->nPoly = 2;
    basis->nBasis = nPts;
    basis->nVtx = nPts;
    basis->delta = deltaSq / (range * range);
    WlzValueCopyDVertexToDVertex(basis->verticies, dPts, nPts);
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aMx + idY) + idX) = 0.0;
      }
      *(bMx + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      tD1 = (dPts + idY)->vtY;
      *(bMx + idY3) = (sPts + idY)->vtX - tD0;
      *(*(aMx + idY3) + 0) = 1.0;
      *(*(aMx + idY3) + 1) = tD0;
      *(*(aMx + idY3) + 2) = tD1;
      *(*(aMx + 0) + idY3) = 1.0;
      *(*(aMx + 1) + idY3) = tD0;
      *(*(aMx + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	tDVx0.vtX = (dPts + idX)->vtX - (dPts + idY)->vtX;
	tDVx0.vtX *= tDVx0.vtX;
	tDVx0.vtY = (dPts + idX)->vtY - (dPts + idY)->vtY;
	tDVx0.vtY *= tDVx0.vtY;
	tD0 = (tDVx0.vtX + tDVx0.vtY) * basis->delta;
	tD1 = (tD0 > DBL_EPSILON)? exp(tD0): 1.0;
	idX3 = idX + 3;
	*(*(aMx + idY3) + idX3) = tD1;
	*(*(aMx + idX3) + idY3) = tD1;
      }
      *(*(aMx + idY3) + idY3) = 1.0;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aMx, nSys, nSys, wMx, vMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    wMax = 0.0;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) > wMax)
      {
	wMax = *(wMx + idN);
      }
    }
    thresh = tol * wMax;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) < thresh)
      {
	*(wMx + idN) = 0.0;
      }
    }
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
    						vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
       and the y polynomial coefficients. */
    WlzBasisFnGaussCoeff(basis, bMx, 1);
    *(bMx + 0) = 0.0;
    *(bMx + 1) = 0.0;
    *(bMx + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      *(bMx + idY3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
    						vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients. */
    WlzBasisFnGaussCoeff(basis, bMx, 0);
  }
  if(bMx)
  {
    AlcFree(bMx);
  }
  if(wMx)
  {
    AlcFree(wMx);
  }
  if(aMx)
  {
    (void )AlcDouble2Free(aMx);
  }
  if(vMx)
  {
    (void )AlcDouble2Free(vMx);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(basis)
    {
      if(basis->verticies)
      {
	AlcFree(basis->verticies);
      }
      if(basis->basis)
      {
	AlcFree(basis->basis);
      }
      if(basis->poly)
      {
	AlcFree(basis->poly);
      }
      AlcFree(basis);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnPolyFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new polynomial basis function transform,	*
*		which will transform an object with the given		*
*		source verticies into an object with the given 		*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		int order:		Order of polynomial.		*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzBasisFnTransform *WlzBasisFnPolyFromCPts(int nPts, int order,
						   WlzDVertex2 *dPts,
						   WlzDVertex2 *sPts,
						   WlzErrorNum *dstErr)
{
  int  		idM,
  		idN,
		idX,
  		idY,
		nCoef;
  double	thresh,
  		wMax;
  double	*bMx = NULL,
  		*wMx = NULL;
  double	**aMx = NULL,
  		**vMx = NULL;
  WlzDVertex2	powVx,
  		sVx;
  WlzBasisFnTransform *basis = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  if((order < 0) || (nPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((basis = (WlzBasisFnTransform *)
	           AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    basis->type = WLZ_TRANSFORM_2D_BASISFN;
    basis->linkcount = 0;
    basis->freeptr = NULL;
    basis->basisFn = WLZ_BASISFN_POLY;
    basis->nPoly = order;
    basis->nBasis = 0;
    basis->nVtx = 0;
    basis->delta = 0.0;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nCoef = (order + 1) * (order + 1);
    if(((wMx = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bMx = (double *)AlcMalloc(sizeof(double) * nPts)) == NULL) ||
       (AlcDouble2Malloc(&vMx, nPts, nCoef) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&aMx, nPts, nCoef) !=  ALC_ER_NONE) ||
       ((basis->poly = (WlzDVertex2 *)
		       AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill matrix A. */
    for(idM = 0; idM < nPts; ++idM)
    {
      idN = 0;
      powVx.vtY = 1.0;
      sVx = *(sPts + idM);
      for(idY = 0; idY <= basis->nPoly; ++idY)
      {
	powVx.vtX = 1.0;
	for(idX = 0; idX <= basis->nPoly; ++idX)
	{
	  *(*(aMx + idM) + idN++) = powVx.vtX * powVx.vtY;
	  powVx.vtX *= sVx.vtX;
	}
	powVx.vtY *= sVx.vtY;
      }
    }
    /* Perform singular value decomposition of matrix A. */
    errNum= WlzErrorFromAlg(AlgMatrixSVDecomp(aMx, nPts, nCoef, wMx, vMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    wMax = 0.0;
    for(idN = 0; idN < nCoef; ++idN)
    {
      if(*(wMx + idN) > wMax)
      {
	wMax = *(wMx + idN);
      }
    }
    thresh = tol * wMax;
    for(idN = 0; idN < nCoef; ++idN)
    {
      if(*(wMx + idN) < thresh)
      {
	*(wMx + idN) = 0.0;
      }
    }
    /* Fill matrix b for x coordinate */
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bMx + idM) = (sPts + idM)->vtX - (dPts + idM)->vtX;
    }
    /* Solve for x polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nPts, nCoef, wMx, vMx,
    					        bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the x polynomial coefficients, fill matrix b for
       y coordinate and re-solve. */
    for(idN = 0; idN < nCoef; ++idN)
    {
      (basis->poly + idN)->vtX = *(bMx + idN);
    }
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bMx + idM) = (sPts + idM)->vtY - (dPts + idM)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nPts, nCoef, wMx, vMx,
    						bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the ypolynomial coefficients. */
    for(idN = 0; idN < nCoef; ++idN)
    {
      (basis->poly + idN)->vtY = *(bMx + idN);
    }
  }
  if(bMx)
  {
    AlcFree(bMx);
  }
  if(wMx)
  {
    AlcFree(wMx);
  }
  if(aMx)
  {
    (void )AlcDouble2Free(aMx);
  }
  if(vMx)
  {
    (void )AlcDouble2Free(vMx);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(basis)
    {
      if(basis->poly)
      {
	AlcFree(basis->poly);
      }
      AlcFree(basis);
      }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnConfFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new conformal basis function transform,	*
*		which will transform an object with the given		*
*		source verticies into an object with the given 		*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		int order:		Order of conformal poly.	*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzBasisFnTransform *WlzBasisFnConfFromCPts(int nPts, int order,
						   WlzDVertex2 *dPts,
						   WlzDVertex2 *sPts,
						   WlzErrorNum *dstErr)
{
  int  		idM,
  		idN,
		idX,
  		idY,
		nCoef;
  double	thresh,
  		wMax;
  double	*bMx = NULL,
  		*wMx = NULL;
  double	**aMx = NULL,
  		**vMx = NULL;
  WlzDVertex2	powVx,
  		sVx;
  ComplexD	z, zPow;
  WlzBasisFnTransform *basis = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  if((order < 0) || (nPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((basis = (WlzBasisFnTransform *)
	           AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    basis->type = WLZ_TRANSFORM_2D_BASISFN;
    basis->linkcount = 0;
    basis->freeptr = NULL;
    basis->basisFn = WLZ_BASISFN_CONF_POLY;
    basis->nPoly = order;
    basis->nBasis = 0;
    basis->nVtx = 0;
    basis->delta = 0.0;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nCoef = (order + 1) + (order + 1);
    if(((wMx = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bMx = (double *)AlcMalloc(sizeof(double) * 2 * nPts)) == NULL) ||
       (AlcDouble2Malloc(&vMx, 2*nPts, nCoef) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&aMx, 2*nPts, nCoef) !=  ALC_ER_NONE) ||
       ((basis->poly = (WlzDVertex2 *)
		       AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Fill matrix A. */
    for(idM = 0; idM < nPts; ++idM)
    {
      sVx = *(sPts + idM);
      z.re = sVx.vtX;
      z.im = sVx.vtY;
      zPow.re = 1.0;
      zPow.im = 0.0;
      for(idY = 0, idX = basis->nPoly + 1; idY <= basis->nPoly; ++idY, idX++)
      {
	aMx[idM][idY] = zPow.re;
	aMx[idM][idX] = -zPow.im;
	aMx[idM + nPts][idY] = zPow.im;
	aMx[idM + nPts][idX] = zPow.re;
	zPow = AlgCMult(zPow, z);
      }
   }
    /* Perform singular value decomposition of matrix A. */
    errNum= WlzErrorFromAlg(AlgMatrixSVDecomp(aMx, nPts, nCoef, wMx, vMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    wMax = 0.0;
    for(idN = 0; idN < nCoef; ++idN)
    {
      if(*(wMx + idN) > wMax)
      {
	wMax = *(wMx + idN);
      }
    }
    thresh = tol * wMax;
    for(idN = 0; idN < nCoef; ++idN)
    {
      if(*(wMx + idN) < thresh)
      {
	*(wMx + idN) = 0.0;
      }
    }
    /* Fill matrix b for x coordinate */
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bMx + idM) = (sPts + idM)->vtX - (dPts + idM)->vtX;
      *(bMx + idM + nPts) = (sPts + idM)->vtY - (dPts + idM)->vtY;
    }
    /* Solve for conformal polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nPts, nCoef, wMx, vMx,
    					        bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the conformal polynomial coefficients */
    for(idN = 0; idN < (order + 1); ++idN)
    {
      (basis->poly + idN)->vtX = *(bMx + idN);
      (basis->poly + idN)->vtY = *(bMx + idN + order + 1);
    }
  }
  if(bMx)
  {
    AlcFree(bMx);
  }
  if(wMx)
  {
    AlcFree(wMx);
  }
  if(aMx)
  {
    (void )AlcDouble2Free(aMx);
  }
  if(vMx)
  {
    (void )AlcDouble2Free(vMx);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(basis)
    {
      if(basis->poly)
      {
	AlcFree(basis->poly);
      }
      AlcFree(basis);
      }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnMQFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new multiquadric basis function transform,	*
*		which will transform an object with the given		*
*		source verticies into an object with the given 		*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		double delta:		Normalized delta value in range	*
*					[> 0.0 , < 1.0 ].		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzBasisFnTransform *WlzBasisFnMQFromCPts(int nPts,
						 WlzDVertex2 *dPts,
						 WlzDVertex2 *sPts,
						 double delta,
						 WlzErrorNum *dstErr)
{
  int		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys;
  double	tD0,
		tD1,
		deltaSq,
		range,
		thresh,
		wMax;
  double	*bMx = NULL,
  		*wMx = NULL;
  double	**aMx = NULL,
  		**vMx = NULL;
  WlzDVertex2	tDVx0;
  WlzDBox2	extentDB;
  WlzBasisFnTransform *basis = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  deltaSq = delta * delta;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basis = (WlzBasisFnTransform *)
	       AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL) ||
     ((basis->poly = (WlzDVertex2 *)
		     AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basis->basis = (WlzDVertex2 *)
		      AlcMalloc(sizeof(WlzDVertex2) * nSys)) == NULL) ||
     ((basis->verticies = (WlzDVertex2 *)
			  AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent(&extentDB, dPts, sPts, nPts);
    tD0 = extentDB.xMax - extentDB.xMin;
    tD1 = extentDB.yMax - extentDB.yMin;
    range = (tD0 > tD1)? tD0: tD1;
    if(range <= 1.0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    basis->type = WLZ_TRANSFORM_2D_BASISFN;
    basis->linkcount = 0;
    basis->freeptr = NULL;
    basis->basisFn = WLZ_BASISFN_MQ;
    basis->nPoly = 2;
    basis->nBasis = nPts;
    basis->nVtx = nPts;
    basis->delta = deltaSq * range * range;
    WlzValueCopyDVertexToDVertex(basis->verticies, dPts, nPts);
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aMx + idY) + idX) = 0.0;
      }
      *(bMx + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      *(bMx + idY3) = (sPts + idY)->vtX - tD0;
      tD0 = (tD0 - extentDB.xMin) / range;
      tD1 = ((dPts + idY)->vtY - extentDB.yMin) / range;
      *(*(aMx + idY3) + 0) = 1.0;
      *(*(aMx + idY3) + 1) = tD0;
      *(*(aMx + idY3) + 2) = tD1;
      *(*(aMx + 0) + idY3) = 1.0;
      *(*(aMx + 1) + idY3) = tD0;
      *(*(aMx + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	tDVx0.vtX = ((dPts + idX)->vtX - (dPts + idY)->vtX) / range;
	tDVx0.vtX *= tDVx0.vtX;
	tDVx0.vtY = ((dPts + idX)->vtY - (dPts + idY)->vtY) / range;
	tDVx0.vtY *= tDVx0.vtY;
	tD0 = tDVx0.vtX + tDVx0.vtY;
	tD1 = (tD0 > DBL_EPSILON)? sqrt(tD0 + deltaSq): delta;
	idX3 = idX + 3;
	*(*(aMx + idY3) + idX3) = tD1;
	*(*(aMx + idX3) + idY3) = tD1;
      }
      *(*(aMx + idY3) + idY3) = delta;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aMx, nSys, nSys, wMx, vMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    wMax = 0.0;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) > wMax)
      {
	wMax = *(wMx + idN);
      }
    }
    thresh = tol * wMax;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) < thresh)
      {
	*(wMx + idN) = 0.0;
      }
    }
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
    						vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
       and the y polynomial coefficients. */
    WlzBasisFnMQCoeff(basis, bMx,  &extentDB, range, 1);
    *(bMx + 0) = 0.0;
    *(bMx + 1) = 0.0;
    *(bMx + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      *(bMx + idY3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
    						vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients. */
    WlzBasisFnMQCoeff(basis, bMx,  &extentDB, range, 0);
  }
  if(bMx)
  {
    AlcFree(bMx);
  }
  if(wMx)
  {
    AlcFree(wMx);
  }
  if(aMx)
  {
    (void )AlcDouble2Free(aMx);
  }
  if(vMx)
  {
    (void )AlcDouble2Free(vMx);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(basis)
    {
      if(basis->verticies)
      {
	AlcFree(basis->verticies);
      }
      if(basis->basis)
      {
	AlcFree(basis->basis);
      }
      if(basis->poly)
      {
	AlcFree(basis->poly);
      }
      AlcFree(basis);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnTPSFromCPts					*
* Returns:	WlzBasisFnTransform *:	New basis function transform.	*
* Purpose:	Creates a new thin plate spline basis function		*
*		transform, which will transform an object with the	*
*		given source verticies into an object with the given 	*
*		destination verticies.					*
* Global refs:	-							*
* Parameters:	int nPts:		Number of control point pairs.	*
*		WlzDVertex2 *dPts:	Destination control points.	*
*		WlzDVertex2 *sPts:	Source control points.		*
*		WlzErrorNum *dstErr:	Destination error pointer,	*
*					may be NULL.			*
************************************************************************/
static WlzBasisFnTransform *WlzBasisFnTPSFromCPts(int nPts,
						  WlzDVertex2 *dPts,
						  WlzDVertex2 *sPts,
						  WlzErrorNum *dstErr)
{
  int		idN,
  		idX,
		idY,
		idX3,
		idY3,
		nSys;
  double	tD0,
		tD1,
		range,
		thresh,
		wMax;
  double	*bMx = NULL,
  		*wMx = NULL;
  double	**aMx = NULL,
  		**vMx = NULL;
  WlzBasisFnTransform *basis = NULL;
  WlzDVertex2	tDVx0;
  WlzDBox2	extentDB;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basis = (WlzBasisFnTransform *)
	       AlcCalloc(sizeof(WlzBasisFnTransform), 1)) == NULL) ||
     ((basis->poly = (WlzDVertex2 *)
		     AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basis->basis = (WlzDVertex2 *)
		      AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basis->verticies = (WlzDVertex2 *)
			  AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent(&extentDB, dPts, sPts, nPts);
    tD0 = extentDB.xMax - extentDB.xMin;
    tD1 = extentDB.yMax - extentDB.yMin;
    range = WLZ_MAX(tD0, tD1);
    if(range <= 1.0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    basis->type = WLZ_TRANSFORM_2D_BASISFN;
    basis->linkcount = 0;
    basis->freeptr = NULL;
    basis->basisFn = WLZ_BASISFN_TPS;
    basis->nPoly = 2;
    basis->nBasis = nPts;
    basis->nVtx = nPts;
    basis->delta = 0.0;
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aMx + idY) + idX) = 0.0;
      }
      *(bMx + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      *(bMx + idY3) = (sPts + idY)->vtX - tD0;
      tD0 = (tD0 - extentDB.xMin) / range;
      tD1 = ((dPts + idY)->vtY - extentDB.yMin) / range;
      *(*(aMx + idY3) + 0) = 1.0;
      *(*(aMx + idY3) + 1) = tD0;
      *(*(aMx + idY3) + 2) = tD1;
      *(*(aMx + 0) + idY3) = 1.0;
      *(*(aMx + 1) + idY3) = tD0;
      *(*(aMx + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	tDVx0.vtX = ((dPts + idX)->vtX - (dPts + idY)->vtX) / range;
	tDVx0.vtX *= tDVx0.vtX;
	tDVx0.vtY = ((dPts + idX)->vtY - (dPts + idY)->vtY) / range;
	tDVx0.vtY *= tDVx0.vtY;
	tD0 = tDVx0.vtX + tDVx0.vtY;
	tD1 = (tD0 > DBL_EPSILON)? tD0 * log(tD0): 0.0;
	idX3 = idX + 3;
	*(*(aMx + idY3) + idX3) = tD1;
	*(*(aMx + idX3) + idY3) = tD1;
      }
      *(*(aMx + idY3) + idY3) = 0.0;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aMx, nSys, nSys, wMx, vMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    wMax = 0.0;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) > wMax)
      {
	wMax = *(wMx + idN);
      }
    }
    thresh = tol * wMax;
    for(idN = 0; idN < nSys; ++idN)
    {
      if(*(wMx + idN) < thresh)
      {
	*(wMx + idN) = 0.0;
      }
    }
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
    						vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValueCopyDVertexToDVertex(basis->verticies, dPts, nPts);
    WlzBasisFnTPSCoeff(basis, bMx,  &extentDB, range, 1);
    *(bMx + 0) = 0.0;
    *(bMx + 1) = 0.0;
    *(bMx + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      *(bMx + idY + 3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aMx, nSys, nSys, wMx,
				    		vMx, bMx));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnTPSCoeff(basis, bMx,  &extentDB, range, 0);
  }
  if(bMx)
  {
    AlcFree(bMx);
  }
  if(wMx)
  {
    AlcFree(wMx);
  }
  if(aMx)
  {
    (void )AlcDouble2Free(aMx);
  }
  if(vMx)
  {
    (void )AlcDouble2Free(vMx);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(basis)
    {
      if(basis->verticies)
      {
	AlcFree(basis->verticies);
      }
      if(basis->basis)
      {
	AlcFree(basis->basis);
      }
      if(basis->poly)
      {
	AlcFree(basis->poly);
      }
      AlcFree(basis);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basis);
}

/************************************************************************
* Function:	WlzBasisFnVxExtent					*
* Returns:	void							*
* Purpose:	Computes the extent (bounding box) of two vectors	*
*		of verticies.						*
* Global refs:	-							*
* Parameters:	WlzDBox2 *extentDB:	Pointer for the extent of the	*
*					verticies.			*
*		WlzDVertex2 *vx0:	First vector of verticies.	*
*		WlzDVertex2 *vx1:	Second vector of verticies.	*
*		int nPts:		Number of verticies in each	*
*					vector.				*
************************************************************************/
static void	WlzBasisFnVxExtent(WlzDBox2 *extentDB,
				   WlzDVertex2 *vx0, WlzDVertex2 *vx1,
				   int nPts)
{
  double	tD0,
		tD1,
		tD2;

  extentDB->xMin = extentDB->xMax = vx0->vtX;
  extentDB->yMin = extentDB->yMax = vx0->vtY;
  while(nPts-- > 0)
  {
    if((tD0 = vx0->vtX) > (tD1 = vx1->vtX))
    {
      tD2 = tD0;
      tD0 = tD1;
      tD1 = tD2;
    }
    if(tD0 < extentDB->xMin)
    {
      extentDB->xMin = tD0;
    }
    if(tD1 > extentDB->xMax)
    {
      extentDB->xMax = tD1;
    }
    if((tD0 = vx0->vtY) > (tD1 = vx1->vtY))
    {
      tD2 = tD0;
      tD0 = tD1;
      tD1 = tD2;
    }
    if(tD0 < extentDB->yMin)
    {
      extentDB->yMin = tD0;
    }
    if(tD1 < extentDB->yMax)
    {
      extentDB->yMax = tD1;
    }
    ++vx0;
    ++vx1;
  }
}

/************************************************************************
* Function:	WlzBasisFnGaussCoeff					*
* Returns:	void							*
* Purpose:	Extracts the Gaussian coefficients from the		*
*		given column vector using the given extent and range	*
*		for transformation to pixel space.			*
* Global refs:	-							*
* Parameters:	WlzBasisFnTransform *basis: Allocated basis function	*
*					transform to be filled in.	*
*					verticies.			*
*		double *vec:		Given column vector.		*
*		int forX:		True if the coefficients are	*
*					for the x coordinate.		*
************************************************************************/
static void	WlzBasisFnGaussCoeff(WlzBasisFnTransform *basis,
				     double *vec, int forX)
{
  int		idN;
  double	*vecP;
  WlzDVertex2	*basisVxP,
  		*polyVxP;

  vecP = vec;
  basisVxP = basis->basis;
  polyVxP = basis->poly;
  if(forX)
  {
    polyVxP++->vtX = *vecP++;
    polyVxP++->vtX = *vecP++;
    polyVxP->vtX = *vecP++;
    for(idN = 0; idN < basis->nBasis; ++idN)
    {
      basisVxP++->vtX = *vecP++;
    }
  }
  else
  {
    polyVxP++->vtY = *vecP++;
    polyVxP++->vtY = *vecP++;
    polyVxP->vtY = *vecP++;
    for(idN = 0; idN < basis->nBasis; ++idN)
    {
      basisVxP++->vtY = *vecP++;
    }
  }
}


/************************************************************************
* Function:	WlzBasisFnMQCoeff					*
* Returns:	void							*
* Purpose:	Extracts the multiquadric coefficients from the		*
*		given column vector using the given extent and range	*
*		for transformation to pixel space.			*
* Global refs:	-							*
* Parameters:	WlzBasisFnTransform *basis: Allocated basis function	*
*					transform to be filled in.	*
*					verticies.			*
*		double *vec:		Given column vector.		*
*		WlzDBox2 *extentDB:	Extent of the verticies.	*
*		double range:		Range of the verticies.		*
*		int forX:		True if the coefficients are	*
*					for the x coordinate.		*
************************************************************************/
static void	WlzBasisFnMQCoeff(WlzBasisFnTransform *basis,
				  double *vec, WlzDBox2 *extentDB,
				  double range, int forX)
{
  int		idN;
  double 	vec0,
  		vec1,
  		vec2;
  double	*vecP;
  WlzDVertex2	*basisVxP,
  		*polyVxP;

  vec0 = *vec;
  vec1 = *(vec + 1);
  vec2 = *(vec + 2);
  vecP = vec + 3;
  basisVxP = basis->basis;
  polyVxP = basis->poly;
  if(forX)
  {
    polyVxP++->vtX = vec0 -
    		     (((vec1 * extentDB->xMin) +
    		       (vec2 * extentDB->yMin)) / range);
    polyVxP++->vtX = vec1 / range;
    polyVxP->vtX = vec2 / range;
    for(idN = 0; idN < basis->nBasis; ++idN)
    {
      basisVxP++->vtX = *vecP++ / range;
    }
  }
  else
  {
    polyVxP++->vtY = vec0 -
    		     (((vec1 * extentDB->xMin) +
    		       (vec2 * extentDB->yMin)) / range);
    polyVxP++->vtY = vec1 / range;
    polyVxP->vtY = vec2 / range;
    for(idN = 0; idN < basis->nBasis; ++idN)
    {
      basisVxP++->vtY = *vecP++ / range;
    }
  }
}

/************************************************************************
* Function:	WlzBasisFnTPSCoeff					*
* Returns:	void							*
* Purpose:	Extracts the thin plate spline coefficients from the	*
*		given column vector using the given extent and range	*
*		for transformation to pixel space.			*
* Global refs:	-							*
* Parameters:	WlzBasisFnTransform *basis: Allocated basis function	*
*					transform to be filled in.	*
*					verticies.			*
*		double *vec:		Given column vector.		*
*		WlzDBox2 *extentDB:	Extent of the verticies.	*
*		double range:		Range of the verticies.		*
*		int forX:		True if the coefficients are	*
*					for the x coordinate.		*
************************************************************************/
static void	WlzBasisFnTPSCoeff(WlzBasisFnTransform *basis,
				   double *vec, WlzDBox2 *extentDB,
				   double range, int forX)
{
  int		idN;
  double	tD0,
		tD1,
		rangeSq,
		sumLogCoeffRSq;
  WlzDVertex2	tDVx0;

  rangeSq = range * range;
  tD0 = 2.0 / rangeSq;
  sumLogCoeffRSq = 0.0;
  for(idN = 0; idN < basis->nBasis; ++idN)
  {
    tD1 = *(vec + idN + 3);
    tDVx0 = *(basis->verticies + idN);
    tDVx0.vtX = tDVx0.vtX - extentDB->xMin;
    tDVx0.vtX *= tDVx0.vtX;
    tDVx0.vtY = tDVx0.vtY - extentDB->yMin;
    tDVx0.vtY *= tDVx0.vtY;
    sumLogCoeffRSq += tD1 * (tDVx0.vtX + tDVx0.vtY);
    if(forX)
    {
      (basis->basis + idN)->vtX = tD1 * tD0;
    }
    else
    {
      (basis->basis + idN)->vtY = tD1 * tD0;
    }
  }
  sumLogCoeffRSq /= rangeSq;
  if(forX)
  {
    (basis->poly + 1)->vtX = *(vec + 1) / range;
    (basis->poly + 2)->vtX = *(vec + 2) / range;
    (basis->poly + 0)->vtX = *(vec + 0) -
			     ((basis->poly + 1)->vtX * extentDB->xMin) -
			     ((basis->poly + 2)->vtX * extentDB->yMin) -
			     (log(rangeSq) * sumLogCoeffRSq);
  }
  else
  {
    (basis->poly + 1)->vtY = *(vec + 1) / range;
    (basis->poly + 2)->vtY = *(vec + 2) / range;
    (basis->poly + 0)->vtY = *(vec + 0) -
			      ((basis->poly + 1)->vtY * extentDB->yMin) -
			      ((basis->poly + 2)->vtY * extentDB->xMin) -
			      (log(rangeSq) * sumLogCoeffRSq);
  }
}
