#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBasisFn.c
* \author       Bill Hill
* \date         January 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz functions for which are themselves the sum of
*		radialy symetric functions.
* \ingroup	WlzFunction
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>

static void			WlzBasisFnVxExtent2D(
				  WlzDBox2 *extentDB,
				  WlzDVertex2 *vx0,
				  WlzDVertex2 *vx1,
				  int nPts);
static void			WlzBasisFnGauss2DCoef(
				  WlzBasisFn *basisFn,
				  double *vec,
				  int forX);
static void			WlzBasisFnMQ2DCoexff(
				  WlzBasisFn *basisFn,
				  double *vec,
				  WlzDBox2 *extentDB,
				  double range,
				  int forX);
static void			WlzBasisFnTPS2DCoef(
				  WlzBasisFn *basisFn,
				  double *vec,
				  WlzDBox2 *extentDB,
				  double range,
				  int forX);
static WlzDVertex2 		WlzBasisFnValueRedPoly2D(
				  WlzDVertex2 *poly,
				  WlzDVertex2 srcVx);

/*!
* \return	Woolz error number.
* \ingroup	WlzFunction
* \brief	Free's the given basis function.
* \param	basisFn			Given basis function, may be NULL.
*/
WlzErrorNum	WlzBasisFnFree(WlzBasisFn *basisFn)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basisFn)
  {
    if(basisFn->poly.v)
    {
      AlcFree(basisFn->poly.v);
    }
    if(basisFn->basis.v)
    {
      AlcFree(basisFn->basis.v);
    }
    if(basisFn->vertices.v)
    {
      AlcFree(basisFn->vertices.v);
    }
    if(basisFn->param)
    {
      AlcFree(basisFn->param);
    }
    AlcFree(basisFn);
  }
  return(errNum);
}

/*!
* \return       New vertex value.
* \ingroup      WlzFunction
* \brief        Calculates the value for the given vertex using
*               a polynomial basis function.
* \param        basis                   Basis function.
* \param        srcVx                   Source vertex.
*/
WlzDVertex2 	WlzBasisFnValuePoly2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           idX,
                idY;
  double        tD0;
  WlzDVertex2   *polyP;
  WlzDVertex2   powVx,
                newVx;


  powVx.vtY = 1.0;
  polyP = basisFn->poly.d2;
  newVx.vtY = newVx.vtX = 0.0;
  for(idY = 0; idY <= basisFn->nPoly; ++idY)
  {
    powVx.vtX = 1.0;
    for(idX = 0; idX <= basisFn->nPoly; ++idX)
    {
      tD0 = powVx.vtX * powVx.vtY;
      newVx.vtX += polyP->vtX * tD0;
      newVx.vtY += polyP->vtY * tD0;
      powVx.vtX *= srcVx.vtX;
      ++polyP;
    }
    powVx.vtY *= srcVx.vtY;
  }
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a Gaussian basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueGauss2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1,
		delta;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;

  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  count = basisFn->nVtx;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    tD1 = (tD0 > DBL_EPSILON)? exp(tD0 * delta): 1.0;
    newVx.vtX += basisCo->vtX * tD1;
    newVx.vtY += basisCo->vtY * tD1;
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnValueRedPoly2D(basisFn->poly.d2, srcVx);
  newVx.vtX = newVx.vtX + polyVx.vtX;
  newVx.vtY = newVx.vtY + polyVx.vtY;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueMQ2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1,
		delta;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;

  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  count = basisFn->nVtx;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    if(tD0 > DBL_EPSILON)
    {
      tD0 = sqrt(tD0 + delta);
      newVx.vtX += basisCo->vtX * tD0;
      newVx.vtY += basisCo->vtY * tD0;
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnValueRedPoly2D(basisFn->poly.d2, srcVx);
  newVx.vtX = newVx.vtX + polyVx.vtX;
  newVx.vtY = newVx.vtY + polyVx.vtY;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a thin plate spline basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueTPS2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           count;
  double        tD0,
		tD1;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;

  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  count = basisFn->nVtx;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  while(count-- > 0)
  {
    tD0 = srcVx.vtX - cPts->vtX;
    tD1 = srcVx.vtY - cPts->vtY;
    tD0 = (tD0 * tD0) + (tD1 * tD1);
    if(tD0 > DBL_EPSILON)
    {
      tD0 *= log(tD0);
      newVx.vtX += basisCo->vtX * tD0;
      newVx.vtY += basisCo->vtY * tD0;
    }
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnValueRedPoly2D(basisFn->poly.d2, srcVx);
  newVx.vtX = (newVx.vtX / 2) + polyVx.vtX;
  newVx.vtY = (newVx.vtY / 2) + polyVx.vtY;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a conformal polynomial basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueConf2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int		i;
  ComplexD	z, w, powW, a, b;
  WlzDVertex2	newVx;
  WlzDVertex2	*polyP;

  polyP = basisFn->poly.d2;
  w.re = srcVx.vtX;
  w.im = srcVx.vtY;
  z.re = polyP[0].vtX;
  z.im = polyP[0].vtY;
  powW.re = 1.0;
  powW.im = 0.0;

  for(i=1; i <= basisFn->nPoly; i++){
    powW = AlgCMult(powW, w);
    a.re = polyP[i].vtX;
    a.im = polyP[i].vtY;
    b = AlgCMult(a, powW);
    z.re += b.re;
    z.im += b.im;
  }
  newVx.vtX = z.re;
  newVx.vtY = z.im;
  return(newVx);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new Gaussian basis function.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range 
*					[> 0.0 , < 1.0 ].
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnGauss2DFromCPts(int nPts,
				    WlzDVertex2 *dPts, WlzDVertex2 *sPts,
				    double delta, WlzErrorNum *dstErr)
{
  int		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys;
  double	tD0,
		tD1,
		deltaRg,
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
  WlzBasisFn 	*basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  deltaSq = delta * delta;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL) ||
     ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basisFn->param = AlcMalloc(sizeof(double))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent2D(&extentDB, dPts, sPts, nPts);
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
    basisFn->type = WLZ_FN_BASIS_2DGAUSS;
    basisFn->nPoly = 2;
    basisFn->nBasis = nPts;
    basisFn->nVtx = nPts;
    deltaRg = deltaSq / (range * range);
    *((double *)(basisFn->param)) = deltaRg;
    WlzValueCopyDVertexToDVertex(basisFn->vertices.d2, dPts, nPts);
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
	tD0 = (tDVx0.vtX + tDVx0.vtY) * deltaRg;
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
    WlzBasisFnGauss2DCoef(basisFn, bMx, 1);
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
    WlzBasisFnGauss2DCoef(basisFn, bMx, 0);
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
    (void )WlzBasisFnFree(basisFn);
    basisFn = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new polynomial basis function.
* \param	nPts			Number of control point pairs.
* \param	order			Order of polynomial.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnPoly2DFromCPts(int nPts, int order,
				   WlzDVertex2 *dPts, WlzDVertex2 *sPts,
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
  WlzBasisFn *basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  if((order < 0) || (nPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((basisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    basisFn->type = WLZ_FN_BASIS_2DPOLY;
    basisFn->nPoly = order;
    basisFn->nBasis = 0;
    basisFn->nVtx = 0;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nCoef = (order + 1) * (order + 1);
    if(((wMx = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bMx = (double *)AlcMalloc(sizeof(double) * nPts)) == NULL) ||
       (AlcDouble2Malloc(&vMx, nPts, nCoef) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&aMx, nPts, nCoef) !=  ALC_ER_NONE) ||
       ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
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
      for(idY = 0; idY <= basisFn->nPoly; ++idY)
      {
	powVx.vtX = 1.0;
	for(idX = 0; idX <= basisFn->nPoly; ++idX)
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
      (basisFn->poly.d2 + idN)->vtX = *(bMx + idN);
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
      (basisFn->poly.d2 + idN)->vtY = *(bMx + idN);
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
    (void )WlzBasisFnFree(basisFn);
    basisFn = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new conformal basis function.
* \param	nPts			Number of control point pairs.
* \param	order			Order of conformal poly.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnConf2DFromCPts(int nPts, int order,
				WlzDVertex2 *dPts, WlzDVertex2 *sPts,
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
  WlzDVertex2	sVx;
  ComplexD	z, zPow;
  WlzBasisFn *basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  if((order < 0) || (nPts <= 0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((basisFn = (WlzBasisFn *)
	           AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    basisFn->type = WLZ_FN_BASIS_2DCONF_POLY;
    basisFn->nPoly = order;
    basisFn->nBasis = 0;
    basisFn->nVtx = 0;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nCoef = (order + 1) + (order + 1);
    if(((wMx = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bMx = (double *)AlcMalloc(sizeof(double) * 2 * nPts)) == NULL) ||
       (AlcDouble2Malloc(&vMx, 2*nPts, nCoef) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&aMx, 2*nPts, nCoef) !=  ALC_ER_NONE) ||
       ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
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
      for(idY = 0, idX = basisFn->nPoly + 1; idY <= basisFn->nPoly; ++idY, idX++)
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
      (basisFn->poly.d2 + idN)->vtX = *(bMx + idN);
      (basisFn->poly.d2 + idN)->vtY = *(bMx + idN + order + 1);
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
    (void )WlzBasisFnFree(basisFn);
    basisFn = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new multiquadric basis function.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range
*					[> 0.0 , < 1.0 ].
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnMQ2DFromCPts(int nPts,
				WlzDVertex2 *dPts, WlzDVertex2 *sPts,
				double delta, WlzErrorNum *dstErr)
{
  int		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys;
  double	tD0,
		tD1,
		deltaRg,
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
  WlzBasisFn *basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  deltaSq = delta * delta;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basisFn = (WlzBasisFn *)
	       AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL) ||
     ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) * nSys)) == NULL) ||
     ((basisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basisFn->param = AlcMalloc(sizeof(double))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent2D(&extentDB, dPts, sPts, nPts);
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
    basisFn->type = WLZ_FN_BASIS_2DMQ;
    basisFn->nPoly = 2;
    basisFn->nBasis = nPts;
    basisFn->nVtx = nPts;
    deltaRg = deltaSq * range * range;
    *((double *)(basisFn->param)) = deltaRg;
    WlzValueCopyDVertexToDVertex(basisFn->vertices.d2, dPts, nPts);
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
    WlzBasisFnMQ2DCoexff(basisFn, bMx,  &extentDB, range, 1);
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
    WlzBasisFnMQ2DCoexff(basisFn, bMx,  &extentDB, range, 0);
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
    if(basisFn)
    {
      (void )WlzBasisFnFree(basisFn);
      basisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new thin plate spline basis function.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnTPS2DFromCPts(int nPts,
				  WlzDVertex2 *dPts, WlzDVertex2 *sPts,
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
  WlzBasisFn *basisFn = NULL;
  WlzDVertex2	tDVx0;
  WlzDBox2	extentDB;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	tol = 1.0e-06;

  nSys = nPts + 3;
  if(((wMx = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
     ((bMx = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
     (AlcDouble2Malloc(&vMx, nSys, nSys) !=  ALC_ER_NONE) ||
     (AlcDouble2Malloc(&aMx, nSys, nSys) !=  ALC_ER_NONE) ||
     ((basisFn = (WlzBasisFn *)
	       AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL) ||
     ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
     ((basisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL) ||
     ((basisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) * nPts)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnVxExtent2D(&extentDB, dPts, sPts, nPts);
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
    basisFn->type = WLZ_FN_BASIS_2DTPS;
    basisFn->nPoly = 2;
    basisFn->nBasis = nPts;
    basisFn->nVtx = nPts;
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
    WlzValueCopyDVertexToDVertex(basisFn->vertices.d2, dPts, nPts);
    WlzBasisFnTPS2DCoef(basisFn, bMx,  &extentDB, range, 1);
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
    WlzBasisFnTPS2DCoef(basisFn, bMx,  &extentDB, range, 0);
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
    if(basisFn)
    {
      (void )WlzBasisFnFree(basisFn);
      basisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	<void>
* \ingroup	WlzFunction
* \brief	Computes the extent (bounding box) of two arrays of
*		2D vertices.
* \param	extentDB		Pointer for the extent of the
* 					vertices.
* \param	vx0			First vector of vertices.
* \param	vx1			Second vector of vertices.
* \param	nPts			Number of vertices in each vector.
*/
static void	WlzBasisFnVxExtent2D(WlzDBox2 *extentDB,
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

/*!
* \return	<void>
* \ingroup	WlzFunction
* \brief	Extracts the Gaussian coefficients from the given column
* 		vector using the given extent and range for function.
* \param	basisFn			Allocated basis function to
* 					be filled in.
* \param	vec			Given column vector.
* \param	forX			True if the coefficients are for the x
* 					coordinate.
*/
static void	WlzBasisFnGauss2DCoef(WlzBasisFn *basisFn, double *vec, int forX)
{
  int		idN;
  double	*vecP;
  WlzDVertex2	*basisVxP,
  		*polyVxP;

  vecP = vec;
  basisVxP = basisFn->basis.d2;
  polyVxP = basisFn->poly.d2;
  if(forX)
  {
    polyVxP++->vtX = *vecP++;
    polyVxP++->vtX = *vecP++;
    polyVxP->vtX = *vecP++;
    for(idN = 0; idN < basisFn->nBasis; ++idN)
    {
      basisVxP++->vtX = *vecP++;
    }
  }
  else
  {
    polyVxP++->vtY = *vecP++;
    polyVxP++->vtY = *vecP++;
    polyVxP->vtY = *vecP++;
    for(idN = 0; idN < basisFn->nBasis; ++idN)
    {
      basisVxP++->vtY = *vecP++;
    }
  }
}

/*!
* \return	<void>
* \ingroup	WlzFunction
* \brief	Extracts the multiquadric coefficients from the given
*		column vector using the given extent and range for
*		function.
* \param	basisFn			Allocated basis function
*					to be filled in.
* \param	vec			Given column vector.
* \param	extentDB		Extent of the vertices.
* \param	range			Range of the vertices.
* \param	forX		 	True if the coefficients are for
*					the x coordinate.
*/
static void	WlzBasisFnMQ2DCoexff(WlzBasisFn *basisFn,
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
  basisVxP = basisFn->basis.d2;
  polyVxP = basisFn->poly.d2;
  if(forX)
  {
    polyVxP++->vtX = vec0 -
    		     (((vec1 * extentDB->xMin) +
    		       (vec2 * extentDB->yMin)) / range);
    polyVxP++->vtX = vec1 / range;
    polyVxP->vtX = vec2 / range;
    for(idN = 0; idN < basisFn->nBasis; ++idN)
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
    for(idN = 0; idN < basisFn->nBasis; ++idN)
    {
      basisVxP++->vtY = *vecP++ / range;
    }
  }
}

/*!
* \return	<void>
* \ingroup	WlzFunction
* \brief	Extracts the thin plate spline coefficients from the given
*		column vector using the given extent and range for
*		function.
* \param	basisFn			Allocated basis function to
* 					be filled in.
* \param	vec			Given column vector.
* \param	extentDB		Extent of the vertices.
* \param	range			Range of the vertices.
* \param	forX			True if the coefficients are for the x
* 					coordinate.
*/
static void	WlzBasisFnTPS2DCoef(WlzBasisFn *basisFn,
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
  for(idN = 0; idN < basisFn->nBasis; ++idN)
  {
    tD1 = *(vec + idN + 3);
    tDVx0 = *(basisFn->vertices.d2 + idN);
    tDVx0.vtX = tDVx0.vtX - extentDB->xMin;
    tDVx0.vtX *= tDVx0.vtX;
    tDVx0.vtY = tDVx0.vtY - extentDB->yMin;
    tDVx0.vtY *= tDVx0.vtY;
    sumLogCoeffRSq += tD1 * (tDVx0.vtX + tDVx0.vtY);
    if(forX)
    {
      (basisFn->basis.d2 + idN)->vtX = tD1 * tD0;
    }
    else
    {
      (basisFn->basis.d2 + idN)->vtY = tD1 * tD0;
    }
  }
  sumLogCoeffRSq /= rangeSq;
  if(forX)
  {
    (basisFn->poly.d2 + 1)->vtX = *(vec + 1) / range;
    (basisFn->poly.d2 + 2)->vtX = *(vec + 2) / range;
    (basisFn->poly.d2 + 0)->vtX = *(vec + 0) -
			     ((basisFn->poly.d2 + 1)->vtX * extentDB->xMin) -
			     ((basisFn->poly.d2 + 2)->vtX * extentDB->yMin) -
			     (log(rangeSq) * sumLogCoeffRSq);
  }
  else
  {
    (basisFn->poly.d2 + 1)->vtY = *(vec + 1) / range;
    (basisFn->poly.d2 + 2)->vtY = *(vec + 2) / range;
    (basisFn->poly.d2 + 0)->vtY = *(vec + 0) -
			      ((basisFn->poly.d2 + 1)->vtY * extentDB->yMin) -
			      ((basisFn->poly.d2 + 2)->vtY * extentDB->xMin) -
			      (log(rangeSq) * sumLogCoeffRSq);
  }
}

/*!
* \return	Displacement due to reduced polynomial.
* \ingroup	WlzFunction
* \brief	Computes the value of the reduced polynomial
* 		used by the TPS, MQ and Gauss basis functions.
* \param	poly			Given polynomial coefficients.
* \param	srcVx			Source vertex.
*/
static WlzDVertex2 WlzBasisFnValueRedPoly2D(WlzDVertex2 *poly,
					WlzDVertex2 srcVx)
{
  WlzDVertex2	newVx;

  newVx.vtX = poly->vtX;
  newVx.vtY = poly->vtY;
  ++poly;
  newVx.vtX += poly->vtX * srcVx.vtX;
  newVx.vtY += poly->vtY * srcVx.vtX;
  ++poly;
  newVx.vtX += poly->vtX * srcVx.vtY;
  newVx.vtY += poly->vtY * srcVx.vtY;
  return(newVx);
}
