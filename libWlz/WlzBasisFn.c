#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBasisFn_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzBasisFn.c
* \author       Bill Hill, Jianguo Rao
* \date         January 2003
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
* \brief	Functions for creating and manipulating basis functions.
* \ingroup	WlzFunction
*/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

static void			WlzBasisFnEditSV(
				  int n,
				  double *vV);
static void			WlzBasisFnVxExtent2D(
				  WlzDBox2 *extentDB,
				  WlzDVertex2 *vx0,
				  WlzDVertex2 *vx1,
				  int nPts);
static void			WlzBasisFnVxExtent3D(
				  WlzDBox3 *extentDB,
				  WlzDVertex3 *vx0,
				  WlzDVertex3 *vx1,
				  int nPts);
static void			WlzBasisFnGauss2DCoef(
				  WlzBasisFn *basisFn,
				  double *vec,
				  int forX);
static void			WlzBasisFnMQCoexff2D(
				  WlzBasisFn *basisFn,
				  double *vec,
				  WlzDBox2 *extentDB,
				  double range,
				  int component,
				  int rescale);
static void	                  WlzBasisFnMQCoeff3D(
                                  WlzBasisFn *basisFn,
				  double     *vec, 
				  WlzDBox3   *extentDB,
				  double range, 
				  int component,
				  int rescale);
static void			WlzBasisFnTPS2DCoef(
				  WlzBasisFn *basisFn,
				  double *vec,
				  WlzDBox2 *extentDB,
				  double range,
				  int forX);
static double   		WlzBasisFnValueMOSPhiPC(
				  double r,
				  double v,
				  double w,
				  double delta,
				  double rv,
				  double rw,
				  double norm);
static double			WlzBasisFnScalarMOS3DEvalFn(
				  void *basisFn,
				  double rad);
static double			WlzBasisFnEucDistFn2D(
				  void *bFnP,
				  int idP,
				  WlzVertex pos);
static double			WlzBasisFnEucDistFn3D(
				  void *bFnP,
				  int idP,
				  WlzVertex pos);
static double			WlzBasisFnMapDistFn2D(
				  void *bFnP,
				  int idP,
				  WlzVertex pos);
static double			WlzBasisFnMapDistFn3D(
				  void *bFnP,
				  int idP,
				  WlzVertex pos);
static WlzDVertex2 		WlzBasisFnValueRedPoly2D(
				  WlzDVertex2 *poly,
				  WlzDVertex2 srcVx);
static WlzDVertex3      	WlzBasisFnValueRedPoly3D(
                                  WlzDVertex3 *poly,
				  WlzDVertex3 srcVx);
static WlzHistogramDomain 	*WlzBasisFnScalarMOS3DEvalTb(
				  int nPts,
				  WlzDVertex3 *cPts,
				  double delta,
				  double tau,
				  WlzErrorNum *dstErr);
/*!
* \return	Woolz error number.
* \ingroup	WlzFunction
* \brief	Free's the given basis function.
* \param	basisFn			Given basis function, may be NULL.
*/
WlzErrorNum	WlzBasisFnFree(WlzBasisFn *basisFn)
{
  int		idx;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(basisFn)
  {
    AlcFree(basisFn->poly.v);
    AlcFree(basisFn->basis.v);
    AlcFree(basisFn->vertices.v);
    AlcFree(basisFn->sVertices.v);
    AlcFree(basisFn->param);
    if((basisFn->evalData != NULL))
    {
      (void )WlzFreeHistogramDomain(basisFn->evalData);
    }
    if(basisFn->distMap != NULL)
    {
      for(idx = 0; idx < basisFn->nVtx; ++idx)
      {
	AlcFree(basisFn->distMap[idx]);
      }
      AlcFree(basisFn->distMap);
    }
    AlcFree(basisFn);
  }
  return(errNum);
}

/*!
* \return       New vertex value.
* \ingroup      WlzFunction
* \brief        Calculates the value for the given vertex using
*               a 2D polynomial basis function. This is not a basis
*		function just a polynomial, but it makes sense for
*		it to be incuded in the basis functions since they
*		have a polynomial component.
* \param        basisFn                 Basis function.
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
*		a 2D Gaussian basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueGauss2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		delta;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;
  WlzVertex	sPt;

  sPt.d2 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD0 = (tD0 * tD0) + (tD1 * tD1);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
    }
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
*		a 2D multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueMQ2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		delta;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;
  WlzVertex	sPt;

  sPt.d2 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD0 = (tD0 * tD0) + (tD1 * tD1);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
    }
    tD0 = sqrt(tD0 + delta);
    newVx.vtX += basisCo->vtX * tD0;
    newVx.vtY += basisCo->vtY * tD0;
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
* \brief	Calculates the displacement value for the given vertex using
*		a 3D multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex3 	WlzBasisFnValueMQ3D(WlzBasisFn *basisFn, WlzDVertex3 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		tD2,
		delta;
  WlzDVertex3   *basisCo,
		*cPts;
  WlzDVertex3    polyVx,
  		 newVx;
  WlzVertex	sPt;

  sPt.d3 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  newVx.vtZ = 0.0;
  cPts    = basisFn->vertices.d3;
  basisCo = basisFn->basis.d3;
  delta = *((double *)(basisFn->param));
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD2 = srcVx.vtZ - cPts->vtZ;
      tD0 = (tD0 * tD0) + (tD1 * tD1) + (tD2 * tD2);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
    }
    tD0 = sqrt(tD0 + delta);
    newVx.vtX += basisCo->vtX * tD0;
    newVx.vtY += basisCo->vtY * tD0;
    newVx.vtZ += basisCo->vtZ * tD0;
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnValueRedPoly3D(basisFn->poly.d3, srcVx);
  newVx.vtX = newVx.vtX + polyVx.vtX;
  newVx.vtY = newVx.vtY + polyVx.vtY;
  newVx.vtZ = newVx.vtZ + polyVx.vtZ;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a 2D inverse multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueIMQ2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		delta;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;
  WlzVertex	sPt;

  sPt.d2 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  delta = *((double *)(basisFn->param));
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD0 = (tD0 * tD0) + (tD1 * tD1);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
    }
    tD0 = 1.0 / sqrt(tD0 + delta);
    newVx.vtX += basisCo->vtX * tD0;
    newVx.vtY += basisCo->vtY * tD0;
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
* \brief	Calculates the displacement value for the given vertex using
*		a 3D inverse multiquadric basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex3 	WlzBasisFnValueIMQ3D(WlzBasisFn *basisFn, WlzDVertex3 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		tD2,
		delta;
  WlzDVertex3   *basisCo,
		*cPts;
  WlzDVertex3    polyVx,
  		 newVx;
  WlzVertex	sPt;

  sPt.d3 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  newVx.vtZ = 0.0;
  cPts    = basisFn->vertices.d3;
  basisCo = basisFn->basis.d3;
  delta = *((double *)(basisFn->param));
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD2 = srcVx.vtZ - cPts->vtZ;
      tD0 = (tD0 * tD0) + (tD1 * tD1) + (tD2 * tD2);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
    }
    tD0 = 1.0 / sqrt(tD0 + delta);
    newVx.vtX += basisCo->vtX * tD0;
    newVx.vtY += basisCo->vtY * tD0;
    newVx.vtZ += basisCo->vtZ * tD0;
    ++cPts;
    ++basisCo;
  }
  polyVx = WlzBasisFnValueRedPoly3D(basisFn->poly.d3, srcVx);
  newVx.vtX = newVx.vtX + polyVx.vtX;
  newVx.vtY = newVx.vtY + polyVx.vtY;
  newVx.vtZ = newVx.vtZ + polyVx.vtZ;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a 2D thin plate spline basis function.
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex2 	WlzBasisFnValueTPS2D(WlzBasisFn *basisFn, WlzDVertex2 srcVx)
{
  int           idx;
  double        tD0,
		tD1;
  WlzDVertex2    *basisCo,
		*cPts;
  WlzDVertex2    polyVx,
  		newVx;
  WlzVertex	sPt;

  sPt.d2 = srcVx;
  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  cPts = basisFn->vertices.d2;
  basisCo = basisFn->basis.d2;
  for(idx = 0; idx < basisFn->nVtx; ++idx)
  {
    if(basisFn->distFn == NULL)
    {
      tD0 = srcVx.vtX - cPts->vtX;
      tD1 = srcVx.vtY - cPts->vtY;
      tD0 = (tD0 * tD0) + (tD1 * tD1);
    }
    else
    {
      tD0 = basisFn->distFn(basisFn, idx, sPt);
      tD0 *= tD0;
    }
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
  newVx.vtX = (newVx.vtX * 0.5) + polyVx.vtX;
  newVx.vtY = (newVx.vtY * 0.5) + polyVx.vtY;
  return(newVx);
}

/*!
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a 2D conformal polynomial basis function.
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
  for(i=1; i <= basisFn->nPoly; i++)
  {
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
* \return	New vertex value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a 3D multiorder basis function:
*		\f[
		\phi(r) = \frac{1}{4 \pi \delta^2 r}
			  (1 -
			   \frac{w}{w - v} e^{- \sqrt{v} r} +
			   \frac{v}{w - v} e^{- \sqrt{w} r})
		\f]
		\f[
		v = \frac{1 + \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		w = \frac{1 - \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
WlzDVertex3 	WlzBasisFnValueMOS3D(WlzBasisFn *basisFn, WlzDVertex3 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		v,
		w,
		rad,
		phi,
		delta,
		tau,
		rv,
		rw,
		norm;
  WlzDVertex3   *basisCo,
		*cPts;
  WlzDVertex3   dispVx,
  		polyVx,
  		newVx;

  newVx.vtX = 0.0;
  newVx.vtY = 0.0;
  newVx.vtZ = 0.0;
  cPts = basisFn->vertices.d3;
  basisCo = basisFn->basis.d3;
  delta = *((double *)(basisFn->param) + 0);
  tau = *((double *)(basisFn->param) + 1);
  if(basisFn->evalFn)
  {
    for(idx = 0; idx < basisFn->nVtx; ++idx)
    {
      WLZ_VTX_3_SUB(dispVx, srcVx, *cPts);
      rad = WLZ_VTX_3_LENGTH(dispVx);
      phi = basisFn->evalFn((void *)basisFn, rad);
      newVx.vtX += basisCo->vtX * phi;
      newVx.vtY += basisCo->vtY * phi;
      newVx.vtZ += basisCo->vtZ * phi;
      ++cPts;
      ++basisCo;
    }
  }
  else
  {
    tD0 = 2.0 * tau * delta;
    tD1 = sqrt(1 - (tD0 * tD0));
    tD0 = 2.0 * tau;
    v = (1 + tD1) / tD0;
    w = (1 - tD1) / tD0;
    rv = sqrt(v);
    rw = sqrt(w);
    norm = 1.0 / (4.0 * ALG_M_PI * delta * delta);
    for(idx = 0; idx < basisFn->nVtx; ++idx)
    {
      WLZ_VTX_3_SUB(dispVx, srcVx, *cPts);
      rad = WLZ_VTX_3_LENGTH(dispVx);
      phi = WlzBasisFnValueMOSPhiPC(rad, v, w, delta, rv, rw, norm);
      newVx.vtX += basisCo->vtX * phi;
      newVx.vtY += basisCo->vtY * phi;
      newVx.vtZ += basisCo->vtZ * phi;
      ++cPts;
      ++basisCo;
    }
  }
  polyVx = WlzBasisFnValueRedPoly3D(basisFn->poly.d3, srcVx);
  WLZ_VTX_3_ADD(newVx, newVx, polyVx);
  return(newVx);
}

/*!
* \return	New scalar value.
* \ingroup	WlzFunction
* \brief	Calculates the value for the given vertex using
*		a scalar 3D multiorder basis function:
*		\f[
		\phi(r) = \frac{1}{4 \pi \delta^2 r}
			  (1 -
			   \frac{w}{w - v} e^{- \sqrt{v} r} +
			   \frac{v}{w - v} e^{- \sqrt{w} r})
		\f]
		\f[
		v = \frac{1 + \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		w = \frac{1 - \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
* \param	basisFn			Basis function.
* \param	srcVx			Source vertex.
*/
double 		WlzBasisFnValueScalarMOS3D(WlzBasisFn *basisFn,
					   WlzDVertex3 srcVx)
{
  int           idx;
  double        tD0,
		tD1,
		v,
		w,
		phi,
		delta,
		tau,
		rv,
		rw,
		rad,
		norm,
		value;
  double	*basisCo;
  WlzDVertex3	*cPts;
  WlzDVertex3   dispVx;

  value = 0.0;
  cPts = basisFn->vertices.d3;
  basisCo = (double *)(basisFn->basis.v);
  if(basisFn->evalFn)
  {
    for(idx = 0; idx < basisFn->nVtx; ++idx)
    {
      WLZ_VTX_3_SUB(dispVx, srcVx, *cPts);
      rad = WLZ_VTX_3_LENGTH(dispVx);
      phi = basisFn->evalFn((void *)basisFn, rad);
      value += *basisCo * phi;
      ++basisCo;
      ++cPts;
    }
  }
  else
  {
    delta = *((double *)(basisFn->param) + 0);
    tau = *((double *)(basisFn->param) + 1);
    tD0 = 2.0 * tau * delta;
    tD1 = sqrt(1 - (tD0 * tD0));
    tD0 = 2.0 * tau * tau;
    v = (1 + tD1) / tD0;
    w = (1 - tD1) / tD0;
    rv = sqrt(v);
    rw = sqrt(w);
    norm = 1.0 / (4.0 * ALG_M_PI * delta * delta);
    for(idx = 0; idx < basisFn->nVtx; ++idx)
    {
      WLZ_VTX_3_SUB(dispVx, srcVx, *cPts);
      rad = WLZ_VTX_3_LENGTH(dispVx);
      phi = WlzBasisFnValueMOSPhiPC(rad, v, w, delta, rv, rw, norm);
      value += *basisCo * phi;
      ++basisCo;
      ++cPts;
    }
  }
  value += *(double *)(basisFn->poly.v);
  return(value);
}

/*!
* \return	The value of a single multiorder radial basis function.
* \ingroup	WlzFunction
* \brief	Computes the value of a single multiorder radial basis
*		function:
*		\f[
		\phi(r) = \frac{1}{4 \pi \delta^2 r}
			  (1 -
			   \frac{w}{w - v} e^{- \sqrt{v} r} +
			   \frac{v}{w - v} e^{- \sqrt{w} r})
		\f]
		\f[
		v = \frac{1 + \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		w = \frac{1 - \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
* \param	r			Radial distance, \f$(r > 0)\f$.
* \param	delta			The 1st order smoothness parameter
*					\f$\delta\f$, \f$(\delta > 0)\f$.
* \param	tau			The 3rd order smoothness parameter
*					 \f$\tau\f$, \f$(\tau > 0)\f$.
*/
double          WlzBasisFnValueMOSPhi(double r, double delta, double tau)
{
  double        tD0,
                tD1,
                v,
                w,
                rv,
                rw,
                norm,
                phi;

  tD0 = 2.0 * tau * delta;
  tD1 = sqrt(1.0 - (tD0 * tD0));
  tD0 = 2.0 * tau * tau;
  v = (1 + tD1) / tD0;
  w = (1 - tD1) / tD0;
  rv = sqrt(v);
  rw = sqrt(w);
  norm = 1.0 / (4.0 * ALG_M_PI * delta * delta);
  phi = WlzBasisFnValueMOSPhiPC(r, v, w, delta, rv, rw, norm);
  return(phi);
}

/*!
* \return	Distance from given position to control point.
* \ingroup	WlzFunction
* \brief	Computes the distance from the given position to the control
*		point with the given index using a simple Euclidean distance.
* \param	bFnP			Used to pass the basis function
*					data structure.
* \param	idP			Index of the control point.
* \param	pos			Given position passed using the
*					vertex pointer union but always
*					either 2D or 3D double..
*/
static double	WlzBasisFnEucDistFn2D(void *bFnP, int idP, WlzVertex pos)
{
  double	dist = DBL_MAX;
  WlzDVertex2	tmp;
  WlzBasisFn	*bFn;

  if((bFn = (WlzBasisFn *)bFnP) != NULL)
  {
    tmp = *(bFn->vertices.d2 + idP);
    WLZ_VTX_2_SUB(tmp, tmp, pos.d2);
    dist = WLZ_VTX_2_SQRLEN(tmp);
  }
  return(dist);
}

/*!
* \return	Distance from given position to control point.
* \ingroup	WlzFunction
* \brief	Computes the distance from the given position to the control
*		point with the given index using a simple Euclidean distance.
* \param	bFnP			Used to pass the basis function
*					data structure.
* \param	idP			Index of the control point.
* \param	pos			Given position passed using the
*					vertex pointer union but always
*					either 2D or 3D double..
*/
static double	WlzBasisFnEucDistFn3D(void *bFnP, int idP, WlzVertex pos)
{
  double	dist = DBL_MAX;
  WlzDVertex3	tmp;
  WlzBasisFn	*bFn;

  if((bFn = (WlzBasisFn *)bFnP) != NULL)
  {
    tmp = *(bFn->vertices.d3 + idP);
    WLZ_VTX_3_SUB(tmp, tmp, pos.d3);
    dist = WLZ_VTX_3_SQRLEN(tmp);
  }
  return(dist);
}

/*!
* \return	Distance from given position to control point.
* \ingroup	WlzFunction
* \brief	Computes the distance from the given position to the control
*		point with the given index for which a distance map has
*		been computed.
* \param	bFnP			Used to pass the basis function
*					data structure.
* \param	idP			Index of the control point.
* \param	pos			Given position passed using the
*					vertex pointer union but always
*					either 2D or 3D double..
*/
static double	WlzBasisFnMapDistFn2D(void *bFnP, int idP, WlzVertex pos)
{
  int		idE;
  double	dist = DBL_MAX;
  double	*map;
  WlzBasisFn	*bFn;
  WlzCMeshElm2D *elm;
  WlzCMeshNod2D *nod[3];

  bFn = (WlzBasisFn *)bFnP;
  map = bFn->distMap[idP];
  if((idE = WlzCMeshElmEnclosingPos2D(bFn->mesh.m2, -1,
                                      pos.d2.vtX, pos.d2.vtY, 0, NULL)) >= 0)
  {
    elm = (WlzCMeshElm2D *)AlcVectorItemGet(bFn->mesh.m2->res.elm.vec, idE);
    nod[0] = elm->edu[0].nod;
    nod[1] = elm->edu[1].nod;
    nod[2] = elm->edu[2].nod;
    dist = WlzGeomInterpolateTri2D(nod[0]->pos, nod[1]->pos, nod[2]->pos,
			 map[nod[0]->idx], map[nod[1]->idx], map[nod[2]->idx],
			 pos.d2);
    dist *= dist;
  }
  return(dist);
}

/*!
* \return	Distance from given position to control point.
* \ingroup	WlzFunction
* \brief	Computes the distance from the given position to the control
*		point with the given index for which a distance map has
*		been computed.
* \param	bFnP			Used to pass the basis function
*					data structure.
* \param	idP			Index of the control point.
* \param	pos			Given position passed using the
*					vertex pointer union but always
*					either 2D or 3D double..
*/
static double	WlzBasisFnMapDistFn3D(void *bFnP, int idP, WlzVertex pos)
{
  int		idE;
  double	dist = DBL_MAX;
  double	*map;
  WlzBasisFn	*bFn;
  WlzCMeshElm3D *elm;
  WlzCMeshNod3D *nod[4];

  bFn = (WlzBasisFn *)bFnP;
  map = bFn->distMap[idP];
  if((idE = WlzCMeshElmEnclosingPos3D(bFn->mesh.m3, -1,
                                      pos.d3.vtX, pos.d3.vtY, pos.d3.vtZ,
				      0, NULL)) >= 0)
  {
    elm = (WlzCMeshElm3D *)AlcVectorItemGet(bFn->mesh.m3->res.elm.vec, idE);
    nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
    nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
    nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
    nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
    dist = WlzGeomInterpolateTet3D(nod[0]->pos, nod[1]->pos,
				   nod[2]->pos, nod[3]->pos,
				   map[nod[0]->idx], map[nod[1]->idx],
				   map[nod[2]->idx], map[nod[3]->idx],
				   pos.d3);
    dist *= dist;
  }
  return(dist);
}

/*!
* \return	The value of a single multiorder radial basis function.
* \ingroup	WlzFunction
* \brief	Computes the value of a single multiorder radial basis
*		function using either an approximation:
*		\f[
		phi'(r) = \frac{\frac{\sqrt{w} - \sqrt{v}}{w - v}
		                (\sqrt{w v} - \frac{w v}{6}) +
		                \frac{v w r^3}{24}}
		               {4 \pi \delta^2}
*		\f]
*		if the radial distance is near to zero or else the
*		function:
*		\f[
		\phi(r) = \frac{1}{4 \pi \delta^2 r}
			  (1 -
			   \frac{w}{w - v} e^{- \sqrt{v} r} +
			   \frac{v}{w - v} e^{- \sqrt{w} r})
		\f]
*		where
		\f[
		v = \frac{1 + \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		w = \frac{1 - \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
* \param	r			Radial distance, \f$(r > 0)\f$.
* \param	v			Precomputed parameter \f$ v \f$.
* \param	w			Precomputed parameter \f$ w \f$.
* \param	delta			The 1st order smoothness parameter
*					\f$ \delta \f$.
* \param	rv			Precomputed parameter \f$ \sqrt{v} \f$.
* \param	rw			Precomputed parameter \f$ \sqrt{w} \f$.
* \param	norm			Precomputed parameter,
*					\f$ 4 \pi \delta^2 \f$.
*/
static double   WlzBasisFnValueMOSPhiPC(double r, double v, double w,
                                        double delta, double rv, double rw,
                                        double norm)
{
  double        phi = 0.0,
                wv,
                rwv,
                rvr,
                rwr;

  wv = w - v;
  if(fabs(wv) > DBL_EPSILON)
  {
    rwv = rw - rv;
    if(((rvr = rv * r) < 0.01) || ((rwr = rw * r) < 0.01))
    {
      phi = (rwv * rw * rv) / wv;
      if(r > DBL_EPSILON)
      {
        phi += -(v * w * r * r) * ((rwv / (6.0 * wv)) + (r / 24.0));
      }
    }
    else
    {
      phi = (1.0  - (w * exp(-rvr) - v * exp(-rwr)) / wv) / r;
    }
    phi *= norm;
  }
  return(phi);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new Gaussian basis function.
*
* 		This function may be used either to compute a new Gaussian
* 		radial basis function, or if given a previous Gaussian
* 		radial basis function, to recompute a Gaussian basis
* 		recycling the buffers (and distance map if used) from
* 		the prevoius basis function when the control points
* 		have changed. If recycling the buffers a new basis
* 		function is created which uses the (reallocated)
* 		buffers of the previous basis function. After this function
* 		has been called the previous basis function should be freed.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range 
*					[> 0.0 , < 1.0 ].
* \param	prvBasisFn		Previous basis function to be recycled,
* 					may be NULL. Arrays will be used for
* 					the new basis function. Must be freed
* 					after recycling.
* \param	mesh			Used to compute constrained distances.
* 					If NULL Euclidean distances will be
* 					used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnGauss2DFromCPts(int nPts, WlzDVertex2 *dPts,
                                      WlzDVertex2 *sPts, double delta,
				      WlzBasisFn *prvBasisFn,
				      WlzCMesh2D *mesh,
				      WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys,
		maxNod,
		newMaxVx;
  double	tD0,
		tD1,
		deltaRg,
		deltaSq,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDBox2	extentDB;
  WlzBasisFn	*newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	stepVx = 10;

  nSys = nPts + 3;
  deltaSq = delta * delta;
  aM.core = NULL;
  vM.core = NULL;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_2D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m2)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_2DGAUSS)
      {
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
	newMaxVx = prvBasisFn->maxVx;
	if(prvBasisFn->maxVx < nPts)
	{
	  newMaxVx = nPts + stepVx;
	  if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
	                                            1)) == NULL) ||
	     ((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->distMap != NULL) &&
	      ((prvBasisFn->distMap = (double **)
	                              AlcRealloc(prvBasisFn->distMap,
				          sizeof(double *) *
					  newMaxVx)) == NULL)))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *newBasisFn = *prvBasisFn;
	  prvBasisFn->poly.v = NULL;
	  prvBasisFn->basis.v = NULL;
	  prvBasisFn->vertices.v = NULL;
	  prvBasisFn->sVertices.v = NULL;
	  prvBasisFn->param = NULL;
	  prvBasisFn->mesh.v = NULL;
	  prvBasisFn->distMap = NULL;
	  if(newBasisFn->distMap != NULL)
	  {
	    tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
	    for(idN = 0; idN < tI0; ++idN)
	    {
	      if((WlzGeomCmpVtx2D(dPts[idN],
	                          newBasisFn->vertices.d2[idN],
				  DBL_EPSILON) == 0) ||
		  (WlzGeomCmpVtx2D(sPts[idN],
				   newBasisFn->sVertices.d2[idN],
				   DBL_EPSILON) == 0))
	      {
		AlcFree(newBasisFn->distMap[idN]);
		newBasisFn->distMap[idN] = NULL;
	      }
	    }
	    for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
	    {
	      AlcFree(newBasisFn->distMap[idN]);
	      newBasisFn->distMap[idN] = NULL;
	    }
	    for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
	    {
	      newBasisFn->distMap[idN] = NULL;
	    }
	  }
	}
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                                       1)) == NULL) ||
	 ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
	 ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) *
	                                   newMaxVx)) == NULL) ||
	 ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	newBasisFn->nVtx = nPts;
        newBasisFn->maxVx = newMaxVx;
	if(newBasisFn->mesh.m2 == NULL)
	{
          newBasisFn->distFn = WlzBasisFnEucDistFn2D;
	}
	else
	{
	  if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
							 newMaxVx)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn2D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m2->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn2D;
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
	if(newBasisFn->distMap[idN] == NULL)
	{
	  if((newBasisFn->distMap[idN] = (double *)
					 AlcMalloc(sizeof(double) *
						   maxNod)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    errNum = WlzCMeshFMarNodes2D(newBasisFn->mesh.m2,
					 newBasisFn->distMap[idN],
					 1, dPts + idN);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    /* Compute range of displacements. */
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
    newBasisFn->nVtx = nPts;
    newBasisFn->maxVx = newMaxVx;
    newBasisFn->type = WLZ_FN_BASIS_2DGAUSS;
    newBasisFn->nVtx = nPts;
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    WlzValueCopyDVertexToDVertex(newBasisFn->vertices.d2, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex(newBasisFn->sVertices.d2, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    deltaRg = deltaSq / (range * range);
    *((double *)(newBasisFn->param)) = deltaRg;
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      tD1 = (dPts + idY)->vtY;
      *(bV + idY3) = (sPts + idY)->vtX - tD0;
      *(*(aA + idY3) + 0) = 1.0;
      *(*(aA + idY3) + 1) = tD0;
      *(*(aA + idY3) + 2) = tD1;
      *(*(aA + 0) + idY3) = 1.0;
      *(*(aA + 1) + idY3) = tD0;
      *(*(aA + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	sPt.d2 = dPts[idX];
	tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	tD0 *= deltaRg;
	tD1 = (tD0 > DBL_EPSILON)? exp(tD0): 1.0;
	idX3 = idX + 3;
	*(*(aA + idY3) + idX3) = tD1;
	*(*(aA + idX3) + idY3) = tD1;
      }
      *(*(aA + idY3) + idY3) = 1.0;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
       and the y polynomial coefficients. */
    WlzBasisFnGauss2DCoef(newBasisFn, bV, 1);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      *(bV + idY3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients. */
    WlzBasisFnGauss2DCoef(newBasisFn, bV, 0);
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzBasisFnFree(newBasisFn);
    newBasisFn = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
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
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzDVertex2	powVx,
  		sVx;
  WlzBasisFn *basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  aM.core = NULL;
  vM.core = NULL;
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
    if(((wV = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nPts)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nPts, nCoef, NULL)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew( nPts, nCoef, NULL)) == NULL) ||
       ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
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
	  *(*(aA + idM) + idN++) = powVx.vtX * powVx.vtY;
	  powVx.vtX *= sVx.vtX;
	}
	powVx.vtY *= sVx.vtY;
      }
    }
    /* Perform singular value decomposition of matrix A. */
    errNum= WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nCoef, wV);
    /* Fill matrix b for x coordinate */
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bV + idM) = (sPts + idM)->vtX - (dPts + idM)->vtX;
    }
    /* Solve for x polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the x polynomial coefficients, fill matrix b for
       y coordinate and re-solve. */
    for(idN = 0; idN < nCoef; ++idN)
    {
      (basisFn->poly.d2 + idN)->vtX = *(bV + idN);
    }
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bV + idM) = (sPts + idM)->vtY - (dPts + idM)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the ypolynomial coefficients. */
    for(idN = 0; idN < nCoef; ++idN)
    {
      (basisFn->poly.d2 + idN)->vtY = *(bV + idN);
    }
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
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
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzDVertex2	sVx;
  ComplexD	z,
  		zPow;
  WlzBasisFn *basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  aM.core = vM.core = NULL;
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
    if(((wV = (double *)AlcCalloc(sizeof(double), nCoef)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * 2 * nPts)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(2 * nPts, nCoef, NULL)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(2 * nPts, nCoef, NULL)) == NULL) ||
       ((basisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * nCoef)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    /* Fill matrix A. */
    for(idM = 0; idM < nPts; ++idM)
    {
      sVx = *(sPts + idM);
      z.re = sVx.vtX;
      z.im = sVx.vtY;
      zPow.re = 1.0;
      zPow.im = 0.0;
      for(idY = 0, idX = basisFn->nPoly + 1; idY <= basisFn->nPoly;
          ++idY, idX++)
      {
	aA[idM][idY] = zPow.re;
	aA[idM][idX] = -zPow.im;
	aA[idM + nPts][idY] = zPow.im;
	aA[idM + nPts][idX] = zPow.re;
	zPow = AlgCMult(zPow, z);
      }
   }
    /* Perform singular value decomposition of matrix A. */
    errNum= WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nCoef, wV);
    /* Fill matrix b for x coordinate */
    for(idM = 0; idM < nPts; ++idM)
    {
      *(bV + idM) = (sPts + idM)->vtX - (dPts + idM)->vtX;
      *(bV + idM + nPts) = (sPts + idM)->vtY - (dPts + idM)->vtY;
    }
    /* Solve for conformal polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Copy out the conformal polynomial coefficients */
    for(idN = 0; idN < (order + 1); ++idN)
    {
      (basisFn->poly.d2 + idN)->vtX = *(bV + idN);
      (basisFn->poly.d2 + idN)->vtY = *(bV + idN + order + 1);
    }
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
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
*
*		To improve the design matix condition number the problem
*		is rescaled when a distance map is not being used.
*
* 		This function may be used either to compute a new MQ
* 		radial basis function, or if given a previous MQ
* 		radial basis function, to recompute a MQ basis
* 		recycling the buffers (and distance map if used) from
* 		the prevoius basis function when the control points
* 		have changed. If recycling the buffers a new basis
* 		function is created which uses the (reallocated)
* 		buffers of the previous basis function. After this function
* 		has been called the previous basis function should be freed.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range
*					[> 0.0 , < 1.0 ].
* \param	prvBasisFn		Previous basis function to be recycled,
* 					may be NULL. Arrays will be used for
* 					the new basis function. Must be freed
* 					after recycling.
* \param	mesh			Used to compute constrained distances.
* 					If NULL Euclidean distances will be
* 					used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnMQ2DFromCPts(int nPts, WlzDVertex2 *dPts,
                                      WlzDVertex2 *sPts, double delta,
				      WlzBasisFn *prvBasisFn,
				      WlzCMesh2D *mesh,
				      WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys,
		maxNod,
		newMaxVx;
  double	tD0,
		tD1,
		deltaRg,
		deltaSq,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDVertex2	tV0;
  WlzDBox2	extentDB;
  WlzBasisFn	*newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	stepVx = 10;

  aM.core = NULL;
  vM.core = NULL;
  nSys = nPts + 3;
  deltaSq = delta * delta;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_2D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m2)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_2DMQ)
      {
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
	newMaxVx = prvBasisFn->maxVx;
	if((newBasisFn == NULL) &&
           ((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
	                                          1)) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if((errNum == WLZ_ERR_NONE) && (prvBasisFn->maxVx < nPts))
	{
	  newMaxVx = nPts + stepVx;
	  if(((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->distMap != NULL) &&
	      ((prvBasisFn->distMap = (double **)
	                              AlcRealloc(prvBasisFn->distMap,
				          sizeof(double *) *
					  newMaxVx)) == NULL)))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *newBasisFn = *prvBasisFn;
	  prvBasisFn->poly.v = NULL;
	  prvBasisFn->basis.v = NULL;
	  prvBasisFn->vertices.v = NULL;
	  prvBasisFn->sVertices.v = NULL;
	  prvBasisFn->param = NULL;
	  prvBasisFn->mesh.v = NULL;
	  prvBasisFn->distMap = NULL;
	  newBasisFn->maxVx = newMaxVx;
	  newBasisFn->nVtx = nPts;
	  if(newBasisFn->distMap != NULL)
	  {
	    tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
	    for(idN = 0; idN < tI0; ++idN)
	    {
	      if((WlzGeomCmpVtx2D(dPts[idN],
	                          newBasisFn->vertices.d2[idN],
				  DBL_EPSILON) != 0) ||
		  (WlzGeomCmpVtx2D(sPts[idN],
				   newBasisFn->sVertices.d2[idN],
				   DBL_EPSILON) != 0))
	      {
		AlcFree(newBasisFn->distMap[idN]);
		newBasisFn->distMap[idN] = NULL;
	      }
	    }
	    for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
	    {
	      AlcFree(newBasisFn->distMap[idN]);
	      newBasisFn->distMap[idN] = NULL;
	    }
	    for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
	    {
	      newBasisFn->distMap[idN] = NULL;
	    }
	  }
	}
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                                       1)) == NULL) ||
	 ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
	 ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) *
	                                   newMaxVx)) == NULL) ||
	 ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	newBasisFn->nVtx = nPts;
        newBasisFn->maxVx = newMaxVx;
	newBasisFn->mesh.m2 = mesh;
	if(newBasisFn->mesh.m2 != NULL)
	{
	  if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
							 newMaxVx)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn2D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m2->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn2D;
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
	if(newBasisFn->distMap[idN] == NULL)
	{
	  if((newBasisFn->distMap[idN] = (double *)
					 AlcMalloc(sizeof(double) *
						   maxNod)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    errNum = WlzCMeshFMarNodes2D(newBasisFn->mesh.m2,
					 newBasisFn->distMap[idN],
					 1, dPts + idN);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function the design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
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
    newBasisFn->type = WLZ_FN_BASIS_2DMQ;
    newBasisFn->nVtx = nPts;
    WlzValueCopyDVertexToDVertex(newBasisFn->vertices.d2, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex(newBasisFn->sVertices.d2, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    deltaRg = deltaSq * range * range;
    *((double *)(newBasisFn->param)) = deltaRg;
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      *(bV + idY3) = (sPts + idY)->vtX - tD0;
      if(newBasisFn->distMap != NULL)
      {
	tD1 = (dPts + idY)->vtY;
      }
      else
      {
	tD0 = (tD0 - extentDB.xMin) / range;
	tD1 = ((dPts + idY)->vtY - extentDB.yMin) / range;
      }
      *(*(aA + idY3) + 0) = 1.0;
      *(*(aA + idY3) + 1) = tD0;
      *(*(aA + idY3) + 2) = tD1;
      *(*(aA + 0) + idY3) = 1.0;
      *(*(aA + 1) + idY3) = tD0;
      *(*(aA + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	if(newBasisFn->distMap)
	{
	  sPt.d2 = dPts[idX];
	  tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	}
	else
	{
	  tV0.vtX = ((dPts + idX)->vtX - (dPts + idY)->vtX) / range;
	  tV0.vtX *= tV0.vtX;
	  tV0.vtY = ((dPts + idX)->vtY - (dPts + idY)->vtY) / range;
	  tV0.vtY *= tV0.vtY;
	  tD0 = tV0.vtX + tV0.vtY;
	}
	tD1 = (tD0 > DBL_EPSILON)? sqrt(tD0 + deltaSq): delta;
	idX3 = idX + 3;
	*(*(aA + idY3) + idX3) = tD1;
	*(*(aA + idX3) + idY3) = tD1;
      }
      *(*(aA + idY3) + idY3) = delta;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
       and the y polynomial coefficients. */
    WlzBasisFnMQCoexff2D(newBasisFn, bV,  &extentDB, range,
    			 0, (newBasisFn->distFn)? 0: 1);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      *(bV + idY3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients. */
    WlzBasisFnMQCoexff2D(newBasisFn, bV,  &extentDB, range,
    			 1, (newBasisFn->distFn)? 0: 1);
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    if(newBasisFn)
    {
      (void )WlzBasisFnFree(newBasisFn);
      newBasisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
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
* \param	prvBasisFn		Previous basis function to be recycled,
* 					may be NULL. Arrays will be used for
* 					the new basis function. Must be freed
* 					after recycling.
* \param	mesh			Used to compute constrained distances.
* 					If NULL Euclidean distances will be
* 					used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnMQ3DFromCPts(int nPts, WlzDVertex3 *dPts, 
				WlzDVertex3 *sPts, double delta, 
				WlzBasisFn *prvBasisFn,
				WlzCMesh3D *mesh,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
		idX,
  		idY,
                idX4,
		idY4,
		nSys,
		maxNod,
		newMaxVx;
  double	tD0,
		tD1,
                tD2,
		deltaRg,
		deltaSq,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDVertex3	tDVx0;
  WlzDBox3	extentDB;
  WlzBasisFn    *newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	stepVx = 10;

  aM.core = NULL;
  vM.core = NULL;
  nSys = nPts + 4;
  deltaSq = delta * delta;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_3D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m3)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_3DMQ)
      {
        errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
        newMaxVx = prvBasisFn->maxVx;
	if((newBasisFn == NULL) &&
	   ((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
	                                          1)) == NULL))
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if((errNum == WLZ_ERR_NONE) && (prvBasisFn->maxVx < nPts))
	{
	  newMaxVx = nPts + stepVx;
	  if(((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->distMap != NULL) &&
	      ((prvBasisFn->distMap = (double **)
				      AlcRealloc(prvBasisFn->distMap,
					  sizeof(double *) *
					  newMaxVx)) == NULL)))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *newBasisFn = *prvBasisFn;
	  prvBasisFn->poly.v = NULL;
	  prvBasisFn->basis.v = NULL;
	  prvBasisFn->vertices.v = NULL;
	  prvBasisFn->sVertices.v = NULL;
	  prvBasisFn->param = NULL;
	  prvBasisFn->mesh.v = NULL;
	  prvBasisFn->distMap = NULL;
	  newBasisFn->maxVx = newMaxVx;
	  newBasisFn->nVtx = nPts;
	  if(newBasisFn->distMap != NULL)
	  {
	    tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
	    for(idN = 0; idN < tI0; ++idN)
	    {
	      if((WlzGeomCmpVtx3D(dPts[idN],
	                          newBasisFn->vertices.d3[idN],
				  DBL_EPSILON) != 0) ||
	         (WlzGeomCmpVtx3D(sPts[idN],
		                  newBasisFn->sVertices.d3[idN],
				  DBL_EPSILON) != 0))
	      {
	        AlcFree(newBasisFn->distMap[idN]);
		newBasisFn->distMap[idN] = NULL;
	      }
	    }
	    for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
	    {
	      AlcFree(newBasisFn->distMap[idN]);
	      newBasisFn->distMap[idN] = NULL;
	    }
	    for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
	    {
	      newBasisFn->distMap[idN] = NULL;
	    }
	  }
	}
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                             1)) == NULL) ||
	 ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex3) * 4)) == NULL) ||
	 ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex3) *
	                                   newMaxVx)) == NULL) ||
	 ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex3) *
	                                      newMaxVx)) == NULL) ||
	 ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex3) *
					       newMaxVx)) == NULL) ||
	 ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        newBasisFn->nVtx = nPts;
	newBasisFn->maxVx = newMaxVx;
	newBasisFn->mesh.m3 = mesh;
	if(newBasisFn->mesh.m3 != NULL)
	{
	  if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
	                                                 newMaxVx)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn3D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m3->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn3D;
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
        if(newBasisFn->distMap[idN] == NULL)
	{
	  if((newBasisFn->distMap[idN] = (double *)
	                                 AlcMalloc(sizeof(double) *
					           maxNod)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    errNum = WlzCMeshFMarNodes3D(newBasisFn->mesh.m3,
	                                 newBasisFn->distMap[idN],
					 1, dPts + idN);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function the design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    WlzBasisFnVxExtent3D(&extentDB, dPts, sPts, nPts);
    tD0 = extentDB.xMax - extentDB.xMin;
    tD1 = extentDB.yMax - extentDB.yMin;
    tD2 = extentDB.zMax - extentDB.zMin;
    range = (tD0 > tD1)? tD0: tD1;
    if(tD2 > range)
    {
      range = tD2;
    }
    if(range <= 1.0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->type = WLZ_FN_BASIS_3DMQ;
    newBasisFn->nVtx = nPts;
    WlzValueCopyDVertexToDVertex3(newBasisFn->vertices.d3, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex3(newBasisFn->sVertices.d3, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    deltaRg = deltaSq * range * range;
    *((double *)(newBasisFn->param)) = deltaRg;
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 4; ++idY)
    {
      for(idX = 0; idX < 4; ++idX)
      {
        *(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      tD0 = (dPts + idY)->vtX;
      *(bV + idY4) = (sPts + idY)->vtX - tD0;
      if(newBasisFn->distMap != NULL)
      {
        tD1 = (dPts + idY)->vtY;
	tD2 = (dPts + idY)->vtZ;
      }
      else
      {
	tD0 = (tD0 - extentDB.xMin) / range;
	tD1 = ((dPts + idY)->vtY - extentDB.yMin) / range;
	tD2 = ((dPts + idY)->vtZ - extentDB.zMin) / range;
      }
      *(*(aA + idY4) + 0) = 1.0;
      *(*(aA + idY4) + 1) = tD0;
      *(*(aA + idY4) + 2) = tD1;
      *(*(aA + idY4) + 3) = tD2;
      *(*(aA + 0) + idY4) = 1.0;
      *(*(aA + 1) + idY4) = tD0;
      *(*(aA + 2) + idY4) = tD1;
      *(*(aA + 3) + idY4) = tD2;
      for(idX = 0; idX < idY; ++idX)
      {
	if(newBasisFn->distMap)
	{
	  sPt.d3 = dPts[idX];
	  tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	}
	else
	{
	  tDVx0.vtX = ((dPts + idX)->vtX - (dPts + idY)->vtX) / range;
	  tDVx0.vtX *= tDVx0.vtX;
	  tDVx0.vtY = ((dPts + idX)->vtY - (dPts + idY)->vtY) / range;
	  tDVx0.vtY *= tDVx0.vtY;
	  tDVx0.vtZ = ((dPts + idX)->vtZ - (dPts + idY)->vtZ) / range;
	  tDVx0.vtZ *= tDVx0.vtZ;
	  tD0 = tDVx0.vtX + tDVx0.vtY + tDVx0.vtZ;
	}
        tD1 = (tD0 > DBL_EPSILON)? sqrt(tD0 + deltaSq): delta;
	idX4 = idX + 4;
	*(*(aA + idY4) + idX4) = tD1;
	*(*(aA + idX4) + idY4) = tD1;
      }
      *(*(aA + idY4) + idY4) = delta;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
     * and the y polynomial coefficients. */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range,
    			0, (newBasisFn->distFn)? 0: 1);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    *(bV + 3) = 0.0;
     for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      *(bV + idY4) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }

  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients, then set up for nu
     * and the z polynomial coefficients */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range,
    			1, (newBasisFn->distFn)? 0: 1);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    *(bV + 3) = 0.0;
     for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      *(bV + idY4) = (sPts + idY)->vtZ - (dPts + idY)->vtZ;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover nu and the z polynomial coefficients. */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range,
    			2, (newBasisFn->distFn)? 0: 1);
  }
  AlcFree(bV);
  AlcFree(wV);
  (void )AlgMatrixFree(aM);
  (void )AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    if(newBasisFn)
    {
      (void )WlzBasisFnFree(newBasisFn);
      newBasisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new inverse-multiquadric basis function.
*
*		To improve the design matix condition number the problem
*		is rescaled when a distance map is not being used.
*
* 		This function may be used either to compute a new IMQ
* 		radial basis function, or if given a previous IMQ
* 		radial basis function, to recompute a IMQ basis
* 		recycling the buffers (and distance map if used) from
* 		the prevoius basis function when the control points
* 		have changed. If recycling the buffers a new basis
* 		function is created which uses the (reallocated)
* 		buffers of the previous basis function. After this function
* 		has been called the previous basis function should be freed.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range
*					[> 0.0 , < 1.0 ].
* \param	prvBasisFn		Previous basis function to be recycled,
* 					may be NULL. Arrays will be used for
* 					the new basis function. Must be freed
* 					after recycling.
* \param	mesh			Used to compute constrained distances.
* 					If NULL Euclidean distances will be
* 					used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnIMQ2DFromCPts(int nPts, WlzDVertex2 *dPts,
                                      WlzDVertex2 *sPts, double delta,
				      WlzBasisFn *prvBasisFn,
				      WlzCMesh2D *mesh,
				      WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
  		idX,
  		idY,
		idX3,
		idY3,
		nSys,
		maxNod,
		newMaxVx;
  double	tD0,
		tD1,
		tD2,
		deltaSq,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDVertex2	tV0;
  WlzDBox2	extentDB;
  WlzBasisFn	*newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	stepVx = 10;

#ifdef WLZ_BASISFN_DELTA_ENV
  double	val;
  char		*envStr;

  if(((envStr = getenv("WLZ_BASISFN_DELTA")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    delta = val;
  }
#endif /* WLZ_BASISFN_DELTA_ENV */
  aM.core = NULL;
  vM.core = NULL;
  nSys = nPts + 3;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_2D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m2)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_2DIMQ)
      {
	errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
	newMaxVx = prvBasisFn->maxVx;
	if((newBasisFn == NULL) &&
           ((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
	                                          1)) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if((errNum == WLZ_ERR_NONE) && (prvBasisFn->maxVx < nPts))
	{
	  newMaxVx = nPts + stepVx;
	  if(((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
				       sizeof(WlzDVertex2) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->distMap != NULL) &&
	      ((prvBasisFn->distMap = (double **)
	                              AlcRealloc(prvBasisFn->distMap,
				          sizeof(double *) *
					  newMaxVx)) == NULL)))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *newBasisFn = *prvBasisFn;
	  prvBasisFn->poly.v = NULL;
	  prvBasisFn->basis.v = NULL;
	  prvBasisFn->vertices.v = NULL;
	  prvBasisFn->sVertices.v = NULL;
	  prvBasisFn->param = NULL;
	  prvBasisFn->mesh.v = NULL;
	  prvBasisFn->distMap = NULL;
	  newBasisFn->maxVx = newMaxVx;
	  newBasisFn->nVtx = nPts;
	  if(newBasisFn->distMap != NULL)
	  {
	    tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
	    for(idN = 0; idN < tI0; ++idN)
	    {
	      if((WlzGeomCmpVtx2D(dPts[idN],
	                          newBasisFn->vertices.d2[idN],
				  DBL_EPSILON) != 0) ||
		  (WlzGeomCmpVtx2D(sPts[idN],
				   newBasisFn->sVertices.d2[idN],
				   DBL_EPSILON) != 0))
	      {
		AlcFree(newBasisFn->distMap[idN]);
		newBasisFn->distMap[idN] = NULL;
	      }
	    }
	    for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
	    {
	      AlcFree(newBasisFn->distMap[idN]);
	      newBasisFn->distMap[idN] = NULL;
	    }
	    for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
	    {
	      newBasisFn->distMap[idN] = NULL;
	    }
	  }
	}
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                                       1)) == NULL) ||
	 ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
	 ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) *
	                                   newMaxVx)) == NULL) ||
	 ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex2) *
					      newMaxVx)) == NULL) ||
	 ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	newBasisFn->nVtx = nPts;
        newBasisFn->maxVx = newMaxVx;
	newBasisFn->mesh.m2 = mesh;
	if(newBasisFn->mesh.m2 != NULL)
	{
	  if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
							 newMaxVx)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn2D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m2->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn2D;
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
	if(newBasisFn->distMap[idN] == NULL)
	{
	  if((newBasisFn->distMap[idN] = (double *)
					 AlcMalloc(sizeof(double) *
						   maxNod)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    errNum = WlzCMeshFMarNodes2D(newBasisFn->mesh.m2,
					 newBasisFn->distMap[idN],
					 1, dPts + idN);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function the design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
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
    delta *= range;
    deltaSq = delta * delta;
    newBasisFn->type = WLZ_FN_BASIS_2DIMQ;
    newBasisFn->nVtx = nPts;
    WlzValueCopyDVertexToDVertex(newBasisFn->vertices.d2, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex(newBasisFn->sVertices.d2, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    *((double *)(newBasisFn->param)) = deltaSq;
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      *(bV + idY3) = (sPts + idY)->vtX - tD0;
      tD1 = (dPts + idY)->vtY;
      *(*(aA + idY3) + 0) = 1.0;
      *(*(aA + idY3) + 1) = tD0;
      *(*(aA + idY3) + 2) = tD1;
      *(*(aA + 0) + idY3) = 1.0;
      *(*(aA + 1) + idY3) = tD0;
      *(*(aA + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	if(newBasisFn->distMap)
	{
	  sPt.d2 = dPts[idX];
	  tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	}
	else
	{
	  tV0.vtX = (dPts + idX)->vtX - (dPts + idY)->vtX;
	  tV0.vtX *= tV0.vtX;
	  tV0.vtY = (dPts + idX)->vtY - (dPts + idY)->vtY;
	  tV0.vtY *= tV0.vtY;
	  tD0 = tV0.vtX + tV0.vtY;
	}
	tD1 = (tD0 > DBL_EPSILON)? sqrt(tD0 + deltaSq): delta;
	tD2 = 1.0 / tD1;
	idX3 = idX + 3;
	*(*(aA + idY3) + idX3) = tD2;
	*(*(aA + idX3) + idY3) = tD2;
      }
      *(*(aA + idY3) + idY3) = 1.0 / delta;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
       and the y polynomial coefficients. */
    WlzBasisFnMQCoexff2D(newBasisFn, bV,  &extentDB, range, 0, 0);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      *(bV + idY3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients. */
    WlzBasisFnMQCoexff2D(newBasisFn, bV,  &extentDB, range, 1, 0);
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    if(newBasisFn)
    {
      (void )WlzBasisFnFree(newBasisFn);
      newBasisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new inverse-multiquadric basis function.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	delta			Normalized delta value in range
*					[> 0.0 , < 1.0 ].
* \param	prvBasisFn		Previous basis function to be recycled,
* 					may be NULL. Arrays will be used for
* 					the new basis function. Must be freed
* 					after recycling.
* \param	mesh			Used to compute constrained distances.
* 					If NULL Euclidean distances will be
* 					used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnIMQ3DFromCPts(int nPts, WlzDVertex3 *dPts, 
				WlzDVertex3 *sPts, double delta, 
				WlzBasisFn *prvBasisFn,
				WlzCMesh3D *mesh,
				WlzErrorNum *dstErr)
{
  int		tI0,
  		idN,
		idX,
  		idY,
                idX4,
		idY4,
		nSys,
		maxNod,
		newMaxVx;
  double	tD0,
		tD1,
                tD2,
		deltaSq,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDVertex3	tDVx0;
  WlzDBox3	extentDB;
  WlzBasisFn    *newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	stepVx = 10;

#ifdef WLZ_BASISFN_DELTA_ENV
  double	val;
  char		*envStr;

  if(((envStr = getenv("WLZ_BASISFN_DELTA")) != NULL) &&
     (sscanf(envStr, "%lg", &val) == 1))
  {
    delta = val;
  }
#endif /* WLZ_BASISFN_DELTA_ENV */
  aM.core = NULL;
  vM.core = NULL;
  nSys = nPts + 4;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_3D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m3)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_3DIMQ)
      {
        errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
        newMaxVx = prvBasisFn->maxVx;
	if((newBasisFn == NULL) &&
	   ((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
	                                          1)) == NULL))
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	if((errNum == WLZ_ERR_NONE) && (prvBasisFn->maxVx < nPts))
	{
	  newMaxVx = nPts + stepVx;
	  if(((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
				       sizeof(WlzDVertex3) *
				       newMaxVx)) == NULL) ||
	     ((prvBasisFn->distMap != NULL) &&
	      ((prvBasisFn->distMap = (double **)
				      AlcRealloc(prvBasisFn->distMap,
					  sizeof(double *) *
					  newMaxVx)) == NULL)))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  *newBasisFn = *prvBasisFn;
	  prvBasisFn->poly.v = NULL;
	  prvBasisFn->basis.v = NULL;
	  prvBasisFn->vertices.v = NULL;
	  prvBasisFn->sVertices.v = NULL;
	  prvBasisFn->param = NULL;
	  prvBasisFn->mesh.v = NULL;
	  prvBasisFn->distMap = NULL;
	  newBasisFn->maxVx = newMaxVx;
	  newBasisFn->nVtx = nPts;
	  if(newBasisFn->distMap != NULL)
	  {
	    tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
	    for(idN = 0; idN < tI0; ++idN)
	    {
	      if((WlzGeomCmpVtx3D(dPts[idN],
	                          newBasisFn->vertices.d3[idN],
				  DBL_EPSILON) != 0) ||
	         (WlzGeomCmpVtx3D(sPts[idN],
		                  newBasisFn->sVertices.d3[idN],
				  DBL_EPSILON) != 0))
	      {
	        AlcFree(newBasisFn->distMap[idN]);
		newBasisFn->distMap[idN] = NULL;
	      }
	    }
	    for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
	    {
	      AlcFree(newBasisFn->distMap[idN]);
	      newBasisFn->distMap[idN] = NULL;
	    }
	    for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
	    {
	      newBasisFn->distMap[idN] = NULL;
	    }
	  }
	}
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                             1)) == NULL) ||
	 ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex3) * 4)) == NULL) ||
	 ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex3) *
	                                   newMaxVx)) == NULL) ||
	 ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex3) *
	                                      newMaxVx)) == NULL) ||
	 ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex3) *
					       newMaxVx)) == NULL) ||
	 ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        newBasisFn->nVtx = nPts;
	newBasisFn->maxVx = newMaxVx;
	newBasisFn->mesh.m3 = mesh;
	if(newBasisFn->mesh.m3 != NULL)
	{
	  if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
	                                                 newMaxVx)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn3D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m3->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn3D;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  WlzErrorNum	errNum1 = WLZ_ERR_NONE;

	  if(newBasisFn->distMap[idN] == NULL)
	  {
	    if((newBasisFn->distMap[idN] = (double *)
					   AlcMalloc(sizeof(double) *
						     maxNod)) == NULL)
	    {
	      errNum1 = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      errNum1 = WlzCMeshFMarNodes3D(newBasisFn->mesh.m3,
					    newBasisFn->distMap[idN],
					    1, dPts + idN);
	    }
	  }
#ifdef _OPENMP
#pragma omp critical
          {
#endif
            if((errNum == WLZ_ERR_NONE) && (errNum1 != WLZ_ERR_NONE))
	    {
	      errNum = errNum1;
	    }
#ifdef _OPENMP
	  }
#endif
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function the design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
    WlzBasisFnVxExtent3D(&extentDB, dPts, sPts, nPts);
    tD0 = extentDB.xMax - extentDB.xMin;
    tD1 = extentDB.yMax - extentDB.yMin;
    tD2 = extentDB.zMax - extentDB.zMin;
    range = (tD0 > tD1)? tD0: tD1;
    if(tD2 > range)
    {
      range = tD2;
    }
    if(range <= 1.0)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    delta *= range;
    deltaSq = delta * delta;
    newBasisFn->type = WLZ_FN_BASIS_3DIMQ;
    newBasisFn->nVtx = nPts;
    WlzValueCopyDVertexToDVertex3(newBasisFn->vertices.d3, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex3(newBasisFn->sVertices.d3, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    *((double *)(newBasisFn->param)) = deltaSq;
    /* Fill matrix A and matrix b for the x component. */
    for(idY = 0; idY < 4; ++idY)
    {
      for(idX = 0; idX < 4; ++idX)
      {
        *(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      tD0 = (dPts + idY)->vtX;
      tD1 = (dPts + idY)->vtY;
      tD2 = (dPts + idY)->vtZ;
      *(bV + idY4) = (sPts + idY)->vtX - tD0;
      *(*(aA + idY4) + 0) = 1.0;
      *(*(aA + idY4) + 1) = tD0;
      *(*(aA + idY4) + 2) = tD1;
      *(*(aA + idY4) + 3) = tD2;
      *(*(aA + 0) + idY4) = 1.0;
      *(*(aA + 1) + idY4) = tD0;
      *(*(aA + 2) + idY4) = tD1;
      *(*(aA + 3) + idY4) = tD2;
      for(idX = 0; idX < idY; ++idX)
      {
	if(newBasisFn->distMap)
	{
	  sPt.d3 = dPts[idX];
	  tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	}
	else
	{
	  tDVx0.vtX = (dPts + idX)->vtX - (dPts + idY)->vtX;
	  tDVx0.vtX *= tDVx0.vtX;
	  tDVx0.vtY = (dPts + idX)->vtY - (dPts + idY)->vtY;
	  tDVx0.vtY *= tDVx0.vtY;
	  tDVx0.vtZ = (dPts + idX)->vtZ - (dPts + idY)->vtZ;
	  tDVx0.vtZ *= tDVx0.vtZ;
	  tD0 = tDVx0.vtX + tDVx0.vtY + tDVx0.vtZ;
	}
        tD1 = (tD0 > DBL_EPSILON)? sqrt(tD0 + deltaSq): delta;
	tD2 = 1.0 / tD1;
	idX4 = idX + 4;
	*(*(aA + idY4) + idX4) = tD2;
	*(*(aA + idX4) + idY4) = tD2;
      }
      *(*(aA + idY4) + idY4) = 1.0 / delta;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover lambda and the x polynomial coefficients, then set up for mu
     * and the y polynomial coefficients. */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range, 0, 0);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    *(bV + 3) = 0.0;
     for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      *(bV + idY4) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }

  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover mu and the y polynomial coefficients, then set up for nu
     * and the z polynomial coefficients */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range, 1, 0);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    *(bV + 3) = 0.0;
     for(idY = 0; idY < nPts; ++idY)
    {
      idY4 = idY + 4;
      *(bV + idY4) = (sPts + idY)->vtZ - (dPts + idY)->vtZ;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Recover nu and the z polynomial coefficients. */
    WlzBasisFnMQCoeff3D(newBasisFn, bV,  &extentDB, range, 2, 0);
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    if(newBasisFn)
    {
      (void )WlzBasisFnFree(newBasisFn);
      newBasisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new thin plate spline basis function.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
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
*                                       indexed using the node index.
* \param	prevNVx			Previous number of control points -
* 					the number of one dimensional arrays
* 					in the distance map.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnTPS2DFromCPts(int nPts,
				  WlzDVertex2 *dPts, WlzDVertex2 *sPts,
				  WlzBasisFn *prvBasisFn,
				  WlzCMesh2D *mesh,
				  WlzErrorNum *dstErr)
{
  int		tI0,
		idN,
		idX,
		idY,
		idX3,
		idY3,
		nSys,
		maxNod,
                newMaxVx;
  double	tD0,
		tD1,
		range;
  double	*bV = NULL,
  		*wV = NULL;
  double	**aA;
  AlgMatrix	aM,
  		vM;
  WlzVertex	sPt;
  WlzDVertex2	tDVx0;
  WlzDBox2	extentDB;
  WlzBasisFn    *newBasisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int     stepVx = 10;

  aM.core = NULL;
  vM.core = NULL;
  nSys = nPts + 3;
  if(mesh != NULL)
  {
    if((mesh->type != WLZ_CMESH_2D) ||
       ((prvBasisFn != NULL) && (mesh != prvBasisFn->mesh.m2)))
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvBasisFn)
    {
      /* Create a new basis function, reallocate the buffers of the previous
       * basis function as required and then use them for the new basis
       * function. */
      if(prvBasisFn->type != WLZ_FN_BASIS_2DTPS)
      {
        errNum = WLZ_ERR_DOMAIN_TYPE;
      }
      else
      {
        newMaxVx = prvBasisFn->maxVx;
        if(prvBasisFn->maxVx < nPts)
        {
          newMaxVx = nPts + stepVx;
          if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                                    1)) == NULL) ||
             ((prvBasisFn->basis.v = AlcRealloc(prvBasisFn->basis.v,
                                       sizeof(WlzDVertex2) *
                                       newMaxVx)) == NULL) ||
             ((prvBasisFn->vertices.v = AlcRealloc(prvBasisFn->vertices.v,
                                       sizeof(WlzDVertex2) *
                                       newMaxVx)) == NULL) ||
             ((prvBasisFn->sVertices.v = AlcRealloc(prvBasisFn->sVertices.v,
                                       sizeof(WlzDVertex2) *
                                       newMaxVx)) == NULL) ||
             ((prvBasisFn->distMap != NULL) &&
              ((prvBasisFn->distMap = (double **)
                                      AlcRealloc(prvBasisFn->distMap,
                                          sizeof(double *) *
                                          newMaxVx)) == NULL)))
          {
            errNum = WLZ_ERR_MEM_ALLOC;
          }
        }
        if(errNum == WLZ_ERR_NONE)
        {
          *newBasisFn = *prvBasisFn;
          prvBasisFn->poly.v = NULL;
          prvBasisFn->basis.v = NULL;
          prvBasisFn->vertices.v = NULL;
          prvBasisFn->sVertices.v = NULL;
          prvBasisFn->param = NULL;
          prvBasisFn->mesh.v = NULL;
          prvBasisFn->distMap = NULL;
          if(newBasisFn->distMap != NULL)
          {
            tI0 = ALG_MIN(newBasisFn->nVtx, prvBasisFn->nVtx);
            for(idN = 0; idN < tI0; ++idN)
            {
              if((WlzGeomCmpVtx2D(dPts[idN],
                                  newBasisFn->vertices.d2[idN],
                                  DBL_EPSILON) == 0) ||
                  (WlzGeomCmpVtx2D(sPts[idN],
                                   newBasisFn->sVertices.d2[idN],
                                   DBL_EPSILON) == 0))
              {
                AlcFree(newBasisFn->distMap[idN]);
                newBasisFn->distMap[idN] = NULL;
              }
            }
            for(idN = tI0; idN < prvBasisFn->nVtx; ++idN)
            {
              AlcFree(newBasisFn->distMap[idN]);
              newBasisFn->distMap[idN] = NULL;
            }
            for(idN = prvBasisFn->nVtx; idN < newBasisFn->maxVx; ++idN)
            {
              newBasisFn->distMap[idN] = NULL;
            }
          }
        }
      }
    }
    else
    {
      /* Create a new basis function, allocate the buffers of the new
       * basis function as required. */
      newMaxVx = nPts + stepVx;
      if(((newBasisFn = (WlzBasisFn *)AlcCalloc(sizeof(WlzBasisFn),
                                                       1)) == NULL) ||
         ((newBasisFn->poly.v = AlcMalloc(sizeof(WlzDVertex2) * 3)) == NULL) ||
         ((newBasisFn->basis.v = AlcMalloc(sizeof(WlzDVertex2) *
                                           newMaxVx)) == NULL) ||
         ((newBasisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex2) *
                                              newMaxVx)) == NULL) ||
         ((newBasisFn->sVertices.v = AlcMalloc(sizeof(WlzDVertex2) *
                                              newMaxVx)) == NULL) ||
         ((newBasisFn->param = AlcMalloc(sizeof(double))) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        newBasisFn->nVtx = nPts;
        newBasisFn->maxVx = newMaxVx;
        newBasisFn->mesh.m2 = mesh;
        if(newBasisFn->mesh.m2 != NULL)
        {
          if((newBasisFn->distMap = (double **)AlcCalloc(sizeof(double *),
                                                         newMaxVx)) == NULL)
          {
            errNum = WLZ_ERR_MEM_ALLOC;
          }
        }
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && (mesh != NULL))
  {
    if(newBasisFn->distMap == NULL)
    {
      newBasisFn->distFn = WlzBasisFnEucDistFn2D;
    }
    else
    {
      /* Allocate and compute the distance maps. */
      maxNod = newBasisFn->mesh.m2->res.nod.maxEnt;
      newBasisFn->distFn = WlzBasisFnMapDistFn2D;
      for(idN = 0; idN < newBasisFn->nVtx; ++idN)
      {
        if(newBasisFn->distMap[idN] == NULL)
        {
          if((newBasisFn->distMap[idN] = (double *)
                                         AlcMalloc(sizeof(double) *
                                                   maxNod)) == NULL)
          {
            errNum = WLZ_ERR_MEM_ALLOC;
          }
          else
          {
            errNum = WlzCMeshFMarNodes2D(newBasisFn->mesh.m2,
                                         newBasisFn->distMap[idN],
                                         1, dPts + idN);
          }
          if(errNum != WLZ_ERR_NONE)
          {
            break;
          }
        }
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate matrices for solving basis function the design equation. */
    if(((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
       ((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
       ((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
       ((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    aA = aM.rect->array;
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
    newBasisFn->type = WLZ_FN_BASIS_2DTPS;
    newBasisFn->nVtx = nPts;
    WlzValueCopyDVertexToDVertex(newBasisFn->vertices.d2, dPts, nPts);
    if(newBasisFn->mesh.v != NULL)
    {
      WlzValueCopyDVertexToDVertex(newBasisFn->sVertices.d2, sPts, nPts);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    newBasisFn->nPoly = 2;
    newBasisFn->nBasis = nPts;
    for(idY = 0; idY < 3; ++idY)
    {
      for(idX = 0; idX < 3; ++idX)
      {
	*(*(aA + idY) + idX) = 0.0;
      }
      *(bV + idY) = 0.0;
    }
    for(idY = 0; idY < nPts; ++idY)
    {
      idY3 = idY + 3;
      tD0 = (dPts + idY)->vtX;
      *(bV + idY3) = (sPts + idY)->vtX - tD0;
      tD0 = (tD0 - extentDB.xMin) / range;
      tD1 = ((dPts + idY)->vtY - extentDB.yMin) / range;
      *(*(aA + idY3) + 0) = 1.0;
      *(*(aA + idY3) + 1) = tD0;
      *(*(aA + idY3) + 2) = tD1;
      *(*(aA + 0) + idY3) = 1.0;
      *(*(aA + 1) + idY3) = tD0;
      *(*(aA + 2) + idY3) = tD1;
      for(idX = 0; idX < idY; ++idX)
      {
	if(newBasisFn->distFn)
	{
	  sPt.d2 = dPts[idX];
	  tD0 = newBasisFn->distFn(newBasisFn, idY, sPt);
	  tD0 *= tD0;
	}
	else
	{
	  tDVx0.vtX = ((dPts + idX)->vtX - (dPts + idY)->vtX) / range;
	  tDVx0.vtX *= tDVx0.vtX;
	  tDVx0.vtY = ((dPts + idX)->vtY - (dPts + idY)->vtY) / range;
	  tDVx0.vtY *= tDVx0.vtY;
	  tD0 = tDVx0.vtX + tDVx0.vtY;
	}
	tD1 = (tD0 > DBL_EPSILON)? tD0 * log(tD0): 0.0;
	idX3 = idX + 3;
	*(*(aA + idY3) + idX3) = tD1;
	*(*(aA + idX3) + idY3) = tD1;
      }
      *(*(aA + idY3) + idY3) = 0.0;
    }
    /* Perform singular value decomposition of matrix A. */
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the X polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnTPS2DCoef(newBasisFn, bV,  &extentDB, range, 1);
    *(bV + 0) = 0.0;
    *(bV + 1) = 0.0;
    *(bV + 2) = 0.0;
    for(idY = 0; idY < nPts; ++idY)
    {
      *(bV + idY + 3) = (sPts + idY)->vtY - (dPts + idY)->vtY;
    }
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzBasisFnTPS2DCoef(newBasisFn, bV,  &extentDB, range, 0);
  }
  AlcFree(bV);
  AlcFree(wV);
  AlgMatrixFree(aM);
  AlgMatrixFree(vM);
  if(errNum != WLZ_ERR_NONE)
  {
    if(newBasisFn)
    {
      (void )WlzBasisFnFree(newBasisFn);
      newBasisFn = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newBasisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Creates a new 3D multi order spline basis function:
*		\f[
		\phi(r) = \frac{1}{4 \pi \delta^2 r}
			  (1 -
			   \frac{w}{w - v} e^{- \sqrt{v} r} +
			   \frac{v}{w - v} e^{- \sqrt{w} r})
		\f]
		\f[
		v = \frac{1 + \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		w = \frac{1 - \sqrt{1 - 4 \pi \tau^2 \delta^2}}{2 \tau^2}
		\f]
		\f[
		f(\mathbf{x}) = P(\mathbf{x}) +
			        \sum_{i=1}^{n}{\lambda_i \phi(r_i)}
		\f]
		\f[
		r_i = |\mathbf{x} - \mathbf{x'_i}|
		\f]
*		The multi order spline may either be an exact interpolating
*		function or one which approximates the given control points.
*		In the approximating case values for the regularization
*		parameters should be given, with the regularization
*		parameters \f$\alpha_i\f$ given by:
*		\f[
		H[f] = \beta[f] +
		       \sum_{i=1}^{n}{\frac{(y_i - f(x))^2}{\alpha_i}}
		\f]
		Giving the design equation:
		\f[
		\left(
		\begin{array}{cccccccc}
		0 & 0 & 0 & 0 & 1 & 1 & \cdots & 1 \\
		0 & 0 & 0 & 0 & x_1 & x_2 & \cdots & x_n \\
		0 & 0 & 0 & 0 & y_1 & y_2 & \cdots & y_n \\
		0 & 0 & 0 & 0 & z_1 & z_2 & \cdots & z_n \\
		1 & x_1 & y_1 & z_1 & \phi(r_{11}) + \alpha_1 & \phi(r_{12}) &
		          \cdots & \phi(r_{1n}) \\
		\dotfill \\
		1 & x_n & y_n & z_n & \phi(r_{n1}) & \phi(r_{n2}) &
		          \cdots & \phi(r_{nn} + \alpha_n)
		\end{array}
		\right)

		\left(
		\begin{array}{c}
		a_o \\
		a_1 \\
		a_2 \\
		a_3 \\
		\lambda_1 \\
		\lambda_1 \\
		\vdots \\
		\lambda_n
		\end{array}
		\right)
		
		=

		\left(
		\begin{array}{c}
		0 \\
		0 \\
		0 \\
		0 \\
		f(\mathbf{x_1}) \\
		f(\mathbf{x_1}) \\
		\vdots \\
		f(\mathbf{x_n})
		\end{array}
		\right)
		\f]
*		The vertex values with components  \f$x\f$, \f$y \f$ and
*		\f$z\f$ are used to solve for the four polynomial
*		coefficients and the three basis function weights 
*		\f$\lambda\f$, \f$\mu\f$ and \f$\nu\f$.
*		The given values of \f$\lambda\f$ and \f$\tau\f$ are
*		constrained by \f$\lambda > 0, \tau > 0\f$ and
*		\f$\lambda \tau < 0.5\f$.
* \param	nPts			Number of control point pairs.
* \param	dPts			Destination control points.
* \param	sPts			Source control points.
* \param	alpha			Regularization parameter for the
*					control points. If NULL all
*					regularization parameters are set to
*					zero for exact interpolation.
* \param	param			Smoothness parameters \f$\delta\f$
*					and \f$\tau\f$ in that order,
*					but may be NULL in which case default
*					parameter values will be used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnMOS3DFromCPts(int nPts,
				  WlzDVertex2 *dPts, WlzDVertex2 *sPts,
				  double *alpha, double *param,
				  WlzErrorNum *dstErr)
{
  WlzBasisFn	*basisFn = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  /* TODO */
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(basisFn);
}

/*!
* \return	New basis function.
* \ingroup	WlzFunction
* \brief	Computes a new 3D multi order spline basis function
*		which either interpolates or approximates the given
*		scalar values. See WlzBasisFnMOS3DFromCPts() for details
*		of the multi order spline.
*		The given values of \f$\lambda\f$ and \f$\tau\f$ are
*		constrained by \f$\lambda > 0, \tau > 0\f$ and
*		\f$\lambda \tau < 0.5\f$.
* \param	nPts			Number of values.
* \param	cPts			Positions of the control points.
* \param	cVal			Values at the control points.
* \param	alpha			Regularization parameter for the
*					control points. If NULL all
*					regularization parameters are set to
*					zero for exact interpolation.
* \param	param			Smoothness parameters \f$\delta\f$
*					and \f$\tau\f$ in that order,
*					but may be NULL in which case default
*					parameter values will be used.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzBasisFn *WlzBasisFnScalarMOS3DFromCPts(int nPts,
				  WlzDVertex3 *cPts, double *cVal,
				  double *alpha, double *param,
				  WlzErrorNum *dstErr)
{
  int		idX,
		idY,
		idX1,
		idY1,
		nSys;
  double	tD0,
		rad,
		phi0,
		delta,
		tau;
  double	*bV = NULL;
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
  double	wMax;
  double  	*wV = NULL;
  double  	**vA;
  AlgMatrix	vM;
#else
  double	*xV = NULL;
  AlgMatrix	wM;
#endif
  double	**aA = NULL;
  AlgMatrix	aM;
  WlzBasisFn	*basisFn = NULL;
  WlzDVertex3	tV0,
  		tV1,
		tV2;
  WlzHistogramDomain *evalData = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  aM.core = NULL;
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
  vM.core = NULL;
#else
  wM.core = NULL;
#endif
  delta = *(param + 0);
  tau = *(param + 1);
  if((delta < DBL_EPSILON) || (tau < DBL_EPSILON) ||
     (delta * tau > (0.5 - DBL_EPSILON)))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    evalData = WlzBasisFnScalarMOS3DEvalTb(nPts, cPts, delta, tau, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    nSys = nPts + 1;
    if(
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
        ((wV = (double *)AlcCalloc(sizeof(double), nSys)) == NULL) ||
	((vM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
#else
	((wM.rect = AlgMatrixRectNew(4, nSys, NULL)) == NULL) ||
	((xV = (double *)AlcCalloc(nSys, sizeof(double))) == NULL) ||
#endif
	((bV = (double *)AlcMalloc(sizeof(double) * nSys)) == NULL) ||
	((aM.rect = AlgMatrixRectNew(nSys, nSys, NULL)) == NULL) ||
	((basisFn = (WlzBasisFn *)
	  AlcCalloc(sizeof(WlzBasisFn), 1)) == NULL) ||
	((basisFn->poly.v = AlcMalloc(sizeof(double) * 1)) == NULL) ||
	((basisFn->basis.v = AlcMalloc(sizeof(double) * nPts)) == NULL) ||
	((basisFn->vertices.v = AlcMalloc(sizeof(WlzDVertex3) *
	                                  nPts)) == NULL) ||
	((basisFn->param = AlcMalloc(sizeof(double) * 2)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
    vA = vM.rect->array;
#endif
    aA = aM.rect->array;
    *((double *)(basisFn->param) + 0) = delta;
    *((double *)(basisFn->param) + 1) = tau;
    basisFn->type = WLZ_FN_BASIS_SCALAR_3DMOS;
    basisFn->nPoly = 1;
    basisFn->nBasis = nPts;
    basisFn->nVtx = nPts;
    basisFn->evalFn = WlzBasisFnScalarMOS3DEvalFn; /* Don't set to avoid LUT */
    basisFn->evalData = evalData;
    WlzValueCopyDVertexToDVertex3(basisFn->vertices.d3, cPts, nPts);
    *(*(aA + 0) + 0) = 0.0;
    for(idX = 0; idX < nPts; ++idX)
    {
      idX1 = idX + 1;
      *(*(aA + 0) + idX1) = 1.0;
    }
    *(bV + 0) = 0.0;
    phi0 = basisFn->evalFn?
    	   basisFn->evalFn((void *)basisFn, 0.0):
	   WlzBasisFnValueMOSPhi(0.0, delta, tau);
    for(idY = 0; idY < nPts; ++idY)
    {
      idY1 = idY + 1;
      tV0 = *(cPts + idY);
      *(bV + idY1) = *(cVal + idY);
      *(*(aA + idY1) + 0) = 1.0;
      for(idX = 0; idX < idY; ++idX)
      {
	idX1 = idX + 1;
	tV1 = *(cPts + idX);
	WLZ_VTX_3_SUB(tV2, tV1, tV0);
	rad = WLZ_VTX_3_LENGTH(tV2);
        tD0 = basisFn->evalFn? 
	      basisFn->evalFn((void *)basisFn, rad):
	      WlzBasisFnValueMOSPhi(rad, delta, tau);
	*(*(aA + idY1) + idX1) = *(*(aA + idX1) + idY1) = tD0;
      }
      *(*(aA + idY1) + idY1) = phi0 + *(alpha + idY);
    }
#ifdef WLZ_BASISFNSCALARMOS3DFROMCPTS_DEBUG
    {
      FILE	*fP;
      fP = fopen("DEBUG_aM.num", "w");
      AlgMatrixWriteAscii(aM, fP);
      fclose(fP);
      fP = fopen("DEBUG_bV.num", "w");
      AlcDouble1WriteAsci(fP, bV, nSys);
      fclose(fP);
    }
#endif
    /* Solve the design equation. */
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
    errNum = WlzErrorFromAlg(AlgMatrixSVDecomp(aM, wV, vM));
#else
    errNum = WlzErrorFromAlg(
    	     AlgMatrixCGSolve(aM, xV, bV, wM, NULL, NULL, 0.000001,
	                      1000, NULL, NULL));
#endif
  }
#ifdef WLZ_BASISFN_MOS_SOLVER_SVD
  if(errNum == WLZ_ERR_NONE)
  {
    /* Edit the singular values. */
    WlzBasisFnEditSV(nSys, wV);
    /* Solve for lambda and the polynomial coefficients. */
    errNum = WlzErrorFromAlg(AlgMatrixSVBackSub(aM, wV, vM, bV));
  }
#ifdef WLZ_BASISFNSCALARMOS3DFROMCPTS_DEBUG
    {
      FILE	*fP;
      fP = fopen("DEBUG_bMxBS.num", "w");
      AlcDouble1WriteAsci(fP, bV, nSys);
      fclose(fP);
    }
#endif
  if(errNum == WLZ_ERR_NONE)
  {
    *(double *)(basisFn->poly.v) = *bV;	              /* Only constant used. */
    for(idY = 0; idY < nPts; ++idY)
    {
      idY1 = idY + 1;
      *((double *)(basisFn->basis.v) + idY) = *(bV + idY1);
    }
  }
  if(wV)
  {
    AlcFree(wV);
  }
  AlgMatrixFree(vM);
#else
  if(errNum == WLZ_ERR_NONE)
  {
    *(double *)(basisFn->poly.v) = *xV;	              /* Only constant used. */
    for(idY = 0; idY < nPts; ++idY)
    {
      idY1 = idY + 1;
      *((double *)(basisFn->basis.v) + idY) = *(xV + idY1);
    }
  }
  if(xV)
  {
    AlcFree(xV);
  }
  AlgMatrixFree(wM);
#endif
  AlcFree(bV);
  AlgMatrixFree(aM);
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
* \return	New histogram domain look up table.
* \ingroup	WlzFunction
* \brief	Computes a look up table for the multiorder spline
*		over the range of distances between the given vertices.
* \param	nPts			Number of vertices.
* \param	cPts			The given vertices.
* \param	delta			Multiorder spline \f$\delta\f$.
* \param	tau			Multiorder spline \f$\tau\f$.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzHistogramDomain *WlzBasisFnScalarMOS3DEvalTb(int nPts,
				WlzDVertex3 *cPts, double delta, double tau,
				WlzErrorNum *dstErr)
{
  int		idx,
  		lutMax;
  double	tD0,
  		tD1,
		v,
		w,
		rv,
		rw,
		norm;
  WlzDVertex3	tV0;
  WlzDBox3	extentDB;
  WlzHistogramDomain *lut = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double lutStep = 0.1;

  extentDB = WlzBoundingBoxVtx3D(nPts, cPts, NULL);
  tV0.vtX = extentDB.xMax - extentDB.xMin;
  tV0.vtY = extentDB.yMax - extentDB.yMin;
  tV0.vtZ = extentDB.zMax - extentDB.zMin;
  lutMax = (int )(ceil(WLZ_VTX_3_LENGTH(tV0)) / lutStep) + 2;
  lut = WlzMakeHistogramDomain(WLZ_HISTOGRAMDOMAIN_FLOAT, lutMax, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    lut->nBins = lutMax;
    lut->origin = 0.0;
    lut->binSize = lutStep;
    tD0 = 2.0 * tau * delta;
    tD1 = sqrt(1 - (tD0 * tD0));
    tD0 = 2.0 * tau * tau;
    v = (1 + tD1) / tD0;
    w = (1 - tD1) / tD0;
    rv = sqrt(v);
    rw = sqrt(w);
    norm = 1.0 / (4.0 * ALG_M_PI * delta * delta);
    for(idx = 0; idx < lutMax; ++idx)
    {
      *(lut->binValues.dbp + idx) = WlzBasisFnValueMOSPhiPC(idx * lutStep,
      					v, w, delta, rv, rw, norm);
    }
  }
  return(lut);
}

/*!
* \return
* \brief	Computes a multiorder spline basis function using a look up
* 		table.
* \param	basisFn			Given basis function with non NULL
*					look up table.
* \param	rad			Radial distance to use.
*/
static double	WlzBasisFnScalarMOS3DEvalFn(void *basisFn, double rad)
{
  int		idx0,
  		idx1;
  double	delta,
  		tau,
		radb,
  		phi,
		binSizeR;
  WlzHistogramDomain *lut;

  lut = ((WlzBasisFn *)basisFn)->evalData;
  binSizeR = 1.0 / lut->binSize;
  /* lut->origin == 0.0 */
  if(rad <= DBL_EPSILON)
  {
    phi = *(lut->binValues.dbp);
  }
  else
  {
    radb = rad * binSizeR;
    if(radb >= lut->nBins)
    {
      delta = *((double *)(((WlzBasisFn *)basisFn)->param) + 0);
      tau = *((double *)(((WlzBasisFn *)basisFn)->param) + 1);
      phi = WlzBasisFnValueMOSPhi(rad, delta, tau);
    }
    else
    {
      idx0 = (int )floor(radb);
      idx1 = (int )ceil(radb);
      if(idx0 == idx1)
      {
	phi = *(lut->binValues.dbp + idx0);
      }
      else
      {
	/* Use linear interpolation. */
	phi = (*(lut->binValues.dbp + idx0) * ((double )idx1 - radb) +
	       *(lut->binValues.dbp + idx1) * (radb - (double )idx0));
      }
    }
  }
  return(phi);
}

/*!
* \return	void
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
* \ingroup      WlzFunction
* \brief        Computes the extent (bounding volume) of two vectors
*		of verticies.
*
* \param    extentDB	Pointer for the extent of the verticies.
* \param    vx0	First vector of verticies.
* \param    vx1	Second vector of verticies.
* \param    nPts	Number of verticies in each vector.
* \par      Source:
*                WlzBasisFn.c
*/
static void	WlzBasisFnVxExtent3D(WlzDBox3 *extentDB,
				   WlzDVertex3 *vx0, WlzDVertex3 *vx1,
				   int nPts)
{
  double	tD0,
		tD1,
		tD2;

  extentDB->xMin = extentDB->xMax = vx0->vtX;
  extentDB->yMin = extentDB->yMax = vx0->vtY;
  extentDB->zMin = extentDB->zMax = vx0->vtZ;
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
    if(tD1 > extentDB->yMax)
    {
      extentDB->yMax = tD1;
    }
    
    if((tD0 = vx0->vtZ) > (tD1 = vx1->vtZ))
    {
      tD2 = tD0;
      tD0 = tD1;
      tD1 = tD2;
    }
    if(tD0 < extentDB->zMin)
    {
      extentDB->zMin = tD0;
    }
    if(tD1 > extentDB->zMax)
    {
      extentDB->zMax = tD1;
    }
     
    ++vx0;
    ++vx1;
  }
}

/*!
* \return	void
* \ingroup	WlzFunction
* \brief	Extracts the Gaussian coefficients from the given column
* 		vector using the given extent and range for function.
* \param	basisFn			Allocated basis function to
* 					be filled in.
* \param	vec			Given column vector.
* \param	forX			True if the coefficients are for the x
* 					coordinate.
*/
static void	WlzBasisFnGauss2DCoef(WlzBasisFn *basisFn, double *vec,
				      int forX)
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
* \return	void
* \ingroup	WlzFunction
* \brief	Extracts the multiquadric coefficients from the given
*		column vector and optionaly rescales them using the given
*		extent and range for the function.
* \param	basisFn			Allocated basis function
*					to be filled in.
* \param	vec			Given column vector.
* \param	extentDB		Extent of the vertices, not used if
* 					rescale == 0.
* \param	range			Range of the vertices, not used if
* 					rescale == 0.
* \param	component		Value indicating which coefficients
* 					are to be recovered: 0 - x or 1 - y.
* \param	rescale			Rescale the coeeficients using
*				        the range and extent if non-zero.
*/
static void	WlzBasisFnMQCoexff2D(WlzBasisFn *basisFn,
				     double *vec, WlzDBox2 *extentDB,
				     double range, int component,
				     int rescale)
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
  if(component == 0)
  {
    if(rescale)
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
      polyVxP++->vtX = vec0;
      polyVxP++->vtX = vec1;
      polyVxP->vtX = vec2;
      for(idN = 0; idN < basisFn->nBasis; ++idN)
      {
        basisVxP++->vtX = *vecP++;
      }
    }
  }
  else
  {
    if(rescale)
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
    else
    {
      polyVxP++->vtY = vec0;
      polyVxP++->vtY = vec1;
      polyVxP->vtY = vec2;
      for(idN = 0; idN < basisFn->nBasis; ++idN)
      {
        basisVxP++->vtY = *vecP++;
      }
    }
  }
}

/*! 
* \ingroup      WlzFunction
* \brief        Extracts the multiquadric coefficients from the
*		given column vector using the given extent and range
*		for transformation to pixel space.
*
* \param    	basis			Allocated basis function transform to
* 					be filled in.
* \param    	vec			Given column vector.
* \param    	extentDB		Extent of the verticies, not used if
* 					rescale == 0.
* \param    	range			Range of the verticies, not used if
* 					rescale == 0.
* \param    	component		Value indicating which coefficients
* 					are to be recovered: 0 - x, 1 - y,
* 					2 - z.
* \param	rescale			Rescale the coeeficients using
*				        the range and extent if non-zero.
*/
static void	WlzBasisFnMQCoeff3D(WlzBasisFn *basis,
				  double *vec, WlzDBox3 *extentDB,
				  double range, int component, int rescale)
{
  int		idN;
  double 	vec0,
  		vec1,
  		vec2,
                vec3;
  double	*vecP;
  WlzDVertex3	*basisVxP,
  		*polyVxP;

  vec0 = *vec;
  vec1 = *(vec + 1);
  vec2 = *(vec + 2);
  vec3 = *(vec + 3);
  vecP = vec + 4;
  basisVxP = basis->basis.d3;
  polyVxP  = basis->poly.d3;

  /* Change back to normal coordinate (the range is the C constant). */
  switch(component)
  {
    case 0: /* x */
      if(rescale)
      {
	polyVxP++->vtX = vec0 -
			 (((vec1 * extentDB->xMin) +
			   (vec2 * extentDB->yMin) +
			   (vec3 * extentDB->zMin)) / range);
	polyVxP++->vtX = vec1 / range;
	polyVxP++->vtX = vec2 / range;
	polyVxP->vtX   = vec3 / range;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtX = *vecP++ / range;
	}
      }
      else
      {
	polyVxP++->vtX = vec0;
	polyVxP++->vtX = vec1;
	polyVxP++->vtX = vec2;
	polyVxP->vtX   = vec3;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtX = *vecP++;
	}
      }
      break;
    case 1: /* y */
      if(rescale)
      {
	polyVxP++->vtY = vec0 -
			 (((vec1 * extentDB->xMin) +
			   (vec2 * extentDB->yMin) +
			   (vec3 * extentDB->zMin)) / range);
	polyVxP++->vtY = vec1 / range;
	polyVxP++->vtY = vec2 / range;
	polyVxP->vtY   = vec3 / range;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtY = *vecP++ / range;
	}
      }
      else
      {
	polyVxP++->vtY = vec0;
	polyVxP++->vtY = vec1;
	polyVxP++->vtY = vec2;
	polyVxP->vtY   = vec3;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtY = *vecP++;
	}
      }
      break;
    case 2: /* z */
      if(rescale)
      {
	polyVxP++->vtZ = vec0 -
			 (((vec1 * extentDB->xMin) +
			   (vec2 * extentDB->yMin) +
			   (vec3 * extentDB->zMin)) / range);
	polyVxP++->vtZ = vec1 / range;
	polyVxP++->vtZ = vec2 / range;
	polyVxP->vtZ   = vec3 / range;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtZ = *vecP++ / range;
	}
      }
      else
      {
	polyVxP++->vtZ = vec0;
	polyVxP++->vtZ = vec1;
	polyVxP++->vtZ = vec2;
	polyVxP->vtZ   = vec3;
	for(idN = 0; idN < basis->nBasis; ++idN)
	{
	  basisVxP++->vtZ = *vecP++;
	}
      }
      break;
    default:
      break;
  }
}

/*!
* \return	void
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
* \return	Displacement due to reduced 2D polynomial.
* \ingroup	WlzFunction
* \brief	Computes the value of the reduced 2D polynomial
* 		used by the TPS, MQ, IMQ and Gauss basis functions.
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


/*!
* \return	Displacement due to reduced polynomial.
* \ingroup	WlzFunction
* \brief	Computes the value of the reduced polynomial
* 		used by the TPS, MQ, IMQ and Gauss basis functions.
* \param	poly			Given polynomial coefficients.
* \param	srcVx			Source vertex.
*/
static WlzDVertex3 WlzBasisFnValueRedPoly3D(WlzDVertex3 *poly,
					WlzDVertex3 srcVx)
{
  WlzDVertex3	dspVx;
  /* a_0 , b_0 , c_0 */
  dspVx.vtX = poly->vtX;
  dspVx.vtY = poly->vtY;
  dspVx.vtZ = poly->vtZ;
   ++poly;
  /* a_1, b_1 , c_1 */
  dspVx.vtX += poly->vtX * srcVx.vtX;
  dspVx.vtY += poly->vtY * srcVx.vtX;
  dspVx.vtZ += poly->vtZ * srcVx.vtX;
   ++poly;
  /* a_2, b_2, c_2 */
  dspVx.vtX += poly->vtX * srcVx.vtY;
  dspVx.vtY += poly->vtY * srcVx.vtY;
  dspVx.vtZ += poly->vtZ * srcVx.vtY;
   ++poly;
  /* a_3, b_3, c_3 */
  dspVx.vtX += poly->vtX * srcVx.vtZ;
  dspVx.vtY += poly->vtY * srcVx.vtZ;
  dspVx.vtZ += poly->vtZ * srcVx.vtZ;
  return(dspVx);
}

/*!
* \return	void
* \ingroup	WlzFunction
* \brief	Edits the 1D array of potentialy singular values,
*		maxing those that are close to zero zero.
* \param	n			Number of values.
* \param	vV			Array of values.
*/
static void	WlzBasisFnEditSV(int n, double *vV)
{
  int		idN;
  double	thresh,
  		vMax = 0.0;
  const double	tol = 1.0e-09;

  for(idN = 0; idN < n; ++idN)
  {
    if(*(vV + idN) > vMax)
    {
      vMax = *(vV + idN);
    }
  }
  thresh = tol * vMax;
  for(idN = 0; idN < n; ++idN)
  {
    if(*(vV + idN) < thresh)
    {
      *(vV + idN) = 0.0;
    }
  }
}
