#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzBSpline.c
* \author       Bill Hill
* \date         June 2020
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Functions to create and compute B-spline domain objects.
* \ingroup	WlzFeatures
*/

#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

static WlzBSpline		*WlzBSplineFromVertices(
				  WlzVertexType vType,
				  int nV,
				  WlzVertexP vtx,
				  int k,
				  int periodic,
				  double sm,
				  WlzErrorNum *dstErr);
/*!
* \return	New Woolz B-spline domain or NULL on error.
* \ingroup	WlzAllocation
* \brief	Creates a new B-spline domain. This can be freed using
* 		using WlzFreeDomain().
* \param	type		Type of B-spline domain: Must be either
* 				WLZ_BSPLINE_C2D or WLZ_BSPLINE_C3D.
* \param	order		Must be in the range [1-5].
* \param	maxKnots	Maximum number of knots per dimension.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzBSpline			*WlzMakeBSpline(
				  WlzObjectType type,
				  int order,
				  int maxKnots,
				  WlzErrorNum *dstErr)
{
  WlzBSpline	*bs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((type != WLZ_BSPLINE_C2D) && (type != WLZ_BSPLINE_C3D))
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((order < WLZ_BSPLINE_ORDER_MIN) || (order > WLZ_BSPLINE_ORDER_MAX) ||
          (maxKnots < 2 * (order + 1)))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  else
  {
    int		dim;

    dim = (type == WLZ_BSPLINE_C2D)? 2: 3;
    if((bs = (WlzBSpline *)
    	     AlcCalloc(sizeof(WlzBSpline) +
	               (maxKnots * (dim + 1)) * sizeof(double), 1)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      bs->type = type;
      bs->order = order;
      bs->maxKnots = maxKnots;
      bs->knots = (double *)(bs + 1);
      bs->coefficients = bs->knots + maxKnots;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bs);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzAllocation
* \brief	Frees a B-spline domain.
* \param	bs		Given B-spline domain.
*/
WlzErrorNum			WlzFreeBSpline(
				  WlzBSpline *bs)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bs)
  {
    if(WlzUnlink(&(bs->linkcount), &errNum))
    {
      AlcFree((void *)bs);
    }
  }
  return(errNum);
}

/*!
* \return	New B-spline domain or NULL on error.
* \ingroup	WlzAllocation
* \brief	Copies a B-spline domain.
* \param	srcBS		Given source B-spline domain.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzBSpline			*WlzBSplineCopy(
				  WlzBSpline *srcBS,
				  WlzErrorNum *dstErr)
{
  int		dim = 0;
  WlzBSpline	*cpyBS = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(srcBS == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    dim = (srcBS->type == WLZ_BSPLINE_C2D)? 2: 3;
    cpyBS = WlzMakeBSpline(srcBS->type, srcBS->order, srcBS->maxKnots,
    			   &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cpyBS->nKnots = srcBS->nKnots;
    WlzValueCopyDoubleToDouble(cpyBS->knots, srcBS->knots, srcBS->nKnots);
    WlzValueCopyDoubleToDouble(cpyBS->coefficients, srcBS->coefficients,
    			       srcBS->maxKnots * dim);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cpyBS);
}

/*!
* \return	New Woolz B-spline domain or NULL on error.
* \ingroup	WlzFeatures
* \brief	Creates a new B-spline domain from the given object.
* \param	gObj		Given object which can give a list
* 				of vertices via WlzVerticesFromObj().
* 				The B-spline is fitted to all vertices.
* \param	order		Must be in the range [1-5].
* \param	closed		If true a periodic B-spline will be computed.
* \param	sm		Smoothing parameter, with value 0.0 for no
* 				smoothing.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzBSpline			*WlzBSplineFromObj(
				 WlzObject *gObj,
				 int order,
				 int closed,
				 double sm,
				 WlzErrorNum *dstErr)
{
  int		nVtx = 0;
  WlzVertexType vType = WLZ_VERTEX_ERROR;
  WlzVertexP	vtx = {0};
  WlzBSpline	*bs = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(gObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_POINTS:
	{
	  WlzPoints	*pts;

	  pts = gObj->domain.pts;
	  vType = WlzPointsVertexType(pts->type, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    bs = WlzBSplineFromVertices(vType, pts->nPoints, pts->points,
	        order, closed, sm, &errNum);
	  }
        }
      default:
	vtx = WlzVerticesFromObj(gObj, NULL, &nVtx, &vType, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  bs = WlzBSplineFromVertices(vType, nVtx, vtx, order, closed,
	      sm, &errNum);
	}
	break;
    }
  }
  AlcFree(vtx.v);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bs);
}

/*!
* \return	New points domain or NULL on error.
* \ingroup	WlzFeatures
* \brief	Creates a new Woolz points domain bt evaluating the given
* 		B-spline at equal intervals along the parametrised curve.
* \param	bs		Given B-spline domain.
* \param	n		Number of points to be evaluated.
* \param	dstErr		Destination error pointer, may be NULL.
*/
WlzPoints			*WlzBSplineEvalPoints(
				  WlzBSpline *bs,
				  int n,
				  WlzErrorNum *dstErr)
{
  int		nn = 0;
  WlzObjectType pType = WLZ_NULL;
  WlzPoints	*pts = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bs == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(bs->type)
    {
      case WLZ_BSPLINE_C2D:
        pType = WLZ_POINTS_2D;
	break;
      case WLZ_BSPLINE_C3D:
        pType = WLZ_POINTS_3D;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzVertexP	nullVtx = {0};

    /* Take care here passing n = 0 indicates use knots avoiding the initial
     * and final padding. */
    nn = (n == 0)? bs->nKnots - 2 * (bs->order): n;
    pts = WlzMakePoints(pType, 0, nullVtx, nn, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzBSplineEval(bs, n, NULL, 0, pts->points);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    pts->nPoints = nn;
  }
  else if(pts)
  {
    WlzDomain	dom;

    dom.pts = pts;
    WlzFreeDomain(dom);
    pts = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pts);
}

/*!
* \return	New Woolz B-spline domain or NULL on error.
* \ingroup	WlzFeatures
* \brief	Creates a new B-spline domain by fitting a B-spline to
* 		the given vertices.
* \param	vType		Type of given vertices.
* \param	nV		Number of given vertices.
* \param	vtx		Given vertices.
* \param	k		Spline order.
* \param	periodic	If true a periodic B-spline will be computed.
* \param	sm		Smoothing parameter, with value 0.0 for no
* 				smoothing.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static WlzBSpline		*WlzBSplineFromVertices(
				  WlzVertexType vType,
				  int nV,
				  WlzVertexP vtx,
				  int k,
				  int periodic,
				  double sm,
				  WlzErrorNum *dstErr)
{
  int		dim,
  		nest,
		nC = 0,
		nT = 0;
  WlzBSpline	*bs = NULL;
  int		*iWrk = NULL;
  double	*c = NULL,
		*t = NULL,
		*u = NULL,
		*w = NULL,
		*x = NULL,
		*dWrk = NULL;
		
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dim = WlzVertexDim(vType, &errNum);
  nest = nV + 2 * (k + 1);
  nC = dim * nest;
  if(errNum == WLZ_ERR_NONE)
  {
    int		lwrk;

    lwrk = (nV * (k + 1)) + (nest * (7 + dim + (5 * k)));
    if(((c = (double *)AlcCalloc(nC, sizeof(double))) == NULL) ||
       ((t = (double *)AlcCalloc(nest, sizeof(double))) == NULL) ||
       ((u = (double *)AlcCalloc(nV, sizeof(double))) == NULL) ||
       ((w = (double *)AlcMalloc(nV * sizeof(double))) == NULL) ||
       ((x = (double *)AlcMalloc(nV * dim * sizeof(double))) == NULL) || 
       ((iWrk = (int *)AlcCalloc(nest, sizeof(int))) == NULL) ||
       ((dWrk = (double *)AlcCalloc(lwrk, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int	 	i,
    		j = 0;
    double 	fp = 0.0;
    AlgError algErr = ALG_ERR_NONE;

    WlzValueSetDouble(w, 1.0, nV);
    switch(vType)
    {
      case WLZ_VERTEX_I2:
	{
	  WlzIVertex2 *v;

	  v = vtx.i2;
	  for(i = 0; i < nV; ++i)
	  {
	    x[j + 0] = v[i].vtX;
	    x[j + 1] = v[i].vtY;
	    j += 2;
	  }
	}
	break;
      case WLZ_VERTEX_L2:
	{
	  WlzLVertex2 *v;

	  v = vtx.l2;
	  for(i = 0; i < nV; ++i)
	  {
	    x[j + 0] = v[i].vtX;
	    x[j + 1] = v[i].vtY;
	    j += 2;
	  }
	}
	break;
      case WLZ_VERTEX_F2:
	{
	  WlzFVertex2 *v;

	  v = vtx.f2;
	  for(i = 0; i < nV; ++i)
	  {
	    x[j + 0] = v[i].vtX;
	    x[j + 1] = v[i].vtY;
	    j += 2;
	  }
	}
	break;
      case WLZ_VERTEX_D2:
	{
	  WlzDVertex2 *v;

	  v = vtx.d2;
	  for(i = 0; i < nV; ++i)
	  {
	    x[j + 0] = v[i].vtX;
	    x[j + 1] = v[i].vtY;
	    j += 2;
	  }
	}
	break;
	case WLZ_VERTEX_I3:
	  {
	    WlzIVertex3 *v;

	    v = vtx.i3;
	    for(i = 0, j = 0; i < nV; ++i)
	    {
	      x[j + 0] = v[i].vtX;
	      x[j + 1] = v[i].vtY;
	      x[j + 2] = v[i].vtZ;
	      j += 3;
	    }
	  }
	  break;
	case WLZ_VERTEX_L3:
	  {
	    WlzLVertex3 *v;

	    v = vtx.l3;
	    for(i = 0, j = 0; i < nV; ++i)
	    {
	      x[j + 0] = v[i].vtX;
	      x[j + 1] = v[i].vtY;
	      x[j + 2] = v[i].vtZ;
	      j += 3;
	    }
	  }
	  break;
	case WLZ_VERTEX_F3:
	  {
	    WlzFVertex3 *v;

	    v = vtx.f3;
	    for(i = 0, j = 0; i < nV; ++i)
	    {
	      x[j + 0] = v[i].vtX;
	      x[j + 1] = v[i].vtY;
	      x[j + 2] = v[i].vtZ;
	      j += 3;
	    }
	  }
	  break;
	case WLZ_VERTEX_D3:
	  {
	    WlzDVertex3 *v;

	    v = vtx.d3;
	    for(i = 0, j = 0; i < nV; ++i)
	    {
	      x[j + 0] = v[i].vtX;
	      x[j + 1] = v[i].vtY;
	      x[j + 2] = v[i].vtZ;
	      j += 3;
	    }
	  }
	  break;
	default:
	  break;
    }
    if(periodic)
    {
      algErr = ALG_ERR_UNIMPLEMENTED;
    }
    else
    {
      algErr = AlgBSplineNDFit(0, 0, dim, nV, u, dim * nV, x, w,
	  0.0, 0.0, k, sm, nest, &nT, t, &nC, c, &fp, dWrk, iWrk);
    }
    errNum = WlzErrorFromAlg(algErr);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bs = WlzMakeBSpline((dim == 2)? WLZ_BSPLINE_C2D: WLZ_BSPLINE_C3D,
         k, nT, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    bs->nKnots = nT;
    WlzValueCopyDoubleToDouble(bs->knots, t, nT);
    WlzValueCopyDoubleToDouble(bs->coefficients, c, dim * nT);
  }
  AlcFree(c);
  AlcFree(t);
  AlcFree(u);
  AlcFree(w);
  AlcFree(x);
  AlcFree(iWrk);
  AlcFree(dWrk);
  if(errNum != WLZ_ERR_NONE)
  {
    WlzDomain	dom;

    dom.bs = bs;
    (void )WlzFreeDomain(dom);
    bs = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(bs);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Evaluates a B-spline at a specified number of points.
* \param	bs		Given B-spline domain.
* \param	n		Number of points at which to evaluate the
* 				B-spline. If zero evaluations are at the
* 				knots.
* \param	x		If NULL then the evaluations will either
* 				be at the knots (n == 0) or at n equaly
* 				spaced points along the parametrised curve.
* 				If not NULL then x must be a pointer to n
* 				points along the parametrised curve (with
* 				a single number per points).
* \param	deriv		Order of derivative, range [0-
* 				(WLZ_BSPLINE_ORDER_MAX - 1)].
* \param	eval		An array of n WlzDVertex2 or WlzDVertex3
* 				vertices for the evaluation.
*/
WlzErrorNum			WlzBSplineEval(
				  WlzBSpline *bs,
				  int n,
				  double *x,
				  int deriv,
				  WlzVertexP eval)
{
  double	*buf = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bs == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((n < 0) || (deriv < 0) || (deriv >= WLZ_BSPLINE_ORDER_MAX))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(eval.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    
    int		nb,
    		nn;

    nb = (deriv == 0)? 0: bs->nKnots;
    nn = 2 * ((n == 0)? bs->nKnots: n);
    if((buf = (double *)AlcMalloc((nn + nb) * sizeof(double))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int 	i,
    		dim;
    double	*y,
    		*w;
    
    if(n == 0)
    {
      y = bs->knots + bs->order;
      n = bs->nKnots - (2 * bs->order);
      WlzValueCopyDoubleToDouble(buf, y, n);
    }
    else
    {
      if(x)
      {
        WlzValueCopyDoubleToDouble(buf, x, n);
      }
      else
      {
	double	n1;

	n1 = n - 1.0;
	for(i = 0; i < n; ++i)
	{
	  buf[i] = (double )i / n1;
	}
      }
    }
    x = buf;
    y = buf + n;
    dim = (bs->type == WLZ_BSPLINE_C2D)? 2: 3;
    for(i = 0; i < dim; ++i)
    {
      double	*coeff;
      AlgError	algErr;

      coeff = bs->coefficients + (i * bs->nKnots);
      if(deriv == 0)
      {
	algErr = AlgBSplineEval(bs->knots, bs->nKnots, coeff, bs->order,
	    x, y, n);
      }
      else
      {
        w = y + n;
        algErr = AlgBSplineDer(bs->knots, bs->nKnots, coeff, bs->order,
	    deriv, x, y, n, w);
      }
      if((errNum = WlzErrorFromAlg(algErr)) != WLZ_ERR_NONE)
      {
        break;
      }
      if(dim == 2)
      {
        int	j;
	WlzDVertex2 *ep;

	ep = eval.d2;
	if(i == 0)
	{
	  for(j = 0; j < n; ++j)
	  {
	    ep[j].vtX = y[j];
	  }
	}
	else
	{
	  for(j = 0; j < n; ++j)
	  {
	    ep[j].vtY = y[j];
	  }
	}
      }
      else /* dim == 3 */
      {
        int	j;
	WlzDVertex3 *ep;

	ep = eval.d3;
	switch(i)
	{
	  case 0:
	    for(j = 0; j < n; ++j)
	    {
	      ep[j].vtX = y[j];
	    }
	    break;
	  case 1:
	    for(j = 0; j < n; ++j)
	    {
	      ep[j].vtY = y[j];
	    }
	    break;
	  case 2:
	    for(j = 0; j < n; ++j)
	    {
	      ep[j].vtZ = y[j];
	    }
	    break;
	  default:
	    break;
	}
      }
    }
  }
  AlcFree(buf);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Evaluates a B-spline at a single specified point.
* \todo		Currently WlzBSplineEvalSP() is just a wrapper for
* 		WlzBSplineEval() however this is quite inefficient and
* 		may change.
* \param	bs		Given B-spline domain.
* \param	x		Parametric coordinate of the point.
* \param	deriv		Order of derivative, range [0-
* 				(WLZ_BSPLINE_ORDER_MAX - 1)].
* \param	eval		Destination pointer for a single WlzDVertex2
* 				or WlzDVertex3 vertex following evaluation.
*/
WlzErrorNum			WlzBSplineEvalSP(
				  WlzBSpline *bs,
				  double x,
				  int deriv,
				  WlzVertexP eval)
{
  WlzErrorNum	errNum;

  errNum = WlzBSplineEval(bs, 1, &x, deriv, eval);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzFeatures
* \brief	Compute the (optional) position and (optional) tangent to
* 		the B-spline at points along the path of it's curve.
* \param	bs		Given B-spline domain.
* \param	n		Number of parametric coordinates given.
* \param	x		Given parametric coordinates at points along
* 				the B-spline curve.
* \param	dstPos		Destination pointer for the n positions of
* 				the points, may be NULL.
* \param	dstTnt		Destination pointer for the n unit tangent
* 				vector at the points, may be NULL.
*/
WlzErrorNum			WlzBSplineTangent(
				  WlzBSpline *bs,
				  int n,
				  double *x,
				  WlzVertexP dstPos,
				  WlzVertexP dstTnt)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bs == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    if(dstPos.v)
    {
      errNum = WlzBSplineEval(bs, n, x, 0, dstPos);
    }
    if((errNum == WLZ_ERR_NONE) && dstTnt.v)
    {
      errNum = WlzBSplineEval(bs, n, x, 1, dstTnt);
      if(errNum == WLZ_ERR_NONE)
      {
	int	i,
		dim;

	dim = (bs->type == WLZ_BSPLINE_C2D)? 2: 3;
	if(dim == 2)
	{
	  for(i = 0; i < n; ++i)
	  {
	    double s;
	    WlzDVertex2 *p;

	    p = dstTnt.d2 + i;
	    s = WLZ_VTX_2_LENGTH(*p);
	    if(s < DBL_EPSILON)
	    {
	      errNum = WLZ_ERR_DOUBLE_DATA;
	      break;
	    }
	    s = 1.0 / s;
	    WLZ_VTX_2_SCALE(*p, *p, s);
	  }
	}
	else /* dim == 3 */
	{
	  for(i = 0; i < n; ++i)
	  {
	    double s;
	    WlzDVertex3 *p;

	    p = dstTnt.d3 + i;
	    s = WLZ_VTX_3_LENGTH(*p);
	    if(s < DBL_EPSILON)
	    {
	      errNum = WLZ_ERR_DOUBLE_DATA;
	      break;
	    }
	    s = 1.0 / s;
	    WLZ_VTX_3_SCALE(*p, *p, s);
	  }
	}
      }
    }
  }
  return(errNum);
}
