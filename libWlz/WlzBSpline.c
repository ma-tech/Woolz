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

static WlzObject               	*WlzBSplineCutOtg(
                                  WlzObject *iObj,
                                  WlzBSpline *bs,
                                  int noGrey,
                                  int radius,
                                  double tB,
                                  double tE,
				  WlzInterpolationType interp,
                                  WlzErrorNum *dstErr);
static WlzObject               	*WlzBSplineCutPar(
                                  WlzObject *iObj,
                                  WlzBSpline *bs,
                                  int noGrey,
                                  int radius,
                                  double tB,
                                  double tE,
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
WlzBSpline			*WlzBSplineFromVertices(
				  WlzVertexType vType,
				  int nV,
				  WlzVertexP vtx,
				  int k,
				  int periodic,
				  double sm,
				  WlzErrorNum *dstErr)
{
  int		nest,
		nC = 0,
		nT = 0,
		dim = 0;
  WlzBSpline	*bs = NULL;
  int		*iWrk = NULL;
  double	*c = NULL,
		*t = NULL,
		*u = NULL,
		*w = NULL,
		*x = NULL,
		*dWrk = NULL;
		
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((nV <= 0) || (vtx.v == NULL) ||
     (k < WLZ_BSPLINE_ORDER_MIN) || (k > WLZ_BSPLINE_ORDER_MAX))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    dim = WlzVertexDim(vType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		lwrk;

    nest = nV + 2 * (k + 1);
    nC = dim * nest;
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
* \return	New spatial domain object or NULL on error.
* \ingroup	WlzFeatures
* \brief	Makes a new interval or plane domain object corresponding
* 		to the evaluation of the given B-spline in the given
* 		parametric coordinate range. The returned object may have
* 		pixel/voxels which are disconnected unless dilated by 1.
* \param	bs			Given B-spline domain.
* \param	tB			Begining of the parametric range.
* \param	TE			End of the parametric range.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzBSplineToDomain(
				  WlzBSpline *bs,
				  double tB,
				  double tE,
				  WlzErrorNum *dstErr)
{
  int		dim = 0,
  		len = 0;
  size_t	vsz = 0;
  double	*buf = NULL;
  WlzVertexP	pos = {0};
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(bs == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(tB > tE)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(bs->type)
    {
      case WLZ_BSPLINE_C2D:
        dim = 2;
	vsz = sizeof(WlzDVertex2);
	break;
      case WLZ_BSPLINE_C3D:
        dim = 3;
	vsz = sizeof(WlzDVertex3);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    len = (int )ceil(WlzBSplineLength(bs, tB, tE, &errNum));
  }
  /* Allocate buffers for evaluating the spline. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((buf = (double *)AlcMalloc(sizeof(double) * len)) == NULL) ||
       ((pos.v = AlcMalloc(vsz * len)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute point positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    double	t = 0.0;

    if(len > 1)
    {
      t = (tE - tB) / (len - 1);
    }
    for(i = 0; i < len; ++i)
    {
      buf[i] = tB + (i * t);
    }
    errNum = WlzBSplineEval(bs, len, buf, 0, pos);
  }
  AlcFree(buf);
  /* Compute bounding box and create domain object. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    WlzPixelV	bgd = {0};

    bgd.type = WLZ_GREY_UBYTE;
    if(dim == 2)
    {
      WlzDVertex2 *p;
      WlzDBox2	b;

      p = pos.d2;
      b.xMin = b.xMax = p[0].vtX;
      b.yMin = b.yMax = p[0].vtY;
      for(i = 1; i < len; ++i)
      {
	if(p[i].vtX < b.xMin)
	{
	  b.xMin = p[i].vtX;
	}
        else if(p[i].vtX > b.xMax)
	{
	  b.xMax = p[i].vtX;
	}
	if(p[i].vtY < b.yMin)
	{
	  b.yMin = p[i].vtY;
	}
        else if(p[i].vtY > b.yMax)
	{
	  b.yMax = p[i].vtY;
	}
      }
      obj = WlzAssignObject(
      	    WlzMakeRect((int )floor(b.yMin), (int )ceil(b.yMax),
			(int )floor(b.xMin), (int )ceil(b.xMax),
			bgd.type, NULL, bgd, NULL, NULL, &errNum), NULL);
    }
    else /* dim == 3 */
    {
      WlzDVertex3 *p;
      WlzDBox3	b;

      p = pos.d3;
      b.xMin = b.xMax = p[0].vtX;
      b.yMin = b.yMax = p[0].vtY;
      b.zMin = b.zMax = p[0].vtZ;
      for(i = 1; i < len; ++i)
      {
	if(p[i].vtX < b.xMin)
	{
	  b.xMin = p[i].vtX;
	}
        else if(p[i].vtX > b.xMax)
	{
	  b.xMax = p[i].vtX;
	}
	if(p[i].vtY < b.yMin)
	{
	  b.yMin = p[i].vtY;
	}
        else if(p[i].vtY > b.yMax)
	{
	  b.yMax = p[i].vtY;
	}
	if(p[i].vtZ < b.zMin)
	{
	  b.zMin = p[i].vtZ;
	}
        else if(p[i].vtZ > b.zMax)
	{
	  b.zMax = p[i].vtZ;
	}
      }
      obj = WlzAssignObject(
            WlzMakeCuboid((int )floor(b.zMin), (int )ceil(b.zMax),
			  (int )floor(b.yMin), (int )ceil(b.yMax),
			  (int )floor(b.xMin), (int )ceil(b.xMax),
			  bgd.type, bgd, NULL, NULL, &errNum), NULL);
    }
  }
  /* Set pixels/voxels in the domain object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzGreyValueWSpace *gVWSp;

    gVWSp = WlzGreyValueMakeWSp(obj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      int	i;
      WlzPixelV	thr;
      WlzObject	*tObj = NULL;

      thr.v.ubv = 1;
      thr.type = WLZ_GREY_UBYTE;
      if(dim == 2)
      {
        WlzDVertex2 *p;

	p = pos.d2;
        for(i = 0; i < len; ++i)
	{
	  WlzGreyValueGet(gVWSp, 0, p[i].vtY, p[i].vtX);
	  *(gVWSp->gPtr[0].ubp) = 1;
	}
      }
      else /* dim == 3 */
      {
        WlzDVertex3 *p;

	p = pos.d3;
        for(i = 0; i < len; ++i)
	{
	  WlzGreyValueGet(gVWSp, p[i].vtZ, p[i].vtY, p[i].vtX);
	  *(gVWSp->gPtr[0].ubp) = 1;
	}
      }
      WlzGreyValueFreeWSp(gVWSp);
      tObj = WlzAssignObject(
             WlzThreshold(obj, thr, WLZ_THRESH_HIGH, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	WlzValues val = {0};

        (void )WlzFreeObj(obj);
	obj = WlzMakeMain(tObj->type, tObj->domain, val, NULL, NULL, &errNum);
      }
      (void )WlzFreeObj(tObj);
    }
  }
  AlcFree(pos.v);
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
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

/*!
* \return	Woolz object cut from the input object.
* \ingroup	WlzFeatures
* \brief	Cuts regions from a spatial domain using a B-spline to define
* 		the region cut. The region may be the domain of the B-spline
* 		or planes orthogonal to the B-spline.
* \param	iObj			Object to be cut.
* \param	bs			B-spline domain.
* \param	cutOrthog		If non-zero cut orthogonal planes
* 					rather than the B-spline dilated by
* 					given radius. Cutting othogonal
* 					planes is only allowed for objects
* 					and B-splines in 3D.
* \param	noGrey			If non-zero then don't fill the
* 					returned objects grey values (if they
* 					exist).
* \param	radius			Radius of region cut, the orthogonal
* 					distance from the B-spline. If
* 					cutting the region of the dilated
* 					B-spline then this is it's dilation,
* 					otherwise if cutting orthogonal planes
* 					this is the radius of the sections.
* 					If zero then either just the
* 					pixels/voxels intersecting the B-spline
* 					or the entire (unclipped) orthogonal
* 					planes will be cut.
* \param	tB			B-spline parametric coordinate at which
* 					to start the cut.
* \param	tE			B-spline parametric coordinate at which
* 					to end the cut.
* \param	interp			Interpolation for cutting sections.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject                	*WlzBSplineCut(
                                  WlzObject *iObj,
                                  WlzBSpline *bs,
                                  int cutOrthog,
                                  int noGrey,
                                  int radius,
                                  double tB,
                                  double tE,
				  WlzInterpolationType interp,
                                  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(iObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((iObj->domain.core == NULL) || (bs == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((radius < 0) || (tB < 0.0) || (tE > 1.0) || (tB > tE))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if(iObj->values.core == NULL)
    {
      noGrey = 1;
    }
    else if(iObj->values.core->type != WLZ_VOXELVALUETABLE_GREY)
    {
      errNum = WLZ_ERR_VALUES_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(iObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if((bs->type != WLZ_BSPLINE_C2D) || cutOrthog)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	if(bs->type != WLZ_BSPLINE_C3D)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
        break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = (cutOrthog)?
        WlzBSplineCutOtg(iObj, bs, noGrey, radius, tB, tE, interp, &errNum):
	WlzBSplineCutPar(iObj, bs, noGrey, radius, tB, tE, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Length along the spline's parametric curve.
* \ingroup	WlzFeatures
* \brief	Computes the length of a path along a spline's curve between
* 		a pair of parametric coordinates.
* 		Given a parametric curve \f$c(t)\f$ then the length \f$l\f$
* 		of a segment on that curve between \f$t_a\f$ and \f$t_b\f$
* 		is given by
* 		\f[
* 		l = \int_{t_a}^{t_b}\sqrt{
  		     {\acute{x}}^2 +{\acute{y}}^2 + \ldots}dt
 		\f]
* 		where \f$\acute{x}, \acute{y}, \ldots\f$ are the first
* 		derivatives of coordinates \f$x, y, \ldots\f$ with respect
* 		to \f$t\f$.
* 		Because the integral is an elliptic integral, then Legendre-
* 		Gauss quadrature is used for numerical integration with
* 		each spline segment integrated using it's own points and
* 		weights.
* \param	bs			Given B-spline domain.
* \param	tB			Parametric coordinate of start.
* \param	tE			Parametric coordinate of end.
* \param	dstErr			Destination error pointer, may be NULL.
*/
double				WlzBSplineLength(
				  WlzBSpline *bs,
				  double tB,
				  double tE,
				  WlzErrorNum *dstErr)
{
  int		n,
		iB = -1,
		iE = -1,
  		dim = 0,
		epk = 0; 		    /* Evaluations per knot interval. */
  double	len = 0.0;
  size_t	vsz = 0;
  double	*x = NULL;
  WlzVertexP	vtx = {0};
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(bs == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((tB < 0.0) || (tE > 1.0) || (tE < tB))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    switch(bs->type)
    {
      case WLZ_BSPLINE_C2D:
        dim = 2;
	vsz = sizeof(WlzDVertex2);
	break;
      case WLZ_BSPLINE_C3D:
        dim = 3;
	vsz = sizeof(WlzDVertex3);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  /* Find knot indices (iB and iE) for which the knots enclose the begining
   * and end parametric coordinates (tB and tE). */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    epk = bs->order + 1;
    /* Find \f$ iB = i st t_i \leq tB < t_{i+1} \f$ */
    for(i = 0; i < bs->nKnots; ++i)
    {
      if(bs->knots[i] >= tB)
      {
	iB = (i == 0)? 0: i - 1;
        break;
      }
    }
    /* Find \f$ iE = i st t_{i-1} < tE \leq t_i \f$ */
    for(i = bs->nKnots - 1; i > iB; --i)
    {
      if(bs->knots[i] <= tE)
      {
	iE = (i < bs->nKnots - 1)? i + 1: bs->nKnots - 1;
        break;
      }
    }
    if(iB >= iE)
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Allocate room for the parametric coordinates and vertices used to
   * compute 1st derivatives. */
  if(errNum == WLZ_ERR_NONE)
  {
    n = (iE - iB + 1) * epk;
    if(((vtx.v = AlcMalloc(n * vsz)) == NULL) ||
       ((x = (double *)AlcMalloc(n * sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute parametic coordinates which cover the knot intervals
   * and then compute the partial derivatives at these coordinates. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
    		j = 0;

    for(i = iB; i < iE; ++i)
    {
      int	k;
      double	t0,
		t1,
		ta,
		ts;

      t0 = (i     == iB)? tB: bs->knots[i];
      t1 = (i + 1 == iE)? tE: bs->knots[i + 1];
      ta = t1 + t0;
      ts = t1 - t0;
      for(k = 0; k < epk; ++k)
      {
	double	p;

	p = AlgGaussLegendrePoints(epk, k);
        x[j] = (ts * p) + ta;
	++j;
      }
    }
    n = j;
    errNum = WlzBSplineEval(bs, n, x, 1, vtx);
  }
  /* Compute the integral of the modulus of the partial derivatives using
   * the Gauss-Legendre weights. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
    		j = 0;

    for(i = iB; i < iE; ++i)
    {
      int	k;
      double	l,
      		t0,
		t1,
		ts;

      l = 0;
      t0 = (i     == iB)? tB: bs->knots[i];
      t1 = (i + 1 == iE)? tE: bs->knots[i + 1];
      ts = t1 - t0;
      for(k = 0; k < epk; ++k)
      {
	double	w;

	w = AlgGaussLegendreWeights(epk, k);
	l += w * ((dim == 2)? WLZ_VTX_2_LENGTH(vtx.d2[j]):
	                      WLZ_VTX_3_LENGTH(vtx.d3[j]));
	++j;
      }
      len += ts * l;
    }
    len *= 0.5;
  }
  AlcFree(x);
  AlcFree(vtx.v);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(len);
}

/*!
* \return	Woolz object cut from the input object.
* \ingroup	WlzFeatures
* \brief	Cuts orthogonal planes from a spatial domain using a B-spline
* 		to define orthogonality of the cutting plane cut.
* 		This may only be called for objects and splines in 3D but
* 		this is not tested for because this is a static function and
* 		the calling function tests for this.
* \param	iObj			Object to be cut.
* \param	bs			B-spline domain.
* \param	noGrey			If non-zero then don't fill the
* 					returned objects grey values (if they
* 					exist).
* \param	radius			Radius of cut domain (centred on the
* 					B-spline). If zero then entire
* 					(unclipped) orthogonal planes will be
* 					cut.
* \param	tB			B-spline parametric coordinate at which
* 					to start the cut.
* \param	tE			B-spline parametric coordinate at which
* 					to end the cut.
* \param	interp			Interpolation for cutting sections.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject               	*WlzBSplineCutOtg(
                                  WlzObject *iObj,
                                  WlzBSpline *bs,
                                  int noGrey,
                                  int radius,
                                  double tB,
                                  double tE,
				  WlzInterpolationType interp,
                                  WlzErrorNum *dstErr)
{
  int		len = 1;
  double	*buf = NULL;
  WlzVertexP	nrm,
  		pos = {0};
  WlzObject	*sObj = NULL,
  		*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-06;

  /* Compute buffer sizes. */
  if(tB + eps < tE)
  {
    len = (int )ceil(WlzBSplineLength(bs, tB, tE, &errNum)) + 1;
  }
  /* Create subdomain in 2D of given radius. */
  if((errNum == WLZ_ERR_NONE) && (radius > 0))
  {
   sObj = WlzAssignObject(
	  WlzMakeSphereObject(WLZ_2D_DOMAINOBJ, radius, 0.0, 0.0, 0.0,
			      &errNum), NULL);
  }
  /* Allocate buffers for evaluating the spline and it's derivative. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((buf = (double *)AlcMalloc(sizeof(double) * len)) == NULL) ||
       ((pos.v = AlcMalloc(sizeof(WlzDVertex3) * 2 * len)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Evaluate the spline and it's derivative. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    double	t = 0.0;

    nrm.d3 = pos.d3 + len;
    if(len > 1)
    {
      t = (tE - tB) / (len - 1);
    }
    for(i = 0; i < len; ++i)
    {
      buf[i] = tB + (i * t);
    }
    errNum = WlzBSplineEval(bs, len, buf, 0, pos);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzBSplineEval(bs, len, buf, 1, nrm);
    }
  }
  AlcFree(buf);
  /* Create 3D domain object with len planes. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	rDom = {0};
    WlzValues	rVal = {0};

    rDom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN, 0, len - 1,
                                0, 1, 0, 1, &errNum);
    if((errNum == WLZ_ERR_NONE) && !noGrey)
    {
      rVal.vox = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY, 0, len - 1,
				     iObj->values.vox->bckgrnd, NULL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, rDom, rVal, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreePlaneDomain(rDom.p);
      (void )WlzFreeVoxelValueTb(rVal.vox);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		p;
    WlzDVertex3 z = {0};

    for(p = 0; p < len; ++p)
    {
      WlzObject	*obj2 = NULL;
      WlzThreeDViewStruct *vs;

      vs = Wlz3DViewStructFromNormal(nrm.d3[p], pos.d3[p], z, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        errNum = WlzInit3DViewStruct(vs, iObj);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	obj2 = WlzGetSubSectionFromObject(iObj, sObj, vs, interp,
		                          NULL, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
        rObj->domain.p->domains[p] = WlzAssignDomain(obj2->domain, NULL);
	if(rObj->values.core)
	{
	  rObj->values.vox->values[p] = WlzAssignValues(obj2->values, NULL);
	}
      }
      (void )WlzFree3DViewStruct(vs);
      (void )WlzFreeObj(obj2);
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  (void )WlzFreeObj(sObj);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzStandardPlaneDomain(rObj->domain.p, rObj->values.vox);
  }
  AlcFree(pos.v);
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeObj(rObj);
    rObj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Woolz object cut from the input object.
* \ingroup	WlzFeatures
* \brief	Cuts a spatial domain parallel to a B-spline within the given
* 		object.
* \param	iObj			Object to be cut.
* \param	bs			B-spline domain.
* \param	noGrey			If non-zero then don't fill the
* 					returned objects grey values (if they
* 					exist).
* \param	radius			Radius of cut domain (centred on the
* 					B-spline). If zero then entire
* 					(unclipped) orthogonal planes will be
* 					cut.
* \param	tB			B-spline parametric coordinate at which
* 					to start the cut.
* \param	tE			B-spline parametric coordinate at which
* 					to end the cut.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject               	*WlzBSplineCutPar(
                                  WlzObject *iObj,
                                  WlzBSpline *bs,
                                  int noGrey,
                                  int radius,
                                  double tB,
                                  double tE,
                                  WlzErrorNum *dstErr)
{
  int		dim;
  WlzObject	*o1,
  		*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  dim = (iObj->type == WLZ_2D_DOMAINOBJ)? 2: 3;
  /* Get domain corresponding to the B-spline and dilate if required. */
  o1 = WlzAssignObject(
       WlzBSplineToDomain(bs, tB, tE, &errNum), NULL);
  if((errNum == WLZ_ERR_NONE) && (radius > 0))
  {
    WlzObject	*o2 = NULL;

    if(radius > 1)
    {
      WlzObject	*os = NULL;

      os = WlzAssignObject(
	   WlzMakeSphereObject(iObj->type, radius, 0.0, 0.0, 0.0,
			       &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
        o2 = WlzAssignObject(
	     WlzStructDilation(o1, os, &errNum), NULL);
      }
      (void )WlzFreeObj(os);
    }
    else
    {
      o2 = WlzAssignObject(
           WlzDilation(o1,
	               (dim == 2)? WLZ_8_CONNECTED: WLZ_26_CONNECTED,
		       &errNum), NULL);
    }
    (void )WlzFreeObj(o1);
    o1 = o2;
  }
  /* Take intersection of the domain with the given object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject   *o2;
    WlzObject	*objs[2];

    objs[0] = iObj;
    objs[1] = o1;
    o2 = WlzAssignObject(                 
         WlzIntersectN(2, objs, 0, &errNum), NULL);
    (void )WlzFreeObj(o1);
    o1 = o2;
  }
  /* Transfer grey values if required. */
  if((errNum == WLZ_ERR_NONE) && !noGrey)
  {
    WlzObject   *o2 = NULL;      
    
    o2 = WlzAssignObject(
         WlzGreyTransfer(o1, iObj, 0, &errNum), NULL);
    (void )WlzFreeObj(o1);       
    o1 = o2;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(o1->type, o1->domain, o1->values, NULL, NULL, &errNum);
  }
  (void )WlzFreeObj(o1);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

