#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvexHullClarkson_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzConvexHullClarkson.c
* \author       Bill Hill
* \date         July 2011
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
* \brief	Functions to compute convex hulls using the Clarkson's
* 		algorithm.
* \ingroup	WlzConvexHull
*/

#include <stdlib.h>
#include <stdio.h>
#include <Wlz.h>

static int			WlzConvHullClarksonCCWI(
				  WlzIVertex2 **v,
				  int i,
				  int j,
				  int k);
static int			WlzConvHullClarksonCCWD(
				  WlzDVertex2 **v,
				  int i,
				  int j,
				  int k);
static int			WlzConvHullClarksonMakeChain2I(
				  WlzIVertex2 **v,
				  int n,
                                  int (*cmp)(const void *, const void *));
static int			WlzConvHullClarksonMakeChain2D(
				  WlzDVertex2 **v,
				  int n,
                                  int (*cmp)(const void *, const void *));
static int			WlzConvHullClarksonCmpIL(
				  const void *a,
				  const void *b);
static int			WlzConvHullClarksonCmpDL(
				  const void *a,
				  const void *b);
static int			WlzConvHullClarksonCmpIH(
				  const void *a,
				  const void *b);
static int			WlzConvHullClarksonCmpDH(
				  const void *a,
				  const void *b);

/*!
* \def		WLZ_CONVHULL_CLARKSON_SM_2D
* \brief	Used to avoid dynamic allocation for small vertex arrays.
* 		The stack will be used for internal workspace arrays smaller
* 		than this.
*/
#define		WLZ_CONVHULL_CLARKSON_SM_2D	(16)

/*!
* \return       New 2D convex hull domain.
* \ingroup      WlzConvexHull
* \brief        Creates a new 2D convex hull domain which encloses the
*               given vertices using Clarkson's algorithm, see
		WlzConvHullClarkson2I() and WlzConvHullClarkson2D().
* \param        pType                   Type of vertex given, must be either
*                                       WLZ_VERTEX_I3 or WLZ_VERTEX_D3.
* \param        nPnt                    Number of given vertices.
* \param        pnt                     The given vertices.
* \param        dstErr                  Destination error pointer, may be NULL.
*                                       If the volume of the tetrahedron
*                                       with maximum/minimum z coordinate
*                                       and minimum x and y coordiantes is
*                                       zero the erro code will be
*                                       WLZ_ERR_DEGENERATE.
*/
WlzConvHullDomain2              *WlzConvexHullFromVtx2(
                                  WlzVertexType pType,
                                  int nPnt,
                                  WlzVertexP pnt,
                                  WlzErrorNum *dstErr)
{
  int		nC = 0;
  int		*idx = NULL;
  WlzConvHullDomain2 *cvh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(nPnt < 3)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if(pnt.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(pType)
    {
      case WLZ_VERTEX_I2:
	nC = WlzConvHullClarkson2I(pnt.i2, nPnt, &idx, &errNum);
	break;
      case WLZ_VERTEX_D2:
	nC = WlzConvHullClarkson2D(pnt.d2, nPnt, &idx, &errNum);
	break;
      default:
	errNum = WLZ_ERR_PARAM_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cvh = WlzMakeConvexHullDomain2(nC, pType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    WlzDVertex2 c;

    WLZ_VTX_2_ZERO(c);
    if(pType == WLZ_VERTEX_I2)
    {
      for(i = 0; i < nC; ++i)
      {
	WlzIVertex2 p;

	p = pnt.i2[idx[i]];
	WLZ_VTX_2_ADD(c, c, p);
	cvh->vertices.i2[i] = p;
      }
      WLZ_VTX_2_SCALE(c, c, (1.0 / nC));
      WLZ_VTX_2_NINT(cvh->centroid.i2, c);
    }
    else /* pType == WLZ_VERTEX_D2 */
    {
      for(i = 0; i < nC; ++i)
      {
	WlzDVertex2 p;

	p = pnt.d2[idx[i]];
	WLZ_VTX_2_ADD(c, c, p);
	cvh->vertices.d2[i] = p;
      }
      WLZ_VTX_2_SCALE(c, c, (1.0 / nC));
      cvh->centroid.d2 = c;
    }
    cvh->nVertices = nC;
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreeConvexHullDomain2(cvh);
  }
  AlcFree(idx);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}

/*!
* \return	Number of vertices in the convex hull or zero on error.
* \ingroup	WlzConvexHull
* \brief	Computes the convex hull of a given array of 2D integer
* 		vertices which is returned as an array of indices into
* 		the array of vertices.
* 		This index array should be freed using AlcFree().
* 		The vertex array is not changed.
* \param	vtx			Given array of vertices.
* \param	n			Number of vertices in the given array.
* \param	dstIdx			Destination pointer for array of
* 					indices to the convex hull vertices.
* 					The convex hull vertices are ordered
* 					counter-clockwise.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzConvHullClarkson2I(WlzIVertex2 *vtx, int n, int **dstIdx,
				      WlzErrorNum *dstErr)
{
  int		i,
  		u = 0;
  int		*idx = NULL;
  WlzIVertex2	**v = NULL;
  WlzIVertex2	*vSmall[WLZ_CONVHULL_CLARKSON_SM_2D];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(n < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((vtx == NULL) || (dstIdx == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if(n <= WLZ_CONVHULL_CLARKSON_SM_2D)
    {
      v = vSmall; 
    }
    else if((v = (WlzIVertex2 **)
		 AlcMalloc(sizeof(WlzIVertex2 *) * (n + 1))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(i = 0; i < n; ++i)
    {
      v[i] = &(vtx[i]);
    }
    u = WlzConvHullClarksonMakeChain2I(v, n, WlzConvHullClarksonCmpIL);
    v[n] = v[0];
    u += WlzConvHullClarksonMakeChain2I(v + u, n - u + 1,
				      WlzConvHullClarksonCmpIH);
    if((idx = (int *)AlcMalloc(sizeof(int) * u)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(idx != NULL)
  {
    for(i = 0; i < u; ++i) 
    {
      idx[i] = v[i] - &(vtx[0]);
    }
  }
  if(n > WLZ_CONVHULL_CLARKSON_SM_2D)
  {
    AlcFree(v);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstIdx = idx;
  }
  else
  {
    u = 0;
  }
  return(u);
}

/*!
* \return	Number of vertices in the convex hull or zero on error.
* \ingroup	WlzConvexHull
* \brief	Computes the convex hull of a given array of 2D double vertices
* 		which is returned as an array of indices into the array of
* 		vertices.
* 		This index array should be freed using AlcFree().
* 		The vertex array is not changed.
* \param	vtx			Given array of vertices.
* \param	n			Number of vertices in the given array.
* \param	dstIdx			Destination pointer for array of
* 					indices to the convex hull vertices.
* 					The convex hull vertices are ordered
* 					counter-clockwise.
* \param	dstErr			Destination error pointer, may be NULL.
*/
int		WlzConvHullClarkson2D(WlzDVertex2 *vtx, int n, int **dstIdx,
				      WlzErrorNum *dstErr)
{
  int		i,
  		u = 0;
  int		*idx = NULL;
  WlzDVertex2	**v = NULL;
  WlzDVertex2	*vSmall[WLZ_CONVHULL_CLARKSON_SM_2D];
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(n < 1)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else if((vtx == NULL) || (dstIdx == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if(n <= WLZ_CONVHULL_CLARKSON_SM_2D)
    {
      v = vSmall; 
    }
    else if((v = (WlzDVertex2 **)
		 AlcMalloc(sizeof(WlzDVertex2 *) * (n + 1))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(i = 0; i < n; ++i)
    {
      v[i] = &(vtx[i]);
    }
    u = WlzConvHullClarksonMakeChain2D(v, n, WlzConvHullClarksonCmpDL);
    v[n] = v[0];
    u += WlzConvHullClarksonMakeChain2D(v + u, n - u + 1,
				      WlzConvHullClarksonCmpDH);
    if((idx = (int *)AlcMalloc(sizeof(int) * u)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(idx != NULL)
  {
    for(i = 0; i < u; ++i) 
    {
      idx[i] = v[i] - &(vtx[0]);
    }
  }
  if(n > WLZ_CONVHULL_CLARKSON_SM_2D)
  {
    AlcFree(v);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstIdx = idx;
  }
  else
  {
    u = 0;
  }
  return(u);
}

static int	WlzConvHullClarksonMakeChain2I(WlzIVertex2 **v, int n,
                                      int (*cmp)(const void*, const void*))
{
  int		i,
  		j,
		s = 1;
  WlzIVertex2 	*t;

  qsort(v, n, sizeof(WlzIVertex2 *), cmp);
  for(i=2; i < n; ++i)
  {
    j = s;
    while((j >= 1) && WlzConvHullClarksonCCWI(v, i, j, j - 1))
    {
      --j;
    }
    s = j + 1;
    t = v[s]; v[s] = v[i]; v[i] = t;
  }
  return(s);
}

static int	WlzConvHullClarksonMakeChain2D(WlzDVertex2 **v, int n,
                                      int (*cmp)(const void*, const void*))
{
  int		i,
  		j,
		s = 1;
  WlzDVertex2 	*t;

  qsort(v, n, sizeof(WlzDVertex2 *), cmp);
  for(i=2; i < n; ++i)
  {
    j = s;
    while((j >= 1) && WlzConvHullClarksonCCWD(v, i, j, j - 1))
    {
      --j;
    }
    s = j + 1;
    t = v[s]; v[s] = v[i]; v[i] = t;
  }
  return(s);
}

/*!
* \return	Returns non-zero if the vertices indexed by i, j and k are
* 		counter-clockwise.
* \ingroup	WlzConvexHull
* \brief	Tests for counter-clockwise order of the vertices indexed
* 		by i, j and k in the given array of vertex pointers.
* \param	v			Given array of vertex pointers.
* \param	i			Index i into the array of vertex
* 					pointers.
* \param	j			Index j into the array of vertex
* 					pointers.
* \param	k			Index k into the array of vertex
* 					pointers.
*/
static int	WlzConvHullClarksonCCWI(WlzIVertex2 **v, int i, int j, int k)
{
  int		ccw;
  WlzLong	a,
  		b,
		c,
		d;

  a = v[i]->vtX - v[j]->vtX,
  b = v[i]->vtY - v[j]->vtY,
  c = v[k]->vtX - v[j]->vtX,
  d = v[k]->vtY - v[j]->vtY;
  ccw = (a * d) <= (b * c);
  return(ccw);
}

/*!
* \return	Returns non-zero if the vertices indexed by i, j and k are
* 		counter-clockwise.
* \ingroup	WlzConvexHull
* \brief	Tests for counter-clockwise order of the vertices indexed
* 		by i, j and k in the given array of vertex pointers.
* \param	v			Given array of vertex pointers.
* \param	i			Index i into the array of vertex
* 					pointers.
* \param	j			Index j into the array of vertex
* 					pointers.
* \param	k			Index k into the array of vertex
* 					pointers.
*/
static int	WlzConvHullClarksonCCWD(WlzDVertex2 **v, int i, int j, int k)
{
  int		ccw;
  double	a,
  		b,
		c,
		d;

  a = v[i]->vtX - v[j]->vtX,
  b = v[i]->vtY - v[j]->vtY,
  c = v[k]->vtX - v[j]->vtX,
  d = v[k]->vtY - v[j]->vtY;
  ccw = (a * d) - (b * c) <= 0.0;
  return(ccw);
}

static int	WlzConvHullClarksonCmpIL(const void *a, const void *b)
{
  int		cmp = 0; 

  if(((*(WlzIVertex2 **)a)->vtX - (*(WlzIVertex2 **)b)->vtX) > 0)
  {
    cmp = 1;
  }
  else if(((*(WlzIVertex2 **)b)->vtY - (*(WlzIVertex2 **)a)->vtY) < 0)
  {
    cmp = -1;
  }
  return(cmp);
}

static int	WlzConvHullClarksonCmpDL(const void *a, const void *b)
{
  int		cmp = 0; 

  if(((*(WlzDVertex2 **)a)->vtX - (*(WlzDVertex2 **)b)->vtX) > 0.0)
  {
    cmp = 1;
  }
  else if(((*(WlzDVertex2 **)b)->vtY - (*(WlzDVertex2 **)a)->vtY) < 0.0)
  {
    cmp = -1;
  }
  return(cmp);
}

static int	WlzConvHullClarksonCmpIH(const void *a, const void *b)
{
  return(WlzConvHullClarksonCmpIL(b,a));
}

static int	WlzConvHullClarksonCmpDH(const void *a, const void *b)
{
  return(WlzConvHullClarksonCmpDL(b,a));
}
