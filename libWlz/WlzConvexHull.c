#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvexHull_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzConvexHull.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
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
* \brief	Functions for computing the convex hull of objects.
* \ingroup	WlzConvexHull
*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

static WlzIntervalDomain 	*WlzConvexVtxPolyToItvDom(
				  int nVtx,
				  WlzIVertex2 *vtx,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzConvexHullToContour2D(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzConvexHullToContour3D(
				  WlzConvHullDomain3 *cvh,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzConvexHullToPixDomObj(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzConvexHullToVoxDomObj(
				  WlzConvHullDomain3 *cvh,
				  WlzErrorNum *dstErr);
static WlzObject 		*WlzObjToConvexPolygon3d(
				  WlzObject *obj,
				  WlzErrorNum *dstErr);
static WlzErrorNum		WlzConvexHullPlaneSweep(
				  WlzPlaneDomain *pDom,
				  WlzConvHullDomain3 *cvh);
static int			WlzConvexHullVtxDZToIZ(
				  double d,
				  double c);
static int 			WlzConvHullFcePlnIsn(
				  WlzConvHullDomain3 *cvh,
				  int z,
				  WlzDVertex2 *isn,
				  int f);

/*!
* \return       New 2D convex hull domain.
* \ingroup      WlzConvexHull
* \brief        Allocates a new 2D convex hull domain with space for at
*               least the requested number of vertices.
* \param        maxVtx                  Number of vertices space allocated for.
* \param        vtxType                 Vertex type, must be either
*                                       WLZ_VERTEX_2I or WLZ_VERTEX_2D.
* \param        dstErr                  Desination error pointer, may be NULL.
*/
WlzConvHullDomain2       	*WlzMakeConvexHullDomain2(
                                  int maxVtx,
                                  WlzVertexType vtxType,
                                  WlzErrorNum *dstErr)
{
  size_t        vtxSz;
  WlzConvHullDomain2 *cvh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((cvh = (WlzConvHullDomain2 *)
            AlcCalloc(1, sizeof(WlzConvHullDomain2))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    switch(vtxType)
    {
      case WLZ_VERTEX_I2:
        vtxSz = sizeof(WlzIVertex2);
        break;
      case WLZ_VERTEX_D2:
        vtxSz = sizeof(WlzDVertex2);
        break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((cvh->vertices.v = AlcMalloc(vtxSz * maxVtx)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(cvh)
    {
      AlcFree(cvh->vertices.v);
      AlcFree(cvh);
    }
  }
  else
  {
    cvh->type = WLZ_CONVHULL_DOMAIN_2D;
    cvh->vtxType = vtxType;
    cvh->maxVertices = maxVtx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}

/*!
* \return       New 3D convex hull domain.
* \ingroup      WlzConvexHull
* \brief        Allocates a new 3D convex hull domain with space for at
*               least the requested number of vertices and faces.
* \param        maxVtx                  Number of vertices space allocated for.
* \param        maxFce                  Number of faces space allocated for.
* \param        vtxType                 Vertex type, must be either
*                                       WLZ_VERTEX_3I or WLZ_VERTEX_3D.
* \param        dstErr                  Desination error pointer, may be NULL.
*/
WlzConvHullDomain3       	*WlzMakeConvexHullDomain3(
                                  int maxVtx,
                                  int maxFce,
                                  WlzVertexType vtxType,
                                  WlzErrorNum *dstErr)
{
  size_t        vtxSz;
  WlzConvHullDomain3 *cvh = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((cvh = (WlzConvHullDomain3 *)
            AlcCalloc(1, sizeof(WlzConvHullDomain3))) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    switch(vtxType)
    {
      case WLZ_VERTEX_I3:
        vtxSz = sizeof(WlzIVertex3);
        break;
      case WLZ_VERTEX_D3:
        vtxSz = sizeof(WlzDVertex3);
        break;
      default:
        errNum = WLZ_ERR_PARAM_TYPE;
        break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((cvh->vertices.v = AlcMalloc(vtxSz * maxVtx)) == NULL) ||
       ((cvh->faces = (int *)AlcMalloc(sizeof(int) * 3 * maxFce)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(cvh)
    {
      AlcFree(cvh->vertices.v);
      AlcFree(cvh->faces);
      AlcFree(cvh);
    }
  }
  else
  {
    cvh->type = WLZ_CONVHULL_DOMAIN_3D;
    cvh->vtxType = vtxType;
    cvh->maxVertices = maxVtx;
    cvh->maxFaces = maxFce;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cvh);
}


/*!
* \return       Woolz error code.
* \ingroup      WlzConvexHull
* \brief        Frees the given 2D convex hull domain.
* \param        cvh                     Given 2D convex hull domain.
*/
WlzErrorNum              	WlzFreeConvexHullDomain2(
                                  WlzConvHullDomain2 *cvh)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(cvh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cvh->type != WLZ_CONVHULL_DOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(WlzUnlink(&(cvh->linkcount), &errNum))
  {
    AlcFree(cvh->vertices.v);
    AlcFree(cvh);
  }
  return(errNum);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzConvexHull
* \brief        Frees the given 3D convex hull domain.
* \param        cvh                     Given 3D convex hull domain.
*/
WlzErrorNum              	WlzFreeConvexHullDomain3(
                                  WlzConvHullDomain3 *cvh)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(cvh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cvh->type != WLZ_CONVHULL_DOMAIN_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(WlzUnlink(&(cvh->linkcount), &errNum))
  {
    AlcFree(cvh->vertices.v);
    AlcFree(cvh->faces);
    AlcFree(cvh);
  }
  return(errNum);
}

/*!
* \return	New 2D convex hull domain.
* \ingroup	WlzConvexHull
* \brief	Copies the given 2D convex hull domain which has the
* 		minimum required number of vertices, with verticies of
* 		the required type.
* \param	gCVH			Given convex hull domain.
* \param	rVtxType		Required vertex type, must be
* 					either WLZ_VERTEX_I2 or WLZ_VERTEX_D2.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzConvHullDomain2		*WlzConvexHullCopy2(
				  WlzConvHullDomain2 *gCVH,
				  WlzVertexType rVtxType,
				  WlzErrorNum *dstErr)
{
  WlzConvHullDomain2 *rCVH = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gCVH == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gCVH->type != WLZ_CONVHULL_DOMAIN_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(((rVtxType != WLZ_VERTEX_I2) &&
           (rVtxType != WLZ_VERTEX_D2)) ||
          ((gCVH->vtxType != WLZ_VERTEX_I2) &&
	   (gCVH->vtxType != WLZ_VERTEX_D2)))
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else
  {
    rCVH = WlzMakeConvexHullDomain2(gCVH->nVertices, rVtxType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rCVH->nVertices = gCVH->nVertices;
    if(gCVH->vtxType == rCVH->vtxType)
    {
      size_t	sz;

      sz = (rCVH->vtxType == WLZ_VERTEX_I2)? sizeof(WlzIVertex2):
                                             sizeof(WlzDVertex2);
      (void )memcpy(rCVH->vertices.v, gCVH->vertices.v, sz * gCVH->nVertices);
    }
    else if(rCVH->vtxType == WLZ_VERTEX_I2)
    {
      WlzValueCopyDVertexToIVertex(rCVH->vertices.i2, gCVH->vertices.d2,
      				   gCVH->nVertices);
    }
    else /* rCVH->vtxType == WLZ_VERTEX_D2 */
    {
      WlzValueCopyIVertexToDVertex(rCVH->vertices.d2, gCVH->vertices.i2,
      				   gCVH->nVertices);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rCVH);
}

/*!
* \return	New 3D convex hull domain.
* \ingroup	WlzConvexHull
* \brief	Copies the given 3D convex hull domain which has the
* 		minimum required number of vertices and faces, with
* 		verticies of the required type.
* \param	gCVH			Given convex hull domain.
* \param	rVtxType		Required vertex type, must be
* 					either WLZ_VERTEX_I3 or WLZ_VERTEX_D3.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzConvHullDomain3		*WlzConvexHullCopy3(
				  WlzConvHullDomain3 *gCVH,
				  WlzVertexType rVtxType,
				  WlzErrorNum *dstErr)
{
  WlzConvHullDomain3 *rCVH = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(gCVH == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gCVH->type != WLZ_CONVHULL_DOMAIN_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(((rVtxType != WLZ_VERTEX_I3) &&
           (rVtxType != WLZ_VERTEX_D3)) ||
          ((gCVH->vtxType != WLZ_VERTEX_I3) &&
	   (gCVH->vtxType != WLZ_VERTEX_D3)))
  {
    errNum = WLZ_ERR_PARAM_TYPE;
  }
  else
  {
    rCVH = WlzMakeConvexHullDomain3(gCVH->nVertices, gCVH->nFaces,
                                    rVtxType, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rCVH->nFaces = gCVH->nFaces;
    rCVH->nVertices = gCVH->nVertices;
    if(gCVH->vtxType == rCVH->vtxType)
    {
      size_t	sz;

      sz = (rCVH->vtxType == WLZ_VERTEX_I3)? sizeof(WlzIVertex3):
                                             sizeof(WlzDVertex3);
      (void )memcpy(rCVH->vertices.v, gCVH->vertices.v, sz * gCVH->nVertices);
    }
    else if(rCVH->vtxType == WLZ_VERTEX_I3)
    {
      WlzValueCopyDVertexToIVertex3(rCVH->vertices.i3, gCVH->vertices.d3,
      				    gCVH->nVertices);
    }
    else /* rCVH->vtxType == WLZ_VERTEX_D3 */
    {
      WlzValueCopyIVertexToDVertex3(rCVH->vertices.d3, gCVH->vertices.i3,
      				    gCVH->nVertices);
    }
    (void )memcpy(rCVH->faces, gCVH->faces, sizeof(int) * gCVH->nFaces);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rCVH);
}

/*!
* \return	New object or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new object of the requested type enclosing the same
* 		region of space as the given convex hull object.
* \param	cObj			Given object.
* \param	rType			Requested object types. Valid
* 					types for a 2D convex hull are:
* 					WLZ_2D_DOMAINOBJ and WLZ_CONTOUR.
* 					Valid types for a 3D convex hull are:
* 					WLZ_3D_DOMAINOBJ and WLZ_CONTOUR.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject			*WlzConvexHullToObj(
				  WlzObject *cObj,
				  WlzObjectType rType,
				  WlzErrorNum *dstErr)
{
  WlzObject	*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  if(cObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(cObj->type != WLZ_CONV_HULL)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(cObj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(cObj->domain.core->type == WLZ_CONVHULL_DOMAIN_2D)
  {
    WlzConvHullDomain2 *cvh;

    cvh = cObj->domain.cvh2;
    if((cvh->vtxType != WLZ_VERTEX_I2) && (cvh->vtxType != WLZ_VERTEX_D2))
    {
      errNum = WLZ_ERR_PARAM_TYPE;
    }
    else
    {
      switch(rType)
      {
	case WLZ_2D_DOMAINOBJ:
          rObj = WlzConvexHullToPixDomObj(cvh, &errNum);
	  break;
	case WLZ_CONTOUR:
          rObj = WlzConvexHullToContour2D(cvh, &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
  }
  else if(cObj->domain.core->type == WLZ_CONVHULL_DOMAIN_3D)
  {
    WlzConvHullDomain3 *cvh;

    cvh = cObj->domain.cvh3;
    if((cvh->vtxType != WLZ_VERTEX_I3) && (cvh->vtxType != WLZ_VERTEX_D3))
    {
      errNum = WLZ_ERR_PARAM_TYPE;
    }
    else
    {
      switch(rType)
      {
	case WLZ_3D_DOMAINOBJ:
          rObj = WlzConvexHullToVoxDomObj(cvh, &errNum);
	  break;
	case WLZ_CONTOUR:
          rObj = WlzConvexHullToContour3D(cvh, &errNum);
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
  }
  else
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Convex hull object, NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes the convex hull of the given object.
* 		If a degenerate convex hull is returned then the returned
* 		error will be WLZ_ERR_DEGENERATE.
* \param	gObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzObjToConvexHull(
  				  WlzObject *gObj,
  				  WlzErrorNum *dstErr)
{
  int		nVtx = 0;
  WlzDomain	dom;
  WlzValues	nullVal;
  WlzObject	*cObj = NULL;
  WlzVertexP	vtx;
  WlzVertexType vType;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  vtx.v = NULL;
  dom.core = NULL;
  nullVal.core = NULL;
  if(gObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	{
	  WlzObject *bObj = NULL;

	  if(gObj->domain.core == NULL)
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  else
	  {
	    bObj = WlzObjToBoundary(gObj, 0, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    cObj = WlzObjToConvexHull(bObj, &errNum);
	  }
	  (void )WlzFreeObj(bObj);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	{
	  int 	nPln;
	  WlzPlaneDomain *pDom;
	  WlzObject	**allObj = NULL;

	  if(gObj->domain.core == NULL)
	  {
	    errNum = WLZ_ERR_DOMAIN_NULL;
	  }
	  else
	  {
	    pDom = gObj->domain.p;
	    /* Allocate an array for planar convex hulls. */
	    if((nPln = pDom->lastpl - pDom->plane1 + 1) <= 0)
	    {
	      errNum = WLZ_ERR_DOMAIN_DATA;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if((allObj = (WlzObject **)
	                 AlcCalloc(nPln, sizeof(WlzObject *))) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  /* For each plane compute the planar convex hull. */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int		p;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	    for(p = 0; p < nPln; ++p)
	    {
	      WlzDomain	*dom2;
	      WlzObject	*obj2D = NULL;
	      WlzErrorNum errNum2 = WLZ_ERR_NONE;

	      if(errNum == WLZ_ERR_NONE)
	      {
		dom2 = pDom->domains + p;
		if((*dom2).core != NULL)
		{
		  obj2D = WlzAssignObject(
			  WlzMakeMain(WLZ_2D_DOMAINOBJ, *dom2, nullVal,
				      NULL, NULL, &errNum2), NULL);
		  if(errNum2 == WLZ_ERR_NONE)
		  {
		    allObj[p] = WlzObjToConvexHull(obj2D, &errNum2);
		  }
		  if((errNum2 != WLZ_ERR_NONE) &&
		     (errNum2 != WLZ_ERR_DEGENERATE))
		  {
#ifdef _OPENMP
#pragma omp critical
		    {
#endif
		      if(errNum == WLZ_ERR_NONE)
		      {
		        errNum = errNum2;
		      }
#ifdef _OPENMP
		    }
#endif
		  }
		  (void )WlzFreeObj(obj2D);
		}
	      }
	    }
	  }
	  /* Create an array of 3D integer vertices which are the vertices
	   * of the planar convex hulls. */
	  if(errNum == WLZ_ERR_NONE)           
	  {
	    int		p;                     

	    for(p = 0; p < nPln; ++p)          
	    {
	      if(allObj[p])
	      {
	        nVtx += allObj[p]->domain.cvh2->nVertices;
	      }
	    }
	    if((vtx.i3 = (WlzIVertex3 *)
	                 AlcMalloc(sizeof(WlzIVertex3) * nVtx)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int		p;         
	    WlzIVertex3 *v0;

	    v0 = vtx.i3;
	    for(p = 0; p < nPln; ++p)
	    {
	      if(allObj[p])
	      {
		int	j,
			z;
                WlzIVertex2 *v1;
		WlzConvHullDomain2 *cvh;
		
		z = pDom->plane1 + p;
		cvh = allObj[p]->domain.cvh2;
		v1 = cvh->vertices.i2;
		for(j = 0; j < cvh->nVertices; ++j)
		{
		  WLZ_VTX_3_SET((*v0), v1->vtX, v1->vtY, z); 
		  ++v0;
		  ++v1;
		}
	      }
	    }
	  }
	  /* Free all planar convex hulls. */
	  if(allObj)
	  {
	    int		p;         

	    for(p = 0; p < nPln; ++p)
	    {
	      (void )WlzFreeObj(allObj[p]);
	    }
	    AlcFree(allObj);
	  }
	  /* Compute the 3D convex hull of all the planar convex hulls. */
	  if(errNum == WLZ_ERR_NONE)           
	  {
	    vType = WLZ_VERTEX_I3;
	    dom.cvh3 = WlzConvexHullFromVtx3(vType, nVtx, vtx, &errNum);
	  }
	}
        break;
      case WLZ_2D_POLYGON:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
        {
	  nVtx = gObj->domain.poly->nvertices;
	  switch(gObj->domain.poly->type)
	  {
	    case WLZ_POLYGON_INT:
	      vType = WLZ_VERTEX_I2;
	      vtx.i2 = gObj->domain.poly->vtx;
	      break;
	    case WLZ_POLYGON_DOUBLE:
	      vType = WLZ_VERTEX_D2;
	      vtx.d2 = (WlzDVertex2 *)(gObj->domain.poly->vtx);
	      break;
	    default:
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    dom.cvh2 = WlzConvexHullFromVtx2(vType, nVtx, vtx, &errNum);
	  }
	}
        break;
      case WLZ_BOUNDLIST:
	/* Here we only collect those vertices which are at the top level of
	 * the boundary list, ie next but not up/down traversal. */
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  WlzBoundList *b0,
	  	       *b1,
		       *b2;

	  b1 = b0 = gObj->domain.b;
	  nVtx = b0->poly->nvertices;
	  while((b2 = b1->next) != NULL)
	  {
	    b1 = b2;
	    nVtx += b1->poly->nvertices;
	    if(b1->poly->type != b0->poly->type)
	    {
	      errNum = WLZ_ERR_DOMAIN_TYPE;
	      break;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    switch(b0->poly->type)
	    {
	      case WLZ_POLYGON_INT:
		vType = WLZ_VERTEX_I2;
	        if((vtx.i2 = (WlzIVertex2 *)
		             AlcMalloc(sizeof(WlzIVertex2) * nVtx)) == NULL)
	        {
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      case WLZ_POLYGON_DOUBLE:
		vType = WLZ_VERTEX_D2;
	        if((vtx.d2 = (WlzDVertex2 *)
		             AlcMalloc(sizeof(WlzDVertex2) * nVtx)) == NULL)
	        {
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		break;
	      default:
	        errNum = WLZ_ERR_DOMAIN_TYPE;
		break;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int		n;

	    b1 = b0 = gObj->domain.b;
	    n = b0->poly->nvertices;
	    if(vType == WLZ_VERTEX_I2)
	    {
	      WlzIVertex2 *p0;

	      p0 = vtx.i2;
	      (void )memcpy(p0, b0->poly->vtx, n * sizeof(WlzIVertex2));
	      n = b0->poly->nvertices;
	      while((b1 = b1->next) != NULL)
	      {
	        p0 += n;
		n = b1->poly->nvertices;
	        (void )memcpy(p0, b1->poly->vtx, n * sizeof(WlzIVertex2));
	      }
	    }
	    else /* vType == WLZ_VERTEX_D2 */
	    {
	      WlzDVertex2 *p0,
	      		  *p1;

	      p0 = vtx.d2;
	      p1 = (WlzDVertex2 *)(b0->poly->vtx);
	      (void )memcpy(p0, p1, n * sizeof(WlzIVertex2));
	      while((b1 = b1->next) != NULL)
	      {
	        p0 += n;
		n = b1->poly->nvertices;
		p1 = (WlzDVertex2 *)(b1->poly->vtx);
	        (void )memcpy(p0, p1, n * sizeof(WlzIVertex2));
	      }
	    }
	    dom.cvh2 = WlzConvexHullFromVtx2(vType, nVtx, vtx, &errNum);
	  }
	}
        break;
      case WLZ_CONTOUR:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  vtx = WlzVerticesFromObj(gObj, NULL, &nVtx, &vType, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    switch(vType)
	    {
	      case WLZ_VERTEX_I2: /* FALLTHROUGH */
	      case WLZ_VERTEX_D2:
		dom.cvh2 = WlzConvexHullFromVtx2(vType, nVtx, vtx, &errNum);
		break;
	      case WLZ_VERTEX_I3: /* FALLTHROUGH */
	      case WLZ_VERTEX_D3:
		dom.cvh3 = WlzConvexHullFromVtx3(vType, nVtx, vtx, &errNum);
		break;
	      default:
	        errNum = WLZ_ERR_DOMAIN_TYPE;
		break;
	    }
	  }
	}
        break;
      case WLZ_CMESH_2D:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  vtx = WlzVerticesFromObj(gObj, NULL, &nVtx, &vType, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dom.cvh2 = WlzConvexHullFromVtx2(vType, nVtx, vtx, &errNum);
	}
        break;
      case WLZ_CMESH_2D5: /* FALLTHROUGH */
      case WLZ_CMESH_3D:
	if(gObj->domain.core == NULL)
	{
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else
	{
	  vtx = WlzVerticesFromObj(gObj, NULL, &nVtx, &vType, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dom.cvh3 = WlzConvexHullFromVtx3(vType, nVtx, vtx, &errNum);
	}
        break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  /* For some object types we have a convex hull object for others just
   * the domain. */
  AlcFree(vtx.v);
  if(cObj == NULL) 
  {
    if((dom.core != NULL) &&
       ((errNum == WLZ_ERR_NONE) || (errNum == WLZ_ERR_DEGENERATE)))
    {
      WlzErrorNum saveErrNum;

      saveErrNum = errNum;
      cObj = WlzMakeMain(WLZ_CONV_HULL, dom, nullVal, NULL, NULL, &errNum);
      if((errNum == WLZ_ERR_NONE) && (saveErrNum == WLZ_ERR_DEGENERATE))
      {
        errNum = WLZ_ERR_DEGENERATE;
      }
    }
    if((errNum != WLZ_ERR_NONE) && (errNum != WLZ_ERR_DEGENERATE) &&
       (dom.core != NULL))
    {
      switch(dom.core->type)
      {
	case WLZ_CONVHULL_DOMAIN_2D:
	  (void )WlzFreeConvexHullDomain2(dom.cvh2);
	  break;
	case WLZ_CONVHULL_DOMAIN_3D:
	  (void )WlzFreeConvexHullDomain3(dom.cvh3);
	  break;
	default:
	  break;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(cObj);
}

/*!
* \return	Convex hull object.
* \ingroup      WlzConvexHull
* \brief	Construct the minimal convex polygonal cover from interval
* 		domain.
* \param	obj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject 			*WlzObjToConvexPolygon(
  				  WlzObject *obj,
  				  WlzErrorNum *dstErr)
{
  WlzObject 		*cvh = NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzIVertex2 		*wtx = NULL, *w1tx = NULL, *w2tx = NULL;
  WlzIntervalWSpace 	iwsp;
  int 			nhalfway;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* check object */
  if( obj ){
    switch( obj->type )
    {
      
    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* check the planedomain type */
      if( obj->domain.p ){
	switch( obj->domain.p->type ){
	case WLZ_PLANEDOMAIN_DOMAIN:
	  return WlzObjToConvexPolygon3d(obj, dstErr);

	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
	}
      }
      else {
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_TRANS_OBJ:
      if((obj->values.core) && 
	 (values.obj = WlzObjToConvexPolygon(obj->values.obj, &errNum))){
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }
  else {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  /* now the real algorithm, 2D only, definitely a domain at this point.
     Make a polygon object with the maximum number of vertices for the
     convex polygon = (2 * num_lines + 1) the extra is to allow the polygon
     to be closed i.e. first = last */
  if( errNum == WLZ_ERR_NONE ){
    if((domain.poly = WlzMakePolygonDomain(WLZ_POLYGON_INT, 0, NULL,
				      3+2*(obj->domain.i->lastln -
					   obj->domain.i->line1),
				      1, &errNum))){
      values.core = NULL;
      cvh = WlzMakeMain(WLZ_2D_POLYGON, domain, values, NULL, NULL, &errNum);
    }
  }
  
  if( errNum == WLZ_ERR_NONE ){
    wtx = cvh->domain.poly->vtx;
    /*
     * proceed down right hand side of object
     */
    if( (errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_ILIC))
       == WLZ_ERR_NONE ){
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	/*
	 * set up first chord
	 */
	if (iwsp.linpos == obj->domain.i->line1) {
	  if (iwsp.nwlpos == 1) {
	    wtx->vtX = iwsp.lftpos;
	    wtx->vtY = iwsp.linpos;
	    wtx++;
	  }
	  if (iwsp.intrmn == 0) {
	    wtx->vtX = iwsp.rgtpos;
	    wtx->vtY = iwsp.linpos;
	    wtx++;
	    cvh->domain.poly->nvertices = 2;
	    w1tx = wtx-1;
	    w2tx = wtx-2;
	  }
	} else {
	  /*
	     * add extra chords, checking concavity condition
	     */
	  if (iwsp.intrmn == 0) {
	    wtx->vtX = iwsp.rgtpos;
	    wtx->vtY = iwsp.linpos;
	    cvh->domain.poly->nvertices++;
	    /*
	     * Concavity condition (may propagate backwards).
	     * Also deals satisfactorily with the case that first
	     * line consists of a single interval, itself a single point.
	     */
	    while ((cvh->domain.poly->nvertices >= 3) &&
		   (wtx->vtY-w2tx->vtY)*(w1tx->vtX-w2tx->vtX) <=
		   (w1tx->vtY-w2tx->vtY)*(wtx->vtX-w2tx->vtX)) {
	      w1tx->vtX = wtx->vtX;
	      w1tx->vtY = wtx->vtY;
	      wtx--;
	      w1tx--;
	      w2tx--;
	      cvh->domain.poly->nvertices--;
	    }
	    wtx++;
	    w1tx++;
	    w2tx++;
	  }
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /*
     * now proceed up left hand side of object
     */
    if( (errNum = WlzInitRasterScan(obj, &iwsp, WLZ_RASTERDIR_DLDC))
       == WLZ_ERR_NONE ){
      nhalfway = cvh->domain.poly->nvertices + 2;
      while((errNum = WlzNextInterval(&iwsp)) == WLZ_ERR_NONE ){
	if (iwsp.intrmn == 0) {
	  wtx->vtX = iwsp.lftpos;
	  wtx->vtY = iwsp.linpos;
	  cvh->domain.poly->nvertices++;
	  /*
	   * Concavity condition (may propagate backwards).
	   * Also deals satisfactorily with the case that last
	   * line consists of a single interval, itself a single point.
	   */
	  while ((cvh->domain.poly->nvertices >= nhalfway) &&
		 (wtx->vtY-w2tx->vtY)*(w1tx->vtX-w2tx->vtX) <=
		 (w1tx->vtY-w2tx->vtY)*(wtx->vtX-w2tx->vtX)) {
	    w1tx->vtX = wtx->vtX;
	    w1tx->vtY = wtx->vtY;
	    wtx--;
	    w1tx--;
	    w2tx--;
	    cvh->domain.poly->nvertices--;
	  }
	  wtx++;
	  w1tx++;
	  w2tx++;
	}
      }
      if( errNum == WLZ_ERR_EOO ){
	errNum = WLZ_ERR_NONE;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return cvh;
}

/*!
* \return	New object or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new 2D contour object enclosing the same
* 		region of space as the given convex hull object.
* 		This function assumes that the given object is a 2D
* 		convex hull and the required object type for return is
* 		a 2D contour.
* \param	cObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzConvexHullToContour2D(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr)
{
  int		vHTSz;
  WlzGMModelType mType;
  WlzObject	*rObj = NULL;
  WlzGMModel	*model = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  dom.core = NULL;
  val.core = NULL;
  mType = (cvh->vtxType == WLZ_VERTEX_I2)?  WLZ_GMMOD_2I: WLZ_GMMOD_2D;
  vHTSz = (cvh->nVertices < 1024)? 1024: cvh->nVertices;
  model = WlzGMModelNew(mType, 0, vHTSz, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    int 	i,
		n;
    WlzDVertex2 buf[2];

    n = cvh->nVertices;
    if(cvh->vtxType == WLZ_VERTEX_I2)
    {
      buf[1].vtX = cvh->vertices.i2[0].vtX;
      buf[1].vtY = cvh->vertices.i2[0].vtY;
      for(i = 0; i < n; ++i)
      {
	buf[0] = buf[1];
	buf[1].vtX = cvh->vertices.i2[i].vtX;
	buf[1].vtY = cvh->vertices.i2[i].vtY;
	if((errNum = WlzGMModelConstructSimplex2D(model, buf)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        buf[0] = buf[1];
	buf[1].vtX = cvh->vertices.i2[0].vtX;
	buf[1].vtY = cvh->vertices.i2[0].vtY;
	errNum = WlzGMModelConstructSimplex2D(model, buf);
      }
    }
    else /* cvh->vtxType == WLZ_VERTEX_D2 */
    {
      buf[1] = cvh->vertices.d2[0];
      for(i = 0; i < n; ++i)
      {
	buf[0] = buf[1];
	buf[1] = cvh->vertices.d2[i];
	if((errNum = WlzGMModelConstructSimplex2D(model, buf)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        buf[0] = buf[1];
	buf[1] = cvh->vertices.d2[0];
	errNum = WlzGMModelConstructSimplex2D(model, buf);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr = WlzMakeContour(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr->model = WlzAssignGMModel(model, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
  }
  else
  {
    if(dom.ctr)
    {
      (void )WlzFreeContour(dom.ctr);
    }
    else if(model)
    {
      (void )WlzGMModelFree(model);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new 3D contour object enclosing the same
* 		region of space as the given convex hull object.
* 		This function assumes that the given object is a 3D
* 		convex hull and the required object type for return is
* 		a 3D contour.
* \param	cObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzConvexHullToContour3D(
				  WlzConvHullDomain3 *cvh,
				  WlzErrorNum *dstErr)
{
  int		vHTSz;
  WlzGMModelType mType;
  WlzObject	*rObj = NULL;
  WlzGMModel	*model = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  dom.core = NULL;
  val.core = NULL;

  mType = WLZ_GMMOD_3D;
  vHTSz = (cvh->nVertices < 1024)? 1024: cvh->nVertices;
  model = WlzGMModelNew(mType, 0, vHTSz, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    int 	i,
		n;
    WlzDVertex3 buf[3];

    n = cvh->nFaces;
    if(cvh->vtxType == WLZ_VERTEX_I3)
    {
      for(i = 0; i < n; ++i)
      {
	int	i3,
		j;

	i3 = 3 * i;
	for(j = 0; j < 3; ++j)
	{
	  WlzIVertex3 v;

	  v = cvh->vertices.i3[cvh->faces[i3 + j]];
	  buf[j].vtX = v.vtX;
	  buf[j].vtY = v.vtY;
	  buf[j].vtZ = v.vtZ;
	}
	if((errNum = WlzGMModelConstructSimplex3D(model, buf)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
    else /* cvh->vtxType == WLZ_VERTEX_D3 */
    {
      for(i = 0; i < n; ++i)
      {
	int	i3,
		j;

	i3 = 3 * i;
	for(j = 0; j < 3; ++j)
	{
	  buf[j] = cvh->vertices.d3[cvh->faces[i3 + j]];
	}
	if((errNum = WlzGMModelConstructSimplex3D(model, buf)) != WLZ_ERR_NONE)
	{
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr = WlzMakeContour(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr->model = WlzAssignGMModel(model, NULL);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
  }
  else
  {
    if(dom.ctr)
    {
      (void )WlzFreeContour(dom.ctr);
    }
    else if(model)
    {
      (void )WlzGMModelFree(model);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	New object or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new 2D spatial (pixel) domain object which
* 		occupies the region of space enclosdd by the given convex
* 		hull object.
* 		This function assumes that the given object is a 2D
* 		convex hull and the required object type for return is
* 		a 2D spatial domain object.
* \param	cObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzConvexHullToPixDomObj(
				  WlzConvHullDomain2 *cvh,
				  WlzErrorNum *dstErr)
{
  int		freeVtx = 0;
  WlzObject	*rObj = NULL;
  WlzIVertex2	*vtx = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if(cvh->vtxType == WLZ_VERTEX_I2)
  {
    /* Just set the vertex pointer, no copy needed. */
    vtx = cvh->vertices.i2;
  }
  else /* cvh->vtxType == WLZ_VERTEX_D2 */
  {
    /* Need to create an array of integer vertices from those of the convex
     * hull. Do this making sure that they enclose the vertices of the given
     * convex hull. */
    if((vtx = AlcMalloc(sizeof(WlzIVertex2) * cvh->nVertices)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int	  i;
      WlzDVertex2 c;

      freeVtx = 1;
      c = cvh->centroid.d2;
      for(i = 0; i < cvh->nVertices; ++i)
      {
	vtx[i] = WlzConvexHullVtxD2ToI2(cvh->vertices.d2[i], c);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dom.i = WlzConvexVtxPolyToItvDom(cvh->nVertices, vtx, &errNum); 
  }
  if(errNum == WLZ_ERR_NONE)
  {
    rObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreeIntervalDomain(dom.i);
    }
  }
  if(freeVtx)
  {
    AlcFree(vtx);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return	Integral vertex.
* \ingroup	WlzConvexHull
* \brief	Convert the given double vertex to an integer vertex
* 		such that it will enclose the given convex hull double
* 		vertex.
* \param	d			Given convex hull double vertex.
* \param	c			Centroid of the convex hull.
*/
WlzIVertex2			WlzConvexHullVtxD2ToI2(
				  WlzDVertex2 d,
				  WlzDVertex2 c)
{
  WlzIVertex2 i;

  i.vtX = (int)((d.vtX < c.vtX)? floor(d.vtX): ceil(d.vtX));
  i.vtY = (int)((d.vtY < c.vtY)? floor(d.vtY): ceil(d.vtY));
  return(i);
}

/*!
* \return	Integral vertex.
* \ingroup	WlzConvexHull
* \brief	Convert the given double vertex to an integer vertex
* 		such that it will enclose the given convex hull double
* 		vertex.
* \param	d			Given convex hull double vertex.
* \param	c			Centroid of the convex hull.
*/
WlzIVertex3			WlzConvexHullVtxD3ToI3(
				  WlzDVertex3 d,
				  WlzDVertex3 c)
{
  WlzIVertex3 i;

  i.vtX = (int)((d.vtX < c.vtX)? floor(d.vtX): ceil(d.vtX));
  i.vtY = (int)((d.vtY < c.vtY)? floor(d.vtY): ceil(d.vtY));
  i.vtZ = (int)((d.vtZ < c.vtZ)? floor(d.vtZ): ceil(d.vtZ));
  return(i);
}

/*!
* \return	Integral z coordinate.
* \ingroup	WlzConvexHull
* \brief	Convert the given double vertex z coordinate to an
* 		integer value which will enclose the given convex hull
* 		double vertex.
* \param	d			Given convex hull double vertex
* 					z coordinate.
* \param	c			Convex hull centroid z coordinate.
*/
static int			WlzConvexHullVtxDZToIZ(
				  double d,
				  double c)
{
  int	 i;

  i = (int)((d < c)? floor(d): ceil(d));
  return(i);
}

/*!
* \return	New object or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new 3D spatial (voxel) domain object which
* 		occupies the region of space enclosdd by the given convex
* 		hull object.
* 		This function assumes that the given object is a 3D
* 		convex hull and the required object type for return is
* 		a 3D spatial domain object.
* \param	cObj			Given object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzConvexHullToVoxDomObj(
				  WlzConvHullDomain3 *cvh,
				  WlzErrorNum *dstErr)
{
  int		plMin,
  		plMax;
  WlzObject	*rObj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  dom.core = NULL;
  val.core = NULL;
  /* Find vertices with minimum and maximum plane coordinate. */
  if(cvh->vtxType == WLZ_VERTEX_I3)
  {
    int		i;
    WlzIVertex3	*pZMin,
    		*pZMax;

    pZMin = pZMax = cvh->vertices.i3;
    for(i = 1; i < cvh->nVertices; ++i)
    {
      WlzIVertex3 *p;

      p = cvh->vertices.i3 + i;
      if(p->vtZ < pZMin->vtZ)
      {
        pZMin = p;
      }
      else if(p->vtZ > pZMax->vtZ)
      {
        pZMax = p;
      }
    }
    plMin = pZMin->vtZ;
    plMax = pZMax->vtZ;
  }
  else /* cvh->vtxType == WLZ_VERTEX_D3 */
  {
    int		i;
    WlzDVertex3	*pZMin,
    		*pZMax;

    pZMin = pZMax = cvh->vertices.d3;
    for(i = 1; i < cvh->nVertices; ++i)
    {
      WlzDVertex3 *p;

      p = cvh->vertices.d3 + i;
      if(p->vtZ < pZMin->vtZ)
      {
        pZMin = p;
      }
      else if(p->vtZ > pZMax->vtZ)
      {
        pZMax = p;
      }
    }
    plMin = (int )floor(pZMin->vtZ);
    plMax = (int )ceil(pZMax->vtZ);
  }
  /* Make a planedomain for the new domain using the minimum and maximum
   * plane coordinates, we will set the line and column bounds latter. */
  dom.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN, plMin, plMax,
		   	     0, 0, 0, 0, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzConvexHullPlaneSweep(dom.p, cvh);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i,
    		nPln;
    WlzDomain	*doms;

    /* Set the planedomain line and column bounds, here every plane is known
     * to have a domain. */
    doms = dom.p->domains;
    nPln = plMax - plMin + 1;
    dom.p->kol1   = (*doms).i->kol1;
    dom.p->lastkl = (*doms).i->lastkl;
    dom.p->line1  = (*doms).i->line1;
    dom.p->lastln = (*doms).i->lastln;
    for(i = 1; i < nPln; ++i)
    {
      if(dom.p->kol1 > (*doms).i->kol1)
      {
        dom.p->kol1 = (*doms).i->kol1;
      }
      if(dom.p->lastkl < (*doms).i->lastkl)
      {
        dom.p->lastkl = (*doms).i->lastkl;
      }
      if(dom.p->line1 > (*doms).i->line1)
      {
        dom.p->line1 = (*doms).i->line1;
      }
      if(dom.p->lastln < (*doms).i->lastln)
      {
        dom.p->lastln = (*doms).i->lastln;
      }
      ++doms;
    }
    /* Make the object. */
    rObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, dom, val, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzFreePlaneDomain(dom.p);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(rObj);
}

/*!
* \return       New interval domain or NULL on error.
* \ingroup	WlzConvexHull
* \brief	Computes a new interval domain from the given convex vertices
* 		which must be in convex polygon order.
* \param	nVtx			Number of vertices.
* \param	vtx			The convex vertices, may be modified
* 					on return..
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzIntervalDomain 	*WlzConvexVtxPolyToItvDom(
				  int nVtx,
				  WlzIVertex2 *vtx,
				  WlzErrorNum *dstErr)
{
  int		i;
  WlzIVertex2	boxL,
  		boxU;
  WlzInterval	*intervals = NULL;
  WlzIntervalDomain *iDom = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Find bounding box. */
  boxL.vtX = boxU.vtX = vtx[0].vtX;
  boxL.vtY = boxU.vtY = vtx[0].vtY;
  for(i = 1; i < nVtx; ++i)
  {
    if(boxU.vtX < vtx[i].vtX)
    {
      boxU.vtX = vtx[i].vtX;
    }
    else if(boxL.vtX > vtx[i].vtX)
    {
      boxL.vtX = vtx[i].vtX;
    }
    if(boxU.vtY < vtx[i].vtY)
    {
      boxU.vtY = vtx[i].vtY;
    }
    else if(boxL.vtY > vtx[i].vtY)
    {
      boxL.vtY = vtx[i].vtY;
    }
  }
  iDom = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
  			       boxL.vtY, boxU.vtY, boxL.vtX, boxU.vtX,
			       &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate intervals only one per line because convex. */
    if((intervals = (WlzInterval *)
                    AlcMalloc((boxU.vtY - boxL.vtY + 1) *
                              sizeof(WlzInterval))) == NULL)
    {
      (void )WlzFreeIntervalDomain(iDom);
      errNum = WLZ_ERR_NONE;
      iDom = NULL;
    }
    else
    {
      iDom->freeptr = AlcFreeStackPush(iDom->freeptr, (void *)intervals, NULL);
      for(i = boxL.vtY; i <= boxU.vtY; ++i)
      {
	int	j;

	/* Make intervals, one for each line and set left point to -ve
	 * value. */
	j = i - boxL.vtY;
	intervals[j].ileft = -1;
        (void )WlzMakeInterval(i, iDom, 1, intervals + j);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzIVertex2 p0,
    		p1;

    /* Rasterise the convex polygon by setting intervals using the
     * polygon edges, traversing the edges using Bresnham's line
     * algorithm. */
    WLZ_VTX_2_SUB(p1, vtx[nVtx - 1], boxL);
    for(i = 0; i < nVtx; ++i)
    {
      int	  e;
      WlzIVertex2 d,
      		  s;
      WlzInterval *itv;

      p0 = p1;
      WLZ_VTX_2_SUB(p1, vtx[i], boxL);
      WLZ_VTX_2_SUB(d, p1, p0);
      s.vtX = ((d.vtX > 0) << 1) - 1;
      s.vtY = ((d.vtY > 0) << 1) - 1;
      WLZ_VTX_2_ABS(d,d);
      e = ((d.vtX > d.vtY)? d.vtX: -d.vtY) / 2;
      itv = intervals + p0.vtY;
      while(1)
      {
        int	e2;

	if(itv->ileft < 0)
	{
	  itv->ileft = itv->iright = p0.vtX;
	}
	else
	{
	  if(p0.vtX < itv->ileft)
	  {
	    itv->ileft = p0.vtX;
	  }
	  else if(p0.vtX > itv->iright)
	  {
	    itv->iright = p0.vtX;
	  }
	}
	if((p0.vtX == p1.vtX) && (p0.vtY == p1.vtY))
	{
	  break;
	}
	e2 = e;
	if(e2 > -d.vtX)
	{
	  e      -= d.vtY;
	  p0.vtX += s.vtX;
	}
	if(e2 < d.vtY)
	{
	  e      += d.vtX;
	  p0.vtY += s.vtY;
	  itv = intervals + p0.vtY;
	}
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(iDom);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzConvexHull
* \brief	Creates new 2D domains for the given 3D plane domain so
* 		that they enclose the given convex hull.
* 		This function assumes that the given plane domain and
* 		convex hull are both valid.
* \param	pDom			Given plane domain.
* \param	cvh			Given convex hull with integer
* 					vertices.
*/
static WlzErrorNum		WlzConvexHullPlaneSweep(
				  WlzPlaneDomain *pDom,
				  WlzConvHullDomain3 *cvh)
{
  int		iMin,      	/* Index into minimum face plane coordinate
			   	 * permutation table for the first face
		           	 * which intersects the plane. */
		nPln,
		nFceBuf = 0,
     		fceBufMax = 0;
  int		*fcePrm,
		*fceMin,
		*allBuf = NULL,
		*fceBuf = NULL;
  WlzDVertex2	*isnBuf = NULL;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  nPln = pDom->lastpl - pDom->plane1 + 1;
  /* Allocate buffers for the face permutation index and face minimum plane
   * coordinate arrays. */
  if((allBuf = AlcMalloc(sizeof(int) * 2 * cvh->nFaces)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    int		i,
    		p,
		n;

    fcePrm = allBuf;
    fceMin = allBuf + cvh->nFaces;
    /* Fill permutation order and minimum plane coordinate arrays for
     * the faces such that the permutation array indexes the faces
     * by increasing minimum plane coordinate. */
    for(i = 0; i < cvh->nFaces; ++i)
    {
      int	j,
		z,
		i3;

      i3 = 3 * i;
      fcePrm[i] = i;
      if(cvh->vtxType == WLZ_VERTEX_I3)
      {
	z = cvh->vertices.i3[cvh->faces[i3 + 0]].vtZ;
	fceMin[i] = z;
	for(j = 1; j < 3; ++j)
	{
	  z = cvh->vertices.i3[cvh->faces[i3 + j]].vtZ;
	  if(z < fceMin[i])
	  {
	    fceMin[i] = z;
	  }
	}
      }
      else
      {
	z = WlzConvexHullVtxDZToIZ(cvh->vertices.d3[cvh->faces[i3 + 0]].vtZ,
	                           cvh->centroid.d3.vtZ);
	fceMin[i] = z;
	for(j = 1; j < 3; ++j)
	{
	  z = WlzConvexHullVtxDZToIZ(cvh->vertices.d3[cvh->faces[i3 + j]].vtZ,
	                             cvh->centroid.d3.vtZ);
	  if(z < fceMin[i])
	  {
	    fceMin[i] = z;
	  }
	}
      }
    }
    AlgHeapSortIdx(fceMin, fcePrm, cvh->nFaces, AlgHeapSortCmpIdxIFn);
    /* For each plane of the planedomain, find the intersection of the
     * plane with the convex hull and construct an enclosing interval
     * domain for that plane. */
    iMin = 0;
    nFceBuf = 0;
    fceBufMax = 0;
    for(p = 0; p < nPln; ++p)
    {
      int	z,
      		bufIdx = 0,
		nIsn = 0;
      WlzConvHullDomain2 *cvh2 = NULL;

      z = p + pDom->plane1;
      /* Compute intersections of faces in the face buffer with this plane,
       * while updateing the face buffer. */
      for(i = 0; i < nFceBuf; ++i)
      {
	n = WlzConvHullFcePlnIsn(cvh, z, isnBuf + nIsn, fcePrm[fceBuf[i]]);
	if(n > 0)
	{
	  iMin = fceBuf[i] + 1;
	  fceBuf[bufIdx++] = fceBuf[i];
	  nIsn += n;
	}
      }
      for(i = iMin; (i < cvh->nFaces) && (fceMin[fcePrm[i]] <= z); ++i)
      {
	if(fceBufMax <= bufIdx)
	{
	  /* Need to reallocate the face and intersection buffers if there's
	   * not enough room. */
	  fceBufMax = ALG_MAX(2 * bufIdx, 256);
	  if(((fceBuf = (int *)AlcRealloc(fceBuf,
		   fceBufMax * sizeof(int))) == NULL) ||
	      ((isnBuf = (WlzDVertex2 *)AlcRealloc(isnBuf,
		   3 * fceBufMax * sizeof(WlzDVertex2))) == NULL))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	}
	/* Add face and intersections to the buffers. */
	n = WlzConvHullFcePlnIsn(cvh, z, isnBuf + nIsn, fcePrm[i]);
	if(n > 0)
	{
	  fceBuf[bufIdx++] = i;
	  nIsn += n;
	}
      }
      nFceBuf = bufIdx;
      /* Compute planar convex hull. */
      if(errNum == WLZ_ERR_NONE)
      {
	WlzVertexP isnP;

	isnP.d2 = isnBuf;
	cvh2 = WlzConvexHullFromVtx2(WLZ_VERTEX_D2, nIsn, isnP, &errNum);
      }
      /* Create interval domain from the planar convex hull. */
      if((errNum == WLZ_ERR_NONE) || (errNum == WLZ_ERR_DEGENERATE))
      {
        WlzObject *cObj2 = NULL;
        
        cObj2 = WlzConvexHullToPixDomObj(cvh2, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
          WlzDomain *doms;

          doms = pDom->domains + p;
	  (*doms) = WlzAssignDomain(cObj2->domain, NULL);
	}
	(void )WlzFreeObj(cObj2);
      }
      (void )WlzFreeConvexHullDomain2(cvh2);
    }
  }
  AlcFree(allBuf);
  AlcFree(fceBuf);
  AlcFree(isnBuf);
  return(errNum);
}

/*!
* \return	Number of intersection vertices added to the intersection
* 		buffer.
* \ingroup	WlzConvexHull
* \brief	Computes intersections between the face edges and a constant
* 		z plane adding these to the 2D double vertex intersection
* 		buffer.
* \param	cvh		Given 3D convex hull domain.
* \param	z		Z (plane) coordinate.
* \param	isn		Intersection buffer.
* \param	f		Index of face in the convex hull.
*/
static int 			WlzConvHullFcePlnIsn(
				  WlzConvHullDomain3 *cvh,
				  int z,
				  WlzDVertex2 *isn,
				  int f)
{
  int		i,
		j,
  		f3,
  		n = 0;
  const double	eps = 1.0e-6;

  f3 = f * 3;
  if(cvh->vtxType == WLZ_VERTEX_I3)
  {
    WlzIVertex3 v[3];

    for(i = 0; i < 3; ++i)
    {
      v[i] = cvh->vertices.i3[cvh->faces[f3 + i]];
    }
    for(i = 0; i < 3; ++i)
    {
      if(v[i].vtZ == z)
      {
	isn[n].vtX = v[i].vtX;
	isn[n].vtY = v[i].vtY;
	++n;
      }
      else
      {
	j = (i + 1) % 3;
	if(v[j].vtZ != z)
	{
	  int	 d;

	  d = v[j].vtZ - v[i].vtZ;
	  if(d != 0)
	  {
	    double f;

	    f = (double )(z - v[i].vtZ) / (double )(d);
	    if((f > -(eps)) && (f < (1.0 + eps)))
	    {
	      isn[n].vtX = (f * (v[j].vtX - v[i].vtX)) + v[i].vtX;
	      isn[n].vtY = (f * (v[j].vtY - v[i].vtY)) + v[i].vtY;
	      ++n;
	    }
	  }
	}
      }
    }
  }
  else /* cvh->vtxType == WLZ_VERTEX_D3 */
  {
    WlzDVertex3	v[3];

    for(i = 0; i < 3; ++i)
    {
      v[i] = cvh->vertices.d3[cvh->faces[f3 + i]];
    }
    for(i = 0; i < 3; ++i)
    {
      if(fabs(v[i].vtZ - z) < eps)
      {
	isn[n].vtX = v[i].vtX;
	isn[n].vtY = v[i].vtY;
	++n;
      }
      else
      {
	j = (i + 1) % 3;
	if(fabs(v[j].vtZ - z) > eps)
	{
	  double d;

	  d = v[j].vtZ - v[i].vtZ;
	  if(fabs(d) > eps)
	  {
	    double f;

	    f = (double )(z - v[i].vtZ) / (double )(d);
	    if((f > -(eps)) && (f < (1.0 + eps)))
	    {
	      isn[n].vtX = (f * (v[j].vtX - v[i].vtX)) + v[i].vtX;
	      isn[n].vtY = (f * (v[j].vtY - v[i].vtY)) + v[i].vtY;
	      ++n;
	    }
	  }
	}
      }
    }
  }
  return(n);
}

/*!
* \return	New object with a 3D polyon domain.
* \ingroup	WlzConvexHull
* \brief	Creates a new object with a 3D polyon domain, a single
* 		convex polygon per plane. The given object and it's
* 		domain are assumed non-NULL and valid.
* \param	obj		Given object which must be a 3D spatial
* 				domain object.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static WlzObject 		*WlzObjToConvexPolygon3d(
  				  WlzObject *obj,
  				  WlzErrorNum *dstErr)
{
  int		p;
  WlzObject	*obj1,
  		*obj2,
		*polygon = NULL;
  WlzDomain	domain;
  WlzDomain	*domains,
  		*new_domains;
  WlzValues	values;
  WlzErrorNum	errNum = WLZ_ERR_NONE;


  if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_POLYGON,
				    obj->domain.p->plane1,
				    obj->domain.p->lastpl,
				    obj->domain.p->line1,
				    obj->domain.p->lastln,
				    obj->domain.p->kol1,
				    obj->domain.p->lastkl,
				    &errNum)) != NULL){
    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
    values.core = NULL;
    polygon = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values, NULL, NULL,
			  &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    domains = obj->domain.p->domains;
    new_domains = domain.p->domains;
    values.core = NULL;
    for(p = obj->domain.p->plane1; p <= obj->domain.p->lastpl; p++)
    {
      if((*domains).core)
      {
	obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, *domains, values,
			   NULL, NULL, NULL);
	if((obj2 = WlzObjToConvexPolygon(obj1, &errNum)) != NULL)
	{
	  *new_domains = WlzAssignDomain(obj2->domain, NULL);
	  WlzFreeObj(obj2);
	}
	else
	{
	  WlzFreeObj(polygon);
	  polygon = NULL;
	  break;
	}
	WlzFreeObj(obj1);
      }
      else
      {
	(*new_domains).core = NULL;
      }
      domains++;
      new_domains++;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(polygon);
}

