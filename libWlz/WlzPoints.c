#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPoints_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzPoints.c
* \author       Bill Hill
* \date         December 2005
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
* \brief	Functions for handling point domains.
* \ingroup	WlzFeatures
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <Wlz.h>

/*!
* \return	New domain object which coresponds to the union of
* 		the given points.
* \ingroup	WlzFeatures
* \brief	Creates a domain object which coresponds to the union of
* 		the given points.
* \param	pnt			Point domain.
* \param	scale			Scale, which if greater than zero
* 					is used as the diameter of a circle
* 					or sphere centred on each of the
* 					points vertices and a multiplier
* 					for the point position.
* \param	dstErr			Destination error poiter, may be NULL.
*/
WlzObject	*WlzPointsToDomObj(WlzPoints *pnt, double scale,
    				   WlzErrorNum *dstErr)
{
  int		idP;
  WlzObjectType	dType;
  WlzObject	*tObj0 = NULL,
		*tObj1 = NULL,
		*tObj2 = NULL,
		*dObj = NULL;
  WlzVertex	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  idP = 0;
  if(scale > DBL_EPSILON)
  {
    pos.d3.vtZ = 0.0;
  }
  else
  {
    pos.i3.vtZ = 0;
  }
  if(pnt == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pnt->type)
    {
      case WLZ_POINTS_2I: /* FALLTHROUGH */
      case WLZ_POINTS_2D:
	dType = WLZ_2D_DOMAINOBJ;
	break;
      case WLZ_POINTS_3I: /* FALLTHROUGH */
      case WLZ_POINTS_3D:
	dType = WLZ_3D_DOMAINOBJ;
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tObj0 = WlzMakeEmpty(&errNum);
  }
  while((errNum == WLZ_ERR_NONE) && (idP < pnt->nPoints))
  {
    if(scale > DBL_EPSILON)
    {
      switch(pnt->type)
      {
	case WLZ_POINTS_2I:
	  pos.d3.vtX = pnt->points.i2[idP].vtX;
	  pos.d3.vtY = pnt->points.i2[idP].vtY;
	  break;
	case WLZ_POINTS_2D:
	  pos.d3.vtX = pnt->points.d2[idP].vtX;
	  pos.d3.vtY = pnt->points.d2[idP].vtY;
	  pos.d3.vtZ = 0.0;
	  break;
	case WLZ_POINTS_3I:
	  pos.d3.vtX = pnt->points.i3[idP].vtX;
	  pos.d3.vtY = pnt->points.i3[idP].vtY;
	  pos.d3.vtZ = pnt->points.i3[idP].vtZ;
	  break;
	case WLZ_POINTS_3D:
	  pos.d3 = pnt->points.d3[idP];
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj1 = WlzMakeSphereObject(dType, scale / 2.0,
				    scale * pos.d3.vtX,
				    scale * pos.d3.vtY,
				    scale * pos.d3.vtZ,
				    &errNum);
      }
    }
    else
    {
      switch(pnt->type)
      {
	case WLZ_POINTS_2I:
	  pos.i3.vtX = pnt->points.i2[idP].vtX;
	  pos.i3.vtY = pnt->points.i2[idP].vtY;
	  break;
	case WLZ_POINTS_2D:
	  pos.i3.vtX = WLZ_NINT(pnt->points.d2[idP].vtX);
	  pos.i3.vtY = WLZ_NINT(pnt->points.d2[idP].vtY);
	  break;
	case WLZ_POINTS_3I:
	  pos.i3 = pnt->points.i3[idP];
	  break;
	case WLZ_POINTS_3D:
	  pos.i3.vtX = WLZ_NINT(pnt->points.d3[idP].vtX);
	  pos.i3.vtY = WLZ_NINT(pnt->points.d3[idP].vtY);
	  pos.i3.vtZ = WLZ_NINT(pnt->points.d3[idP].vtZ);
	  break;
        default:
	  errNum = WLZ_ERR_PARAM_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	tObj1 = WlzMakeSinglePixelObject(dType,
					 pos.i3.vtX, pos.i3.vtY, pos.i3.vtZ,
				         &errNum);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      tObj2 = WlzUnion2(tObj0, tObj1, &errNum);
    }
    (void )WlzFreeObj(tObj0);
    (void )WlzFreeObj(tObj1);
    tObj0 = tObj2;
    ++idP;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dObj = tObj0;
  }
  else
  {
    (void )WlzFreeObj(tObj0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dObj);
}

/*!
* \return	Pointer to the points value.
* \ingroup	WlzValueUtils
* \brief	Computes the pointer to the individual value corresponding
* 		to the given index.
* \param	pts			Given points values.
* \param	idx			Given index.
*/
void		*WlzPointValueGet(WlzPointValues *pts, int idx)
{
  void		*val;

  val = (char *)pts->values.v + (idx * pts->pSz);
  return(val);
}

/*!
* \return	New points domain or NULL on error.
* \ingroup	WlzFeatures
* \brief	Finds points which are within the given opjects domain
* 		and are seperated by at least the given minimum distance.
* \param	gvnObj			Given spatial domain object.
* \param	dMin			Given minimum distance (if less than
* 					1.0 then will be set to 1.0)..
* \param	voxelScaling		Use voxel scaling if non-zero.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPoints			*WlzPointsFromDomObj(
				  WlzObject *gvnObj,
				  double dMin,
				  int voxelScaling,
				  WlzErrorNum *dstErr)
{
  int		dim = 0,
		vecIdx = 0,
  		nShfBuf = 0;
  int		*shfBuf = NULL;
  AlcKDTTree	*tree = NULL;
  WlzPoints	*pts = NULL;
  WlzObject	*curObj = NULL,
  		*strObj = NULL;;
  WlzDVertex3	voxSz;
  AlcVector	*vec = NULL;
  size_t	vtxSz = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  voxSz.vtX = voxSz.vtY = voxSz.vtZ = 1.0;
  if(gvnObj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else
  {
    switch(gvnObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        dim = 2;
	vtxSz = sizeof(WlzIVertex2);
	break;
      case WLZ_3D_DOMAINOBJ:
        dim = 3;
	vtxSz = sizeof(WlzIVertex3);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if(dim)
    {
      if(gvnObj->domain.core == NULL)
      {
        errNum = WLZ_ERR_DOMAIN_NULL;
      }
      if(dim == 2)
      {
        switch(gvnObj->domain.core->type)
	{
	  case WLZ_INTERVALDOMAIN_INTVL:
	  case WLZ_INTERVALDOMAIN_RECT:
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
      }
      else /* dim == 3 */
      {
        switch(gvnObj->domain.core->type)
	{
	  case WLZ_PLANEDOMAIN_DOMAIN:
	    if(voxelScaling)
	    {
	      voxSz.vtX = gvnObj->domain.p->voxel_size[0];
	      voxSz.vtY = gvnObj->domain.p->voxel_size[1];
	      voxSz.vtZ = gvnObj->domain.p->voxel_size[2];
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
      }
    }
  }
  /* Create a vector in which to collect vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    vec = AlcVectorNew(1024, vtxSz, 1024, NULL);
    if(vec == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Create a KD-tree for proximity queries. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(dMin < 1.0)
    {
      dMin = 1.0;
    }
    tree = AlcKDTTreeNew(ALC_POINTTYPE_INT, dim, 0.0, 0, NULL);
    if(tree == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Create a structuring element with which to successively erode the given
   * object. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		rad;

    rad = (dMin < 2.0)? 1: ALG_NINT(dMin - 1.0);
    if(dim == 2)
    {
      strObj = WlzAssignObject(
      	       WlzMakeRectangleObject(rad, rad, 0.0, 0.0, &errNum), NULL);
    }
    else
    {
      strObj = WlzAssignObject(
      	       WlzMakeCuboidObject(WLZ_3D_DOMAINOBJ, rad, rad, rad,
				   0.0, 0.0, 0.0, &errNum), NULL);
    }
  }
  /* Create a current object, initially this is just the given object's
   * domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzValues	nullVal;

    nullVal.core = NULL;
    curObj = WlzAssignObject(
    	     WlzMakeMain(gvnObj->type, gvnObj->domain, nullVal,
	                 NULL, NULL, &errNum), NULL);
  }
  /* Erode the current object while adding points from the boundary until
   * the current object has been eroded away. */
  while((errNum == WLZ_ERR_NONE) && (curObj != NULL))
  {
    int		nVtx = 0;
    WlzVertexP	vtx;
    WlzObject	*shlObj;

    vtx.v = NULL;
    /* Get the vertices that lie on the boundry (shell) domain of the
     * current object. */
    shlObj = WlzAssignObject(
             WlzBoundaryDomain(curObj, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      if(dim == 2)
      {
        errNum = WlzVerticesFromObj2I(shlObj, &nVtx, &(vtx.i2));
      }
      else
      {
        errNum = WlzVerticesFromObj3I(shlObj, &nVtx, &(vtx.i3));
      }
    }
    (void )WlzFreeObj(shlObj);
    /* Update the shuffle buffer, used to make randomise the insertions into
     * the KD-tree. */
    if((errNum == WLZ_ERR_NONE) && (nShfBuf < nVtx))
    {
      nShfBuf = nVtx + 1024;
      if((shfBuf = (int *)AlcRealloc(shfBuf, sizeof(int) * nShfBuf)) == NULL)
      {
	nShfBuf = 0;
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )AlgShuffleIdx(nVtx, shfBuf, 0);
      if(AlcVectorExtend(vec, vecIdx + nVtx + 1) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    /* Add vertices to the KD-tree as long as their seperation distance is
     * >= dMin. */
    if(errNum == WLZ_ERR_NONE)
    {
      int	i;

      for(i = 0; i < nVtx; ++i)
      {
        int 	j;
	double	d;
	int	pos[3];
	AlcKDTNode *nod;

	j = shfBuf[i];
	if(dim == 2)
	{
	  pos[0] = vtx.i2[j].vtX;
	  pos[1] = vtx.i2[j].vtY;
	}
	else
	{
	  pos[0] = vtx.i3[j].vtX;
	  pos[1] = vtx.i3[j].vtY;
	  pos[2] = vtx.i3[j].vtZ;
	}
	nod = AlcKDTGetNN(tree, pos, DBL_MAX, &d, NULL);
	if((nod == NULL) || (d > dMin))
	{
	  AlcErrno alcErr = ALC_ER_NONE;

	  (void )AlcKDTInsert(tree, pos, NULL, &alcErr);
	  if(alcErr == ALC_ER_NONE)
	  {
	    WlzVertexP p;

	    p.v = AlcVectorItemGet(vec, vecIdx);
	    if(dim == 2)                               
	    {
	      *(p.i2) = vtx.i2[j];
	    }
	    else
	    {
	      *(p.i3) = vtx.i3[j];
	    }
	    ++vecIdx;
	  }
	  else
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	}
      }
    }
    AlcFree(vtx.v);
    /* Erode the current object. */
    if(errNum == WLZ_ERR_NONE)
    {
      WlzObject	*tObj;

      tObj = WlzAssignObject(WlzStructErosion(curObj, strObj, &errNum), NULL);
      (void )WlzFreeObj(curObj);
      curObj = tObj;
    }
    /* Check for completion when the current object will have been eroded
     * away. */
    if(errNum == WLZ_ERR_NONE)
    {
      if(WlzIsEmpty(curObj, &errNum))
      {
        (void )WlzFreeObj(curObj);
	curObj = NULL;
      }
    }
  }
  AlcFree(shfBuf);
  (void )WlzFreeObj(curObj);
  (void )WlzFreeObj(strObj);
  (void )AlcKDTTreeFree(tree);
  /* Construct the points using the vector of vertex locations. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzVertexP nullP;
    WlzObjectType pntType;

    nullP.v = NULL;
    if(dim == 2)
    {
      pntType = WLZ_POINTS_2I;
    }
    else 
    {
      if(voxelScaling)
      {
        pntType = WLZ_POINTS_3D;
      }
      else
      {
        pntType = WLZ_POINTS_3I;
      }
    }
    pts = WlzMakePoints(pntType, 0, nullP, vecIdx, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;
    WlzVertexP 	p;

    pts->nPoints = vecIdx;
    if(dim == 2)
    {
      for(i = 0; i < vecIdx; ++i)
      {
        p.v = AlcVectorItemGet(vec, i);
        pts->points.i2[i] = *(p.i2);
      }
    }
    else
    {
      if(voxelScaling)
      {
	for(i = 0; i < vecIdx; ++i)
	{
	  p.v = AlcVectorItemGet(vec, i);
	  pts->points.d3[i].vtX = p.i3->vtX * voxSz.vtX;
	  pts->points.d3[i].vtY = p.i3->vtY * voxSz.vtY;
	  pts->points.d3[i].vtZ = p.i3->vtZ * voxSz.vtZ;
	}
      }
      else
      {
	for(i = 0; i < vecIdx; ++i)
	{
	  p.v = AlcVectorItemGet(vec, i);
	  pts->points.i3[i] = *(p.i3);
	}
      }
    }
  }
  (void )AlcVectorFree(vec);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pts);
}

/*!
* \return       New domain object without values.
* \ingroup	WlzFeatures
* \brief        Constructs a domain from the union of marker domains with
*               a marker domain at each of the given point positions.
* \param	pts			Given points.
* \param        mType                   Marker type.
* \param        mSz                     Marker size. This is the radius of a
*					sphere marker. The marker size is
*					ignored for point markers.
* \param        dstErr                  Destination error pointer, may be NULL.
*/
WlzObject			*WlzPointsToMarkers(
				  WlzPoints *pts,
				  WlzMarkerType mType,
				  int mSz,
				  WlzErrorNum *dstErr)
{
  WlzVertexType vType = WLZ_VERTEX_ERROR;
  WlzObject	*mObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(pts == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pts->type)
    {
      case WLZ_POINTS_2I:
        vType = WLZ_VERTEX_I2;
	break;
      case WLZ_POINTS_2D:
        vType = WLZ_VERTEX_D2;
	break;
      case WLZ_POINTS_3I:
        vType = WLZ_VERTEX_I3;
	break;
      case WLZ_POINTS_3D:
        vType = WLZ_VERTEX_D3;
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mObj = WlzMakeMarkers(vType, pts->nPoints, pts->points, mType, mSz,
    			  &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(mObj);
}

/*!
* \return	New points domain.
* \ingroup	WlzFeatures
* \brief	Dithers the vertices of the given points by adding a
* 		random displacement within the range -d to +d, where
* 		d is the given dither size.
* \param	gPts			Given points.
* \param	dSz			Dither size (z component is not
* 					used if points are 2D).
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzPoints			*WlzPointsDither(
				  WlzPoints *gPts,
				  WlzDVertex3 dSz,
				  WlzErrorNum *dstErr)
{
  WlzPoints	*dPts = NULL;
  WlzObjectType dPntType = WLZ_POINTS_2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gPts == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(gPts->type)
    {
      case WLZ_POINTS_2I:        
      case WLZ_POINTS_2D:        
        dPntType = WLZ_POINTS_2D;
	break;
      case WLZ_POINTS_3I:        
      case WLZ_POINTS_3D:        
        dPntType = WLZ_POINTS_3D;
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzVertexP nullP;

    nullP.v = NULL;
    dPts = WlzMakePoints(dPntType, 0, nullP, gPts->nPoints, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)                
  {
    int		idx;

    dPts->nPoints = gPts->nPoints;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idx = 0; idx < gPts->nPoints; ++idx)
    {
      WlzDVertex3 d;

      d.vtX = dSz.vtX * ((AlgRandUniform() * 2.0) - 1.0);
      d.vtY = dSz.vtY * ((AlgRandUniform() * 2.0) - 1.0);
      switch(gPts->type)
      {
	case WLZ_POINTS_2I:        
	  dPts->points.d2[idx].vtX = gPts->points.i2[idx].vtX + d.vtX;
	  dPts->points.d2[idx].vtY = gPts->points.i2[idx].vtY + d.vtY;
	  break;
	case WLZ_POINTS_2D:        
	  dPts->points.d2[idx].vtX = gPts->points.d2[idx].vtX + d.vtX;
	  dPts->points.d2[idx].vtY = gPts->points.d2[idx].vtY + d.vtY;
	  break;
	case WLZ_POINTS_3I:        
          d.vtZ = dSz.vtZ * ((AlgRandUniform() * 2.0) - 1.0);
	  dPts->points.d3[idx].vtX = gPts->points.i3[idx].vtX + d.vtX;
	  dPts->points.d3[idx].vtY = gPts->points.i3[idx].vtY + d.vtY;
	  dPts->points.d3[idx].vtZ = gPts->points.i3[idx].vtZ + d.vtZ;
	  break;
	case WLZ_POINTS_3D:        
          d.vtZ = dSz.vtZ * ((AlgRandUniform() * 2.0) - 1.0);
	  dPts->points.d3[idx].vtX = gPts->points.d3[idx].vtX + d.vtX;
	  dPts->points.d3[idx].vtY = gPts->points.d3[idx].vtY + d.vtY;
	  dPts->points.d3[idx].vtZ = gPts->points.d3[idx].vtZ + d.vtZ;
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
  return(dPts);
}

