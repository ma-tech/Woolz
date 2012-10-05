#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSnapFit_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzSnapFit.c
* \author       Bill Hill
* \date         July 2004
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
* \brief	Functions to compute correspondences between a pair of
* 		objects using closest points.
* \ingroup	WlzTransform
*/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Computes correspondences between the given target and source
* 		objects, based only on closest points. These may be used
*		to snap-fit one object to another when the alignment of
*		the objects is known to be close. If the optional transform
*		is given then all distances are with respect to the transformed
*		source.
*
*		The algorithm used is:
*		<ol>
*		  <li>
*		  Get vertices from target and source object.
*		  </li>
* 		  <li>
*		  Select those target vertices such that for any pair of
*		  target vertices \f$\mathbf{v_i}\f$ and \f$\mathbf{v_j}\f$
*		  the separation distance
*		  \f$|\mathbf{v_i} - \mathbf{v_j}| < d_t \forall i,j\f$
*		  </li>
*		  <li>
*		  For each source vertex find the closest target vertex
*		  and add those source vertices for which
*		  \f$|\mathbf{v_{ti}} - \mathbf{v_{si}}| \leq d_c\f$ to
*		  a list of correspondences.
*		  </li>
*		  <li>
*		  Reject all correspondences for which
*		  \f$|\mathbf{v_i} - \mathbf{v_j}| < d_s \forall i,j\f$,
*		  where \f$\mathbf{v_i}\f$ and \f$\mathbf{v_j}\f$ are
*		  the corresponding source vertices, keeping those with
*		  the minimum target-source separation by preference.
*		  </li>
*		</ol>

* \param	tObj			Target object.
* \param	sObj			Source object.
* \param	tr			Initial affine transform for source
*					object.
* \param	vType			Type of vertices returned, which is
*					always either WLZ_VERTEX_D2 or
*					WLZ_VERTEX_D3.
* \param	dstNVtx			Destination pointer for the number of
*					target vertces.
* \param	dstTVtxP		Destination pointer for the target
* 					vertces.
* \param	dstSVtxP		Destination pointer for the target
* 					vertces.
* \param	maxCDist		Maximum distance between target and
*					source vertex pairs \f$d_c\f$.
* \param	minTDist		Minimum distance between target
* 					vertces, \f$d_t\f$.
* \param	minSDist		Minimum distance between source
*					vertces \f$d_s\f$.
*/
WlzErrorNum	WlzSnapFit(WlzObject *tObj, WlzObject *sObj,
			   WlzAffineTransform *tr,
			   WlzVertexType *vType,
			   int *dstNVtx,
			   WlzVertexP *dstTVtxP, WlzVertexP *dstSVtxP,
			   double maxCDist, double minTDist, double minSDist)
{
  int		idM,
  		idN,
  		idT,
		idS = 0;
  double	tD0;
  WlzVertexP	tVP0,
  		tVtxP,
		sVtxP;
  int		*idxBuf = NULL;
  double	*dist = NULL;
  AlcKDTTree	*tree = NULL;
  AlcKDTNode	*node = NULL;
  int		nVtx[2],			/* Number of extracted
  						 * vertices with target 0 and
						 * source 1. */
		nCor[4],			/* Number of correspondences
						 * with: Target 0, closest
						 * source 1 and minimum
						 * distance closest source
						 * 2. */
  		dim[2];				/* Dimension of target 0 and
						 * source 1. */
  int		*idx[3];			/* Indices of corresponding
  						 * vertices in the extracted
						 * vertices target 0,
						 * closest source 1,
						 * minimum distance closest
						 * source 2 and minimum
						 * distance closest source with
						 * unique targets 3. */
  double	pos[3];
  WlzVertexP	vtxP[2];			/* Extracted target 0 and
  						 * source 1 vertices. */
  WlzVertexType	vtxType[2];			/* Extracted target 0 and
  						 * source 1 vertex type. */
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  tVtxP.v = sVtxP.v = NULL;
  vtxP[0].v = vtxP[1].v = NULL;
  idx[0] = idx[1] = idx[2] = NULL;
  nCor[0] = nCor[1] = nCor[2] = nCor[3] = 0;
  if((tObj == NULL) || (sObj == NULL))
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if((tObj->domain.core == NULL) || (sObj->domain.core == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((maxCDist < 0.0) || (minTDist < 0.0) || (minSDist < 0.0))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  /* Extract vertices from the objects. */
  if(errNum == WLZ_ERR_NONE)
  {
    vtxP[0] = WlzVerticesFromObj(tObj, NULL, nVtx + 0, vtxType + 0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vtxP[1] = WlzVerticesFromObj(sObj, NULL, nVtx + 1, vtxType + 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((idx[0] = (int *)AlcMalloc(nVtx[0] * sizeof(int))) == NULL) ||
       ((idx[1] = (int *)AlcMalloc(nVtx[1] * sizeof(int))) == NULL) ||
       ((idx[2] = (int *)AlcMalloc(nVtx[1] * sizeof(int))) == NULL) ||
       ((dist = (double *)AlcMalloc(nVtx[1] * sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Make sure that all the vertices are either WLZ_VERTEX_D2 or
   * WLZ_VERTEX_D3. */
  idN = 0;
  while((errNum == WLZ_ERR_NONE) && (idN < 2))
  {
    switch(vtxType[idN])
    {
      case WLZ_VERTEX_I2:
	if((tVP0.v = AlcMalloc(sizeof(WlzDVertex2) * nVtx[idN])) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  WlzValueCopyIVertexToDVertex(tVP0.d2, vtxP[idN].i2, nVtx[0]);
	  AlcFree(vtxP[0].v); vtxP[0].v = tVP0.v; tVP0.v = NULL;
	  vtxType[idN] = WLZ_VERTEX_D2;
	  dim[idN] = 2;
	}
	break;
      case WLZ_VERTEX_F2:
	if((tVP0.v = AlcMalloc(sizeof(WlzDVertex2) * nVtx[idN])) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  WlzValueCopyFVertexToDVertex(tVP0.d2, vtxP[idN].f2, nVtx[0]);
	  AlcFree(vtxP[0].v); vtxP[0].v = tVP0.v; tVP0.v = NULL;
	  vtxType[idN] = WLZ_VERTEX_D2;
	  dim[idN] = 2;
	}
        break;
      case WLZ_VERTEX_D2:
	dim[idN] = 3;
        break;
      case WLZ_VERTEX_I3:
	if((tVP0.v = AlcMalloc(sizeof(WlzDVertex3) * nVtx[idN])) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  WlzValueCopyIVertexToDVertex3(tVP0.d3, vtxP[idN].i3, nVtx[0]);
	  AlcFree(vtxP[0].v); vtxP[0].v = tVP0.v; tVP0.v = NULL;
	  vtxType[idN] = WLZ_VERTEX_D3;
	  dim[idN] = 3;
	}
	break;
      case WLZ_VERTEX_F3:
	if((tVP0.v = AlcMalloc(sizeof(WlzDVertex3) * nVtx[idN])) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  WlzValueCopyIVertexToDVertex3(tVP0.d3, vtxP[idN].i3, nVtx[0]);
	  AlcFree(vtxP[0].v); vtxP[0].v = tVP0.v; tVP0.v = NULL;
	  vtxType[idN] = WLZ_VERTEX_D3;
	  dim[idN] = 3;
	}
        break;
      case WLZ_VERTEX_D3:
	dim[idN] = 3;
        break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
    ++idN;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((dim[0] != dim[1]) ||
       (tr && (WlzAffineTransformDimension(tr, &errNum) != dim[0]) &&
        (errNum == WLZ_ERR_NONE)))
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  /* Create a kD-tree and populate it with the target vertices such that the
   * minimum target vertex seperation is >= minTDist. Index of nodes in the
   * kD-tree is that of the extracted target vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((idxBuf = (int *)AlcMalloc(sizeof(int) *
    				  WLZ_MAX(nVtx[0], nVtx[1]))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      (void )AlgShuffleIdx(nVtx[0], idxBuf, 0);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((tree = AlcKDTTreeNew(ALC_POINTTYPE_DBL, dim[0], 1.0E-06, nVtx[0],
    			     NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idT = idN = 0;
    while((alcErr == ALC_ER_NONE) && (idN < nVtx[0]))
    {
      idM = idxBuf[idN];
      if(dim[0] == 2)
      {
	tVP0.d2 = vtxP[0].d2 + idM;
	pos[0] = tVP0.d2->vtX;
	pos[1] = tVP0.d2->vtY;
      }
      else /* dim[0] == 3 */
      {
	tVP0.d3 = vtxP[0].d3 + idM;
	pos[0] = tVP0.d3->vtX;
	pos[1] = tVP0.d3->vtY;
	pos[2] = tVP0.d3->vtZ;
      }
      node = AlcKDTGetNN(tree, pos, DBL_MAX, &tD0, NULL);
      if((node == NULL) || (tD0 > minTDist))
      {
	if((node = AlcKDTInsert(tree, pos, NULL, &alcErr)) != NULL)
	{
	  *(idx[0] + idT++) = node->idx = idM;
	}
      }
      ++idN;
    }
    nCor[0] = idT;
    if(alcErr != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Foreach source vertex transform it and then find the closest vertex
   * in target kD-tree clollecting those which are closer than the minimum
   * distance. Then sort these correspondences by distance. */
  if(errNum == WLZ_ERR_NONE)
  {
    idS = idN = 0;
    while(idN < nVtx[1])
    {
      if(dim[1] == 2)
      {
	tVP0.d2 = vtxP[1].d2 + idN;
	*(tVP0.d2) = WlzAffineTransformVertexD2(tr, *(tVP0.d2), NULL);
	pos[0] = tVP0.d2->vtX;
	pos[1] = tVP0.d2->vtY;
      }
      else /* dim[1] == 3 */
      {
        tVP0.d3 = vtxP[1].d3 + idN;
	if(tr)
	{
	  *(tVP0.d3) = WlzAffineTransformVertexD3(tr, *(tVP0.d3), NULL);
	}
	pos[0] = tVP0.d3->vtX;
	pos[1] = tVP0.d3->vtY;
	pos[2] = tVP0.d3->vtZ;
      }
      node = AlcKDTGetNN(tree, pos, DBL_MAX, &tD0, NULL);
      if(node && (tD0 < maxCDist))
      {
        *(idx[1] + idS) = node->idx;
	idxBuf[idS] = idS;
	dist[idS++] = tD0;
      }
      ++idN;
    }
    nCor[1] = idS;
    (void )AlgHeapSortIdx(dist, idxBuf, (unsigned )(nCor[1]),
    			  AlgHeapSortCmpIdxDFn);
    /* Destroy the target kD-tree and then rebuild to ensure that no
     * source vertices are closer than minSDist. */
    (void )AlcKDTTreeFree(tree);
    if((tree = AlcKDTTreeNew(ALC_POINTTYPE_DBL, dim[1], 1.0E-06, nCor[1],
    			     NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Ensure that all source vertices have a seperation distance which
   * is >= minSDist, keeping those with the minimum target-source distance
   * by preference. */
  if(errNum == WLZ_ERR_NONE)
  {
    idN = idS = 0;
    while((alcErr == ALC_ER_NONE) && (idN < nCor[1]))
    {
      idM = idxBuf[idN];
      if(dim[1] == 2)
      {
	tVP0.d2 = vtxP[1].d2 + idM;
	if(tr)
	{
	  *(tVP0.d2) = WlzAffineTransformVertexD2(tr, *(tVP0.d2), NULL);
	}
	pos[0] = tVP0.d2->vtX;
	pos[1] = tVP0.d2->vtY;
      }
      else /* dim[1] == 3 */
      {
        tVP0.d3 = vtxP[1].d3 + idM;
	if(tr)
	{
	  *(tVP0.d3) = WlzAffineTransformVertexD3(tr, *(tVP0.d3), NULL);
	}
	pos[0] = tVP0.d3->vtX;
	pos[1] = tVP0.d3->vtY;
	pos[2] = tVP0.d3->vtZ;
      }
      node = AlcKDTGetNN(tree, pos, DBL_MAX, &tD0, NULL);
      if((node == NULL) || (tD0 > minSDist))
      {
	if(AlcKDTInsert(tree, pos, NULL, &alcErr) != NULL)
	{
	  *(idx[2] + idS) = idM;
	  ++idS;
	}
      }
      ++idN;
    }
    if(alcErr != ALC_ER_NONE)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    nCor[2] = idS;
  }
  (void )AlcKDTTreeFree(tree); tree = NULL;
  /* Allocate and set correspondences while making sure that target vertices
   * aren't multiply included. */
  if(errNum == WLZ_ERR_NONE)
  {
    idM = 0;
    WlzValueSetInt(idx[0], 0, nVtx[0]);
    if(dim[0] == 2)
    {
      if(((tVtxP.d2 = AlcMalloc(sizeof(WlzDVertex2) * nCor[2])) == NULL) ||
         ((sVtxP.d2 = AlcMalloc(sizeof(WlzDVertex2) * nCor[2])) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        for(idN = 0; idN < nCor[2]; ++idN)
	{
	  idT = *(idx[1] + idS);
	  idS = *(idx[2] + idN);
	  if(*(idx[0] + idT) == 0)
	  {
	    *(idx[0] + idT) = 1;
	    *(tVtxP.d2 + idM) = *(vtxP[0].d2 + idT);
	    *(sVtxP.d2 + idM) = *(vtxP[1].d2 + idS);
	    ++idM;
	  }
	}
      }
    }
    else /* dim[0] == 3 */
    {
      if(((tVtxP.d3 = AlcMalloc(sizeof(WlzDVertex3) * nCor[2])) == NULL) ||
         ((sVtxP.d3 = AlcMalloc(sizeof(WlzDVertex3) * nCor[2])) == NULL))
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
      else
      {
        for(idN = 0; idN < nCor[2]; ++idN)
	{
	  idS = *(idx[2] + idN);
	  idT = *(idx[1] + idS);
	  if(*(idx[0] + idT) == 0)
	  {
	    *(idx[0] + idT) = 1;
	    *(tVtxP.d3 + idM) = *(vtxP[0].d3 + idT);
	    *(sVtxP.d3 + idM) = *(vtxP[1].d3 + idS);
	    ++idM;
	  }
        }
      }
    }
    nCor[3] = idM;
  }
  /* Free stuff. */
  AlcFree(idxBuf);
  AlcFree(dist);
  AlcFree(idx[0]); AlcFree(idx[1]); AlcFree(idx[2]);
  AlcFree(vtxP[0].v); AlcFree(vtxP[1].v);
  /* Set return values. */
  if(errNum == WLZ_ERR_NONE)
  {
    *vType = vtxType[0];
    *dstNVtx = nCor[3];
    *dstTVtxP = tVtxP;
    *dstSVtxP = sVtxP;
  }
  return(errNum);
}
