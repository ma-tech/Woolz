#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshUtils_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzCMeshUtils.c
* \author       Bill Hill
* \date         June 2003
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
* \brief	Utility functions for 2D and 3D conforming simplical meshes.
*
* \ingroup	WlzMesh
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

static void			WlzCMeshFilterLPL2D(
				  WlzCMesh2D *mesh,
				  WlzDVertex2 *vGIn,
				  WlzDVertex2 *vGOut,
				  double lambda,
				  int doBnd);
static void			WlzCMeshFilterLPL3D(
				  WlzCMesh3D *mesh,
				  WlzDVertex3 *vGIn,
				  WlzDVertex3 *vGOut,
				  double lambda,
				  int doBnd);
static WlzDVertex2 		WlzCMeshFilterLPLDelta2D(
				  WlzCMesh2D *mesh,
				  WlzCMeshNod2D *nod,
				  WlzDVertex2 *vBuf,
				  int doBnd);
static WlzDVertex3 		WlzCMeshFilterLPLDelta3D(
				  WlzCMesh3D *mesh,
				  WlzCMeshNod3D *nod,
				  WlzDVertex3 *vBuf,
				  int doBnd);

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Computes the mesh maximum edge length which is used to
*               terminate vertex location. This should not be allowed
*               to become less than the actual maximum edge length or
*               vertex location may fail, also if it is far larger than
*               the actual maximum edge length then vertex location
*               will be inefficient when vertices are outside the mesh.
* \param        mesh                    The mesh.
*/
void     	WlzCMeshUpdateMaxSqEdgLen2D(WlzCMesh2D *mesh)
{
  int           idE;
  double        dSq;
  WlzCMeshElm2D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    mesh->maxSqEdgLen = 0.0;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec,
      				              idE);
      if(elm->idx >= 0)
      {
	dSq = WlzGeomDistSq2D(elm->edg[0].nod->pos,
	                      elm->edg[1].nod->pos);
	if(dSq > mesh->maxSqEdgLen)
	{
	  mesh->maxSqEdgLen = dSq;
	}
	dSq = WlzGeomDistSq2D(elm->edg[1].nod->pos,
			      elm->edg[2].nod->pos);
	if(dSq > mesh->maxSqEdgLen)
	{
	  mesh->maxSqEdgLen = dSq;
	}
	dSq = WlzGeomDistSq2D(elm->edg[2].nod->pos,
	                      elm->edg[0].nod->pos);
	if(dSq > mesh->maxSqEdgLen)
	{
	  mesh->maxSqEdgLen = dSq;
	}
      }
    }
  }
}

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Computes the mesh maximum edge length which is used to
*               terminate vertex location. This should not be allowed
*               to become less than the actual maximum edge length or
*               vertex location may fail, also if it is far larger than
*               the actual maximum edge length then vertex location
*               will be inefficient when vertices are outside the mesh.
* \param        mesh                    The mesh.
*/
void            WlzCMeshUpdateMaxSqEdgLen3D(WlzCMesh3D *mesh)
{
  int           idE,
                idF;
  double        dSq;
  WlzCMeshElm3D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    mesh->maxSqEdgLen = 0.0;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec,
      					      idE);
      if(elm->idx >= 0)
      {
	for(idF = 0; idF < 3; ++idF) /* Only need to check 3 faces
	                                for all edges. */
	{
	  dSq = WlzGeomDistSq3D(elm->face[idF].edg[0].nod->pos,
				elm->face[idF].edg[1].nod->pos);
	  if(dSq > mesh->maxSqEdgLen)
	  {
	    mesh->maxSqEdgLen = dSq;
	  }
	  dSq = WlzGeomDistSq3D(elm->face[idF].edg[1].nod->pos,
				elm->face[idF].edg[2].nod->pos);
	  if(dSq > mesh->maxSqEdgLen)
	  {
	    mesh->maxSqEdgLen = dSq;
	  }
	  dSq = WlzGeomDistSq3D(elm->face[idF].edg[2].nod->pos,
				elm->face[idF].edg[0].nod->pos);
	  if(dSq > mesh->maxSqEdgLen)
	  {
	    mesh->maxSqEdgLen = dSq;
	  }
	}
      }
    }
  }
}

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Updates the bounding box of the 2D conforming mesh.
* \param        mesh                    The mesh.
*/
void     	WlzCMeshUpdateBBox2D(WlzCMesh2D *mesh)
{
  int           idN,
		firstNod;
  WlzCMeshNod2D *nod;
  WlzDBox2      bBox;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    /* Update the bounding box. */
    firstNod = 1;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	if(firstNod)
	{
	  firstNod = 0;
	  bBox.xMin = bBox.xMax = nod->pos.vtX;
	  bBox.yMin = bBox.yMax = nod->pos.vtY;
	}
	else
	{
	  if(nod->pos.vtX < bBox.xMin)
	  {
	    bBox.xMin = nod->pos.vtX;
	  }
	  else if(nod->pos.vtX > bBox.xMax)
	  {
	    bBox.xMax = nod->pos.vtX;
	  }
	  if(nod->pos.vtY < bBox.yMin)
	  {
	    bBox.yMin = nod->pos.vtY;
	  }
	  else if(nod->pos.vtY > bBox.yMax)
	  {
	    bBox.yMax = nod->pos.vtY;
	  }
	}
      }
    }
    if(firstNod == 0)
    {
      mesh->bBox = bBox;
    }
  }
}

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Updates the bounding box of the 3D conforming mesh.
* \param        mesh                    The mesh.
*/
void            WlzCMeshUpdateBBox3D(WlzCMesh3D *mesh)
{
  int           idN,
                firstNod;
  WlzCMeshNod3D *nod;
  WlzDBox3      bBox;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    /* Update the bounding box. */
    firstNod = 1;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	if(firstNod)
	{
	  firstNod = 0;
	  bBox.xMin = bBox.xMax = nod->pos.vtX;
	  bBox.yMin = bBox.yMax = nod->pos.vtY;
	  bBox.zMin = bBox.zMax = nod->pos.vtZ;
	}
	else
	{
	  if(nod->pos.vtX < bBox.xMin)
	  {
	    bBox.xMin = nod->pos.vtX;
	  }
	  else if(nod->pos.vtX > bBox.xMax)
	  {
	    bBox.xMax = nod->pos.vtX;
	  }
	  if(nod->pos.vtY < bBox.yMin)
	  {
	    bBox.yMin = nod->pos.vtY;
	  }
	  else if(nod->pos.vtY > bBox.yMax)
	  {
	    bBox.yMax = nod->pos.vtY;
	  }
	  if(nod->pos.vtZ < bBox.zMin)
	  {
	    bBox.zMin = nod->pos.vtZ;
	  }
	  else if(nod->pos.vtZ > bBox.zMax)
	  {
	    bBox.zMax = nod->pos.vtZ;
	  }
	}
      }
    }
    if(firstNod == 0)
    {
      mesh->bBox = bBox;
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the mesh.
* \param	mesh			Given mesh.
*/
void		WlzCMeshSetBoundNodFlags(WlzCMeshP mesh)
{
  int		idN;
  WlzCMeshNod2D *nod;

  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        WlzCMeshSetBoundNodFlags2D(mesh.m2);
	break;
      case WLZ_CMESH_TET3D:
        WlzCMeshSetBoundNodFlags3D(mesh.m3);
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the 2D mesh.
* \param	mesh			Given mesh.
*/
void		WlzCMeshSetBoundNodFlags2D(WlzCMesh2D *mesh)
{
  int		idN;
  WlzCMeshNod2D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	nod->flags &= ~(WLZ_CMESH_NOD_FLAG_BOUNDARY);
	if(WlzCMeshNodIsBoundary2D(nod))
	{
	  nod->flags |= WLZ_CMESH_NOD_FLAG_BOUNDARY;
	}
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the 3D mesh.
* \param	mesh			Given mesh.
*/
void		WlzCMeshSetBoundNodFlags3D(WlzCMesh3D *mesh)
{
  int		idN;
  WlzCMeshNod3D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	nod->flags &= ~(WLZ_CMESH_NOD_FLAG_BOUNDARY);
	if(WlzCMeshNodIsBoundary3D(nod))
	{
	  nod->flags |= WLZ_CMESH_NOD_FLAG_BOUNDARY;
	}
      }
    }
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief
* \param	mesh			Given 3D mesh.
* \param	dstNBndElm		Destination pointer for the number
*					of boundary elements.
* \param	dstBndVec		Destination pointer for a vector
*					of boundary element indices.
* \param	trustNndFlags		If non-zero the element boundary
*					flags can be trusted.
*/
WlzErrorNum	WlzCMeshGetBndElm3D(WlzCMesh3D *mesh,
				    int *dstNBndElm,  AlcVector **dstBndVec,
				    int trustNndFlags)
{
  int		idB,
  		idE,
		idF,
		isBnd;
  int		*idxP;
  AlcVector	*bndElm = NULL;
  WlzCMeshElm3D	*elm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if (mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((bndElm = AlcVectorNew(1, sizeof(int),
                              mesh->res.elm.vec->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idE = 0;
    idB = 0;
    while(idE < mesh->res.elm.maxEnt)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	/* Establish whether the element is a boundary element. */
	isBnd = 0;
	if(trustNndFlags)
	{
	  isBnd = elm->flags & WLZ_CMESH_ELM_FLAG_BOUNDARY;
	}
	else
	{
	  for(idF = 0; idF < 4; ++idF)
	  {
            if((elm->face[idF].opp == NULL) ||
	       (elm->face[idF].opp == &(elm->face[idF])))
	    {
	      isBnd = 1;
	      break;
	    }
	  }
	}
	/* If element is a boundary element then record it's index in the
	 * vector of boundary elements. */
	if(isBnd)
	{
	  if((idxP = (int *)AlcVectorExtendAndGet(bndElm, idB)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	  *idxP = elm->idx;
	  ++idB;
	}
      }
      ++idE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNBndElm = idB;
    *dstBndVec = bndElm;
  }
  else
  {
    (void )AlcVectorFree(bndElm);
  }
  return(errNum);
}

/*!
* \return				Non-zero if the node is a boundary
* 					node.
* \ingroup	WlzMesh
* \brief	Checks whether the node is a boundary node by examining
* 		the edges which use the node. If one of these edges does
*		not have an opposite edge (other than itself) the node is
*		a boundary node.
* \param	nod			Given node of mesh.
*/
int		WlzCMeshNodIsBoundary2D(WlzCMeshNod2D *nod)
{
  int		isBnd = 0;
  WlzCMeshEdg2D	*edg0,
  		*edg1;

  if(nod && (nod->idx >= 0))
  {
    edg0 = edg1 = nod->edg;
    do
    {
      if((edg1->opp == NULL) || (edg1->opp == edg1))
      {
	isBnd = 1;
	break;
      }
      edg1 = edg1->nnxt;
    } while(edg1 != edg0);
  }
  return(isBnd);
}

/*!
* \return				Non-zero if the node is a boundary
* 					node.
* \ingroup	WlzMesh
* \brief	Checks whether the node is a boundary node by examining
* 		the faces which use the node. If one of these faces does
*		not have an opposite face (other than itself) the node is
*		a boundary node.
* \param	nod			Given node of mesh.
*/
int		WlzCMeshNodIsBoundary3D(WlzCMeshNod3D *nod)
{
  int		isBnd = 0;
  WlzCMeshFace	*fce;
  WlzCMeshEdg3D	*edg0,
  		*edg1;

  if(nod && (nod->idx >= 0))
  {
    edg0 = edg1 = nod->edg;
    do
    {
      fce = edg1->face;
      if((fce->opp == NULL) || (fce->opp == fce))
      {
	isBnd = 1;
	break;
      }
      edg1 = edg1->nnxt;
    } while(edg1 != edg0);
  }
  return(isBnd);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Applies a Laplacian smoothing to the 2D mesh in which
*		nodes are iteratively moved to the centroid of their
*		imediate neighbours. If a node is on a boundary of the
*		mesh then it is moved on the boundary, ie to the centroid
*		of it's neighboring boundary nodes.
*		Before calling this function all nodes must have had the
*		boundary node flag bit set or cleared appropriately.
*		This function will shrink meshes when the boundary
*		parameter is set.
*
*		Each node at position \f$p_i\f$ is moved to \f$p'_i\f$:
*		\f[
                    p'_i = (1 - \alpha)p_i +
		           \frac{\alpha}{n}\sum_{j}^{n}{p_{ij}}
		\f]
*		where \f$\alpha\f$ is the weight factor.
* \param	mesh			Given mesh.
* \param	itr			Number of iterations.
* \param	alpha			Weight factor.
* \param	doBnd			Apply smoothing to boundary nodes
*					if non-zero.
* \param	update			Update the mesh bucket grid and
*					maximum edge length.
*/
WlzErrorNum	WlzCMeshLaplacianSmooth(WlzCMeshP mesh,
					  int itr, double alpha,
					  int doBnd, int update)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        errNum = WlzCMeshLaplacianSmooth2D(mesh.m2, itr, alpha, doBnd, update);
	break;
      case WLZ_CMESH_TET3D:
        errNum = WlzCMeshLaplacianSmooth3D(mesh.m3, itr, alpha, doBnd, update);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Applies a Laplacian smoothing to the 2D mesh in which
*		nodes are iteratively moved to the centroid of their
*		imediate neighbours. See WlzCMeshLaplacianSmooth().
* \param	mesh			Given mesh.
* \param	itr			Number of iterations.
* \param	alpha			Weight factor.
* \param	doBnd			Apply smoothing to boundary nodes
*					if non-zero.
* \param	update			Update the mesh bucket grid and
*					maximum edge length.
*/
WlzErrorNum	WlzCMeshLaplacianSmooth2D(WlzCMesh2D *mesh,
					  int itr, double alpha,
					  int doBnd, int update)
{
  int		idI,
  		idN,
  		nCnt;
  WlzDVertex2	nPos;
  WlzCMeshNod2D *nod,
  		*oNod;
  WlzCMeshEdg2D	*edg0,
  		*edg1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    for(idI = 0; idI < itr; ++idI)
    {
      for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
      {
	nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	if(nod->idx >= 0)
	{
          if(doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) == 0))
	  {
	    nCnt = 0;
	    nPos.vtX = nPos.vtY = 0.0;
	    edg0 = edg1 = nod->edg;
	    do
	    {
	      oNod = edg1->next->nod;
	      nPos.vtX += oNod->pos.vtX;
	      nPos.vtY += oNod->pos.vtY;
	      ++nCnt;
	      edg1 = edg1->nnxt;
	    } while(edg0 != edg1);
	    if(nCnt > 0)
	    {
	      nod->pos.vtX = (1.0 - alpha) * nod->pos.vtX + 
			     alpha * nPos.vtX / nCnt;
	      nod->pos.vtY = (1.0 - alpha) * nod->pos.vtY + 
			     alpha * nPos.vtY / nCnt;
	    }
	  }
	}
      }
    }
    if(update)
    {
      WlzCMeshUpdateBBox2D(mesh);
      WlzCMeshUpdateMaxSqEdgLen2D(mesh);
      errNum = WlzCMeshReassignBuckets2D(mesh, 0);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Applies a Laplacian smoothing to the 3D mesh in which
*		nodes are iteratively moved to the centroid of their
*		imediate neighbours. See WlzCMeshLaplacianSmooth().
* \param	mesh			Given mesh.
* \param	itr			Number of iterations.
* \param	alpha			Weight factor.
* \param	doBnd			Apply smoothing to boundary nodes
*					if non-zero.
* \param	update			Update the mesh bucket grid and
*					maximum edge length.
*/
WlzErrorNum	WlzCMeshLaplacianSmooth3D(WlzCMesh3D *mesh,
					  int itr, double alpha,
					  int doBnd, int update)
{
  int		idI,
  		idN,
  		nCnt;
  WlzDVertex3	nPos;
  WlzCMeshNod3D *nod,
  		*oNod;
  WlzCMeshEdg3D	*edg0,
  		*edg1;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    for(idI = 0; idI < itr; ++idI)
    {
      for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
      {
	nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
	if(nod->idx >= 0)
	{
          if(doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) == 0))
	  {
	    nCnt = 0;
	    nPos.vtX = nPos.vtY = nPos.vtZ = 0.0;
	    edg0 = edg1 = nod->edg;
	    do
	    {
	      oNod = edg1->next->nod;
	      nPos.vtX += oNod->pos.vtX;
	      nPos.vtY += oNod->pos.vtY;
	      nPos.vtZ += oNod->pos.vtZ;
	      ++nCnt;
	      edg1 = edg1->nnxt;
	    } while(edg0 != edg1);
	    if(nCnt > 0)
	    {
	      nod->pos.vtX = (1.0 - alpha) * nod->pos.vtX + 
			     alpha * nPos.vtX / nCnt;
	      nod->pos.vtY = (1.0 - alpha) * nod->pos.vtY + 
			     alpha * nPos.vtY / nCnt;
	      nod->pos.vtZ = (1.0 - alpha) * nod->pos.vtZ + 
			     alpha * nPos.vtZ / nCnt;
	    }
	  }
	}
      }
    }
    if(update)
    {
      WlzCMeshUpdateBBox3D(mesh);
      WlzCMeshUpdateMaxSqEdgLen3D(mesh);
      errNum = WlzCMeshReassignBuckets3D(mesh, 0);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Applies a low pass filter to the mesh in which
*		nodes are moved. If a node is on a boundary of the
*		mesh then it is moved on the boundary.
*		Before calling this function all nodes must have had the
*		boundary node flag bit set or cleared appropriately.
*		This function should not significantly shrink meshes
*		because it applies a low pass filter.
* \param	mesh			Given mesh.
* \param	kPB			The band pass frequency parameter.
* \param	kSB			The band stop frequency parameter.
* \param	dPB			The pass band maximum deviation.
* \param	dSB			The stop band maximum deviation.
* \param	maxItr			Maximum number of iterations.
* \param	doBnd			Apply smoothing to boundary nodes
*					if non-zero.
* \param	update			Update the mesh bucket grid and
*					maximum edge length.
*/
WlzErrorNum	WlzCMeshLPFilter(WlzCMeshP mesh,
                                 double kPB, double kSB,
				 double dPB, double dSB,
				 int maxItr, int doBnd,
				 int update)
{
  int		nItr;
  double	lambda,
  		mu;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_NONE;
  }
  else
  {
    errNum = WlzGMFilterGeomLPParam(&lambda, &mu, &nItr, kPB, kSB, dPB, dSB);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nItr > maxItr)
    {
      nItr = maxItr;
    }
    errNum = WlzCMeshLPFilterLM(mesh, lambda, mu, nItr, doBnd, update);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(update)
    {
      switch(mesh.m2->type)
      {
        case WLZ_CMESH_TRI2D:
	  WlzCMeshUpdateBBox2D(mesh.m2);
	  WlzCMeshUpdateMaxSqEdgLen2D(mesh.m2);
	  errNum = WlzCMeshReassignBuckets2D(mesh.m2, 0);
	  break;
        case WLZ_CMESH_TET3D:
	  WlzCMeshUpdateBBox3D(mesh.m3);
	  WlzCMeshUpdateMaxSqEdgLen3D(mesh.m3);
	  errNum = WlzCMeshReassignBuckets3D(mesh.m3, 0);
	  break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Applies a low pass filter to the geometry of the given
*		mesh. See WlzGMFilterGeomLPLM().
* \param	mesh			Given mesh.
* \param	lambda			Positive filter parameter.
* \param	mu			Negative filter parameter.
* \param	nItr			Number of itterations.
* \param	doBnd			Filter boundary nodes in non-zero.
* \param	update			Update the mesh bucket grid,
*					bounding box and maximum edge length.
*/
WlzErrorNum	WlzCMeshLPFilterLM(WlzCMeshP mesh,
				   double lambda, double mu,
				   int nItr, int doBnd, int update)
{
  int		idI,
  		nVtx = 0;
  WlzVertexP	vtxBuf[2];
  WlzVertexType	vtxType;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_NONE;
  }
  else if((lambda < DBL_EPSILON) || (lambda > (1.0 - DBL_EPSILON)) ||
          (-mu < DBL_EPSILON) || (-mu > (1.0 - DBL_EPSILON)) ||
	  (nItr < 1))
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vtxBuf[0] = WlzDVerticesFromCMesh(mesh, &nVtx, &vtxType, 0, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        if((vtxBuf[1].d2 = (WlzDVertex2 *)
	                   AlcMalloc(sizeof(WlzDVertex2) * nVtx)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idI = 0; idI < nItr; ++idI)
	  {
	    WlzCMeshFilterLPL2D(mesh.m2, vtxBuf[0].d2, vtxBuf[1].d2,
	                        lambda, doBnd);
	    WlzCMeshFilterLPL2D(mesh.m2, vtxBuf[1].d2, vtxBuf[0].d2,
	                        lambda, doBnd);
	  }
	}
	break;
      case WLZ_CMESH_TET3D:
        if((vtxBuf[1].d3 = (WlzDVertex3 *)
	                   AlcMalloc(sizeof(WlzDVertex3) * nVtx)) == NULL)
        {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idI = 0; idI < nItr; ++idI)
	  {
	    WlzCMeshFilterLPL3D(mesh.m3, vtxBuf[0].d3, vtxBuf[1].d3,
	                        lambda, doBnd);
	    WlzCMeshFilterLPL3D(mesh.m3, vtxBuf[1].d3, vtxBuf[0].d3,
	                        lambda, doBnd);
	  }
	}
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshSetVertices(mesh, vtxBuf[0], update);
  }
  return(errNum); 
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief        Sets the position of all valid nodes in the mesh.
*		The mesh will still need to have it's bounding box,
*		maximum edge length and grid box set after this function
*		has repositioned the nodes.
*		Vertex type must be correct for mesh type, ie
*		WlzDVertex2 for WlzCMesh2D and WlzDVertex3 for
*		WlzCMesh3D.
* \param        mesh			The given mesh.
* \param        vtxBuf                  The buffer with vnode positions
*                                       that are to be set in the model.
* \param	update			Update the mesh bucket grid,
*					bounding box and maximum edge length.
*/
WlzErrorNum	WlzCMeshSetVertices(WlzCMeshP mesh, WlzVertexP vtxBuf,
				    int update)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  
  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(vtxBuf.v == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        WlzCMeshSetVertices2D(mesh.m2, vtxBuf.d2, update);
	break;
      case WLZ_CMESH_TET3D:
        WlzCMeshSetVertices3D(mesh.m3, vtxBuf.d3, update);
	break;
      default:
	errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \ingroup	WlzMesh
* \brief        Sets the position of all valid nodes in the mesh.
*		The mesh will still need to have it's bounding box,
*		maximum edge length and grid box set after this function
*		has repositioned the nodes.
*		The parameters are assumed vaild.
* \param        mesh			The given mesh.
* \param        vtxBuf                  The buffer with vnode positions
*                                       that are to be set in the model.
* \param	update			Update the mesh bucket grid,
*					bounding box and maximum edge length.
*/
void		WlzCMeshSetVertices2D(WlzCMesh2D *mesh, WlzDVertex2 *vtxBuf,
				      int update)
{
  int		idx,
  		cnt;
  AlcVector	*vec;
  WlzCMeshNod2D *nod;

  cnt = mesh->res.elm.maxEnt;
  vec = mesh->res.nod.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    nod = (WlzCMeshNod2D *)AlcVectorItemGet(vec, idx);
    if(nod && (nod->idx >= 0))
    {
      nod->pos = *(vtxBuf + idx);
    }
  }
  if(update)
  {
    WlzCMeshUpdateBBox2D(mesh);
    WlzCMeshUpdateMaxSqEdgLen2D(mesh);
    (void )WlzCMeshReassignBuckets2D(mesh, 0);
  }
}

/*!
* \ingroup	WlzMesh
* \brief        Sets the position of all valid nodes in the mesh.
*		The mesh will still need to have it's bounding box,
*		maximum edge length and grid box set after this function
*		has repositioned the nodes.
*		The parameters are assumed vaild.
* \param        mesh			The given mesh.
* \param        vtxBuf                  The buffer with vnode positions
*                                       that are to be set in the model.
* \param	update			Update the mesh bucket grid,
*					bounding box and maximum edge length.
*/
void		WlzCMeshSetVertices3D(WlzCMesh3D *mesh, WlzDVertex3 *vtxBuf,
				      int update)
{
  int		idx,
  		cnt;
  AlcVector	*vec;
  WlzCMeshNod3D *nod;

  cnt = mesh->res.nod.maxEnt;
  vec = mesh->res.nod.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    nod = (WlzCMeshNod3D *)AlcVectorItemGet(vec, idx);
    if(nod && (nod->idx >= 0))
    {
      nod->pos = *(vtxBuf + idx);
    }
  }
  if(update)
  {
    WlzCMeshUpdateBBox3D(mesh);
    WlzCMeshUpdateMaxSqEdgLen3D(mesh);
    (void )WlzCMeshReassignBuckets3D(mesh, 0);
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Checks that the 2D or 3D mesh has valid connectivities.
*		This function is slow and should only be used when
*		debugging mesh connectivities - it is not intended for
*		routine use. With an invalid mesh this checking function
*		may provoke NULL pointer access or segmentation faults.
* \param	mesh			Given mesh.
* \param	dstElm			Destination mesh element pointer
*					for last mesh element, may be NULL.
* \param	allErr			If non zero the checking conmtinues
*					after an error has been found, else if
*					zero the checking stops after the first
*					error has been found.
* \param	fP			Stream for diagnostic output
*					statements - may be NULL in which case
*					there will be no diagnostic output.
*/
WlzErrorNum 	WlzCMeshVerify(WlzCMeshP mesh, void **dstElm,
				 int allErr, FILE *fP)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        errNum = WlzCMeshVerify2D(mesh.m2, (WlzCMeshElm2D **)dstElm,
				  allErr, fP);
        break;
      case WLZ_CMESH_TET3D:
        errNum = WlzCMeshVerify3D(mesh.m3, (WlzCMeshElm3D **)dstElm,
				  allErr, fP);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Checks that the 2D mesh has valid connectivities.
*		This function is slow and should only be used when
*		debugging mesh connectivities - it is not intended for
*		routine use. With an invalid mesh this checking function
*		may provoke NULL pointer access or segmentation faults.
* \param	mesh			Given mesh.
* \param	dstElm			Destination mesh element pointer
*					for last mesh element, may be NULL.
* \param	allErr			If non zero the checking conmtinues
*					after an error has been found, else if
*					zero the checking stops after the first
*					error has been found.
* \param	fP			Stream for diagnostic output
*					statements - may be NULL in which case
*					there will be no diagnostic output.
*/
WlzErrorNum 	WlzCMeshVerify2D(WlzCMesh2D *mesh, WlzCMeshElm2D **dstElm,
				 int allErr, FILE *fP)
{
  int		cnt,
  		idE,
  		idN;
  WlzCMeshEdg2D *edg0,
  		*edg1;
  WlzCMeshElm2D	*elm;
  WlzErrorNum	errNum0,
  		errNum1 = WLZ_ERR_NONE;
  const		nnxtLimit = 1000;
  char		msgBuf[1000];

  if(mesh == NULL)
  {
    errNum1 = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum1 = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    idE = 0;
    while((idE < mesh->res.elm.maxEnt) &&
          ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
    {
      /* Verify elements of mesh. */
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	idN = 0;
	while((idN < 3) &&
	      ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
	{
	  errNum0 = WLZ_ERR_NONE;
	  /* Verify each edge of element. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edg[idN].next != &(elm->edg[(idN + 1) % 3])))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	                   "elm[%d]->edg[%d].next != &(elm[%d]->edg[%d])",
			   idE, idN, idE, (idN + 1) % 3);
	  }
	  /* Verify that each edge is directed from a node. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edg[idN].nod == NULL))
	  {
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edg[%d].nod == NULL",
			   idE, idN);
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	  }
	  /* Verify that each edge's node has not been deleted. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edg[idN].nod->idx < 0))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edg[%d].nod->idx < 0",
			   idE, idN);
	  }
	  /* Verify that the each edge's node is the node. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edg[idN].nod->edg->nod != elm->edg[idN].nod))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
		"elm[%d]->edg[%d].nod->edg->nod != elm[%d]->edg[%d].nod",
		idE, idN, idE, idN);
	  }
	  /* Verify that an opposite opposite edge is the edge. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     ((elm->edg[idN].opp != NULL) &&
	     (elm->edg[idN].opp->opp != &(elm->edg[idN]))))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edg[%d].opp->opp != &(elm[%d]->edg[%d])",
			   idE, idN, idE, idN);
	  }
	  /* Check the number of edges directed from a node is reasonable. */
	  if((allErr == 0)  || (errNum0 == WLZ_ERR_NONE))
	  {
	    cnt = 0;
	    edg1 = edg0 = elm->edg[idN].nod->edg;
	    do
	    {
	      edg1 = edg1->nnxt;
	    }
	    while((cnt++ < nnxtLimit) && (edg1 != edg0));
	    if(cnt >= nnxtLimit)
	    {
	      errNum0 = WLZ_ERR_DOMAIN_DATA;
	      (void )sprintf(msgBuf,
			     "elm[%d]->edg[%d].nod->edg->nnxt cycle > %d",
			     idE, idN, nnxtLimit);
	    }
	  }
	  if(errNum1 == WLZ_ERR_NONE)
	  {
	    errNum1 = errNum0;
	  }
	  ++idN;
	}
	/* Check element areas are positive). */
	if((allErr == 0)  || (errNum1 == WLZ_ERR_NONE))
	{
	  if(WlzCMeshElmSnArea22D(elm) < WLZ_MESH_TOLERANCE_SQ)
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	    		   "WlzCMeshElmSnArea22D(elm[%d]) < %g",
			   idE, WLZ_MESH_TOLERANCE_SQ);
	  }
	}
      }
      ++idE;
    }
  }
  if(dstElm)
  {
    *dstElm = elm;
  }
  return(errNum1);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Checks that the 3D mesh has valid connectivities.
*		This function is slow and should only be used when
*		debugging mesh connectivities - it is not intended for
*		routine use. With an invalid mesh this checking function
*		may provoke NULL pointer access or segmentation faults.
* \param	mesh			Given mesh.
* \param	dstElm			Destination mesh element pointer
*					for last mesh element, may be NULL.
* \param	allErr			If non zero the checking conmtinues
*					after an error has been found, else if
*					zero the checking stops after the first
*					error has been found.
* \param	fP			Stream for diagnostic output
*					statements - may be NULL in which case
*					there will be no diagnostic output.
*/
WlzErrorNum 	WlzCMeshVerify3D(WlzCMesh3D *mesh, WlzCMeshElm3D **dstElm,
				 int allErr, FILE *fP)
{
  int		idE,
  		idF,
		idN,
		cnt;
  WlzCMeshNod3D	*nod;
  WlzCMeshEdg3D	*edg,
  		*edg0,
		*edg1;
  WlzCMeshElm3D	*elm;
  WlzCMeshFace	*fce;
  WlzErrorNum	errNum0,
  		errNum1 = WLZ_ERR_NONE;
  const		nnxtLimit = 1000;
  char		msgBuf[1000];


  if(mesh == NULL)
  {
    errNum1 = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum1 = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    idE = 0;
    while((idE < mesh->res.elm.maxEnt) &&
          ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
    {
      /* Verify elements of mesh. */
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	/* Verify faces of element. */
	idF = 0;
	while((idF < 4) &&
	      ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
	{
	  errNum0 = WLZ_ERR_NONE;
	  fce = elm->face + idF;
	  if(fce->elm != elm)
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	                   "elm[%d]->face[%d].elm != &(elm[%d])",
			   idE, idF, idE);
	  }
	  if(errNum0 == WLZ_ERR_NONE)
	  {
	    for(idN = 0; idN < 3; ++idN)
	    {
	      edg = fce->edg + idN;
	      if(edg->face != fce)
	      {
		errNum0 = WLZ_ERR_DOMAIN_DATA;
	        (void )sprintf(msgBuf,
		  "elm[%d]->face[%d].edg[%d]->face != elm[%d]->face[%d]\n",
		  idE, idF, idN, idE, idF);
	      }
	      if(edg->next->next->next != edg)
	      {
	        errNum0 = WLZ_ERR_DOMAIN_DATA;
		(void )sprintf(msgBuf,
                               "elm[%d]->face[%d].edg[%d].next->next->next |= "
			       "elm[%d]->face[%d].edg[%d]\n",
			       idE, idF, idN, idE, idF, idN);
	      }
	      nod = edg->nod;
	      if(nod == NULL)
	      {
	        errNum0 = WLZ_ERR_DOMAIN_DATA;
		(void )sprintf(msgBuf,
                               "elm[%d]->face[%d].edg[%d].nod == NULL\n",
			       idE, idF, idN);
	      }
	      else
	      {
	        if(nod->idx < 0)
		{
		  errNum0 = WLZ_ERR_DOMAIN_DATA;
		  (void )sprintf(msgBuf,
		                 "elm[%d]->face[%d].edg[%d].nod->idx < 0\n",
				 idE, idF, idN);

		}
		else if(nod->edg == NULL)
		{
		  errNum0 = WLZ_ERR_DOMAIN_DATA;
		  (void )sprintf(msgBuf,
			       "elm[%d]->face[%d].edg[%d].nod->edg == NULL\n",
			       idE, idF, idN);
		}
		if((allErr == 0)  || (errNum0 == WLZ_ERR_NONE))
		{
		  cnt = 0;
		  edg1 = edg0 = nod->edg;
		  do
		  {
		    edg1 = edg1->nnxt;
		  }
		  while((cnt++ < nnxtLimit) && (edg1 != edg0));
		  if(cnt >= nnxtLimit)
		  {
		    errNum0 = WLZ_ERR_DOMAIN_DATA;
		    (void )sprintf(msgBuf,
			 "elm[%d]->face[%d].edg[%d].nod->edg->nnxt cycle > %d",
				   idE, idF, idN, nnxtLimit);
		  }
		}
	      }
	    }
	  }
	  ++idF;
	}
	if(errNum1 == WLZ_ERR_NONE)
	{
	  errNum1 = errNum0;
	}
	/* Check element volumes are positive). */
	if((allErr == 0)  || (errNum1 == WLZ_ERR_NONE))
	{
	  if(WlzCMeshElmSnVolume63D(elm) < WLZ_MESH_TOLERANCE_SQ)
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
			   "WlzCMeshElmSnVolume6(elm[%d]) < %g",
			   idE, WLZ_MESH_TOLERANCE_SQ);
	  }
	}
      }
      ++idE;
    }
  }
  if(dstElm)
  {
    *dstElm = elm;
  }
  return(errNum1);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes some simple geometric features of all valid
*		elements in a mesh.
* \param	mesh			Given mesh.
* \param	dstNElm			Destination pointer for the number
*					of mesh elements. May be NULL.
* \param	dstIdx			Destination pointer for element
* 					indices. May be NULL.
* \param	dstVol			Destination pointer for the area
*					or volume of the elements. May be
*					NULL.
* \param	dstMinLen		Destination pointer for the minimum
*					edge length of the elements.
* \param	dstMaxLen		Destination pointer for the minimum
*					edge length of the elements. May be
*					NULL.
*/
WlzErrorNum 	WlzCMeshCmpElmFeat(WlzCMeshP mesh, int *dstNElm,
				   int **dstIdx, double **dstVol,
				   double **dstMinLen, double **dstMaxLen)
{
  int		nElm = 0;
  int		*idx = NULL;
  double	*vol = NULL,
  		*minLen = NULL,
		*maxLen = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        errNum = WlzCMeshCmpElmFeat2D(mesh.m2, &nElm, &idx, &vol,
	                              &minLen, &maxLen);
        break;
      case WLZ_CMESH_TET3D:
        errNum = WlzCMeshCmpElmFeat3D(mesh.m3, &nElm, &idx, &vol,
	                              &minLen, &maxLen);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNElm)
    {
      *dstNElm = nElm;
    }
    if(dstIdx)
    {
      *dstIdx = idx;
    }
    if(dstVol)
    {
      *dstVol = vol;
    }
    if(dstMinLen)
    {
      *dstMinLen = minLen;
    }
    if(dstMaxLen)
    {
      *dstMaxLen = maxLen;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes some simple geometric features of all valid
*		elements in a 2D mesh.
* \param	mesh			Given 2D mesh.
* \param	dstNElm			Destination pointer for the number
*					of mesh elements. May be NULL.
* \param	dstIdx			Destination pointer for element
* 					indices. May be NULL.
* \param	dstVol			Destination pointer for the area
*					or volume of the elements. May be
*					NULL.
* \param	dstMinLen		Destination pointer for the minimum
*					edge length of the elements.
* \param	dstMaxLen		Destination pointer for the minimum
*					edge length of the elements. May be
*					NULL.
*/
WlzErrorNum 	WlzCMeshCmpElmFeat2D(WlzCMesh2D *mesh, int *dstNElm,
				   int **dstIdx, double **dstVol,
				   double **dstMinLen, double **dstMaxLen)
{
  int		idE,
		idV,
  		nElm = 0;
  int		*idx = NULL;
  double	*vol = NULL,
  		*minLen = NULL,
		*maxLen = NULL;
  double	len[3];
  WlzDVertex2	tV0;
  WlzCMeshElm2D	*elm;
  WlzCMeshNod2D	*nod[3];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(mesh->res.elm.maxEnt > 0)
  {
    if(dstIdx &&
       ((idx = AlcCalloc(mesh->res.elm.maxEnt, sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstVol &&
       ((vol = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstMinLen &&
       ((minLen = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstMaxLen &&
       ((maxLen = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idV = 0;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	nod[0] = elm->edg[0].nod;
	nod[1] = elm->edg[1].nod;
	nod[2] = elm->edg[2].nod;
	if(idx)
	{
	  *(idx + idV) = elm->idx;
	}
	if(vol)
	{
	  *(vol + idV) = WlzGeomTriangleSnArea2(nod[0]->pos,
                                nod[1]->pos, nod[2]->pos) / 2.0;
	}
	if(minLen || maxLen)
	{
	  WLZ_VTX_2_SUB(tV0, nod[0]->pos, nod[1]->pos);
	  len[0] = WLZ_VTX_2_SQRLEN(tV0);
	  WLZ_VTX_2_SUB(tV0, nod[0]->pos, nod[2]->pos);
	  len[1] = WLZ_VTX_2_SQRLEN(tV0);
	  WLZ_VTX_2_SUB(tV0, nod[1]->pos, nod[2]->pos);
	  len[2] = WLZ_VTX_2_SQRLEN(tV0);
	}
	if(minLen)
	{
	  AlgRankSelectD(len, 3, 0);
	  *(minLen + idV) = (len[0] < DBL_EPSILON)? 0.0: sqrt(len[0]);
	}
	if(maxLen)
	{
	  AlgRankSelectD(len, 3, 2);
	  *(maxLen + idV) = (len[2] < DBL_EPSILON)? 0.0: sqrt(len[2]);
	}
        ++idV;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNElm)
    {
      *dstNElm = idV;
    }
    if(dstIdx)
    {
      *dstIdx = idx;
    }
    if(dstVol)
    {
      *dstVol = vol;
    }
    if(dstMinLen)
    {
      *dstMinLen = minLen;
    }
    if(dstMaxLen)
    {
      *dstMaxLen = maxLen;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzMesh
* \brief	Computes some simple geometric features of all valid
*		elements in a 3D mesh.
* \param	mesh			Given 3D mesh.
* \param	dstNElm			Destination pointer for the number
*					of mesh elements. May be NULL.
* \param	dstIdx			Destination pointer for element
* 					indices. May be NULL.
* \param	dstVol			Destination pointer for the area
*					or volume of the elements. May be
*					NULL.
* \param	dstMinLen		Destination pointer for the minimum
*					edge length of the elements.
* \param	dstMaxLen		Destination pointer for the minimum
*					edge length of the elements. May be
*					NULL.
*/
WlzErrorNum 	WlzCMeshCmpElmFeat3D(WlzCMesh3D *mesh, int *dstNElm,
				   int **dstIdx, double **dstVol,
				   double **dstMinLen, double **dstMaxLen)
{
  int		idE,
		idV,
  		nElm = 0;
  int		*idx = NULL;
  double	*vol = NULL,
  		*minLen = NULL,
		*maxLen = NULL;
  double	len[6];
  WlzDVertex3	tV0;
  WlzCMeshElm3D	*elm;
  WlzCMeshNod3D	*nod[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(mesh->res.elm.maxEnt > 0)
  {
    if(dstIdx &&
       ((idx = AlcCalloc(mesh->res.elm.maxEnt, sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstVol &&
       ((vol = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstMinLen &&
       ((minLen = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    if(dstMaxLen &&
       ((maxLen = AlcCalloc(mesh->res.elm.maxEnt, sizeof(double))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idV = 0;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	WlzCMeshElmGetNodes3D(elm, nod + 0, nod + 1, nod + 2, nod + 3);
	if(idx)
	{
	  *(idx + idV) = elm->idx;
	}
	if(vol)
	{
	  *(vol + idV) = WlzGeomTetraSnVolume6(nod[0]->pos,
                                nod[1]->pos, nod[2]->pos, nod[3]->pos) / 6.0;
	}
	if(minLen || maxLen)
	{
	  WLZ_VTX_3_SUB(tV0, nod[0]->pos, nod[1]->pos);
	  len[0] = WLZ_VTX_3_SQRLEN(tV0);
	  WLZ_VTX_3_SUB(tV0, nod[0]->pos, nod[2]->pos);
	  len[1] = WLZ_VTX_3_SQRLEN(tV0);
	  WLZ_VTX_3_SUB(tV0, nod[0]->pos, nod[3]->pos);
	  len[2] = WLZ_VTX_3_SQRLEN(tV0);
	  WLZ_VTX_3_SUB(tV0, nod[1]->pos, nod[2]->pos);
	  len[3] = WLZ_VTX_3_SQRLEN(tV0);
	  WLZ_VTX_3_SUB(tV0, nod[1]->pos, nod[3]->pos);
	  len[4] = WLZ_VTX_3_SQRLEN(tV0);
	  WLZ_VTX_3_SUB(tV0, nod[2]->pos, nod[3]->pos);
	  len[5] = WLZ_VTX_3_SQRLEN(tV0);
	}
	if(minLen)
	{
	  AlgRankSelectD(len, 6, 0);
	  *(minLen + idV) = (len[0] < DBL_EPSILON)? 0.0: sqrt(len[0]);
	}
	if(maxLen)
	{
	  AlgRankSelectD(len, 6, 5);
	  *(maxLen + idV) = (len[5] < DBL_EPSILON)? 0.0: sqrt(len[5]);
	}
        ++idV;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(dstNElm)
    {
      *dstNElm = idV;
    }
    if(dstIdx)
    {
      *dstIdx = idx;
    }
    if(dstVol)
    {
      *dstVol = vol;
    }
    if(dstMinLen)
    {
      *dstMinLen = minLen;
    }
    if(dstMaxLen)
    {
      *dstMaxLen = maxLen;
    }
  }
  return(errNum);
}

/*!
* \return       Twice the signed area of the 2D mesh element.
* \ingroup      WlzMesh
* \brief        Computes twice the signed area of the 2D mesh element.
* \param        elm                     Given mesh element.
*/
double          WlzCMeshElmSnArea22D(WlzCMeshElm2D *elm)
{
  double        area;

  area = WlzGeomTriangleSnArea2(elm->edg[0].nod->pos,
                                elm->edg[1].nod->pos,
                                elm->edg[2].nod->pos);
  return(area);
}

/*!
* \return       Siz times the signed volume of the 3D mesh element.
* \ingroup      WlzMesh
* \brief        Computes six times the signed volume of the 3D mesh
*		element.
* \param        elm                     Given mesh element.
*/
double          WlzCMeshElmSnVolume63D(WlzCMeshElm3D *elm)
{
  double        vol;
  WlzCMeshNod3D	*nod[4];

  WlzCMeshElmGetNodes3D(elm, nod + 0, nod + 1, nod + 2, nod + 3);
  vol = WlzGeomTetraSnVolume6(nod[0]->pos,
                              nod[1]->pos,
                              nod[2]->pos,
                              nod[3]->pos);
  return(vol);
}

/*!
* \ingroup	WlzMesh
* \brief	Gets the four nodes of a 3D element.
* \param	elm			Given mesh element.
* \param	dstNod0			First destination pointer for node.
* \param	dstNod1			Second destination pointer for node.
* \param	dstNod2			Third destination pointer for node.
* \param	dstNod3			Forth destination pointer for node.
*/
void		WlzCMeshElmGetNodes3D(WlzCMeshElm3D *elm,
				      WlzCMeshNod3D **dstNod0,
				      WlzCMeshNod3D **dstNod1,
				      WlzCMeshNod3D **dstNod2,
				      WlzCMeshNod3D **dstNod3)
{
  int		idN;
  WlzCMeshNod3D	*nod;

  *dstNod0 = elm->face[0].edg[0].nod;
  *dstNod1 = elm->face[0].edg[1].nod;
  *dstNod2 = elm->face[0].edg[2].nod;
  for(idN = 0; idN < 3; ++idN)
  {
    nod = elm->face[1].edg[idN].nod;
    if((nod != *dstNod0) &&
       (nod != *dstNod1) &&
       (nod != *dstNod2))
    {
      *dstNod3  = nod;
      break;
    }
  }
}

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Filters the geometry of the verticies in a 2D mesh using
*               the given input and output buffers for the mesh node
*		positions.
* \note         See WlzGMFilterGeomLPLM().
* \param        model                   The given model.
* \param        vGIn                    Input vertex geometries.
* \param        vGOut                   Output vertex geometries.
* \param        lambda                  The filter parameter.
* \param        nonMan                  If non-zero allows non manifold
*/
static void	WlzCMeshFilterLPL2D(WlzCMesh2D *mesh,
				    WlzDVertex2 *vGIn, WlzDVertex2 *vGOut,
				    double lambda, int doBnd)
{
  int           idx,
                cnt;
  WlzDVertex2   tV0,
                tV1;
  WlzCMeshNod2D  *cN;
  AlcVector     *vec;

  cnt = mesh->res.nod.maxEnt;
  vec = mesh->res.nod.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    cN = (WlzCMeshNod2D *)AlcVectorItemGet(vec, idx);
    if(cN->idx >= 0)
    {
      tV0 = *(vGIn + idx);
      if(doBnd || ((cN->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) == 0))
      {
	tV1 = WlzCMeshFilterLPLDelta2D(mesh, cN, vGIn, doBnd);
	tV0.vtX += lambda * tV1.vtX;
	tV0.vtY += lambda * tV1.vtY;
      }
      *(vGOut + idx) = tV0;
    }
  }
}

/*!
* \return       void
* \ingroup      WlzMesh
* \brief        Filters the geometry of the verticies in a 3D mesh using
*               the given input and output buffers for the mesh node
*		positions.
* \note         See WlzGMFilterGeomLPLM().
* \param        model                   The given model.
* \param        vGIn                    Input vertex geometries.
* \param        vGOut                   Output vertex geometries.
* \param        lambda                  The filter parameter.
* \param        nonMan                  If non-zero allows non manifold
*/
static void	WlzCMeshFilterLPL3D(WlzCMesh3D *mesh,
				    WlzDVertex3 *vGIn, WlzDVertex3 *vGOut,
				    double lambda, int doBnd)
{
  int           idx,
                cnt;
  WlzDVertex3   tV0,
                tV1;
  WlzCMeshNod3D  *cN;
  AlcVector     *vec;

  cnt = mesh->res.nod.maxEnt;
  vec = mesh->res.nod.vec;
  for(idx = 0; idx < cnt; ++idx)
  {
    cN = (WlzCMeshNod3D *)AlcVectorItemGet(vec, idx);
    if(cN->idx >= 0)
    {
      tV0 = *(vGIn + idx);
      if(doBnd || ((cN->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) == 0))
      {
	tV1 = WlzCMeshFilterLPLDelta3D(mesh, cN, vGIn, doBnd);
	tV0.vtX += lambda * tV1.vtX;
	tV0.vtY += lambda * tV1.vtY;
	tV0.vtZ += lambda * tV1.vtZ;
      }
      *(vGOut + idx) = tV0;
    }
  }
}

/*!
* \return	Vertex displacement.
* \ingroup	WlzMesh
* \brief	Computes the displacement of the given node, for use
*		by WlzCMeshFilterLPL2D(). This is just the displacement
*		from the nodes position to the mean of the directly
*		connected neighbour node positions. All positions are
*		taken from the vertex buffer using the nodes index value.
* \param	mesh			The mesh.
* \param	nod			The current node.
* \param	vBuf			Buffer of node positions.
* \param	doBnd			Filter boundary nodes if non-zero.
*/
static WlzDVertex2 WlzCMeshFilterLPLDelta2D(WlzCMesh2D *mesh,
					WlzCMeshNod2D *nod, WlzDVertex2 *vBuf,
					int doBnd)
{
  int           nN = 0;
  double        tD0;
  WlzCMeshNod2D	*oNod;
  WlzCMeshEdg2D	*edg0,
  		*edg1;
  WlzDVertex2   nP,
                sP;

  sP.vtX = sP.vtY = 0.0;
  if((nod->idx >= 0) &&
     (doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0)))
  {
    edg0 = edg1 = nod->edg;
    do
    {
      oNod = edg1->next->nod;
      nP = *(vBuf + oNod->idx);
      WLZ_VTX_2_ADD(sP, sP, nP);
      ++nN;
      edg1 = edg1->nnxt;
    } while(edg0 != edg1);
    if(nN > 0)
    {
      tD0 = 1.0 / nN;
      nP = *(vBuf + nod->idx);
      sP.vtX = (tD0 * sP.vtX) - nP.vtX;
      sP.vtY = (tD0 * sP.vtY) - nP.vtY;
    }
  }
  return(sP);
}

/*!
* \return	Vertex displacement.
* \ingroup	WlzMesh
* \brief	Computes the displacement of the given node, for use
*		by WlzCMeshFilterLPL3D(). This is just the displacement
*		from the nodes position to the mean of the directly
*		connected neighbour node positions. All positions are
*		taken from the vertex buffer using the nodes index value.
* \param	mesh			The mesh.
* \param	nod			The current node.
* \param	vBuf			Buffer of node positions.
* \param	doBnd			Filter boundary nodes if non-zero.
*/
static WlzDVertex3 WlzCMeshFilterLPLDelta3D(WlzCMesh3D *mesh,
					WlzCMeshNod3D *nod, WlzDVertex3 *vBuf,
					int doBnd)
{
  int           nN = 0;
  double        tD0;
  WlzCMeshNod3D	*oNod;
  WlzCMeshEdg3D	*edg0,
  		*edg1;
  WlzDVertex3   nP,
                sP;

  sP.vtX = sP.vtY = sP.vtZ = 0.0;
  if((nod->idx >= 0) &&
     (doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0)))
  {
    edg0 = edg1 = nod->edg;
    do
    {
      oNod = edg1->next->nod;
      nP = *(vBuf + oNod->idx);
      WLZ_VTX_3_ADD(sP, sP, nP);
      ++nN;
      edg1 = edg1->nnxt;
    } while(edg0 != edg1);
    if(nN > 0)
    {
      tD0 = 1.0 / nN;
      nP = *(vBuf + nod->idx);
      sP.vtX = (tD0 * sP.vtX) - nP.vtX;
      sP.vtY = (tD0 * sP.vtY) - nP.vtY;
      sP.vtZ = (tD0 * sP.vtZ) - nP.vtZ;
    }
  }
  return(sP);
}
