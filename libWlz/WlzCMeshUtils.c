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
* \brief	Utility functions for 2D and 3D graph based conforming
* 		simplical meshes.
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
	dSq = WlzGeomDistSq2D(elm->edu[0].nod->pos,
	                      elm->edu[1].nod->pos);
	if(dSq > mesh->maxSqEdgLen)
	{
	  mesh->maxSqEdgLen = dSq;
	}
	dSq = WlzGeomDistSq2D(elm->edu[1].nod->pos,
			      elm->edu[2].nod->pos);
	if(dSq > mesh->maxSqEdgLen)
	{
	  mesh->maxSqEdgLen = dSq;
	}
	dSq = WlzGeomDistSq2D(elm->edu[2].nod->pos,
	                      elm->edu[0].nod->pos);
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
	  dSq = WlzGeomDistSq3D(elm->face[idF].edu[0].nod->pos,
				elm->face[idF].edu[1].nod->pos);
	  if(dSq > mesh->maxSqEdgLen)
	  {
	    mesh->maxSqEdgLen = dSq;
	  }
	  dSq = WlzGeomDistSq3D(elm->face[idF].edu[1].nod->pos,
				elm->face[idF].edu[2].nod->pos);
	  if(dSq > mesh->maxSqEdgLen)
	  {
	    mesh->maxSqEdgLen = dSq;
	  }
	  dSq = WlzGeomDistSq3D(elm->face[idF].edu[2].nod->pos,
				elm->face[idF].edu[0].nod->pos);
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
* \brief	Sets the node flags for all valid nodes of the mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to set.
*/
void		WlzCMeshSetNodFlags(WlzCMeshP mesh, unsigned int flags)
{
  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        WlzCMeshSetNodFlags2D(mesh.m2, flags);
	break;
      case WLZ_CMESH_TET3D:
        WlzCMeshSetNodFlags3D(mesh.m3, flags);
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Sets node flags for all valid nodes of the 2D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to set.
*/
void		WlzCMeshSetNodFlags2D(WlzCMesh2D *mesh, unsigned int flags)
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
	nod->flags |= flags;
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Sets node flags for all valid nodes of the 3D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to set.
*/
void		WlzCMeshSetNodFlags3D(WlzCMesh3D *mesh, unsigned int flags)
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
	nod->flags |= flags;
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears the node flags for all valid nodes of the mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearNodFlags(WlzCMeshP mesh, unsigned int flags)
{
  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        WlzCMeshClearNodFlags2D(mesh.m2, flags);
	break;
      case WLZ_CMESH_TET3D:
        WlzCMeshClearNodFlags3D(mesh.m3, flags);
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears node flags for all valid nodes of the 2D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearNodFlags2D(WlzCMesh2D *mesh, unsigned int flags)
{
  int		idN;
  WlzCMeshNod2D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    flags = ~flags;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	nod->flags &= flags;
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears node flags for all valid nodes of the 3D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearNodFlags3D(WlzCMesh3D *mesh, unsigned int flags)
{
  int		idN;
  WlzCMeshNod3D *nod;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    flags = ~flags;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	nod->flags &= flags;
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears the element flags for all valid elements of the mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearElmFlags(WlzCMeshP mesh, unsigned int flags)
{
  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        WlzCMeshClearElmFlags2D(mesh.m2, flags);
	break;
      case WLZ_CMESH_TET3D:
        WlzCMeshClearElmFlags3D(mesh.m3, flags);
	break;
      default:
        break;
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears element flags for all valid elements of the 2D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearElmFlags2D(WlzCMesh2D *mesh, unsigned int flags)
{
  int		idN;
  WlzCMeshElm2D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    flags = ~flags;
    for(idN = 0; idN < mesh->res.elm.maxEnt; ++idN)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idN);
      if(elm->idx >= 0)
      {
	elm->flags &= flags;
      }
    }
  }
}

/*!
* \ingroup	WlzMesh
* \brief	Clears element flags for all valid elements of the 3D mesh.
* \param	mesh			Given mesh.
* \param	flags			Flags to clear.
*/
void		WlzCMeshClearElmFlags3D(WlzCMesh3D *mesh, unsigned int flags)
{
  int		idN;
  WlzCMeshElm3D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    flags = ~flags;
    for(idN = 0; idN < mesh->res.elm.maxEnt; ++idN)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idN);
      if(elm->idx >= 0)
      {
	elm->flags &= flags;
      }
    }
  }
}

/*!
* \return	Number of boundary nodes.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundNodFlags(WlzCMeshP mesh)
{
  int		nBnd = 0;
  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        nBnd = WlzCMeshSetBoundNodFlags2D(mesh.m2);
	break;
      case WLZ_CMESH_TET3D:
        nBnd = WlzCMeshSetBoundNodFlags3D(mesh.m3);
	break;
      default:
        break;
    }
  }
  return(nBnd);
}

/*!
* \return	Number of boundary nodes.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the 2D mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundNodFlags2D(WlzCMesh2D *mesh)
{
  int		idN,
  		nBnd = 0;
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
	  ++nBnd;
	  nod->flags |= WLZ_CMESH_NOD_FLAG_BOUNDARY;
	}
      }
    }
  }
  return(nBnd);
}

/*!
* \return	Number of boundary nodes.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary node flag bit for all nodes
*		of the 3D mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundNodFlags3D(WlzCMesh3D *mesh)
{
  int		idN,
  		nBnd = 0;
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
	  ++nBnd;
	  nod->flags |= WLZ_CMESH_NOD_FLAG_BOUNDARY;
	}
      }
    }
  }
  return(nBnd);
}

/*!
* \return	Number of boundary elements.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary element flag bit for all elements
*		of the mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundElmFlags(WlzCMeshP mesh)
{
  int		nBnd = 0;
  if(mesh.v)
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        nBnd = WlzCMeshSetBoundElmFlags2D(mesh.m2);
	break;
      case WLZ_CMESH_TET3D:
        nBnd = WlzCMeshSetBoundElmFlags3D(mesh.m3);
	break;
      default:
        break;
    }
  }
  return(nBnd);
}

/*!
* \return	Number of boundary elements.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary element flag bit for all elements
*		of the 2D mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundElmFlags2D(WlzCMesh2D *mesh)
{
  int		idN,
  		nBnd = 0;
  WlzCMeshElm2D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TRI2D))
  {
    for(idN = 0; idN < mesh->res.elm.maxEnt; ++idN)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idN);
      if(elm->idx >= 0)
      {
	elm->flags &= ~(WLZ_CMESH_ELM_FLAG_BOUNDARY);
	if(WlzCMeshElmIsBoundary2D(elm))
	{
	  ++nBnd;
	  elm->flags |= WLZ_CMESH_ELM_FLAG_BOUNDARY;
	}
      }
    }
  }
  return(nBnd);
}

/*!
* \return	Number of boundary elements.
* \ingroup	WlzMesh
* \brief	Sets or clears the boundary element flag bit for all elements
*		of the 3D mesh.
* \param	mesh			Given mesh.
*/
int		WlzCMeshSetBoundElmFlags3D(WlzCMesh3D *mesh)
{
  int		idN,
  		nBnd = 0;
  WlzCMeshElm3D *elm;

  if(mesh && (mesh->type == WLZ_CMESH_TET3D))
  {
    for(idN = 0; idN < mesh->res.elm.maxEnt; ++idN)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idN);
      if(elm->idx >= 0)
      {
	elm->flags &= ~(WLZ_CMESH_ELM_FLAG_BOUNDARY);
	if(WlzCMeshElmIsBoundary3D(elm))
	{
	  ++nBnd;
	  elm->flags |= WLZ_CMESH_ELM_FLAG_BOUNDARY;
	}
      }
    }
  }
  return(nBnd);
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
  WlzCMeshEdgU2D *edu0,
  		*edu1;

  if(nod && (nod->idx >= 0))
  {
    edu0 = edu1 = nod->edu;
    do
    {
      if((edu1->opp == NULL) || (edu1->opp == edu1))
      {
	isBnd = 1;
	break;
      }
      edu1 = edu1->nnxt;
    } while(edu1 != edu0);
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
  WlzCMeshEdgU3D *edu0,
  		*edu1;

  if(nod && (nod->idx >= 0))
  {
    edu0 = edu1 = nod->edu;
    do
    {
      fce = edu1->face;
      if((fce->opp == NULL) || (fce->opp == fce))
      {
	isBnd = 1;
	break;
      }
      edu1 = edu1->nnxt;
    } while(edu1 != edu0);
  }
  return(isBnd);
}

/*!
* \return				Non-zero if the element is a boundary
* 					element.
* \ingroup	WlzMesh
* \brief	Checks whether the element is a boundary node by examining
* 		the edges for opposite edges.
* \param	elm			Given element of mesh.
*/
int		WlzCMeshElmIsBoundary2D(WlzCMeshElm2D *elm)
{
  int		isBnd = 0;

  isBnd = (elm->edu[0].opp == NULL) || (elm->edu[0].opp == &(elm->edu[0])) ||
          (elm->edu[1].opp == NULL) || (elm->edu[1].opp == &(elm->edu[1])) ||
          (elm->edu[2].opp == NULL) || (elm->edu[2].opp == &(elm->edu[2]));
  return(isBnd);
}

/*!
* \return				Non-zero if the element is a boundary
* 					element.
* \ingroup	WlzMesh
* \brief	Checks whether the element is a boundary node by examining
* 		the faces for opposite faces.
* \param	elm			Given element of mesh.
*/
int		WlzCMeshElmIsBoundary3D(WlzCMeshElm3D *elm)
{
  int		isBnd = 0;

  isBnd = (elm->face[0].opp == NULL) || (elm->face[0].opp == &(elm->face[0])) ||
          (elm->face[1].opp == NULL) || (elm->face[1].opp == &(elm->face[1])) ||
          (elm->face[2].opp == NULL) || (elm->face[2].opp == &(elm->face[2])) ||
          (elm->face[3].opp == NULL) || (elm->face[3].opp == &(elm->face[3]));
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
  WlzCMeshEdgU2D *edu0,
  		*edu1;
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
	    edu0 = edu1 = nod->edu;
	    do
	    {
	      oNod = edu1->next->nod;
	      nPos.vtX += oNod->pos.vtX;
	      nPos.vtY += oNod->pos.vtY;
	      ++nCnt;
	      edu1 = edu1->nnxt;
	    } while(edu0 != edu1);
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
  WlzCMeshEdgU3D *edu0,
  		*edu1;
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
	    edu0 = edu1 = nod->edu;
	    do
	    {
	      oNod = edu1->next->nod;
	      nPos.vtX += oNod->pos.vtX;
	      nPos.vtY += oNod->pos.vtY;
	      nPos.vtZ += oNod->pos.vtZ;
	      ++nCnt;
	      edu1 = edu1->nnxt;
	    } while(edu0 != edu1);
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
  int		cnt0,
		cnt1,
  		idE0,
		idE1,
  		idN;
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  WlzCMeshElm2D	*elm;
  WlzCMeshNod2D	*nod;
  WlzErrorNum	errNum0,
  		errNum1 = WLZ_ERR_NONE;
  const int	nnxtLimit = 1000;
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
    idE0 = 0;
    while((idE0 < mesh->res.elm.maxEnt) &&
          ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
    {
      /* Verify elements of mesh. */
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE0);
      if(elm->idx >= 0)
      {
	idN = 0;
	while((idN < 3) &&
	      ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
	{
	  errNum0 = WLZ_ERR_NONE;
	  /* Verify each edge of element. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edu[idN].next != &(elm->edu[(idN + 1) % 3])))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	                   "elm[%d]->edu[%d].next != &(elm[%d]->edu[%d])",
			   idE0, idN, idE0, (idN + 1) % 3);
	  }
	  /* Verify that each edge is directed from a node. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edu[idN].nod == NULL))
	  {
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edu[%d].nod == NULL",
			   idE0, idN);
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	  }
	  /* Verify that each edge's node has not been deleted. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edu[idN].nod->idx < 0))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edu[%d].nod->idx < 0",
			   idE0, idN);
	  }
	  /* Verify that the each edge's node is the node. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     (elm->edu[idN].nod->edu->nod != elm->edu[idN].nod))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
		"elm[%d]->edu[%d].nod->edu->nod != elm[%d]->edu[%d].nod",
		idE0, idN, idE0, idN);
	  }
	  /* Verify that an opposite opposite edge is the edge. */
	  if(((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)) &&
	     ((elm->edu[idN].opp != NULL) &&
	     (elm->edu[idN].opp->opp != &(elm->edu[idN]))))
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	    		   "elm[%d]->edu[%d].opp->opp != &(elm[%d]->edu[%d])",
			   idE0, idN, idE0, idN);
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
			   idE0, WLZ_MESH_TOLERANCE_SQ);
	  }
	}
      }
      ++idE0;
    }
  }
  if(errNum1 == WLZ_ERR_NONE)
  {
    idN = 0;
    while((idN < mesh->res.nod.maxEnt) &&
          ((allErr == 0)  || (errNum1 == WLZ_ERR_NONE)))
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	/* Check that the number of edge uses of each node is reasonable. */
	cnt0 = 0;
        edu0 = edu1 = nod->edu;
	do
	{
	  edu1 = edu1->nnxt;
	}
	while((cnt0++ < nnxtLimit) && (edu1 != edu0));
	if(cnt0 >= nnxtLimit)
	{
	  errNum0 = WLZ_ERR_DOMAIN_DATA;
	  (void )sprintf(msgBuf,
			 "elm[%d]->edu[%d].nod->edu->nnxt cycle > %d",
			 idE0, idN, nnxtLimit);
	}
	/* Check that the number of edge uses of each node is correct. */
	if((allErr == 0)  || (errNum0 == WLZ_ERR_NONE))
	{
	  cnt1 = 0;
	  idE0 = 0;
	  while((idE0 < mesh->res.elm.maxEnt) &&
		((allErr == 0)  || (errNum0 == WLZ_ERR_NONE)))
	  {
	    elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE0);
	    if(elm->idx >= 0)
	    {
	      for(idE1 = 0; idE1 < 3; ++idE1)
	      {
	        if(elm->edu[idE1].nod == nod)
		{
		  ++cnt1;
		}
	      }
	    }
	    ++idE0;
	  }
	  if(cnt0 != cnt1)
	  {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	                   "node %d edu->nnxt loop is inconsistant",
			   idN);
	  }
	}
	if(errNum1 == WLZ_ERR_NONE)
	{
	  errNum1 = errNum0;
	}
      }
      ++idN;
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
  WlzCMeshEdgU3D *edu,
  		*edu0,
		*edu1;
  WlzCMeshElm3D	*elm;
  WlzCMeshFace	*fce;
  WlzErrorNum	errNum0,
  		errNum1 = WLZ_ERR_NONE;
  const int	nnxtLimit = 1000;
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
	      edu = fce->edu + idN;
	      if(edu->face != fce)
	      {
		errNum0 = WLZ_ERR_DOMAIN_DATA;
	        (void )sprintf(msgBuf,
		  "elm[%d]->face[%d].edu[%d]->face != elm[%d]->face[%d]\n",
		  idE, idF, idN, idE, idF);
	      }
	      if(edu->next->next->next != edu)
	      {
	        errNum0 = WLZ_ERR_DOMAIN_DATA;
		(void )sprintf(msgBuf,
                               "elm[%d]->face[%d].edu[%d].next->next->next |= "
			       "elm[%d]->face[%d].edu[%d]\n",
			       idE, idF, idN, idE, idF, idN);
	      }
	      nod = edu->nod;
	      if(nod == NULL)
	      {
	        errNum0 = WLZ_ERR_DOMAIN_DATA;
		(void )sprintf(msgBuf,
                               "elm[%d]->face[%d].edu[%d].nod == NULL\n",
			       idE, idF, idN);
	      }
	      else
	      {
	        if(nod->idx < 0)
		{
		  errNum0 = WLZ_ERR_DOMAIN_DATA;
		  (void )sprintf(msgBuf,
		                 "elm[%d]->face[%d].edu[%d].nod->idx < 0\n",
				 idE, idF, idN);

		}
		else if(nod->edu == NULL)
		{
		  errNum0 = WLZ_ERR_DOMAIN_DATA;
		  (void )sprintf(msgBuf,
			       "elm[%d]->face[%d].edu[%d].nod->edu == NULL\n",
			       idE, idF, idN);
		}
		if((allErr == 0)  || (errNum0 == WLZ_ERR_NONE))
		{
		  cnt = 0;
		  edu1 = edu0 = nod->edu;
		  do
		  {
		    edu1 = edu1->nnxt;
		  }
		  while((cnt++ < nnxtLimit) && (edu1 != edu0));
		  if(cnt >= nnxtLimit)
		  {
		    errNum0 = WLZ_ERR_DOMAIN_DATA;
		    (void )sprintf(msgBuf,
			 "elm[%d]->face[%d].edu[%d].nod->edu->nnxt cycle > %d",
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
	if(errNum1 == WLZ_ERR_NONE)
	{
	  errNum1 = errNum0;
	}
	/* Check the node ordering is consistent. */
	if((allErr == 0)  || (errNum1 == WLZ_ERR_NONE))
	{
	  if((elm->face[1].edu[0].nod != elm->face[0].edu[0].nod) ||
	     (elm->face[1].edu[2].nod != elm->face[0].edu[1].nod) ||
	     (elm->face[2].edu[0].nod != elm->face[0].edu[0].nod) ||
	     (elm->face[2].edu[1].nod != elm->face[0].edu[2].nod) ||
	     (elm->face[2].edu[2].nod != elm->face[1].edu[1].nod) ||
	     (elm->face[3].edu[0].nod != elm->face[0].edu[2].nod) ||
	     (elm->face[3].edu[1].nod != elm->face[0].edu[1].nod) ||
	     (elm->face[3].edu[2].nod != elm->face[1].edu[1].nod))
          {
	    errNum0 = WLZ_ERR_DOMAIN_DATA;
	    (void )sprintf(msgBuf,
	                   "Node ordering for elm[%d] is inconsistent",
			   idE);
	  }
	}
	if(errNum1 == WLZ_ERR_NONE)
	{
	  errNum1 = errNum0;
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
		idV;
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
	nod[0] = elm->edu[0].nod;
	nod[1] = elm->edu[1].nod;
	nod[2] = elm->edu[2].nod;
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
		idV;
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

  area = WlzGeomTriangleSnArea2(elm->edu[0].nod->pos,
                                elm->edu[1].nod->pos,
                                elm->edu[2].nod->pos);
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
  *dstNod0 = elm->face[0].edu[0].nod;
  *dstNod1 = elm->face[0].edu[1].nod;
  *dstNod2 = elm->face[0].edu[2].nod;
  *dstNod3 = elm->face[1].edu[1].nod;
}

/*!
* \ingroup	WlzMesh
* \brief	Gets the three nodes of a 2D element.
* \param	elm			Given mesh element.
* \param	dstNod0			First destination pointer for node.
* \param	dstNod1			Second destination pointer for node.
* \param	dstNod2			Third destination pointer for node.
*/
void		WlzCMeshElmGetNodes2D(WlzCMeshElm2D *elm,
				      WlzCMeshNod2D **dstNod0,
				      WlzCMeshNod2D **dstNod1,
				      WlzCMeshNod2D **dstNod2)
{
  *dstNod0 = elm->edu[0].nod;
  *dstNod1 = elm->edu[1].nod;
  *dstNod2 = elm->edu[2].nod;
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
  WlzCMeshEdgU2D *edu0,
  		*edu1;
  WlzDVertex2   nP,
                sP;

  sP.vtX = sP.vtY = 0.0;
  if((nod->idx >= 0) &&
     (doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0)))
  {
    edu0 = edu1 = nod->edu;
    do
    {
      oNod = edu1->next->nod;
      nP = *(vBuf + oNod->idx);
      WLZ_VTX_2_ADD(sP, sP, nP);
      ++nN;
      edu1 = edu1->nnxt;
    } while(edu0 != edu1);
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
  WlzCMeshEdgU3D *edu0,
  		*edu1;
  WlzDVertex3   nP,
                sP;

  sP.vtX = sP.vtY = sP.vtZ = 0.0;
  if((nod->idx >= 0) &&
     (doBnd || ((nod->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0)))
  {
    edu0 = edu1 = nod->edu;
    do
    {
      oNod = edu1->next->nod;
      nP = *(vBuf + oNod->idx);
      WLZ_VTX_3_ADD(sP, sP, nP);
      ++nN;
      edu1 = edu1->nnxt;
    } while(edu0 != edu1);
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

/*!
* \return	Copy of the given mesh.
* \ingroup	WlzMesh
* \brief	Creates a copy of the given constrained mesh in which
* 		all deleted entities have been squeezed out. In the
* 		mew mesh the valid mesh entities are contiguous and start
* 		with index 0.
* 		As well as copying the given mesh this function will
* 		also copy data associated with the nodes of the mesh.
* 		The associated data are held in AlcVector data structures
* 		and are indexed using the node index.
* 		If (datSz == 0) || (newDat == NULL) || (gvnDat == NULL)
* 		then the associated data will be ignored.
* \param	gvnMesh			Given mesh.
* \param	datSz			Size of associated datum.
* \param	newDat			Destination pointer for the copied
* 					associated data, may be NULL.
* \param	gvnDat			Given associated data, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMeshP	WlzCMeshCopy(WlzCMeshP gvnMesh, size_t datSz,
				AlcVector **newDat, AlcVector *gvnDat,
				WlzErrorNum *dstErr)
{
  WlzCMeshP	newMesh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  newMesh.v = NULL;
  if(gvnMesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(gvnMesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        newMesh.m2 = WlzCMeshCopy2D(gvnMesh.m2, datSz, newDat, gvnDat,
				    &errNum);
        break;
      case WLZ_CMESH_TET3D:
        newMesh.m3 = WlzCMeshCopy3D(gvnMesh.m3, datSz, newDat, gvnDat,
				    &errNum);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(newMesh);
}

/*!
* \return	Copy of the given mesh.
* \ingroup	WlzMesh
* \brief	Creates a copy of the given constrained mesh in which
* 		all deleted entities have been squeezed out. In the
* 		mew mesh the valid mesh entities are contiguous and start
* 		with index 0.
* 		As well as copying the given mesh this function will
* 		also copy data associated with the nodes of the mesh.
* 		The associated data are held in AlcVector data structures
* 		and are indexed using the node index.
* 		If (datSz == 0) || (newDat == NULL) || (gvnDat == NULL)
* 		then the associated data will be ignored.
* \param	gvnMesh			Given mesh.
* \param	datSz			Size of associated datum.
* \param	newDat			Destination pointer for the copied
* 					associated data, may be NULL.
* \param	gvnDat			Given associated data, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMesh2D	*WlzCMeshCopy2D(WlzCMesh2D *gvnMesh, size_t datSz,
				AlcVector **newDat, AlcVector *gvnDat,
				WlzErrorNum *dstErr)
{
  int		idE,
  		idN;
  void		*ascP;
  AlcVector	*prvDat = NULL;
  WlzIVertex2	dumGrdPos;
  WlzCMeshNod2D	*dumNod;
  WlzCMeshNod2D	*gvnNodes[3],
  		*newNodes[3];
  WlzCMeshElm2D	*gvnElm;
  WlzCMesh2D	*newMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gvnMesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gvnMesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((datSz > 0) && (newDat != NULL) && (gvnDat != NULL))
  {
    if((prvDat = AlcVectorNew(1, datSz, gvnDat->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) &&
     (gvnMesh->res.nod.numEnt > 0) && (gvnMesh->res.elm.numEnt > 0))
  {
    idE = 0;
    newMesh = WlzCMeshNew2D(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      newMesh->bBox = gvnMesh->bBox;
      errNum = WlzCMeshReassignBuckets2D(newMesh, gvnMesh->res.nod.numEnt);
    }
    while((errNum == WLZ_ERR_NONE) && (idE < gvnMesh->res.elm.maxEnt))
    {
      gvnElm = (WlzCMeshElm2D *)AlcVectorItemGet(gvnMesh->res.elm.vec, idE);
      if(gvnElm->idx >= 0)
      {
	/* Copy Nodes. */
        WlzCMeshElmGetNodes2D(gvnElm, gvnNodes + 0, gvnNodes + 1,
	                      gvnNodes + 2);
	idN = 0;
	do
	{
	  if(WlzCMeshLocateNod2D(newMesh, gvnNodes[idN]->pos,
	                         &dumGrdPos, &dumNod, newNodes + idN) == 0)

	  {
	    newNodes[idN] = WlzCMeshNewNod2D(newMesh, gvnNodes[idN]->pos,
	                                     &errNum);
	    /* Copy associated data. */
	    if((errNum == WLZ_ERR_NONE) && (prvDat != NULL))
	    {
	      if((ascP = AlcVectorExtendAndGet(prvDat,
	                                       newNodes[idN]->idx)) != NULL)
	      {
	        memcpy(ascP, AlcVectorItemGet(prvDat, gvnNodes[idN]->idx),
		       datSz);
	      }
	    }
	  }
	} while((errNum == WLZ_ERR_NONE) && (++idN < 3));
	/* Copy element. */
	if(errNum == WLZ_ERR_NONE)
	{
	  (void )WlzCMeshNewElm2D(newMesh,
	                          newNodes[0], newNodes[1], newNodes[2],
				  &errNum);
	}
      }
      ++idE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvDat != NULL)
    {
      *newDat = prvDat;
    }
  }
  else
  {
    if(newMesh != NULL)
    {
      WlzCMeshFree2D(newMesh);
      newMesh = NULL;
    }
    if(prvDat != NULL)
    {
      (void )AlcVectorFree(prvDat);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(newMesh);
}

/*!
* \return	Copy of the given mesh.
* \ingroup	WlzMesh
* \brief	Creates a copy of the given constrained mesh in which
* 		all deleted entities have been squeezed out. In the
* 		mew mesh the valid mesh entities are contiguous and start
* 		with index 0.
*		As well as copying the given mesh this function will
*               also copy data associated with the nodes of the mesh.
*               The associated data are held in AlcVector data structures
*               and are indexed using the node index.
*               If (datSz == 0) || (newDat == NULL) || (gvnDat == NULL)
*               then the associated data will be ignored.
* \param	gvnMesh			Given mesh.
* \param        datSz                   Size of associated datum.
* \param        newDat                  Destination pointer for the copied
*                                       associated data, may be NULL.
* \param        gvnDat                  Given associated data, may be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzCMesh3D	*WlzCMeshCopy3D(WlzCMesh3D *gvnMesh, size_t datSz,
                                AlcVector **newDat, AlcVector *gvnDat,
                                WlzErrorNum *dstErr)
{
  int		idE,
  		idN;
  void          *ascP;
  AlcVector     *prvDat = NULL;
  WlzIVertex3	dumGrdPos;
  WlzCMeshNod3D	*dumNod;
  WlzCMeshNod3D	*gvnNodes[4],
  		*newNodes[4];
  WlzCMeshElm3D	*gvnElm;
  WlzCMesh3D	*newMesh = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gvnMesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(gvnMesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if((datSz > 0) && (newDat != NULL) && (gvnDat != NULL))
  {
    if((prvDat = AlcVectorNew(1, datSz, gvnDat->blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) &&
     (gvnMesh->res.nod.numEnt > 0) && (gvnMesh->res.elm.numEnt > 0))
  {
    idE = 0;
    newMesh = WlzCMeshNew3D(&errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      newMesh->bBox = gvnMesh->bBox;
      errNum = WlzCMeshReassignBuckets3D(newMesh, gvnMesh->res.nod.numEnt);
    }
    while((errNum == WLZ_ERR_NONE) && (idE < gvnMesh->res.elm.maxEnt))
    {
      gvnElm = (WlzCMeshElm3D *)AlcVectorItemGet(gvnMesh->res.elm.vec, idE);
      if(gvnElm->idx >= 0)
      {
	/* Copy Nodes. */
        WlzCMeshElmGetNodes3D(gvnElm, gvnNodes + 0, gvnNodes + 1,
	                      gvnNodes + 2, gvnNodes + 3);
	idN = 0;
	do
	{
	  if(WlzCMeshLocateNod3D(newMesh, gvnNodes[idN]->pos,
	                         &dumGrdPos, &dumNod, newNodes + idN) == 0)

	  {
	    newNodes[idN] = WlzCMeshNewNod3D(newMesh, gvnNodes[idN]->pos,
	                                     &errNum);
            /* Copy associated data. */
            if((errNum == WLZ_ERR_NONE) && (prvDat != NULL))
            {
              if((ascP = AlcVectorExtendAndGet(prvDat,
                                               newNodes[idN]->idx)) != NULL)
              {
                memcpy(ascP, AlcVectorItemGet(prvDat, gvnNodes[idN]->idx),
		       datSz);
              }
            }
	  }
	} while((errNum == WLZ_ERR_NONE) && (++idN < 4));
	/* Copy element. */
	if(errNum == WLZ_ERR_NONE)
	{
	  (void )WlzCMeshNewElm3D(newMesh,
	                          newNodes[0], newNodes[1], newNodes[2],
				  newNodes[3], &errNum);
	}
      }
      ++idE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(prvDat != NULL)
    {
      *newDat = prvDat;
    }
  }
  else
  {
    if(newMesh != NULL)
    {
      WlzCMeshFree3D(newMesh);
      newMesh = NULL;
    }
    if(prvDat != NULL)
    {
      (void )AlcVectorFree(prvDat);
    }
  }
  if(dstErr != NULL)
  {
    *dstErr = errNum;
  }
  return(newMesh);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzMesh
* \brief        Reorders nodes in any elements which have negative
		area.
* \param        mesh                     The given constrained mesh.
*/
WlzErrorNum     WlzCMeshFixNegativeElms(WlzCMeshP mesh)
{
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if(mesh.v == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_TRI2D:
        errNum = WlzCMeshFixNegativeElms2D(mesh.m2);
        break;
      case WLZ_CMESH_TET3D:
        errNum = WlzCMeshFixNegativeElms3D(mesh.m3);
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
        break;
    }
  }
  return(errNum);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzMesh
* \brief        Reorders nodes in any elements which have negative
		area.
* \param        mesh                     Given 2D constrained mesh.
*/
WlzErrorNum     WlzCMeshFixNegativeElms2D(WlzCMesh2D *mesh)
{
  int		idE,
  		idG;
  double	sA2;
  WlzCMeshElm2D *elm;
  WlzCMeshEdgU2D *edu0,
  		*edu1,
		*edu2;
  WlzCMeshNod2D	*nod[3];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* TODO Check this function. */
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TRI2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	WlzCMeshElmGetNodes2D(elm, nod + 0, nod + 1, nod + 2);
	sA2 = WlzGeomTriangleSnArea2(nod[0]->pos, nod[1]->pos, nod[2]->pos);
	if(fabs(sA2) < WLZ_MESH_TOLERANCE_SQ)
	{
	  /* Might need to remove the element and repair the mesh here? */
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else if(sA2 < WLZ_MESH_TOLERANCE_SQ)
	{
	  /* Unlink the edge uses from the rest of the mesh. */
	  for(idG = 0; idG < 3; ++idG)
	  {
	    edu0 = elm->edu + idG;
	    if((edu0->opp != NULL) && (edu0->opp->opp != NULL) &&
	       (edu0->opp->opp->elm == elm))
            {
	      edu0->opp->opp = NULL;
	    }
	    if(edu0 == edu0->nnxt)
	    {
	      edu0->nnxt = NULL;
	    }
	    else
	    {
	      edu1 = edu0;
	      while((edu2 = edu1->nnxt) != edu0)
	      {
	        edu1 = edu2;
	      }
	      edu1->nnxt = edu0->nnxt;
	      if(edu0->nod->edu == edu0)
	      {
	        edu0->nod->edu = edu1;
	      }
	    }
	  }
	  /* Create connnectivities with the element and to the rest of the
	   * mesh. */
	  errNum = WlzCMeshSetElm2D(mesh, elm, nod[1], nod[0], nod[2]);
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return       Woolz error code.
* \ingroup      WlzMesh
* \brief        Reorders nodes in any elements which have negative
		area.
* \param        mesh                     Given 3D constrained mesh.
*/
WlzErrorNum     WlzCMeshFixNegativeElms3D(WlzCMesh3D *mesh)
{
  int		idE,
		idF,
  		idG;
  double	sV6;
  WlzCMeshElm3D *elm;
  WlzCMeshFace	*fce0;
  WlzCMeshEdgU3D *edu0,
  		*edu1,
		*edu2;
  WlzCMeshNod3D	*nod[4];
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  /* TODO Check this function. */
  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_TET3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	WlzCMeshElmGetNodes3D(elm, nod + 0, nod + 1, nod + 2, nod + 3);
	sV6 = WlzGeomTetraSnVolume6(nod[0]->pos, nod[1]->pos, nod[2]->pos,
				    nod[3]->pos);
	if(fabs(sV6) < WLZ_MESH_TOLERANCE_SQ)
	{
	  /* Might need to remove the element and repair the mesh here? */
	  errNum = WLZ_ERR_DOMAIN_DATA;
	}
	else if(sV6 < WLZ_MESH_TOLERANCE_SQ)
	{
	  for(idF = 0; idF < 4; ++idF)
	  {
	    fce0 = elm->face + idF;
	    /* Unlink edge uses from the nodes. */
	    for(idG = 0; idG < 3; ++idG)
	    {
	      edu0 = fce0->edu + idG;
	      if(edu0 == edu0->nnxt)
	      {
		edu0->nnxt = NULL;
	      }
	      else
	      {
		edu1 = edu0;
		edu2 = edu1->nnxt;
		while(edu2 != edu0)
		{
		  edu1 = edu2;
		  edu2 = edu2->nnxt;
		}
		edu1->nnxt = edu2->nnxt;
		edu1->nod->edu = edu1;
	      }
	    }
	    /* Unlink face. Need to make sure that the opp - opp link is back
	     * to this element and not some other that will replace it. */
	    if((fce0->opp != NULL) && (fce0->opp->opp != NULL) &&
	       (fce0->opp->opp->elm == elm))
	    {
	      fce0->opp->opp = NULL;
	    }
	  }
	  /* Create connnectivities with the element and to the rest of the
	   * mesh. */
	  errNum = WlzCMeshSetElm3D(mesh, elm, nod[1], nod[0], nod[2],
	  			    nod[3]);
	}
      }
    }
  }
  return(errNum);
}
