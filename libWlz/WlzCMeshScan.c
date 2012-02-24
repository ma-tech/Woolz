#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshScan_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzCMeshScan.c
* \author       Bill Hill
* \date         October 2009
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
* \brief	Iterators for scanning through CMesh based objects and
* 		domains.
* \ingroup	WlzMesh
*/

#include <Wlz.h>

/*!
* \return       Next node in mesh or NULL one error of if there are no more
* 		nodes in the mesh.
* \ingroup      WlzMesh
* \brief        Returns the next node in the 2D mesh. Repeatedly calling
* 		this function will itterate through all nodes of the mesh.
* 		Before starting to itterate through the nodes of a mesh
* 		the index should be initialised to zero (or the index
* 		value from which to start). On return the index will have
* 		been incremented.
* \param        mesh                    The given 2D mesh, must be valid.
* \param        idx                     Pointer to the index of the next
*                                       node which should be set to zero
*                                       at the start of an itteration,
*                                       must not be NULL.
* \param        all			Non-zero if invalid nodes are to
* 					be returned. If zero non-valid
* 					nodes (eg deleted nodes) are
* 					skiped.
*/
WlzCMeshNod2D   *WlzCMeshNextNod2D(WlzCMesh2D *mesh, int *idx, int all)
{
  int           maxIdx;
  AlcVector     *vec;
  WlzCMeshNod2D *nod = NULL;

  if((mesh != NULL) && (idx != NULL))
  {
    vec = mesh->res.nod.vec;
    maxIdx = mesh->res.nod.maxEnt;
    if(all == 0)
    {
      if(*idx < maxIdx)
      {
        nod = (WlzCMeshNod2D *)AlcVectorItemGet(vec, (*idx)++);
      }
    }
    else
    {
      while(*idx < maxIdx)
      {
        nod = (WlzCMeshNod2D *)AlcVectorItemGet(vec, (*idx)++);
        if(nod->idx >= 0)
        {
          break;
        }
      }
    }
  }
  return(nod);
}

/*!
* \return       Next node in mesh or NULL on error or if there are no more
* 		nodes in the mesh.
* \ingroup      WlzMesh
* \brief        Returns the next node in the 3D mesh. Repeatedly calling
* 		this function will itterate through all nodes of the mesh.
* 		Before starting to itterate through the nodes of a mesh
* 		the index should be initialised to zero (or the index
* 		value from which to start). On return the index will have
* 		been incremented.
* \param        mesh                    The given 3D mesh, must be valid.
* \param        idx                     Pointer to the index of the next
*                                       node which should be set to zero
*                                       at the start of an itteration,
*                                       must not be NULL.
* \param        all			Non-zero if invalid nodes are to
* 					be returned. If zero non-valid
* 					nodes (eg deleted nodes) are
* 					skiped.
*/
WlzCMeshNod3D   *WlzCMeshNextNod3D(WlzCMesh3D *mesh, int *idx, int all)
{
  int           maxIdx;
  AlcVector     *vec;
  WlzCMeshNod3D *nod = NULL;

  if((mesh != NULL) && (idx != NULL))
  {
    vec = mesh->res.nod.vec;
    maxIdx = mesh->res.nod.maxEnt;
    if(all == 0)
    {
      if(*idx < maxIdx)
      {
        nod = (WlzCMeshNod3D *)AlcVectorItemGet(vec, (*idx)++);
      }
    }
    else
    {
      while(*idx < maxIdx)
      {
        nod = (WlzCMeshNod3D *)AlcVectorItemGet(vec, (*idx)++);
        if(nod->idx >= 0)
        {
          break;
        }
      }
    }
  }
  return(nod);
}

/*!
* \return       Next element in mesh or NULL one error of if there are no
* 		more elements in the mesh.
* \ingroup      WlzMesh
* \brief        Returns the next element in the 2D mesh. Repeatedly calling
* 		this function will itterate through all elements of the mesh.
* 		Before starting to itterate through the elements of a mesh
* 		the index should be initialised to zero (or the index
* 		value from which to start). On return the index will have
* 		been incremented.
* \param        mesh                    The given 2D mesh, must be valid.
* \param        idx                     Pointer to the index of the next
*                                       element which should be set to zero
*                                       at the start of an itteration,
*                                       must not be NULL.
* \param        all			Non-zero if invalid elements are to
* 					be returned. If zero non-valid
* 					elements (eg deleted elements) are
* 					skiped.
*/
WlzCMeshElm2D   *WlzCMeshNextElm2D(WlzCMesh2D *mesh, int *idx, int all)
{
  int           maxIdx;
  AlcVector     *vec;
  WlzCMeshElm2D *elm = NULL;

  if((mesh != NULL) && (idx != NULL))
  {
    vec = mesh->res.elm.vec;
    maxIdx = mesh->res.elm.maxEnt;
    if(all == 0)
    {
      if(*idx < maxIdx)
      {
        elm = (WlzCMeshElm2D *)AlcVectorItemGet(vec, (*idx)++);
      }
    }
    else
    {
      while(*idx < maxIdx)
      {
        elm = (WlzCMeshElm2D *)AlcVectorItemGet(vec, (*idx)++);
        if(elm->idx >= 0)
        {
          break;
        }
      }
    }
  }
  return(elm);
}

/*!
* \return       Next element in mesh or NULL on error or if there are no
* 		more elements in the mesh.
* \ingroup      WlzMesh
* \brief        Returns the next element in the 3D mesh. Repeatedly calling
* 		this function will itterate through all elements of the mesh.
* 		Before starting to itterate through the elements of a mesh
* 		the index should be initialised to zero (or the index
* 		value from which to start). On return the index will have
* 		been incremented.
* \param        mesh                    The given 3D mesh, must be valid.
* \param        idx                     Pointer to the index of the next
*                                       element which should be set to zero
*                                       at the start of an itteration,
*                                       must not be NULL.
* \param        all			Non-zero if invalid elements are to
* 					be returned. If zero non-valid
* 					elements (eg deleted elements) are
* 					skiped.
*/
WlzCMeshElm3D   *WlzCMeshNextElm3D(WlzCMesh3D *mesh, int *idx, int all)
{
  int           maxIdx;
  AlcVector     *vec;
  WlzCMeshElm3D *elm = NULL;

  if((mesh != NULL) && (idx != NULL))
  {
    vec = mesh->res.elm.vec;
    maxIdx = mesh->res.elm.maxEnt;
    if(all == 0)
    {
      if(*idx < maxIdx)
      {
        elm = (WlzCMeshElm3D *)AlcVectorItemGet(vec, (*idx)++);
      }
    }
    else
    {
      while(*idx < maxIdx)
      {
        elm = (WlzCMeshElm3D *)AlcVectorItemGet(vec, (*idx)++);
        if(elm->idx >= 0)
        {
          break;
        }
      }
    }
  }
  return(elm);
}
