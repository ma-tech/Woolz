#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFMesh_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFMesh.c
* \author       Bill Hill
* \date         September 2009
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
* \brief	Functions for reading and writting Woolz objects to and from
* 		the the NETGEN tetrahedral mesh file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjMeshRec(
				  FILE *fP,
				  char *buf,
				  int bufLen);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		NETGEN mesh file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjMesh(FILE *fP, WlzErrorNum *dstErr)
{
  int		cnt,
  		idE,
  		idN,
  		nElm = 0,
		nNod = 0;
  double	vol;
  WlzDBox3	bBox;
  char		*str;
  WlzDVertex3   *vBuf = NULL;
  WlzCMesh3D	*mesh = NULL;
  WlzObject	*obj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		fBuf[5];
  char		cBuf[256];
  WlzCMeshNod3D	*nBuf[4];

  dom.core = NULL;
  val.core = NULL;
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    /* Optional records may be comments starting with a '#' but apart
     * from these the first line should have a single integer specifying
     * the number of nodes. */
    cnt = 0;
    do
    {
      if((str = WlzEffReadObjMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((cnt = sscanf(str, "%d", &nNod)) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (cnt == 0));
  }
  /* Check for reasonable number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(nNod <= 0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Read in the node positions into a temporary buffer computing their
   * bounding box and then create the nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((vBuf = (WlzDVertex3 *)AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idN = 0;
    while(idN < nNod)
    {
      if((str = WlzEffReadObjMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      else if(*str != '#')
      {
	if(sscanf(str, "%lg %lg %lg",
	          &(vBuf[idN].vtX), &(vBuf[idN].vtY), &(vBuf[idN].vtZ)) == 3)
	{
	  if(idN == 0)
	  {
	    bBox.xMin = bBox.xMax = vBuf[idN].vtX;
	    bBox.yMin = bBox.yMax = vBuf[idN].vtY;
	    bBox.zMin = bBox.zMax = vBuf[idN].vtZ;
	  }
	  else
	  {
	    if(vBuf[idN].vtX < bBox.xMin)
	    {
	      bBox.xMin = vBuf[idN].vtX;
	    }
	    else if(vBuf[idN].vtX > bBox.xMax)
	    {
	      bBox.xMax = vBuf[idN].vtX;
	    }
	    if(vBuf[idN].vtY < bBox.yMin)
	    {
	      bBox.yMin = vBuf[idN].vtY;
	    }
	    else if(vBuf[idN].vtY > bBox.yMax)
	    {
	      bBox.yMax = vBuf[idN].vtY;
	    }
	    if(vBuf[idN].vtZ < bBox.zMin)
	    {
	      bBox.zMin = vBuf[idN].vtZ;
	    }
	    else if(vBuf[idN].vtZ > bBox.zMax)
	    {
	      bBox.zMax = vBuf[idN].vtZ;
	    }
	  }
          ++idN;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
      }
    }
  }
  /* Read the number of elements and check that it's reasonable. */
  /* Create a new 3D constrained mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    do
    {
      if((str = WlzEffReadObjMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((cnt = sscanf(str, "%d", &nElm)) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (cnt == 0));
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(nElm <= 0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Make a new mesh then allocate it's nodes and elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew3D(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mesh->bBox = bBox;
    errNum = WlzCMeshReassignGridCells3D(mesh, nNod);
  }
  /* Create the mesh nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idN = 0; idN < nNod; ++idN)
    {
       nBuf[0] = WlzCMeshNewNod3D(mesh, vBuf[idN], NULL);
       nBuf[0]->pos = vBuf[idN];
       nBuf[0]->flags = 0;
    }
  }
  AlcFree(vBuf);
  /* Read the node indices for the elements and create the elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    idE = 0;
    while((errNum == WLZ_ERR_NONE) && (idE < nElm))
    {
      if((str = WlzEffReadObjMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if(sscanf(str, "%d %d %d %d %d",
	          fBuf + 0, fBuf + 1, fBuf + 2, fBuf + 3, fBuf + 4) != 5)
	{
	  break;
	}
	else if((fBuf[0] != 1) ||
	        (fBuf[2] < 2) || (fBuf[2] > nNod + 1) ||
		(fBuf[3] < 2) || (fBuf[3] > nNod + 1) ||
		(fBuf[4] < 2) || (fBuf[4] > nNod + 1))

	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  for(idN = 0; idN < 4; ++idN)
	  {
	    nBuf[idN] = (WlzCMeshNod3D *)
			AlcVectorItemGet(mesh->res.nod.vec, fBuf[idN + 1] - 2);
	  }
	  vol = WlzGeomTetraSnVolume6(nBuf[0]->pos, nBuf[1]->pos,
	                              nBuf[2]->pos, nBuf[3]->pos);
	  if(vol < 0)
	  {
	    (void )WlzCMeshNewElm3D(mesh, nBuf[0], nBuf[1], nBuf[3], nBuf[2],
		                    &errNum);
	  }
	  else
	  {
	    (void )WlzCMeshNewElm3D(mesh, nBuf[0], nBuf[1], nBuf[2], nBuf[3],
		                    &errNum);
	  }
	}
        ++idE;
      }
    }
  }
  /* Compute maximum edge length and then create the Woolz object. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzCMeshUpdateMaxSqEdgLen3D(mesh);
    dom.cm3 = mesh;
    obj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzCMeshFree3D(mesh);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Start of string having skipped past possible white space
* 		characters.
* \ingroup	WlzExtFF
* \brief	Gets the next record from a file.
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
*/
static char	*WlzEffReadObjMeshRec(FILE *fP, char *buf, int bufLen)
{
  char		*str = NULL;

  if((str = fgets(buf, bufLen, fP)) != NULL)
  {
    buf[bufLen - 1] = '\0';
    while((*str != '\0') && isspace(*str))
    {
      ++str;
    }
  }
  return(str);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given stream using the
* 		NETGEN tetrahedral mesh file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjMesh(FILE *fP, WlzObject *obj)
{
  int		cnt,
  		idE,
  		idN,
  		nElm = 0,
		nNod = 0;
  int		*nodTbl = NULL;
  WlzCMesh3D	*mesh;
  WlzCMeshElm3D	*elm;
  WlzCMeshNod3D	*nod;
  WlzCMeshNod3D	*nBuf[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    mesh = obj->domain.cm3;
    if(mesh->type != WLZ_CMESH_3D)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  /* Output the number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    nElm = mesh->res.elm.numEnt;
    nNod = mesh->res.nod.numEnt;

    if(fprintf(fP, "%d\n", nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the node positions while building a node table to avoid deleted
   * nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodTbl = (int *)AlcMalloc(sizeof(int) * nNod)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
        nodTbl[idN] = cnt++;
	if(fprintf(fP, " %lg %lg %lg\n",
	           nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  /* Output the number of elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "%d\n", nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the elements defined by the indices of their nodes. note that
   * the node indices have an offset of 2 for some reason. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	nBuf[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
	nBuf[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
	nBuf[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	nBuf[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
	if(fprintf(fP, " %d %d %d %d %d\n",
		   1,
		   nodTbl[nBuf[0]->idx] + 2, nodTbl[nBuf[1]->idx] + 2,
		   nodTbl[nBuf[2]->idx] + 2, nodTbl[nBuf[3]->idx] + 2) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  AlcFree(nodTbl);
  return(errNum);
}
