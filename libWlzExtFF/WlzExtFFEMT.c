#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFEMT_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFEMT.c
* \author       Bill Hill
* \date         March 2011
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
* \brief	Functions for reading and writting Woolz objects to and from
* 		the the Netgen neutral (tetrahedral) mesh file format.
* \ingroup	WlzExtFF
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjEMTRec(
				  FILE *fP,
				  char *buf,
				  int bufLen);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		GRUMMP vmesh tetrahedral mesh file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjEMT(FILE *fP, WlzErrorNum *dstErr)
{
  int		nElm = 0,
		nNod = 0;
  char		*str;
  WlzDVertex3   *vBuf = NULL;
  WlzCMesh3D	*mesh = NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		cBuf[256];

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    int		idC = 0;

    /* Optional records may be comments starting with a '#' but apart
     * from these the first line should have a single integer specifying
     * the number of nodes. */
    do
    {
      if((str = WlzEffReadObjEMTRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((idC = sscanf(str, "%d", &nNod)) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (idC != 1));
  }
  /* Check for reasonable number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(nNod <= 0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Create a new 3D constrained mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew3D(&errNum);
  }
  /* Read in the node positions into a temporary buffer computing their
   * bounding box and then create the nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((vBuf = AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN = 0;

    while(idN < nNod)
    {
      if((str = WlzEffReadObjEMTRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      else if(*str != '#')
      {
	if(sscanf(str, "%lg %lg %lg",
		  &(vBuf[idN].vtX), &(vBuf[idN].vtY), &(vBuf[idN].vtZ)) != 3)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	++idN;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mesh->bBox = WlzBoundingBoxVtx3D(nNod, vBuf, NULL);
    if(AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshReassignGridCells3D(mesh, nNod);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN;

    for(idN = 0; idN < nNod; ++idN)
    {
       (void )WlzCMeshNewNod3D(mesh, vBuf[idN], NULL);
    }
  }
  AlcFree(vBuf);
  /* Read the number of elements in the mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idC = 0;

    /* Optional records may be comments starting with a '#' but apart
     * from these the first line should have a single integer specifying
     * the number of nodes. */
    do
    {
      if((str = WlzEffReadObjEMTRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((idC = sscanf(str, "%d", &nElm)) != 1)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (idC != 1));
  }
  /* Check for reasonable number of elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(nElm <= 0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Allocate room for the elements in the mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Read the elements adding them to the mesh. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE = 0;

    while(idE < nElm)
    {
      if((str = WlzEffReadObjEMTRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      else if(*str != '#')
      {
	int	dummyI;
	int	nIdx[4];
	WlzCMeshNod3D *nBuf[4];

	if((sscanf(str, "%d %d %d %d %d",
		   &dummyI,
		   &(nIdx[0]), &(nIdx[1]), &(nIdx[2]), &(nIdx[3])) != 5) ||
	   (nIdx[0] <= 0) || (nIdx[0] > nNod) ||
	   (nIdx[1] <= 0) || (nIdx[1] > nNod) ||
	   (nIdx[2] <= 0) || (nIdx[2] > nNod) ||
	   (nIdx[3] <= 0) || (nIdx[3] > nNod))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	else
	{
	  int	idN;

	  for(idN = 0; idN < 4; ++idN)
	  {
	    nBuf[idN] = (WlzCMeshNod3D *)
			AlcVectorItemGet(mesh->res.nod.vec, nIdx[idN] - 1);
	  }
	  (void )WlzCMeshNewElm3D(mesh, nBuf[0], nBuf[1], nBuf[3], nBuf[2], 1,
				  &errNum);
	  if(errNum != WLZ_ERR_NONE)
	  {
	    break;
	  }
	  ++idE;
	}
      }
    }
  }
  /* Ignore the boundary faces in the file, they're redundant information. */
  if(errNum == WLZ_ERR_NONE)
  {
    WlzDomain	dom;
    WlzValues	val;

    val.core = NULL;
    dom.cm3 = mesh;
    WlzCMeshDelUnusedNodes3D(mesh);
    WlzCMeshUpdateMaxSqEdgLen3D(mesh);
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
static char	*WlzEffReadObjEMTRec(FILE *fP, char *buf, int bufLen)
{
  char		*str = NULL;

  if((str = fgets(buf, bufLen, fP)) != NULL)
  {
    buf[bufLen - 1] = '\0';
    while((*str != '\0') && isspace(*str))
    {
      ++str;
    }
    if(*str == '\0')
    {
      buf[0] = '#';
      str = buf;
    }
  }
  return(str);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given stream using the
* 		Netgen neutral mesh file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjEMT(FILE *fP, WlzObject *obj)
{
  int		nBFce = 0,
     		nElm = 0,
		nNod = 0;
  int		*nodTbl = NULL;
  WlzCMesh3D	*mesh;
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
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodTbl = (int *)
		 AlcMalloc(mesh->res.nod.maxEnt * sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute the number of boundary faces while building an element table to
   * avoid deleted elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idC = 0,
    		idE;

    nNod = mesh->res.nod.numEnt;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	int	idF;

	for(idF = 0; idF < 4; ++idF)
	{
	  if((elm->face[idF].opp == NULL) ||
	     (elm->face[idF].opp == &(elm->face[idF])))
	  
	  {
	    ++nBFce;
	  }
	}
        ++nElm;
      }
    }
  }
  /* Output the number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "%d\n", nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the node positions while building a node table to avoid deleted
   * nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idC = 0,
    		idN;

    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      WlzCMeshNod3D *nod;

      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
        nodTbl[idN] = ++idC;
	if(fprintf(fP, "  %lg %lg %lg\n",
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
  /* Output the elements using the node indices. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;
    WlzCMeshNod3D *nBuf[4];

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	WlzCMeshElmGetNodes3D(elm,
			      &(nBuf[0]), &(nBuf[1]), &(nBuf[2]), &(nBuf[3]));
        if(fprintf(fP, "  1 %d %d %d %d\n",
	           nodTbl[(nBuf[0])->idx],
		   nodTbl[(nBuf[1])->idx],
	           nodTbl[(nBuf[2])->idx],
		   nodTbl[(nBuf[3])->idx]) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  /* Output the number of boundary faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "%d\n", nBFce) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the boundary faces defined by the node indices. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      WlzCMeshElm3D *elm;

      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	int	idF;

        for(idF = 0; idF < 4; ++idF)
	{
	  WlzCMeshFace *fce;

	  fce = elm->face + idF;
	  if((fce->opp == NULL) || (fce->opp == fce))
	  {
	    WlzCMeshEdgU3D *edu;

	    edu = fce->edu;
            if(fprintf(fP, "  1 %d %d %d\n",
		       nodTbl[edu->nod->idx], 
		       nodTbl[edu->next->nod->idx], 
		       nodTbl[edu->next->next->nod->idx]) <= 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      break;
	    }
	  }
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  AlcFree(nodTbl);
  return(errNum);
}
