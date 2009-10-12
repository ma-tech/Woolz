#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFVMesh_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFVMesh.c
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
* 		the the GRUMMP VMESH tetrahedral mesh file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjVMeshRec(
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
WlzObject	*WlzEffReadObjVMesh(FILE *fP, WlzErrorNum *dstErr)
{
  int		cnt,
		id0,
		id1,
  		idE,
  		idF,
  		idN,
  		nElm = 0,
  		nFce = 0,
		nBFce = 0,
		nNod = 0;
  double	vol;
  char		*str;
  int		*eBufP,
  		*eBuf = NULL;
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
     * from these the first line should have four integers specifying
     * the number of elements, faces, boundary faces and nodes. */
    cnt = 0;
    do
    {
      if((str = WlzEffReadObjVMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else if(*str != '#')
      {
	if((cnt = sscanf(str, "%d %d %d %d",
	                 &nElm, &nFce, &nBFce, &nNod)) != 4)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (cnt != 4));
  }
  /* Check for reasonable number of elements, faces, etc.... */
  if(errNum == WLZ_ERR_NONE)
  {
    if((nElm <= 0) || (nFce <= 0) || (nBFce <= 0) || (nNod <= 0))
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
    if(((vBuf = AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.nod.vec, nNod) == NULL) ||
       (AlcVectorExtendAndGet(mesh->res.elm.vec, nElm) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idN = 0;
    while(idN < nNod)
    {
      if((str = WlzEffReadObjVMeshRec(fP, cBuf, 256)) == NULL)
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
	    mesh->bBox.xMin = mesh->bBox.xMax = vBuf[idN].vtX;
	    mesh->bBox.yMin = mesh->bBox.yMax = vBuf[idN].vtY;
	    mesh->bBox.zMin = mesh->bBox.zMax = vBuf[idN].vtZ;
	  }
	  else
	  {
	    if(vBuf[idN].vtX < mesh->bBox.xMin)
	    {
	      mesh->bBox.xMin = vBuf[idN].vtX;
	    }
	    else if(vBuf[idN].vtX > mesh->bBox.xMax)
	    {
	      mesh->bBox.xMax = vBuf[idN].vtX;
	    }
	    if(vBuf[idN].vtY < mesh->bBox.yMin)
	    {
	      mesh->bBox.yMin = vBuf[idN].vtY;
	    }
	    else if(vBuf[idN].vtY > mesh->bBox.yMax)
	    {
	      mesh->bBox.yMax = vBuf[idN].vtY;
	    }
	    if(vBuf[idN].vtZ < mesh->bBox.zMin)
	    {
	      mesh->bBox.zMin = vBuf[idN].vtZ;
	    }
	    else if(vBuf[idN].vtZ > mesh->bBox.zMax)
	    {
	      mesh->bBox.zMax = vBuf[idN].vtZ;
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
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzCMeshReassignGridCells3D(mesh, nNod);
  }
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
  /* Build a table of nodes for each element. Each element to have 7 ints:
   * 0   - number of nodes found for this element so far
   * 1-4 - indices of the four nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((eBuf = AlcCalloc(5 * nElm, sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idF = 0;
    while(idF < nFce)
    {
      if((str = WlzEffReadObjVMeshRec(fP, cBuf, 256)) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      else if(*str != '#')
      {
	if(sscanf(str, "%d %d %d %d %d",
	          fBuf + 0, fBuf + 1, fBuf + 2, fBuf + 3, fBuf + 4) == 5)
	{
	  if((fBuf[0] >= nElm) ||
	     (fBuf[1] >= nElm) ||
	     (fBuf[2] < 0) || (fBuf[2] >= nNod) ||
	     (fBuf[3] < 0) || (fBuf[3] >= nNod) ||
	     (fBuf[4] < 0) || (fBuf[4] >= nNod))
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	    break;
	  }
	  for(idE = 0; idE < 2; ++idE)
	  {
	    if(fBuf[idE] >= 0)
	    {
	      eBufP = eBuf + (5 *  fBuf[idE]);
	      if(eBufP[0] < 4)
	      {
		if(eBufP[0] == 0)
		{
		  eBufP[0] = 3;
		  eBufP[1] = fBuf[2];
		  eBufP[2] = fBuf[3];
		  eBufP[3] = fBuf[4];
		}
		else
		{
		  for(id0 = 2; id0 < 5; ++id0)    /* Input face field index. */
		  {
		    cnt = 0;
		    for(id1 = 1; id1 < 4; ++id1)   /* Elm buffer node index. */
		    {
		      if(fBuf[id0] == eBufP[id1])
		      {
		        ++cnt;
		      }
		    }
		    if(cnt == 0)
		    {
		      eBufP[0] = 4;
		      eBufP[4] = fBuf[id0];
		      break;
		    }
		  }
		}
	      }
	    }
	  }
          ++idF;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      eBufP = eBuf + (5 * idE);
      if(eBufP[0] != 4)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      for(idN = 0; idN < 4; ++idN)
      {
        nBuf[idN] = (WlzCMeshNod3D *)
	            AlcVectorItemGet(mesh->res.nod.vec, eBufP[1 + idN]);
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
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  AlcFree(eBuf);
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
static char	*WlzEffReadObjVMeshRec(FILE *fP, char *buf, int bufLen)
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
* 		GRUMP VMESH tetrahedral mesh file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjVMesh(FILE *fP, WlzObject *obj)
{
  int		cnt,
  		idE,
		idF,
  		idN,
		idE0,
		idE1,
  		nElm = 0,
  		nFce = 0,
		nBFce = 0,
		nNod = 0;
  int		*elmTbl = NULL,
  		*nodTbl = NULL;
  WlzCMesh3D	*mesh;
  WlzCMeshEdgU3D *edu;
  WlzCMeshElm3D	*elm;
  WlzCMeshFace	*fce;
  WlzCMeshNod3D	*nod;
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
    if(mesh->type != WLZ_CMESH_TET3D)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((elmTbl = (int *)
                  AlcMalloc(mesh->res.elm.maxEnt * sizeof(int))) == NULL) ||
       ((nodTbl = (int *)
                  AlcMalloc(mesh->res.nod.maxEnt * sizeof(int))) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Compute the number of faces and boundary faces while building an
   * element table to avoid deleted elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    nNod = mesh->res.nod.numEnt;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	++nElm;
	elmTbl[idE] = cnt++;
	for(idF = 0; idF < 4; ++idF)
	{
	  if((elm->face[idF].opp == NULL) ||
	     (elm->face[idF].opp == &(elm->face[idF])))
	  
	  {
	    ++nFce;
	    ++nBFce;
	  }
	  else if(elm->idx > elm->face[idF].opp->elm->idx)
	  {
	    ++nFce;
	  }
	}
        ++nElm;
      }
    }
  }
  /* Print number of elements, faces, boundary faces and nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "%d %d %d %d\n", nElm / 2, nFce, nBFce, nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Print node positions while building a node table to avoid deleted
   * nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    cnt = 0;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
        nodTbl[idN] = cnt++;
	if(fprintf(fP, "%lg %lg %lg\n",
	           nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  /* Print faces defined by element indices and then node indices. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        for(idF = 0; idF < 4; ++idF)
	{
	  idE0 = idE1 = -1;
	  fce = elm->face + idF;
	  if((fce->opp == NULL) || (fce->opp == fce))
	  {
	    idE1 = elmTbl[elm->idx];
	  }
	  else if(fce->opp->elm->idx > elm->idx)
	  {
	    idE0 = elmTbl[fce->opp->elm->idx];
	    idE1 = elmTbl[elm->idx];
	  }
	  if(idE1 >= 0)
	  {
	    edu = fce->edu;
            if(fprintf(fP, "%d %d %d %d %d\n",
	               idE0, idE1,
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
  /* Print boundary faces again defined by the element index and then node
   * indices. */
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        for(idF = 0; idF < 4; ++idF)
	{
	  fce = elm->face + idF;
	  if((fce->opp == NULL) || (fce->opp == fce))
	  {
	    edu = fce->edu;
            if(fprintf(fP, "%d %d %d %d %d\n",
	               elmTbl[elm->idx], 1,
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
  AlcFree(elmTbl);
  AlcFree(nodTbl);
  return(errNum);
}
