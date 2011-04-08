#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFSMesh_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFSMesh.c
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
* 		the the GRUMMP surface mesh (smesh) file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static int			WlzEffReadObjSMeshRec(
				  FILE *fP, 
				  char *buf,
				  int bufLen,
				  int nFld);
static WlzErrorNum		WlzEffWriteObjCM2SMesh(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteObjCtrSMesh(
				  FILE *fP,
				  WlzObject *obj);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		GRUMMP smesh file format. This format is:
* 		\verbatim
    # comment for rest of record
    <n nodes> <dim (always 3)> <n attrib> <n boundary markers>
    <node idx> <x> <y> <z> [<attributes>] [<boundary mark>]
    ...
    <n facets> <n boundary markers>
    <n nodes> <node idx> <node idx> <node idx> <facet idx> [<boundary mark>]
    ...
    <n holes>
    <hole idx> <x> <y> <z>
    ...
    <n region>
    <region idx> <x> <y> <z> <region n> <region attributes>
    ...
	 	\endverbatim
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjSMesh(FILE *fP, WlzErrorNum *dstErr)
{
  int		dim = 0,
		nAtr = 0,
		nBMk = 0,
  		nFct = 0,
		nVtx = 0;
  WlzGMModel	*model = NULL;
  WlzObject	*obj = NULL;
  WlzDVertex3	*vBuf = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		cBuf[256];

  dom.core = NULL;
  val.core = NULL;
  /* Read number of nodes and elements. */
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if((WlzEffReadObjSMeshRec(fP, cBuf, 256, 4) != 4) ||
       (sscanf(cBuf, "%d %d %d %d", &nVtx, &dim, &nAtr, &nBMk) != 4) ||
       (nVtx < 3) || (dim != 3) || (nAtr < 0) || (nBMk < 0))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Create contour. */
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr = WlzMakeContour(&errNum);
  }
  /* Create a new geometric model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		vHTSz;

    vHTSz = nVtx / 16;
    if(vHTSz < 1024)
    {
      vHTSz = 1024;
    }
    model = WlzGMModelNew(WLZ_GMMOD_3D, 0, vHTSz, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create a vertex buffer. */
    if((vBuf = AlcMalloc(sizeof(WlzDVertex3) * nVtx)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Read the vertex positions into the buffer. */
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idN,
    		idX,
    		nFld;

    nFld = 4 + nAtr + nBMk;
    for(idN = 0; idN < nVtx; ++idN)
    {
      WlzDVertex3 p;
      if((WlzEffReadObjSMeshRec(fP, cBuf, 256, nFld) != nFld) ||
	 (sscanf(cBuf, "%d %lg %lg %lg",
	         &idX, &(p.vtX), &(p.vtY), &(p.vtZ)) != 4) ||
	 (idX < 1) || (idX > nVtx))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      *(vBuf + idX - 1) = p;
    }
  }
  /* Read number of facets and boundary markers. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((WlzEffReadObjSMeshRec(fP, cBuf, 256, 2) != 2) ||
       (sscanf(cBuf, "%d %d", &nFct, &nBMk) != 2) ||
       (nFct < 1) || (nBMk < 0))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  /* Read the facets indicies and build the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idE,
    		nFld;
    int		idx[3];
    WlzDVertex3	pos[3];

    for(idE = 0; idE < nFct; ++idE)
    {
      int	nFVx;

      if((WlzEffReadObjSMeshRec(fP, cBuf, 256, nFld) < 4) ||
	 (sscanf(cBuf, "%d %d %d %d",
	         &nFVx, idx + 0, idx + 1, idx + 2) != 4) ||
	 (nFVx != 3) ||
	 (idx[0] < 1)     || (idx[1] < 1)     || (idx[2] < 1) ||
	 (idx[0] > nVtx) || (idx[1] > nVtx) || (idx[2] > nVtx))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
	pos[0] = *(vBuf + idx[0] - 1);
	pos[1] = *(vBuf + idx[1] - 1);
	pos[2] = *(vBuf + idx[2] - 1);
	errNum = WlzGMModelConstructSimplex3D(model, pos);
      }
      if(errNum != WLZ_ERR_NONE)
      {
        break;
      }
    }
  }
  AlcFree(vBuf);
  /* Compute maximum edge length and then create the Woolz object. */
  if(errNum == WLZ_ERR_NONE)
  {
    dom.ctr->model = WlzAssignGMModel(model, NULL);
    obj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(obj != NULL)
    {
      (void )WlzFreeObj(obj);
    }
    else if(dom.ctr != NULL)
    {
      (void )WlzFreeContour(dom.ctr);
    }
    else
    {
      (void )WlzGMModelFree(model);
    }
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given file stream
* 		using the GRUMMP smesh file format, see
* 		WlzEffReadObjSMesh().
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjSMesh(FILE *fP, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  {
    switch(obj->type)
    {
      case WLZ_CMESH_2D5:
        errNum = WlzEffWriteObjCM2SMesh(fP, obj);
        break;
      case WLZ_CONTOUR:
        errNum = WlzEffWriteObjCtrSMesh(fP, obj);
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object (which is known to be a
* 		WLZ_CMESH_2D5) object to the given file stream using the
* 		GRUMMP smesh file format, see WlzEffReadObjSMesh().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCM2SMesh(FILE *fP, WlzObject *obj)
{
  int		nElm = 0,
		nNod = 0;
  int		*nodTbl = NULL;
  WlzCMesh2D5	*mesh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    mesh = obj->domain.cm2d5;
    if(mesh->type != WLZ_CMESH_2D5)
    {
      errNum = WLZ_ERR_DOMAIN_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((nodTbl = (int *)AlcMalloc(mesh->res.nod.maxEnt * sizeof(int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Output the number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    nNod = mesh->res.nod.numEnt;
    nElm = mesh->res.elm.numEnt;
    if(fprintf(fP,
               "# Exported to smesh format via Woolz.\n"
	       "%d 3 0 0\n", nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Create a valid node LUT while outputing the node positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idN,
    		nCnt;
    AlcVector	*vec;
    WlzCMeshNod2D5 *nod;

    nCnt = 0;
    vec = mesh->res.nod.vec;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(vec, idN);
      if(nod->idx >= 0)
      {
        nodTbl[idN] = nCnt++;
	if(fprintf(fP, "%d  %g %g %g\n",
	           nCnt, nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ) <=0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  /* Output the number of elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
	       "%d 1\n", nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE,
    		eCnt;
    WlzCMeshElm2D5 *elm;
    WlzCMeshNod2D5 *nod[3];

    eCnt = 0;
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	++eCnt;
        WlzCMeshElmGetNodes2D5(elm, nod + 0, nod + 1, nod + 2);
	if(fprintf(fP, "3 %d %d %d %d\n",
		   nodTbl[nod[0]->idx] + 1, 
		   nodTbl[nod[1]->idx] + 1, 
		   nodTbl[nod[2]->idx] + 1,
		   eCnt) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  AlcFree(nodTbl);
  /* Output the number of holes and regions (both zero). */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "0\n0\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}


/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object (which is known to be a
* 		WLZ_CONTOUR) object to the given file stream using the
* 		GRUMMP smesh file format, see WlzEffReadObjSMesh().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCtrSMesh(FILE *fP, WlzObject *obj)
{
  int		nFce,
  		nVtx;
  WlzGMModel	*model;
  WlzGMResIdxTb *resIdxTb = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((obj->domain.core == NULL) || ((model = obj->domain.ctr->model) == NULL))
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(model->type)
    {
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Index the vertices. */
    resIdxTb = WlzGMModelResIdx(model, WLZ_GMELMFLG_VERTEX, &errNum);
  }
  /* Check there are verticies and faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(((nVtx = resIdxTb->vertex.idxCnt) != model->res.vertex.numElm) ||
       (model->res.vertex.numElm < 3) || ((nFce = model->res.face.numElm) < 1))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  /* Output the number of vertices. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
               "# Exported to smesh format via Woolz.\n"
	       "%d 3 0 0\n", nVtx) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the vertex positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idV,
    		vCnt;
    AlcVector	*vec;
    WlzDVertex3	pos;

    vCnt = 0;
    vec = model->res.vertex.vec;
    for(idV = 0; idV < model->res.vertex.numIdx; ++idV)
    {
      WlzGMVertex *vtx;

      vtx = (WlzGMVertex *)AlcVectorItemGet(vec, idV);
      if(vtx->idx >= 0)
      {
	++vCnt;
	(void )WlzGMVertexGetG3D(vtx, &pos);
	if(fprintf(fP, "%d  %g %g %g\n",
	           vCnt, pos.vtX, pos.vtY, pos.vtZ) <=0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  /* Output the number of faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
	       "%d 1\n", nFce) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the vertex indices for the facets. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idF,
    		fCnt;
    int		*lut;
    WlzGMFace	*fce;
    AlcVector	*vec;
    int		idx[3];

    fCnt = 0;
    vec = model->res.face.vec;
    lut = resIdxTb->vertex.idxLut;
    for(idF = 0; idF < model->res.face.numIdx; ++idF)
    {
      fce = (WlzGMFace *)AlcVectorItemGet(vec, idF);
      if(fce->idx >= 0)
      {
        WlzGMEdgeT	*tET;

	++fCnt;
	tET = fce->loopT->edgeT;
	idx[0] = *(lut + tET->vertexT->diskT->vertex->idx);
	idx[1] = *(lut + tET->next->vertexT->diskT->vertex->idx);
	idx[2] = *(lut + tET->prev->vertexT->diskT->vertex->idx);

	if(fprintf(fP, "3 %d %d %d %d\n",
		   idx[0] + 1, idx[1] + 1, idx[2] + 1, fCnt) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  WlzGMModelResIdxFree(resIdxTb);
  /* Output the number of holes and regions (both zero). */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "0\n0\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Number of fields found.
* \ingroup	WlzExtFF
* \brief	Attempts to get the requeted number of fields from a single
* 		record of the file.
*               
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
* \param	nf			Number of fields required.
*/
static int	WlzEffReadObjSMeshRec(FILE *fP, char *buf, int bufLen,
				      int nf)
{
  int		f = 0;

  while(!f)
  {
    int		ws;
    char	*s0,
  		*s1,
		*s2; 

    /* Read record. */
    if(fgets(buf, bufLen, fP) == NULL)
    {
      break;
    }
    buf[bufLen - 1] = '\0';
    /* Squeeze out possible leading space and all multiple spaces from the
     * record. Skip remainder if comment character encountered. */
    ws = 1;
    s0 = s1 = s2 = buf;
    while(*s2 && (*s2 != '#') && (*s2 != '\n'))
    {
      int	sp;

      sp = isspace(*s2);
      if(sp)
      {
        if(!ws)
	{
	  *s0++ = ' ';
	}
      }
      else
      {
	if(ws || (s0 == buf))
	{
	  ++f;
	}
        *s0++ = *s2;;
      }
      ws = sp;
      s1 = s2;
      ++s2;
    }
    *s0 = '\0';
  }
  return(f);
}
