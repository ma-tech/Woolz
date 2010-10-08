#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFPly2_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFPly2.c
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
* 		the the Riken PLY2 file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjPly2Rec(
				  FILE *fP, 
				  char *buf,
				  int bufLen,
				  int nFld);
static WlzErrorNum		WlzEffWriteObjCM2Ply2(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteObjCtrPly2(
				  FILE *fP,
				  WlzObject *obj);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		Riken PLY2 triangular mesh file format. This format
* 		simply consists of:
* 		\verbatim
	        <n vertices (int)>
	        <n triangles (int)>
	        <vertex (double double double)>
	        ...
	        <triangle (int int int)>
	        ...
	 	\endverbatim
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjPly2(FILE *fP, WlzErrorNum *dstErr)
{
  int		nFce = 0,
		nVtx = 0;
  char		*str;
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
    if(((str = WlzEffReadObjPly2Rec(fP, cBuf, 256, 1)) == NULL) ||
       (sscanf(str, "%d", &nVtx) != 1))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(((str = WlzEffReadObjPly2Rec(fP, cBuf, 256, 1)) == NULL) ||
       (sscanf(str, "%d", &nFce) != 1))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((nVtx < 3) || (nFce < 1))
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
    int 	idN;
    WlzDVertex3 *pos;

    for(idN = 0; idN < nVtx; ++idN)
    {
      pos = vBuf + idN;
      if(((str = WlzEffReadObjPly2Rec(fP, cBuf, 256, 3)) == NULL) ||
	 (sscanf(str, "%lg %lg %lg", &(pos->vtX), &(pos->vtY), &(pos->vtZ)) != 3))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
    }
  }
  /* Read the vertex indicies and build the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idE;
    int		idx[3];
    WlzDVertex3	pos[3];

    for(idE = 0; idE < nFce; ++idE)
    {
      int	ck;

      if(((str = WlzEffReadObjPly2Rec(fP, cBuf, 256, 4)) == NULL) ||
	 (sscanf(str, "%d %d %d %d", &ck, idx + 0, idx + 1, idx + 2) != 4) ||
	 (ck != 3) || (idx[0] < 0) || (idx[1] < 0) || (idx[2] < 0) ||
	 (idx[0] >= nVtx) || (idx[1] >= nVtx) || (idx[2] >= nVtx))
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
	pos[0] = *(vBuf + idx[0]);
	pos[1] = *(vBuf + idx[1]);
	pos[2] = *(vBuf + idx[2]);
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
* 		using the Riken PLY2 triangular mesh file format, see
* 		WlzEffReadObjPly2().
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjPly2(FILE *fP, WlzObject *obj)
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
        errNum = WlzEffWriteObjCM2Ply2(fP, obj);
        break;
      case WLZ_CONTOUR:
        errNum = WlzEffWriteObjCtrPly2(fP, obj);
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
* \brief	Writes the given Woolz object (which is known to be a WLZ_CMESH_2D5)
* 		object to the given file stream using the Riken PLY2 triangular mesh
* 		file format, see WlzEffReadObjPly2().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCM2Ply2(FILE *fP, WlzObject *obj)
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
  /* Output the number of nodes and elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    nNod = mesh->res.nod.numEnt;
    nElm = mesh->res.elm.numEnt;
    if(fprintf(fP, "%d\n%d\n", nNod, nElm) <= 0)
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
	if(fprintf(fP, "%lg\n%lg\n%lg\n",
	           nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ) <=0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  /* Output the elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;
    WlzCMeshElm2D5 *elm;
    WlzCMeshNod2D5 *nod[3];

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        WlzCMeshElmGetNodes2D5(elm, nod + 0, nod + 1, nod + 2);
	if(fprintf(fP, "3\n%d\n%d\n%d\n",
		   nodTbl[nod[0]->idx], 
		   nodTbl[nod[1]->idx], 
		   nodTbl[nod[2]->idx]) <= 0)
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


/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object (which is known to be a WLZ_CONTOUR)
* 		object to the given file stream using the Riken PLY2 triangular mesh
* 		file format, see WlzEffReadObjPly2().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCtrPly2(FILE *fP, WlzObject *obj)
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
  /* Output the number of vertices and faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "%d\n%d\n", nVtx, nFce) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the vertex positions. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idV;
    AlcVector	*vec;
    WlzDVertex3	pos;

    vec = model->res.vertex.vec;
    for(idV = 0; idV < model->res.vertex.numIdx; ++idV)
    {
      WlzGMVertex *vtx;

      vtx = (WlzGMVertex *)AlcVectorItemGet(vec, idV);
      if(vtx->idx >= 0)
      {
	(void )WlzGMVertexGetG3D(vtx, &pos);
	if(fprintf(fP, "%lg\n%lg\n%lg\n", pos.vtX, pos.vtY, pos.vtZ) <=0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  /* Output the vertex indices for the faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idF;
    int		*lut;
    WlzGMFace	*fce;
    AlcVector	*vec;
    int		idx[3];

    vec = model->res.face.vec;
    lut = resIdxTb->vertex.idxLut;
    for(idF = 0; idF < model->res.face.numIdx; ++idF)
    {
      fce = (WlzGMFace *)AlcVectorItemGet(vec, idF);
      if(fce->idx >= 0)
      {
        WlzGMEdgeT	*tET;

	tET = fce->loopT->edgeT;
	idx[0] = *(lut + tET->vertexT->diskT->vertex->idx);
	idx[1] = *(lut + tET->next->vertexT->diskT->vertex->idx);
	idx[2] = *(lut + tET->prev->vertexT->diskT->vertex->idx);
	if(fprintf(fP, "3\n%d\n%d\n%d\n", idx[0], idx[1], idx[2]) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  WlzGMModelResIdxFree(resIdxTb);
  return(errNum);
}

/*!
* \return	String with requested number of space seperated fields.
* \ingroup	WlzExtFF
* \brief	Gets the requeted number of fields from the file.
*               
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
* \param	nFld			Number of fields required.
*/
static char	*WlzEffReadObjPly2Rec(FILE *fP, char *buf, int bufLen,
				      int nFld)
{
  int		c0,
  		c1,
		f;
  char		*s,
  		*str = NULL;
  const char	space = ' ';

  f = 0;
  s = buf;
  c0 = space;
  while(1)
  {
    if(s - buf > bufLen - 2)
    {
      break;
    }
    c1 = fgetc(fP);
    if(c1 == EOF)
    {
      break;
    }
    if(isspace(c1))
    {
      if(c0 != space)
      {
	if(++f >= nFld)
	{
	  str = buf;
	  break;
	}
	else
	{
	  c1 = space;
	  *s++ = space;
	}
      }
    }
    else
    {
      *s++ = c1;
    }
    c0 = c1;
  }
  *s = '\0';
  return(str);
}
