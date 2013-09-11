#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFObj.c
* \author       Bill Hill
* \date         September 2009
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
* 		the the OBJ Wavefront geometry file format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjObjRec(
				  FILE *fP, 
				  char *buf,
				  int bufLen);
static WlzErrorNum		WlzEffWriteObjCM2Obj(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteObjCtrObj(
				  FILE *fP,
				  WlzObject *obj);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		Wavefront OBJ file format. This format consists of the
* 		following possible fields (one per line):
* 		Vertex data
* 		   -v          vertex geometry
* 		   -vt         ignored
* 		   -vn         vertex normal
* 		   -vp         ignored
* 		   -cstype     ignored
* 		   -deg        ignored
* 		   -bmat       ignored
* 		   -step       ignored
*               Elements
*                  -p          ignored
*                  -l          ignored
*                  -f	       face 
*                  -curv       ignored
*                  -curv2      ignored
*                  -surf       ignored
*		Freeform elements
*		   -param      ignored
*		   -trim       ignored
*		   -hole       ignored
*		   -scrv       ignored
*		   -sp         ignored
*		   -end        ignored
*		Freeform surface connectivity
*		   -con        ignored
*		Grouping
*		   -g          ignored
*		   -o          ignored
*		   -s          ignored
*		   -mg         ignored
*		Rendering attributes
*		   -bevel      ignored
*		   -c_interp   ignored
*		   -d_interp   ignored
*		   -lod        ignored
*		   -usemtl     ignored
*		   -mtllib     ignored
*		   -shadow_obj ignored
*		   -trace_obj  ignored
*		   -ctech      ignored
*		   -stech      ignored
* 		Comments
* 		   -#          ignored
* 		Only triangulated surface models can be read. All groupings
* 		will be amalgamated into one.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjObj(FILE *fP, WlzErrorNum *dstErr)
{
  int		nF = 0,
		nN = 0,
		nV = 0;
  char		*sav,
  		*str,
  		*tok;
  int		*iP;
  WlzDVertex3	*dP;
  WlzGMModel	*model = NULL;
  WlzObject	*obj = NULL;
  AlcVector	*fVec = NULL,
  		*nVec = NULL,
  		*vVec = NULL;
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
  else if(((vVec = AlcVectorNew(1, sizeof(WlzDVertex3), 4096, NULL)) == NULL) ||
          ((nVec = AlcVectorNew(1, sizeof(WlzDVertex3), 4096, NULL)) == NULL) ||
          ((fVec = AlcVectorNew(1, 3 * sizeof(int), 4096, NULL)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  while((errNum == WLZ_ERR_NONE) &&
        ((str = WlzEffReadObjObjRec(fP, cBuf, 256)) != NULL))
  {
    if((tok = ALC_STRTOK_R(str, " \t", &sav)) != NULL)
    {
      if(strcmp(tok, "v") == 0)
      {
	if((dP = (WlzDVertex3 *)AlcVectorExtendAndGet(vVec, nV++)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else if(((str = ALC_STRTOK_R(NULL, "\n", &sav)) == NULL) ||
		(sscanf(str, "%lg %lg %lg",
		        &(dP->vtX), &(dP->vtY), &(dP->vtZ)) != 3))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strcmp(tok, "vn") == 0)
      {
	if((dP = (WlzDVertex3 *)AlcVectorExtendAndGet(nVec, nN++)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else if(((str = ALC_STRTOK_R(NULL, "\n", &sav)) == NULL) ||
		(sscanf(str, "%lg %lg %lg",
		        &(dP->vtX), &(dP->vtY), &(dP->vtZ)) != 3))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strcmp(tok, "f") == 0)
      {
	if((iP = (int *)AlcVectorExtendAndGet(fVec, nF++)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else if(((str = ALC_STRTOK_R(NULL, "\n", &sav)) == NULL) ||
		(sscanf(str, "%d %d %d", iP + 0, iP + 1, iP + 2) != 3))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  --iP[0];
	  --iP[1];
	  --iP[2];
	}
      }
      /* All other tokens are ignored. */
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((nV < 3) || (nF < 1) || ((nN > 0) && (nN != nV)))
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

    vHTSz = nV / 8;
    if(vHTSz < 1024)
    {
      vHTSz = 1024;
    }
    if(nN > 0)
    {
      model = WlzGMModelNew(WLZ_GMMOD_3N, 0, vHTSz, &errNum);
    }
    else
    {
      model = WlzGMModelNew(WLZ_GMMOD_3D, 0, vHTSz, &errNum);
    }
  }
  /* Build the model. */
  if(errNum == WLZ_ERR_NONE)
  {
    int 	idF;
    int	 	*idx;
    WlzDVertex3 pos[3],
		nrm[3];

    for(idF = 0; idF < nF; ++idF)
    {
      idx = (int *)AlcVectorItemGet(fVec, idF);
      if((idx[0] < 0) || (idx[1] < 0) || (idx[2] < 0) ||
	 (idx[0] >= nV) || (idx[1] >= nV) || (idx[2] >= nV))
      {
	errNum = WLZ_ERR_DOMAIN_DATA;
      }
      else
      {
	pos[0] = *(WlzDVertex3 *)AlcVectorItemGet(vVec, idx[0]);
	pos[1] = *(WlzDVertex3 *)AlcVectorItemGet(vVec, idx[1]);
	pos[2] = *(WlzDVertex3 *)AlcVectorItemGet(vVec, idx[2]);
	if(nN > 0)
	{
	  nrm[0] = *(WlzDVertex3 *)AlcVectorItemGet(nVec, idx[0]);
	  nrm[1] = *(WlzDVertex3 *)AlcVectorItemGet(nVec, idx[1]);
	  nrm[2] = *(WlzDVertex3 *)AlcVectorItemGet(nVec, idx[2]);
	  errNum = WlzGMModelConstructSimplex3N(model, pos, nrm);
	}
	else
	{
	  errNum = WlzGMModelConstructSimplex3D(model, pos);
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
	break;
      }
    }
  }
  /* Free termporary vectors. */
  (void )AlcVectorFree(vVec);
  (void )AlcVectorFree(nVec);
  (void )AlcVectorFree(fVec);
  /* Create the Woolz object. */
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
* 		using the Wavefront OBJ file format, see WlzEffReadObjObj().
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjObj(FILE *fP, WlzObject *obj)
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
        errNum = WlzEffWriteObjCM2Obj(fP, obj);
        break;
      case WLZ_CONTOUR:
        errNum = WlzEffWriteObjCtrObj(fP, obj);
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
* 		wavefront obj file format, see WlzEffReadObjObj().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCM2Obj(FILE *fP, WlzObject *obj)
{
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
  /* Output comment header. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
               "# wavefront obj file written by WlzEffWriteObjCM2Obj()\n") <= 0)
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
        nodTbl[idN] = ++nCnt;     /* ++nCnt because indices start at 1 not 0. */
	if(fprintf(fP, "v %lg %lg %lg\n",
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
	if(fprintf(fP, "v %d %d %d\n",
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
* \brief	Writes the given Woolz object (which is known to be a
* 		WLZ_CONTOUR) to the given file stream using the Wavefront
* 		obj file format, see WlzEffReadObjObj().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCtrObj(FILE *fP, WlzObject *obj)
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
  /* Output comment header. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
               "# wavefront obj file written by WlzEffWriteObjCtrObj()\n") <= 0)
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
	if(fprintf(fP, "v %lg %lg %lg\n", pos.vtX, pos.vtY, pos.vtZ) <=0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
    }
  }
  /* Output vertex normals if the model has them. */
  if((errNum == WLZ_ERR_NONE) && (model->type == WLZ_GMMOD_3N))
  {
    int		idV;
    AlcVector	*vec;
    WlzDVertex3	nrm;

    vec = model->res.vertex.vec;
    for(idV = 0; idV < model->res.vertex.numIdx; ++idV)
    {
      WlzGMVertex *vtx;

      vtx = (WlzGMVertex *)AlcVectorItemGet(vec, idV);
      if(vtx->idx >= 0)
      {
	(void )WlzGMVertexGetG3N(vtx, NULL, &nrm);
	if(fprintf(fP, "vn %lg %lg %lg\n", nrm.vtX, nrm.vtY, nrm.vtZ) <=0)
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
	idx[0] = *(lut + tET->vertexT->diskT->vertex->idx) + 1;
	idx[1] = *(lut + tET->next->vertexT->diskT->vertex->idx) + 1;
	idx[2] = *(lut + tET->prev->vertexT->diskT->vertex->idx) + 1;
	if(fprintf(fP, "f %d %d %d\n", idx[0], idx[1], idx[2]) <= 0)
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
static char	*WlzEffReadObjObjRec(FILE *fP, char *buf, int bufLen)
{
  char		*str;

  str = fgets(buf, bufLen, fP);
  buf[bufLen - 1] = '\0';
  return(str);
}
