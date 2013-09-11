#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFStl_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFStl.c
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
* 		the the stereolithography stl surface file sormat.
* \ingroup	WlzExtFF
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static char			*WlzEffReadObjStlRec(
				  FILE *fP, 
				  char *buf,
				  int bufLen);
static WlzErrorNum		WlzEffWriteObjCM2D5Stl(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum		WlzEffWriteObjCtrStl(
				  FILE *fP,
				  WlzObject *obj);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		stereolithography stl surface file format. This format
* 		has surface definitions defined by a name and a collection
* 		of facet normals and loops. Each facet loop has the
* 		cooordinates of it's vertices.
* 		Only triangulated surface models can be read and all solids
* 		will be treated as one.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjStl(FILE *fP, WlzErrorNum *dstErr)
{
  int		vCnt = 0,
  		inSolid = 0,
  		inFacet = 0,
		inLoop = 0;
  char		*sav,
  		*str,
  		*tok;
  WlzDVertex3	vBuf[3];
  WlzGMModel	*model = NULL;
  WlzObject	*obj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		cBuf[256];

  dom.core = NULL;
  val.core = NULL;
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    model = WlzGMModelNew(WLZ_GMMOD_3D, 0, 0, &errNum);
  }
  while((errNum == WLZ_ERR_NONE) &&
        ((str = WlzEffReadObjStlRec(fP, cBuf, 256)) != NULL))
  {
    if((tok = ALC_STRTOK_R(str, " \t", &sav)) != NULL)
    {
      if(strncmp(tok, "solid", 5) == 0)
      {
        if(inSolid == 0)
	{
	  inSolid = 1;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strncmp(tok, "facet", 5) == 0)
      {
        if((inSolid == 1) && (inFacet == 0))
	{
	  inFacet = 1;
	  /* Normal vector is ignored. */
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strncmp(tok, "outer", 5) == 0)
      {
        if(((tok = ALC_STRTOK_R(NULL, " \t", &sav)) == NULL) ||
	   (strncmp(tok, "loop", 4) != 0) ||
           (inSolid == 0) || (inFacet == 0) || (inLoop != 0))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  vCnt = 0;
	  inLoop = 1;
	}
      }
      else if(strncmp(tok, "vertex", 6) == 0)
      {
	char *pTok[3];

        if((vCnt < 3) &&
	   ((pTok[0] = ALC_STRTOK_R(NULL, " \t", &sav)) != NULL) &&
	   ((pTok[1] = ALC_STRTOK_R(NULL, " \t", &sav)) != NULL) &&
	   ((pTok[2] = ALC_STRTOK_R(NULL, " \t", &sav)) != NULL) &&
	   (sscanf(pTok[0], "%lg", &(vBuf[vCnt].vtX)) == 1) &&
	   (sscanf(pTok[1], "%lg", &(vBuf[vCnt].vtY)) == 1) &&
	   (sscanf(pTok[2], "%lg", &(vBuf[vCnt].vtZ)) == 1))
	{
	  ++vCnt;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strncmp(tok, "endloop", 7) == 0)
      {
        if(inLoop == 1)
	{
	  inLoop = 0;
	  if(vCnt == 3)
	  {
	    errNum = WlzGMModelConstructSimplex3D(model, vBuf);
	  }
	  else
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strncmp(tok, "endfacet", 8) == 0)
      {
        if(inFacet == 1)
	{
	  inFacet = 0;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
      else if(strncmp(tok, "endsolid", 8) == 0)
      {
        if(inSolid == 1)
	{
	  inSolid = 0;
	}
	else
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
      }
    }
  }
  /* Create the Woolz object. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((dom.ctr = WlzMakeContour(&errNum)) != NULL)
    {
      dom.ctr->model = WlzAssignGMModel(model, NULL);
      obj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
    }
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
* 		using the stereolithography stl file format, see
* 		WlzEffReadObjStl().
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjStl(FILE *fP, WlzObject *obj)
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
        errNum = WlzEffWriteObjCM2D5Stl(fP, obj);
        break;
      case WLZ_CONTOUR:
        errNum = WlzEffWriteObjCtrStl(fP, obj);
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
* 		stereolithography stl file format, see WlzEffReadObjStl().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCM2D5Stl(FILE *fP, WlzObject *obj)
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
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "solid ascii\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the elements. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		idE;
    WlzDVertex3	nrm;
    WlzCMeshElm2D5 *elm;
    WlzCMeshNod2D5 *nod[3];

    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
        WlzCMeshElmGetNodes2D5(elm, nod + 0, nod + 1, nod + 2);
        nrm = WlzGeomTriangleNormal(nod[0]->pos, nod[1]->pos, nod[2]->pos);
	if(fprintf(fP,
	           "  facet normal %g %g %g\n"
		   "    outer loop\n"
	           "      vertex %g %g %g\n"
	           "      vertex %g %g %g\n"
	           "      vertex %g %g %g\n"
		   "    endloop\n"
	           "  endfacet\n",
		   nrm.vtX, nrm.vtY, nrm.vtZ,
		   nod[0]->pos.vtX, nod[0]->pos.vtY, nod[0]->pos.vtZ,
		   nod[1]->pos.vtX, nod[1]->pos.vtY, nod[1]->pos.vtZ,
		   nod[2]->pos.vtX, nod[2]->pos.vtY, nod[2]->pos.vtZ) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "endsolid\n") <= 0)
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
* 		WLZ_CONTOUR) to the given file stream using the
* 		stereolithography stl file format, see WlzEffReadObjStl().
* \param	fP			Output file stream.
* \param	obj			Given woolz object (must not be NULL).
*/
static WlzErrorNum WlzEffWriteObjCtrStl(FILE *fP, WlzObject *obj)
{
  int		nFce;
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
  /* Check there are verticies and faces. */
  if(errNum == WLZ_ERR_NONE)
  {
    if((model->res.vertex.numElm < 3) || ((nFce = model->res.face.numElm) < 1))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "solid ascii\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		idF;
    WlzGMFace	*fce;
    AlcVector	*vec;
    WlzDVertex3 nrm;
    WlzDVertex3	vBuf[3];

    vec = model->res.face.vec;
    for(idF = 0; idF < model->res.face.numIdx; ++idF)
    {
      fce = (WlzGMFace *)AlcVectorItemGet(vec, idF);
      if(fce->idx >= 0)
      {
        errNum = WlzGMFaceGetG3D(fce, vBuf + 0, vBuf + 1, vBuf + 2);
	if(errNum == WLZ_ERR_NONE)
	{
	  nrm = WlzGeomTriangleNormal(vBuf[0], vBuf[1], vBuf[2]);
	  if(fprintf(fP,
		     "  facet normal %g %g %g\n"
		     "    outer loop\n"
		     "      vertex %g %g %g\n"
		     "      vertex %g %g %g\n"
		     "      vertex %g %g %g\n"
		     "    endloop\n"
		     "  endfacet\n",
		     nrm.vtX, nrm.vtY, nrm.vtZ,
		     vBuf[0].vtX, vBuf[0].vtY, vBuf[0].vtZ,
		     vBuf[1].vtX, vBuf[1].vtY, vBuf[1].vtZ,
		     vBuf[2].vtX, vBuf[2].vtY, vBuf[2].vtZ) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
      }
    }
  }
  WlzGMModelResIdxFree(resIdxTb);
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "endsolid\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
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
static char	*WlzEffReadObjStlRec(FILE *fP, char *buf, int bufLen)
{
  char		*str;

  str = fgets(buf, bufLen, fP);
  buf[bufLen - 1] = '\0';
  return(str);
}
