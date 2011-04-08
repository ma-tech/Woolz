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
* 		the the tetrahedral mesh file format used by Pascal Frey's
* 		medit. See the INRIA technical report 0253 for a detailed
* 		explaination of the fromat. This file format is also used by
* 		the tetrhedral mesh generator tetgen.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static int			WlzEffReadObjMeshRec(
				  FILE *fP,
				  char *buf,
				  int bufLen,
				  int nFld);
static WlzErrorNum		WlzEffReadObjMeshSkip(
				  FILE *fP,
				  char *buf,
				  int bufLen,
				  int nFld);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using Pascal Frey's
* 		tetrahedral mesh file format.
* 		The file must start with the records:
* 		\verbatim
  MeshVersionFormatted 1
  Dimension 3
	 	\endverbatim
*		Following this a number of keywords followed by a cardinality
*		and then that number of further data fields follow:
*		\verbatim
  Vertices                       np     xi yi zi ri    {i in 1-np}
  Edges                          ne     e1i e2i ri     {i in 1-ne}
  Triangles                      nt     vij ri         {i in 1-nt},{j in 1-3}
  Quadrilaterals                 nq     vij ri         {i in 1-nq},{j in 1-4}
  Tetrahedra                     ntet   vij ri         {i in 1-ntet},{j in 1-4}
  Hexaedra                       nh     vij ri         {i in 1-nh},{j in 1-8}
  Corners                        nc     vi             {i in 1-nc}
  RequiredVertices               nrv    vi             {i in 1-nrv}
  Ridges                         nr     ei             {i in 1-nr}
  RequiredEdges                  nre    vi             {i in 1-nre}
  Normals                        nn     xi yi zi       {i in 1-nn}
  Tangents                       nnt    xi yi zi       {i in 1-nnt}
  NormalAtVertices               nv     vi ni          {i in 1-nv}
  NormalAtTriangleVertices       ntv    ti vj ni       {i in 1-ntv}
  NormalAtQuadrilateralVertices  nqv    qi vj ni       {i in 1-nqv}
  TangentAtEdges                 te     ei vj ni       {i in 1-te}
	 	\endverbatim
* 		At any point outside of a data field a comment character ('#')
* 		indicates that the rest of the line is to be ignored.
* 		The file terminates with the keyword: End.
*
* 		Only tetrahedral meshes can be read by this function and
* 		all fields except for the start, Vertices and Tetrahedra
* 		fields are ignored.
*
* 		The resulting object will be of type WLZ_CMESH_3D.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjMesh(FILE *fP, WlzErrorNum *dstErr)
{

  int		tI0,
  		tI1,
		nT = 0,
		nV = 0;
  WlzDVertex3   *vBuf = NULL;
  WlzCMesh3D	*mesh = NULL;
  WlzObject	*obj = NULL;
  WlzDomain	dom;
  WlzValues	val;
  char		inBuf[256],
  		sBuf0[252],
		sBuf1[256];
  const int	bufLen = 256;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  dom.core = NULL;
  val.core = NULL;
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((WlzEffReadObjMeshRec(fP, inBuf, bufLen, 4) != 4) ||
       (sscanf(inBuf, "%s %d %s %d",
               sBuf0, &tI0, sBuf1, &tI1) != 4) ||
       strcmp(sBuf0, "MeshVersionFormatted") || (tI0 != 1) ||
       strcmp(sBuf1, "Dimension") || (tI1 != 3))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    mesh = WlzCMeshNew3D(&errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		allKeywordsDone = 0;

    /* Read keywords and then thier associated cardinality followed by their
       data. */
    do
    {
      if(WlzEffReadObjMeshRec(fP, inBuf, bufLen, 1) != 1)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(strcmp(inBuf, "Vertices") == 0)
	{
	  /* Read vertices into buffer. */
	  if((WlzEffReadObjMeshRec(fP, inBuf, bufLen, 1) != 1) ||
	     (sscanf(inBuf, "%d", &nV) != 1) || (nV < 4))
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  else if((vBuf = (WlzDVertex3 *)
	                  AlcMalloc(sizeof(WlzDVertex3) * nV)) == NULL)
	  {
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  else
	  {
	    int		idV;

	    for(idV = 0; idV < nV; ++idV)
	    {
	      double	dummy;
	      WlzDVertex3 p;

	      if((WlzEffReadObjMeshRec(fP, inBuf, bufLen, 4) != 4) ||
	        (sscanf(inBuf, "%lg %lg %lg %lg",
		        &(p.vtX), &(p.vtY), &(p.vtZ), &dummy) != 4))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
		break;
	      }
	      vBuf[idV] = p;
	    }
	  }
	  /* Add vertices to the mesh. */
	  if(errNum == WLZ_ERR_NONE)
	  {
	    mesh->bBox = WlzBoundingBoxVtx3D(nV, vBuf, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    errNum = WlzCMeshReassignGridCells3D(mesh, nV);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(AlcVectorExtendAndGet(mesh->res.nod.vec, nV) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int         idV;

	    for(idV = 0; idV < nV; ++idV)
	    {
	      WlzCMeshNod3D *nod;

	      nod = WlzCMeshNewNod3D(mesh, vBuf[idV], &errNum);
	      nod->flags = 0;
	    }
	  }
	}
        else if(strcmp(inBuf, "Edges") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "Triangles") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 4);
	}
        else if(strcmp(inBuf, "Quadrilaterals") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 5);
	}
        else if(strcmp(inBuf, "Tetrahedra") == 0)
	{
	  if((nV < 4) ||
	     (WlzEffReadObjMeshRec(fP, inBuf, bufLen, 1) != 1) ||
	     (sscanf(inBuf, "%d", &nT) != 1) || (nT < 1))
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  else
	  {
	    if(AlcVectorExtendAndGet(mesh->res.elm.vec, nT) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    int		idT;
	    WlzCMeshNod3D *nod[4];

	    for(idT = 0; idT < nT; ++idT)
	    {
	      int	dummy;
	      int	idV[4];

	      if((WlzEffReadObjMeshRec(fP, inBuf, bufLen, 5) != 5) ||
	        (sscanf(inBuf, "%d %d %d %d %d",
		        idV + 0, idV + 1, idV + 2, idV + 3, &dummy) != 5) ||
	        (idV[0] < 1) || (idV[0] > nV) ||
	        (idV[1] < 1) || (idV[1] > nV) ||
	        (idV[2] < 1) || (idV[2] > nV) ||
	        (idV[3] < 1) || (idV[3] > nV))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
	        int 	idN;
		WlzCMeshElm3D *elm;

	        for(idN = 0; idN < 4; ++idN)
		{
		  nod[idN] = (WlzCMeshNod3D *)
		             AlcVectorItemGet(mesh->res.nod.vec, idV[idN] - 1);
		}
		elm = WlzCMeshNewElm3D(mesh,
		                       nod[0], nod[1], nod[2], nod[3], 1,
				       &errNum);
	        if(elm != NULL)
		{
		  elm->flags = 0;
		}
	      }
	      if(errNum != WLZ_ERR_NONE)
	      {
	        break;
	      }
	    }
	  }
	}
        else if(strcmp(inBuf, "Hexaedra") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 9);
	}
        else if(strcmp(inBuf, "Corners") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 1);
	}
        else if(strcmp(inBuf, "RequiredVertices") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 1);
	}
        else if(strcmp(inBuf, "Ridges") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 1);
	}
        else if(strcmp(inBuf, "RequiredEdges") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 1);
	}
        else if(strcmp(inBuf, "Normals") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "Tangents") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "NormalAtVertices") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 2);
	}
        else if(strcmp(inBuf, "NormalAtTriangleVertices") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "NormalAtQuadrilateralVertices") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "TangentAtEdges") == 0)
	{
	  errNum = WlzEffReadObjMeshSkip(fP, inBuf, bufLen, 3);
	}
        else if(strcmp(inBuf, "End") == 0)
	{
	  allKeywordsDone = 1;
	}
      }
    } while((errNum == WLZ_ERR_NONE) && (allKeywordsDone == 0));
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
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given stream using
* 		Pascal Frey's tetrahedral mesh file format. See
* 		WlzEffReadObjMesh() for details of the format.
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
  /* Output file identification and the number of nodes. */
  if(errNum == WLZ_ERR_NONE)
  {
    nElm = mesh->res.elm.numEnt;
    nNod = mesh->res.nod.numEnt;

    if(fprintf(fP,
               "MeshVersionFormatted 1\n"
	       "Dimension\n"
	       "3\n"
	       "Vertices\n"
	       "%d\n", nNod) <= 0)
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
	if(fprintf(fP, " %lg %lg %lg 0\n",
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
    if(fprintf(fP, "Tetrahedra\n"
                   "%d\n", nElm) <= 0)
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
	if(fprintf(fP, "% 8d % 8d % 8d % 8d 0\n",
		   nodTbl[nBuf[0]->idx] + 1, nodTbl[nBuf[1]->idx] + 1,
		   nodTbl[nBuf[2]->idx] + 1, nodTbl[nBuf[3]->idx] + 1) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  /* Output the end of file marker. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP, "End\n") <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  AlcFree(nodTbl);
  return(errNum);
}

/*!
* \return	Number of fields read.
* \ingroup	WlzExtFF
* \brief	Reads the next record with the required number of fields
* 		into the given buffer. On return the fields will be
* 		separated by a single space character.
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
* \param	nFld			Requested number of fields.
*/
static int	WlzEffReadObjMeshRec(FILE *fP, char *buf, int bufLen,
				      int nFld)
{
  int		c0,
  		c1,
		f,
		skip;
  char		*s;
  const char	comment = '#',
  		eol = '\n',
  		space = ' ';

  f = 0;
  s = buf;
  skip = 0;
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
    else if(c1 == comment)
    {
      skip = 1;
    }
    else if(c1 == eol)
    {
      skip = 0;
    }
    if(skip == 0)
    {
      if(isspace(c1))
      {
	c1 = space;
	if(c0 != space)
	{
	  if(++f >= nFld)
	  {
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
  }
  *s = '\0';
  return(f);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads the next record and parses it for an integer number
* 		of the entity to skip. Then skips that many records
* 		with the required number of fields per record.
* \param	fP			Input file.
* \param	buf			Buffer for input record.
* \param	bufLen			Buffer length.
* \param	nFld			Requested number of fields.
*/
static WlzErrorNum	WlzEffReadObjMeshSkip(FILE *fP, char *buf, int bufLen,
					      int nFld)
{
  int		idx,
		nE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((WlzEffReadObjMeshRec(fP, buf, bufLen, 1) != 1) ||
     (sscanf(buf, "%d", &nE) != 1) || (nE < 0))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    for(idx = 0; idx < nE; ++idx)
    {
      if(WlzEffReadObjMeshRec(fP, buf, bufLen, nFld) != nFld)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}
