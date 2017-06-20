#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFVtk_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFVtk.c
* \author       Bill Hill
* \date         March 1999
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
* 		the VTK '.vtk' data format.
* \ingroup	WlzExtFF
*/

#include <errno.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <string.h>
#include <float.h>

static WlzErrorNum 		WlzEffWriteImgVtk(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteCtrVtk(
				  FILE *fP,
				  WlzContour *ctr);
static WlzErrorNum 		WlzEffWriteGMModelVtk(
				  FILE *fP,
				  WlzGMModel *model);
static WlzErrorNum 		WlzEffWriteCMesh2DVtk(
				  FILE *fP,
				  WlzCMesh2D *mesh);
static WlzErrorNum 		WlzEffWriteCMesh2D5Vtk(
				  FILE *fP,
				  WlzCMesh2D5 *mesh);
static WlzErrorNum 		WlzEffWriteCMesh3DVtk(
				  FILE *fP,
				  WlzCMesh3D *mesh);
static WlzErrorNum		WlzEffHeadReadVtk(
				  WlzEffVtkHeader *header,
				  FILE *fP);
static WlzObject		*WlzEffReadCMeshVtk(
				  FILE *fP,
				  WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzEffReadPolyVtk(
				  FILE *fP,
				  WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr);
static WlzObject		*WlzEffReadImgVtk(
				  FILE *fP,
				  WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		Visualization Toolkit file format. There is no check on the
* 		major and minor version numbers since all files seem
* 		compatable.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjVtk(FILE *fP, WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzEffVtkHeader header;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    (void )memset((void *)&header, 0, sizeof(WlzEffVtkHeader));
    errNum = WlzEffHeadReadVtk(&header, fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(header.type)
    {
      case WLZEFF_VTK_TYPE_STRUCTURED_POINTS:
	obj = WlzEffReadImgVtk(fP, &header, &errNum);
        break;
      case WLZEFF_VTK_TYPE_POLYDATA:
	obj = WlzEffReadPolyVtk(fP, &header, &errNum);
        break;
      case WLZEFF_VTK_TYPE_UNSTRUCTURED_GRID:
        obj = WlzEffReadCMeshVtk(fP, &header, &errNum);
	break;
      case WLZEFF_VTK_TYPE_STRUCTURED_GRID:   /* FALLTHROUGH */
      case WLZEFF_VTK_TYPE_RECTILNEAR_GRID:   /* FALLTHROUGH */
      default:
        errNum = WLZ_ERR_READ_INCOMPLETE;
        break;
    }
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
* \brief	Writes the given Woolz object to the given stream using the
* 		Visualization Toolkit file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjVtk(FILE *fP, WlzObject *obj)
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
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_3D_DOMAINOBJ:
        errNum = WlzEffWriteImgVtk(fP, obj);
	break;
      case WLZ_CONTOUR:
        errNum = WlzEffWriteCtrVtk(fP, obj->domain.ctr);
	break;
      case WLZ_CMESH_2D:
        errNum = WlzEffWriteCMesh2DVtk(fP, obj->domain.cm2);
        break;
      case WLZ_CMESH_2D5:
        errNum = WlzEffWriteCMesh2D5Vtk(fP, obj->domain.cm2d5);
        break;
      case WLZ_CMESH_3D:
        errNum = WlzEffWriteCMesh3DVtk(fP, obj->domain.cm3);
        break;
      case WLZ_POINTS:
        errNum = WlzEffWritePointsVtk(fP, obj, 0);
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
* \brief	Writes the given Woolz 3D domain object with values object to
* 		the given stream using the Visualization Toolkit (structured
* 		points) file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
static WlzErrorNum WlzEffWriteImgVtk(FILE *fP, WlzObject *obj)
{
  unsigned char	***data = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  int		dataCount;
  WlzFVertex3	aspect;
  WlzIVertex3	origin,
  		size;

  if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(obj->values.core == NULL)
  {
     errNum = WLZ_ERR_VALUES_NULL;
  }
  else if(obj->values.core->type != WLZ_VOXELVALUETABLE_GREY)
  {
    errNum = WLZ_ERR_VALUES_TYPE;
  }
  else
  {
    origin.vtX = obj->domain.p->kol1;
    origin.vtY = obj->domain.p->line1;
    origin.vtZ = obj->domain.p->plane1;
    aspect.vtX = obj->domain.p->voxel_size[0];
    aspect.vtY = obj->domain.p->voxel_size[1];
    aspect.vtZ = obj->domain.p->voxel_size[2];
    size.vtX = obj->domain.p->lastkl - origin.vtX + 1;
    size.vtY = obj->domain.p->lastln - origin.vtY + 1;
    size.vtZ = obj->domain.p->lastpl - origin.vtZ + 1;
    if((size.vtX <= 0) || (size.vtY <= 0) || (size.vtZ <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      errNum = WlzToArray3D((void ****)&data, obj, size, origin,
			    0, WLZ_GREY_UBYTE);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    dataCount = size.vtX * size.vtY * size.vtZ;
    if((fprintf(fP, 
		"# vtk DataFile Version %d.%d\n"
		"\n"
		"BINARY\n"
		"DATASET STRUCTURED_POINTS\n"
		"DIMENSIONS %d %d %d\n"
		"ORIGIN %d %d %d\n"
		"ASPECT_RATIO %g %g %g\n"
		"POINT_DATA %d\n"
		"SCALARS scalars unsigned_char\n"
		"LOOKUP_TABLE default\n",
		WLZEFF_VTK_VERSION_MAJOR, WLZEFF_VTK_VERSION_MINOR,
		size.vtX, size.vtY, size.vtZ,
		origin.vtX, origin.vtY, origin.vtZ,
		aspect.vtX, aspect.vtY, aspect.vtZ,
		dataCount) <= 0) ||
	(fwrite(**data, sizeof(unsigned char), dataCount, fP) != dataCount))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(data)
  {
    (void )AlcUnchar3Free(data);
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz contour (2D or 3D model) to the
*		given stream using the Visualization Toolkit (polydata) file
*		format.
* \param	fP			Output file stream.
* \param	ctr			Given woolz contour.
*/
static WlzErrorNum WlzEffWriteCtrVtk(FILE *fP, WlzContour *ctr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(ctr == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if(ctr->type != WLZ_CONTOUR)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else
  {
    errNum = WlzEffWriteGMModelVtk(fP, ctr->model);
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz gemetric model (2D or 3D) to the
*		given stream using the Visualization Toolkit (polydata) file
*		format.
* \param	fP			Output file stream.
* \param	model			Given gemetric model.
*/
static WlzErrorNum WlzEffWriteGMModelVtk(FILE *fP, WlzGMModel *model)
{

  int		idI,
  		iCnt;
  int		bufI[3];
  AlcVector	*vec;
  WlzGMEdgeT	*tET;
  WlzGMElemP	eP;
  WlzGMResIdxTb	*resIdxTb = NULL;
  WlzDVertex3	vtx,
  		nrm;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    /* Check the model type. */
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
	break;
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
    /* Index the verticies. */
    resIdxTb = WlzGMModelResIdx(model, WLZ_GMELMFLG_VERTEX, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Check there are verticies! */
    if(resIdxTb->vertex.idxCnt < 1)
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    (void )fprintf(fP,
    		   "# vtk DataFile Version 1.0\n"
    		   "WlzGeoModel test output\n"
		   "ASCII\n"
		   "DATASET POLYDATA\n"
		   "POINTS %d float\n",
		   resIdxTb->vertex.idxCnt);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the vertex geometries. */
    idI = 0;
    vec = model->res.vertex.vec;
    iCnt = model->res.vertex.numIdx;
    while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
    {
      eP.vertex = (WlzGMVertex *)AlcVectorItemGet(vec, idI++);
      if(eP.vertex->idx >= 0)
      {
	(void )WlzGMVertexGetG3D(eP.vertex, &vtx);
	(void )fprintf(fP, "%g %g %g\n", vtx.vtX, vtx.vtY, vtx.vtZ);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the vertex indicies of the simplicies in the order that they
     * are found in the model resouces. */
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D: /* FALLTHROUGH */
      case WLZ_GMMOD_2N:
	(void )fprintf(fP,
		        "LINES %d %d\n",
			model->res.edge.numElm, 3 * model->res.edge.numElm);
	idI = 0;
	vec = model->res.edge.vec;
	iCnt = model->res.edge.numIdx;
	while(iCnt-- > 0)
	{
	  eP.edge = (WlzGMEdge *)AlcVectorItemGet(vec, idI++);
	  if(eP.edge->idx >= 0)
	  {
	    tET = eP.edge->edgeT;
	    bufI[0] = *(resIdxTb->vertex.idxLut +
	    		tET->vertexT->diskT->vertex->idx);
	    bufI[1] = *(resIdxTb->vertex.idxLut +
	    		tET->opp->vertexT->diskT->vertex->idx);
	    (void )fprintf(fP, "2 %d %d\n",
	    		   bufI[0], bufI[1]);
	  }
	}
        break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	(void )fprintf(fP,
		       "POLYGONS %d %d\n",
		       model->res.face.numElm, 4 * model->res.face.numElm);
	idI = 0;
	vec = model->res.face.vec;
	iCnt = model->res.face.numIdx;
	while(iCnt-- > 0)
	{
	  eP.face = (WlzGMFace *)AlcVectorItemGet(vec, idI++);
	  if(eP.face->idx >= 0)
	  {
	    /* Face IS a triangle, in 3D nothing else is allowed. */
	    tET = eP.face->loopT->edgeT;
	    bufI[0] = *(resIdxTb->vertex.idxLut +
	    		tET->vertexT->diskT->vertex->idx);
	    bufI[1] = *(resIdxTb->vertex.idxLut +
	    		tET->next->vertexT->diskT->vertex->idx);
	    bufI[2] = *(resIdxTb->vertex.idxLut +
	    		tET->prev->vertexT->diskT->vertex->idx);
	    (void )fprintf(fP, "3 %d %d %d\n",
	    	           bufI[0], bufI[1], bufI[2]);
	  }
	}
        break;
    }
  }
  /* Output the normals if they are in the model, i.e. the model is
   * either a WLZ_GMMOD_2N or WLZ_GMMOD_3N. */
  if((errNum == WLZ_ERR_NONE)  &&
     ((model->type == WLZ_GMMOD_2N) || (model->type == WLZ_GMMOD_3N)))
  {
    (void )fprintf(fP,
		   "\n"
		   "POINT_DATA %d\n"
		   "NORMALS normals float\n",
		   resIdxTb->vertex.idxCnt);
    /* Output the vertex normals. */
    idI = 0;
    vec = model->res.vertex.vec;
    iCnt = model->res.vertex.numIdx;
    while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
    {
      eP.vertex = (WlzGMVertex *)AlcVectorItemGet(vec, idI++);
      if(eP.vertex->idx >= 0)
      {
	(void )WlzGMVertexGetG3N(eP.vertex, &vtx, &nrm);
	(void )fprintf(fP, "%g %g %g\n", nrm.vtX, nrm.vtY, nrm.vtZ);
      }
    }
  }
  if(resIdxTb)
  {
    WlzGMModelResIdxFree(resIdxTb);
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz 2D constrained mesh to the
*		given stream using the Visualization Toolkit
*		unstructured grid file format with triangular elements.
* \param	fP			Output file stream.
* \param	model			Given gemetric model.
*/
static WlzErrorNum WlzEffWriteCMesh2DVtk(FILE *fP, WlzCMesh2D *mesh)
{

  int		cnt,
		idE,
		idN,
  		nElm,
  		nNod;
  WlzCMeshElm2D	*elm;
  WlzCMeshNod2D	*nod[3];
  int		*nodTbl = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(((nNod = mesh->res.nod.numEnt) < 3) ||
          ((nElm = mesh->res.elm.numEnt) < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate a node table to avoid deleted nodes. */
    if((nodTbl = (int *)AlcMalloc(sizeof(int) *
                                  mesh->res.nod.maxEnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    if(fprintf(fP,
	       "# vtk DataFile Version 1.0\n"
	       "Written by WlzEffWriteCMesh2DVtk().\n"
	       "ASCII\n"
	       "DATASET UNSTRUCTURED_GRID\n"
	       "POINTS %d float\n",
	       nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the node positions while building a table of valid nodes. */
    cnt = 0;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod[0] = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod[0]->idx >= 0)
      {
        if(fprintf(fP, "%g %g 0.0\n",
	               nod[0]->pos.vtX, nod[0]->pos.vtY) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
	nodTbl[idN] = cnt++;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element node indices. */
    if(fprintf(fP,
               "CELLS %d %d\n", nElm, 4 * nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	nod[0] = WLZ_CMESH_ELM2D_GET_NODE_0(elm);
	nod[1] = WLZ_CMESH_ELM2D_GET_NODE_1(elm);
	nod[2] = WLZ_CMESH_ELM2D_GET_NODE_2(elm);
        if(fprintf(fP, "3 %d %d %d\n",
	           nodTbl[nod[0]->idx], nodTbl[nod[2]->idx],
	           nodTbl[nod[1]->idx]) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  AlcFree(nodTbl);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element cell types (all triangles, type 5). */
    if(fprintf(fP,
               "CELL_TYPES %d\n", nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      if(fprintf(fP, "5\n") <= 0)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz points object to the
*		given stream using the Visualization Toolkit
*		polydata format, but without polygons.
* \param	fP			Output file stream.
* \param	model			Given points object.
* \param        onlyDom			Only write the points domain if
* 					non-zero.
*/
WlzErrorNum 	WlzEffWritePointsVtk(FILE *fP, WlzObject *obj, int onlyDom)
{

  WlzPoints	 *pdm;
  WlzPointValues *pvl;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((pdm = obj->domain.pts) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pdm->type)
    {
      case WLZ_POINTS_2I: /* FALLTHROUGH */
      case WLZ_POINTS_2D: /* FALLTHROUGH */
      case WLZ_POINTS_3I: /* FALLTHROUGH */
      case WLZ_POINTS_3D:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    if(fprintf(fP,
	       "# vtk DataFile Version 1.0\n"
	       "Written by WlzEffWritePointsVtk().\n"
	       "ASCII\n"
	       "DATASET POLYDATA\n"
	       "POINTS %d float\n",
	       pdm->nPoints) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  /* Output the point domain. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    switch(pdm->type)
    {
      case WLZ_POINTS_2I:
	for(i = 0; i < pdm->nPoints; ++i)
        {
	  WlzIVertex2 *p;

	  p = pdm->points.i2 + i;
	  if(fprintf(fP, "%g %g 0\n",
	             (double )(p->vtX), (double )(p->vtY)) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
	break;
      case WLZ_POINTS_2D:
	for(i = 0; i < pdm->nPoints; ++i)
        {
	  WlzDVertex2 *p;

	  p = pdm->points.d2 + i;
	  if(fprintf(fP, "%g %g 0\n", p->vtX, p->vtY) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
	break;
      case WLZ_POINTS_3I:
	for(i = 0; i < pdm->nPoints; ++i)
        {
	  WlzIVertex3 *p;

	  p = pdm->points.i3 + i;
	  if(fprintf(fP, "%g %g %g\n",
	             (double )(p->vtX), (double )(p->vtY),
		     (double )(p->vtZ)) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
	break;
      case WLZ_POINTS_3D:
	for(i = 0; i < pdm->nPoints; ++i)
        {
	  WlzDVertex3 *p;

	  p = pdm->points.d3 + i;
	  if(fprintf(fP, "%g %g %g\n", p->vtX, p->vtY, p->vtZ) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    break;
	  }
	}
        break;
      default:
	break;
    }
  }
  /* Output point values if they exist. */
  if((errNum == WLZ_ERR_NONE) && (onlyDom == 0) &&
     ((pvl = obj->values.pts) != NULL))
  {
    if(pvl->rank == 0)
    {
      errNum = WlzEffWritePointsVtkScalarValues(fP, obj);
    }
    else
    {
      errNum = WlzEffWritePointsVtkFieldValues(fP, obj);
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the point scalar values of the given object to the given
* 		file. It is assumed that the points domain has already been
* 		written by WlzEffWritePointsVtk().
* \param	fP			Output file stream.
* \param	obj			Object vith points values for output.
*/
WlzErrorNum			WlzEffWritePointsVtkScalarValues(
				  FILE *fP, 
				  WlzObject *obj)
{
  WlzPoints	 *pdm;
  WlzPointValues *pvl;
  char		 *dStr = "scalars",
  		 *nStr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((pdm = obj->domain.pts) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pdm->type)
    {
      case WLZ_POINTS_2I: /* FALLTHROUGH */
      case WLZ_POINTS_2D: /* FALLTHROUGH */
      case WLZ_POINTS_3I: /* FALLTHROUGH */
      case WLZ_POINTS_3D:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && obj->plist)
  {
    WlzNameProperty *np;

    if((np = WlzPropertyListContainsName(obj->plist, NULL)) != NULL)
    {
      dStr = np->name;
    }
    if((nStr = WlzStringCopyReplace(dStr, " \t\n", '_', 0)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && ((pvl = obj->values.pts) != NULL))
  {
    if(pvl->rank == 0)
    {
      const char *gStr;
      int	 unsup = 1;

      switch(pvl->vType)
      {
        case WLZ_GREY_UBYTE:
	  unsup = 0;
	  gStr = "char";
	  break;
        case WLZ_GREY_INT:   /* FALLTHROUGH */
	case WLZ_GREY_SHORT: /* FALLTHROUGH */
	case WLZ_GREY_FLOAT: /* FALLTHROUGH */
	case WLZ_GREY_DOUBLE:
	  unsup = 0;
	  gStr = "float";
	  break;
	default:
	  break;
      }
      if(unsup == 0)
      {
	if(fprintf(fP,
		   "POINT_DATA %d\n"
		   "SCALARS %s %s\n"
		   "LOOKUP_TABLE default\n",
		   pdm->nPoints, nStr, gStr) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	else
	{
          int	i;

	  switch(pvl->vType)
	  {
	    case WLZ_GREY_UBYTE:
	      for(i = 0; i < pdm->nPoints; ++i)
	      {
	        if(fprintf(fP, "%d\n", pvl->values.ubp[i]) <= 0)
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		  break;
		}
	      }
	      break;
	    case WLZ_GREY_INT:
	      for(i = 0; i < pdm->nPoints; ++i)
	      {
	        if(fprintf(fP, "%g\n", (float )(pvl->values.inp[i])) <= 0)
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		  break;
		}
	      }
	      break;
	    case WLZ_GREY_SHORT:
	      for(i = 0; i < pdm->nPoints; ++i)
	      {
	        if(fprintf(fP, "%g\n", (float )(pvl->values.shp[i])) <= 0)
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		  break;
		}
	      }
	      break;
	    case WLZ_GREY_FLOAT:
	      for(i = 0; i < pdm->nPoints; ++i)
	      {
	        if(fprintf(fP, "%g\n", pvl->values.flp[i]) <= 0)
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		  break;
		}
	      }
	      break;
	    case WLZ_GREY_DOUBLE:
	      for(i = 0; i < pdm->nPoints; ++i)
	      {
	        if(fprintf(fP, "%g\n", (float )(pvl->values.dbp[i])) <= 0)
		{
		  errNum = WLZ_ERR_WRITE_INCOMPLETE;
		  break;
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_GREY_TYPE;
	      break;
	  }
	}
      }
    }
  }
  AlcFree(nStr);
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the point field values of the given object to the given
* 		file. It is assumed that the points domain has already been
* 		written by WlzEffWritePointsVtk().
*
* 		The VTK file specification does not include an appropriate
* 		data type for arbitrary dimensional fields so here one is
* 		defined for point data:
* 		\verbatim
	POINT_DATA <n>
	FIELD <field name> <data type> <rank> [<dim_0>,[...]]
		\endverbatim
*		with the following:
*		<ul>
*		<li>n - Number of points.
*		<li>field name - Non-whitespace ASCII string naming the field.
*		<li>rank - Rank of the field, 0 for scalar data, 1 for vector
*		           data, ....
*		<li>dim0 - Comma seperated dimensions of the field.
*		</ul>
*
*		Example for a rank two (3x3) tensor:
*		\verbatim
	POINT_DATA 2
	FIELD my_tensor float 2 3x3
	1.019 -0.002 0.079  0.003 1.055  0.012 0.002   0.0085 1.019
	1.012 -0.008 0.015 -0.023 1.013 -0.018 7.9e-05 0.0961 1.005
		\endverbatim
* \param	fP			Output file stream.
* \param	obj			Object vith points values for output.
*/
WlzErrorNum			WlzEffWritePointsVtkFieldValues(
				  FILE *fP, 
				  WlzObject *obj)
{
  WlzPoints	 *pdm;
  WlzPointValues *pvl;
  char		 *dStr = "field",
  		 *nStr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((pdm = obj->domain.pts) == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    switch(pdm->type)
    {
      case WLZ_POINTS_2I: /* FALLTHROUGH */
      case WLZ_POINTS_2D: /* FALLTHROUGH */
      case WLZ_POINTS_3I: /* FALLTHROUGH */
      case WLZ_POINTS_3D:
        break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  if((errNum == WLZ_ERR_NONE) && obj->plist)
  {
    WlzNameProperty *np;

    if((np = WlzPropertyListContainsName(obj->plist, NULL)) != NULL)
    {
      dStr = np->name;
    }
    if((nStr = WlzStringCopyReplace(dStr, " \t\n", '_', 0)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if((errNum == WLZ_ERR_NONE) && ((pvl = obj->values.pts) != NULL))
  {
    if(pvl->rank >= 0)
    {
      int	 vpe = 1,
      		 unsup = 1;
      const char *gStr;

      switch(pvl->vType)
      {
        case WLZ_GREY_UBYTE:
	  unsup = 0;
	  gStr = "char";
	  break;
        case WLZ_GREY_INT:   /* FALLTHROUGH */
	case WLZ_GREY_SHORT: /* FALLTHROUGH */
	case WLZ_GREY_FLOAT: /* FALLTHROUGH */
	case WLZ_GREY_DOUBLE:
	  unsup = 0;
	  gStr = "float";
	  break;
	default:
	  break;
      }
      if(unsup == 0)
      {
	if(fprintf(fP,
		   "POINT_DATA %d\n"
		   "FIELD %s %s %d ",
		   pdm->nPoints, nStr, gStr, pvl->rank) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	else if(pvl->rank >= 0)
	{
	  int	i,
	  	penult;

	  penult = pvl->rank - 1;
	  for(i = 0; i < penult; ++i)
	  {
	    vpe *= pvl->dim[i];
	    (void )fprintf(fP, "%dx", pvl->dim[i]);
	  }
	  vpe *= pvl->dim[i];
	  if(vpe < 1)
	  {
	    errNum = WLZ_ERR_VALUES_DATA;
	  }
	  else if(fprintf(fP, "%d\n", pvl->dim[i]) <= 0)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
          int		i;
	  size_t	cnt = 0;
	  WlzGreyP	p;

	  for(i = 0; i < pdm->nPoints; ++i)
	  {
	    int		j;

	    switch(pvl->vType)
	    {
	      case WLZ_GREY_UBYTE:
		p.ubp = pvl->values.ubp + cnt;
		for(j = 0; j < (vpe - 1); ++j)
		{
		  (void )fprintf(fP, "%d ", p.ubp[j]);
		}
		(void )fprintf(fP, "%d", p.ubp[j]);
		break;
	      case WLZ_GREY_INT:
		p.inp = pvl->values.inp + cnt;
		for(j = 0; j < (vpe - 1); ++j)
		{
		  (void )fprintf(fP, "%d ", p.inp[j]);
		}
		(void )fprintf(fP, "%d", p.inp[j]);
		break;
	      case WLZ_GREY_SHORT:
		p.shp = pvl->values.shp + cnt;
		for(j = 0; j < (vpe - 1); ++j)
		{
		  (void )fprintf(fP, "%d ", p.shp[j]);
		}
		(void )fprintf(fP, "%d", p.shp[j]);
		break;
	      case WLZ_GREY_FLOAT:
		p.flp = pvl->values.flp + cnt;
		for(j = 0; j < (vpe - 1); ++j)
		{
		  (void )fprintf(fP, "%g ", p.flp[j]);
		}
		(void )fprintf(fP, "%g", p.flp[j]);
		break;
	      case WLZ_GREY_DOUBLE:
		p.dbp = pvl->values.dbp + cnt;
		for(j = 0; j < (vpe - 1); ++j)
		{
		  (void )fprintf(fP, "%lg ", p.dbp[j]);
		}
		(void )fprintf(fP, "%lg", p.dbp[j]);
		break;
	      default:
		errNum = WLZ_ERR_GREY_TYPE;
		break;
	    }
	    cnt += vpe;
	    if(fprintf(fP, "\n") <= 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	      break;
	    }
	  }
        }
      }
    }
  }
  AlcFree(nStr);
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz 2D5 constrained mesh to the
*		given stream using the Visualization Toolkit
*		unstructured grid file format with triangular elements.
* \param	fP			Output file stream.
* \param	model			Given gemetric model.
*/
static WlzErrorNum WlzEffWriteCMesh2D5Vtk(FILE *fP, WlzCMesh2D5 *mesh)
{

  int		cnt,
		idE,
		idN,
  		nElm,
  		nNod;
  WlzCMeshElm2D5 *elm;
  WlzCMeshNod2D5 *nod[3];
  int		*nodTbl = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_2D5)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(((nNod = mesh->res.nod.numEnt) < 3) ||
          ((nElm = mesh->res.elm.numEnt) < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate a node table to avoid deleted nodes. */
    if((nodTbl = (int *)AlcMalloc(sizeof(int) *
                                  mesh->res.nod.maxEnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    if(fprintf(fP,
	       "# vtk DataFile Version 1.0\n"
	       "Written by WlzEffWriteCMesh2D5Vtk().\n"
	       "ASCII\n"
	       "DATASET UNSTRUCTURED_GRID\n"
	       "POINTS %d float\n",
	       nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the node positions while building a table of valid nodes. */
    cnt = 0;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod[0] = (WlzCMeshNod2D5 *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod[0]->idx >= 0)
      {
        if(fprintf(fP, "%g %g %g\n",
	               nod[0]->pos.vtX, nod[0]->pos.vtY, nod[0]->pos.vtZ) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
	nodTbl[idN] = cnt++;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element node indices. */
    if(fprintf(fP,
               "CELLS %d %d\n", nElm, 4 * nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm2D5 *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	nod[0] = WLZ_CMESH_ELM2D5_GET_NODE_0(elm);
	nod[1] = WLZ_CMESH_ELM2D5_GET_NODE_1(elm);
	nod[2] = WLZ_CMESH_ELM2D5_GET_NODE_2(elm);
        if(fprintf(fP, "3 %d %d %d\n",
	           nodTbl[nod[0]->idx], nodTbl[nod[2]->idx],
	           nodTbl[nod[1]->idx]) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  AlcFree(nodTbl);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element cell types (all triangles, type 5). */
    if(fprintf(fP,
               "CELL_TYPES %d\n", nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      if(fprintf(fP, "5\n") <= 0)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz 3D constrained mesh to the
*		given stream using the Visualization Toolkit
*		unstructured grid file format with tetrahedral elements.
* \param	fP			Output file stream.
* \param	model			Given gemetric model.
*/
static WlzErrorNum WlzEffWriteCMesh3DVtk(FILE *fP, WlzCMesh3D *mesh)
{

  int		cnt,
		idE,
		idN,
  		nElm,
  		nNod;
  WlzCMeshElm3D	*elm;
  WlzCMeshNod3D	*nod[4];
  int		*nodTbl = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(mesh == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(mesh->type != WLZ_CMESH_3D)
  {
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }
  else if(((nNod = mesh->res.nod.numEnt) < 4) ||
          ((nElm = mesh->res.elm.numEnt) < 1))
  {
    errNum = WLZ_ERR_DOMAIN_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Allocate a node table to avoid deleted nodes. */
    if((nodTbl = (int *)AlcMalloc(sizeof(int) *
                                  mesh->res.nod.maxEnt)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the file header. */
    if(fprintf(fP,
	       "# vtk DataFile Version 1.0\n"
	       "Written by WlzEffWriteCMesh3DVtk().\n"
	       "ASCII\n"
	       "DATASET UNSTRUCTURED_GRID\n"
	       "POINTS %d float\n",
	       nNod) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the node positions while building a table of valid nodes. */
    cnt = 0;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod[0] = (WlzCMeshNod3D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod[0]->idx >= 0)
      {
        if(fprintf(fP, "%g %g %g\n",
	               nod[0]->pos.vtX, nod[0]->pos.vtY, nod[0]->pos.vtZ) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
	nodTbl[idN] = cnt++;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element node indices. */
    if(fprintf(fP,
               "CELLS %d %d\n", nElm, 5 * nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
    {
      elm = (WlzCMeshElm3D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
      if(elm->idx >= 0)
      {
	nod[0] = WLZ_CMESH_ELM3D_GET_NODE_0(elm);
	nod[1] = WLZ_CMESH_ELM3D_GET_NODE_1(elm);
	nod[2] = WLZ_CMESH_ELM3D_GET_NODE_2(elm);
	nod[3] = WLZ_CMESH_ELM3D_GET_NODE_3(elm);
        if(fprintf(fP, "4 %d %d %d %d\n",
	           nodTbl[nod[0]->idx], nodTbl[nod[1]->idx],
	           nodTbl[nod[3]->idx], nodTbl[nod[2]->idx]) <= 0)
        {
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  break;
	}
      }
    }
  }
  AlcFree(nodTbl);
  if(errNum == WLZ_ERR_NONE)
  {
    /* Output the element cell types (all tetrahedra, type 10). */
    if(fprintf(fP,
               "CELL_TYPES %d\n", nElm) <= 0)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    for(idE = 0; idE < nElm; ++idE)
    {
      if(fprintf(fP, "10\n") <= 0)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
	break;
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Reads the VTK file header which is:
*		<ul>
*		  <li> ident\n
*		  <li> title\n
*		  <li> data type\n
*		  <li> type\n
*		</ul>
*		where
*		<ul>
*		  <li> ident = "vtk DataFile Version MAJOR.MINOR"
*		  <li> title = less than 256 chars followed by a \n
*		  <li> data type = "ASCII"|"BINARY"
*		  <li> type = DATASET "STRUCTURED_POINTS"|"STRUCTURED_GRID"|
*		              "UNSTRUCTURED_GRID"|"POLYDATA"|
*		              "RECTILNEAR_GRID"
*		</ul>
* \param	header			Given header data structure be filled.
* \param	fP			Output file stream.
*/
static WlzErrorNum WlzEffHeadReadVtk(WlzEffVtkHeader *header, FILE *fP)
{
  int		valI;
  char		*valS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		buf[256];

  /* Read and check the identifier, major and minor version numbers are
   * assumed to be ok. */
  if(fgets(buf, 256, fP) == NULL)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    buf[255] = '\0';
    if(sscanf(buf, "# vtk DataFile Version %d.%d",
	      &(header->versionMajor), &(header->versionMinor)) != 2)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    /* Discard the version because we can probably read it anyway. */
  }
  /* Read the title. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fgets(header->title, 256, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      buf[255] = '\0';
    }
  }
  /* Read and parse the data type. */
  if(errNum == WLZ_ERR_NONE)
  {
    valI = WLZEFF_INVALID;
    do
    {
      if(fgets(buf, 256, fP) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      buf[255] = '\0';
      valS = strtok(buf, " \t\n\r\f\v");
      if((valS != NULL) && (*valS != '#'))
      {
	if(WlzStringMatchValue(&valI, valS,
			       "ASCII", WLZEFF_VTK_DATATYPE_ASCII,
			       "BINARY", WLZEFF_VTK_DATATYPE_BINARY,
			       NULL) == 0)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	header->dataType = valI;
      }
    } while(valI == WLZEFF_INVALID);
  }
  /* Read and parse the type. */
  if(errNum == WLZ_ERR_NONE)
  {
    valI = WLZEFF_INVALID;
    do
    {
      if(fgets(buf, 256, fP) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      buf[255] = '\0';
      valS = strtok(buf, " \t\n\r\f\v");
      if((valS != NULL) && (*valS != '#'))
      {
	if(strcmp(valS, "DATASET"))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	valS = strtok(NULL, " \t\n\r\f\v");
	if((valS == NULL) ||
	   (WlzStringMatchValue(&valI, valS,
		"STRUCTURED_POINTS", WLZEFF_VTK_TYPE_STRUCTURED_POINTS,
		"STRUCTURED_GRID", WLZEFF_VTK_TYPE_STRUCTURED_GRID,
		"UNSTRUCTURED_GRID", WLZEFF_VTK_TYPE_UNSTRUCTURED_GRID,
		"POLYDATA", WLZEFF_VTK_TYPE_POLYDATA,
		"RECTILNEAR_GRID", WLZEFF_VTK_TYPE_RECTILNEAR_GRID,
		NULL) == 0))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
	}
	header->type = valI;
      }
    } while(valI == WLZEFF_INVALID);
  }
  return(errNum);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz CMESH object from the given stream using
*		the Visualization Toolkit (unstructured grid) file format.
* \param	fP			Input file stream.
* \param	header			Header data structure.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadCMeshVtk(FILE *fP, WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr)
{
  int		cnt,
		idE,
		idN,
		dim,
  		valI,
		nElm = 0,
		nNod = 0;
  char 		*valS0,
  		*valS1;
  WlzDBox3	bBox;
  WlzCMeshP	mesh;
  WlzEffVtkUnstructuredGridType prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDVertex2	pos2D;
  WlzDVertex3	*vBuf = NULL;
  WlzCMeshNod2D	*nBuf2[3];
  WlzCMeshNod3D	*nBuf3[4];
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*obj = NULL;
  int		eBuf[4];
  char		buf[256];

  mesh.v = NULL;
  dom.core = NULL;
  val.core = NULL;
  if(header->dataType != WLZEFF_VTK_DATATYPE_ASCII)
  {
    /* Can only read ascii meshes. */
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    do
    {
      /* Read line containing token. */
      valS0 = NULL;
      if(fgets(buf, 256, fP) == NULL)
      {
	if(mesh.v == NULL)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
        break;
      }
      buf[255] = '\0';
      valS0 = strtok(buf, " \t\n\r\f\v");
      if(valS0)
      {
	/* Any other tokens apart from these are ignored. */
	if((WlzStringMatchValue(&valI, valS0,
		"POINTS", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINTS,
		"CELLS", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELLS,
		"CELL_TYPES", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELL_TYPES,
		"POINT_DATA", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINT_DATA,
		"VECTORS", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_VECTORS,
		"SCALARS", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_SCALARS,
		"LOOKUP_TABLE", WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_LOOKUP_TABLE,
		NULL) != 0))
	{
	  prim = valI;
	  switch(prim)
	  {
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINTS:
	      if(nNod != 0)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		valS0 = strtok(NULL, " \t\n\r\f\v");
		if((valS0 == NULL) || (sscanf(valS0, "%d", &nNod) != 1) ||
		    (nNod <= 0))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
		else
		{
		  valS0 = strtok(NULL, " \t\n\r\f\v");
		  if((valS0 == NULL) || strcmp(valS0, "float"))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  /* Allocate a 3D vertex buffer. 3D vertices are used even
		   * for 2D meshes. */
		  if((vBuf = AlcMalloc(sizeof(WlzDVertex3) * nNod)) == NULL)
		  {
		    errNum = WLZ_ERR_MEM_ALLOC;
		  }
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  /* Read nNod 3D verticies into the 3D vertex buffer. */
		  for(idN = 0; idN < nNod; ++idN)
		  {
		    if((fscanf(fP, "%lg", &(vBuf[idN].vtX)) != 1) ||
		       (fscanf(fP, "%lg", &(vBuf[idN].vtY)) != 1) ||
		       (fscanf(fP, "%lg", &(vBuf[idN].vtZ)) != 1))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		      break;
		    }
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
		  }
		}
		break;
	      case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELLS:
		if(nNod <= 0)
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  valS0 = strtok(NULL, " \t\n\r\f\v");
		  valS1 = strtok(NULL, " \t\n\r\f\v");
		  if((valS0 == NULL) || (valS1 == NULL) ||
		      (sscanf(valS0, "%d", &nElm) != 1) ||
		      (sscanf(valS1, "%d", &cnt) != 1) ||
		      (nElm < 0) || (cnt <= nElm))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  /* Assume that there will only be triangular or tetrhedral
		   * mesh elements with either (3 + 1) or (4 + 1) fields
		   * per record for 2D and 3D. If this isn't the case can
		   * return an error later. */
		  dim = (cnt / nElm) - 2;
		  switch(dim)
		  {
		    case 2:
		      if(errNum == WLZ_ERR_NONE)
		      {
			mesh.m2 = WlzCMeshNew2D(&errNum);
		      }
		      if(AlcVectorExtendAndGet(mesh.m2->res.nod.vec,
		                               nNod) == NULL)
		      {
		        errNum = WLZ_ERR_MEM_ALLOC;
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			mesh.m2->bBox.xMin = bBox.xMin;
			mesh.m2->bBox.yMin = bBox.yMin;
			mesh.m2->bBox.xMax = bBox.xMax;
			mesh.m2->bBox.yMax = bBox.yMax;
			errNum = WlzCMeshReassignGridCells2D(mesh.m2, nNod);
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			for(idN = 0; idN < nNod; ++idN)
			{
			  pos2D.vtX = vBuf[idN].vtX;
			  pos2D.vtY = vBuf[idN].vtY;
			  nBuf2[0] = WlzCMeshNewNod2D(mesh.m2, pos2D, NULL);
			  nBuf2[0]->flags = 0;
			}
			for(idE = 0; idE < nElm; ++idE)
			{
			  if((fscanf(fP, "%d", &valI) != 1) || (valI != 3) ||
			      (fscanf(fP, "%d %d %d",
				      eBuf + 0, eBuf + 1, eBuf + 2) != 3) ||
			      (eBuf[0] < 0) || (eBuf[0] >= nNod) ||
			      (eBuf[1] < 0) || (eBuf[1] >= nNod) ||
			      (eBuf[2] < 0) || (eBuf[2] >= nNod))
			  {
			    errNum = WLZ_ERR_READ_INCOMPLETE;
			    break;
			  }
			  for(idN = 0; idN < 3; ++idN)
			  {
			    nBuf2[idN] = (WlzCMeshNod2D *)
			                 AlcVectorItemGet(mesh.m2->res.nod.vec,
				                          eBuf[idN]);
			  }
			  /* Add triangle to the mesh. It's possible that the
			   * orienation will be wrong so try both. */
			  (void )WlzCMeshNewElm2D(mesh.m2, nBuf2[0], nBuf2[2],
						  nBuf2[1], 1, &errNum);
			}
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			WlzCMeshUpdateMaxSqEdgLen2D(mesh.m2);
		      }
		      break;
		    case 3:
		      if(errNum == WLZ_ERR_NONE)
		      {
			mesh.m3 = WlzCMeshNew3D(&errNum);
		      }
		      if(AlcVectorExtendAndGet(mesh.m3->res.nod.vec,
		                               nNod) == NULL)
		      {
		        errNum = WLZ_ERR_MEM_ALLOC;
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			mesh.m3->bBox = bBox;
			errNum = WlzCMeshReassignGridCells3D(mesh.m3, nNod);
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			for(idN = 0; idN < nNod; ++idN)
			{
			  nBuf3[0] = WlzCMeshNewNod3D(mesh.m3, vBuf[idN],
			                              NULL);
			  nBuf3[0]->flags = 0;
			}
			for(idE = 0; idE < nElm; ++idE)
			{
			  if((fscanf(fP, "%d", &valI) != 1) || (valI != 4) ||
			      (fscanf(fP, "%d %d %d %d",
				      eBuf + 0, eBuf + 1,
				      eBuf + 2, eBuf + 3) != 4) ||
			      (eBuf[0] < 0) || (eBuf[0] >= nNod) ||
			      (eBuf[1] < 0) || (eBuf[1] >= nNod) ||
			      (eBuf[2] < 0) || (eBuf[2] >= nNod) ||
			      (eBuf[3] < 0) || (eBuf[3] >= nNod))
			  {
			    errNum = WLZ_ERR_READ_INCOMPLETE;
			    break;
			  }
			  for(idN = 0; idN < 4; ++idN)
			  {
			    nBuf3[idN] = (WlzCMeshNod3D *)
			                 AlcVectorItemGet(mesh.m3->res.nod.vec,
				                          eBuf[idN]);
			  }
			  /* Add tetrahedron to the mesh. It's possible that
			   * the orienation will be wrong so try both. */
			  (void )WlzCMeshNewElm3D(mesh.m3, nBuf3[0], nBuf3[1],
						  nBuf3[3], nBuf3[2], 1, &errNum);
			}
		      }
		      if(errNum == WLZ_ERR_NONE)
		      {
			WlzCMeshUpdateMaxSqEdgLen3D(mesh.m3);
		      }
		      break;
		    default:
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		      break;
		  }
		}
	      }
	      break;
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_CELL_TYPES:
	      if(mesh.v == NULL)
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		valS0 = strtok(NULL, " \t\n\r\f\v");
		if((valS0 == NULL) || (sscanf(valS0, "%d", &cnt) != 1))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        switch(mesh.m2->type)
		{
		  case WLZ_CMESH_2D:
		    /* Read cell types, just make sure that they are all
		     * triangles (type 5).  */
		    for(idE = 0; idE < cnt; ++idE)
		    {
		      if((fscanf(fP, "%d", &valI) != 1) || (valI != 5))
		      {
			errNum = WLZ_ERR_READ_INCOMPLETE;
			break;
		      }
		    }
		    break;
		  case WLZ_CMESH_3D:
		    /* Read cell types, just make sure that they are all
		     * tetrahedra (type 10).  */
		    for(idE = 0; idE < cnt; ++idE)
		    {
		      if((fscanf(fP, "%d", &valI) != 1) || (valI != 10))
		      {
			errNum = WLZ_ERR_READ_INCOMPLETE;
			break;
		      }
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		    break;
		}
	      }
	      break;
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_POINT_DATA:
	      /* TODO: Point data not handled yet. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_VECTORS:
	      /* TODO: Point vector data not handled yet. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_SCALARS:
	      /* TODO: Point scalar data not handled yet. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_UNSTRUCTUREDGRIDTYPE_LOOKUP_TABLE:
	      /* TODO: Point LUT data not handled yet. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    default:
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	  }
	}
      }
    }
    while(errNum == WLZ_ERR_NONE);
  }
  AlcFree(vBuf);
  if((errNum == WLZ_ERR_NONE) && (mesh.v != NULL))
  {
    switch(mesh.m2->type)
    {
      case WLZ_CMESH_2D:
	WlzCMeshUpdateMaxSqEdgLen2D(mesh.m2);
	dom.cm2 = mesh.m2;
	obj = WlzMakeMain(WLZ_CMESH_2D, dom, val, NULL, NULL, &errNum);
	break;
      case WLZ_CMESH_3D:
	WlzCMeshUpdateMaxSqEdgLen3D(mesh.m3);
	dom.cm3 = mesh.m3;
	obj = WlzMakeMain(WLZ_CMESH_3D, dom, val, NULL, NULL, &errNum);
	break;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(mesh.v != NULL)
    {
      switch(mesh.m2->type)
      {
	case  WLZ_CMESH_2D:
	  (void )WlzCMeshFree2D(mesh.m2);
	  break;
	case WLZ_CMESH_3D:
	  (void )WlzCMeshFree3D(mesh.m3);
	  break;
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a WlzGMModel from the given stream using the
*		Visualization Toolkit (polydata) file format.
* \param	fP			Input file stream.
* \param	header			Header data structure.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadPolyVtk(FILE *fP, WlzEffVtkHeader *header,
				   WlzErrorNum *dstErr)
{
  int		valI,
		pIdx,
		nLine = 0,
		nPoly = 0,
		nPoints = 0,
		nPointData = 0,
		sumPoints = 0,
		endOfData = 0,
		vHTSz;
  char 		*valS = NULL;
  WlzGMModel	*model = NULL;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDVertex3	*pointBuf = NULL;
  float		*pointDataBuf = NULL;
  char		buf[256];

  dom.core = NULL;
  val.core = NULL;
  if(header->dataType == WLZEFF_VTK_DATATYPE_BINARY)
  {
    /* Can only read ascii polydata. */
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    do
    {
      WlzEffVtkPolyDataType prim;

      /* Read line containing token. */
      if(fgets(buf, 256, fP) == NULL)
      {
        endOfData = 1;
      }
      else
      {
	buf[255] = '\0';
	valS = strtok(buf, " \t\n\r\f\v");
      }
      if(valS)
      {
	/* Any other tokens apart from these are ignored. */
	if((WlzStringMatchValue(&valI, valS,
		"POINTS", WLZEFF_VTK_POLYDATATYPE_POINTS,
		"VERTICIES", WLZEFF_VTK_POLYDATATYPE_VERTICIES,
		"LINES", WLZEFF_VTK_POLYDATATYPE_LINES,
		"POLYGONS", WLZEFF_VTK_POLYDATATYPE_POLYGONS,
		"TRIANGLE_STRIPS", WLZEFF_VTK_POLYDATATYPE_TRIANGLE_STRIPS,
		"POINT_DATA", WLZEFF_VTK_POLYDATATYPE_POINT_DATA,
		"VECTORS", WLZEFF_VTK_POLYDATATYPE_VECTORS,
		"SCALARS", WLZEFF_VTK_POLYDATATYPE_SCALARS,
		"LOOKUP_TABLE", WLZEFF_VTK_POLYDATATYPE_LOOKUP_TABLE,
		NULL) != 0))
	{
	  prim = valI;
	  switch(prim)
	  {
	    case WLZEFF_VTK_POLYDATATYPE_POINTS:
	      valS = strtok(NULL, " \t\n\r\f\v");
	      if((valS == NULL) || (sscanf(valS, "%d", &nPoints) != 1) ||
		 (nPoints <= 0))
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		  valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || strcmp(valS, "float"))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        if((pointBuf = (WlzDVertex3 *)AlcRealloc(pointBuf,
			              sizeof(WlzDVertex3) * nPoints)) == NULL)
	        {
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        /* Read nPoints 3D verticies. */
		pIdx = 0;
		while((errNum == WLZ_ERR_NONE) && (pIdx < nPoints))
		{
		  if((fscanf(fP, "%lg", &((pointBuf + pIdx)->vtX)) != 1) ||
		     (fscanf(fP, "%lg", &((pointBuf + pIdx)->vtY)) != 1) ||
		     (fscanf(fP, "%lg", &((pointBuf + pIdx)->vtZ)) != 1))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  ++pIdx;
		}
		sumPoints += nPoints;
	      }
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_POLYGONS:
	      if(obj == NULL)
	      {
		vHTSz = (sumPoints < 1024)? 1024: sumPoints / 3;
		model = WlzGMModelNew(WLZ_GMMOD_3D, 0, vHTSz, &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  dom.ctr = WlzMakeContour(&errNum);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    dom.ctr->model = WlzAssignGMModel(model, NULL);
	            obj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL,
		                      &errNum);
		    if(obj == NULL)
		    {
		      (void )WlzFreeContour(dom.ctr);
		    }
		  }
		  else
		  {
		    (void )WlzGMModelFree(model);
		  }
		}
	      }
	      else if(model && (model->type == WLZ_GMMOD_3D))
	      {
		vHTSz = (sumPoints < 1024)? 1024: sumPoints / 3;
		if(model->vertexHTSz < vHTSz)
		{
		  errNum = WlzGMModelRehashVHT(model, vHTSz);
		}
	      }
	      else
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || (sscanf(valS, "%d", &nPoly) != 1) ||
		   (nPoly < 0))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || (sscanf(valS, "%d", &valI) != 1) ||
		   (valI < 0) || ((valI / nPoly) != 4))
		{
		  /* Can only use triangles ((valI / nPoly) != (1 + 3)). */
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
                int	polyBuf[3];

		pIdx = 0;
		while((errNum == WLZ_ERR_NONE) && (pIdx < nPoly))
		{
		  if((fscanf(fP, "%d", &valI) != 1) || (valI != 3) ||
		     (fscanf(fP, "%d %d %d",
			     polyBuf + 0, polyBuf + 1, polyBuf + 2) != 3) ||
		     (*(polyBuf + 0) < 0) || (*(polyBuf + 0) >= sumPoints) ||
		     (*(polyBuf + 1) < 0) || (*(polyBuf + 1) >= sumPoints) ||
		     (*(polyBuf + 2) < 0) || (*(polyBuf + 2) >= sumPoints))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
                    WlzDVertex3	triBuf[3];

		    /* Add triangle to GM. */
		    *(triBuf + 0) = *(pointBuf + *(polyBuf + 0));
		    *(triBuf + 1) = *(pointBuf + *(polyBuf + 1));
		    *(triBuf + 2) = *(pointBuf + *(polyBuf + 2));
		    errNum = WlzGMModelConstructSimplex3D(model, triBuf);
		    ++pIdx;
		  }
		}
	      }
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_TRIANGLE_STRIPS:
	      /* Read triangle strip polydata. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_VERTICIES: /* FALLTHROUGH */
	      /* Read triangle strip polydata. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_LINES:     /* FALLTHROUGH */
	      if(obj == NULL)
	      {
		vHTSz = (sumPoints < 1024)? 1024: sumPoints / 2;
		model = WlzGMModelNew(WLZ_GMMOD_2D, 0, vHTSz, &errNum);
		if(errNum == WLZ_ERR_NONE)
		{
		  dom.ctr = WlzMakeContour(&errNum);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    dom.ctr->model = WlzAssignGMModel(model, NULL);
	            obj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL,
		                      &errNum);
		  }
		}
	      }
	      else if(model && (model->type == WLZ_GMMOD_2D))
	      {
	        vHTSz = (sumPoints < 1024)? 1024: sumPoints / 2;
		if(model->vertexHTSz < vHTSz)
		{
		  errNum = WlzGMModelRehashVHT(model, vHTSz);
		}
	      }
	      else
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || (sscanf(valS, "%d", &nLine) != 1) ||
		   (nLine < 0))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || (sscanf(valS, "%d", &valI) != 1) ||
		   (valI < 0) || ((valI / nLine) != 3))
	        {
		  /* Can only use line formed by two verticies at it's
		   * ends ((valI / nLine) != (1 + 2)). */
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
                int	polyBuf[2];

	        pIdx = 0;
		while((errNum == WLZ_ERR_NONE) && (pIdx < nLine))
		{
		  if((fscanf(fP, "%d", &valI) != 1) || (valI != 2) ||
		     (fscanf(fP, "%d %d",
		             polyBuf + 0, polyBuf + 1) != 2) ||
		     (*(polyBuf + 0) < 0) || (*(polyBuf + 0) >= sumPoints) ||
		     (*(polyBuf + 1) < 0) || (*(polyBuf + 1) >= sumPoints))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
                    WlzDVertex2	linBuf[3];

		    /* Add line to GM. */
		    (linBuf + 0)->vtX = (pointBuf + *(polyBuf + 0))->vtX;
		    (linBuf + 0)->vtY = (pointBuf + *(polyBuf + 0))->vtY;
		    (linBuf + 1)->vtX = (pointBuf + *(polyBuf + 1))->vtX;
		    (linBuf + 1)->vtY = (pointBuf + *(polyBuf + 1))->vtY;
		    errNum = WlzGMModelConstructSimplex2D(model, linBuf);
		    ++pIdx;
		  }
		}
	      }
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_POINT_DATA:
	      if(obj || nPointData)
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {

		valS = strtok(NULL, " \t\n\r\f\v");
		if((valS == NULL) || (sscanf(valS, "%d", &nPointData) != 1) ||
		   (nPointData != sumPoints))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_VECTORS:
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_SCALARS:
	      if(nPointData <= 0)
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		if((strtok(NULL, " \t\n\r\f\v") == NULL) ||
		   ((valS = strtok(NULL, " \t\n\r\f\v")) == NULL) ||
		   (strcmp(valS, "float") && strcmp(valS, "char")))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_LOOKUP_TABLE:
	      if((nPointData <= 0) |
	         ((valS = strtok(NULL, " \t\n\r\f\v")) == NULL) ||
	         strcmp(valS, "default"))
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
		WlzGreyType vType = WLZ_GREY_UBYTE;

		if((pointDataBuf = (float *)
		    AlcMalloc(sizeof(float) * nPointData)) == NULL)
		{
		  errNum = WLZ_ERR_MEM_ALLOC;
		}
		else
		{
		  float pVal;

		  pIdx = 0;
		  while((errNum == WLZ_ERR_NONE) && (pIdx < nPointData))
		  {
		    if(fscanf(fP, "%g", &pVal) != 1)
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      int iVal;

		      iVal = WLZ_NINT(pVal);
		      if(vType != WLZ_GREY_FLOAT)
		      {
			if(fabs(pVal - iVal) > FLT_EPSILON)
			{
			  vType = WLZ_GREY_FLOAT;
			}
			else if((unsigned int )iVal > 255)
			{
			  vType = WLZ_GREY_INT;
			}
		      }
		      pointDataBuf[pIdx] = pVal;
		    }
		    ++pIdx;
		  }
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  WlzVertexP pBuf;

		  pBuf.d3 = pointBuf;
		  dom.pts = WlzMakePoints(WLZ_POINTS_3D, nPointData,  pBuf,
					  nPointData, &errNum);
		  if(errNum == WLZ_ERR_NONE)
		  {
		    val.pts = WlzMakePointValues(nPointData, 0, NULL, vType,
						 &errNum);
		  }
		  if(errNum == WLZ_ERR_NONE)
		  {
		    obj = WlzMakeMain(WLZ_POINTS, dom, val, NULL, NULL,
		                      &errNum);
		  }
		}
		if(errNum == WLZ_ERR_NONE)
		{
		  switch(vType)
		  {
		    case WLZ_GREY_INT:
		      for(pIdx = 0; pIdx < nPointData; ++pIdx)
		      {
			val.pts->values.inp[pIdx] = 
			    WLZ_NINT(pointDataBuf[pIdx]);
		      }
		      break;
		    case WLZ_GREY_UBYTE:
		      for(pIdx = 0; pIdx < nPointData; ++pIdx)
		      {
			val.pts->values.ubp[pIdx] = 
			    WLZ_NINT(pointDataBuf[pIdx]);
		      }
		      break;
		    case WLZ_GREY_FLOAT:
		      for(pIdx = 0; pIdx < nPointData; ++pIdx)
		      {
			val.pts->values.flp[pIdx] = pointDataBuf[pIdx];
		      }
		      break;
		    default:
		      errNum = WLZ_ERR_GREY_TYPE;
		      break;
		  }
		}
	      }
	      break;
	    default:
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	      break;
	  }
	}
      }
    }
    while((endOfData == 0) && (errNum == WLZ_ERR_NONE));
  }
  if((errNum == WLZ_ERR_NONE) && (obj == NULL) && (sumPoints > 0))
  {
    WlzVertexP pBuf;

    pBuf.d3 = pointBuf;
    dom.pts = WlzMakePoints(WLZ_POINTS_3D, sumPoints,  pBuf, sumPoints,
                            &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      obj = WlzMakeMain(WLZ_POINTS, dom, val, NULL, NULL, &errNum);
    }
  }
  AlcFree(pointBuf);
  AlcFree(pointDataBuf);
  if(errNum != WLZ_ERR_NONE)
  {
    WlzFreeObj(obj);
    obj = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz domain object from the given stream using
*		the Visualization Toolkit (structured points) file format.
* \todo		WlzEffReadImgVtk() is unimplemented.
* \param	fP			Input file stream.
* \param	header			Header data structure.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadImgVtk(FILE *fP, WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO: Read VTK images. */
  errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}
