#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFVtk.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2005 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for reading and writting Woolz objects to and from
* 		the VTK '.vtk' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>
#include <WlzExtFF.h>
#include <string.h>

static WlzErrorNum 		WlzEffWriteImgVtk(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteCtrVtk(
				  FILE *fP,
				  WlzContour *ctr);
static WlzErrorNum 		WlzEffWriteGMModelVtk(
				  FILE *fP,
				  WlzGMModel *model);
static WlzErrorNum		WlzEffHeadReadVtk(
				  WlzEffVtkHeader *header,
				  FILE *fP);
static WlzObject		*WlzEffReadCtrVtk(
				  FILE *fP,
				  WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr);
static WlzGMModel		*WlzEffReadGMVtk(
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
	obj = WlzEffReadCtrVtk(fP, &header, &errNum);
        break;
      case WLZEFF_VTK_TYPE_STRUCTURED_GRID:   /* FALLTHROUGH */
      case WLZEFF_VTK_TYPE_UNSTRUCTURED_GRID: /* FALLTHROUGH */
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
  int		bufCount,
  		bufIdx,
  		dataCount;
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

  int		dim,
  		idI,
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
        dim = 2;
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D: /* FALLTHROUGH */
      case WLZ_GMMOD_3N:
	dim = 3;
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
    if(fgets(buf, 256, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      buf[255] = '\0';
      valS = strtok(buf, " \t\n\r\f\v");
      if((valS == NULL) ||
         (WlzStringMatchValue(&valI, valS,
      			      "ASCII", WLZEFF_VTK_DATATYPE_ASCII,
			      "BINARY", WLZEFF_VTK_DATATYPE_BINARY,
			      NULL) == 0))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
        header->dataType = valI;
      }
    }
  }
  /* Read and parse the type. */
  if(errNum == WLZ_ERR_NONE)
  {
    if(fgets(buf, 256, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      buf[255] = '\0';
      valS = strtok(buf, " \t\n\r\f\v");
      if((valS == NULL) || strcmp(valS, "DATASET"))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
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
	}
	else
	{
          header->type = valI;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz contour object from the given stream using
*		the Visualization Toolkit (polydata) file format.
* \param	fP			Input file stream.
* \param	header			Header data structure.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadCtrVtk(FILE *fP, WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr)
{
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*obj = NULL;
  WlzContour	*ctr = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((ctr = WlzMakeContour(&errNum)) != NULL)
  {
    ctr->model = WlzAssignGMModel(WlzEffReadGMVtk(fP, header, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      dom.ctr = ctr;
      val.core = NULL;
      obj = WlzMakeMain(WLZ_CONTOUR, dom, val, NULL, NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      WlzFreeContour(ctr);
      ctr = NULL;
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
WlzGMModel	*WlzEffReadGMVtk(FILE *fP, WlzEffVtkHeader *header,
				 WlzErrorNum *dstErr)
{
  int		valI,
		pIdx,
		nLine,
		nPoly,
		nPoints,
		sumPoints,
		endOfData,
		vHTSz;
  char *valS;
  WlzGMModel	*model = NULL;
  WlzObject	*obj = NULL;
  WlzEffVtkPolyDataType prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDVertex3	*pointBuf = NULL;
  WlzDVertex3	triBuf[3];
  WlzDVertex2	linBuf[3];
  int		polyBuf[3];
  char		buf[256];

  if(header->dataType == WLZEFF_VTK_DATATYPE_BINARY)
  {
    /* Can only read ascii polydata. */
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sumPoints = 0;
    endOfData = 0;
    do
    {
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
		  errNum == WLZ_ERR_MEM_ALLOC;
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
	      if(model == NULL)
	      {
		vHTSz = (nPoints < 1024)? 1024: nPoints / 3;
    		model = WlzGMModelNew(WLZ_GMMOD_3D, 0, vHTSz, &errNum);
	      }
	      else
	      {
	        vHTSz = (sumPoints < 1024)? 1024: sumPoints / 3;
		if(model->vertexHTSz < vHTSz)
		{
		  errNum = WlzGMModelRehashVHT(model, vHTSz);
		}
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
	        pIdx = 0;
		while((errNum == WLZ_ERR_NONE) && (pIdx < nPoly))
		{
		  if((fscanf(fP, "%d", &valI) != 1) || (valI != 3) ||
		     (fscanf(fP, "%d %d %d",
		             polyBuf + 0, polyBuf + 1, polyBuf + 2) != 3) ||
		     (*(polyBuf + 0) < 0) || (*(polyBuf + 0) >= nPoints) ||
		     (*(polyBuf + 1) < 0) || (*(polyBuf + 1) >= nPoints) ||
		     (*(polyBuf + 2) < 0) || (*(polyBuf + 2) >= nPoints))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
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
	      /* TODO Read triangle strip polydata. */
	      errNum == WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_VERTICIES: /* FALLTHROUGH */
	      /* TODO Read triangle strip polydata. */
	      errNum == WLZ_ERR_READ_INCOMPLETE;
	      break;
	    case WLZEFF_VTK_POLYDATATYPE_LINES:     /* FALLTHROUGH */
	      if(model == NULL)
	      {
		vHTSz = (nPoints < 1024)? 1024: nPoints / 2;
    		model = WlzGMModelNew(WLZ_GMMOD_2D, 0, vHTSz, &errNum);
	      }
	      else
	      {
	        vHTSz = (sumPoints < 1024)? 1024: sumPoints / 2;
		if(model->vertexHTSz < vHTSz)
		{
		  errNum = WlzGMModelRehashVHT(model, vHTSz);
		}
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
	        pIdx = 0;
		while((errNum == WLZ_ERR_NONE) && (pIdx < nLine))
		{
		  if((fscanf(fP, "%d", &valI) != 1) || (valI != 2) ||
		     (fscanf(fP, "%d %d",
		             polyBuf + 0, polyBuf + 1) != 2) ||
		     (*(polyBuf + 0) < 0) || (*(polyBuf + 0) >= nPoints) ||
		     (*(polyBuf + 1) < 0) || (*(polyBuf + 1) >= nPoints))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		  else
		  {
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
	    default:
	      errNum == WLZ_ERR_READ_INCOMPLETE;
	      break;
	  }
	}
      }
    }
    while((endOfData == 0) && (errNum == WLZ_ERR_NONE));
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(model)
    {
      (void )WlzGMModelFree(model);
      model = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(model);
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

  /* TODO Read VTK images. */
  errNum = WLZ_ERR_UNIMPLEMENTED;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}
