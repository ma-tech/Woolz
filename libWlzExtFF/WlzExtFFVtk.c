#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFVtk.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for reading and writting Woolz objects to
*		and from the VTK '.vtk' data format.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 10-01-01 bill Add WlzEffHeadReadVtk(), WlzEffReadCtrVtk(),
*		WlzEffReadGMVtk() nad WlzEffReadImgVtk().
* 16-08-00 bill	Add WlzEffWriteCtrVtk() and WlzEffWriteGMModelVtk().
************************************************************************/
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

/************************************************************************
* Function:	WlzEffReadObjVtk
* Returns:	WlzObject *:		Object read from file.
* Purpose:	Reads a Woolz object from the given stream using
*		the Visualization Toolkit (structured points) file
*		format.
*		So far ther is no check on the major and minor version
*		numbers since all files seem compatable.
* Global refs:	-
* Parameters:	FILE *fP:		Input file stream.
* 		WlzErrorNum *dstErr:	Destination error number ptr,
*					may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzEffWriteObjVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Writes the given Woolz object to the given stream
*		using the Visualization Toolkit (structured points)
*		file format.
* Global refs:	-
* Parameters:	FILE *fP:		Output file stream.
*		WlzObject *obj:		Given woolz object.
************************************************************************/
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

/************************************************************************
* Function:	WlzEffWriteImgVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Writes the given Woolz 3D domain object with values
*		object to the given stream using the Visualization
*		Toolkit (structured points) file format.
* Global refs:	-
* Parameters:	FILE *fP:		Output file stream.
*		WlzObject *obj:		Given woolz object.
************************************************************************/
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

/************************************************************************
* Function:	WlzEffWriteCtrVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Writes the given Woolz contour (2D or 3D model) to the
*		given stream using the Visualization Toolkit (polydata)
*		file format.
* Global refs:	-
* Parameters:	FILE *fP:		Output file stream.
*		WlzContour *ctr:	Given woolz contour.
************************************************************************/
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

/************************************************************************
* Function:	WlzEffWriteGMModelVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Writes the given Woolz gemetric model (2D or 3D) to the
*		given stream using the Visualization Toolkit (polydata)
*		file format.
* Global refs:	-
* Parameters:	FILE *fP:		Output file stream.
*		WlzGMModel *model:	Given gemetric model.
************************************************************************/
static WlzErrorNum WlzEffWriteGMModelVtk(FILE *fP, WlzGMModel *model)
{

  int		idI,
  		iCnt;
  int		bufI[3];
  AlcVector	*vec;
  WlzGMEdgeT	*tET;
  WlzGMElemP	eP;
  WlzGMResIdxTb	*resIdxTb = NULL;
  WlzDVertex3	vtx;
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
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
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
    /* Output the vertex indicies of the simplicies. */
    switch(model->type)
    {
      case WLZ_GMMOD_2I:
      case WLZ_GMMOD_2D:
	/* TODO Not sure how to output a 2D model for VTK. */
        break;
      case WLZ_GMMOD_3I:
      case WLZ_GMMOD_3D:
	(void )fprintf(fP,
		       "POLYGONS %d %d\n",
		       model->res.loop.numElm, 4 * model->res.loop.numElm);
	idI = 0;
	vec = model->res.loop.vec;
	iCnt = model->res.loop.numIdx;
	while((errNum == WLZ_ERR_NONE) && (iCnt-- > 0))
	{
	  eP.loop = (WlzGMLoop *)AlcVectorItemGet(vec, idI++);
	  if(eP.loop->idx >= 0)
	  {
	    /* Loop IS a triangle, in 3D nothing else is allowed. */
	    tET = eP.loop->loopT->edgeT;
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
  if(resIdxTb)
  {
    WlzGMModelResIdxFree(resIdxTb);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffHeadReadVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Reads the VTK file header which is
*		  ident\n
*		  title\n
*		  data type\n
*		  type\n
*		where
*		  ident = "vtk DataFile Version MAJOR.MINOR"
*		  title = less than 256 chars followed by a \n
*		  data type = "ASCII"|"BINARY"
*		  type = DATASET "STRUCTURED_POINTS"|"STRUCTURED_GRID"|
*		         "UNSTRUCTURED_GRID"|"POLYDATA"|
*		         "RECTILNEAR_GRID"
* Global refs:	-
* Parameters:	WlzEffVtkHeader *header: Given header data structure
*					be filled.
*		FILE *fP:		Output file stream.
************************************************************************/
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
      valS = strtok(buf, " \t\n");
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
      valS = strtok(buf, " \t\n");
      if((valS == NULL) || strcmp(valS, "DATASET"))
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
        valS = strtok(NULL, " \t\n");
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

/************************************************************************
* Function:	WlzEffReadCtrVtk
* Returns:	WlzObject *:		Object read from file.
* Purpose:	Reads a Woolz contour object from the given stream
*		using the Visualization Toolkit (polydata) file
*		format.
* Global refs:	-
* Parameters:	FILE *fP:		Input file stream.
*		WlzEffVtkHeader *header: Header data structure.
* 		WlzErrorNum *dstErr:	Destination error number ptr,
*					may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	WlzEffReadGMVtk
* Returns:	WlzGMModel *:		Object read from file.
* Purpose:	Reads a WlzGMModel from the given stream
*		using the Visualization Toolkit (polydata) file
*		format.
* Global refs:	-
* Parameters:	FILE *fP:		Input file stream.
*		WlzEffVtkHeader *header: Header data structure.
* 		WlzErrorNum *dstErr:	Destination error number ptr,
*					may be NULL.
************************************************************************/
WlzGMModel	*WlzEffReadGMVtk(FILE *fP, WlzEffVtkHeader *header,
				 WlzErrorNum *dstErr)
{
  int		valI,
		pIdx,
		nPoly,
		nPoints,
		sumPoints,
		endOfData,
		vHTSz;
  char		*valS;
  WlzGMModel	*model = NULL;
  WlzObject	*obj = NULL;
  WlzEffVtkPolyDataType prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDVertex3	*pointBuf = NULL;
  WlzDVertex3	triBuf[3];
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
	valS = strtok(buf, " \t\n");
      }
      if(valS)
      {
	if((WlzStringMatchValue(&valI, valS,
		"POINTS", WLZEFF_VTK_POLYDATATYPE_POINTS,
		"VERTICIES", WLZEFF_VTK_POLYDATATYPE_VERTICIES,
		"LINES", WLZEFF_VTK_POLYDATATYPE_LINES,
		"POLYGONS", WLZEFF_VTK_POLYDATATYPE_POLYGONS,
		"TRIANGLE_STRIPS", WLZEFF_VTK_POLYDATATYPE_TRIANGLE_STRIPS,
		NULL) == 0))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  prim = valI;
	  switch(prim)
	  {
	    case WLZEFF_VTK_POLYDATATYPE_POINTS:
	      valS = strtok(NULL, " \t\n");
	      if((valS == NULL) || (sscanf(valS, "%d", &nPoints) != 1) ||
		 (nPoints <= 0))
	      {
	        errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
	        valS = strtok(NULL, " \t\n");
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
		     (fscanf(fP, "%lg", &((pointBuf + pIdx++)->vtZ)) != 1))
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
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
		valS = strtok(NULL, " \t\n");
		if((valS == NULL) || (sscanf(valS, "%d", &nPoly) != 1) ||
		   (nPoly < 0))
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
	      }
	      if(errNum == WLZ_ERR_NONE)
	      {
	        valS = strtok(NULL, " \t\n");
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
	    case WLZEFF_VTK_POLYDATATYPE_LINES:     /* FALLTHROUGH */
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

/************************************************************************
* Function:	WlzEffReadImgVtk
* Returns:	WlzObject *:		Object read from file.
* Purpose:	Reads a Woolz domain object from the given stream
*		using the Visualization Toolkit (structured points)
*		file format.
* Global refs:	-
* Parameters:	FILE *fP:		Input file stream.
*		WlzEffVtkHeader *header: Header data structure.
* 		WlzErrorNum *dstErr:	Destination error number ptr,
*					may be NULL.
************************************************************************/
WlzObject	*WlzEffReadImgVtk(FILE *fP, WlzEffVtkHeader *header,
				  WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* TODO Read VTK images. */
  errNum = WLZ_ERR_READ_INCOMPLETE;

  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}
