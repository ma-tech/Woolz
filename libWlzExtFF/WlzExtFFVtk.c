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
* 16-08-00 bill	Add WlzEffWriteCtrVtk() and WlzEffWriteGMModelVtk().
************************************************************************/
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzErrorNum 		WlzEffWriteVoxVtk(
				  FILE *fP,
				  WlzObject *obj);
static WlzErrorNum 		WlzEffWriteCtrVtk(
				  FILE *fP,
				  WlzContour *ctr);
static WlzErrorNum 		WlzEffWriteGMModelVtk(
				  FILE *fP,
				  WlzGMModel *model);

/************************************************************************
* Function:	WlzEffReadObjVtk
* Returns:	WlzObject *:		Object read from file.
* Purpose:	Reads a Woolz object from the given stream using
*		the Visualization Toolkit (structured points) file
*		format.
* Global refs:	-
* Parameters:	FILE *fP:		Input file stream.
* 		WlzErrorNum *dstErr:	Destination error number ptr,
*					may be NULL.
************************************************************************/
WlzObject	*WlzEffReadObjVtk(FILE *fP, WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_READ_INCOMPLETE;
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
        errNum = WlzEffWriteVoxVtk(fP, obj);
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
* Function:	WlzEffWriteVoxVtk
* Returns:	WlzErrorNum		Woolz error number.
* Purpose:	Writes the given Woolz 3D domain object with values
*		object to the given stream using the Visualization
*		Toolkit (structured points) file format.
* Global refs:	-
* Parameters:	FILE *fP:		Output file stream.
*		WlzObject *obj:		Given woolz object.
************************************************************************/
static WlzErrorNum WlzEffWriteVoxVtk(FILE *fP, WlzObject *obj)
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
	/* TODO */
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
