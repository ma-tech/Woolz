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
************************************************************************/
#include <Wlz.h>
#include <WlzExtFF.h>

/************************************************************************
* Function:	WlzEffReadObjVtk					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a Woolz object from the given stream using 	*
*		the Visualization Toolkit (structured points) file	*
*		format.							*
* Global refs:	-							*
* Parameters:	FILE *fP:		Input file stream.		*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
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
* Function:	WlzEffWriteObjVtk					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given Woolz object to the given stream 	*
*		using the Visualization Toolkit (structured points)	*
*		file format.						*
* Global refs:	-							*
* Parameters:	FILE *fP:		Output file stream.		*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
WlzErrorNum	WlzEffWriteObjVtk(FILE *fP, WlzObject *obj)
{
  unsigned char	***data = NULL;
  WlzErrorNum	errNum = WLZ_ERR_WRITE_INCOMPLETE;
  int		bufCount,
  		bufIdx,
  		dataCount;
  WlzFVertex3	aspect;
  WlzIVertex3	origin,
  		size;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->type != WLZ_3D_DOMAINOBJ)
  {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
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
