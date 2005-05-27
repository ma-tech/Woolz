#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFDen.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Functions for reading and writting Woolz objects to
*		and from the Stanford University '.den' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>
#include <WlzExtFF.h>


/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
*		Stanford density file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjDen(FILE *fP, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_READ_INCOMPLETE;
  WlzIVertex3	origin,
  		size;
  unsigned char	***data = NULL;
  WlzObject	*obj = NULL;
  WlzEffDenHeader header;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if ((fread(&(header.version), sizeof(short), 1, fP) != 1) ||
	   (header.version != WLZEFF_DEN_VERSION) ||
	   (fread(header.orgMin, sizeof(short), 3, fP) != 3) ||
	   (fread(header.orgMax, sizeof(short), 3, fP) != 3) ||
	   (fread(header.orgLen, sizeof(short), 3, fP) != 3) ||
	   (fread(header.extrMin, sizeof(short), 3, fP) != 3) ||
	   (fread(header.extrMax, sizeof(short), 3, fP) != 3) ||
	   (fread(header.extrLen, sizeof(short), 3, fP) != 3) ||
	   (fread(header.mapMin, sizeof(short), 3, fP) != 3) ||
	   (fread(header.mapMax, sizeof(short), 3, fP) != 3) ||
	   (fread(header.mapLen, sizeof(short), 3, fP) != 3) ||
	   (header.mapLen[0] <= 0) ||
	   (header.mapLen[1] <= 0) ||
	   (header.mapLen[2] <= 0) ||
	   (fread(&(header.mapWarps), sizeof(short), 1, fP) != 1) ||
	   (fread(&(header.mapLength), sizeof(unsigned int), 1, fP) != 1) ||
	   (header.mapLength <= 0) ||
	   (header.mapLength != (header.mapLen[0] * header.mapLen[1] *
				 header.mapLen[2])))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else if((AlcUnchar3Malloc(&data, header.mapLen[2],
    			 header.mapLen[1], header.mapLen[0]) != ALC_ER_NONE) ||
          (data == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else if(fread(**data, sizeof(unsigned char),
		header.mapLength, fP) != header.mapLength)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    size.vtX = header.mapLen[0];
    size.vtY = header.mapLen[1];
    size.vtZ = header.mapLen[2];
    origin.vtX = header.mapMin[0];
    origin.vtY = header.mapMin[1];
    origin.vtZ = header.mapMin[2];
    obj = WlzFromArray3D((void ***)data, size, origin,
			 WLZ_GREY_UBYTE, WLZ_GREY_UBYTE,
			 0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = 1.0;
    obj->domain.p->voxel_size[1] = 1.0;
    obj->domain.p->voxel_size[2] = 1.0;
  }
  else
  {
    if(data)
    {
      (void )AlcUnchar3Free(data);
      obj = NULL;
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
* \brief	Writes the given Woolz object to the given stream using
*		the Stanford density file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjDen(FILE *fP, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_WRITE_INCOMPLETE;
  WlzIVertex3	origin,
		size;
  unsigned char	***data = NULL;
  WlzEffDenHeader header;

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
    header.version = WLZEFF_DEN_VERSION;
    header.orgMin[0] = origin.vtX;
    header.orgMin[1] = origin.vtY;
    header.orgMin[2] = origin.vtZ;
    header.orgMax[0] = origin.vtX + size.vtX - 1;
    header.orgMax[1] = origin.vtY + size.vtY - 1;
    header.orgMax[2] = origin.vtZ + size.vtZ - 1;
    header.orgLen[0] = size.vtX;
    header.orgLen[1] = size.vtY;
    header.orgLen[2] = size.vtZ;
    header.extrMin[0] = origin.vtX;
    header.extrMin[1] = origin.vtY;
    header.extrMin[2] = origin.vtZ;
    header.extrMax[0] = origin.vtX + size.vtX - 1;
    header.extrMax[1] = origin.vtY + size.vtY - 1;
    header.extrMax[2] = origin.vtZ + size.vtZ - 1;
    header.extrLen[0] = size.vtX;
    header.extrLen[1] = size.vtY;
    header.extrLen[2] = size.vtZ;
    header.mapMin[0] = origin.vtX;
    header.mapMin[1] = origin.vtY;
    header.mapMin[2] = origin.vtZ;
    header.mapMax[0] = origin.vtX + size.vtX - 1;
    header.mapMax[1] = origin.vtY + size.vtY - 1;
    header.mapMax[2] = origin.vtZ + size.vtZ - 1;
    header.mapLen[0] = size.vtX;
    header.mapLen[1] = size.vtY;
    header.mapLen[2] = size.vtZ;
    header.mapLength = size.vtX * size.vtY * size.vtZ;
    if((fwrite(&(header.version), sizeof(short), 1, fP) == 1) &&
	(fwrite(header.orgMin, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.orgMax, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.orgLen, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.extrMin, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.extrMax, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.extrLen, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.mapMin, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.mapMax, sizeof(short), 3, fP) == 3) &&
	(fwrite(header.mapLen, sizeof(short), 3, fP) == 3) &&
	(fwrite(&(header.mapWarps), sizeof(short), 1, fP) == 1) &&
	(fwrite(&(header.mapLength), sizeof(unsigned int), 1, fP) == 1) &&
	(fwrite(**data, sizeof(unsigned char),
		header.mapLength, fP) == header.mapLength))
    {
      errNum = WLZ_ERR_NONE;
    }
    (void )AlcUnchar3Free(data);
  }
  return(errNum);
}
