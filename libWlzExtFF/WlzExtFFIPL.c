#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFIPL_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFIPL.c
* \author       Richard Baldock
* \date         November 1999
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
* \brief	Functions for reading and writing IPLab image formats.
* \ingroup	WlzExtFF
*/

#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzErrorNum WlzEffHeadReadIPL(WlzEffIPLHeader *header, FILE *fP),
		WlzEffHeadWriteIPL(WlzEffIPLHeader *header, FILE *fP);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the '.ipl'
* 		file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjIPL(FILE *fP, WlzErrorNum *dstErr)
{
  int		idxZ, i,
  		planeSz, pixelSz;
  void 		***data;
  unsigned char	***ubData = NULL;
  short		***shData = NULL;
  int		***inData = NULL;
  float		***flData = NULL;
  unsigned char	*planeBuf = NULL,
		*planePtr;
  WlzIVertex3	origin, size;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_READ_INCOMPLETE;
  WlzEffIPLHeader header;
  WlzGreyType	gType;

  origin.vtX = origin.vtY = origin.vtZ = 0;
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    (void )memset((void *)&header, 0, sizeof(WlzEffIPLHeader));
    errNum = WlzEffHeadReadIPL(&header, fP);
    size.vtX = header.nWidth;
    size.vtY = header.nHeight;
    size.vtZ = header.nFrames;
  }

  if(errNum == WLZ_ERR_NONE)
  {
    switch( header.dType ){
    case WLZEFF_IPL_TYPE_UBYTE:
      gType = WLZ_GREY_UBYTE;
      pixelSz = sizeof(char);
      if((AlcUnchar3Malloc(&ubData, size.vtZ, size.vtY,
			   size.vtX) != ALC_ER_NONE) || (ubData == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      data = (void ***) ubData;
      break;

    case WLZEFF_IPL_TYPE_SHORT:
    case WLZEFF_IPL_TYPE_U_16:
    case WLZEFF_IPL_TYPE_COL_16:
      gType = WLZ_GREY_SHORT;
      pixelSz = sizeof(short);
      if((AlcShort3Malloc(&shData, size.vtZ, size.vtY,
			  size.vtX) != ALC_ER_NONE) || (shData == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      data = (void ***) shData;
      break;

    case WLZEFF_IPL_TYPE_INT:
    case WLZEFF_IPL_TYPE_COL_24:
      gType = WLZ_GREY_INT;
      pixelSz = sizeof(int);
      if((AlcInt3Malloc(&inData, size.vtZ, size.vtY,
			size.vtX) != ALC_ER_NONE) || (inData == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      data = (void ***) inData;
      break;

    case WLZEFF_IPL_TYPE_FLOAT:
      gType = WLZ_GREY_FLOAT;
      pixelSz = sizeof(float);
      if((AlcFloat3Malloc(&flData, size.vtZ, size.vtY,
			  size.vtX) != ALC_ER_NONE) || (flData == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
      data = (void ***) flData;
      break;

    default:
      errNum = WLZ_ERR_PARAM_TYPE;
      break;
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    planeSz = size.vtX * size.vtY;
    idxZ = 0;
    while((idxZ < size.vtZ) && (errNum == WLZ_ERR_NONE))
    {
      /* Read in each slice into the buffer */
      planePtr = (unsigned char *) (**(data + idxZ));

      switch( header.dType ){
      case WLZEFF_IPL_TYPE_UBYTE:
      case WLZEFF_IPL_TYPE_SHORT:
      case WLZEFF_IPL_TYPE_U_16:
      case WLZEFF_IPL_TYPE_COL_16:
      case WLZEFF_IPL_TYPE_INT:
      case WLZEFF_IPL_TYPE_FLOAT:
	if(fread(planePtr, pixelSz, planeSz, fP) != planeSz)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	break;

      case WLZEFF_IPL_TYPE_COL_24:
	for(i=0; (i < planeSz) && (errNum == WLZ_ERR_NONE); i++, planePtr++){
	  if(fread(planePtr, sizeof(char)*3, 1, fP) != 1)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  planePtr[3] = 0;
	}
	break;
      default:
        break;
      }

      ++idxZ;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzFromArray3D((void ***)data, size, origin,
    			 gType, gType,
			 0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = 1.0;
    obj->domain.p->voxel_size[1] = 1.0;
    obj->domain.p->voxel_size[2] = 1.0;
  }

  if(planeBuf)
  {
    AlcFree(planeBuf);
  }

  if(errNum != WLZ_ERR_NONE)
  {
    if(obj)
    {
      WlzFreeObj(obj);
    }
    else if(data)
    {
      Alc3Free((void ***)data);
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
* 		'.ipl' file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjIPL(FILE *fP, WlzObject *obj)
{
  WlzIVertex3   size;
  WlzErrorNum	errNum = WLZ_ERR_WRITE_INCOMPLETE;
  WlzIVertex3	origin;
  unsigned char ***data = NULL;
  WlzEffIPLHeader header;

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

    if((size.vtX <= 0) || (size.vtY <= 0) ||
       (size.vtZ <= 0))
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
    errNum = WlzEffHeadWriteIPL(&header, fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* write the data */
  }
  if(data)
  {
    Alc3Free((void ***)data);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an IPL header from the given file stream.
* \param	header			Given header data structure to be
*					filled in.
* \param	fP			Given input file stream.
*/
static WlzErrorNum WlzEffHeadReadIPL(WlzEffIPLHeader *header, FILE *fP)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  unsigned char	ubBuf[16];

  /* version, format and datatype */
  if( fread(&(ubBuf[0]), sizeof(char), 6, fP) != 6 ){
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else {
    /* should test this version */
    strncpy(&(header->version[0]), (const char *) &(ubBuf[0]), 4);
    if( (header->format = ubBuf[4]) != 0 ){
      errNum = WLZ_ERR_PARAM_TYPE;
    }
    switch( (int) ubBuf[5] ){
    case 0:
      header->dType = WLZEFF_IPL_TYPE_UBYTE;
      break;
    case 1:
      header->dType = WLZEFF_IPL_TYPE_SHORT;
      break;
    case 2:
      header->dType = WLZEFF_IPL_TYPE_INT;
      break;
    case 3:
      header->dType = WLZEFF_IPL_TYPE_FLOAT;
      break;
    case 4:
      header->dType = WLZEFF_IPL_TYPE_COL_16;
      break;
    case 5:
      header->dType = WLZEFF_IPL_TYPE_COL_24;
      break;
    case 6:
      header->dType = WLZEFF_IPL_TYPE_U_16;
      break;
    default:
      errNum = WLZ_ERR_PARAM_TYPE;
      break;
    }
  }

  /* height and width */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->nWidth = *((int *) &(ubBuf[0]));
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->nHeight = *((int *) &(ubBuf[0]));
  }

  /* various id's and modes */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 6, fP) != 6 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->fileClutID = ubBuf[0];
    header->overlayInFile = ubBuf[1];
    header->viewMode = ubBuf[4];
  }

  /* number of sections */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 2, fP) != 2 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->nFrames = *((short *) &(ubBuf[0]));
  }

  /* delta? plus reserved space */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 8, fP) != 8 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->delta = *((double *) &(ubBuf[0]));
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* units string */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 10, fP) != 10 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    strncpy(&(header->units[0]), (const char *) &(ubBuf[0]), 10);
  }

  /* normalisation data */
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->normType = ubBuf[1];
    header->normSource = ubBuf[2];
    header->numRegNarks = ubBuf[3];
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 8, fP) != 8 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->normMin = *((double *) &(ubBuf[0]));
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 8, fP) != 8 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    header->normMax = *((double *) &(ubBuf[0]));
  }
  if( errNum == WLZ_ERR_NONE ){
    if( fread(&(ubBuf[0]), sizeof(char), 4, fP) != 4 ){
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }

  /* now read the color table */
  if( errNum == WLZ_ERR_NONE ){
    if( (header->colorTable = (void *) AlcCalloc(2048, sizeof(char))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      if( fread(header->colorTable, sizeof(char), 2048, fP) != 2048 ){
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
  }

  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes an IPL header to the given file stream.
* \todo		Just a stub, still to be written.
* \param	header			Given header data structure to
*					be written.
* \param	fP			Given output file stream.
*/
static WlzErrorNum WlzEffHeadWriteIPL(WlzEffIPLHeader *header, FILE *fP)
{
  WlzErrorNum	errNum = WLZ_ERR_UNIMPLEMENTED;

  return(errNum);
}

