#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFVff_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFVff.c
* \author       Bill Hill
* \date         March 1999
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
* 		the Sunvision '.vff' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static int			WlzEffHeadRecReadVff(
				  char *recPtr,
				  int recMax,
				  FILE *fP);
static int			WlzEffHeadRecParseVff(
				  WlzEffVffHeader *header,
				  char *recStr);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the Sunvision
* 		vff (rank 3) file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjVff(FILE *fP, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_READ_INCOMPLETE;
  int		datumSize,
  		dataCount,
  		recC;
  WlzEffVffHeader header;
  WlzIVertex3	origin,
  		size;
  char		recBuf[WLZEFF_VFF_REC_LEN_MAX];
  WlzObject	*obj = NULL;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  void		***data = NULL;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    (void )memset((void *)&header, 0, sizeof(WlzEffVffHeader));
    header.aspect.vtX = 1.0;
    header.aspect.vtY = 1.0;
    header.aspect.vtZ = 1.0;
    do						     /* Parse the VFF header */
    {
      recC = WlzEffHeadRecReadVff(recBuf, WLZEFF_VFF_REC_LEN_MAX, fP);
      if((recC == '\n') && (WlzEffHeadRecParseVff(&header, recBuf) < 0))
      {
        recC = EOF;			  /* Used to flag an error condition */
      }
    } while(header.ncaa && (recC == '\n'));
    if((recC == '\f') && 			   /* Only binary data files */
       header.ncaa &&
       (header.type == WLZEFF_VFF_TYPE_RASTER) &&
       (header.format == WLZEFF_VFF_FORMAT_SLICE) &&
       (header.rank == 3) &&
       (header.bands == 1) &&
       ((header.bits == 8) || (header.bits == 16) || (header.bits == 32)) &&
       (header.size.vtX > 0) &&
       (header.size.vtY > 0) &&
       (header.size.vtZ > 0))
    {
      errNum = WLZ_ERR_NONE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    origin.vtX = WLZ_NINT(header.origin.vtX);
    origin.vtY = WLZ_NINT(header.origin.vtY);
    origin.vtZ = WLZ_NINT(header.origin.vtZ);
    size.vtX = header.size.vtX;
    size.vtY = header.size.vtY;
    size.vtZ = header.size.vtZ;
    dataCount = size.vtX * size.vtY * size.vtZ;
    switch(header.bits)
    {
      case 8:
	gType = WLZ_GREY_UBYTE;
	if((AlcUnchar3Malloc((unsigned char ****)&data, size.vtZ, size.vtY,
			     size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  data = NULL;
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(unsigned char);
	}
	break;
      case 16:
	gType = WLZ_GREY_SHORT;
	if((AlcShort3Malloc((short ****)&data, size.vtZ, size.vtY,
			    size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  data = NULL;
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(short);
	}
	break;
      case 32:
	gType = WLZ_GREY_INT;
	if((AlcInt3Malloc((int ****)&data, size.vtZ, size.vtY,
			  size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  data = NULL;
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(int);
	}
	break;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if(fread(**data, datumSize, dataCount, fP) != dataCount)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj = WlzFromArray3D((void ***)data, size, origin, gType, gType,
			   0.0, 1.0, 0, 1, &errNum);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = header.aspect.vtX;
    obj->domain.p->voxel_size[1] = header.aspect.vtY;
    obj->domain.p->voxel_size[2] = header.aspect.vtZ;
  }
  else
  {
    if(data)
    {
      Alc3Free(data);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	Last char of record, ie
*		 * EOF - end of file or error.
*		 * '\f' or '{' - end of header.
*		 * '\n' - successfuly read record.
* \ingroup	WlzExtFF
* \brief	Reads a single SunVision VFF file format header record.
* \param	recPtr			Record string buffer.
* \param	recMax			Max number of chars in record.
* \param	fP			Input file stream.
*/
static int	WlzEffHeadRecReadVff(char *recPtr, int recMax, FILE *fP)
{
  int		tI0,
  		recC = 0,
  		recIdx = 0;
  
  if(fP && recPtr)
  {
    while(((recC = getc(fP)) != '\n') &&
          (recC != '\f') && (recC != '{') && (recC != EOF) &&
	  (recIdx < recMax))
    {
      *recPtr++ = (char )recC;
    }
    *recPtr = '\0';
    if(recC == '\f')
    {
      if((tI0 = getc(fP)) != '\n')
      {
        (void )ungetc(tI0, fP);
      }
    }
  }
  return(recC);
}

/*!
* \return	Negative on error.
* \ingroup	WlzExtFF
* \brief	Parses the given string for a single SunVision VFF file format
* 		header record.
* \param	header			Header to be filled in.
* \param	recStr			Record string.
*/
static int	WlzEffHeadRecParseVff(WlzEffVffHeader *header, char *recStr)
{
  int		ok = -1;
  unsigned int	field1,
  		recId = (unsigned int )WLZEFF_VFF_REC_NONE;
  char		*token;

  if(header && recStr)
  {
    if((token = strtok(recStr, " =\t")) == NULL)
    {
      ok = 0;
    }
    else
    {
      if((WlzStringMatchValue((int *)&recId, token,
			"ncaa", WLZEFF_VFF_REC_NCAA,
			"type", WLZEFF_VFF_REC_TYPE,
			"format", WLZEFF_VFF_REC_FORMAT,
			"rank", WLZEFF_VFF_REC_RANK,
			"bands", WLZEFF_VFF_REC_BANDS,
			"bits", WLZEFF_VFF_REC_BITS,
			"rawsize", WLZEFF_VFF_REC_RAWSIZE,
			"size", WLZEFF_VFF_REC_SIZE,
			"origin", WLZEFF_VFF_REC_ORIGIN,
			"extent", WLZEFF_VFF_REC_EXTENT,
			"aspect", WLZEFF_VFF_REC_ASPECT,
			NULL) != 0) &&
         (recId != WLZEFF_VFF_REC_NONE))
      {
        switch(recId)
	{
	  case WLZEFF_VFF_REC_NCAA:
	    header->ncaa = 1;
	    ok = recId;
	    break;
	  case WLZEFF_VFF_REC_TYPE:
	    if((token = strtok(NULL, " \t;")) != NULL)
	    {
	      if(WlzStringMatchValue((int *)&field1, token,
			"connectivity", WLZEFF_VFF_TYPE_CONNECTIVITY,
			"include", WLZEFF_VFF_TYPE_INCLUDE,
			"nurbs", WLZEFF_VFF_TYPE_NURBPATCH,
			"raster", WLZEFF_VFF_TYPE_RASTER,
			"verticies", WLZEFF_VFF_TYPE_VERTICIES,
			NULL) != 0)
	      {
	        header->type = (WlzEffVffType )field1;
		ok = recId;
	      }
	    }
	    break;
	  case WLZEFF_VFF_REC_FORMAT:
	    if((token = strtok(NULL, " \t;")) != NULL)
	    {
	      if(WlzStringMatchValue((int *)&field1, token,
			"base", WLZEFF_VFF_FORMAT_BASE,
			"slice", WLZEFF_VFF_FORMAT_SLICE,
			NULL) != 0)
	      {
	        header->format = (WlzEffVffFormat )field1;
		ok = recId;
	      }
	    }
	    break;
	  case WLZEFF_VFF_REC_RANK:
	    if(((token = strtok(NULL, " \t;")) != NULL) &&
	       (sscanf(token, "%d", &(header->rank)) == 1))
	    {
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_BANDS:
	    if(((token = strtok(NULL, " \t;")) != NULL) &&
	       (sscanf(token, "%d", &(header->bands)) == 1))
	    {
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_BITS:
	    if(((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%d", &(header->bits)) == 1))
	    {
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_RAWSIZE:
	    if(((token = strtok(NULL, " \t;")) != NULL) &&
	       (sscanf(token, "%d", &(header->rawsize)) == 1))
	    {
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_SIZE:
	    if(((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%d", &(header->size.vtX)) == 1) &&
	       ((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%d", &(header->size.vtY)) == 1))
	    {
	      if(((token = strtok(NULL, " \t;")) == NULL) ||
	         (sscanf(token, "%d", &(header->size.vtZ)) != 1))
	      {
	        header->size.vtZ = 1;
	      }
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_ORIGIN:
	    if(((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->origin.vtX)) == 1) &&
	       ((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->origin.vtY)) == 1))
	    {
	      if(((token = strtok(NULL, " \t;")) == NULL) ||
	         (sscanf(token, "%lg", &(header->origin.vtZ)) != 1))
	      {
	        header->size.vtZ = 1;
	      }
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_EXTENT:
	    if(((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->extent.vtX)) == 1) &&
	       ((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->extent.vtY)) == 1))
	    {
	      if(((token = strtok(NULL, " \t;")) == NULL) ||
	         (sscanf(token, "%lg", &(header->extent.vtZ)) != 1))
	      {
	        header->size.vtZ = 1;
	      }
	      ok = recId;
	    }
	    break;
	  case WLZEFF_VFF_REC_ASPECT:
	    if(((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->aspect.vtX)) == 1) &&
	       ((token = strtok(NULL, " \t,;")) != NULL) &&
	       (sscanf(token, "%lg", &(header->aspect.vtY)) == 1))
	    {
	      if(((token = strtok(NULL, " \t;")) == NULL) ||
	         (sscanf(token, "%lg", &(header->aspect.vtZ)) != 1))
	      {
	        header->size.vtZ = 1;
	      }
	      ok = recId;
	    }
	    break;
	}
      }
    }

  }
  return(ok);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Writes the given Woolz object to the given stream using the
* 		Sunvision vff (rank 3) file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjVff(FILE *fP, WlzObject *obj)
{
  unsigned char	***data = NULL;
  int		byteCount;
  WlzIVertex3	origin,
  		size;
  WlzErrorNum	errNum = WLZ_ERR_WRITE_INCOMPLETE;
  WlzEffVffHeader header;

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
    header.origin.vtX = origin.vtX;
    header.origin.vtY = origin.vtY;
    header.origin.vtZ = origin.vtZ;
    header.size.vtX = size.vtX;
    header.size.vtY = size.vtY;
    header.size.vtZ = size.vtZ;
    header.aspect.vtX = obj->domain.p->voxel_size[0];
    header.aspect.vtY = obj->domain.p->voxel_size[1];
    header.aspect.vtZ = obj->domain.p->voxel_size[2];
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
    header.rank = 3;
    header.bands = 1;
    header.bits = 8;
    header.rawsize = size.vtX * size.vtY * size.vtZ;
    header.extent.vtX = size.vtX;
    header.extent.vtY = size.vtY;
    header.extent.vtZ = size.vtZ;
    if((byteCount = fprintf(fP,
			    "ncaa\n"
			    "type=raster;\n"
			    "format=slice;\n"
			    "rank=%d;\n"
			    "bands=%d;\n"
			    "bits=%d;\n"
			    "rawsize=%d;\n"
			    "size=%d %d %d;\n"
			    "origin=%f %f %f;\n"
			    "extent=%f %f %f;\n"
			    "aspect=%f %f %f;\n",
			    header.rank,
			    header.bands,
			    header.bits,
			    sizeof(char) * header.rawsize,
			    header.size.vtX,
			    header.size.vtY,
			    header.size.vtZ,
			    header.origin.vtX,
			    header.origin.vtY,
			    header.origin.vtZ,
			    header.extent.vtX,
			    header.extent.vtY,
			    header.extent.vtZ,
			    header.aspect.vtX,
			    header.aspect.vtY,
			    header.aspect.vtZ)) > 0)
    {
      while((byteCount % 8) != 6)   /* Align VFF data on dbl word boundary */
      {
	byteCount += fprintf(fP, " ");
      }
      byteCount += fprintf(fP, "\f\n");
      if(fwrite(**data, sizeof(unsigned char),
		header.rawsize, fP) != header.rawsize)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;;
      }
    }
  }
  if(data)
  {
    (void )AlcUnchar3Free(data);
  }
  return(errNum);
}
