#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlzExtFF/WlzExtFFSlc.c
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
* \brief	Functions for reading and writting Woolz objects to
*		and from the '.slc' data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzErrorNum 		WlzEffHeadReadSlc(
				  WlzEffSlcHeader *header,
				  FILE *fP);
static WlzErrorNum		WlzEffHeadWriteSlc(
				  WlzEffSlcHeader *header,
				  FILE *fP);
static WlzErrorNum		WlzEffSlcDecode8(
				  unsigned char *outBuf,
			   	  unsigned char *inBuf,
				  int outSize);
static WlzErrorNum	        WlzEffSlcEncode8(
				  unsigned char **outBuffer,
				  int *outBufSz,
				  int *outBytes,
				  unsigned char *inBuffer,
				  int inBytes);
static void			WlzEffSlcCreateIcon(
				  WlzEffSlcHeader *header,
				  unsigned char ***data);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the '.slc'
*		file format.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjSlc(FILE *fP, WlzErrorNum *dstErr)
{
  int		idxZ,
  		compPlaneSz,
  		planeBufSz,
  		planeSz;
  unsigned char	***data = NULL;
  unsigned char	*planeBuf = NULL,
		*planePtr;
  WlzIVertex3	origin;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_READ_INCOMPLETE;
  WlzEffSlcHeader header;

  header.icon = NULL;
  origin.vtX = origin.vtY = origin.vtZ = 0;
  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    (void )memset((void *)&header, 0, sizeof(WlzEffSlcHeader));
    errNum = WlzEffHeadReadSlc(&header, fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    planeBufSz = planeSz = header.size.vtX * header.size.vtY;
    if((planeBuf = AlcMalloc(planeSz * sizeof(unsigned char))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((AlcUnchar3Malloc(&data, header.size.vtZ, header.size.vtY,
			 header.size.vtX) != ALC_ER_NONE) || (data == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idxZ = 0;
    while((idxZ < header.size.vtZ) && (errNum == WLZ_ERR_NONE))
    {
      /* Read in each slice into the buffer */
      planePtr = **(data + idxZ);
      switch(header.compression)
      {
	case 0:
	  if(fread(planePtr, sizeof(unsigned char), planeSz, fP) != planeSz)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  break;
	case 1:
	  if(fscanf(fP, "%d X", &compPlaneSz) != 1)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  else
	  {
	    if(compPlaneSz > planeBufSz)
	    {
	      planeBufSz = compPlaneSz;
	      if((planeBuf = AlcRealloc(planeBuf,
		  planeSz * sizeof(unsigned char))) == NULL)
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    if(fread(planeBuf, 1, compPlaneSz, fP) != compPlaneSz)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      errNum = WlzEffSlcDecode8(planePtr, planeBuf, planeSz);
	    }
	  }
	  break;
	  default:
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
      }
      ++idxZ;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzFromArray3D((void ***)data, header.size, origin,
    			 WLZ_GREY_UBYTE, WLZ_GREY_UBYTE,
			 0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = header.spacing.vtX;
    obj->domain.p->voxel_size[1] = header.spacing.vtY;
    obj->domain.p->voxel_size[2] = header.spacing.vtZ;
  }
  if(header.icon)
  {
    AlcFree(header.icon);
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
* 		'.slc' file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjSlc(FILE *fP, WlzObject *obj)
{
  int		idxZ,
  		planeSz,
		compPlaneSz,
		compBufSz;
  WlzErrorNum	errNum = WLZ_ERR_WRITE_INCOMPLETE;
  WlzIVertex3	origin;
  unsigned char	*compBuf = NULL;
  unsigned char ***data = NULL;
  WlzEffSlcHeader header;

  header.icon = NULL;
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
    header.size.vtX = obj->domain.p->lastkl - origin.vtX + 1;
    header.size.vtY = obj->domain.p->lastln - origin.vtY + 1;
    header.size.vtZ = obj->domain.p->lastpl - origin.vtZ + 1;
    header.spacing.vtX = obj->domain.p->voxel_size[0];
    header.spacing.vtY = obj->domain.p->voxel_size[1];
    header.spacing.vtZ = obj->domain.p->voxel_size[2];
    if((header.size.vtX <= 0) || (header.size.vtY <= 0) ||
       (header.size.vtZ <= 0))
    {
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else
    {
      errNum = WlzToArray3D((void ****)&data, obj, header.size, origin,
      			    0, WLZ_GREY_UBYTE);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    header.magic = WLZEFF_SLC_MAGIC;
    header.bits = 8;
    header.units = WLZEFF_SLC_DATAUNITS_MICROMETER;
    header.source = WLZEFF_SLC_DATASRC_OTHER;
    header.modification = WLZEFF_SLC_DATAMOD_OTHER;
    header.compression = 1;
    header.iconSize.vtX = 100;
    header.iconSize.vtY = 100;
    if((header.icon = AlcMalloc(header.iconSize.vtX * header.iconSize.vtY *
    			        sizeof(unsigned char) * 3)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzEffSlcCreateIcon(&header, data);
    errNum = WlzEffHeadWriteSlc(&header, fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    planeSz = header.size.vtX * header.size.vtY;
    compPlaneSz = planeSz + 64;
    compBufSz = compPlaneSz;
    if((compBuf = AlcMalloc(sizeof(unsigned char) * compBufSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    idxZ = 0;
    while((idxZ < header.size.vtZ) && (errNum == WLZ_ERR_NONE))
    {
      /* Encode each plane into the buffer */
      errNum = WlzEffSlcEncode8(&compBuf, &compBufSz, &compPlaneSz,
      				**(data + idxZ), planeSz);
      if(errNum == WLZ_ERR_NONE)
      {
        if((fprintf(fP, "%d X", compPlaneSz) < 3) ||
	   (fwrite(compBuf, sizeof(unsigned char), compPlaneSz,
	           fP) != compPlaneSz))
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
      }
      ++idxZ;
    }
  }
  if(header.icon)
  {
    AlcFree(header.icon);
  }
  if(data)
  {
    Alc3Free((void ***)data);
  }
  if(compBuf)
  {
    AlcFree(compBuf);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Reads an SLC header from the given file stream.
* \param	header			Given header data structure to
*					be filled in.
* \param	fP			Given input file stream.
*/
static WlzErrorNum WlzEffHeadReadSlc(WlzEffSlcHeader *header, FILE *fP)
{
  int		iconCnSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fscanf(fP, "%d", &(header->magic)) != 1) ||
     (header->magic != WLZEFF_SLC_MAGIC))
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  else if((fscanf(fP, "%d", &(header->size.vtX)) != 1) ||
	  (fscanf(fP, "%d", &(header->size.vtY)) != 1) ||
	  (fscanf(fP, "%d", &(header->size.vtZ)) != 1) ||
	  (fscanf(fP, "%d", &(header->bits)) != 1) ||
	  (header->size.vtX <= 0) ||
	  (header->size.vtY <= 0) ||
	  (header->size.vtZ <= 0) ||
	  (header->bits != 8) ||
	  (fscanf(fP, "%g", &(header->spacing.vtX)) != 1) ||
	  (fscanf(fP, "%g", &(header->spacing.vtY)) != 1) ||
	  (fscanf(fP, "%g", &(header->spacing.vtZ)) != 1) ||
	  (fscanf(fP, "%d", &(header->units)) != 1) ||
	  (fscanf(fP, "%d", &(header->source)) != 1) ||
	  (fscanf(fP, "%d", &(header->modification)) != 1) ||
	  (fscanf(fP, "%d\n", &(header->compression)) != 1) ||
	  (fscanf(fP, "%d %d X",
	          &(header->iconSize.vtX), &(header->iconSize.vtY)) != 2) ||
	  (header->iconSize.vtX <= 0) || (header->iconSize.vtY <= 0))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    iconCnSz = header->iconSize.vtX * header->iconSize.vtY;
    if((header->icon = (unsigned char *)AlcMalloc(sizeof(unsigned char) *
		                		  iconCnSz * 3)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((fread(header->icon, 1, iconCnSz, fP) != iconCnSz) ||
       (fread(header->icon + iconCnSz, 1, iconCnSz, fP) != iconCnSz) ||
       (fread(header->icon + (2 * iconCnSz), 1, iconCnSz, fP) != iconCnSz))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Writes an SLC header to the given file stream.
* \param	header			Given header data structure to
*					be written.
* \param	fP			Given output file stream.
*/
static WlzErrorNum WlzEffHeadWriteSlc(WlzEffSlcHeader *header, FILE *fP)
{
  int		iconCnSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  iconCnSz = header->iconSize.vtX * header->iconSize.vtY;
  if(fprintf(fP, "%d\n", header->magic) <= 0)
  {
    errNum = WLZ_ERR_WRITE_EOF;
  }
  else if((fprintf(fP, "%d ", header->size.vtX) <= 0) ||
	  (fprintf(fP, "%d ", header->size.vtY) <= 0) ||
	  (fprintf(fP, "%d ", header->size.vtZ) <= 0) ||
	  (fprintf(fP, "%d\n", header->bits) <= 0) ||
	  (fprintf(fP, "%g ", header->spacing.vtX) <= 0) ||
	  (fprintf(fP, "%g ", header->spacing.vtY) <= 0) ||
	  (fprintf(fP, "%g\n", header->spacing.vtZ) <= 0) ||
	  (fprintf(fP, "%d ", header->units) <= 0) ||
	  (fprintf(fP, "%d ", header->source) <= 0) ||
	  (fprintf(fP, "%d ", header->modification) <= 0) ||
	  (fprintf(fP, "%d\n", header->compression) <= 0) ||
	  (fprintf(fP, "%d %d X",
	           header->iconSize.vtX, header->iconSize.vtY) < 5) ||
          (fwrite(header->icon, 1, iconCnSz * 3, fP) != (3 * iconCnSz)))
  {
    errNum = WLZ_ERR_WRITE_INCOMPLETE;
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Decodes the given input buffer into the output buffer.
* \param	outBuf			Output buffer.
* \param	nBuf			Input buffer to be decoded.
* \param	outSize			Size of output buffer.
*/
static WlzErrorNum WlzEffSlcDecode8(unsigned char *outBuf,
				    unsigned char *inBuf, int outSize)
{
  int		count;
  unsigned char val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  while(((count = (unsigned int )((val = *(inBuf++)) & 0x7f)) != 0) &&
        (errNum == WLZ_ERR_NONE))
  {
    if(count > outSize)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      outSize -= count;
      if(val & 0x80)
      {
	while(count-- > 0)
	{
	  *(outBuf++) = *(inBuf++);
	}
      }
      else
      {
	val = *(inBuf++);
	while(count-- > 0)
	{
	  *(outBuf++) = val;
	}
      }
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzExtFF
* \brief	Encodes the given input buffer into the output buffer.
* \param	outBuffer		Output buffer which may be reallocated
*					using AlcRealloc().
* \param	outBufSz		Source and destination ptr for the size
* 					of the output buffer (ie the number of
* 					bytes allocated).
* \param	outBytes		Destination ptr for the number of bytes
* 					in the output buffer.
* \param	inBuffer		Input buffer to be encoded.
* \param	inBytes			Number of bytes in the input buffer.
*/
static WlzErrorNum WlzEffSlcEncode8(unsigned char **outBuffer, int *outBufSz,
				    int *outBytes, unsigned char *inBuffer,
				    int inBytes)
{
  int 		tI0,
  		processCount,
		inByteCount,
		outByteCount;
  unsigned char	*inBuf0,
  		*inBuf1,
  		*inBufEnd,
		*outBuf;
  unsigned char	val;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  outByteCount = 0;
  outBuf = *outBuffer;
  inBuf1 = inBuffer;
  inBufEnd = inBuffer + inBytes - 1;
  inByteCount = inBytes;
  while((inBuf1 < inBufEnd) && (errNum == WLZ_ERR_NONE))
  {
    inBuf0 = inBuf1;
    inBuf1 += 2;
    while((inBuf1 < inBufEnd) &&
	  ((*(inBuf1 - 2) != *(inBuf1 - 1)) ||
	   (*(inBuf1 - 1) != *(inBuf1))))
    {
      inBuf1++;
    }
    inBuf1 -= 2;
    inByteCount = inBuf1 - inBuf0;
    while(inByteCount && (errNum == WLZ_ERR_NONE))
    {
      processCount = (inByteCount > 127)? 127: inByteCount;
      inByteCount  -= processCount;
      *(outBuf++) = 0x80 | processCount;
      outByteCount += processCount + 1;
      if(outByteCount >= *outBufSz) 
      {
        tI0 = (*outBufSz > inBytes)? *outBufSz: inBytes;
	*outBufSz = tI0 + (tI0 / 4) + 1024;
	tI0 = outBuf - *outBuffer;
	if((*outBuffer = AlcRealloc(*outBuffer, *outBufSz)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  outBuf = *outBuffer + tI0;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	while(processCount--)
	{
	  *(outBuf++) = *(inBuf0++);
	}
      }
    }
    val = *(inBuf1++);
    while((inBuf1 < inBufEnd) && (*inBuf1 == val) && (errNum == WLZ_ERR_NONE))
    {
      inBuf1++;
    }
    inByteCount = inBuf1 - inBuf0;
    while(inByteCount && (errNum == WLZ_ERR_NONE))
    {
      processCount = (inByteCount > 127)? 127: inByteCount;
      inByteCount  -= processCount;
      outByteCount += 2;
      if(outByteCount  >= *outBufSz ) 
      {
        tI0 = (*outBufSz > inBytes)? *outBufSz: inBytes;
	*outBufSz = tI0 + (tI0 / 4) + 1024;
	tI0 = outBuf - *outBuffer;
	if((*outBuffer = AlcRealloc(*outBuffer, *outBufSz)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  outBuf = *outBuffer + tI0;
	}
      }
      else
      {
	*(outBuf++) = processCount;
	*(outBuf++) = val;
      }
    }
  }
  ++outByteCount;
  *(outBuf++) = 0;
  *outBytes = (errNum == WLZ_ERR_NONE)? outByteCount: 0;
  return(errNum);
}

/*!
* \return	void
* \ingroup	WlzExtFF
* \brief	Creates an icon for the SLC header from the given cuboid
*		of data by sampling the centre plane.
* \param	header			SLC header data structure.
* \param	data			The cuboid of data.
*/
static void	WlzEffSlcCreateIcon(WlzEffSlcHeader *header,
				    unsigned char ***data)
{
  int		idxX,
		idxY,
		dstOffY,
		iconSz;
  WlzDVertex2	samStep;
  unsigned char	*dst1D,
  		*src1D;
  unsigned char	**src2D;

  dst1D = header->icon;
  src2D = *(data + (header->size.vtZ / 2));
  samStep.vtX = (double )( header->size.vtX) / (double )(header->iconSize.vtX);
  samStep.vtY = (double )( header->size.vtY) / (double )(header->iconSize.vtY);
  dstOffY = 0;
  for(idxY = 0; idxY < header->iconSize.vtY; ++idxY)
  {
    src1D = *(src2D + (int )(idxY * samStep.vtY));
    for(idxX = 0; idxX < header->iconSize.vtX; ++idxX)
    {
      *(dst1D + dstOffY + idxX) = *(src1D + (int )(idxX * samStep.vtX));
    }
    dstOffY += header->iconSize.vtX;
  }
  iconSz = header->iconSize.vtX * header->iconSize.vtY;
  (void )memcpy(header->icon + iconSz, header->icon, iconSz);
  (void )memcpy(header->icon + (2 * iconSz), header->icon, iconSz);
}
