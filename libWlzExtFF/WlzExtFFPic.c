#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFPic.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for reading and writting Woolz objects to
*		and from the Biorad '.pic' data format.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>
#include <WlzExtFF.h>

#ifdef SUNOS5
#ifndef __x86
extern void swab(const void *, void *, size_t);
#endif /* __x86 */
#endif /* SUNOS5 */

/************************************************************************
* Function:	WlzEffReadObjPic					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a Woolz object from the given stream using 	*
*		the BioRad Confocal '.pic' format.			*
* Global refs:	-							*
* Parameters:	FILE *fP:		Input file stream.		*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
WlzObject 	*WlzEffReadObjPic(FILE *fP, WlzErrorNum *dstErr)
{
  int		tI0,
  		byteData,
		magic;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIVertex3	origin,
  		size;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzObject	*obj = NULL;
  void		***data = NULL;
  unsigned char	header[WLZEFF_PIC_HEADBYTES];

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(fread(header, sizeof(char), WLZEFF_PIC_HEADBYTES,
  		fP) != WLZEFF_PIC_HEADBYTES)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  else
  {
    origin.vtX = 0;
    origin.vtY = 0;
    origin.vtZ = 0;
    WLZEFF_PIC_HEAD_WORD_GET(size.vtX, WLZEFF_PIC_OFF_NX, header);
    WLZEFF_PIC_HEAD_WORD_GET(size.vtY, WLZEFF_PIC_OFF_NY, header);
    WLZEFF_PIC_HEAD_WORD_GET(size.vtZ, WLZEFF_PIC_OFF_NPIC, header);
    WLZEFF_PIC_HEAD_WORD_GET(byteData, WLZEFF_PIC_OFF_BYTEPIX, header);
    WLZEFF_PIC_HEAD_WORD_GET(magic, WLZEFF_PIC_OFF_MAGIC, header);
    if((size.vtX <= 0) || (size.vtY <= 0) || (size.vtZ <= 0) ||
       ((byteData != 0) && (byteData != 1)) ||
       (magic != WLZEFF_PIC_MAGIC))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(byteData)
    {
      gType = WLZ_GREY_UBYTE;
      tI0 = sizeof(unsigned char) * size.vtZ * size.vtY * size.vtX;
      if((AlcUnchar3Malloc((unsigned char ****)&data, size.vtZ, size.vtY,
			   size.vtX) != ALC_ER_NONE) ||
	 (data == NULL))
      {
	data = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    else				  /* Only 8 and 16 bit data are used */
    {
      gType = WLZ_GREY_SHORT;
      tI0 = sizeof(short) * size.vtZ * size.vtY * size.vtX;
      if((AlcShort3Malloc((short ****)&data, size.vtZ, size.vtY,
			  size.vtX) != ALC_ER_NONE) ||
	 (data == NULL))
      {
	data = NULL;
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(data)
    {
      if(fread((unsigned char *)**data, sizeof(unsigned char),
	       tI0, fP) == tI0)
      {
	if(!byteData)
	{
	  swab((const void *)**data, (void *)**data, tI0);
	}
	obj = WlzFromArray3D((void ***)data, size, origin, gType, gType,
			     0.0, 1.0, 0, 1, &errNum);
      }
    }
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
      Alc3Free(data);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
* Function:	WlzEffWriteObjPic					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given Woolz object to the given stream 	*
*		using the BioRad Confocal '.pic' format. 		*
* Global refs:	-							*
* Global refs:	-							*
* Parameters:	FILE *fP:		Output file stream.		*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
WlzErrorNum	WlzEffWriteObjPic(FILE *fP, WlzObject *obj)
{
  int		tI0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzIVertex3	origin,
  		size;
  unsigned char	***data = NULL;
  unsigned char	header[WLZEFF_PIC_HEADBYTES];

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
    tI0 = size.vtX * size.vtY * size.vtZ;
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_NX, size.vtX);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_NY, size.vtY);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_NPIC, size.vtZ);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_BLACKVAL, 0);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_WHITEVAL, 255);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_BYTEPIX, 1);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_MERGED, 0);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_RGB, 7);
    WLZEFF_PIC_HEAD_WORD_SET(header, WLZEFF_PIC_OFF_MAGIC,
			       WLZEFF_PIC_MAGIC);
    if((fwrite(header, sizeof(char), WLZEFF_PIC_HEADBYTES,
	       fP) != WLZEFF_PIC_HEADBYTES) ||
       (fwrite(**data, sizeof(unsigned char), tI0, fP) != tI0))
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
