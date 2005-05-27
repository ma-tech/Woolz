#pragma ident "MRC HGU $Id$"
#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzExtFFIcs.c
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
*		and from the International Cytometry Standard '.ics'/'.ids'
*		data format.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Builds the ics file names from the given file name.
*		These strings should be free'd using AlcFree() when
*		no longer required.
* \param	fileBody		Dest ptr for the file body.
* \param	icsFileName		Dest ptr for the '.ics' file name.
* \param	idsFileName		Dest ptr for the '.ids' file name.
* \param	gvnFileName	Given file name with .ics, .ids or no
* 					extension.
*/
WlzErrorNum	WlzEffIcsFileNames(char **fileBody,
				   char **icsFileName,
				   char **idsFileName,
				   const char *gvnFileName)
{
  int		tI0;
  WlzErrorNum	errFlag = WLZ_ERR_MEM_ALLOC;

  tI0 = ((int )strlen(gvnFileName) + 5) * sizeof(char);
  if(((*fileBody = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*icsFileName = (char *)AlcMalloc(tI0)) != NULL) &&
     ((*idsFileName = (char *)AlcMalloc(tI0)) != NULL))
  {
    (void )strcpy(*fileBody, gvnFileName);
    if((tI0 = (int )strlen(*fileBody) - 4) >= 0)
    {
      if((strcmp(*fileBody +  tI0, ".ics") == 0) ||
	 (strcmp(*fileBody +  tI0, ".ids") == 0))
      {
	*(*fileBody +  tI0) = '\0';
      }
    }
    (void )sprintf(*icsFileName, "%s.ics", *fileBody);
    (void )sprintf(*idsFileName, "%s.ids", *fileBody);
    errFlag = WLZ_ERR_NONE;
  }
  return(errFlag);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the ICS
*		format. The given file name is used to generate the '.ics' and
*		'.ids' filenames.
* \param	gvnFileName		Given file name.
* \param	dstErr			Destination error number ptr,
*					may be NULL.
*/
WlzObject	*WlzEffReadObjIcs(const char *gvnFileName, WlzErrorNum *dstErr)
{
  int		bits,
		dataCount,
		datumSize,
  		idx;
  void		***data = NULL;
  FILE		*fP = NULL;
  char		*fileName = NULL,
  		*icsFileName = NULL,
		*idsFileName = NULL,
		*tokenStr;
  unsigned int	tokenId;
  WlzIVertex3	origin,
  		size;
  WlzEffIcsHeader header;
  WlzObject	*obj = NULL;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		recBuf[WLZEFF_ICS_REC_LEN_MAX];

  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzEffIcsFileNames(&fileName, &icsFileName, &idsFileName,
				gvnFileName);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((icsFileName == NULL) || (*icsFileName == '\0') ||
       ((fP = fopen(icsFileName, "r")) == NULL))
    {
      errNum = WLZ_ERR_READ_EOF;
    }

	  #ifdef _WIN32
  if (fP != NULL){
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
  #endif
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )memset((void *)&header, 0, sizeof(WlzEffIcsHeader));
    header.coords = WLZEFF_ICS_TKN_VIDEO;
    while((errNum == WLZ_ERR_NONE) &&
	  (fgets(recBuf, WLZEFF_ICS_REC_LEN_MAX - 1, fP)))
    {
      recBuf[WLZEFF_ICS_REC_LEN_MAX - 1] = '\0';
      if((tokenStr = strtok(recBuf, " \t\n")) != NULL)
      {
	tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
	if((WlzStringMatchValue((int *)&tokenId, tokenStr,
		      "ics_version", WLZEFF_ICS_TKN_VERSION,
		      "filename", WLZEFF_ICS_TKN_FILENAME,
		      "layout", WLZEFF_ICS_TKN_LAYOUT,
		      "representation", WLZEFF_ICS_TKN_REPRESENTATION,
		      NULL) == 0) ||
	   (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE) ||
	   ((tokenStr = strtok(NULL, " \t\n")) == NULL))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  switch(tokenId)
	  {
	    case WLZEFF_ICS_TKN_VERSION:
	      if((sscanf(tokenStr, "%d.%d", &(header.versionMajor),
		  &(header.versionMinor)) != 2) ||
		 (header.versionMajor != WLZEFF_ICS_VERSION_MAJOR) ||
		 (header.versionMinor != WLZEFF_ICS_VERSION_MINOR))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      break;
	    case WLZEFF_ICS_TKN_FILENAME:
	      break;	    /* Ignore the filename record of the '.ics' file */
	    case WLZEFF_ICS_TKN_LAYOUT:
	      tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
	      if((WlzStringMatchValue((int *)&tokenId, tokenStr,
		      "parameters", WLZEFF_ICS_TKN_PARAMETERS,
		      "order", WLZEFF_ICS_TKN_ORDER,
		      "sizes", WLZEFF_ICS_TKN_SIZES,
		      "coordinates", WLZEFF_ICS_TKN_COORDS,
		      "significant_bits", WLZEFF_ICS_TKN_SIGBITS,
		      NULL) == 0) ||
		 ((tokenStr = strtok(NULL, " \t\n")) == NULL))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		switch(tokenId)
		{
		  case WLZEFF_ICS_TKN_PARAMETERS:
		    if((sscanf(tokenStr, "%d", &(header.parameters)) != 1) ||
		       (header.parameters < 1) ||
		       (header.parameters > WLZEFF_ICS_PARAMETERS_MAX))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    break;
		  case WLZEFF_ICS_TKN_ORDER:
		    if(header.parameters < 1)
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      idx = 0;
		      while((idx < header.parameters) &&
			    (errNum == WLZ_ERR_NONE))
		      {
			tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
			if((tokenStr == NULL) ||
			   (WlzStringMatchValue((int *)&tokenId, tokenStr,
			      "bits", WLZEFF_ICS_TKN_BITS,
			      "x", WLZEFF_ICS_TKN_X,
			      "y", WLZEFF_ICS_TKN_Y,
			      "z", WLZEFF_ICS_TKN_Z,
			      NULL) == 0) ||
			   (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE))
			{
			  errNum = WLZ_ERR_READ_INCOMPLETE;
			}
			else
			{
			  header.order[idx] = (WlzEffIcsToken )tokenId;
			  tokenStr = strtok(NULL, " \t\n");
			  ++idx;
			}
		      }
		    }
		    break;
		  case WLZEFF_ICS_TKN_SIZES:
		    if(header.parameters < 1)
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      idx = 0;
		      while((idx < header.parameters) &&
			    (errNum == WLZ_ERR_NONE))
		      {
			if((tokenStr == NULL) ||
			   (sscanf(tokenStr, "%d",
				   &(header.sizes[idx])) != 1))
			{
			  errNum = WLZ_ERR_READ_INCOMPLETE;
			}
			else
			{
			  tokenStr = strtok(NULL, " \t\n");
			  ++idx;
			}
		      }
		    }
		    break;
		  case WLZEFF_ICS_TKN_COORDS:
		    tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
		    if((WlzStringMatchValue((int *)&tokenId, tokenStr,
			      "video", WLZEFF_ICS_TKN_VIDEO,
			      NULL) == 0) ||
		       (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      header.coords = (WlzEffIcsToken )tokenId;
		    }
		    break;
		  case WLZEFF_ICS_TKN_SIGBITS:
		    if((sscanf(tokenStr, "%d", &(header.sigBits)) != 1) ||
		       (header.sigBits < 1) ||
		       (header.sigBits > 32))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    break;
		  default:
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		    break;
		}
	      }
	      break;
	    case WLZEFF_ICS_TKN_REPRESENTATION:
	      tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
	      if((WlzStringMatchValue((int *)&tokenId, tokenStr,
		      "format", WLZEFF_ICS_TKN_FORMAT,
		      "sign", WLZEFF_ICS_TKN_SIGN,
		      "SCIL_TYPE", WLZEFF_ICS_TKN_SCIL,
		      NULL) == 0) ||
		 ((tokenStr = strtok(NULL, " \t\n")) == NULL))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
		switch(tokenId)
		{
		  case WLZEFF_ICS_TKN_FORMAT:
		    tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
		    if((WlzStringMatchValue((int *)&tokenId, tokenStr,
			      "integer", WLZEFF_ICS_TKN_INT,
			      "float", WLZEFF_ICS_TKN_FLOAT,
			      NULL) == 0) ||
		       (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      header.format = (WlzEffIcsToken )tokenId;
		    }
		    break;
		  case WLZEFF_ICS_TKN_SIGN:
		    tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
		    if((WlzStringMatchValue((int *)&tokenId, tokenStr,
			      "signed", WLZEFF_ICS_TKN_SIGNED,
			      "unsigned", WLZEFF_ICS_TKN_UNSIGNED,
			      NULL) == 0) ||
		       (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      header.sign = (WlzEffIcsToken )tokenId;
		    }
		    break;
		  case WLZEFF_ICS_TKN_SCIL:
		    tokenId = (unsigned int )WLZEFF_ICS_TKN_NONE;
		    if((WlzStringMatchValue((int *)&tokenId, tokenStr,
			      "g3d", WLZEFF_ICS_TKN_G3D,
			      "l3d", WLZEFF_ICS_TKN_L3D,
			      NULL) == 0) ||
		       (tokenId == (unsigned int )WLZEFF_ICS_TKN_NONE))
		    {
		      errNum = WLZ_ERR_READ_INCOMPLETE;
		    }
		    else
		    {
		      header.scil = (WlzEffIcsToken )tokenId;
		    }
		    break;
		  default:
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
      }
    }
    fclose(fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((header.versionMajor != WLZEFF_ICS_VERSION_MAJOR) ||
       (header.versionMinor != WLZEFF_ICS_VERSION_MINOR) ||
       (header.parameters != 4) || 			    /* 3D files only */
       (header.order[0] == WLZEFF_ICS_TKN_NONE) ||
       (header.sizes[0] == WLZEFF_ICS_TKN_NONE) ||
       (header.coords != WLZEFF_ICS_TKN_VIDEO) ||
       (header.format != WLZEFF_ICS_TKN_INT) ||
       ((header.scil != WLZEFF_ICS_TKN_G3D) &&
	(header.scil != WLZEFF_ICS_TKN_G3D)))
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      origin.vtX = 0;
      origin.vtY = 0;
      origin.vtZ = 0;
      for(idx = 0; idx < header.parameters; ++idx)
      {
	switch(header.order[idx])
	{
	  case WLZEFF_ICS_TKN_BITS:
	    bits = header.sizes[idx];
	    break;
	  case WLZEFF_ICS_TKN_X:
	    size.vtX = header.sizes[idx];
	    break;
	  case WLZEFF_ICS_TKN_Y:
	    size.vtY = header.sizes[idx];
	    break;
	  case WLZEFF_ICS_TKN_Z:
	    size.vtZ = header.sizes[idx];
	    break;
	}
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) &&
     (size.vtX > 0) && (size.vtY > 0) && (size.vtZ > 0) &&
     (header.format == WLZEFF_ICS_TKN_INT) &&
     idsFileName && ((fP = fopen(idsFileName, "r")) != NULL))
  {
      #ifdef _WIN32
  if (fP != NULL){
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
  #endif

  dataCount = size.vtX * size.vtY * size.vtZ;
    switch(bits)
    {
      case 8:
	if((AlcUnchar3Malloc((unsigned char ****)&data, size.vtZ, size.vtY,
			     size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(unsigned char);
	  gType = WLZ_GREY_UBYTE;
	}
	break;
      case 16:
	if((AlcShort3Malloc((short ****)&data, size.vtZ, size.vtY,
			   size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(short);
	  gType = WLZ_GREY_SHORT;
	}
	break;
      case 32:
	if((AlcInt3Malloc((int ****)&data, size.vtZ, size.vtY,
			   size.vtX) != ALC_ER_NONE) ||
	   (data == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  datumSize = sizeof(int);
	  gType = WLZ_GREY_INT;
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
    obj->domain.p->voxel_size[0] = 1.0;
    obj->domain.p->voxel_size[1] = 1.0;
    obj->domain.p->voxel_size[2] = 1.0;
  }
  if(fileName)
  {
    AlcFree(fileName);
  }
  if(icsFileName)
  {
    AlcFree(icsFileName);
  }
  if(idsFileName)
  {
    AlcFree(idsFileName);
  }
  if((errNum != WLZ_ERR_NONE) && data)
  {
    Alc3Free(data);
  }
  if(fP)
  {
    fclose(fP);
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
* \brief	Writes the given Woolz object to the given stream
*		using the ICS format. The given file name is used to generate
*		the '.ics' and '.ids' filenames.
* \param	gvnFileName		Given file name with .ics, .ids or
*					no extension.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjIcs(const char *gvnFileName, WlzObject *obj)
{
  int  		tI0;
  FILE		*fP = NULL;
  unsigned char	 ***data = NULL;
  WlzIVertex3	origin,
  		size;
  char		*tCP0,
  		*fileName = NULL,
  		*icsFileName = NULL,
		*idsFileName = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((gvnFileName == NULL) || (*gvnFileName == '\0'))
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
    errNum = WlzEffIcsFileNames(&fileName, &icsFileName, &idsFileName,
				gvnFileName);
  }
  if(errNum == WLZ_ERR_NONE)
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
    if((errNum == WLZ_ERR_NONE) && fileName && icsFileName)
    {
      if((fP = fopen(icsFileName, "w")) == NULL)
      {
        errNum = WLZ_ERR_WRITE_EOF;
      }
      else
      {
		    #ifdef _WIN32
  if (fP != NULL){
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
  #endif

	if((tCP0 = strrchr(fileName, '/')) == NULL)
	{
	  tCP0 = fileName;
	}
	else
	{
	  ++tCP0;
	}
	if(fprintf(fP, "\t\nics_version\t%d.%d\n"
		       "filename\t%s\n"
		       "layout\tparameters\t4\n"
		       "layout\torder\tbits\tx\ty\tz\n"
		       "layout\tsizes\t8\t%d\t%d\t%d\n"
		       "layout\tcoordinates\tvideo\n"
		       "layout\tsignificant_bits\t8\n"
		       "representation\tformat\tinteger\n"
		       "representation\tsign\tunsigned\n"
		       "representation\tSCIL_TYPE\tg3d\n",
		       WLZEFF_ICS_VERSION_MAJOR,
		       WLZEFF_ICS_VERSION_MINOR,
		       tCP0,
		       size.vtX, size.vtY, size.vtZ) <= 0)
	{
	  errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	fclose(fP);
      }
    }
  }
  if((errNum == WLZ_ERR_NONE) && idsFileName)
  {
    if((fP = fopen(idsFileName, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    else
    {
		  #ifdef _WIN32
  if (fP != NULL){
	if(_setmode(_fileno(fP), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
  #endif

      tI0 = size.vtX * size.vtY * size.vtZ;
      if(fwrite(**data, sizeof(unsigned char), tI0, fP) != tI0)
      {
	errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      fclose(fP);
    }
  }
  if(fileName)
  {
    AlcFree(fileName);
  }
  if(icsFileName)
  {
    AlcFree(icsFileName);
  }
  if(idsFileName)
  {
    AlcFree(idsFileName);
  }
  if(data)
  {
    (void )AlcUnchar3Free(data);
  }
  return(errNum);
}

