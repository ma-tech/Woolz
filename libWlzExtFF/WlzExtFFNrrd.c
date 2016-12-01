#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFNrrd_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFNrrd.c
* \author       Bill Hill
* \date         July 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Functions for reading and writting Woolz objects to and from
* 		the Utah nearly raw raster data (NRRD) '.nrrd' single file
* 		data format.
* \ingroup	WlzExtFF
*/

#include <string.h>
#include <float.h>
#include <time.h>
#if HAVE_ZLIB != 0
#include <zlib.h>
#endif /* HAVE_ZLIB */
#if HAVE_BZLIB != 0
#include <bzlib.h>
#endif /* HAVE_BZLIB */
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzErrorNum		WlzEffHeadReadNrrd(
				  WlzEffNrrdHeader *header,
				  FILE *fP);
static WlzObject		*WlzEffReadImgNrrd(
				  FILE *fP,
				  WlzEffNrrdHeader *header,
				  WlzErrorNum *dstErr);
static WlzErrorNum 		WlzEffWriteImgNrrd(
				  FILE *fP,
				  WlzObject *obj);

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz object from the given stream using the
* 		NRRD file format. Not all NRRD files and not all encodings
* 		have been implemented. This function will currently only read
* 		pixel/voxel data, where this data is either raw or encoded
* 		using bzip2 or gzip compression. Use of these compression
* 		methods depends on their availability at compile time.
* 		If the input NRRD pixel/voxel type is either signed char or
* 		unsigned short then it is promoted to either signed short
* 		or signed int in the output Woolz object.
* 		Unused NRRD header fields are encoded as a text property of
* 		the output object.
* \param	fP			Input file stream.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObjNrrd(FILE *fP, WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzEffNrrdHeader header;

  if(fP == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    errNum = WlzEffHeadReadNrrd(&header, fP);
    if(errNum == WLZ_ERR_NONE)
    {
      switch(header.dimension)
      {
        case 2: /* FALLTHROUGH */
	case 3: /* FALLTHROUGH */
	case 4:
	  obj = WlzEffReadImgNrrd(fP, &header, &errNum);
	  break;
        default:
          errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
      }
    }
    AlcFree(header.space);
    AlcFree(header.spaceDirections);
    AlcFree(header.spaceUnits);
    AlcFree(header.thickness);
    AlcFree(header.centers);
    AlcFree(header.centering);
    AlcFree(header.labels);
    AlcFree(header.kinds);
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
* 		NRRD file format. This function will currently only write
* 		2 or 3D domain objects. If the object has no values then
* 		unsigned byte values are written with background = 0 and
* 		forground = 255. If gzip compression is available at compile
* 		time the output file uses gzip compression, otherwise raw
* 		data are output.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
WlzErrorNum	WlzEffWriteObjNrrd(FILE *fP, WlzObject *obj)
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
      case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
      case WLZ_3D_DOMAINOBJ:
        errNum = WlzEffWriteImgNrrd(fP, obj);
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
* \brief	Writes the given Woolz 2 or 3D domain object (with or
* 		without values) to the given stream using the NRRD
* 		file format.
* \param	fP			Output file stream.
* \param	obj			Given woolz object.
*/
static WlzErrorNum WlzEffWriteImgNrrd(FILE *fP, WlzObject *obj)
{
  int		ln,
  		greySz = 1;
  size_t	nBytes = 0,
  		nData = 0;
  time_t	tm;
  char		*dateS,
		*dimS = NULL,
  		*typeS = NULL;
  char		orgS[64],
                sizeS[64],
		spaceS[64];
  WlzVertex	sz,
  		org;
  WlzIBox3	bBox;
  WlzGreyP	data;
  WlzGreyType gType = WLZ_GREY_UBYTE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
#if HAVE_ZLIB != 0
  const char    *encodingS = "gzip";
#else
  const char    *encodingS = "raw";
#endif /* HAVE_ZLIB */

  data.v = NULL;
  tm = time(NULL);
  dateS = ctime(&tm);
  ln = strlen(dateS);
  while((--ln > 0) && ((dateS[ln] == '\n') || (dateS[ln] == '\r')))
  {
    dateS[ln] = '\0';		 /* Remove trailing \n and/or \r from dateS. */
  }
  bBox = WlzBoundingBox3I(obj, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    if(obj->values.core)
    {
      gType = WlzGreyTypeFromObj(obj, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      switch(gType)
      {
        case WLZ_GREY_INT:
	  greySz = 4;
	  typeS = "int32";
	  break;
        case WLZ_GREY_SHORT:
	  greySz = 2;
	  typeS = "int16";
	  break;
        case WLZ_GREY_UBYTE:
	  greySz = 1;
	  typeS = "uint8";
	  break;
        case WLZ_GREY_FLOAT:
	  greySz = 4;
	  typeS = "float";
	  break;
        case WLZ_GREY_DOUBLE:
	  greySz = 8;
	  typeS = "double";
	  break;
        case WLZ_GREY_RGBA:
	  greySz = 4;
	  typeS = "uint32";
	  break;
        default:
	  errNum = WLZ_ERR_GREY_TYPE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
        dimS = "2";
	sz.i2.vtX = bBox.xMax - bBox.xMin + 1;
	sz.i2.vtY = bBox.yMax - bBox.yMin + 1;
	org.i2.vtX = bBox.xMin;
	org.i2.vtY = bBox.yMin;
	(void )sprintf(sizeS, "%d %d", sz.i2.vtX, sz.i2.vtY);
	(void )sprintf(orgS, "(%d,%d)", org.i2.vtX, org.i2.vtY);
	(void )sprintf(spaceS, "nan nan");
	break;
      case WLZ_3D_DOMAINOBJ:
        dimS = "3";
	sz.i3.vtX = bBox.xMax - bBox.xMin + 1;
	sz.i3.vtY = bBox.yMax - bBox.yMin + 1;
	sz.i3.vtZ = bBox.zMax - bBox.zMin + 1;
	org.i3.vtX = bBox.xMin;
	org.i3.vtY = bBox.yMin;
	org.i3.vtZ = bBox.zMin;
	(void )sprintf(sizeS, "%d %d %d", sz.i3.vtX, sz.i3.vtY, sz.i3.vtZ);
	(void )sprintf(orgS, "(%d,%d,%d)", org.i3.vtX, org.i3.vtY, org.i3.vtZ);
	(void )sprintf(spaceS, "%f %f %f",
	               obj->domain.p->voxel_size[0],
	               obj->domain.p->voxel_size[1],
	               obj->domain.p->voxel_size[2]);
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fprintf(fP,
	       "NRRD0004\n"
	       "# This NRRD file was generated from Woolz\n"
	       "# on %s.\n"
	       "# Complete NRRD file format specification at:\n"
	       "# http://teem.sourceforge.net/nrrd/format.html\n"
	       "type: %s\n"
	       "dimension: %s\n"
	       "sizes: %s\n"
	       "kinds: domain domain domain\n"
	       "endian: little\n"
	       "encoding: %s\n"
	       "space origin: %s\n"
	       "spacings: %s\n\n",
	       dateS, typeS, dimS, sizeS, encodingS, orgS, spaceS) < 256)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzObject *nObj = NULL;

    if(obj->values.core == NULL)
    {
      WlzObjectType gTT;
      WlzPixelV bV,
      		fV;

      bV.type = fV.type = WLZ_GREY_UBYTE;
      bV.v.ubv = 0;
      fV.v.ubv = 255;
      gTT = WlzGreyValueTableType(0, WLZ_GREY_TAB_RAGR, WLZ_GREY_UBYTE, NULL);
      nObj = WlzNewObjectValues(obj, gTT, bV, 1, fV, &errNum);
    }
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	{
	  void 	**ary = NULL;

	  nData = sz.i2.vtX * sz.i2.vtY;
          nBytes = greySz * nData;
	  errNum = WlzToArray2D(&ary, (nObj)? nObj: obj, sz.i2, org.i2,
				0, gType);
	  if(ary)
	  {
	    data.v = *ary;
	    AlcFree(ary);
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
	{
	  void  ***ary = NULL;

	  nData = sz.i3.vtX * sz.i3.vtY * sz.i3.vtZ;
          nBytes = greySz * nData;
	  errNum = WlzToArray3D(&ary, (nObj)? nObj: obj, sz.i3, org.i3,
				0, gType);
	  if(ary)
	  {
	    data.v = **ary;
	    AlcFree(*ary);
	    AlcFree(ary);
	  }
	}
	break;
      default:
        errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    (void )WlzFreeObj(nObj);
  }
  if(errNum == WLZ_ERR_NONE)
#if HAVE_ZLIB != 0
  {
    size_t 	maxSz;
    void	*cBuf = NULL;


    maxSz = 1024 + ((nBytes * 110) / 100); /* Add 1k + 10% and hope that this
    				              copes with bad compression. */
    if((cBuf = AlcMalloc(maxSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      int	zer;
      z_stream	zs;

      zs.zalloc = NULL;
      zs.zfree = NULL;
      zs.next_in = data.v;
      zs.avail_in = nBytes;
      zs.next_out = cBuf;
      zs.avail_out = maxSz;
      if((zer = deflateInit(&zs, Z_BEST_COMPRESSION)) == Z_OK)
      {
        if((zer = deflate(&zs, Z_FINISH)) == Z_STREAM_END)
	{
	  nBytes = zs.total_out;
	  if(fwrite(cBuf, sizeof(char), nBytes, fP) != nBytes)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	}
	else
	{
          errNum = WLZ_ERR_WRITE_INCOMPLETE;
	}
	deflateEnd(&zs);
      }
      else
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      AlcFree(cBuf);
    }
  }
#else
  {
    if(fwrite(data.v, greySz, nData, fP) != nBytes)
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
#endif /* HAVE_ZLIB */
  AlcFree(data.v);
  return(errNum);
}

/*!
* \return	Woolz error number.
* \ingroup	WlzExtFF
* \brief	Reads the NRRD file header.
* \param	header			Given header data structure be filled.
* \param	fP			Output file stream.
*/
static WlzErrorNum WlzEffHeadReadNrrd(WlzEffNrrdHeader *header, FILE *fP)
{
  int		i,
  		fldE;
  char		*fldS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		buf[256];

  /* Set default header values. */
  (void )memset((void *)header, 0, sizeof(WlzEffNrrdHeader));
  for(i = 0; i < WLZEFF_NRRD_MAX_DIMENSION; ++i)
  {
    header->spacings[i] = 1.0;
  }
  /* Read and check the identifier and version number. */
  if(fgets(buf, 256, fP) == NULL)
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    buf[255] = '\0';
    if(sscanf(buf, "NRRD000%d", &(header->version)) != 1)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    /* Ignore the version number because we can probably read it anyway. */
  }
  /* Read and parse all remaining records. */
  if(errNum == WLZ_ERR_NONE)
  {
    int		finished = 0;
    do
    {
      if(fgets(buf, 256, fP) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
      }
      buf[255] = '\0';
      fldS = strtok(buf, ":\n\r");
      if((fldS == NULL) || (*fldS == '\0'))
      {
        finished = 1;
      }
      else if(*fldS != '#')
      {
	int	dscE;
	char	*dscS;

        (void )WlzStringWhiteSpSkip(fldS);
        (void )WlzStringToLower(fldS);
	if(WlzStringMatchValue(&fldE, fldS,
			"axismaxs",          WLZEFF_NRRD_FIELD_AXIS_MAX,        
			"axismins",          WLZEFF_NRRD_FIELD_AXIS_MIN,        
			"blocksize",         WLZEFF_NRRD_FIELD_BLOCK_SIZE,      
			"byteskip",          WLZEFF_NRRD_FIELD_BYTE_SKIP,       
			"centerings",        WLZEFF_NRRD_FIELD_CENTERING,       
			"center",            WLZEFF_NRRD_FIELD_CENTERS,         
			"datafile",          WLZEFF_NRRD_FIELD_DATA_FILE,       
			"dimension",         WLZEFF_NRRD_FIELD_DIMENSION,       
			"encoding",          WLZEFF_NRRD_FIELD_ENCODING,        
			"endian",            WLZEFF_NRRD_FIELD_ENDIAN,          
			"kinds",             WLZEFF_NRRD_FIELD_KINDS,           
			"labels",            WLZEFF_NRRD_FIELD_LABELS,          
			"lineskip",          WLZEFF_NRRD_FIELD_LINE_SKIP,       
			"measurementframe",  WLZEFF_NRRD_FIELD_MEASUREMENT_FRM, 
			"number",            WLZEFF_NRRD_FIELD_NUMBER,          
			"sizes",             WLZEFF_NRRD_FIELD_SIZES,           
			"spacedimension",    WLZEFF_NRRD_FIELD_SPACE_DIM,       
			"spacedirections",   WLZEFF_NRRD_FIELD_SPACE_DIR,       
			"spaceorigin",       WLZEFF_NRRD_FIELD_SPACE_ORG,       
			"spaceunits",        WLZEFF_NRRD_FIELD_SPACE_UNITS,     
			"space",             WLZEFF_NRRD_FIELD_SPACE,           
			"spacings",          WLZEFF_NRRD_FIELD_SPACINGS,        
			"thickness",         WLZEFF_NRRD_FIELD_THICKNESS,       
			"type",              WLZEFF_NRRD_FIELD_TYPE,            
		  NULL) == 0)
	{
	  fldE = WLZEFF_NRRD_FIELD_UNSUP;
	}
	switch(fldE)
	{
          case WLZEFF_NRRD_FIELD_AXIS_MAX:
	    /* Ignore. */
	    break;
          case WLZEFF_NRRD_FIELD_AXIS_MIN:
	    /* Ignore. */
	    break;
          case WLZEFF_NRRD_FIELD_BLOCK_SIZE:
	    if(((dscS = strtok(NULL, " \n\r")) == NULL) ||
	       (sscanf(dscS, "%zd", &(header->blockSz)) != 1) ||
	       (header->blockSz <= 0))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_BYTE_SKIP:
	    if(((dscS = strtok(NULL, " \n\r")) == NULL) ||
	       (sscanf(dscS, "%zd", &(header->byteSkip)) != 1))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_CENTERING:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->centering = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_CENTERS:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->centers = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_DATA_FILE:
	    /* Separate data files are not supported. */
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	    break;
          case WLZEFF_NRRD_FIELD_DIMENSION:
	    dscS = strtok(NULL, " \n\r");
	    if((dscS == NULL) || (*dscS == '\0') ||
	       (header->dimension != 0) ||
	       (sscanf(dscS, "%d", &(header->dimension)) != 1) ||
	       (header->dimension < 1) ||
	       (header->dimension > WLZEFF_NRRD_MAX_DIMENSION))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_ENCODING:
	    dscS = strtok(NULL, " \n\r");
	    if((dscS == NULL) || (*dscS == '\0') ||
	       (WlzStringMatchValue(&dscE, dscS,
			"bz2",   WLZEFF_NRRD_ENCODE_BZ2,
			"bzip2", WLZEFF_NRRD_ENCODE_BZ2,
			"gz",    WLZEFF_NRRD_ENCODE_GZ,
			"gzip",  WLZEFF_NRRD_ENCODE_GZ,
			"hex",   WLZEFF_NRRD_ENCODE_HEX,
			"raw",   WLZEFF_NRRD_ENCODE_RAW,
			"txt",   WLZEFF_NRRD_ENCODE_TXT,
			"text",  WLZEFF_NRRD_ENCODE_TXT,
			"ascii", WLZEFF_NRRD_ENCODE_TXT,
			NULL) == 0))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      header->encoding = (WlzEffNrrdDataType )dscE;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_ENDIAN:
	    dscS = strtok(NULL, " \n\r");
	    if((dscS == NULL) || (*dscS == '\0') ||
	       (WlzStringMatchValue(&dscE, dscS,
			"big",    WLZEFF_NRRD_ENDIAN_BIG,
			"little", WLZEFF_NRRD_ENDIAN_LITTLE,
			NULL) == 0))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      header->endian = (WlzEffNrrdEndian )dscE;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_KINDS:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->kinds = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_LABELS:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->labels = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_LINE_SKIP:
	    if(((dscS = strtok(NULL, " \n\r")) == NULL) ||
	       (sscanf(dscS, "%zd", &(header->lineSkip)) != 1))
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    break;
          case WLZEFF_NRRD_FIELD_MEASUREMENT_FRM:
	    /* Ignore. */
	    break;
          case WLZEFF_NRRD_FIELD_NUMBER:
	    /* Ignore. */
	    break;
          case WLZEFF_NRRD_FIELD_SIZES:
	    if(header->dimension < 1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      for(i = 0; i < header->dimension; ++i)
	      {
	        if(((dscS = strtok(NULL, " \n\r")) == NULL) ||
		   (sscanf(dscS, "%d", &(header->sizes[i])) != 1))

	        {
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		  break;
		}
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_SPACE:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->space = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_SPACE_DIM:
	    /* Ignore. */
	    break;
          case WLZEFF_NRRD_FIELD_SPACE_DIR:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->spaceDirections = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_SPACE_ORG:
	    if(header->dimension < 1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      int	n;

	      n = (header->dimension > 3)? 3: header->dimension;
	      for(i = 0; i < n; ++i)
	      {
	        if(((dscS = strtok(NULL, " ,()\n\r")) == NULL) ||
		   (sscanf(dscS, "%lg", &(header->origin[i])) != 1))
	        
	        {
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		  break;
		}
	      }
            }
	    break;
          case WLZEFF_NRRD_FIELD_SPACE_UNITS:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->spaceUnits = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_SPACINGS:
	    if(header->dimension < 1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      for(i = 0; i < header->dimension; ++i)
	      {
	        if(((dscS = strtok(NULL, " \n\r")) == NULL) ||
		   (sscanf(dscS, "%lg", &(header->spacings[i])) != 1))

	        {
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		  break;
		}
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_THICKNESS:
	    if((dscS = strtok(NULL, "\n\r")) != NULL)
	    {
	      dscS = WlzStringWhiteSpSkipLeading(dscS);
	      if((header->thickness = AlcStrDup(dscS)) == NULL)
	      {
	        errNum = WLZ_ERR_MEM_ALLOC;
	      }
	    }
	    break;
          case WLZEFF_NRRD_FIELD_TYPE:
	    dscS = strtok(NULL, "\n\r");
	    dscS = WlzStringWhiteSpSkipLeading(dscS);
	    if((dscS == NULL) || (*dscS == '\0'))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      if(WlzStringMatchValue(&dscE, dscS,
			"int8",                    WLZEFF_NRRD_TYPE_INT8,
			"int8_t",                  WLZEFF_NRRD_TYPE_INT8,
			"signed char",             WLZEFF_NRRD_TYPE_INT8,
			"uchar",                   WLZEFF_NRRD_TYPE_UINT8,
			"uint8",                   WLZEFF_NRRD_TYPE_UINT8,
			"uint8_t",                 WLZEFF_NRRD_TYPE_UINT8,
			"unsigned char",           WLZEFF_NRRD_TYPE_UINT8,
			"int16",                   WLZEFF_NRRD_TYPE_INT16,
			"int16_t",                 WLZEFF_NRRD_TYPE_INT16,
			"short",                   WLZEFF_NRRD_TYPE_INT16,
			"short int",               WLZEFF_NRRD_TYPE_INT16,
			"signed short",            WLZEFF_NRRD_TYPE_INT16,
			"signed short int",        WLZEFF_NRRD_TYPE_INT16,
			"uint16",                  WLZEFF_NRRD_TYPE_UINT16,
			"uint16_t",                WLZEFF_NRRD_TYPE_UINT16,
			"ushort",                  WLZEFF_NRRD_TYPE_UINT16,
			"unsigned short",          WLZEFF_NRRD_TYPE_UINT16,
			"unsigned short int",      WLZEFF_NRRD_TYPE_UINT16,
			"int",                     WLZEFF_NRRD_TYPE_INT32,
			"int32",                   WLZEFF_NRRD_TYPE_INT32,
			"int32_t",                 WLZEFF_NRRD_TYPE_INT32,
			"signed int",              WLZEFF_NRRD_TYPE_INT32,
			"uint",                    WLZEFF_NRRD_TYPE_UINT32,
			"uint32",                  WLZEFF_NRRD_TYPE_UINT32,
			"uint32_t",                WLZEFF_NRRD_TYPE_UINT32,
			"unsigned int",            WLZEFF_NRRD_TYPE_UINT32,
			"int64",                   WLZEFF_NRRD_TYPE_INT64,
			"int64_t",                 WLZEFF_NRRD_TYPE_INT64,
			"longlong",                WLZEFF_NRRD_TYPE_INT64,
			"long long",               WLZEFF_NRRD_TYPE_INT64,
			"long long int",           WLZEFF_NRRD_TYPE_INT64,
			"signed long long",        WLZEFF_NRRD_TYPE_INT64,
			"signed long long int",    WLZEFF_NRRD_TYPE_INT64,
			"uint64",                  WLZEFF_NRRD_TYPE_UINT64,
			"uint64_t",                WLZEFF_NRRD_TYPE_UINT64,
			"ulonglong",               WLZEFF_NRRD_TYPE_UINT64,
			"unsigned long long",      WLZEFF_NRRD_TYPE_UINT64,
			"unsigned long long int",  WLZEFF_NRRD_TYPE_UINT64,
			"float",                   WLZEFF_NRRD_TYPE_FLOAT,
			"double",                  WLZEFF_NRRD_TYPE_DOUBLE,
			"block",                   WLZEFF_NRRD_TYPE_BLOCK,
		     NULL) == 0)
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
	        header->dataType = (WlzEffNrrdDataType )dscE;
	      }
	    }
	    break;
	  default:
	    /* Ignore unknown fields. */
	    break;
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  finished = 1;
	}
      }
    } while(!finished);
  }
  return(errNum);
}

/*!
* \return	Object read from file.
* \ingroup	WlzExtFF
* \brief	Reads a Woolz domain object from the given stream using
*		the NRRD file format.
* \param	fP			Input file stream.
* \param	header			Header data structure.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
static WlzObject *WlzEffReadImgNrrd(FILE *fP, WlzEffNrrdHeader *header,
				    WlzErrorNum *dstErr)
{
  int		vLen = 1;
  size_t	nData = 0,
  		nrrdSz = 0,
  		wlzSz = 0;
  WlzGreyP	buf;
  WlzGreyType	gType = WLZ_GREY_ERROR;
  WlzIVertex3	objSz,
	        objOrg;
  WlzObjectType objType = WLZ_NULL;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* First check the header is valid while computing the number of data,
   * their size and the appropriate Woolz grey type. */
  buf.v = NULL;
  if((header->version < 1) ||
     (header->encoding == WLZEFF_NRRD_ENCODE_UNSUP) ||
     (header->dataType == WLZEFF_NRRD_TYPE_UNSUP))
  {
    errNum = WLZ_ERR_READ_INCOMPLETE;
  }
  else
  {
    switch(header->dimension)
    {
      case 2: /* FALLTHROUGH */
      case 3:
        break;
      case 4:
	if(strncmp("vector", header->kinds, 6) != 0)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
        break;
      default:
        errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    int		i;

    nData = 1;
    for(i = 0; i < header->dimension; ++i)
    {
      if(header->sizes[i] < 1)
      {
        errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      nData *= header->sizes[i];
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(header->dataType)
    {
      case WLZEFF_NRRD_TYPE_INT8:
	nrrdSz = 1;
        gType = WLZ_GREY_SHORT;
	wlzSz  = sizeof(short);
        break;
      case WLZEFF_NRRD_TYPE_UINT8:
	nrrdSz = 1;
        gType = WLZ_GREY_UBYTE;
	wlzSz  = sizeof(WlzUByte);
        break;
      case WLZEFF_NRRD_TYPE_INT16:
	nrrdSz = 2;
        gType = WLZ_GREY_SHORT;
	wlzSz  = sizeof(short);
        break;
      case WLZEFF_NRRD_TYPE_UINT16:
	nrrdSz = 2;
        gType = WLZ_GREY_INT;
	wlzSz  = sizeof(int);
        break;
      case WLZEFF_NRRD_TYPE_INT32:
	nrrdSz = 4;
        gType = WLZ_GREY_INT;
	wlzSz  = sizeof(int);
        break;
      case WLZEFF_NRRD_TYPE_FLOAT:
	nrrdSz = 4;
        gType = WLZ_GREY_FLOAT;
	wlzSz  = sizeof(float);
        break;
      case WLZEFF_NRRD_TYPE_DOUBLE:
	nrrdSz = 8;
        gType = WLZ_GREY_DOUBLE;
	wlzSz  = sizeof(double);
        break;
      default:
        /* Don't attempt to handle anything else. */
	errNum = WLZ_ERR_READ_INCOMPLETE;
	break;
    }
  }
  /* Read the NRRD data. */
  if(errNum == WLZ_ERR_NONE)
  {
    size_t	maxSz,
	  	nNrrdBytes;

    maxSz = nData * ((nrrdSz > wlzSz)? nrrdSz: wlzSz);
    maxSz = 1024 + ((maxSz * 110) / 100); /* Add 1k + 10% and hope that this
    				             copes with bad compression. */
    nNrrdBytes = nrrdSz * nData;
    if((buf.v = AlcMalloc(maxSz)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      switch(header->encoding)
      {
        case WLZEFF_NRRD_ENCODE_BZ2:
#if HAVE_BZLIB != 0
	  {
	    size_t	nRead;
	    void	*cBuf = NULL;

	    /* This BZ2 code seems OK, but it's hard to find NRRD datasets
	     * to test it with. */
	    if((nRead = fread(buf.v, 1, nNrrdBytes, fP)) < 1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else if((cBuf = AlcMalloc(nRead)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      int 	bze;
	      bz_stream bzs;

	      (void )memset(&bzs, 0, sizeof(bz_stream));
	      (void )memcpy(cBuf, buf.v, nRead);
	      bzs.next_in  = (char *)cBuf;
	      bzs.next_out = (char *)(buf.v);
	      bzs.avail_in = nRead;
	      bzs.avail_out = nNrrdBytes;
	      if((bze = BZ2_bzDecompressInit(&bzs, 0, 0)) == BZ_OK)
	      {
	        if((bze = BZ2_bzDecompress(&bzs)) != BZ_STREAM_END)
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
		BZ2_bzDecompressEnd(&bzs);
	      }
	      else
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      AlcFree(cBuf);
	    }
	  }
#else
          errNum = WLZ_ERR_READ_INCOMPLETE;
#endif /* HAVE_BZLIB */
	  break;
	case WLZEFF_NRRD_ENCODE_GZ:
#if HAVE_ZLIB != 0
	  {
	    size_t	nRead;
	    void	*cBuf = NULL;

	    if((nRead = fread(buf.v, 1, nNrrdBytes, fP)) < 1)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else if((cBuf = AlcMalloc(nRead)) == NULL)
	    {
	      errNum = WLZ_ERR_MEM_ALLOC;
	    }
	    else
	    {
	      int	zer;
	      z_stream 	zs;

	      /* Here we need to use zlib's inflate() rather than uncompress()
	       * as the NRRD compression uses gzip formats. */
	      (void )memset(&zs, 0, sizeof(z_stream));
	      (void )memcpy(cBuf, buf.v, nRead);
	      zs.total_in  = zs.avail_in  = nRead;
	      zs.total_out = zs.avail_out = nNrrdBytes;
	      zs.next_in   = (Bytef *)cBuf;
	      zs.next_out  = (Bytef *)(buf.v);
	      if((zer = inflateInit2(&zs, (15 + 32))) == Z_OK)
	      {
		if((zer = inflate(&zs, Z_FINISH)) == Z_STREAM_END)
		{
		  if(zs.total_out != nNrrdBytes)
		  {
		    errNum = WLZ_ERR_READ_INCOMPLETE;
		  }
		}
		else
		{
		  errNum = WLZ_ERR_READ_INCOMPLETE;
		}
		inflateEnd(&zs);
	      }
	      else
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      AlcFree(cBuf);
	    }
	  }
#else
          errNum = WLZ_ERR_READ_INCOMPLETE;
#endif /* HAVE_ZLIB */
	  break;
	case WLZEFF_NRRD_ENCODE_RAW:
	  if(fread(buf.v, nrrdSz, nData, fP) != nData)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  break;
	case WLZEFF_NRRD_ENCODE_HEX: /* FALLTHROUGH */
	case WLZEFF_NRRD_ENCODE_TXT: /* FALLTHROUGH */
	default:
	  /* All other encodings including text/asci and hex are
	   * unsupported. */
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Promote values if needed. */
    if(wlzSz != nrrdSz)
    {
      size_t	i;

      switch(gType)
      {
	case WLZ_GREY_SHORT:
	  for(i = nData - 1; i >= 0; --i)
	  {
	    buf.shp[i] = ((char *)(buf.ubp))[i];
	  }
	  break;
	case WLZ_GREY_INT:
	  for(i = nData - 1; i > 0; --i)
	  {
	    buf.inp[i] = ((unsigned short *)(buf.shp))[i];
	  }
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;         /* This should never happen. */
	  break;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Create a Woolz object with values allocated. */
    switch(header->dimension)
    {
      case 2:
        objType = WLZ_2D_DOMAINOBJ;
	WLZ_VTX_3_SET(objSz,
		      header->sizes[0],
		      header->sizes[1],
		      0);
	WLZ_VTX_3_SET(objOrg,
		      ALG_NINT(header->origin[0]),
		      ALG_NINT(header->origin[1]),
		      0);
	break;
      case 3:
        objType = WLZ_3D_DOMAINOBJ;
	WLZ_VTX_3_SET(objSz, header->sizes[0], header->sizes[1],
	              header->sizes[2]);
	WLZ_VTX_3_SET(objOrg,
		      ALG_NINT(header->origin[0]),
		      ALG_NINT(header->origin[1]),
		      ALG_NINT(header->origin[2]));
	break;
      case 4:
        objType = WLZ_3D_DOMAINOBJ;
	vLen = header->sizes[0];
	WLZ_VTX_3_SET(objSz, header->sizes[1], header->sizes[2],
	              header->sizes[3]);
	WLZ_VTX_3_SET(objOrg,
		      ALG_NINT(header->origin[1]),
		      ALG_NINT(header->origin[2]),
		      ALG_NINT(header->origin[3]));
	break;
      default:
        /* Handled on entry. */
	break;
    }
    obj = WlzFromArray1D(objType, objSz, objOrg, gType, vLen, buf, 0, &errNum);
  }
  AlcFree(buf.v);
  if(errNum == WLZ_ERR_NONE)
  {
    size_t	len;
    char	*infS = NULL;
    WlzProperty p;
    WlzPropertyList *pLst = NULL;

    p.core = NULL;
    len = 128 +
          ((header->space)? strlen(header->space): 0) +
	  ((header->spaceDirections)? strlen(header->spaceDirections): 0) +
	  ((header->spaceUnits)? strlen(header->spaceUnits): 0) +
	  ((header->thickness)? strlen(header->thickness): 0) +
	  ((header->centers)? strlen(header->centers): 0) +
	  ((header->centering)? strlen(header->centering): 0) +
	  ((header->labels)? strlen(header->labels): 0);
    if(objType == WLZ_3D_DOMAINOBJ)
    {
      obj->domain.p->voxel_size[0] = header->spacings[0];
      obj->domain.p->voxel_size[1] = header->spacings[1];
      obj->domain.p->voxel_size[2] = header->spacings[2];
    }
    if(((infS = (char *)AlcMalloc(len)) == NULL) ||
       ((pLst = WlzMakePropertyList(NULL)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      (void )sprintf(infS,
		     "space: %s\n"
		     "space directions: %s\n"
		     "space units: %s\n"
		     "thickness: %s\n"
		     "centers: %s\n"
		     "centering: %s\n"
		     "labels: %s\n"
		     "kinds: %s\n",
		     (header->space)? header->space: "",
		     (header->spaceDirections)? header->spaceDirections: "",
		     (header->spaceUnits)? header->spaceUnits: "",
		     (header->thickness)? header->thickness: "",
		     (header->centers)? header->centers: "",
		     (header->centering)? header->centering: "",
		     (header->labels)? header->labels: "",
		     (header->kinds)? header->kinds: "");
      p.text = WlzMakeTextProperty("Nrrd Info", infS, &errNum);
    }
    AlcFree(infS);
    if(errNum == WLZ_ERR_NONE)
    {
      if(AlcDLPListEntryAppend(pLst->list, NULL, (void *)(p.core),
                               WlzFreePropertyListEntry) != ALC_ER_NONE)
      {
        errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      obj->plist = WlzAssignPropertyList(pLst, NULL);
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(obj)
    {
      (void )WlzFreeObj(obj);
      obj = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}
