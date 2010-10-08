#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFF_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFF.c
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
* 		external data formats.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/*!
* \return	File format, WLZEFF_FORMAT_NONE if string not matched.
* \ingroup	WlzExtFF
* \brief	Returns a file format for the given file extension string.
* \param	extStr			Given file extension string.
*/
WlzEffFormat	WlzEffStringExtToFormat(const char *extStr)
{
  unsigned int	fileFmt;

  if(WlzStringMatchValue((int *)&fileFmt, extStr,
			 "am",    WLZEFF_FORMAT_AM,
			 "bmp",   WLZEFF_FORMAT_BMP,
			 "den",   WLZEFF_FORMAT_DEN,
			 "gif",   WLZEFF_FORMAT_GIF,
			 "hdr",   WLZEFF_FORMAT_ANL,
			 "ics",   WLZEFF_FORMAT_ICS,
			 "ids",   WLZEFF_FORMAT_ICS,
			 "img",   WLZEFF_FORMAT_ANL,
			 "ipl",   WLZEFF_FORMAT_IPL,
			 "jpeg",  WLZEFF_FORMAT_JPEG,
			 "jpg",   WLZEFF_FORMAT_JPEG,
			 "mesh",  WLZEFF_FORMAT_MESH,
			 "node",  WLZEFF_FORMAT_NODEELE,
			 "obj",   WLZEFF_FORMAT_OBJ,
			 "pgm",   WLZEFF_FORMAT_PNM,
			 "pic",   WLZEFF_FORMAT_PIC,
			 "ply2",  WLZEFF_FORMAT_PLY2,
			 "pnm",   WLZEFF_FORMAT_PNM,
			 "raw",   WLZEFF_FORMAT_RAW,
			 "slc",   WLZEFF_FORMAT_SLC,
			 "tif",   WLZEFF_FORMAT_TIFF,
			 "txt",    WLZEFF_FORMAT_TXT,
			 "vmesh", WLZEFF_FORMAT_VMESH,
			 "vff",   WLZEFF_FORMAT_VFF,
			 "vtk",   WLZEFF_FORMAT_VTK,
			 "wlz",   WLZEFF_FORMAT_WLZ,
			 NULL) == 0)
  {
    fileFmt = (unsigned int )WLZEFF_FORMAT_NONE;
  }
  return((WlzEffFormat )fileFmt);
}

/*!
* \return	File format, WLZEFF_FORMAT_NONE if string not matched.
* \ingroup	WlzExtFF
* \brief	Returns a file format for the given file format string.
* \param	fmtStr			Given file format string.
*/
WlzEffFormat	WlzEffStringToFormat(const char *fmtStr)
{
  unsigned int	fileFmt;

  if(WlzStringMatchValue((int *)&fileFmt, fmtStr,
		     "Amira Lattice", WLZEFF_FORMAT_AM,
		     "ANALYZE 7.5", WLZEFF_FORMAT_ANL,
		     "BioRad Confocal", WLZEFF_FORMAT_PIC,
		     "Graphics Interchange Format", WLZEFF_FORMAT_GIF,
		     "GRUMMP VMESH", WLZEFF_FORMAT_VMESH,
		     "Image Cytometry Standard", WLZEFF_FORMAT_ICS,
		     "IPLab", WLZEFF_FORMAT_IPL,
		     "JPEG", WLZEFF_FORMAT_JPEG,
		     "Microsoft Bitmap", WLZEFF_FORMAT_BMP,
		     "NETGEN tetrahedral mesh", WLZEFF_FORMAT_MESH,
		     "Jonathan Shewchuk's mesh format", WLZEFF_FORMAT_NODEELE,
		     "OBJ", WLZEFF_FORMAT_OBJ,
		     "PNM", WLZEFF_FORMAT_PNM,
		     "Raw", WLZEFF_FORMAT_RAW,
		     "Riken PLY2", WLZEFF_FORMAT_PLY2,
		     "SLC", WLZEFF_FORMAT_SLC,
		     "Stanford Density", WLZEFF_FORMAT_DEN,
		     "Sunvision VFF", WLZEFF_FORMAT_VFF,
		     "Text", WLZEFF_FORMAT_TXT,
		     "Tiff", WLZEFF_FORMAT_TIFF,
		     "Visualization Toolkit", WLZEFF_FORMAT_VTK,
		     "Woolz", WLZEFF_FORMAT_WLZ,
		     NULL) == 0)
  {
    fileFmt = (unsigned int )WLZEFF_FORMAT_NONE;
  }
  return((WlzEffFormat )fileFmt);
}

/*!
* \return	String appropriate to format, or NULL on error.
* \ingroup	WlzExtFF
* \brief	Returns a string which is suitable for displaying the given
*		format and sets a destination pointer for the formats file
*		extension.
* \param	fileFmt			Given file format.
* \param	dstExtStr		Destination ptr for extension string.
*/
const char	*WlzEffStringFromFormat(WlzEffFormat fileFmt,
					const char **dstExtStr)
{
  const char	*extStr = NULL,
  		*fmtStr = NULL;
  const char	*extAmStr   = "am",
                *fmtAmStr   = "Amira Lattice",
		*extAnlStr  = "hdr",
		*fmtAnlStr  = "ANALYZE HDR",
  		*extBmpStr  = "bmp",
		*fmtBmpStr  = "Microsoft Bitmap",
  		*extDenStr  = "den",
  		*fmtDenStr  = "Stanford Density",
		*extGifStr  = "gif",
		*fmtGifStr  = "Graphics Interchange Format",
  		*extIcsStr  = "ics",
  		*fmtIcsStr  = "Image Cytometry Standard",
		*extIPLStr  = "ipl",
		*fmtIPLStr  = "IPLab",
		*extJpegStr = "jpg",
		*fmtJpegStr = "JPEG",
		*extMeshStr = "mesh",
		*fmtMeshStr = "NETGEN tetrahedral mesh",
		*extNodeEleStr = "node",
		*fmtNodeEleStr = "Jonathan Shewchuk's mesh format",
  		*extObjStr  = "obj",
  		*fmtObjStr  = "Wavefront",
  		*extPnmStr  = "pnm",
  		*fmtPnmStr  = "PNM",
  		*extPicStr  = "pic",
  		*fmtPicStr  = "BioRad Confocal",
  		*extPly2Str = "ply2",
  		*fmtPly2Str = "Riken PLY2",
		*extRawStr  = "raw",
		*fmtRawStr  = "Raw",
  		*extSlcStr  = "slc",
		*fmtSlcStr  = "SLC",
		*extTxtStr  = "txt",
		*fmtTxtStr  = "Text",
		*extTiffStr = "tif",
		*fmtTiffStr = "Tiff",
  		*extVffStr  = "vff",
  		*fmtVffStr  = "Sunvision VFF",
		*extVMeshStr = "vmesh",
		*fmtVMeshStr = "GRUMMP VMESH",
  		*extVtkStr  = "vtk",
  		*fmtVtkStr  = "Visualization Toolkit VTK",
		*extWlzStr  = "wlz",
		*fmtWlzStr  = "Woolz";

  switch(fileFmt)
  {
    case WLZEFF_FORMAT_BMP:
      fmtStr = fmtBmpStr;
      extStr = extBmpStr;
      break;
    case WLZEFF_FORMAT_DEN:
      fmtStr = fmtDenStr;
      extStr = extDenStr;
      break;
    case WLZEFF_FORMAT_GIF:
      fmtStr = fmtGifStr;
      extStr = extGifStr;
      break;
    case WLZEFF_FORMAT_ICS:
      fmtStr = fmtIcsStr;
      extStr = extIcsStr;
      break;
    case WLZEFF_FORMAT_OBJ:
      fmtStr = fmtObjStr;
      extStr = extObjStr;
      break;
    case WLZEFF_FORMAT_PNM:
      fmtStr = fmtPnmStr;
      extStr = extPnmStr;
      break;
    case WLZEFF_FORMAT_PIC:
      fmtStr = fmtPicStr;
      extStr = extPicStr;
      break;
    case WLZEFF_FORMAT_PLY2:
      fmtStr = fmtPly2Str;
      extStr = extPly2Str;
      break;
    case WLZEFF_FORMAT_SLC:
      fmtStr = fmtSlcStr;
      extStr = extSlcStr;
      break;
    case WLZEFF_FORMAT_VFF:
      fmtStr = fmtVffStr;
      extStr = extVffStr;
      break;
    case WLZEFF_FORMAT_VTK:
      fmtStr = fmtVtkStr;
      extStr = extVtkStr;
      break;
    case WLZEFF_FORMAT_WLZ:
      fmtStr = fmtWlzStr;
      extStr = extWlzStr;
      break;
    case WLZEFF_FORMAT_IPL:
      fmtStr = fmtIPLStr;
      extStr = extIPLStr;
      break;
    case WLZEFF_FORMAT_TIFF:
      fmtStr = fmtTiffStr;
      extStr = extTiffStr;
      break;
    case WLZEFF_FORMAT_RAW:
      fmtStr = fmtRawStr;
      extStr = extRawStr;
      break;
    case WLZEFF_FORMAT_AM:
      fmtStr = fmtAmStr;
      extStr = extAmStr;
      break;
    case WLZEFF_FORMAT_JPEG:
      fmtStr = fmtJpegStr;
      extStr = extJpegStr;
      break;
    case WLZEFF_FORMAT_ANL:
      fmtStr = fmtAnlStr;
      extStr = extAnlStr;
      break;
    case WLZEFF_FORMAT_NODEELE:
      fmtStr = fmtNodeEleStr;
      extStr = extNodeEleStr;
      break;
    case WLZEFF_FORMAT_MESH:
      fmtStr = fmtMeshStr;
      extStr = extMeshStr;
      break;
    case WLZEFF_FORMAT_VMESH:
      fmtStr = fmtVMeshStr;
      extStr = extVMeshStr;
      break;
    case WLZEFF_FORMAT_TXT:
      fmtStr = fmtTxtStr;
      extStr = extTxtStr;
      break;
    default:
      break;
  }
  if(dstExtStr)
  {
    *dstExtStr = extStr;
  }
  return(fmtStr);
}

/*!
* \return	File format.
* \ingroup	WlzExtFF
* \brief	Returns a the file format implied by the file extension used
*		in the given file name. If no match is found
*		WlzEffFormat::WLZEFF_FORMAT_NONE is returned.
* \param	fNameStr		Given file name.
*/
WlzEffFormat 	WlzEffStringFormatFromFileName(const char *fNameStr)
{
  int           len;
  char          *dot,
                *ext,
                *sep;
  WlzEffFormat  fmt = WLZEFF_FORMAT_NONE;

  if(((len = strlen(fNameStr)) >= 3) &&  /* Minimum of 3 chars. */
     (*(fNameStr + len - 1) != '/'))     /* Directory not plain file. */
  {
    sep = strrchr(fNameStr, '/');
    if(((dot = strrchr(fNameStr, '.')) != NULL) &&
       ((sep == NULL) || ((dot - sep) > 1)))
    {
      ext = ++dot;
      if(ext != NULL)
      {
        fmt = WlzEffStringExtToFormat(ext);
      }
    }
  }
  return(fmt);
}

/*!
* \return	Woolz object read from file, or NULL on error.
* \ingroup	WlzExtFF
* \brief	Reads a woolz object from the the given file names(s)
*		or opened file, where the file format is that given as a
*		parameter.
*		The file ANL, ICS and PNM formats require a file name, for all
*		other formats either the file name or a file pointer may be
*		given.
* \param	fP			Input file.
* \param	fName			File name (and path).
* \param	fFmt			File format.
* \param	split			If non zero and the file contains a
*					index labeled image then the image
*					is split into seperate domains.
* \param	dstErr			Destination error number ptr, may be
* 					NULL.
*/
WlzObject	*WlzEffReadObj(FILE *fP, const char *fName, WlzEffFormat fFmt,
			       int split, WlzErrorNum *dstErr)
{
  int		openFileFlag = 0;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

#ifdef _WIN32
  if(fP != NULL)
  {
    if(_setmode(_fileno(fP), 0x8000) == -1)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
  }
#endif
  if((fP == NULL) && (fName == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(fP == NULL)
  {
    switch(fFmt)
    {
      case WLZEFF_FORMAT_ANL:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_BMP:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_ICS:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_NODEELE: /* FALLTHROUGH */
      case WLZEFF_FORMAT_PNM:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_TIFF:
        break;
      default:
	if((fP = fopen(fName, "rb")) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  openFileFlag = 1;
	}
	break;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    switch(fFmt)
    {
      case WLZEFF_FORMAT_BMP:
	obj = WlzEffReadObjBmp(fName, &errNum);
	break;
      case WLZEFF_FORMAT_DEN:
	obj = WlzEffReadObjDen(fP, &errNum);
	break;
      case WLZEFF_FORMAT_GIF:
        obj = WlzEffReadObjGif(fP, &errNum);
	break;
      case WLZEFF_FORMAT_ICS:
	obj = WlzEffReadObjIcs(fName, &errNum);
	break;
      case WLZEFF_FORMAT_OBJ:
	obj = WlzEffReadObjObj(fP, &errNum);
	break;
      case WLZEFF_FORMAT_PNM:
	obj = WlzEffReadObjPnm(fName, &errNum);
	break;
      case WLZEFF_FORMAT_PIC:
	obj = WlzEffReadObjPic(fP, &errNum);
	break;
      case WLZEFF_FORMAT_PLY2:
	obj = WlzEffReadObjPly2(fP, &errNum);
	break;
      case WLZEFF_FORMAT_SLC:
	obj = WlzEffReadObjSlc(fP, &errNum);
	break;
      case WLZEFF_FORMAT_VFF:
	obj = WlzEffReadObjVff(fP, &errNum);
	break;
      case WLZEFF_FORMAT_NODEELE:
	obj = WlzEffReadObjNodeEle(fName, &errNum);
	break;
      case WLZEFF_FORMAT_MESH:
	obj = WlzEffReadObjMesh(fP, &errNum);
	break;
      case WLZEFF_FORMAT_VMESH:
	obj = WlzEffReadObjVMesh(fP, &errNum);
	break;
      case WLZEFF_FORMAT_VTK:
	obj = WlzEffReadObjVtk(fP, &errNum);
	break;
      case WLZEFF_FORMAT_WLZ:
        obj = WlzReadObj(fP, &errNum);
	break;
      case WLZEFF_FORMAT_IPL:
        obj = WlzEffReadObjIPL(fP, &errNum);
	break;
      case WLZEFF_FORMAT_TIFF:
	obj = WlzEffReadObjTiff(fName, split, &errNum);
	break;
      case WLZEFF_FORMAT_AM:
        obj = WlzEffReadObjAm(fP, split, &errNum);
	break;
      case WLZEFF_FORMAT_JPEG:
        obj = WlzEffReadObjJpeg(fP, &errNum);
	break;
      case WLZEFF_FORMAT_ANL:
        obj = WlzEffReadObjAnl(fName, &errNum);
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(openFileFlag)
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
* \brief	Writes a woolz object to the the given file stream, with the
* 		given file format.
* \param	fP			Output file stream.
* \param	fName			File name (and path).
* \param	obj			Given Woolz object.
* \param	fFmt			File format.
*/
WlzErrorNum	WlzEffWriteObj(FILE *fP, const char *fName, WlzObject *obj,
			       WlzEffFormat fFmt)
{
  int		openFileFlag = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP == NULL) && (fName == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(fP == NULL)
  {
    switch(fFmt)
    {
      case WLZEFF_FORMAT_ANL:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_BMP:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_ICS:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_PNM:     /* FALLTHROUGH */
      case WLZEFF_FORMAT_NODEELE: /* FALLTHROUGH */
      case WLZEFF_FORMAT_TIFF:
        break;
      default:
	if((fP = fopen(fName, "wb")) == NULL)
	{
	  errNum = WLZ_ERR_WRITE_EOF;
	}
	else
	{
	  openFileFlag = 1;
	}
	break;
    }
  }
#ifdef _WIN32
  if(fP != NULL)
  {
    if(_setmode(_fileno(fP), 0x8000) == -1)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
  }
#endif
  if(errNum == WLZ_ERR_NONE)
  {
    switch(fFmt)
    {
      case WLZEFF_FORMAT_BMP:
	errNum = WlzEffWriteObjBmp(fName, obj);
	break;
      case WLZEFF_FORMAT_DEN:
	errNum = WlzEffWriteObjDen(fP, obj);
	break;
      case WLZEFF_FORMAT_GIF:
        errNum = WLZ_ERR_UNIMPLEMENTED;
	break;
      case WLZEFF_FORMAT_ICS:
	errNum = WlzEffWriteObjIcs(fName, obj);
	break;
      case WLZEFF_FORMAT_OBJ:
	errNum = WlzEffWriteObjObj(fP, obj);
	break;
      case WLZEFF_FORMAT_PNM:
	errNum = WlzEffWriteObjPnm(fName, obj);
	break;
      case WLZEFF_FORMAT_PIC:
	errNum = WlzEffWriteObjPic(fP, obj);
	break;
      case WLZEFF_FORMAT_PLY2:
	errNum = WlzEffWriteObjPly2(fP, obj);
	break;
      case WLZEFF_FORMAT_SLC:
	errNum = WlzEffWriteObjSlc(fP, obj);
	break;
      case WLZEFF_FORMAT_VFF:
	errNum = WlzEffWriteObjVff(fP, obj);
	break;
      case WLZEFF_FORMAT_NODEELE:
	errNum = WlzEffWriteObjNodeEle(fName, obj);
	break;
      case WLZEFF_FORMAT_MESH:
	errNum = WlzEffWriteObjMesh(fP, obj);
	break;
      case WLZEFF_FORMAT_VMESH:
	errNum = WlzEffWriteObjVMesh(fP, obj);
	break;
      case WLZEFF_FORMAT_VTK:
	errNum = WlzEffWriteObjVtk(fP, obj);
	break;
      case WLZEFF_FORMAT_WLZ:
        errNum = WlzWriteObj(fP, obj);
	break;
      case WLZEFF_FORMAT_IPL:
	errNum = WlzEffWriteObjIPL(fP, obj);
	break;
      case WLZEFF_FORMAT_TIFF:
	errNum = WlzEffWriteObjTiff(fName, obj);
	break;
      case WLZEFF_FORMAT_AM:
        errNum = WlzEffWriteObjAm(fP, obj);
	break;
      case WLZEFF_FORMAT_JPEG:
        errNum = WlzEffWriteObjJpeg(fP, obj, "");
	break;
      case WLZEFF_FORMAT_TXT:
        errNum = WlzEffWriteObjTxt(fP, obj, "");
	break;
      case WLZEFF_FORMAT_ANL:
	errNum = WlzEffWriteObjAnl(fName, obj);
	break;
      default:
        errNum = WLZ_ERR_PARAM_DATA;
	break;
    }
  }
  if(openFileFlag)
  {
    fclose(fP);
  }
  return(errNum);
}

/*!
* \return	Number of file formats.
* \ingroup	WlzExtFF
* \brief	Retuns the number of file formats, ie WLZEFF_FORMAT_COUNT - 1.
*/
int 		WlzEffNumberOfFormats(void)
{
  return WLZEFF_FORMAT_COUNT - 1;
}

/*!
* \return	Ascii table of valid file formats.
* \ingroup	WlzExtFF
* \brief	Creates an ascii table of the valid formats supported by
* 		WlzExtFF. This may be useful for giving usage information
* 		in binaries.
* \param        indWth			Indent width (no default).
* \param	desWth			Maximum description width (if 0 a
* 					default field width is used).
* \param	extWth			Maximum extension width (if 0 a default
* 					field width is used).
* \param	dstErr			Destination error pointer, may be NULL.
*/
char		*WlzEffFormatTable(unsigned indWth,
				   unsigned desWth, unsigned extWth,
				   WlzErrorNum *dstErr)
{
  int		fmtTblSz,
		maxFmtRecSz;
  int		*idxTbl;
  char		*fmtTbl = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const int	defExtWth = 10,
  		defDesWth = 50;
  
  if(desWth == 0)
  {
    desWth = defDesWth;
  }
  if(extWth == 0)
  {
    extWth = defExtWth;
  }
  maxFmtRecSz = indWth + desWth + extWth + 2;
  fmtTblSz = WLZEFF_FORMAT_COUNT * maxFmtRecSz;
  if(((idxTbl = (int *)AlcMalloc(sizeof(int) * WLZEFF_FORMAT_COUNT)) == NULL) ||
     ((fmtTbl = (char *)AlcMalloc(sizeof(char) * fmtTblSz)) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else
  {
    int		id0,
    		id1;
    char	*s0;

    for(id0 = 0; id0 < WLZEFF_FORMAT_COUNT - 1; ++id0)
    {
      idxTbl[id0] = id0;
    }
    /* Insertion sort of the index table. */
    for(id0 = 0; id0 < WLZEFF_FORMAT_COUNT - 2; ++id0)
    {
      const char *e0,
		 *e1;

      (void )WlzEffStringFromFormat((WlzEffFormat)idxTbl[id0] + 1, &e0);
      for(id1 = id0 + 1; id1 < WLZEFF_FORMAT_COUNT - 1; ++id1)
      {
        (void )WlzEffStringFromFormat((WlzEffFormat)idxTbl[id1] + 1, &e1);
	if(strcmp(e0, e1) > 0)
	{
	  int	t;

	  t = idxTbl[id0]; idxTbl[id0] = idxTbl[id1]; idxTbl[id1] = t;
	  e0 = e1;
	}
      }
    }
    /* Fill in the format table string. */
    s0 = fmtTbl;
    for(id0 = 0; id0 < WLZEFF_FORMAT_COUNT - 2; ++id0)
    {
      int 	len;
      const char *ext,
      		 *des;

      des = WlzEffStringFromFormat((WlzEffFormat )(idxTbl[id0] + 1), &ext);
      for(id1 = 0; id1 < indWth; ++id1)
      {
        *s0++ = ' ';
      }
      (void )strncpy(s0, des, desWth);
      if((len = strlen(des)) > desWth)
      {
        len = desWth;
      }
      s0 += len;
      for(id1 = len; id1 < desWth; ++id1)
      {
        *s0++ = ' ';
      }
      (void )strncpy(s0, ext, extWth);
      if((len = strlen(ext)) > extWth)
      {
        len = extWth;
      }
      s0 += len;
      *s0++ = '\n';
    }
    *s0 = '\0';
  }
  AlcFree(idxTbl);
  return(fmtTbl);
}
