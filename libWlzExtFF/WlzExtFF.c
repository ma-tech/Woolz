#pragma ident "MRC HGU $Id$"
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
			 "bmp", WLZEFF_FORMAT_BMP,
			 "den", WLZEFF_FORMAT_DEN,
			 "ics", WLZEFF_FORMAT_ICS,
			 "ids", WLZEFF_FORMAT_ICS,
			 "pgm", WLZEFF_FORMAT_PNM,
			 "pnm", WLZEFF_FORMAT_PNM,
			 "pic", WLZEFF_FORMAT_PIC,
			 "slc", WLZEFF_FORMAT_SLC,
			 "vff", WLZEFF_FORMAT_VFF,
			 "vtk", WLZEFF_FORMAT_VTK,
			 "wlz", WLZEFF_FORMAT_WLZ,
			 "ipl", WLZEFF_FORMAT_IPL,
			 "tif", WLZEFF_FORMAT_TIFF,
			 "raw", WLZEFF_FORMAT_RAW,
			 "am",  WLZEFF_FORMAT_AM,
			 "jpg", WLZEFF_FORMAT_JPEG,
			 "jpeg", WLZEFF_FORMAT_JPEG,
			 "hdr", WLZEFF_FORMAT_ANL,
			 "img", WLZEFF_FORMAT_ANL,
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
			 "Microsoft Bitmap", WLZEFF_FORMAT_BMP,
			 "Stanford Density", WLZEFF_FORMAT_DEN,
			 "ICS", WLZEFF_FORMAT_ICS,
			 "PNM", WLZEFF_FORMAT_PNM,
			 "BioRad Confocal", WLZEFF_FORMAT_PIC,
			 "SLC", WLZEFF_FORMAT_SLC,
			 "Sunvision VFF", WLZEFF_FORMAT_VFF,
			 "Visualization Toolkit VTK", WLZEFF_FORMAT_VTK,
			 "Woolz", WLZEFF_FORMAT_WLZ,
			 "IPLab", WLZEFF_FORMAT_IPL,
			 "Tiff", WLZEFF_FORMAT_TIFF,
			 "Raw", WLZEFF_FORMAT_RAW,
			 "Amira Lattice", WLZEFF_FORMAT_AM,
			 "JPEG", WLZEFF_FORMAT_JPEG,
			 "ANALYZE 7.5", WLZEFF_FORMAT_ANL,
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
  const char	*extBmpStr = "bmp",
		*fmtBmpStr = "Microsoft Bitmap",
  		*extDenStr = "den",
  		*fmtDenStr = "Stanford Density",
  		*extIcsStr = "ics",
  		*fmtIcsStr = "ICS",
  		*extPnmStr = "pnm",
  		*fmtPnmStr = "PNM",
  		*extPicStr = "pic",
  		*fmtPicStr = "BioRad Confocal",
  		*extSlcStr = "slc",
		*fmtSlcStr = "SLC",
  		*extVffStr = "vff",
  		*fmtVffStr = "Sunvision VFF",
  		*extVtkStr = "vtk",
  		*fmtVtkStr = "Visualization Toolkit VTK",
		*extWlzStr = "wlz",
		*fmtWlzStr = "Woolz",
		*extIPLStr = "ipl",
		*fmtIPLStr = "IPLab",
		*extTiffStr = "tif",
		*fmtTiffStr = "Tiff",
		*extRawStr = "raw",
		*fmtRawStr = "Raw",
		*extAmStr  = "am",
                *fmtAmStr  = "Amira Lattice",
		*extJpegStr  = "jpg",
		*fmtJpegStr  = "JPEG",
		*extAnlStr = "hdr",
		*fmtAnlStr = "HDR";

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
    case WLZEFF_FORMAT_ICS:
      fmtStr = fmtIcsStr;
      extStr = extIcsStr;
      break;
    case WLZEFF_FORMAT_PNM:
      fmtStr = fmtPnmStr;
      extStr = extPnmStr;
      break;
    case WLZEFF_FORMAT_PIC:
      fmtStr = fmtPicStr;
      extStr = extPicStr;
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
      case WLZEFF_FORMAT_ANL: /* FALLTHROUGH */
      case WLZEFF_FORMAT_BMP: /* FALLTHROUGH */
      case WLZEFF_FORMAT_ICS: /* FALLTHROUGH */
      case WLZEFF_FORMAT_PNM: /* FALLTHROUGH */
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
      case WLZEFF_FORMAT_ICS:
	obj = WlzEffReadObjIcs(fName, &errNum);
	break;
      case WLZEFF_FORMAT_PNM:
	obj = WlzEffReadObjPnm(fName, &errNum);
	break;
      case WLZEFF_FORMAT_PIC:
	obj = WlzEffReadObjPic(fP, &errNum);
	break;
      case WLZEFF_FORMAT_SLC:
	obj = WlzEffReadObjSlc(fP, &errNum);
	break;
      case WLZEFF_FORMAT_VFF:
	obj = WlzEffReadObjVff(fP, &errNum);
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
	obj = WlzEffReadObjTiff(fName, &errNum);
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
      case WLZEFF_FORMAT_ANL: /* FALLTHROUGH */
      case WLZEFF_FORMAT_BMP: /* FALLTHROUGH */
      case WLZEFF_FORMAT_ICS: /* FALLTHROUGH */
      case WLZEFF_FORMAT_PNM: /* FALLTHROUGH */
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
      case WLZEFF_FORMAT_ICS:
	errNum = WlzEffWriteObjIcs(fName, obj);
	break;
      case WLZEFF_FORMAT_PNM:
	errNum = WlzEffWriteObjPnm(fName, obj);
	break;
      case WLZEFF_FORMAT_PIC:
	errNum = WlzEffWriteObjPic(fP, obj);
	break;
      case WLZEFF_FORMAT_SLC:
	errNum = WlzEffWriteObjSlc(fP, obj);
	break;
      case WLZEFF_FORMAT_VFF:
	errNum = WlzEffWriteObjVff(fP, obj);
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
