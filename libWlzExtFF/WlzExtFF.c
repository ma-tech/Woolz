#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFF.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for reading and writting Woolz objects to
*		and from external data formats.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <Wlz.h>
#include <WlzExtFF.h>

/************************************************************************
* Function:	WlzEffStringExtToFormat					*
* Returns:	WlzEffFormat:		File format, WLZEFF_FORMAT_NONE	*
*					if string not matched.		*
* Purpose:	Returns a file format for the given file extension 	*
*		string.							*
* Global refs:	-							*
* Parameters:	const char *extStr:	Given file extension string.	*
************************************************************************/
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
			 NULL) == 0)
  {
    fileFmt = (unsigned int )WLZEFF_FORMAT_NONE;
  }
  return((WlzEffFormat )fileFmt);
}

/************************************************************************
* Function:	WlzEffStringToFormat					*
* Returns:	WlzEffFormat:		File format, WLZEFF_FORMAT_NONE	*
*					if string not matched.		*
* Purpose:	Returns a file format for the given file format string.	*
* Global refs:	-							*
* Parameters:	const char *fmtStr:	Given file format string.	*
************************************************************************/
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
			 NULL) == 0)
  {
    fileFmt = (unsigned int )WLZEFF_FORMAT_NONE;
  }
  return((WlzEffFormat )fileFmt);
}

/************************************************************************
* Function:	WlzEffStringFromFormat					*
* Returns:	const char *:		String appropriate to format,	*
*					or NULL on error.		*
* Purpose:	Returns a string which is suitable for displaying	*
*		the given format and sets a destination pointer for the	*
*		formats file extension.					*
* Global refs:	-							*
* Parameters:	WlzEffFormat fileFmt:	Given file format.		*
*		const char **dstExtStr:	Destination ptr for extension	*
*					string.				*
************************************************************************/
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
		*fmtWlzStr = "Woolz";

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
  }
  if(dstExtStr)
  {
    *dstExtStr = extStr;
  }
  return(fmtStr);
}

/************************************************************************
* Function:	WlzEffReadObj						*
* Returns:	WlzObject:		Woolz object read from file, or	*
*					NULL on error.			*
* Purpose:	Reads a woolz object from the the given file names(s)	*
*		or opened file stream, where the file format is that	*
*		given as a parameter.					*
*		The file ICS and PNM formats require a file name, for	*
*		all other formats either the file name or a file stream	*
*		pointer may be given.					*
* Global refs:	-							*
* Parameters:	FILE *fP:		Input file stream.		*
*		const char *fName:	File name (and path).		*
*		WlzEffFormat fFmt:	File format.			*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzEffReadObj(FILE *fP, const char *fName, WlzEffFormat fFmt,
			       WlzErrorNum *dstErr)
{
  int		openFileFlag = 0;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP == NULL) && (fName == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else if(fP == NULL)
  {
    if((fFmt != WLZEFF_FORMAT_ICS) && (fFmt != WLZEFF_FORMAT_PNM) &&
       (fFmt != WLZEFF_FORMAT_BMP))
    {
      if((fP = fopen(fName, "r")) == NULL)
      {
	errNum = WLZ_ERR_READ_EOF;
      }
      else
      {
	openFileFlag = 1;
      }
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


/************************************************************************
* Function:	WlzEffWriteObj						*
* Returns:	WlzErrorNum:		Woolz error number.		*
* Purpose:	Writes a woolz object to the the given file stream,	*
*		with the given file format.				*
* Global refs:	-							*
* Parameters:	FILE *fP:		Output file stream.		*
*		const char *fName:	File name (and path).		*
*		WlzObject *obj:		Given woolz object.		*
*		WlzEffFormat fFmt:	File format.			*
************************************************************************/
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
    if((fFmt != WLZEFF_FORMAT_ICS) && (fFmt != WLZEFF_FORMAT_PNM) &&
       (fFmt != WLZEFF_FORMAT_BMP))
    {
      if((fP = fopen(fName, "w")) == NULL)
      {
	errNum = WLZ_ERR_WRITE_EOF;
      }
      else
      {
	openFileFlag = 1;
      }
    }
  }
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
