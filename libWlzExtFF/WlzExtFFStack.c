#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFStack.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for reading and writting Woolz objects to
*		and from a stack of 2D files.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* 16-06-2000 bill Fixed freeing unallocated memory bug in
*		WlzEffReadObjStack3D.
* GFeng add more delimited to the token.
************************************************************************/
 #include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

static WlzObject *WlzEffReadObjStack3D(const char *, WlzEffFormat,
				WlzErrorNum *),
		*WlzEffReadObjStack2D(FILE *, WlzEffFormat,
				WlzErrorNum *);
static WlzErrorNum WlzEffHeadParseStackCtr(WlzEffStackCtrHeader *,
				FILE *),
		WlzEffStackFileNameParse(const char *, WlzEffFormat,
				 char **, char **, char **, char **),
		WlzEffReadObjStackData2D(FILE *, WlzEffFormat,
				WlzIVertex2 *, unsigned char ***),
		WlzEffWriteObjStack3D(const char *, const char *,
				const char *, const char *,
				WlzEffFormat, WlzObject *),
		WlzEffWriteObjStack2D(const char *, WlzEffFormat,
				WlzObject *, WlzIVertex2, WlzIVertex2,
				unsigned char *, unsigned char);

/************************************************************************
* Function:	WlzEffReadObjStack					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a Woolz object from the given file(s) using	*
*		the given (2D) file format.				*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzEffReadObjStack(const char *gvnFileName, WlzEffFormat fFmt,
				    WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  FILE		*fP = NULL;

  if(gvnFileName == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    if((fP = fopen(gvnFileName, "r")) == NULL)
    {
      obj = WlzEffReadObjStack3D(gvnFileName, fFmt, &errNum);
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

      obj = WlzEffReadObjStack2D(fP, fFmt, &errNum);
    }
    if(fP)
    {
      fclose(fP);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
* Function:	WlzEffWriteObjStack					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given Woolz object to the given file(s)  	*
*		using the given (2D) file format.			*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
WlzErrorNum	WlzEffWriteObjStack(const char *gvnFileName, WlzEffFormat fFmt,
				    WlzObject *obj)
{
  int		tI0;
  char		*fPathStr = NULL,
		*fBodyStr = NULL,
  		*fExtStr = NULL,
		*fCtrStr = NULL;
  unsigned char	*data = NULL;
  WlzPixelV	bgdV;
  WlzIVertex2	imgSz2D,
  		imgOrg2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(gvnFileName == NULL)
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
  else if(obj->values.core == NULL)
  {
     errNum = WLZ_ERR_VALUES_NULL;
  }
  else
  {
    switch(obj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	if((obj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
	   (obj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else
	{
	  imgOrg2D.vtX = obj->domain.i->kol1;
	  imgOrg2D.vtY = obj->domain.i->line1;
	  imgSz2D.vtX = obj->domain.i->lastkl - imgOrg2D.vtX + 1;
	  imgSz2D.vtY = obj->domain.i->lastln - imgOrg2D.vtY + 1;
	  tI0 = imgSz2D.vtX * imgSz2D.vtY;
	  if((data = (unsigned char *)AlcMalloc(tI0 *
					       sizeof(unsigned char))) == NULL)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  bgdV = WlzGetBackground(obj, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
          (void )WlzValueConvertPixel(&bgdV, bgdV, WLZ_GREY_UBYTE);
	  errNum = WlzEffWriteObjStack2D(gvnFileName, fFmt,
	  				 obj, imgSz2D, imgOrg2D,
	  			         data, bgdV.v.ubv);
	}
        break;
      case WLZ_3D_DOMAINOBJ:
	if(obj->domain.core->type != WLZ_PLANEDOMAIN_DOMAIN)
	{
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	}
	else if(obj->values.core->type != WLZ_VOXELVALUETABLE_GREY)
	{
	  errNum = WLZ_ERR_VALUES_TYPE;
	}
	else
	{
	  errNum = WlzEffStackFileNameParse(gvnFileName, fFmt,
	  				    &fPathStr, &fBodyStr,
					    &fExtStr, &fCtrStr);
        }
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzEffWriteObjStack3D(fPathStr, fBodyStr, fExtStr, fCtrStr,
	  			         fFmt, obj);
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
  }
  if(data)
  {
    AlcFree(data);
  }
  if(fPathStr)
  {
    AlcFree(fPathStr);
  }
  if(fBodyStr)
  {
    AlcFree(fBodyStr);
  }
  if(fExtStr)
  {
    AlcFree(fExtStr);
  }
  if(fCtrStr)
  {
    AlcFree(fCtrStr);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffWriteObjStack2D					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given 2D Woolz object to the given file  	*
*		using the given (2D) file format.			*
* Global refs:	-							*
* Parameters:	const char *fNameStr:	Given file name.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
*		WlzObject *obj:		Given woolz object.		*
*		WlzIVertex2 imgSz:	Required image size.		*
*		WlzIVertex2 imgOrg:	Required image origin.		*
*		unsigned char *data:	Buffer of imgSz bytes.		*
*		unsigned char bgd:	Background value.		*
************************************************************************/
static WlzErrorNum WlzEffWriteObjStack2D(const char *fNameStr,
					 WlzEffFormat fFmt, WlzObject *obj,
					 WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
					 unsigned char *data,
					 unsigned char bgd)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(fFmt)
  {
    case WLZEFF_FORMAT_BMP:
      errNum = WlzEffWriteObjBmp2D(fNameStr, obj, imgSz, imgOrg, data, bgd);
      break;
    case WLZEFF_FORMAT_PNM:
      errNum = WlzEffWriteObjPnm2D(fNameStr, obj, imgSz, imgOrg, data, bgd);
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffWriteObjStack3D					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given 3D Woolz object to the given file(s)  	*
*		using the given (2D) file format.			*
* Global refs:	-							*
* Parameters:	const char *fPathStr:	File name base.			*
*		const char *fBodyStr:	File body.			*
*		const char *fExtStr:	File extension.			*
*		const char *fCtrStr:	File control string (ctr).	*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
static WlzErrorNum WlzEffWriteObjStack3D(const char *fPathStr,
					 const char *fBodyStr,
					 const char *fExtStr,
					 const char *fCtrStr,
					 WlzEffFormat fFmt,
					 WlzObject *obj)
{
  int		tI0,
		tI1,
		planeIdx,
  		planeOff;
  WlzObject	*obj2D = NULL;
  WlzIVertex2	imgSz2D,
  		imgOrg2D;
  WlzDomain	dom2D;
  WlzValues	val2D;
  WlzPixelV	bgdV;
  unsigned char	*data = NULL;
  char		*fNameStr = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  imgOrg2D.vtX = obj->domain.p->kol1;
  imgOrg2D.vtY = obj->domain.p->line1;
  imgSz2D.vtX = obj->domain.p->lastkl - imgOrg2D.vtX + 1;
  imgSz2D.vtY = obj->domain.p->lastln - imgOrg2D.vtY + 1;
  tI0 = strlen(fPathStr) + strlen(fBodyStr) + strlen(fExtStr) + 16;
  tI1 = imgSz2D.vtX * imgSz2D.vtY;
  if(((fNameStr = AlcMalloc(tI0 * sizeof(char))) == NULL) ||
     ((data = (unsigned char *)AlcMalloc(tI1 *
     				         sizeof(unsigned char))) == NULL))
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    sprintf(fNameStr, "%s%s.%s", fPathStr, fBodyStr, fCtrStr);
    if((fP = fopen(fNameStr, "w")) == NULL)
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

		if(fprintf(fP, "%s%s%s%s%s\n"
                 "%s%s%d %d %d\n"
		 "%s%s%d %d %d\n"
		 "%s%s%g %g %g\n"
		 "%s\n",
		 WLZEFF_STACK_CTR_IDENT, WLZEFF_STACK_CTR_FIELDSEP,
		 WLZEFF_STACK_CTR_IDENTSTR, WLZEFF_STACK_CTR_FIELDSEP,
		 WLZEFF_STACK_CTR_IDENTVERSION,
      		 WLZEFF_STACK_CTR_VOLORIGIN, WLZEFF_STACK_CTR_FIELDSEP,
		 obj->domain.p->kol1, obj->domain.p->line1,
		 obj->domain.p->plane1,
		 WLZEFF_STACK_CTR_VOLSIZE, WLZEFF_STACK_CTR_FIELDSEP,
		 obj->domain.p->lastkl - obj->domain.p->kol1 + 1,
		 obj->domain.p->lastln - obj->domain.p->line1 + 1,
		 obj->domain.p->lastpl - obj->domain.p->plane1 + 1,
		 WLZEFF_STACK_CTR_VOXSIZE, WLZEFF_STACK_CTR_FIELDSEP,
		 obj->domain.p->voxel_size[0],
		 obj->domain.p->voxel_size[1],
		 obj->domain.p->voxel_size[2],
		 WLZEFF_STACK_CTR_FILES) < 16)
      {
        errNum = WLZ_ERR_WRITE_INCOMPLETE;
      }
      else
      {
        planeOff = 0;
	planeIdx = obj->domain.p->plane1;
	while((planeIdx <= obj->domain.p->lastpl) && (errNum == WLZ_ERR_NONE))
	{
	  if(fprintf(fP, "%*d %s%0*d.%s\n",
		     WLZEFF_STACK_NAMEDIGITS, planeIdx,
	  	     fBodyStr, WLZEFF_STACK_NAMEDIGITS, planeIdx,
		     fExtStr) < WLZEFF_STACK_NAMEDIGITS)
	  {
	    errNum = WLZ_ERR_WRITE_INCOMPLETE;
	  }
	  ++planeIdx;
	}
      }
      fclose(fP);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    planeOff = 0;
    planeIdx = obj->domain.p->plane1;
    bgdV = WlzGetBackground(obj, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    (void )WlzValueConvertPixel(&(bgdV), bgdV, WLZ_GREY_UBYTE);
    while((planeIdx <= obj->domain.p->lastpl) && (errNum == WLZ_ERR_NONE))
    {
      dom2D = *(obj->domain.p->domains + planeOff);
      val2D = *(obj->values.vox->values + planeOff);
      sprintf(fNameStr, "%s%s%0*d.%s",
      	      fPathStr, fBodyStr, WLZEFF_STACK_NAMEDIGITS, planeIdx, fExtStr);
      if(dom2D.core && val2D.core)
      {
	obj2D = WlzAssignObject(WlzMakeMain(WLZ_2D_DOMAINOBJ, dom2D, val2D,
					 NULL, NULL, &errNum), NULL);
      }
      else
      {
	obj2D = NULL;
      }
      if(errNum == WLZ_ERR_NONE)
      {
	errNum = WlzEffWriteObjStack2D(fNameStr, fFmt,
				       obj2D, imgSz2D, imgOrg2D,
				       data, bgdV.v.ubv);
      }
      if(obj2D)
      {
	WlzFreeObj(obj2D);
	obj2D = NULL;
      }
      ++planeOff;
      ++planeIdx;
    }
  }
  if(fNameStr)
  {
    AlcFree(fNameStr);
  }
  if(data)
  {
    AlcFree(data);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffReadObjStack2D					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a 2D Woolz object from the given file using	*
*		the given (2D) file format.				*
* Global refs:	-							*
* Parameters:	FILE *fP:		Given file.	 		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
static WlzObject *WlzEffReadObjStack2D(FILE *fP, WlzEffFormat fFmt,
				       WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  unsigned char **data = NULL;
  WlzIVertex2	imgSz,
  		imgOrg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  imgSz.vtX = 0;
  imgSz.vtY = 0;
  imgOrg.vtX = 0;
  imgOrg.vtY = 0;
  errNum = WlzEffReadObjStackData2D(fP, fFmt, &imgSz, &data);
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzFromArray2D((void **)data, imgSz, imgOrg,
    			 WLZ_GREY_UBYTE, WLZ_GREY_UBYTE,
			 0.0, 1.0, 0, 1, &errNum);
  }
  if(data)
  {
    if(obj)
    {
      AlcFree(data); 			    /* **data used by object */
    }
    else
    {
      Alc2Free((void **)data);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
* Function:	WlzEffReadObjStackData2D				*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Reads the  data from a file with the given (2D) format	*
*		into the data array given.				*
* Global refs:	-							*
* Parameters:	FILE *fP:		Given file stream.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
*		WlzIVertex2 *imgSz:	Dst ptr for image size which	*
*					is assumed valid if height and	*
*					width are non zero.		*
*		unsigned char ***data:	Ptr to 2D array for data.	*
************************************************************************/
static WlzErrorNum WlzEffReadObjStackData2D(FILE *fP, WlzEffFormat fFmt,
					    WlzIVertex2 *imgSz,
					    unsigned char ***data)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(fFmt)
  {
    case WLZEFF_FORMAT_BMP:
      errNum = WlzEffReadObjBmpData2D(fP, imgSz, data);
      break;
    case WLZEFF_FORMAT_PNM:
      errNum = WlzEffReadObjPnmData2D(fP, imgSz, data);
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffReadObjStack3D					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a 3D Woolz object from the given file(s) using	*
*		the given (2D) file format.				*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
static WlzObject *WlzEffReadObjStack3D(const char *gvnFileName,
				       WlzEffFormat fFmt,
				       WlzErrorNum *dstErr)
{
  int		tI0,
		planeIdx,
		planeOff;
  char		*tCP0,
  		*recTok,
  		*fNameStr = NULL,
  		*fPathStr = NULL,
  		*fBodyStr = NULL,
  		*fExtStr = NULL,
		*fCtrStr = NULL;
  FILE		*fP = NULL,
  		*fP2D = NULL;
  WlzIVertex2	imgSz2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  unsigned char	***data = NULL;
  char		fRecord[WLZEFF_STACK_CTR_RECORDMAX];
  WlzEffStackCtrHeader header;

  errNum = WlzEffStackFileNameParse(gvnFileName, fFmt, &fPathStr, &fBodyStr,
				    &fExtStr, &fCtrStr);
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = strlen(fPathStr) + strlen(fBodyStr) + strlen(fExtStr) + 16;
    if((fNameStr = AlcMalloc(tI0 * sizeof(char))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      sprintf(fNameStr, "%s%s.%s", fPathStr, fBodyStr, fCtrStr);
      if((fP = fopen(fNameStr, "r")) == NULL)
      {
        errNum = WLZ_ERR_READ_EOF;
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

        errNum = WlzEffHeadParseStackCtr(&header, fP);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    imgSz2D.vtX = header.volSize.vtX;
    imgSz2D.vtY = header.volSize.vtY;
    if((AlcUnchar3Malloc(&data, header.volSize.vtZ, header.volSize.vtY,
			 header.volSize.vtX) != ALC_ER_NONE) ||
       (data == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    planeOff = 0;
    planeIdx = header.volOrigin.vtZ;
    while((planeOff < header.volSize.vtZ) && (errNum == WLZ_ERR_NONE))
    {
      if(fgets(fRecord, WLZEFF_STACK_CTR_RECORDMAX, fP) == NULL)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      else
      {
	fRecord[WLZEFF_STACK_CTR_RECORDMAX - 1] = '\0';
	if((tCP0 = strrchr(fRecord,  '\n')) != NULL)
	{
	  *tCP0 = '\0';
	}
	if(((recTok = strtok(fRecord, " \t")) == NULL) ||
	   (sscanf(recTok, "%d", &tI0) != 1) ||
	   (tI0 != planeIdx) ||
	   ((recTok = strtok(NULL, " \t\n")) == NULL))
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  sprintf(fNameStr, "%s%s", fPathStr, recTok);
	  if((fP2D = fopen(fNameStr, "r")) == NULL)
	  {
	    errNum = WLZ_ERR_READ_EOF;
	  }
	  else
	  {
		    #ifdef _WIN32
  if (fP2D != NULL){
	if(_setmode(_fileno(fP2D), 0x8000) == -1)
	{
		errNum = WLZ_ERR_READ_EOF;
	}
  }
  #endif
	    errNum = WlzEffReadObjStackData2D(fP2D, fFmt, &imgSz2D,
					      (data + planeOff));
	    fclose(fP2D);
	  }
	}
      }
      ++planeOff;
      ++planeIdx;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj = WlzFromArray3D((void ***)data, header.volSize, header.volOrigin,
    			 WLZ_GREY_UBYTE, WLZ_GREY_UBYTE,
			 0.0, 1.0, 0, 1, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj->domain.p->voxel_size[0] = header.voxSize.vtX;
    obj->domain.p->voxel_size[1] = header.voxSize.vtY;
    obj->domain.p->voxel_size[2] = header.voxSize.vtZ;
  }
  if(data)
  {
    if(obj)
    {
      /* **data used in object */
      AlcFree(*data);
      AlcFree(data);
    }
    else
    {
      if(data)
      {
	Alc3Free((void ***)data);
      }
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  if(fNameStr)
  {
    AlcFree(fNameStr);
  }
  if(fPathStr)
  {
    AlcFree(fPathStr);
  }
  if(fBodyStr)
  {
    AlcFree(fBodyStr);
  }
  if(fExtStr)
  {
    AlcFree(fExtStr);
  }
  if(fCtrStr)
  {
    AlcFree(fCtrStr);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
* Function:	WlzEffStackFileNameParse				*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Parses the given file name into it's base, body and	*
*		extension parts. If there is no extension then the	*
*		extension '???' is used. The control file extension	*
*		'.ctr' is also created.					*
*		Eg, gvnFileName = "/some/where/here.???"		*
*		    *fPathStr = "/some/where/"				*
*		    *fBodyStr = "here"					*
*		    *fExtStr = "???"					*
*		    *fCtrStr = "ctr"					*
*		Where the string ??? is: the appropriate extension for	*
*		the given (2D) file format, ie:	bmp, pgm.		*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
*		WlzEffFormat fFmt:	Given file format (must be a 2D	*
*					file format).			*
*		char **fPathStr:	File name base.			*
*		char **fBodyStr:	File name body.			*
*		char **fExtStr:		File extension (pgm).		*
*		char **fCtrStr		Control file extension (ctr).	*
************************************************************************/
static WlzErrorNum WlzEffStackFileNameParse(const char *gvnFileName,
					    WlzEffFormat fFmt,
					    char **fPathStr, char **fBodyStr,
					    char **fExtStr, char **fCtrStr)
{
  int		gvnLen,
		pathLen = 0,
  		bodyLen = 0,
  		extLen = 0;
  const char	*extStr,
  		*sepPath = NULL,
  		*sepExt = NULL,
		*pathStart = NULL,
		*bodyStart = NULL,
		*extStart = NULL;
  const char	extBmp[] = "bmp",
		extCtr[] = "ctr",
  		extPnm[] = "pgm";
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  *fPathStr = NULL;
  *fBodyStr = NULL;
  *fExtStr = NULL;
  *fCtrStr = NULL;
  switch(fFmt)
  {
    case WLZEFF_FORMAT_BMP:
      extStr = extBmp;
      break;
    case WLZEFF_FORMAT_PNM:
      extStr = extPnm;
      break;
    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((gvnLen = strlen(gvnFileName)) > 0)
    {
      sepPath = strrchr(gvnFileName, '/');
      sepExt = strrchr(gvnFileName, '.');
      if(sepExt <= sepPath)
      {
	sepExt = NULL;
      }
      if(sepPath)
      {
	pathStart = gvnFileName;
	pathLen = sepPath - gvnFileName + 1;
      }
      if(sepExt)
      {
	if(sepPath)
	{
	  bodyLen = sepExt - sepPath - 1; 
	  bodyStart = sepPath + 1;
	}
	else
	{
	  bodyLen = sepExt - gvnFileName;
	  bodyStart = gvnFileName;
	}
      }
      else
      {
	if(sepPath)
	{
	  bodyLen = gvnLen - pathLen;
	  bodyStart = sepPath + 1;
	}
	else
	{
	  bodyLen = gvnLen;
	  bodyStart = gvnFileName;
	}
      }
      if(sepExt)
      {
	extLen = gvnLen - (sepExt - gvnFileName + 1);
	if(extLen == 3)
	{
	  extStart = sepExt + 1;
	}
	else
	{
	  extLen = 0;
	}
      }
    }
    if(((*fPathStr = AlcMalloc(sizeof(char) * (pathLen + 1))) == NULL) ||
       ((*fBodyStr = AlcMalloc(sizeof(char) * (bodyLen + 1))) == NULL) ||
       ((*fExtStr = AlcMalloc(sizeof(char) * 4)) == NULL) ||
       ((*fCtrStr = AlcMalloc(sizeof(char) * 4)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(pathStart && (pathLen > 0))
    {
      strncpy(*fPathStr, pathStart, pathLen);
    }
    *(*fPathStr + pathLen) = '\0';
    if(bodyStart && (bodyLen > 0))
    {
      strncpy(*fBodyStr, bodyStart, bodyLen);
    }
    *(*fBodyStr + bodyLen) = '\0';
    if(extStart && (extLen == 3))
    {
      strncpy(*fExtStr, extStart, 3);
    }
    else
    {
      strcpy(*fExtStr, extStr);
    }
    *(*fExtStr + 3) = '\0';
    strcpy(*fCtrStr, extCtr);
    *(*fCtrStr + 3) = '\0';
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(*fBodyStr)
    {
      AlcFree(*fBodyStr);
    }
    if(*fPathStr)
    {
      AlcFree(*fPathStr);
    }
    if(*fExtStr)
    {
      AlcFree(*fExtStr);
    }
    if(*fCtrStr)
    {
      AlcFree(*fCtrStr);
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffHeadParseStackCtr					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Parses the stack control file header until the 'files'	*
*		record.							*
* Global refs:	-							*
* Parameters:	WlzEffStackCtrHeader *header: Header to fill in.	*
*		FILE *fP:		Opened control file stream.	*
************************************************************************/
static WlzErrorNum WlzEffHeadParseStackCtr(WlzEffStackCtrHeader *header,
					   FILE *fP)
{
  int		parsedFlag = 0,
  		volSzFlag = 0,
		recIdx = 0;
  char		*tCP0,
  		*recTok;
  char		fRecord[WLZEFF_STACK_CTR_RECORDMAX];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  header->volOrigin.vtX = 0;
  header->volOrigin.vtY = 0;
  header->volOrigin.vtZ = 0;
  header->voxSize.vtX = 1.0;
  header->voxSize.vtY = 1.0;
  header->voxSize.vtZ = 1.0;
  while((parsedFlag == 0) && (errNum == WLZ_ERR_NONE))
  {
    if(fgets(fRecord, WLZEFF_STACK_CTR_RECORDMAX, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      fRecord[WLZEFF_STACK_CTR_RECORDMAX - 1] = '\0';
      if((tCP0 = strrchr(fRecord,  '\n')) != NULL)
      {
	*tCP0 = '\0';
      }
      if((recTok = strtok(fRecord, WLZEFF_STACK_CTR_FIELDSEP)) != NULL)
      {
	if(recIdx == 0)
	{
	  if(strcmp(recTok, WLZEFF_STACK_CTR_IDENT) ||
	     ((recTok = strtok(NULL,
			       WLZEFF_STACK_CTR_FIELDSEP)) == NULL) ||
	     strcmp(recTok, WLZEFF_STACK_CTR_IDENTSTR))
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	}
	else
	{
	  if(strcmp(recTok, WLZEFF_STACK_CTR_VOLORIGIN) == 0)
	  {
	    if(((recTok = strtok(NULL,
				 WLZEFF_STACK_CTR_FIELDSEP)) == NULL) ||
	       (sscanf(recTok, "%d %d %d",
	               &(header->volOrigin.vtX),
		       &(header->volOrigin.vtY),
		       &(header->volOrigin.vtZ)) != 3))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	  else if(strcmp(recTok, WLZEFF_STACK_CTR_VOLSIZE) == 0)
	  {
	    if(((recTok = strtok(NULL,
				 WLZEFF_STACK_CTR_FIELDSEP)) == NULL) ||
	       (sscanf(recTok, "%d %d %d",
	       	       &(header->volSize.vtX),
		       &(header->volSize.vtY),
		       &(header->volSize.vtZ)) != 3))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      volSzFlag = 1;
	    }
	  }
	  else if(strcmp(recTok, WLZEFF_STACK_CTR_VOXSIZE) == 0)
	  {
	    if(((recTok = strtok(NULL,
				 WLZEFF_STACK_CTR_FIELDSEP)) == NULL) ||
	       (sscanf(recTok, "%g %g %g",
	       	       &(header->voxSize.vtX),
		       &(header->voxSize.vtY),
		       &(header->voxSize.vtZ)) != 3))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	  else if(strcmp(recTok, WLZEFF_STACK_CTR_FILES) == 0)
	  {
	    parsedFlag = 1;
	  }
	  else
	  {
	    while(*recTok && isspace(*recTok))
	    {
	      ++recTok;
	    }
	    if(*recTok && (strcmp(recTok, WLZEFF_STACK_CTR_COMMENT) == 0))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	  }
	}
	++recIdx;
      }
    }
  }
  return(errNum);
}
