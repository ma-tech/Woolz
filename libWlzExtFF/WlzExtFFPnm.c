#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFPnm.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Functions for reading and writting Woolz objects to
*		and from the portable anymap '.pnm' data format.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
* GFeng add more delimited to the token.
************************************************************************/
 #include <ctype.h>
#include <string.h>
#include <Wlz.h>
#include <WlzExtFF.h>

/************************************************************************
* Function:	WlzEffReadObjPnm					*
* Returns:	WlzObject *:		Object read from file.		*
* Purpose:	Reads a Woolz object from the given file(s) using	*
*		the '.pnm' file format.					*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
* 		WlzErrorNum *dstErr:	Destination error number ptr,	*
*					may be NULL.			*
************************************************************************/
WlzObject	*WlzEffReadObjPnm(const char *gvnFileName, WlzErrorNum *dstErr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;

  obj = WlzEffReadObjStack(gvnFileName, WLZEFF_FORMAT_PNM, &errNum);
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/************************************************************************
* Function:	WlzEffWriteObjPnm					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given Woolz object to the given file(s)  	*
*		using the '.pnm' file format.				*
* Global refs:	-							*
* Parameters:	const char *gvnFileName: Given file name.		*
*		WlzObject *obj:		Given woolz object.		*
************************************************************************/
WlzErrorNum	WlzEffWriteObjPnm(const char *gvnFileName, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzEffWriteObjStack(gvnFileName, WLZEFF_FORMAT_PNM, obj);
  return(errNum);
}

/************************************************************************
* Function:	WlzEffWriteObjPnm2D					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Writes the given 2D Woolz object to the given file  	*
*		using the '.pnm' file format.				*
* Global refs:	-							*
* Parameters:	const char *fNameStr:	Given file name.		*
*		WlzObject *obj:		Given woolz object.		*
*		WlzIVertex2 imgSz:	Required image size.		*
*		WlzIVertex2 imgOrg:	Required image origin.		*
*		unsigned char *data:	Buffer of imgSz bytes.		*
*		unsigned char bgd:	Background value.		*
************************************************************************/
WlzErrorNum 	WlzEffWriteObjPnm2D(const char *fNameStr, WlzObject *obj,
				    WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
				    unsigned char *data,
				    unsigned char bgd)
{
  int		tI0;
  FILE		*fP = NULL;
  WlzIBox2	cutBox;
  WlzObject	*cutObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

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

    if((obj == NULL) || (obj->values.core == NULL))
    {
      (void )memset(data, (unsigned int )bgd, imgSz.vtX * imgSz.vtY);
    }
    else
    {
      cutBox.xMin = imgOrg.vtX;
      cutBox.yMin = imgOrg.vtY;
      cutBox.xMax = imgOrg.vtX + imgSz.vtX - 1;
      cutBox.yMax = imgOrg.vtY + imgSz.vtY - 1;
      cutObj = WlzCutObjToValBox2D(obj, cutBox, WLZ_GREY_UBYTE, data,
      			        0, 0.0, 0.0, &errNum);
      if(cutObj)
      {
	if(cutObj->type == WLZ_2D_DOMAINOBJ)
	{
	  /* (void )WlzPopFreePtr(cutObj->values.r->freeptr, NULL, NULL); */
	  AlcFree(cutObj->values.r->freeptr);
	  cutObj->values.r->freeptr = NULL;
	  cutObj->values.r->values.inp = NULL;
	}
	WlzFreeObj(cutObj);
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = imgSz.vtX * imgSz.vtY;
    if((fprintf(fP, "%s\n%d %d\n255\n",
    		WLZEFF_PGM_MAGIC, imgSz.vtX, imgSz.vtY) < 8) ||
       (fwrite(data, sizeof(unsigned char), tI0, fP) !=  tI0))

    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  if(fP)
  {
    fclose(fP);
  }
  return(errNum);
}

/************************************************************************
* Function:	WlzEffReadObjPnmData2D					*
* Returns:	WlzErrorNum		Woolz error number.		*
* Purpose:	Reads the  data from a pgm file into the data array	*
*		given.							*
* Global refs:	-							*
* Parameters:	FILE *fP:		Given file stream.		*
*		WlzIVertex2 gvnImgSz:	Dst ptr for image size which	*
*					is assumed valid if height and	*
*					width are non zero.		*
*		unsigned char ***data:	Ptr to 2D array for data.	*
************************************************************************/
WlzErrorNum 	WlzEffReadObjPnmData2D(FILE *fP, WlzIVertex2 *gvnImgSz,
				       unsigned char ***data)
{
  int		tI0,
  		parsedFlag = 0,
  		recIdx = 0;
  char		*tCP0,
  		*rec;
  WlzIVertex2	imgSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char		fRecord[256];

  while((parsedFlag == 0) && (errNum == WLZ_ERR_NONE))
  {
    if(fgets(fRecord, 256, fP) == NULL)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
    else
    {
      fRecord[255] = '\0';
      if((tCP0 = strrchr(fRecord,  '\n')) != NULL)
      {
	*tCP0 = '\0';
      }
      rec = fRecord;
      while(*rec && isspace(*rec))
      {
	++rec;
      }
      if(*rec && (*rec != '#'))
      {
	switch(recIdx)
	{
	  case 0:
	    if(strcmp(rec, WLZEFF_PGM_MAGIC))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    break;
	  case 1:
	    if(sscanf(rec, "%d %d", &(imgSz.vtX), &(imgSz.vtY)) != 2)
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      if((gvnImgSz->vtX && (imgSz.vtX != gvnImgSz->vtX)) ||
		 (gvnImgSz->vtY && (imgSz.vtY != gvnImgSz->vtY)))
	      {
		errNum = WLZ_ERR_READ_INCOMPLETE;
	      }
	      else
	      {
	        *gvnImgSz = imgSz;
	      }
	    }
	    break;
	  case 2:
	    if((sscanf(rec, "%d", &tI0) != 1) || (tI0 != 255))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    parsedFlag = 1;
	    break;
	}
	++recIdx;
      }
    }
  }
  if(*data == NULL)
  {
    if((AlcUnchar2Malloc(data, imgSz.vtY, imgSz.vtX) != ALC_ER_NONE) ||
       (*data == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    tI0 = gvnImgSz->vtX * gvnImgSz->vtY;
    if(fread(**data, sizeof(unsigned char), tI0, fP) != tI0)
    {
      errNum = WLZ_ERR_READ_INCOMPLETE;
    }
  }
  return(errNum);
}
