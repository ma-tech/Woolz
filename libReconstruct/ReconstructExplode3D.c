#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructExplode3D.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides functions for explodeing a 3D object into
*		a section list and 2D files.			
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>

/************************************************************************
* Function:	RecExplode3DObjToFile				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Given a destination directory and file body and a 3D
*		woolz object, the 3D object is and exploded into a
*		series of 2D section files (in the destination	
*		directory) each with a name formed from the given file
*		body and the source plane index. A section list is
*		also created in the directory using the given file body
*		together with the  '.bib' file extension.	
* Global refs:	-						
* Parameters:	char *dstDirStr:	Destination directory.	
*		char *dstBodyStr:	Destination file body.	
*		WlzObject *srcObj:	Source 3D woolz	object.	
*		char *srcFName:		Source file name, may be NULL.
*		WlzEffFormat srcFFormat: Source file format, may be
*					WLZEFF_FORMAT_NONE.
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
RecError	RecExplode3DObjToFile(char *dstDirStr, char *dstBodyStr,
				      WlzObject *srcObj, char *srcFName,
				      WlzEffFormat srcFFormat, char **eMsg)
{
  int		objIdx,
		planeIdx,
  		objVecCount,
		fileStrSpace;
  char		*fileStr = NULL,
  		*fileTemplateStr = NULL;
  FILE		*fP;
  WlzObject	**objVec = NULL;
  RecSectionList *secList = NULL;
  RecError	errFlag = REC_ERR_NONE;
  struct stat	statBuf;
  const char 	*errDirAccessStr = "Failed to access directory",
		*errDirCreateStr = "Failed to create directory",
  		*errExplodeStr = "Failed to explode 3D object",
  		*errWriteWoolzStr = "Failed to write 2D section file",
		*errListAppendStr = "Failed to build section list",
		*errWriteListStr = "Failed to write section list to file";
  const char	*errStr = NULL;

  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecExplode3DObjToFile FE 0x%lx 0x%lx 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )dstDirStr, (unsigned long )dstBodyStr,
	   (unsigned long )srcObj, (unsigned long )srcFName,
	   srcFFormat, (unsigned long )eMsg));
  if((dstDirStr == NULL) || (strlen(dstDirStr) == 0) ||
     (dstBodyStr == NULL) || (strlen(dstBodyStr) == 0) || 
     (srcObj == NULL) || (srcObj->type != WLZ_3D_DOMAINOBJ) ||
     (srcObj->domain.core == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)       /* Check for/create destination directory */
  {
    if(stat(dstDirStr, &statBuf) == 0)			 /* Directory exists */
    {
#ifdef LINUX2
      if((statBuf.st_mode & (__S_IFDIR|__S_IREAD|__S_IWRITE|__S_IEXEC)) == 0)
#else /* LINUX2 */
      if((statBuf.st_mode & (S_IFDIR | S_IRWXU)) == 0)
#endif /* LINUX2 */
      {
	errStr = errDirAccessStr;
	errFlag = REC_ERR_WRITE;         /* Can't read/write/search directory */
      }
    }
    else
    {
      if(mkdir(dstDirStr, S_IRWXU))
      {
	errStr = errDirCreateStr;
	errFlag = REC_ERR_WRITE;
      }
      
    }
    REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	    ("RecExplode3DObjToFile 01 %d\n",
	     (int )errFlag));
  }
  if(errFlag == REC_ERR_NONE) 	  	      /* Explode the 3D Woolz object */
  {
    if((errFlag = RecErrorFromWlz(
    		  WlzExplode3D(&objVecCount, &objVec, srcObj))) != REC_ERR_NONE)
    {
      errStr = errExplodeStr;
    }
    REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	    ("RecExplode3DObjToFile 02 %d\n",
	     (int )errFlag));
  }
  if(errFlag == REC_ERR_NONE)                       /* Create a section list */
  {
    if(((secList = RecSecNewSectionList(NULL)) == NULL) ||
       ((secList->list = HGUDlpListCreate(NULL)) == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
    else
    {
      secList->reconstruction.fileName = srcFName;
      secList->reconstruction.fileFormat = srcFFormat;
    }
    REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	    ("RecExplode3DObjToFile 03 %d\n",
	     (int )errFlag));
  }
  if(errFlag == REC_ERR_NONE)     /* Write 2D files, build/write section list */
  {
    fP = NULL;
    objIdx = 0;
    planeIdx = srcObj->domain.p->plane1;
    fileStrSpace = (int )strlen(dstDirStr) + (int )strlen(dstBodyStr) + 64;
    if(((fileStr = AlcMalloc(sizeof(char) * fileStrSpace)) == NULL) ||
       ((fileTemplateStr = AlcMalloc(sizeof(char) * fileStrSpace)) == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
    else
    {
      (void )sprintf(fileTemplateStr, "%s/%s", dstDirStr, dstBodyStr);
    }
    while((errFlag == REC_ERR_NONE) && (objIdx < objVecCount))
    {
      (void )sprintf(fileStr, "%s_%06d.wlz", fileTemplateStr, planeIdx);
      if((fP = fopen(fileStr, "w")) == NULL)
      {
	errStr = errWriteWoolzStr;
        errFlag = REC_ERR_WRITE;
      }
      if(errFlag == REC_ERR_NONE)
      {
	if((errFlag = RecErrorFromWlz(
		      WlzWriteObj(fP, *(objVec + objIdx)))) != REC_ERR_NONE)
	{
	  errStr = errWriteWoolzStr;
	}
        fclose(fP);
      }
      if(errFlag == REC_ERR_NONE)
      {
        errFlag = RecSecAppendListFromFiles(secList->list, NULL, &fileStr, 1,
					    planeIdx, 1);
        if(errFlag != REC_ERR_NONE)
	{
	  errStr = errListAppendStr;
	}
      }
      if(errFlag == REC_ERR_NONE)
      {
	++objIdx;
	++planeIdx;
      }
    }
    if(objVec && (objVecCount > 0))
    {
      for(objIdx = 0; objIdx < objVecCount; ++objIdx)
      {
	if(*(objVec + objIdx))
	{
	  (void )WlzFreeObj(*(objVec + objIdx));
	}
      }
      free(objVec);
    }
    if(errFlag == REC_ERR_NONE)
    {
      (void )sprintf(fileStr, "%s.bib", fileTemplateStr);
      if((fP = fopen(fileStr, "w")) == NULL)
      {
        errStr = errWriteListStr;
	errFlag = REC_ERR_WRITE;
      }
      else
      {
	errFlag = RecFileSecListWrite(fP, secList, objVecCount, eMsg);
	fclose(fP);
      }
    }
    if(fileStr)
    {
      free(fileStr);
    }
    if(fileTemplateStr)
    {
      free(fileTemplateStr);
    }
    if(secList)
    {
      HGUDlpListDestroy(secList->list);
      AlcFree(secList);
    }
    REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	    ("RecExplode3DObjToFile 04 %d\n",
	     (int )errFlag));
  }
  if((errFlag != REC_ERR_NONE) && eMsg && (*eMsg == NULL) && (errStr != NULL))
  {
    *eMsg = AlcStrDup(errStr);
  }
  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
  	  ("RecExplode3DObjToFile FX %d\n",
	   (int )errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecExplode3DFileToFile				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Given a destination directory and file body and a 3D
*		image file, the 3D object is and exploded into a
*		series of 2D section files (in the destination	
*		directory) each with a name formed from the given file
*		body and the source plane index. A section list is
*		also created in the directory using the given file body
*		together with the  '.bib' file extension.	
* Global refs:	-						
* Parameters:	char *dstDirStr:	Destination directory.	
*		char *dstBodyStr:	Destination file body.	
*		char *srcFile:		Source 3D image file.	
*		WlzEffFormat srcFmt:	Source 3D image file format.
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
RecError	RecExplode3DFileToFile(char *dstDirStr, char *dstBodyStr,
				      char *srcFile, WlzEffFormat srcFmt,
				      char **eMsg)
{
  RecError	errFlag = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  WlzObject	*srcObj = NULL;

  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecExplode3DFileToFile FE 0x%lx 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )dstDirStr, (unsigned long )dstBodyStr,
	   (unsigned long )srcFile, (int )srcFmt, (unsigned long )eMsg));
  if((dstDirStr == NULL) || (strlen(dstDirStr) == 0) ||
     (dstBodyStr == NULL) || (strlen(dstBodyStr) == 0) || 
     (srcFile == NULL) || (strlen(srcFile) == 0))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    srcObj = WlzEffReadObj(NULL, srcFile, srcFmt, 0, &wlzErr);
    errFlag = RecErrorFromWlz(wlzErr);
  }
  if(srcObj)
  {
    if(errFlag == REC_ERR_NONE)
    {
      errFlag = RecExplode3DObjToFile(dstDirStr, dstBodyStr, srcObj,
      				      srcFile, srcFmt, eMsg);
    }
    (void )WlzFreeObj(srcObj);
  }
  REC_DBG((REC_DBG_3D|REC_DBG_LVL_FN|REC_DBG_LVL_1),
  	  ("RecExplode3DFileToFile FX %d\n",
	   (int )errFlag));
  return(errFlag);
}
