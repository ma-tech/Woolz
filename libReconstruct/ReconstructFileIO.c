#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructFileIO.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides functions for file based I/O for the	
*		MRC Human Genetics Unit reconstruction library.	
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

static RecError RecFileSecHeadRead(int *version, FILE *fP, char **eMsg),
		RecFileHeadRead(int *version, FILE *fP, char **ident,
				int versionMin, int versionMax, char **eMsg),
		RecFileSecBodyRead(HGUDlpList *secList, int *numSec, FILE *fP,
			  int version, char **eMsg),
		RecFileIcsFileNames(char **fileBody,
				    char **icsFileName, char **idsFileName,
				    char *gvnFileName),
		RecFileSecRecordReadAndParse(RecSection *sec, int *doneFlag,
					     int *recNum, FILE *fP,
					     char **eMsg);
static BibFileField *RecFileTransfToField(WlzAffineTransform *transf);
static WlzAffineTransform *RecFileFieldToTransf(BibFileField *field,
						int defaultFlg);

typedef enum	       /* Warning: keep consistant with recFileTransfFieldS! */
{
  REC_FILE_TRANSF_MIN = 0,
  REC_FILE_TRANSF_TYPE = 0,
  REC_FILE_TRANSF_TX,
  REC_FILE_TRANSF_TY,
  REC_FILE_TRANSF_TZ,
  REC_FILE_TRANSF_SCALE,
  REC_FILE_TRANSF_THETA,
  REC_FILE_TRANSF_PHI,
  REC_FILE_TRANSF_ALPHA,
  REC_FILE_TRANSF_PSI,
  REC_FILE_TRANSF_XSI,
  REC_FILE_TRANSF_INVERT,
  REC_FILE_TRANSF_MAX
} RecFileTransfFields;

static char *recFileTransfFieldS[REC_FILE_TRANSF_MAX] =
{
  "TransformType",
  "TransformTx",
  "TransformTy",
  "TransformTz",
  "TransformScale",
  "TransformTheta",
  "TransformPhi",
  "TransformAlpha",
  "TransformPsi",
  "TransformXsi",
  "TransformInvert"
};

/************************************************************************
* Function:	RecFileSecWrite					
* Returns:	RecError:		Non zero if read fails.	
* Purpose:	Write registration sections file.		
* Notes:	Files are written using the restricted bixtex file
*		syntax of the BibFileIO library.		
*		Records have the format				
*		  @Ident					
*		  {						
*		    0,						
*		    Text = {registration section file},		
*		    Version = {1}				
*		  }						
*		  @Comment					
*		  {						
*		    0,						
*		    Text = {...},				
*		    .						
*		    .						
*		    Text = {...}				
*		  }						
*		  @Section					
*		  {						
*		    <section index>,				
*		    Iterations = {<iterations to find best match>},
*		    Correlation = {<correlation value>},	
*		    File = {<section image data file>},		
*		    TransformType = {1},			
*		    TransformTx = {<x>}				
*		    TransformTy = {<y>}				
*		    TransformTz = {<z>}				
*		    TransformScale = {<scale>}			
*		    TransformTheta = {<theta>}			
*		    TransformPhi = {<phi>}			
*		    TransformAlpha = {<alpha>}			
*		    TransformPsi = {<psi>}			
*		    TransformXsi = {<xsi>}			
*		    TransformInvert  = {<invert}		
*		  }						
* Global refs:	-						
* Parameters:	FILE *fP:		Output file stream.	
*		HGUDlpList *secList:	Sections list.		
*		int numSec:		Number of sections.	
*		char **eMsg:		Ptr for error message strings.
************************************************************************/
RecError	RecFileSecWrite(FILE *fP, HGUDlpList *secList, int numSec,
				char **eMsg)
{
  int		idx;
  RecError	errFlag = REC_ERR_NONE;
  time_t	tmpTime;
  char		*tmpS,
		*idxS = NULL,
		*itrS = NULL,
		*corS = NULL,
		*dateS = NULL,
		*hostS = NULL,
		*userS = NULL;
  RecSection 	*sec;
  HGUDlpListItem *secItem;
  BibFileRecord	*record = NULL;
  BibFileField	*field0 = NULL,
		*field1 = NULL,
		*field2 = NULL;
  char		tmpBuf[256];
  static char	unknownS[] = "unknown";

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecWrite FE 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )fP, (unsigned long )secList,
	   numSec, (unsigned long )eMsg));
  if((((field0 = BibFileFieldMakeVa("Text", "registration section file",
				    "Version", "1",
				    NULL)) == NULL) ||
      ((record = BibFileRecordMake("Ident", "0", field0)) == NULL)))
  {
    errFlag = REC_ERR_MALLOC;
  }
  if((errFlag == REC_ERR_NONE) && 
     (BibFileRecordWrite(fP, eMsg, record) != BIBFILE_ER_NONE))
  {
    errFlag = REC_ERR_WRITE;
  }
  if(record)
  {
    BibFileRecordFree(&record);
  }
  else
  {
    BibFileFieldFree(&field0);
  }
  if(errFlag == REC_ERR_NONE)
  {
    tmpS = getenv("USER");
    (void )sprintf(tmpBuf, "User: %s", tmpS?tmpS:unknownS);
    if((userS = AlcStrDup(tmpBuf)) == NULL)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    tmpTime = time(NULL);
    tmpS = ctime(&tmpTime);
    *(tmpS + strlen(tmpS) - 1) = '\0';
    (void )sprintf(tmpBuf, "Date: %s", tmpS?tmpS:unknownS);
    if((dateS = AlcStrDup(tmpBuf)) == NULL)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    tmpS = getenv("HOST");
    (void )sprintf(tmpBuf, "Host: %s", tmpS?tmpS:unknownS);
    if((hostS = AlcStrDup(tmpBuf)) == NULL)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if((errFlag == REC_ERR_NONE) &&
     (((field0 = BibFileFieldMakeVa("Text", userS,
				    "Text", dateS,
				    "Text", hostS,
				    NULL)) == NULL) ||
       ((record = BibFileRecordMake("Comment", "0", field0)) == NULL)))
  {
      errFlag = REC_ERR_MALLOC;
  }
  if((errFlag == REC_ERR_NONE) &&
     (BibFileRecordWrite(fP, eMsg, record) != BIBFILE_ER_NONE))
  {
    errFlag = REC_ERR_WRITE;
  }
  if(dateS)
  {
    AlcFree(dateS);
  }
  if(hostS)
  {
    AlcFree(hostS);
  }
  if(userS)
  {
    AlcFree(userS);
  }
  if(record)
  {
    BibFileRecordFree(&record);
  }
  else
    BibFileFieldFree(&field0);
  idx = 0;
  secItem = HGUDlpListHead(secList);
  while((errFlag == REC_ERR_NONE) && (idx < numSec) && secItem &&
        ((sec = (RecSection *)HGUDlpListEntryGet(secList, secItem)) != NULL))
  {
    if((sprintf(tmpBuf, "%d", sec->index) <= 0) ||
       ((idxS = AlcStrDup(tmpBuf)) == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
    if((errFlag == REC_ERR_NONE) &&
       ((sprintf(tmpBuf, "%d", sec->iterations) <= 0) ||
       ((itrS = AlcStrDup(tmpBuf)) == NULL)))
    {
      errFlag = REC_ERR_MALLOC;
    }
    if((errFlag == REC_ERR_NONE) &&
       ((sprintf(tmpBuf, "%1.7g", sec->correl) <= 0) ||
       ((corS = AlcStrDup(tmpBuf)) == NULL)))
    {
      errFlag = REC_ERR_MALLOC;
    }
    if((errFlag == REC_ERR_NONE) &&
       (((field0 = RecFileTransfToField(sec->transform)) == NULL) ||
        ((field1 = BibFileFieldMakeVa("Iterations", itrS,
				      "Correlation", corS,
			              "File", sec->imageFile,
			              NULL)) == NULL) ||
        ((field2 = BibFileFieldJoin(field1, field0, NULL)) == NULL) ||
        ((record = BibFileRecordMake("Section", idxS,
				     field2)) == NULL) ||
        (BibFileRecordWrite(fP, eMsg, record) != BIBFILE_ER_NONE)))
    {
      errFlag = REC_ERR_WRITE;
    }
    if(idxS)
    {
      AlcFree(idxS);
    }
    if(itrS)
    {
      AlcFree(itrS);
    }
    if(corS)
    {
      AlcFree(corS);
    }
    if(record)
    {
      BibFileRecordFree(&record);
    }
    else if(field2)
    {
      BibFileFieldFree(&field2);
    }
    else
    {
      BibFileFieldFree(&field0);
      BibFileFieldFree(&field1);
    }
    if(errFlag == REC_ERR_NONE)
    {
      secItem = HGUDlpListNext(secList, secItem);
    }
    ++idx;
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecObjRead				
* Returns:	RecError:		Non zero if read fails.	
* Purpose:	Read section image object from file.		
* Global refs:	-						
* Parameters:	RecSection *sec:	Section with obj to be read in.
*		char **eMsg:		Ptr for error message strings.
************************************************************************/
RecError	RecFileSecObjRead(RecSection *sec, char **eMsg)
{
  RecError	errFlag = REC_ERR_NONE;
  char		errBuf[256];
  FILE		*fP = NULL;
  WlzObject	*obj;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjRead FE 0x%lx 0x%lx\n",
	   (unsigned long )sec, (unsigned long )eMsg));
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_1),
  	  ("RecFileSecObjRead 01 0x%lx  0x%lx %d\n",
	   (unsigned long )sec,
	   (sec)? ((unsigned long )(sec->obj)): NULL,
	   (sec)? ((sec->obj)? (sec->obj->linkcount): -999): -999));
  if((sec->obj = WlzAssignObject(sec->obj, NULL)) == NULL)
  {
    if((sec->imageFile == NULL) || (strlen(sec->imageFile) == 0))
    {
      errFlag = REC_ERR_READ;
    }
    else if(RecSecIsEmpty(sec))
    {
      obj = WlzAssignObject(WlzMakeEmpty(&wlzErr), NULL);
      errFlag = RecErrorFromWlz(wlzErr);
    }
    else
    {
      if(((fP = fopen(sec->imageFile, "r")) == NULL) ||
	 ((obj = WlzAssignObject(WlzReadObj(fP, &wlzErr), NULL)) == NULL))
      {
	errFlag = (wlzErr == WLZ_ERR_NONE)? REC_ERR_READ:
				            RecErrorFromWlz(wlzErr);
	(void )sprintf(errBuf, "failed to open section image file %s.",
		       sec->imageFile);
      }
      if(fP)
      {
	(void )fclose(fP);
      }
    }
    if(errFlag == REC_ERR_NONE)
    {
      sec->obj = obj;
    }
    else if(*eMsg == NULL)
    {
      *eMsg = AlcStrDup(errBuf);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecObjFree				
* Returns:	void						
* Purpose:	Free (ie decrement linkcount) section image object.
* Global refs:	-						
* Parameters:	RecSection *sec:	Section with obj to be read in.
************************************************************************/
void		RecFileSecObjFree(RecSection *sec)
{
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjFree FE 0x%lx\n",
	   (unsigned long )sec));
  if(sec)
  {
    if(sec->obj)
    {
      REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_3),
      	      ("RecFileSecObjFree 01 0x%lx %d\n",
	       (unsigned long )(sec->obj), sec->obj->linkcount));
      WlzFreeObj(sec->obj);
      sec->obj = NULL;
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjFree FX\n"));
}

/************************************************************************
* Function:	RecFileSecObjsRead				
* Returns:	RecError:		Non zero if read fails.	
* Purpose:	Read section image objects from files for all sections
*		with indicies in given range.			
* Global refs:	-						
* Parameters:	HGUDlpList *secList:	Given section list.	
*		int minIdx:		Minimum index of range.	
*		int maxIdx:		Maximum index of range.	
*		int allSecFlag:		Flag to read all section's if
*					non zero.		
*		char **eMsg:		Ptr for error message strings.
************************************************************************/
RecError	RecFileSecObjsRead(HGUDlpList *secList, int minIdx, int maxIdx,
				   int allSecFlag, char **eMsg)
{
  RecError	errFlag = REC_ERR_NONE;
  HGUDlpListItem *secItem;
  RecSection	*sec;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjsRead FE 0x%lx %d %d %d 0x%lx\n",
	   (unsigned long )secList, minIdx, maxIdx, allSecFlag,
           (unsigned long )eMsg));

  if((secList == NULL) || ((allSecFlag == 0) && (minIdx > maxIdx)) ||
     ((secItem = HGUDlpListHead(secList)) == NULL) ||
     ((sec = (RecSection *)HGUDlpListEntryGet(secList, secItem)) == NULL))
  {
    errFlag = REC_ERR_READ;
  }
  else
  {
    while(sec && (errFlag == REC_ERR_NONE))
    {
      if(allSecFlag || ((sec->index >= minIdx) && (sec->index <= maxIdx)))
      {
	errFlag = RecFileSecObjRead(sec, eMsg);
      }
      secItem = HGUDlpListNext(secList, secItem);
      sec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjsRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecObjsFree				
* Returns:	void						
* Purpose:	Free section image objects from files for all sections
*		with indicies in given range.			
* Global refs:	-						
* Parameters:	HGUDlpList *secList:	Given section list.	
*		int minIdx:		Minimum index of range.	
*		int maxIdx:		Maximum index of range.	
*		int allSecFlag:		Flag to read all section's if
*					non zero.		
************************************************************************/
void		RecFileSecObjsFree(HGUDlpList *secList, int minIdx, int maxIdx,
				   int allSecFlag)
{
  HGUDlpListItem *secItem;
  RecSection	*sec;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjsFree FE 0x%lx %d %d %d\n",
	   (unsigned long )secList, minIdx, maxIdx, allSecFlag));

  if(secList && (allSecFlag || (minIdx <= maxIdx)) &&
     ((secItem = HGUDlpListHead(secList)) != NULL))
  {
    while((sec = (RecSection *)HGUDlpListEntryGet(secList, secItem)) != NULL)
    {
      if(allSecFlag || ((sec->index >= minIdx) && (sec->index <= maxIdx)))
      {
	RecFileSecObjFree(sec);
      }
      secItem = HGUDlpListNext(secList, secItem);
      sec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecObjsFree FX\n"));
}

/************************************************************************
* Function:	RecFileSecRead					
* Returns:	RecError:		Non zero on error.	
* Purpose:	Reads a bibtex style list of section transform records
*		from the given file stream. 			
* Global refs:	-						
* Parameters:	HGUDlpList *secList:	Section list.		
*		int *numSec:		Ptr for number of sections in 
*					the list.		
*		FILE *fP:		Input file.		
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
RecError	RecFileSecRead(HGUDlpList *secList, int *numSec, FILE *fP,
			       char **eMsg)
{
  int		version;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecRead FE 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )numSec,
	   (unsigned long )fP, (unsigned long)eMsg));
  if((secList == NULL) || (numSec == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  else
  {
    *numSec = 0;
    if(HGUDlpListCount(secList) > 0)   /* Free any existing sections in list */
    {
      HGUDlpListDeleteAll(secList);
    }
    errFlag = RecFileSecHeadRead(&version, fP, eMsg);
    if(errFlag == REC_ERR_NONE)
    {
      errFlag = RecFileSecBodyRead(secList, numSec, fP, version, eMsg);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecHeadRead				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Reads the header of a bibtex style list of section
*		transform records from the given file stream, to verify
*		that this is infact the expected file type and also to
*		establish the file version number.		
* Global refs:	-						
* Parameters:	int *version:		Ptr for file version number.
*		FILE *fP:		Input file.		
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
static RecError	RecFileSecHeadRead(int *version, FILE *fP, char **eMsg)
{
  int		versionMin = 1,
		versionMax = 1;
  RecError	errFlag = REC_ERR_NONE;
  static char   *fileIdent[] =
		{
		  "registration",
		  "section",
		  "file",
		  ""
		};

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecHeadRead FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )version, (unsigned long )fP, (unsigned long)eMsg));
  errFlag = RecFileHeadRead(version, fP, fileIdent, versionMin, versionMax,
			    eMsg);
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecHeadRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileHeadRead					
* Returns:	RecError:		Non zero if read fails.	
* Purpose:	Read the file header to establish that this is a valid 
*		file and its version number.			
* Global refs:	-						
* Parameters:	int *version:		Ptr for version number.	
*		FILE *fP:		Input file stream.	
*		char **ident:		File identification strings.
*		char **eMsg:		Ptr for error message strings.
*		int *version:		Ptr for version number.	
*		FILE *fP:		Opened file pointer.	
************************************************************************/
static RecError	RecFileHeadRead(int *version, FILE *fP, char **ident,
				int versionMin, int versionMax, char **eMsg)
{
  BibFileError	bibErrFlag = BIBFILE_ER_NONE;
  RecError	errFlag = REC_ERR_NONE;
  BibFileField	*idField,
		*vField;
  BibFileRecord	*idRecord = NULL,
		*lastRecord = NULL;
  char		**idWord;
  char		*bibErrMsg = NULL,
		*idTok,
		*tmpPtr;
  char		errBuf[256];

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileHeadRead FE 0x%lx 0x%lx 0x%lx %d %d 0x%lx\n",
	   (unsigned long )version, (unsigned long )fP, (unsigned long )ident,
	   versionMin, versionMax, (unsigned long )eMsg));
  while(((bibErrFlag = BibFileRecordRead(&idRecord, &bibErrMsg,
					 fP)) == BIBFILE_ER_NONE) &&
	idRecord && idRecord->name && idRecord->id &&
	strcmp(idRecord->name, "Ident"))
  {
    if(lastRecord)
    {
      BibFileRecordFree(&lastRecord);
    }
    lastRecord = idRecord;
  }
  if(bibErrMsg)
  {
    AlcFree(bibErrMsg);
  }
  if((bibErrFlag != BIBFILE_ER_NONE) || 
     (idRecord == NULL) ||
     (idRecord->name == NULL) || (idRecord->id == NULL) ||
     ((idField = idRecord->field) == NULL) ||
     (idField->name == NULL) || strcmp(idField->name, "Text") ||
     (idField->value  == NULL) ||
     ((vField = idField->next) == NULL) ||
     (vField->name == NULL) || strcmp(vField->name, "Version") ||
     (idField->value  == NULL))
  {
    errFlag = REC_ERR_READ;
  }
  if(errFlag == REC_ERR_NONE)
  {
    idWord = ident;
    tmpPtr = idField->value;
    while((errFlag == REC_ERR_NONE) &&
	  idWord && *idWord && **idWord &&
	  ((idTok = strtok(tmpPtr, " \t")) != NULL) && *idTok)
    {
      if(strcmp(*idWord, idTok))
      {
	errFlag = REC_ERR_READ;
      }
      else
      {
        ++idWord;
	tmpPtr = NULL;
      }
    }
  }
  if(errFlag != REC_ERR_NONE)
  {
    (void )sprintf(errBuf,
		   "inappropriate file identification record.");
  }
  if(errFlag == REC_ERR_NONE)
  {
    idTok = strtok(vField->value, "\t\n");
    if((version == NULL) || (idTok == NULL) || (*idTok == '\0') ||
       (sscanf(idTok, "%d", version) != 1) ||
       (*version < versionMin) || (*version > versionMax))
    {
      errFlag = REC_ERR_READ;
      (void )sprintf(errBuf,
		     "unsupported file version %d, versions %d to %d only.",
		     *version, versionMin, versionMax);
    }
  }
  if(idRecord)
  {
    BibFileRecordFree(&idRecord);
  }
  if((errFlag != REC_ERR_NONE) && (*eMsg == NULL))
  {
    *eMsg = AlcStrDup(errBuf);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileHeadRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecBodyRead				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Reads the body of a bibtex style list of section
*		transform records from the given file stream, parsing
*		parsing them and building a linked list of sections.
* Global refs:	-						
* Parameters:	HGUDlpList *secList:	Section list.		
*		int *numSec:		Ptr for number of sections in 
*					the list.		
*		FILE *fP:		Input file.		
*		int version:		File version number.	
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
static RecError	RecFileSecBodyRead(HGUDlpList *secList, int *numSec, FILE *fP,
			  int version, char **eMsg)
{
  int		doneFlag = 0,
		secCount = 0,
		recCount = 0;
  RecSection	*sec = NULL;
  RecError	errFlag = REC_ERR_NONE;
  char		errBuf[256];
  
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecBodyRead FE 0x%lx 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )secList, (unsigned long )numSec,
	   (unsigned long )fP, version, (unsigned long )eMsg));
  if((secList == NULL) || (numSec == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(version != 1)
    {
      errFlag = REC_ERR_READ;
      (void )sprintf(errBuf,
		     "version %d section file, only version %d supported.",
		     version, 1);
    }
  }
  while((errFlag == REC_ERR_NONE) && (doneFlag == 0))
  {
    if((sec = (RecSection *)AlcCalloc(1, sizeof(RecSection))) == NULL)
    {
      errFlag = REC_ERR_MALLOC;
    }
    if(errFlag == REC_ERR_NONE)
    {
      errFlag = RecFileSecRecordReadAndParse(sec, &doneFlag, &recCount,
      					     fP, eMsg);
    }
    if((errFlag == REC_ERR_NONE) && (doneFlag == 0))
    {
      ++(sec->linkcount);
      if(HGUDlpListAppend(secList, NULL, (void *)sec,
      			  (void (*)(void *) )RecSecFree) == NULL)
      {
        errFlag = REC_ERR_MALLOC; 			 /* May not be true! */
      }
      else
      {
	sec = NULL;
        ++secCount;
      }
    }
    REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
	    ("RecFileSecBodyRead 01 0x%lx %d %d %d\n",
	     (unsigned long )sec, secCount, doneFlag, errFlag));
  }
  if(((errFlag != REC_ERR_NONE) || doneFlag) && sec)
  {
    RecSecFree(sec);
  }
  if((errFlag == REC_ERR_NONE) && (secCount < 2))
  {
    errFlag = REC_ERR_READ;
    (void )strcpy(errBuf,
		  "at least two sections are required in a section list.");
  }
  if(errFlag == REC_ERR_NONE)
  {
    *numSec = secCount;
  }
  else
  {
    if(*eMsg == NULL)
    {
      *eMsg = AlcStrDup(errBuf);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecBodyRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileObjWlzRead				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Reads a woolz object from the the given file stream.
* Global refs:	-						
* Parameters:	FILE *fP:		Output file stream.	
*		WlzObject **obj:	Given woolz object.	
************************************************************************/
RecError	RecFileObjWlzRead(FILE *fP, WlzObject **obj)
{
  RecError	errFlag = REC_ERR_FUNC;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileObjWlzRead FE 0x%lx 0x%lx\n",
	   (unsigned long )fP, (unsigned long )obj));
  if(fP && obj)
  {
    if((*obj = WlzReadObj(fP, &wlzErr)) != NULL)
    {
      errFlag = (wlzErr == WLZ_ERR_NONE)? REC_ERR_NONE:
      					  RecErrorFromWlz(wlzErr);
    }
    else
    {
      errFlag = REC_ERR_READ;
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileObjWlzRead FX %d\n",
	   errFlag));
  return(errFlag);
}


/************************************************************************
* Function:	RecFileObjWlzWrite				
* Returns:	RecError:		Non zero on error.	
* Purpose:	Writes the given woolz object to the given file stream.
* Global refs:	-						
* Parameters:	FILE *fP:		Output file stream.	
*		WlzObject *obj:		Given woolz object.	
************************************************************************/
RecError	RecFileObjWlzWrite(FILE *fP, WlzObject *obj)
{
  RecError	errFlag = REC_ERR_WRITE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileObjWlzWrite FE 0x%lx 0x%lx\n",
	   (unsigned long )fP, (unsigned long )obj));
  if(fP && obj)
  {
    if((wlzErr = WlzWriteObj(fP, obj)) != WLZ_ERR_NONE)
    {
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileObjWlzWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileSecRecordReadAndParse			
* Returns:	RecError		Non zero on error.	
* Purpose:	Reads and parses the next BibFile section record from
*		the given file stream. The record is parsed into the
*		given section data structure.			
* Global refs:	-						
* Parameters:	RecSection *sec:	Given section (allocated).
*		int *doneFlag:		Non-error, done flag.	
*		int *recNum:		Used to count records parsed.
*		FILE *fP:		Input file stream.	
*		char **eMsg:		Ptr for any error messages.
************************************************************************/
static RecError	RecFileSecRecordReadAndParse(RecSection *sec, int *doneFlag,
					     int *recNum, FILE *fP,
					     char **eMsg)
{
  RecError	errFlag = REC_ERR_NONE;
  BibFileError	bibErrFlag = BIBFILE_ER_NONE;
  int		secIdx,
  		secIter;
  char		*bibErrMsg = NULL,
  		*secFile = NULL;
  double	secCor;
  WlzAffineTransform	*secTransf = NULL;
  BibFileRecord	*secRecord = NULL,
  		*lastRecord = NULL;
  char		errBuf[256];
  static char	intFmt[] = "%d",
  		dblFmt[] = "%lf",
		strFmt[] = "%s";

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
  	  ("RecFileSecRecordReadAndParse FE 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )sec, (unsigned long )doneFlag, (unsigned long )fP,
	   (unsigned long )eMsg));
  if((sec == NULL) || (doneFlag == NULL) || (recNum == NULL) || (fP == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    while(((bibErrFlag = BibFileRecordRead(&secRecord, &bibErrMsg,
					   fP)) == BIBFILE_ER_NONE) &&
	  secRecord && secRecord->name && strcmp(secRecord->name, "Section"))
    {
      REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
	      ("RecFileSecRecordReadAndParse 01 0x%lx >%s< %d\n",
	       (unsigned long )secRecord,
	       (secRecord->name)? (secRecord->name): ("null"),
	       *recNum));
      ++*recNum;
      if(lastRecord)
      {
	BibFileRecordFree(&lastRecord);
      }
      lastRecord = secRecord;
      secRecord = NULL;
    }
    ++*recNum;
    switch(bibErrFlag)
    {
      case BIBFILE_ER_NONE:
	break;
      case BIBFILE_ER_EOF:
	*doneFlag = 1;
	break;
      default:
	errFlag = REC_ERR_READ;
	break;
    }
    REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
	    ("RecFileSecRecordReadAndParse 02 0x%lx >%s< %d %d %d %d\n",
	     (unsigned long )secRecord,
      (secRecord)? ((secRecord->name)? (secRecord->name): ("null")): ("NULL"),
	     *recNum, *doneFlag, bibErrFlag, errFlag));
  }
  if(*doneFlag == 0)
  {
    if((errFlag != REC_ERR_NONE) || (secRecord->id == NULL) ||
       (sscanf(secRecord->id, "%d", &secIdx) != 1))
    {
      errFlag = REC_ERR_READ;
      (void )sprintf(errBuf, "syntax error in section file, record %d%s%s%s.",
		     *recNum,
		     bibErrMsg?" (":"",
		     bibErrMsg?bibErrMsg:"",
		     bibErrMsg?")":"");
      if(bibErrMsg)
      {
        AlcFree(bibErrMsg);
      }
    }
    if(errFlag == REC_ERR_NONE)
    {
      secIter = 0;
      secCor = 0.0;
      (void )BibFileFieldParseFmt(secRecord->field,
                                  &secIter, intFmt, "Iterations",
                                  &secCor, dblFmt, "Correlation",
                                  NULL);
      if((BibFileFieldParseFmt(secRecord->field,
                               &secFile, strFmt, "File",
                               NULL) != 1) ||
         (secFile == NULL) || (*secFile == '\0') ||
         ((secTransf = RecFileFieldToTransf(secRecord->field, 1)) == NULL))
      {
        errFlag = REC_ERR_READ;
        (void )sprintf(errBuf, "syntax error in section file, record %d.",
                       *recNum);
      }
      else
      {
        sec->index = secIdx;
        sec->iterations = secIter;
        sec->correl = secCor;
        sec->imageFile = secFile;
	sec->transform = WlzAssignAffineTransform(secTransf, NULL);
        sec->obj = NULL;
      }
    }
    if(secRecord)
    {
      BibFileRecordFree(&secRecord);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
  	  ("RecFileSecRecordReadAndParse FX %d\n",
	   errFlag));
  return(errFlag);
}

/************************************************************************
* Function:	RecFileTransfToField				
* Returns:      BibFileField *:		New field, NULL on error.
* Purpose:      Creates fields using the bibtex file syntax for the
*		given transform and allocate storage as required.
* Global refs:  char *bibFileTransfFieldS: Transform field names.
* Parameters:   WlzAffineTransform *transf:	Given transform.	
************************************************************************/
BibFileField	*RecFileTransfToField(WlzAffineTransform *transf)
{
  BibFileField	*field = NULL;
  char 		tTypeS[32],
		tTxS[32],
		tTyS[32],
		tTzS[32],
		tScaleS[32],
		tThetaS[32],
		tPhiS[32],
		tAlphaS[32],
		tPsiS[32],
		tXsiS[32],
		tInvertS[32];

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileTransfToField FE 0x%lx\n",
	   (unsigned long )transf));
  if(transf)
  {
    (void )sprintf(tTypeS, "%d", transf->type);
    (void )sprintf(tTxS, "%g", transf->tx);
    (void )sprintf(tTyS, "%g", transf->ty);
    (void )sprintf(tTzS, "%g", transf->tz);
    (void )sprintf(tScaleS, "%g", transf->scale);
    (void )sprintf(tThetaS, "%g", transf->theta);
    (void )sprintf(tPhiS, "%g", transf->phi);
    (void )sprintf(tAlphaS, "%g", transf->alpha);
    (void )sprintf(tPsiS, "%g", transf->psi);
    (void )sprintf(tXsiS, "%g", transf->xsi);
    (void )sprintf(tInvertS, "%d", transf->invert);
    field = BibFileFieldMakeVa(
    			 recFileTransfFieldS[REC_FILE_TRANSF_TYPE], tTypeS,
			 recFileTransfFieldS[REC_FILE_TRANSF_TX], tTxS,
			 recFileTransfFieldS[REC_FILE_TRANSF_TY], tTyS,
			 recFileTransfFieldS[REC_FILE_TRANSF_TZ], tTzS,
			 recFileTransfFieldS[REC_FILE_TRANSF_SCALE], tScaleS,
			 recFileTransfFieldS[REC_FILE_TRANSF_THETA], tThetaS,
			 recFileTransfFieldS[REC_FILE_TRANSF_PHI], tPhiS,
			 recFileTransfFieldS[REC_FILE_TRANSF_ALPHA], tAlphaS,
			 recFileTransfFieldS[REC_FILE_TRANSF_PSI], tPsiS,
			 recFileTransfFieldS[REC_FILE_TRANSF_XSI], tXsiS,
			 recFileTransfFieldS[REC_FILE_TRANSF_INVERT], tInvertS,
			 NULL);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileTransfToField FX 0x%lx\n",
	   (unsigned long )field));
  return(field);
}

/************************************************************************
* Function:	RecFileFieldToTransf				
* Returns:      WlzAffineTransform:	New transform, NULL on error.
* Purpose:      Parses bibtex file syntax fields for the transform
*		fields and allocate storage for the new transform as
*		required.					
* Note:		The transform's linkcount is not set by this function.
* Global refs:  char *bibFileTransfFieldS: Transform field names.
* Parameters:   BibFileField *field:	Given field which is at the top
*					of the transform fields to be
*					parsed.			
*		int defaultFlg:		Allow default transform fields 
*					if non zero.		
************************************************************************/
WlzAffineTransform	*RecFileFieldToTransf(BibFileField *field,
					      int defaultFlg)
{
  int		count,
		tType,
		tInvert;
  double	tX,
		tY,
		tZ,
		tScale,
		tTheta,
		tPhi,
		tAlpha,
		tPsi,
		tXsi;
  WlzAffineTransform	*transf = NULL;
  static char	intFmt[] = "%d",
		dblFmt[] = "%lf";

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileFieldToTransf FE 0x%lx %d\n",
	   (unsigned long )field, defaultFlg));
  if(defaultFlg)
  {
    tType = WLZ_TRANSFORM_2D_AFFINE;
    tX = 0.0;
    tY = 0.0;
    tZ = 0.0;
    tScale = 1.0;
    tTheta = 0.0;
    tPhi = 0.0;
    tAlpha = 0.0;
    tPsi = 0.0;
    tXsi = 0.0;
    tInvert = 0;
  }
  count = BibFileFieldParseFmt(field,
	     &tType,	intFmt,	recFileTransfFieldS[REC_FILE_TRANSF_TYPE],
	     &tX,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_TX],
	     &tY,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_TY],
	     &tZ,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_TZ],
	     &tScale,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_SCALE],
	     &tTheta,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_THETA],
	     &tPhi,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_PHI],
	     &tAlpha,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_ALPHA],
	     &tPsi,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_PSI],
	     &tXsi,	dblFmt,	recFileTransfFieldS[REC_FILE_TRANSF_XSI],
	     &tInvert,	intFmt,	recFileTransfFieldS[REC_FILE_TRANSF_INVERT],
	     NULL);
   REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
   	   ("RecFileFieldToTransf 01 %d %d %f %f %f %f %f %f %f %f %f %d\n",
	    count, tType, tX, tY, tZ, tScale, tTheta, tPhi, tAlpha, tPsi, tXsi,
	    tInvert));
   if((defaultFlg || (count == REC_FILE_TRANSF_MAX)) &&
      (tType == WLZ_TRANSFORM_2D_AFFINE))
   {
     transf = WlzAffineTransformFromPrim(tType, tX, tY, 0.0,
     					 tScale, tTheta, tPhi,
					 tAlpha, tPsi, tXsi,
					 tInvert,
					 NULL);
   }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileFieldToTransf FX 0x%lx\n",
	   (unsigned long )transf));
  return(transf);
}
