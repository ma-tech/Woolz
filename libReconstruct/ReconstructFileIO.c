#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructFileIO_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructFileIO.c
* \author       Bill Hill
* \date         November 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
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
* \brief	File based I/O functions for the Reconstruct library.
* \ingroup	Reconstruct
*/
#include <Reconstruct.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

static RecError			RecFileSecSecListRead(
				  HGUDlpList *secList,
				  int *numSec,
				  int version,
				  FILE *fP,
				  char **eMsg) ;
static RecError			RecFileSecAtrbRead(
				  RecSectionListAtrb *atrb,
				  int version,
				  FILE *fP,
				  char **eMsg);
static RecError  		RecFileSecRecRead(
				  RecReconstruction *rec,
				  int version,
				  FILE *fP,
				  char **eMsg);
static RecError 		RecFileSecIdentRead(
				  int *version,
				  FILE *fP,
				  char **eMsg);
static RecError 		RecFileSecIdentWrite(
				  FILE *fP,
				  char **eMsg);
static RecError 		RecFileSecCommentsWrite(
				  FILE *fP,
				  char **eMsg);
static RecError			RecFileSecAtrbWrite(
				  FILE *fP,
				  RecSectionListAtrb *atrb,
				  char **eMsg);
static RecError 		RecFileSecRecWrite(
				  FILE *fP,
				  RecReconstruction *rec,
				  char **eMsg);
static RecError 		RecFileSecSecListWrite(
				  FILE *fP,
				  HGUDlpList *secList,
				  int numSec,
				  char **eMsg);
static RecError			RecFileSecRecordReadAndParse(
				  RecSection *sec,
				  int *doneFlag,
				  int *recNum,
				  FILE *fP,
				  char **eMsg);
static BibFileField 		*RecFileTransfToField(
				  WlzAffineTransform *transf);
static WlzAffineTransform 	*RecFileFieldToTransf(
				  BibFileField *field,
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

const char *recFileTransfFieldS[REC_FILE_TRANSF_MAX] =
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

typedef enum             /* Warning: keep consistant with recFileAtrbValueS! */
{
  REC_FILE_SECATRB_MIN = 0,
  REC_FILE_SECATRB_REL = 0,
  REC_FILE_SECATRB_ABS,
  REC_FILE_SECATRB_MAX
} RecFileSecAtrbValues;

const char *recFileAtrbValueS[REC_FILE_SECATRB_MAX] =
{
  "rel",
  "abs"
};

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief 	Writes registration sections file.

* 		Files are written using the restricted bixtex file
*		syntax of the BibFileIO library.
 		File format version 1 records have the format:
* \verbatim
 		  @Ident
 		  {
 		    0,
 		    Text = {registration section file},
 		    Version = {1}
 		  }
 		  @Comment
 		  {
 		    0,
 		    Text = {...},
 		    .
 		    .
 		    Text = {...}
 		  }
 		  @Section
 		  {
 		    <section index>,
 		    Iterations = {<iterations to find best match>},
 		    Correlation = {<correlation value>},
 		    File = {<section image data file>},
 		    TransformType = {1},
 		    TransformTx = {<x>}
 		    TransformTy = {<y>}
 		    TransformTz = {<z>}
 		    TransformScale = {<scale>}
 		    TransformTheta = {<theta>}
 		    TransformPhi = {<phi>}
 		    TransformAlpha = {<alpha>}
 		    TransformPsi = {<psi>}
 		    TransformXsi = {<xsi>}
 		    TransformInvert  = {<invert}
 		  }
 		File format version 2 records have the format:
 		  @Ident
 		  {
 		    0,
 		    Text = {registration section file},
 		    Version = {2}
 		  }
 		  @Comment
 		  {
 		    0,
 		    Text = {...},
 		    .
 		    .
 		    Text = {...}
 		  }
 		  @Attributes
 		  {
 		    0,
 		    TransformMode = {"rel"|"abs"}
 		  }
 		  @Reconstruction
 		  {
 		    0,
 		    File = {<complete reconstruction file name string>}
 		    Format = {<file format extension string>}
 		    Filter = {0|1}
 		    FastSampling = {0|1}
 		    IntScale = {0|1}
 		    Greedy = {0|1}
 		    Scale = {<scale x> <scale y> <scale z>}
 		    MatchHistograms = {0|1}
 		    MatchSection = {<section index>},
 		    ClipSrc = {0|1}
 		    SrcBox = {<x min> <y min> <z min> <x max> <y max> <z max>}
 		    ClipDst = {0|1}
 		    DstBox = {<x min> <y min> <z min> <x max> <y max> <z max>}
 		  }
 		  @Section
 		  {
 		    <section index>,
 		    Iterations = {<iterations to find best match>},
 		    Correlation = {<correlation value>},
 		    File = {<section image data file>},
 		    TransformType = {1},
 		    TransformTx = {<x>}
 		    TransformTy = {<y>}
 		    TransformTz = {<z>}
 		    TransformScale = {<scale>}
 		    TransformTheta = {<theta>}
 		    TransformPhi = {<phi>}
 		    TransformAlpha = {<alpha>}
 		    TransformPsi = {<psi>}
 		    TransformXsi = {<xsi>}
 		    TransformInvert  = {<invert}
 		  }
* \endverbatim
* \param	fP			Output file pointer.
* \param	secList			Sections list.
* \param	numSec			Number of sections.
* \param	eMsg			Destination pointer for messages.
*/
RecError	RecFileSecListWrite(FILE *fP, RecSectionList *secList,
				    int numSec, char **eMsg)
{
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecListWrite FE 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )fP, (unsigned long )secList,
	   numSec, (unsigned long )eMsg));
  if((fP == NULL) || (secList == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecIdentWrite(fP, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecCommentsWrite(fP, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecAtrbWrite(fP, &(secList->attributes), eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecRecWrite(fP, &(secList->reconstruction), eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecSecListWrite(fP, secList->list, numSec, eMsg);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecListWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes out file identification.
* \param	fP			Output file.
* \param	eMsg			Destination pointer for error messages.
*/
static RecError	RecFileSecIdentWrite(FILE *fP, char **eMsg)
{
  BibFileRecord	*record = NULL;
  BibFileField	*field = NULL;
  const char	versionS[] = "2";
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
          ("RecFileSecIdentWrite FE 0x%lx 0x%lx\n",
	  (unsigned long )fP, (unsigned long )eMsg));
  if((((field = BibFileFieldMakeVa("Text", "registration section file",
				   "Version", versionS,
				   NULL)) == NULL) ||
      ((record = BibFileRecordMake("Ident", "0", field)) == NULL)))
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
    BibFileFieldFree(&field);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecIdentWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code
* \ingroup	Reconstruct
* \brief	Writes the comments section of an output file.
* \param	fP			Output file.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecCommentsWrite(FILE *fP, char **eMsg)
{
  time_t	tmpTime;
  char		*tmpS,
		*dateS = NULL,
		*hostS = NULL,
		*userS = NULL;
  BibFileRecord	*record = NULL;
  BibFileField	*field = NULL;
  char		tmpBuf[256];
  const char	unknownS[] = "unknown";
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
  	  ("RecFileSecListCommentsWrite FE 0x%lx 0x%lx\n",
	   (unsigned long )fP, (unsigned long )eMsg));
  if(errFlag == REC_ERR_NONE)
  {
    tmpS = getenv("USER");
    (void )sprintf(tmpBuf, "User: %s", tmpS? tmpS: unknownS);
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
     (((field = BibFileFieldMakeVa("Text", userS,
				   "Text", dateS,
				   "Text", hostS,
				   NULL)) == NULL) ||
       ((record = BibFileRecordMake("Comment", "0", field)) == NULL)))
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
  {
    BibFileFieldFree(&field);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
  	  ("RecFileSecListCommentsWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes the section list attributes. Global recFileAtrbValueS
*		is used for the attribute strings.
* \param	fP			Output file.
* \param	atrb			Section list attributes.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecAtrbWrite(FILE *fP, RecSectionListAtrb *atrb,
				   char **eMsg)
{
  BibFileRecord	*record = NULL;
  BibFileField	*field = NULL;
  const char	*trModeS;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
  	  ("RecFileSecAtrbWrite FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )fP, (unsigned long )atrb, (unsigned long )eMsg));
  switch(atrb->trMode)
  {
    case REC_TRMODE_REL:
      trModeS = recFileAtrbValueS[REC_FILE_SECATRB_REL];
      break;
    case REC_TRMODE_ABS:
      trModeS = recFileAtrbValueS[REC_FILE_SECATRB_ABS];
      break;
    default:
      errFlag = REC_ERR_FUNC;
      break;

  }
  if(errFlag == REC_ERR_NONE)
  {
    if((((field = BibFileFieldMake("TransformMode", (char *)trModeS,
    				   NULL)) == NULL) ||
       ((record = BibFileRecordMake("Attributes", "0", field)) == NULL)))
    {
      errFlag = REC_ERR_MALLOC;
    }
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
    BibFileFieldFree(&field);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
          ("RecFileSecAtrbWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes out reconstruction details.
* \param	fP			Output file.
* \param	rec			Reconstruction details.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecRecWrite(FILE *fP, RecReconstruction *rec,
				   char **eMsg)
{
  BibFileRecord	*record = NULL;
  BibFileField	*topField = NULL;
  const char	*fileNameS,
  		*fileFormatS;
  char		scaleS[256],
  		matchSecS[256],
  		srcBoxS[256],
		dstBoxS[256];
  RecError	errFlag = REC_ERR_NONE;
  const char	zeroS[] = "0",
  		oneS[] = "1",
  		nullS[] = "";

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
  	  ("RecFileSecRecWrite FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )fP, (unsigned long )rec, (unsigned long )eMsg));
  fileNameS = (rec->fileName)? rec->fileName: nullS;
  (void )WlzEffStringFromFormat(rec->fileFormat, &fileFormatS);
  if(fileFormatS == NULL)
  {
    fileFormatS = nullS;
  }
  (void )sprintf(scaleS, "%g %g %g",
		 rec->scale.vtX, rec->scale.vtY, rec->scale.vtZ);
  if(rec->matchHst)
  {
    (void )sprintf(matchSecS, "%d", rec->matchSec);
  }
  else
  {
    strcpy(matchSecS, zeroS);
  }
  if(rec->clipSrc)
  {
    (void )sprintf(srcBoxS, "%g %g %g %g %g %g",
		   rec->srcBox.xMin, rec->srcBox.yMin, rec->srcBox.zMin,
		   rec->srcBox.xMax, rec->srcBox.yMax, rec->srcBox.zMax);
  }
  else
  {
    (void )sprintf(srcBoxS, "0.0 0.0 0.0 0.0 0.0 0.0");
  }
  if(rec->clipDst)
  {
    (void )sprintf(dstBoxS, "%g %g %g %g %g %g",
		   rec->dstBox.xMin, rec->dstBox.yMin, rec->dstBox.zMin,
		   rec->dstBox.xMax, rec->dstBox.yMax, rec->dstBox.zMax);
  }
  else
  {
    (void )sprintf(dstBoxS, "0.0 0.0 0.0 0.0 0.0 0.0");
  }
  if(((topField = BibFileFieldMakeVa(
			   "File", (char *)fileNameS,
			   "Format", (char *)fileFormatS,
			   "Filter", (rec->gaussFlt)? oneS: zeroS,
			   "FastSampling", (rec->fastSam)? oneS: zeroS,
			   "IntScale", (rec->intScale)? oneS: zeroS,
			   "Greedy", (rec->greedy)? oneS: zeroS,
			   "Scale", scaleS,
			   "MatchHistograms", (rec->matchHst)? oneS: zeroS,
			   "MatchSection", matchSecS,
			   "ClipSrc", (rec->clipSrc)? oneS: zeroS,
			   "SrcBox", srcBoxS,
			   "ClipDst", (rec->clipDst)? oneS: zeroS,
			   "DstBox", dstBoxS,
			   NULL)) == NULL) ||
     ((record = BibFileRecordMake("Reconstruction", "0",
				  topField)) == NULL))
  {
    errFlag = REC_ERR_MALLOC;
  }
  else if(BibFileRecordWrite(fP, eMsg, record) != BIBFILE_ER_NONE)
  {
    errFlag = REC_ERR_WRITE;
  }
  if(record)
  {
    BibFileRecordFree(&record);
  }
  else if(topField)
  {
    BibFileFieldFree(&topField);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
          ("RecFileSecRecWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes the body of the section list to the given file.
* \param	fP			Output file.
* \param	secList			Section list to write.
* \param	numSec			Number of sections to write.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecSecListWrite(FILE *fP, HGUDlpList *secList,
				       int numSec, char **eMsg)
{
  int		idx;
  char		*idxS = NULL,
		*itrS = NULL,
		*corS = NULL;
  RecSection 	*sec;
  HGUDlpListItem *secItem;
  BibFileRecord	*record = NULL;
  BibFileField	*field0 = NULL,
		*field1 = NULL,
		*field2 = NULL;
  char		tmpBuf[256];
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecListBodyWrite FE 0x%lx 0x%lx %d 0x%lx\n",
	   (unsigned long )fP, (unsigned long )secList,
	   numSec, (unsigned long )eMsg));
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
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecListBodyWrite FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads section image object from file.
* \param	sec			Section with obj to be read.
* \param	eMsg			Destination pointer for messages.
*/
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
	   (sec)? ((unsigned long )(sec->obj)): (unsigned long)0,
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

/*!
* \ingroup	Reconstruct
* \brief	Frees section image object.
* \param	sec			Section with obj.
*/
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

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads section image objects from files for all sections
*		with indicies in given range.
* \param	secList			Given section list.
* \param	minIdx			Minimum index of range.
* \param	maxIdx			Maximum index of range.
* \param	allSecFlag		Flag to read all section's if
*					non zero.
* \param	eMsg			Destination pointer for messages.
*/
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

/*!
* \ingroup	Reconstruct
* \brief	Free section image objects from files for all sections
*		with indicies in given range.
* \param	secList			Given section list.
* \param	minIdx			Minimum index of range.
* \param	maxIdx			Maximum index of range.
* \param	allSecFlag		Flag to read all section's if non zero.
*/
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

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads a bibtex style reconstruction list which contains
*		both details of how a reconstruction was made and the
*		list of section transform records.
* \param	secList			Reconstruction section list.
* \param	numSec			Destination pointer for the number
*					of sections in the list.
* \param	fP			Destination error pointer.
* \param	eMsg			Destination pointer for messages.
*/
RecError	RecFileSecListRead(RecSectionList *secList,
				   int *numSec, FILE *fP, char **eMsg)
{
  int		version;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecListRead FE 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )numSec,
	   (unsigned long )fP, (unsigned long)eMsg));
  if((secList == NULL) || (numSec == NULL) || (fP == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecIdentRead(&version, fP, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    *numSec = 0;
    secList->attributes.currentItem = NULL;
    if((secList->list = HGUDlpListCreate(NULL)) == NULL)
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecAtrbRead(&(secList->attributes), version,
    				 fP, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecRecRead(&(secList->reconstruction), version,
    				fP, eMsg);
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecSecListRead(secList->list, numSec,
    				    version, fP, eMsg);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads the file identification to establish that this is
*		a valid  reconstruction file and that is's version
* 		number is reasonable.
* \param	version			Destination pinter for version number.
* \param	fP			File pointer.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecIdentRead(int *version, FILE *fP, char **eMsg)
/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief
* \param	(REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2)
*/
{
  BibFileField	*idField,
		*vField;
  BibFileRecord	*idRecord = NULL,
		*lastRecord = NULL;
  const char	**idWord;
  char		*bibErrMsg = NULL,
		*idTok,
		*tmpPtr;
  char		errBuf[256];
  BibFileError	bibErrFlag = BIBFILE_ER_NONE;
  RecError	errFlag = REC_ERR_NONE;
  const int	versionMin = 1,
		versionMax = 2;
  const char   *fileIdent[] =
		{
		  "registration",
		  "section",
		  "file",
		  ""
		};

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecIdentRead FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )version, (unsigned long )fP,
	   (unsigned long )eMsg));
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
    idWord = fileIdent;
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
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecIdentRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads the reconstruct section list attributes.
* \param	atrb			Section list attributres data 
* 					attributres to be filled in.
* \param	version			Version number.
* \param	fP			File pointer.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	 RecFileSecAtrbRead(RecSectionListAtrb *atrb, int version,
				    FILE *fP, char **eMsg)
{

  BibFileRecord	*atRecord = NULL,
		*lastRecord = NULL;
  char		*trModeS = NULL,
  		*bibErrMsg = NULL;
  char		errBuf[256];
  BibFileError	bibErrFlag = BIBFILE_ER_NONE;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
  	  ("RecFileSecAtrbRead FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )atrb, version,
	   (unsigned long )fP, (unsigned long)eMsg));
  switch(version)
  {
    case 1:
      atrb->trMode = REC_TRMODE_REL;
      atrb->currentItem = NULL;
      break;
    case 2:
      atrb->currentItem = NULL;
      while(((bibErrFlag = BibFileRecordRead(&atRecord, &bibErrMsg,
					     fP)) == BIBFILE_ER_NONE) &&
	    atRecord && atRecord->name &&
	    strcmp(atRecord->name, "Attributes"))
      {
	if(lastRecord)
	{
	  BibFileRecordFree(&lastRecord);
	}
	lastRecord = atRecord;
      }
      if(bibErrMsg)
      {
	AlcFree(bibErrMsg);
      }
      if((bibErrFlag != BIBFILE_ER_NONE) ||
	 (atRecord == NULL) || (atRecord->name == NULL) ||
	 (atRecord->field == NULL) ||
	 ((BibFileFieldParseFmt(atRecord->field,
	 			&trModeS, "%s", "TransformMode",
				NULL) == 0)) ||
	 (trModeS == NULL) ||
         (WlzStringMatchValue((int *)&(atrb->trMode), trModeS,
		  recFileAtrbValueS[REC_FILE_SECATRB_REL], REC_TRMODE_REL,
		  recFileAtrbValueS[REC_FILE_SECATRB_ABS], REC_TRMODE_ABS,
		  NULL) == 0))
      {
	errFlag = REC_ERR_READ;
	(void )sprintf(errBuf,
		       "failed to read list attributes.");
      }
      if(trModeS)
      {
        AlcFree(trModeS);
      }
      if(atRecord)
      {
	BibFileRecordFree(&atRecord);
      }
      REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
      	      ("RecFileSecAtrbRead 01 %d\n",
	       atrb->trMode));
      break;
    default:
      errFlag = REC_ERR_READ;
      (void )sprintf(errBuf,
		     "unsupported file version %d.",
		     version);
     break;
  }
  if((errFlag != REC_ERR_NONE) && (*eMsg == NULL))
  {
    *eMsg = AlcStrDup(errBuf);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecAtrbRead FX 0x%lx\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads the reconstruct reconstruction information.
* \param	rec			Reconstruct information data structure
* 					to be filled in.
* \param	version			Version number.
* \param	fP			File pointer.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	 RecFileSecRecRead(RecReconstruction *rec, int version,
				   FILE *fP, char **eMsg)
{
  BibFileRecord	*rcRecord = NULL,
  		*lastRecord = NULL;
  char		*fileS = NULL,
  		*formatS = NULL,
		*scaleS = NULL,
		*srcBoxS = NULL,
		*dstBoxS = NULL,
		*bibErrMsg = NULL;
  char		errBuf[256];
  RecError	errFlag = REC_ERR_NONE;
  BibFileError	bibErrFlag = BIBFILE_ER_NONE;

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecRecRead FE 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )rec, version,
	   (unsigned long )fP, (unsigned long)eMsg));
  switch(version)
  {
    case 1:
      RecSecRecSetDefaults(rec);
      break;
    case 2:
      RecSecRecSetDefaults(rec);
      while(((bibErrFlag = BibFileRecordRead(&rcRecord, &bibErrMsg,
					     fP)) == BIBFILE_ER_NONE) &&
	    rcRecord && rcRecord->name &&
	    strcmp(rcRecord->name, "Reconstruction"))
      {
	if(lastRecord)
	{
	  BibFileRecordFree(&lastRecord);
	}
	lastRecord = rcRecord;
      }
      if(bibErrMsg)
      {
	AlcFree(bibErrMsg);
      }
      if(bibErrFlag != BIBFILE_ER_NONE)
      {
        errFlag = REC_ERR_READ;
      }
      else
      {
        if(rcRecord && rcRecord->name && rcRecord->field &&
	   (BibFileFieldParseFmt(rcRecord->field,
				 &fileS, "%s", "File",
				 &formatS, "%s", "Format",
				 &(rec->gaussFlt), "%d", "Filter",
				 &(rec->fastSam), "%d", "FastSampling",
				 &(rec->intScale), "%d", "IntScale",
				 &(rec->greedy), "%d", "Greedy",
				 &scaleS, "%s", "Scale",
				 &(rec->matchHst), "%d", "MatchHistograms",
				 &(rec->matchSec), "%d", "MatchSection",
				 &(rec->clipSrc), "%d", "ClipSrc",
				 &srcBoxS, "%s", "SrcBox",
				 &(rec->clipDst), "%d", "ClipDst",
				 &dstBoxS, "%s", "DstBox",
				 NULL) > 0))
	{
	  if(fileS)
	  {
	    rec->fileName = fileS;
	    fileS = NULL;
	  }
	  if(formatS)
	  {
	    if(errFlag == REC_ERR_NONE)
	    {
	      rec->fileFormat = WlzEffStringExtToFormat(formatS);
	      if(rec->fileFormat == WLZEFF_FORMAT_NONE)
	      {
		errFlag = REC_ERR_READ;
		(void )sprintf(errBuf,
			       "unsupported reconstruction file format %s.",
			       formatS);
	      }
	    }
	    AlcFree(formatS);
	  }
	  if(scaleS)
	  {
	    if(errFlag == REC_ERR_NONE)
	    {
	      if(sscanf(scaleS, "%lg %lg %lg",
			&(rec->scale.vtX), &(rec->scale.vtY),
			&(rec->scale.vtZ)) != 3)
	      {
		errFlag = REC_ERR_READ;
		(void )sprintf(errBuf,
			       "failed to parse reconstruction scale.");
	      }
	    }
	    AlcFree(scaleS);
	  }
	  if(srcBoxS)
	  {
	    if((errFlag == REC_ERR_NONE) && rec->clipSrc)
	    {
	      if(sscanf(srcBoxS, "%lg %lg %lg %lg %lg %lg",
	      		&(rec->srcBox.xMin), &(rec->srcBox.yMin),
			&(rec->srcBox.zMin), &(rec->srcBox.xMax),
			&(rec->srcBox.yMax), &(rec->srcBox.zMax)) != 6)
	      {
	        errFlag = REC_ERR_READ;
		(void )sprintf(errBuf,
			       "failed to parse reconstruction source box.");
	      }
	    }
	    AlcFree(srcBoxS);
	  }
	  if(dstBoxS)
	  {
	    if((errFlag == REC_ERR_NONE) && rec->clipDst)
	    {
	      if(sscanf(dstBoxS, "%lg %lg %lg %lg %lg %lg",
	      		&(rec->dstBox.xMin), &(rec->dstBox.yMin),
			&(rec->dstBox.zMin), &(rec->dstBox.xMax),
			&(rec->dstBox.yMax), &(rec->dstBox.zMax)) != 6)
	      {
	        errFlag = REC_ERR_READ;
		(void )sprintf(errBuf,
			       "failed to parse reconstruction source box.");
	      }
	    }
	    AlcFree(dstBoxS);
	  }
	}
      }
      if(rcRecord)
      {
	BibFileRecordFree(&rcRecord);
      }
      REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
      	      ("RecFileSecRecRead 01 %d 0x%lx %s %d "
	       "%d %d %d %d "
	       "{%g %g %g} "
	       "%d %d "
	       "%d {%g %g %g %g %g %g} "
	       "%d {%g %g %g %g %g %g}\n",
	       errFlag, (unsigned long )(rec->obj),
	       (rec->fileName)? rec->fileName: "(null)", rec->fileFormat,
	       rec->gaussFlt, rec->fastSam, rec->intScale, rec->greedy,
	       rec->scale.vtX, rec->scale.vtY, rec->scale.vtZ,
	       rec->matchHst, rec->matchSec,
	       rec->clipSrc, rec->srcBox.xMin, rec->srcBox.yMin,
	       rec->srcBox.zMin, rec->srcBox.xMax,
	       rec->srcBox.yMax, rec->srcBox.zMax,
	       rec->clipDst, rec->dstBox.xMin, rec->dstBox.yMin,
	       rec->dstBox.zMin, rec->dstBox.xMax,
	       rec->dstBox.yMax, rec->dstBox.zMax));
      break;
    default:
      errFlag = REC_ERR_READ;
      (void )sprintf(errBuf,
		     "unsupported file version %d.",
		     version);
     break;
  }
  if((errFlag != REC_ERR_NONE) && (*eMsg == NULL))
  {
    *eMsg = AlcStrDup(errBuf);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecFileSecRecRead FX 0x%lx\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads the body of a bibtex style list of section
*		transform records from the given file stream, parsing
*		them and building a linked list of sections.
* \param	secList			Section list.
* \param	numSec			Ptr for number of sections in
*					the list.
* \param	version			Version number.
* \param	fP			File pointer.
* \param	eMsg			Destination pointer for messages.
*/
static RecError	RecFileSecSecListRead(HGUDlpList *secList, int *numSec,
				      int version, FILE *fP, char **eMsg)
{
  int		doneFlag = 0,
		secCount = 0,
		recCount = 0;
  RecSection	*sec = NULL;
  RecError	errFlag = REC_ERR_NONE;
  char		errBuf[256];

  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecSecListRead FE 0x%lx 0x%lx %d 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )numSec,
	   version, (unsigned long )fP, (unsigned long )eMsg));
  if((secList == NULL) || (numSec == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    switch(version)
    {
      case 1: /* FALLTHROUGH */
      case 2:
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
	      errFlag = REC_ERR_MALLOC; 		 /* May not be true! */
	    }
	    else
	    {
	      sec = NULL;
	      ++secCount;
	    }
	  }
	  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
		  ("RecFileSecSecListRead 01 0x%lx %d %d %d\n",
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
			"at least two sections required in a section list.");
	}
	if(errFlag == REC_ERR_NONE)
	{
	  *numSec = secCount;
	}
	break;
      default:
	errFlag = REC_ERR_READ;
	(void )sprintf(errBuf,
		       "unsupported file version %d.",
		       version);
        break;
    }
  }
  if(errFlag != REC_ERR_NONE)
  {
    if(*eMsg == NULL)
    {
      *eMsg = AlcStrDup(errBuf);
    }
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileSecSecListRead FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads a woolz object from the the given file.
* \param	fP			File pointer.
* \param	obj			Destination pointer fow Woolz
*					object.
*/
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

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes a woolz object to the given file.
* \param	fP			File pointer.
* \param	obj			Woolz object to be written.
*/
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

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Reads and parses the next BibFile section record from
*               the given file stream. The record is parsed into the
*               given section data structure.
* \param	sec			Section list (already allocated).
* \param	doneFlag		Destination pointer for a non-error,
*					done flag.
* \param	recNum			Destination pointer for a record
*					counter, which may be used to count
*					the number of words parsed.
* \param	fP			File pointer.
* \param	eMsg			Destination pointer for messages.
*/
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
  const char	intFmt[] = "%d",
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
                                  &secIter, (char *)intFmt, "Iterations",
                                  &secCor, (char *)dblFmt, "Correlation",
                                  NULL);
      if((BibFileFieldParseFmt(secRecord->field,
                               &secFile, (char *)strFmt, "File",
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

/*!
* \return
* \ingroup	Reconstruct
* \brief	Creates fields using the bibtex file syntax for the
*               given transform and allocate storage as required.
* 		Global bibFileTransfFieldS is used fro transform field
*		names.
* \param	transf			Given transform.
*/
BibFileField	*RecFileTransfToField(WlzAffineTransform *transf)
{
  BibFileField	*field = NULL;
  WlzAffineTransformPrim prim;
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
    (void )WlzAffineTransformPrimGet(transf, &prim);
    (void )sprintf(tTypeS, "%d", transf->type);
    (void )sprintf(tTxS, "%g", prim.tx);
    (void )sprintf(tTyS, "%g", prim.ty);
    (void )sprintf(tTzS, "%g", prim.tz);
    (void )sprintf(tScaleS, "%g", prim.scale);
    (void )sprintf(tThetaS, "%g", prim.theta);
    (void )sprintf(tPhiS, "%g", prim.phi);
    (void )sprintf(tAlphaS, "%g", prim.alpha);
    (void )sprintf(tPsiS, "%g", prim.psi);
    (void )sprintf(tXsiS, "%g", prim.xsi);
    (void )sprintf(tInvertS, "%d", prim.invert);
    field = BibFileFieldMakeVa(
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_TYPE], tTypeS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_TX], tTxS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_TY], tTyS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_TZ], tTzS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_SCALE], tScaleS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_THETA], tThetaS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_PHI], tPhiS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_ALPHA], tAlphaS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_PSI], tPsiS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_XSI], tXsiS,
		 (char *)recFileTransfFieldS[REC_FILE_TRANSF_INVERT], tInvertS,
		 NULL);
  }
  REC_DBG((REC_DBG_FILE|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecFileTransfToField FX 0x%lx\n",
	   (unsigned long )field));
  return(field);
}

/*!
* \return
* \ingroup	Reconstruct
* \brief	Parses bibtex file syntax fields for the transform
*               fields and allocate storage for the new transform as
*               required.
* 		Note: The transform's linkcount is not set by this function.
* \param	field			Given field which is at the top
*                                       of the transform fields to be
*                                       parsed.
* \param	defaultFlg		Allow default transform fields
*                                       if non zero.
*/
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
  const char	intFmt[] = "%d",
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
			   &tType,   (char *)intFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_TYPE],
			   &tX,      (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_TX],
			   &tY,      (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_TY],
			   &tZ,      (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_TZ],
			   &tScale,  (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_SCALE],
			   &tTheta,  (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_THETA],
			   &tPhi,    (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_PHI],
			   &tAlpha,  (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_ALPHA],
			   &tPsi,    (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_PSI],
			   &tXsi,    (char *)dblFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_XSI],
			   &tInvert, (char *)intFmt,
			   (char *)recFileTransfFieldS[REC_FILE_TRANSF_INVERT],
			   NULL);
   REC_DBG((REC_DBG_FILE|REC_DBG_LVL_2),
   	   ("RecFileFieldToTransf 01 %d %d %f %f %f %f %f %f %f %f %f %d\n",
	    count, tType, tX, tY, tZ, tScale, tTheta, tPhi, tAlpha, tPsi, tXsi,
	    tInvert));
   if((defaultFlg || (count == REC_FILE_TRANSF_MAX)) &&
      (tType == WLZ_TRANSFORM_2D_AFFINE))
   {
     transf = WlzAffineTransformFromPrimVal(tType, tX, tY, 0.0,
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
