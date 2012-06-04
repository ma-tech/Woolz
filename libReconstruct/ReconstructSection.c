#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructSection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructSection.c
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
* \brief	Provides functions for manipulationg sections and section lists
* 		for the Reconstruct library.
* \ingroup	Reconstruct
*/
#include <Reconstruct.h>
#include <string.h>

typedef struct
{
  int		indexCur;
  int		indexInc;
} RecSecListIndexArgs;

typedef struct
{
  unsigned int	checkMask;
  unsigned int	invalidMask;
  RecSection	*prevSec;
} RecSecListInvalidArgs;

static int			RecSecSortFnSignCmpDbl(
				  double d0,
				  double d1);
static int			RecSecSortFnTransf(
				  void *entry0,
				  void *entry1,
				  unsigned int mask);
static int			RecSecSortFnIndex(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnIterations(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnCorrel(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnImageFile(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfTx(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfTy(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfTz(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfScale(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfTheta(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfPhi(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfAlpha(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfPsi(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfXsi(
				  void *entry0,
				  void *entry1);
static int			RecSecSortFnTransfInvert(
				  void *entry0,
				  void *entry1);
static int			RecSecListIndexFn(
				  HGUDlpList *secList,
				  HGUDlpListItem *item,
				  void *data);
static int			RecSecListInvalidFn(
				  HGUDlpList *secList,
				  HGUDlpListItem *item,
				  void *data);
static int			RecSecCumTransfClearItFn(
				  HGUDlpList *secList,
				  HGUDlpListItem *secItem,
				  void *itData);
static int			RecSecCumTransfSetItFn(
				  HGUDlpList *secList,
				  HGUDlpListItem *secItem,
				  void *itData);
static int			RecSecFindItemIndexItFn(
				  HGUDlpList *secList,
				  HGUDlpListItem *secItem,
				  void *data);

/*!
* \return	Section with incremented link count, if non NULL.
* \ingroup	Reconstruct
* \brief	Increments the given sections link count to ease registration
* 		section assignment.
* \param	sec			Section to be assigned.
*/
RecSection	*RecSecAssign(RecSection *sec)
{
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecAssign FE 0x%lx\n",
	   (unsigned long )sec));
  if(sec)
  {
    ++(sec->linkcount);
    REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
	    ("RecSecAssign 01 %d\n",
	     sec->linkcount));
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecAssign FX 0x%lx\n",
	   (unsigned long )sec));
  return(sec);
}

/*!
* \return	Duplicated section, NULL on error.
* \ingroup	Reconstruct
* \brief	Duplicates the given section.
* \param	sec			Section to be duplicated.
*/
RecSection	*RecSecDup(RecSection *sec)
{
  RecSection	*newSec = NULL;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecDup FE 0x%lx\n",
	   (unsigned long )sec));
  if(sec)
  {
    if(((newSec = (RecSection *)AlcCalloc(1, sizeof(RecSection))) != NULL) &&
       sec->imageFile &&
       ((newSec->imageFile = AlcStrDup(sec->imageFile)) != NULL) &&
       sec->transform &&
       ((newSec->transform = WlzAssignAffineTransform(
			 WlzAffineTransformFromMatrix(sec->transform->type,
       					      sec->transform->mat, &wlzErr),
					      NULL)) != NULL) &&
       (wlzErr == WLZ_ERR_NONE))
    {
      newSec->index = sec->index;
      newSec->iterations = sec->iterations;
      newSec->correl = sec->correl;
    }
    else
    {
      if(newSec)
      {
        if(newSec->imageFile)
	{
	  AlcFree(newSec->imageFile);
	}
        AlcFree(newSec);
        newSec = NULL;
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecDup FX 0x%lx\n",
	   (unsigned long )newSec));
  return(newSec);
}

/*!
* \return	New registration section.
* \ingroup	Reconstruct
* \brief	Makes a registration section using the given member values.
* \param	index			Section index.
* \param	iterations		Number of iterations to find
*                                       section transform.
* \param	correlation		Section correlation value.
* \param	imageFile		Image file path, this is duplicated
*					so that the original may be freed
*					The image file path must not be NULL.
* \param	transform		Section transform, if NULL an identity
* 					transform is created.
* \param	obj			Woolz object corresponding to the given
* 					image file. This may be NULL without
* 					causing the object to be read from the
* 					associated file.
*/
RecSection	*RecSecMake(int index, int iterations, double correlation,
			    char *imageFile,
			    WlzAffineTransform *transform, WlzObject *obj)
{
  RecSection	*sec = NULL;
  char		*newImageFile = NULL;
  WlzAffineTransform	*newTransform = NULL;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecMake FE %d %d %g 0x%lx 0x%lx 0x%lx\n",
	   index, iterations, correlation,
	   (unsigned long )imageFile, (unsigned long )transform,
	   (unsigned long )obj));
  if(imageFile)
  {
    newImageFile = AlcStrDup(imageFile);
  }
  if(newImageFile && (transform == NULL))
  {
    newTransform = WlzAffineTransformFromPrimVal(WLZ_TRANSFORM_2D_AFFINE,
    					         0.0, 0.0, 0.0,
					         1.0, 0.0, 0.0,
					         0.0, 0.0, 0.0, 0, &wlzErr);
  }
  if(newImageFile && (newTransform || transform) && (wlzErr == WLZ_ERR_NONE))
  {
    sec = (RecSection *)AlcMalloc(sizeof(RecSection));
  }
  if(sec == NULL)
  {
    if(newImageFile)
    {
      AlcFree(newImageFile);
    }
    if(newTransform)
    {
      WlzFreeAffineTransform(newTransform);
    }
  }
  else
  {
    sec->linkcount = 0;
    sec->index = index;
    sec->iterations = iterations;
    sec->correl = correlation;
    sec->obj = WlzAssignObject(obj, NULL);
    sec->imageFile = newImageFile;
    sec->transform = WlzAssignAffineTransform(transform? transform:
    							 newTransform, NULL);
    sec->transObj = NULL;
    sec->cumTransform = NULL;
    sec->cumTransObj = NULL;
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecMake FX 0x%lx\n",
	   (unsigned long )sec));
  return(sec);
}

/*!
* \return	New registration section.
* \ingroup	Reconstruct
* \brief	Makes an empty registration section.
* \param	index			Section index.
*/
RecSection	*RecSecMakeEmpty(int index)
{
  RecSection	*sec = NULL;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecMakeEmpty FE %d\n",
	   index));
  sec = RecSecMake(index, 0, 0.0, "Empty", NULL, NULL);
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecMakeEmpty FX 0x%lx\n",
	   (unsigned long )sec));
  return(sec);
}

/*!
* \ingroup	Reconstruct
* \brief	Free's the given section.
* \param	sec			Section to free.
*/
void		RecSecFree(RecSection *sec)
{
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecFree FE 0x%lx\n",
	   (unsigned long )sec));
  if(sec)
  {
    REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
	    ("RecSecFree 01 %d\n",
	     sec->linkcount));
    if(--(sec->linkcount) <= 0)
    {
      if(sec->imageFile)
      {
	AlcFree(sec->imageFile);
      }
      if(sec->transform)
      {
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
		("RecSecFree 02 %d\n",
		 sec->transform->linkcount));
	(void )WlzFreeAffineTransform(sec->transform);
      }
      if(sec->obj)
      {
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
		("RecSecFree 03 %d\n",
		 sec->obj->linkcount));
	(void )WlzFreeObj(sec->obj);
      }
      if(sec->transObj)
      {
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
		("RecSecFree 04 %d\n",
		 sec->transObj->linkcount));
        (void )WlzFreeObj(sec->transObj);
      }
      if(sec->cumTransform)
      {
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
		("RecSecFree 05 %d\n",
		 sec->cumTransform->linkcount));
	(void )WlzFreeAffineTransform(sec->cumTransform);
      }
      if(sec->cumTransObj && (sec->cumTransObj->linkcount > 0))
      {
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_3),
		("RecSecFree 05 %d\n",
		 sec->cumTransObj->linkcount));
        (void )WlzFreeObj(sec->cumTransObj);
      }
      AlcFree(sec);
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecFree FX\n"));
}

/*!
* \return	Non zero if section is empty.
* \ingroup	Reconstruct
* \brief	Checks if this is an empty section by comparing the image file
* 		name with the case insensitive string 'empty'.
* \param	sec			Section to be checked.
*/
int		RecSecIsEmpty(RecSection *sec)
{
  int		isEmpty = 0;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecIsEmpty FE 0x%lx\n",
	   (unsigned long )sec));
  if(sec)
  {
    if((sec->imageFile == NULL) || (*(sec->imageFile) == '\0') ||
       (strcasecmp(sec->imageFile, "Empty") == 0))
    {
      isEmpty = 1;
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecIsEmpty FX %d\n",
	   isEmpty));
  return(isEmpty);
}

/*!
* \return	Next registration section, or NULL if at end of list.
* \ingroup	Reconstruct
* \brief	Given a section list and a current item pointer, this
*               function sets the given destination pointer to the next
*               item and returns that list item's section entry.
* \param	secList			Section list.
* \param	cItem			Current list item, NULL to
*                                       get first list entry.
* \param	nItemP			Destination pointer for the
*                                       next list item, will be set to
*                                       NULL at end of list.
* \param	skipEmpty		Skip over empty entries if flag
*                                       is non-zero.
*/
RecSection	*RecSecNext(HGUDlpList *secList, HGUDlpListItem *cItem,
			    HGUDlpListItem **nItemP, int skipEmpty)
{
  RecSection	*sec = NULL;
  HGUDlpListItem *item;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecNext FE 0x%lx 0x%lx 0x%lx %d\n",
	   (unsigned long )secList, (unsigned long )cItem,
	   (unsigned long )nItemP, skipEmpty));
  if(secList)
  {
    item = cItem;
    do
    {
      item = (item)? HGUDlpListNext(secList, item): HGUDlpListHead(secList);
      sec = (item)? (RecSection *)HGUDlpListEntryGet(secList, item): NULL;
    } while(skipEmpty && sec && RecSecIsEmpty(sec));
    if(nItemP)
    {
      *nItemP = item;
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecNext FX 0x%lx\n",
	   (unsigned long )sec));
  return(sec);
}

/*!
* \return	Previous registration section or NULL if at head of list.
* \ingroup	Reconstruct
* \brief	Given a section list and a current item pointer, this
*               function sets the given destination pointer to the
*               previous item and returns that list item's section
*               entry.
* \param	secList			Section list.
* \param	cItem			Current list item, NULL to get the
*					last list entry.
* \param	pItemP			Destination pointer for the
*					previous list item, will be set to NULL
*					at head of list.
* \param	skipEmpty		Skip over empty entries if flag
*                                       is non-zero.
*/
RecSection	*RecSecPrev(HGUDlpList *secList, HGUDlpListItem *cItem,
			    HGUDlpListItem **pItemP, int skipEmpty)
{
  RecSection	*sec = NULL;
  HGUDlpListItem *item;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecPrev FE 0x%lx 0x%lx 0x%lx %d\n",
	   (unsigned long )secList, (unsigned long )cItem,
	   (unsigned long )pItemP, skipEmpty));
  if(secList)
  {
    item = cItem;
    do
    {
      item = (item)? HGUDlpListPrev(secList, item): HGUDlpListTail(secList);
      sec = (item)? (RecSection *)HGUDlpListEntryGet(secList, item): NULL;
    } while(skipEmpty && sec && RecSecIsEmpty(sec));
    if(pItemP)
    {
      *pItemP = item;
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecPrev FX 0x%lx\n",
	   (unsigned long )sec));
  return(sec);
}

/*!
* \return	New section string.
* \ingroup	Reconstruct
* \brief	Creates a string from the given registration serial
*               section and bit mask for the fields required.
* \param	sec			Given section.
* \param	reqFields		Bit mask for fields required
*                                       in the strings.
* \param	eMsg			Destination pointer for messages.
*/
char		*RecSecToStr(RecSection *sec, unsigned int reqFields,
			     char **eMsg)
{
  unsigned long tUL0,
  		filePathLen,
  		strIdx;
  char		*secStr = NULL;
  WlzAffineTransformPrim prim;
  char		strBuf[1024]; 	   /* This must be big enough for any string */

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecToStr FE 0x%lx %d 0x%lx %d\n",
	   (unsigned long )sec, reqFields, (unsigned long )eMsg));
  if(sec)
  {
    strIdx = 0;
    strBuf[0] = '\0';
    if(reqFields & REC_SECMSK_INDEX)
    {
      (void )sprintf(strBuf + strIdx, "Index % 12d", sec->index);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_ITERATIONS)
    {
      (void )sprintf(strBuf + strIdx, ", Iterations % 12d",
		     sec->iterations);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_CORREL)
    {
      (void )sprintf(strBuf + strIdx, ", Correlation % 12.8g",
		     sec->correl);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_IMAGEFILE)
    {
      filePathLen = strlen(sec->imageFile);
      tUL0 = filePathLen % 10; 			       /* Pad out file paths */
      if(tUL0 == 0)
      {
        filePathLen += 10;
      }
      else
      {
        filePathLen += 20 - tUL0;
      }
      filePathLen = strlen(sec->imageFile) + 10;
      filePathLen -= filePathLen % 10;
      (void )sprintf(strBuf + strIdx, ", Imagefile %*s",
		     (int )filePathLen,
		     sec->imageFile);	 	 /* Dynamic string precision */
      strIdx = strlen(strBuf);
    }
    (void )WlzAffineTransformPrimGet(sec->transform, &prim);
    if(reqFields & REC_SECMSK_TRANSF_TX)
    {
      (void )sprintf(strBuf + strIdx, ", Tx % 12.8g",
		     prim.tx);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_TY)
    {
      (void )sprintf(strBuf + strIdx, ", Ty % 12.8g",
		     prim.ty);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_TZ)
    {
      (void )sprintf(strBuf + strIdx, ", Tz % 12.8g",
		     prim.tz);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_SCALE)
    {
      (void )sprintf(strBuf + strIdx, ", Tscale % 12.8g",
		     prim.scale);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_THETA)
    {
      (void )sprintf(strBuf + strIdx, ", Ttheta % 12.8g",
		     prim.theta);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_PHI)
    {
      (void )sprintf(strBuf + strIdx, ", Tphi % 12.8g",
		     prim.phi);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_ALPHA)
    {
      (void )sprintf(strBuf + strIdx, ", Talpha % 12.8g",
		     prim.alpha);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_PSI)
    {
      (void )sprintf(strBuf + strIdx, ", Tpsi % 12.8g",
		     prim.psi);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_XSI)
    {
      (void )sprintf(strBuf + strIdx, ", Txsi % 12.8g",
		     prim.xsi);
      strIdx = strlen(strBuf);
    }
    if(reqFields & REC_SECMSK_TRANSF_INVERT)
    {
      (void )sprintf(strBuf + strIdx, ", Tinvert %d",
		     (prim.invert)? 1: 0);
      strIdx = strlen(strBuf);
    }
    if(strIdx > 0)
    {
      secStr = AlcStrDup(strBuf);
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecToStr FX %s\n",
	   (secStr)? (secStr): ("null")));
  return(secStr);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Creates a list of strings from a registration serial
*               section list and a bit mask for the fields required.
*               This is NOT intended to be used for output to a file
*               (RecFileSecWrite() should be used). This function was
*               written to allow a user to be presented with a simple
*               list of serial sections within a GUI application.
* \param	strList			Destination pointer for list of
*                                       strings.
* \param	secList			Given section list.
* \param	numSec			Number of sections in secList.
* \param	eMsg			Destination pointer for messages.
* \param	reqFields		Bit mask for fields required
*                                       in the strings.
*/
RecError	RecSecListToStrList(char ***strList, HGUDlpList *secList,
				    int numSec, char **eMsg,
				    unsigned int reqFields)
{
  int		secIdx = 0;
  HGUDlpListItem *secItem;
  RecSection	*sec;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListToStrList FE 0x%lx 0x%lx %d 0x%lx %d\n",
	   (unsigned long )strList, (unsigned long )secList, numSec,
	   (unsigned long )eMsg, reqFields));
  if((*strList = (char **)AlcMalloc(sizeof(char **) *
  				 numSec)) != NULL)   /* Allocate string list */
  {
    secItem = HGUDlpListHead(secList);
    while(secItem &&
    	  ((sec = HGUDlpListEntryGet(secList, secItem)) != NULL) &&
    	  (secIdx < numSec) && (errFlag == REC_ERR_NONE))
    {
      if((*(*strList + secIdx) = RecSecToStr(sec, reqFields, eMsg)) == NULL)
      {
        errFlag = REC_ERR_MALLOC;
      }
      else
      {
	secItem = HGUDlpListNext(secList, secItem);
	++secIdx;
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListToStrList FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \ingroup	Reconstruct
* \brief	Sorts the given HGUDlpList section list using the given
*               section field mask.
* \param	secList			Given section list.
* \param	sortMask		Bit mask for the section field
*                                       to sort on.
*/
void		RecSecListSort(HGUDlpList *secList, unsigned int sortMask)
{
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListSort FE 0x%lx %d\n",
	   (unsigned long )secList, sortMask));
  sortMask &= REC_SECMSK_ALL;
  if(secList && sortMask)
  {
    switch(sortMask)
    {
      case REC_SECMSK_INDEX:
        (void )HGUDlpListSort(secList,  RecSecSortFnIndex);
	break;
      case REC_SECMSK_ITERATIONS:
        (void )HGUDlpListSort(secList,  RecSecSortFnIterations);
	break;
      case REC_SECMSK_CORREL:
        (void )HGUDlpListSort(secList,  RecSecSortFnCorrel);
	break;
      case REC_SECMSK_IMAGEFILE:
        (void )HGUDlpListSort(secList,  RecSecSortFnImageFile);
	break;
      case REC_SECMSK_TRANSF_TX:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfTx);
	break;
      case REC_SECMSK_TRANSF_TY:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfTy);
	break;
      case REC_SECMSK_TRANSF_TZ:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfTz);
	break;
      case REC_SECMSK_TRANSF_SCALE:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfScale);
	break;
      case REC_SECMSK_TRANSF_THETA:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfTheta);
	break;
      case REC_SECMSK_TRANSF_PHI:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfPhi);
	break;
      case REC_SECMSK_TRANSF_ALPHA:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfAlpha);
	break;
      case REC_SECMSK_TRANSF_PSI:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfPsi);
	break;
      case REC_SECMSK_TRANSF_XSI:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfXsi);
	break;
      case REC_SECMSK_TRANSF_INVERT:
        (void )HGUDlpListSort(secList,  RecSecSortFnTransfInvert);
	break;
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListSort FX\n"));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given a pair of doubles, returns an integer with the
*               sign of their difference. Ie sign(d0 - d1).
*               Used by some of the sort functions.
* \param	d0			First of the doubles.
* \param	d1			Second of the doubles.
*/
static int	RecSecSortFnSignCmpDbl(double d0, double d1)
{
  int		sgnCmp = 0;

  if(d0 < d1)
  {
    sgnCmp = -1;
  }
  else if(d0 > d1)
  {
    sgnCmp = 1;
  }
  return(sgnCmp);
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               indicies for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnIndex(void *entry0, void *entry1)
{
  int		tI0 = 0;
  RecSection	*sec0,
		*sec1;

  if(entry0 && entry1)
  {
    sec0 = (RecSection *)entry0,
    sec1 = (RecSection *)entry1;
    tI0 = sec1->index - sec0->index;
  }
  return(tI0);
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               iterations for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnIterations(void *entry0, void *entry1)
{
  int		tI0 = 0;
  RecSection	*sec0,
		*sec1;

  if(entry0 && entry1)
  {
    sec0 = (RecSection *)entry0,
    sec1 = (RecSection *)entry1;
    tI0 = sec1->iterations - sec0->iterations;
  }
  return(tI0);
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               cross-correlations for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnCorrel(void *entry0, void *entry1)
{
  int		tI0 = 0;
  RecSection	*sec0,
		*sec1;

  if(entry0 && entry1)
  {
    sec0 = (RecSection *)entry0,
    sec1 = (RecSection *)entry1;
    tI0 = RecSecSortFnSignCmpDbl(sec1->correl, sec0->correl);
  }
  return(tI0);
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               image file names for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnImageFile(void *entry0, void *entry1)
{
  int		tI0 = 0;
  RecSection	*sec0,
		*sec1;

  if(entry0 && entry1)
  {
    sec0 = (RecSection *)entry0,
    sec1 = (RecSection *)entry1;
    if(sec0->imageFile && sec1->imageFile)
    {
      tI0 = strcmp(sec1->imageFile, sec0->imageFile);
    }
  }
  return(tI0);
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be
*               pointers to registration sections, compares the sections
*               transform x components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfTx(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_TX));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform y components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfTy(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_TY));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform z components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfTz(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_TZ));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform scale components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfScale(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_SCALE));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform theta components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfTheta(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_THETA));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform phi components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfPhi(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_PHI));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform alpha components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfAlpha(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_ALPHA));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform psi components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfPsi(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_PSI));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform xsi components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfXsi(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_XSI));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, compares the sections
*               transform invert components for HGUDlpListSort().
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
*/
static int	RecSecSortFnTransfInvert(void *entry0, void *entry1)
{
  return(RecSecSortFnTransf(entry0, entry1, REC_SECMSK_TRANSF_INVERT));
}

/*!
* \return	Signed value as for sort().
* \ingroup	Reconstruct
* \brief	Given two HGUDlpList entries which are known to be pointers
*               to registration sections, and a registration section
*               field mask, compares the sections transform fields
*               for RecSecSortFnTransfXXX() functions.
* \param	entry0			First list entry.
* \param	entry1			Second list entry.
* \param	mask
*/
static int	RecSecSortFnTransf(void *entry0, void *entry1,
				   unsigned int mask)
{
  int		tI0 = 0;
  RecSection	*sec0,
		*sec1;
  WlzAffineTransform	*transf0,
  		*transf1;
  WlzAffineTransformPrim prim0,
  		prim1;

  if(entry0 && entry1)
  {
    sec0 = (RecSection *)entry0,
    sec1 = (RecSection *)entry1;
    transf0 = sec0->transform;
    transf1 = sec1->transform;
    if(transf0 && transf1)
    {
      (void )WlzAffineTransformPrimGet(transf0, &prim0);
      (void )WlzAffineTransformPrimGet(transf1, &prim1);
      switch(mask)
      {
        case REC_SECMSK_TRANSF_TX:
          tI0 = RecSecSortFnSignCmpDbl(prim1.tx, prim0.tx);
	  break;
	case REC_SECMSK_TRANSF_TY:
          tI0 = RecSecSortFnSignCmpDbl(prim1.ty, prim0.ty);
	  break;
	case REC_SECMSK_TRANSF_TZ:
          tI0 = RecSecSortFnSignCmpDbl(prim1.tz, prim0.tz);
	  break;
	case REC_SECMSK_TRANSF_SCALE:
	  tI0 = prim1.scale - prim0.scale;
	  break;
	case REC_SECMSK_TRANSF_THETA:
          tI0 = RecSecSortFnSignCmpDbl(prim1.theta, prim0.theta);
	  break;
	case REC_SECMSK_TRANSF_PHI:
          tI0 = RecSecSortFnSignCmpDbl(prim1.phi, prim0.phi);
	  break;
	case REC_SECMSK_TRANSF_ALPHA:
          tI0 = RecSecSortFnSignCmpDbl(prim1.alpha, prim0.alpha);
	  break;
	case REC_SECMSK_TRANSF_PSI:
          tI0 = RecSecSortFnSignCmpDbl(prim1.psi, prim0.psi);
	  break;
	case REC_SECMSK_TRANSF_XSI:
          tI0 = RecSecSortFnSignCmpDbl(prim1.xsi, prim0.xsi);
	  break;
	case REC_SECMSK_TRANSF_INVERT:
	  tI0 = prim1.invert - prim0.invert;
	  break;
      }
    }
  }
  return(tI0);
}

/*!
* \ingroup	Reconstruct
* \brief	Sets the indicies for the sections in the given section
*               list using the given first index value and the given
*               index increment.
* \param	secList			Given section list.
* \param	indexFirst		First index value.
* \param	indexInc		Index increment between the
*                                       sections.
*/
void		RecSecListIndiciesSet(HGUDlpList *secList,
				      int indexFirst, int indexInc)
{
  HGUDlpListItem *head;
  RecSecListIndexArgs args;

  args.indexCur = indexFirst;
  args.indexInc = indexInc;
  if(secList)
  {
    head = HGUDlpListHead(secList);
    (void )HGUDlpListIterate(secList, head, HGU_DLPLIST_DIR_TOTAIL,
			     RecSecListIndexFn, (void *)&args);
  }
}

/*!
* \return	Non-zero for itteration to continue.
* \ingroup	Reconstruct
* \brief	Iterated function for section lists associated with
*               RecSecListIndiciesSet() which sets the section index
*               of the given item.
* \param	secList			Given section list.
* \param	item			List item with section entry.
* \param	data			Used to pass other args.
*/
static int	RecSecListIndexFn(HGUDlpList *secList, HGUDlpListItem *item,
				  void *data)
{
  int		continueFlag = 0;
  RecSecListIndexArgs *args;
  RecSection	*sec;

  if(secList && item && data)
  {
    args = (RecSecListIndexArgs *)data;
    sec = (RecSection *)HGUDlpListEntryGet(secList, item);
    sec->index = args->indexCur;
    args->indexCur += args->indexInc;
    continueFlag = 1;
  }
  return(continueFlag);
}

/*!
* \return	Registration section mask with any fields that are INVALID set.
* \ingroup	Reconstruct
* \brief	Checks the registration section list for validity.
*               Only section fields which are indicated in the given
*               section mask are checked for validity.
* \param	secList			Given section list.
* \param	invalidItem		Destination pointer for the
*                                       first invalid section list item
*                                       found, maybe NULL if not required.
* \param	mask			Mask specifying which fields
*                                       are to be checked for validity.
*/
unsigned int	 RecSecListInvalid(HGUDlpList *secList,
				   HGUDlpListItem **invalidItem,
				   unsigned int mask)
{
  HGUDlpListItem *item;
  RecSecListInvalidArgs args;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListInvalid FE 0x%lx 0x%lx %d\n",
	   (unsigned long )secList, (unsigned long )invalidItem,
	   (unsigned long )mask));
  args.checkMask = mask;
  args.invalidMask = REC_SECMSK_NONE;
  args.prevSec = NULL;
  if(secList)
  {
    item = HGUDlpListHead(secList);
    item = HGUDlpListIterate(secList, item, HGU_DLPLIST_DIR_TOTAIL,
			     RecSecListInvalidFn, (void *)&args);
    if(args.prevSec)
    {
      RecSecFree(args.prevSec);
    }
    if(invalidItem != NULL)
    {
      if(args.invalidMask == REC_SECMSK_NONE)
      {
        *invalidItem = NULL;
      }
      else
      {
        *invalidItem = item;
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecListInvalid FX %d\n",
	   (unsigned long )(args.invalidMask)));
  return(args.invalidMask);
}

/*!
* \return	Non-zero to continue iteration.
* \ingroup	Reconstruct
* \brief	Iterated function for section lists associated with
*               RecSecListInvalid() which checks the validity of the
*               given section item.
* \param	secList			Given section list.
* \param	item			List item with section entry.
* \param	data			Used to pass other args.
*/
static int	RecSecListInvalidFn(HGUDlpList *secList, HGUDlpListItem *item,
				    void *data)
{
  int		continueFlag = 0;
  RecSecListInvalidArgs *args;
  RecSection	*sec;

  if(secList && item && data)
  {
    args = (RecSecListInvalidArgs *)data;
    sec = (RecSection *)HGUDlpListEntryGet(secList, item);
    if(args->checkMask & REC_SECMSK_INDEX)
    {
      if((sec->index < 0) ||
         ((args->prevSec) && (sec->index <= args->prevSec->index)))
      {
	args->invalidMask |= REC_SECMSK_INDEX;
      }
    }
    if(args->checkMask & REC_SECMSK_ITERATIONS)
    {
      if(sec->iterations < 0)
      {
        args->invalidMask |= REC_SECMSK_ITERATIONS;
      }
    }
    if(args->checkMask & REC_SECMSK_CORREL)
    {
      if(sec->correl < 0.0)
      {
        args->invalidMask |=  REC_SECMSK_CORREL;
      }
    }
    if(args->checkMask & REC_SECMSK_IMAGEFILE)
    {
      if((sec->imageFile == NULL) || (strlen(sec->imageFile) == 0))
      {
        args->invalidMask |= REC_SECMSK_IMAGEFILE;
      }
    }
    if(args->prevSec)
    {
      RecSecFree(args->prevSec);
    }
    args->prevSec = RecSecDup(sec);
    if(args->invalidMask == REC_SECMSK_NONE)
    {
      continueFlag = 1;
    }
  }
  return(continueFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Appends new section records to the given section list
*               using the given list of image file names. The section
*               indicies are set using the given first index and index
*               increment values. The new section records have identity
*               transforms, zero iteration count and cross correlation
*               values. The woolz image object is set to NULL.
* \param	secList			Given section list.
* \param	secMark			List item after which the new
*                                       sections should be appended, if
*                                       this is NULL the the new
*                                       sections are appended at the
*                                       end of the list.
* \param	fileVec			Vector of file name strings.
* \param	nFiles			Number of file names.
* \param	indexFirst		Index value of first section of
*                                       the list.
* \param	indexInc		Index increment value.
*/
RecError	 RecSecAppendListFromFiles(HGUDlpList *secList,
					   HGUDlpListItem *secMark,
					   char **fileVec, int nFiles,
					   int indexFirst, int indexInc)
{
  int		index;
  RecSection	*sec;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCreateListFromFiles FE 0x%lx 0x%lx 0x%lx %d %d %d\n",
	   (unsigned long )secList, (unsigned long )fileVec,
	   (unsigned long )secMark, nFiles, indexFirst, indexInc));
  if((secList == NULL) || (fileVec == NULL) || (nFiles <= 0) ||
     (nFiles > REC_MAX_LIST))
  {
    errFlag = REC_ERR_FUNC;
  }
  else
  {
    index = indexFirst;
    while((errFlag == REC_ERR_NONE) && *fileVec && (nFiles > 0))
    {
      if(((sec = RecSecMake(index, 0, 0.0, *fileVec, NULL, NULL)) != NULL) &&
         ((secMark = HGUDlpListAppend(secList, secMark, sec,
	 			      (void (*)(void *) )RecSecFree)) != NULL))
      {
	--nFiles;
	++fileVec;
        index += indexInc;
      }
      else
      {
        errFlag = REC_ERR_MALLOC;
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCreateListFromFiles FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Sets the cumulative transform field of each section in
*               the section list to NULL, from the given section to the
*               tail of the list. When set to NULL the cumulative
*               transforms then need recalculating prior to section
*               display, ie NULL is used to mark the cumulative
*               transforms as invalid/outdated.
*               If the given section item is NULL then all list items
*               are processed.
* \param	secList			Given section list.
* \param	secItem			Given section item.
*/
RecError	 RecSecCumTransfClear(HGUDlpList *secList,
			    	      HGUDlpListItem *secItem)
{
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCumTransfClear FE 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )secItem));
  if(secList)
  {
    (void )HGUDlpListIterate(secList, secItem, HGU_DLPLIST_DIR_TOTAIL,
    		      	     RecSecCumTransfClearItFn, NULL);
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCumTransfClear FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Non-zero to continue iteration.
* \ingroup	Reconstruct
* \brief	Iterated function which clears the cumulative transform
*               and cumulative transform transformed object fields of
*               the given section list item to be NULL.
* \param	secList			Given section list.
* \param	secItem			Given section item.
* \param	itData			Not used.
*/
static int	RecSecCumTransfClearItFn(HGUDlpList *secList,
					 HGUDlpListItem *secItem,
					void *itData)
{
  int		iterate = 0;
  RecSection	*sec;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_3),
          ("RecSecCumTransfClearItFn FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )secItem,
	   (unsigned long )itData));
  if(secList && secItem &&
     ((sec = (RecSection *)HGUDlpListEntryGet(secList, secItem)) != NULL))
  {
    if(sec->cumTransObj)
    {
      (void )WlzFreeObj(sec->cumTransObj);
      sec->cumTransObj = NULL;
    }
    if(sec->cumTransform)
    {
      (void )WlzFreeAffineTransform(sec->cumTransform);
      sec->cumTransform = NULL;
    }
    iterate = 1;
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_3),
          ("RecSecCumTransfClearItFn FX %d\n",
	   iterate));
  return(iterate);
}

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Calculates and sets the cumulative transform field of
*               each section in the section list from the head of
*               the list to the given section. If a section has a
*               non-NULL cumulative transform then it is assumed valid.
*               This function should be called AFTER the function
*               RecSecCumTransfClear has been called to set invalid
*               entries to NULL.
*               If the given section item is NULL then all entries in
*               the list will be processed.
* \param	secList			Given section list.
* \param	secItem			Given section item.
*/
RecError	 RecSecCumTransfSet(HGUDlpList *secList,
		 		   HGUDlpListItem *secItem)
{
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCumTransfSet FE 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )secItem));
  if(secList)
  {
    (void )HGUDlpListIterate(secList, NULL, HGU_DLPLIST_DIR_TOTAIL,
    		      	     RecSecCumTransfSetItFn, secItem);
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
          ("RecSecCumTransfSet FX %d\n",
	   errFlag));
  return(errFlag);
}

/*!
* \return	Non-zero to continue iteration.
* \ingroup	Reconstruct
* \brief	Iterated function which sets the cumulative transform
*               fields of a section list, working from the head of the
*               list towards the tail and stoping when the specified
*               list item has been processed.
* \param	secList			Given section list.
* \param	secItem			Given section item.
* \param	itData			Used to pass item to stop at.
*/
static int	RecSecCumTransfSetItFn(HGUDlpList *secList,
					 HGUDlpListItem *secItem,
					void *itData)
{
  int		iterate = 0;
  HGUDlpListItem *prevItem,
  		*stopItem;
  WlzAffineTransform	*cumTransf = NULL;
  RecSection	*curSec,
  		*prevSec;
  WlzAffineTransformPrim prim;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
          ("RecSecCumTransfSetItFn FE 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )secList, (unsigned long )secItem,
	   (unsigned long )itData));
  if(secList && secItem && itData)
  {
    stopItem = (HGUDlpListItem *)itData;
    curSec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
    if(curSec && curSec->transform)		/* Should always be true */
    {
      if(HGUDlpListItemIsHead(secList, secItem))
      {
	if(curSec->cumTransform == NULL)
	{
	  curSec->cumTransform = curSec->transform;
	  ++(curSec->cumTransform->linkcount);
	}
	cumTransf = curSec->cumTransform;
      }
      else
      {
	if(curSec->cumTransform == NULL)
	{
	  prevItem = HGUDlpListPrev(secList, secItem);
	  prevSec = (RecSection *)HGUDlpListEntryGet(secList, prevItem);
	  if(prevSec && prevSec->cumTransform)      /* Should always be true */
	  {
	    cumTransf = WlzAffineTransformProduct(curSec->transform,
	    					  prevSec->cumTransform,
						  NULL);
	    curSec->cumTransform = cumTransf;
	    ++(curSec->cumTransform->linkcount);
	  }
	}
	else
	{
	  cumTransf = curSec->cumTransform;
	}
      }
      REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	      ("RecSecCumTransfSetItFn 02 0x%lx\n",
	       (unsigned long )cumTransf));
      if(cumTransf)
      {
	(void )WlzAffineTransformPrimGet(cumTransf, &prim);
	REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_3),
		("RecSecCumTransfSetItFn 02 %d %g %g %g %g %g %g %g %g %g %d\n",
		 cumTransf->type,
		 prim.tx,
		 prim.ty,
		 prim.tz,
		 prim.scale,
		 prim.theta,
		 prim.phi,
		 prim.alpha,
		 prim.psi,
		 prim.xsi,
		 prim.invert));
	if(secItem != stopItem)
	{
	  iterate = 1;
	}
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
          ("RecSecCumTransfSetItFn FX %d\n",
	   iterate));
  return(iterate);
}

/*!
* \return	First item found with matching index or NULL if no matching
*		index (or on error).
* \ingroup	Reconstruct
* \brief	Finds the first item in the given section list with a
*               matching index and returns a pointer to the list item.
* \param	list			Section list to search.
* \param	from			Item to search first and from.
* \param	index			Section index to search for.
* \param	dir			Search direction.
*/
HGUDlpListItem	*RecSecFindItemIndex(HGUDlpList *list, HGUDlpListItem *from,
				     int index, HGUDlpListDirection dir)
{
  HGUDlpListItem *found = NULL;
  RecSection 	*sec;

  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecFindItemIndex FE 0x%lx 0x%lx %d %d\n",
	   (unsigned long )list, (unsigned long )from, index, (int )dir));
  if(list && 
     ((dir == HGU_DLPLIST_DIR_TOTAIL) || (dir == HGU_DLPLIST_DIR_TOHEAD)))
  {
    found = HGUDlpListIterate(list, from, dir, RecSecFindItemIndexItFn,
    			      (void *)&index);
    if(found)			    /* Check for indicies same and not error */
    {
      sec = (RecSection *)HGUDlpListEntryGet(list, found);
      if((sec == NULL) || (sec->index != index))
      {
        found = NULL;
      }
    }
  }
  REC_DBG((REC_DBG_SEC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecSecFindItemIndex FX 0x%lx\n",
	   (unsigned long )found));
  return(found);
}

/*!
* \return	Non zero if the item's index is equal to the given index.
* \ingroup	Reconstruct
* \brief	Iterated function for RecSecFindItemIndex().
* \param	secList			Given section list.
* \param	secItem			The given item.
* \param	data			Used to pass given index.
*/
static int	RecSecFindItemIndexItFn(HGUDlpList *secList,
					HGUDlpListItem *secItem,
				  	void *data)
{
  int		notEqual = 0;
  RecSection	*sec;

  if(secList && secItem)
  {
    sec = (RecSection *)HGUDlpListEntryGet(secList, secItem);
    if(sec)
    {
      notEqual = sec->index != *(int *)data;
    }
  }
  return(notEqual);
}

/*!
* \return	New secton list or NULL on error.
* \ingroup	Reconstruct
* \brief	Creates a new empty section list with default settings.
* \param	dstErr			Destination error pointer, may
*                                       be NULL.
*/
RecSectionList	*RecSecNewSectionList(RecError *dstErr)
{
  RecSectionList *newSecList = NULL;
  RecError	errFlag = REC_ERR_NONE;

  if((newSecList = (RecSectionList *)
  		   AlcCalloc(1, sizeof(RecSectionList))) == NULL)
  {
    errFlag = REC_ERR_MALLOC;
  }
  else
  {
    /* Set fields to default values. */
    newSecList->list = NULL;
    newSecList->attributes.trMode = REC_TRMODE_REL;
    newSecList->attributes.currentItem = NULL;
    RecSecRecSetDefaults(&(newSecList->reconstruction));
  }
  if(dstErr)
  {
    *dstErr = errFlag;
  }
  return(newSecList);
}

/*!
* \ingroup	Reconstruct
* \brief	Sets the reconstruction information to default values.
* \param	rec			Reconstruction information.
*/
void		RecSecRecSetDefaults(RecReconstruction *rec)
{
  if(rec)
  {
    rec->obj = NULL;
    rec->fileName = NULL;
    rec->fileFormat = WLZEFF_FORMAT_WLZ;
    rec->gaussFlt = 0;
    rec->fastSam = 0;
    rec->intScale = 0;
    rec->greedy = 0;
    rec->scale.vtX = 1.0;
    rec->scale.vtY = 1.0;
    rec->scale.vtZ = 1.0;
    rec->matchHst = 0;
    rec->clipSrc = 0;
    rec->clipDst = 0;
  }
}
