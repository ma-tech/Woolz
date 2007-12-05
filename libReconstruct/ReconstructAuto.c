#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _ReconstructAuto_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         ReconstructAuto.c
* \author       Bill Hill
* \date         April 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Provides functions for the automatic registration of
*		serial sections.
* \ingroup	Reconstruct
* \todo         -
* \bug          None known.
*/

#include <Reconstruct.h>
#include <string.h>

/*!
* \return	Non zero if registration fails.
* \ingroup	Reconstruct
* \brief	Performs the automatic registration of serial sections.
* \param	rCtrl			The registration control data
* 					structure.
* \param	ppCtrl			Pre-processing control data
*					structure.
* \param	secList			Section list.
* \param	cancelFlag		Cancel if flag pointed to is non-zero.
* \param	secFn			application supplied section update
*					function. This function is responsible
*					for replacing the section in the list,
*					it may also display it, etc, ....
* \param	secData			Application supplied data for section
* 					update function.
* \param	workFn			Application supplied work function.
* \param	workData		Application supplied data for the
*					work function.
* \param	eMsg			Pointer for error message strings.
*/
RecError	RecAuto(RecControl *rCtrl, RecPPControl *ppCtrl,
			HGUDlpList *secList, int *cancelFlag,
			RecSecUpdateFunction secFn, void *secData,
			RecWorkFunction workFn, void *workData,
			char **eMsg)
{
  RecSection	*oSec0,
		*oSec1,
		*nSec0,
		*nSec1;
  HGUDlpListItem *item;
  RecState	rState;
  static char	errMsgInvalidListStr[] =
	     		"Section list or the registration limits are invalid.",
	     	errMsgMallocStr[] = "Not enough memory available.";
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_AUTO|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecAuto FE 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx 0x%lx\n",
	   (unsigned long )rCtrl, (unsigned long )ppCtrl,
	   (unsigned long )secList, (unsigned long )cancelFlag,
	   (unsigned long )secFn, (unsigned long )secData,
	   (unsigned long )workFn, (unsigned long )workData,
	   (unsigned long )eMsg));
  if((rCtrl == NULL) || (ppCtrl == NULL) || (secList == NULL))
  {
    errFlag = REC_ERR_FUNC;
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(((item = RecSecFindItemIndex(secList, NULL, rCtrl->firstIdx,
     				    HGU_DLPLIST_DIR_TOTAIL)) == NULL) ||
       ((oSec0 = (RecSection *)HGUDlpListEntryGet(secList, item)) == NULL))
    {
      errFlag = REC_ERR_LIST;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(RecSecIsEmpty(oSec0))
    {
      if((oSec0 = RecSecNext(secList, item, &item, 1)) == NULL)
      {
        errFlag = REC_ERR_LIST;
      }
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if((oSec1 = RecSecNext(secList, item, &item, 1)) == NULL)
    {
      errFlag = REC_ERR_LIST;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if((oSec0->index != rCtrl->firstIdx) ||
       (oSec1->index < rCtrl->firstIdx) ||
       (oSec1->index > rCtrl->lastIdx))
    {
      errFlag = REC_ERR_LIST;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    if(((nSec0 = RecSecDup(oSec0)) == NULL) ||
       ((nSec1 = RecSecDup(oSec1)) == NULL))
    {
      errFlag = REC_ERR_MALLOC;
    }
  }
  if(errFlag == REC_ERR_NONE)
  {
    errFlag = RecFileSecObjRead(nSec0, eMsg);
  }
  while((errFlag == REC_ERR_NONE) && (*cancelFlag == 0) && (nSec1 != NULL))
  {
    if(secFn && nSec0)
    {
      (*secFn)(nSec0, secData);	/* Replaces oSec0 with copy of nSec0 in list */
    }
    errFlag = RecFileSecObjRead(nSec1, eMsg);
    if(errFlag == REC_ERR_NONE)
    {
      errFlag = RecRegisterPair(&(nSec1->transform),
				&(nSec1->correl), &(nSec1->iterations),
      				rCtrl, ppCtrl,
				nSec0->obj, nSec1->obj,
      				workFn, workData, eMsg);
    }
    if(errFlag == REC_ERR_NONE)
    {
      RecSecFree(nSec0);
      nSec0 = nSec1;
      nSec1 = NULL;
      if(nSec0->index < rCtrl->lastIdx)
      {
        if((oSec1 = RecSecNext(secList, item, &item, 1)) == NULL)
	{
	  errFlag = REC_ERR_LIST;
	}
	else if(oSec1->index <= rCtrl->lastIdx)
	{
	  if((nSec1 = RecSecDup(oSec1)) == NULL)
	  {
	    errFlag = REC_ERR_MALLOC;
	  }
	}
      }
    }
  }
  if((errFlag == REC_ERR_NONE) && secFn && nSec0)
  {
    (*secFn)(nSec0, secData);   /* Replaces oSec0 with copy of nSec0 in list */
  }
  if(nSec0)
  {
    RecSecFree(nSec0);
  }
  if(nSec1)
  {
    RecSecFree(nSec1);
  }
  if(*cancelFlag && (errFlag == REC_ERR_NONE))
  {
    errFlag = REC_ERR_CANCEL;
  }
  if((errFlag != REC_ERR_NONE) && (*eMsg == NULL))
  {
    switch(errFlag)
    {
      case REC_ERR_MALLOC:
        *eMsg = AlcStrDup(errMsgMallocStr);
        break;
      case REC_ERR_LIST:
        *eMsg = AlcStrDup(errMsgInvalidListStr);
        break;
      default:
        break;
    }
  }
  if(workFn && workData)
  {
    rState.approach = 0;
    rState.iteration = 0;
    rState.lastMethod = REC_MTHD_NONE;
    rState.transform = NULL;
    rState.correl = 0.0;
    rState.errFlag = errFlag;
    (*workFn)(&rState, workData);
  }
  REC_DBG((REC_DBG_AUTO|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecAuto FX %d\n",
	   errFlag));
  return(errFlag);
}
