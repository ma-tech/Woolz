#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructPreProc.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides functions for pre-processing woolz objects
*		to aid registration through cross correlation for
*		the MRC Human Genetics Unit reconstruction library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>
#include <string.h>

/************************************************************************
* Function:	RecPreProcObj					
* Returns:	WlzObject *:		Pre-processed woolz object with
*					link count incremented or
*					NULL on error.		
* Purpose:	Pre-process the given object according to the given
*		pre-processing mask.				
* Note:		The link count of the returned (pre-processed) object
*		IS incremented by this function.		
* Global refs:	-						
* Parameters:	WlzObject *obj:		Given object.		
*		RecPPControl *ppCtrl:	Pre-processing control data
*					structure.		
*		RecError *dstErr:	Destination ptr for error code.
************************************************************************/
WlzObject	*RecPreProcObj(WlzObject *obj, RecPPControl *ppCtrl,
			       RecError *dstErr)
{
  int		tI0;
  WlzObject	*ppObj1 = NULL,
		*ppObj2 = NULL,
		*ppObj3 = NULL,
		*rtnObj = NULL;
  double	tD0;
  WlzPixelV	tV0,
  		tV1;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;
  WlzWindowFnType winFn;
  WlzIVertex2	samFac,
  		winRad,
		winOrg;
  RecError	errFlag = REC_ERR_NONE;

  REC_DBG((REC_DBG_PPROC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecPreProcObj FE 0x%lx 0x%lx\n",
	   (unsigned long )obj, (unsigned long )ppCtrl));
  if((obj == NULL) || (obj->type != WLZ_2D_DOMAINOBJ) ||
     (obj->domain.core == NULL) || (obj->values.core == NULL) ||
     (ppCtrl == NULL))
  {
    errFlag = REC_ERR_WLZ;
  }
  else
  {
    REC_DBG((REC_DBG_PPROC|REC_DBG_LVL_1),
    	    ("RecPreProcObj 01 %d 0x%lx %d\n",
	     ppCtrl->method, (unsigned long )obj, obj->linkcount));
    if(ppCtrl->method & REC_PP_SAMPLE)
    {
      samFac.vtX = ppCtrl->sample.factor;
      samFac.vtY = ppCtrl->sample.factor;
    }
    else
    {
      samFac.vtX = 1;
      samFac.vtY = 1;
    }
    if(((ppObj1 = WlzAssignObject(
                  WlzSampleObj(obj, samFac, ppCtrl->sample.function,
    			       &wlzErr), NULL)) == NULL) ||
       (wlzErr != WLZ_ERR_NONE))
    {
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_INVERT))
  {
    tV0.type = WLZ_GREY_INT;
    tV0.v.inv = 0;
    tV1.type = WLZ_GREY_INT;
    tV1.v.inv = 255;
    if((wlzErr = WlzGreyInvertMinMax(ppObj1, tV0, tV1)) != WLZ_ERR_NONE)
    {
      errFlag = REC_ERR_WLZ;
    }
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_LAPLAC))
  {
    if(((ppObj2 = WlzAssignObject(WlzLaplacian(ppObj1, 3, 1, 1,
    				               &wlzErr), NULL)) != NULL) &&
       (wlzErr == WLZ_ERR_NONE))
    {
      WlzFreeObj(ppObj1);
      ppObj1 = ppObj2;
      ppObj2 = NULL;
    }
    errFlag = RecErrorFromWlz(wlzErr);
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_SOBEL))
  {
    if(((ppObj2 = WlzAssignObject(WlzSobel(ppObj1, 1, 1,
    				           &wlzErr), NULL)) != NULL) &&
       (wlzErr == WLZ_ERR_NONE))
    {
      WlzFreeObj(ppObj1);
      ppObj1 = ppObj2;
      ppObj2 = NULL;
    }
    errFlag = RecErrorFromWlz(wlzErr);
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_THRESH))
  {
    if(((ppObj3 = WlzHistogramObj(ppObj1, 256, 0.0, 1.0, &wlzErr)) != NULL) &&
       (wlzErr == WLZ_ERR_NONE))
    {
      if((wlzErr = WlzCompThreshold(&(tV0.v.dbv), ppObj3,
      				    WLZ_COMPTHRESH_GRADIENT,
				    0.05)) ==  WLZ_ERR_NONE)
      {
        tV0.type = WLZ_GREY_DOUBLE;
      }
    }
    errFlag = RecErrorFromWlz(wlzErr);
    if(ppObj3)
    {
      WlzFreeObj(ppObj3);
      ppObj3 = NULL;
    }
    if(errFlag == REC_ERR_NONE)
    {
      if(((ppObj2 = WlzAssignObject(WlzThreshold(ppObj1, tV0,
      						 WLZ_THRESH_HIGH,
						 &wlzErr), NULL)) != NULL) &&
         (wlzErr == WLZ_ERR_NONE))
      {
        WlzFreeObj(ppObj1);
	ppObj1 = ppObj2;
	ppObj2 = NULL;
      }
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_ERODE))
  {
    if((tI0 = ppCtrl->erode) > 0)
    {
      while((tI0-- > 0) && ppObj1 && (wlzErr == WLZ_ERR_NONE))
      {
        if(((ppObj2 = WlzAssignObject(WlzErosion(ppObj1, WLZ_4_CONNECTED,
				                 &wlzErr), NULL)) != NULL) &&
	   (wlzErr == WLZ_ERR_NONE))
	{
	  if(ppObj2->type == WLZ_2D_DOMAINOBJ)
	  {
	    ppObj3 = WlzAssignObject(
	    	     WlzMakeMain(WLZ_2D_DOMAINOBJ, ppObj2->domain,
		     		 ppObj1->values, NULL, NULL, &wlzErr), NULL);
	    WlzFreeObj(ppObj1);
	    WlzFreeObj(ppObj2);
	    ppObj1 = ppObj3;
	    ppObj2 = NULL;
	    ppObj3 = NULL;
	  }
	  else
	  {
	    ppObj1 = WlzMakeEmpty(&wlzErr);
	    WlzFreeObj(ppObj2);
	    ppObj2 = NULL;
	  }
	}
      }
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_BACKGROUND))
  {
    tV0.type = WLZ_GREY_INT;
    tV0.v.ubv = 0;
    wlzErr = WlzSetBackground(ppObj1, tV0);
    errFlag = RecErrorFromWlz(wlzErr);
  }
  if((errFlag == REC_ERR_NONE) && (ppCtrl->method & REC_PP_WINDOW))
  {
    if(ppObj1->type == WLZ_2D_DOMAINOBJ)
    {
      winFn = ppCtrl->window.function;
      winRad.vtX = ((ppObj1->domain.i->lastkl - ppObj1->domain.i->kol1) *
		    ppCtrl->window.size.vtX) / 200;
      winRad.vtY = ((ppObj1->domain.i->lastln - ppObj1->domain.i->line1) *
		    ppCtrl->window.size.vtY) / 200;
      winOrg.vtX = ((ppObj1->domain.i->lastkl + ppObj1->domain.i->kol1) / 2) +
		    ppCtrl->window.offset.vtX;
      winOrg.vtY = ((ppObj1->domain.i->lastln + ppObj1->domain.i->line1) / 2) +
		    ppCtrl->window.offset.vtY;
      if(((ppObj2 = WlzAssignObject(WlzWindow(ppObj1, winFn, winOrg, winRad,
      					      &wlzErr), NULL)) != NULL) &&
         (wlzErr == WLZ_ERR_NONE))
      {
	WlzFreeObj(ppObj1);
	ppObj1 = ppObj2;
	ppObj2 = NULL;
      }
      errFlag = RecErrorFromWlz(wlzErr);
    }
  }
  if(errFlag != REC_ERR_NONE)
  {
    WlzFreeObj(ppObj1);
    ppObj1 = NULL;
    WlzFreeObj(ppObj2);
  }
  else
  {
    rtnObj = ppObj1;
    REC_DBG((REC_DBG_PPROC|REC_DBG_LVL_1),
    	    ("RecPreProcObj 02 %d\n",
	     rtnObj->linkcount));
  }
  if(dstErr)
  {
    *dstErr = errFlag;
  }
  REC_DBG((REC_DBG_PPROC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecPreProcObj FX 0x%lx\n",
	   (unsigned long )rtnObj));
  return(rtnObj);
}
