#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzDebug.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Provides default debugging text and object output
*		functions for the woolz library flag based debugging
*		system.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <Wlz.h>

WlzErrorNum 	WlzDbgWrite(char *, ...),
		WlzDbgObjWrite(WlzObject *, int);

WlzDbgMask	wlzDbgMask = WLZ_DBG_NONE,
		wlzDbgObjMask = WLZ_DBG_NONE;
void		*wlzDbgData = NULL,
		*wlzDbgObjData = NULL;

WlzDbgFn	wlzDbgOutFn = WlzDbgWrite;
WlzDbgObjFn	wlzDbgOutOutFn = WlzDbgObjWrite;

/************************************************************************
* Function:	WlzDbgWrite						*
* Returns:	WlzErrorNum:		Non zero if fails to output	*
*					given message.			*
* Purpose:	Writes out the given debug message to the debug file.	*
* Global refs:	void *wlzDbgData:	Used to hold debug output file.	*
* Parameters:	char *fmt:		Format for printing message.	*
*		...			Varargs function.		*
************************************************************************/
WlzErrorNum	WlzDbgWrite(char *fmt, ...)
{
  WlzErrorNum	errFlag = WLZ_ERR_NONE;
  va_list	ap;

  if(wlzDbgData)
  {
    va_start(ap, fmt);
    if((vfprintf((FILE *)wlzDbgData, fmt, ap) == EOF) ||
       (fflush((FILE *)wlzDbgData) == EOF))
    {
      errFlag = WLZ_ERR_UNSPECIFIED;
    }
    va_end (ap);
  }
  return(errFlag);
}

/************************************************************************
* Function:	WlzDbgObjWrite						*
* Returns:	WlzErrorNum:		Non zero if fails to output	*
*					given message.			*
* Purpose:	Writes out the given debug Woolz object to the woolz	*
*		debug file.						*
* Global refs:	void *wlzDbgObjData:	Used to hold woolz debug output	*
*					file.				*
* Parameters:	WlzObject *obj:		Woolz debug object for output.	*
*		int freeFlg:		If non zero free the object	*
*					after writting it.		*
************************************************************************/
WlzErrorNum	WlzDbgObjWrite(WlzObject *obj, int freeFlg)
{
  WlzErrorNum	errFlag = WLZ_ERR_NONE;

  if(wlzDbgObjData)
  {
    if(WlzWriteObj((FILE *)wlzDbgObjData, obj) != 0)
    {
      errFlag = WLZ_ERR_UNSPECIFIED;
    }
    (void )fflush((FILE *)wlzDbgObjData);
    if(freeFlg && obj)
    {
      WlzFreeObj(obj);
    }
  }
  return(errFlag);
}
