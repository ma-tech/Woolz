#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructDebug.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides functions for flag based debugging for the
*		MRC Human Genetics Unit	reconstruction library.	
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>
#include <unistd.h>

RecDbgMask	recDbgMask = REC_DBG_NONE,
		recDbgWlzMask = REC_DBG_NONE;
void		*recDbgData = NULL,
		*recDbgWlzData = NULL;

RecDbgFn	recDbgOutFn = RecDbgWrite;
RecDbgWlzFn	recDbgOutWlzFn = RecDbgWlzWrite;

/************************************************************************
* Function:	RecDbgWrite					
* Returns:	RecError:		Non zero if fails to output
*					given message.		
* Purpose:	Writes out the given debug message to the debug file.
* Global refs:	void *recDbgData:	Used to pass debug output file.
* Parameters:	char *fmt:		Format for printing message.
*		...			Varargs function.	
************************************************************************/
RecError	RecDbgWrite(char *fmt, ...)
{
  RecError	errFlag = REC_ERR_NONE;
  va_list	ap;

  if(recDbgData)
  {
    flockfile((FILE *)recDbgData);
    va_start(ap, fmt);
    if((vfprintf((FILE *)recDbgData, fmt, ap) == EOF) ||
       (fflush((FILE *)recDbgData) == EOF))
    {
      errFlag = REC_ERR_WRITE;
    }
    va_end (ap);
    funlockfile((FILE *)recDbgData);
  }
  return(errFlag);
}

/************************************************************************
* Function:	RecDbgWlzWrite					
* Returns:	RecError:		Non zero if fails to output
*					given message.		
* Purpose:	Writes out the given debug Woolz object to the woolz
*		debug file.					
* Global refs:	FILE *recDbgWlzData:	Used to pass woolz debug output
*					file.			
* Parameters:	WlzObject *obj:		Woolz debug object for output.
*		int freeFlg:		If non zero free the object
*					after writting it.	
************************************************************************/
RecError	RecDbgWlzWrite(WlzObject *obj, int freeFlg)
{
  RecError	errFlag = REC_ERR_NONE;

  if(recDbgWlzData)
  {
    flockfile((FILE *)recDbgWlzData);
    errFlag = RecFileObjWlzWrite((FILE *)recDbgWlzData, obj);
    (void )fflush((FILE *)recDbgWlzData);
    if(freeFlg && obj)
    {
      WlzFreeObj(obj);
    }
    funlockfile((FILE *)recDbgWlzData);
    (void )sleep(1);
  }
  return(errFlag);
}
