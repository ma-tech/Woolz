#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgDebug.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Provides functions for debug output from the MRC
*	  	Human Genetics Unit numerical algorithm library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
************************************************************************/
#include <Alg.h>
#include <float.h>

AlgError 	AlgDbgWrite(char *, ...);

AlgDbgMask	algDbgMask = ALG_DBG_NONE;
void		*algDbgData = NULL;
AlgDbgFn	algDbgOutFn = AlgDbgWrite;

/************************************************************************
* Function:	AlgDbgWrite						*
* Returns:	AlgError:		Non zero if fails to output	*
*					given message.			*
* Purpose:	Writes out the given debug message to the debug file.	*
* Global refs:	void *algDbgData:	Used to hold debug output file.	*
* Parameters:	char *fmt:		Format for printing message.	*
*		...			Varargs function.		*
************************************************************************/
AlgError	AlgDbgWrite(char *fmt, ...)
{
  AlgError	errFlag = ALG_ERR_NONE;
  va_list	ap;

  if(algDbgData)
  {
    va_start(ap, fmt);
    (void )vfprintf((FILE *)algDbgData, fmt, ap);
    (void )fflush((FILE *)algDbgData);
    va_end (ap);
  }
  return(errFlag);
}
