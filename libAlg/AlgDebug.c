#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgDebug.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides debug output.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgDebug
* @{
*/

#include <Alg.h>
#include <float.h>

AlgError 	AlgDbgWrite(char *, ...);

AlgDbgMask	algDbgMask = ALG_DBG_NONE;
void		*algDbgData = NULL;
AlgDbgFn	algDbgOutFn = AlgDbgWrite;

/*!
* \return				Non zero if fails to output
*					given message.
* \brief	Writes out the given debug message to the debug file.
* \note		Uses global void *algDbgData to hold debug output file.
* \param	fmt:			Format for printing message.
* \param	...			Varargs function.
*/
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
