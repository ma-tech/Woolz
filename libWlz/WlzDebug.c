#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzDebug.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Default debugging text and object output functions for the
* 		woolz library flag based debugging system.
* \ingroup	WlzDebug
* \todo         -
* \bug          None known.
*/
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

/*!
* \return	Woolz error code.
* \ingroup	WlzDebug
* \brief	Writes out the given debug message to the debug file.
* \param	fmt			Format for printing message.
*/
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

/*!
* \return	Woolz error code.
* \ingroup	WlzDebug
* \brief	Writes out the given debug Woolz object to the woolz debug
* 		file.
* \param	obj			Woolz debug object for output.
* \param	freeFlg			If non zero free the object after
* 					writing it.
*/
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
