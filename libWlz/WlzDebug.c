#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDebug_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzDebug.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Default debugging text and object output functions for
* 		the Woolz library flag based debugging system.
* \ingroup	WlzDebug
*/

#include <stdio.h>
#include <stdarg.h>
#include <Wlz.h>

WlzErrorNum 	WlzDbgWrite(char *fmt, ...),
		WlzDbgObjWrite(WlzObject *obj, int freeFlg);

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

