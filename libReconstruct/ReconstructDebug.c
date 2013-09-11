#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructDebug_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructDebug.c
* \author       Bill Hill
* \date         April 1999
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
* \brief	Provides functions for flag based debugging of the
*		Reconstruct library.
* \ingroup	Reconstruct
*/

#include <Reconstruct.h>
#include <unistd.h>


#if (defined DARWIN || defined _WIN32 || defined __MINGW32__ )
#define flockfile(F)
#define funlockfile(F)
#endif

RecDbgMask	recDbgMask = REC_DBG_NONE,
		recDbgWlzMask = REC_DBG_NONE;
void		*recDbgData = NULL,
		*recDbgWlzData = NULL;

RecDbgFn	recDbgOutFn = RecDbgWrite;
RecDbgWlzFn	recDbgOutWlzFn = RecDbgWlzWrite;

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes out the given debug message to the debug file.
* 		Global recDbgData is used to pass debug output file.
* \param	fmt			Format for printing message.
* \param	 ...			Varargs function.
*/
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

/*!
* \return	Error code.
* \ingroup	Reconstruct
* \brief	Writes out the given debug Woolz object to the woolz
*		debug file.
*		Global recDbgWlzData: is used to pass woolz debug output
*		file.
* \param	obj			Woolz debug object for output.
* \param	freeFlg			If non zero free the object after
*					writting it.
*/
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
#if (defined _WIN32 || defined __MINGW32__ )
    (void )Sleep(1);
#else
    (void )sleep(1);
#endif
  }
  return(errFlag);
}
