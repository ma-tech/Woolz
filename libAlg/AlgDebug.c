#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlg/AlgDebug.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Provides debug output.
* \ingroup	AlgDebug
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <float.h>

AlgError 	AlgDbgWrite(char *, ...);

AlgDbgMask	algDbgMask = ALG_DBG_NONE;
void		*algDbgData = NULL;
AlgDbgFn	algDbgOutFn = AlgDbgWrite;

/*!
* \return	Non zero if fails to output given message.
* \ingroup	AlgDebug
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
