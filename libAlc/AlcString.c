#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcString_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlc/AlcString.c
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
* \brief        Provides functions for string duplication.
* \ingroup	AlcString
* \todo		-
* \bug          None found.
*/
#include <stdio.h>
#include <string.h>
#include <Alc.h>

/*!
* \return	Duplicated string or NULL on error.
* \ingroup	AlcString
* \brief	Allocates space for and duplicates the given NULL
*		terminated character string.
* \param	srcStr			 Given string.
*/
char		*AlcStrDup(const char *srcStr)
{
  char		*dstStr = NULL;

  if(srcStr &&
     ((dstStr = AlcMalloc((strlen(srcStr) + 1) * sizeof(char))) != NULL))
  {
    (void )strcpy(dstStr, srcStr);
  }
  return(dstStr);
}
