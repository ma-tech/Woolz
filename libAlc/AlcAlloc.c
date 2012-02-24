#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcAlloc_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcAlloc.c
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
* \brief        Provides functions for basic storage allocation.
*               In their most basic form are simple wrappers for the
*               ANSI functions malloc(3), calloc(3), realloc(3) and
*               free(3) but they may be used to encapsulate more
*               complex allocation such as for persistant storage.
* \ingroup	AlcAlloc
*/

#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/*!
* \return	Allocated storage or NULL on error.
* \ingroup	AlcAlloc
* \brief	Allocates space for the given number of elements with
*		each element set to zero. At it's most basic this
*		function is a wrapper for calloc(3).
* \param	elCount 		Number of elements.
* \param	elSz 			Size of an element.
*/
void		*AlcCalloc(size_t elCount, size_t elSz)
{
  void		*data = NULL;

  if((elCount > 0) && (elSz > 0))
  {
    data = calloc(elCount, elSz);
  }
  return(data);
}

/*!
* \return	Allocated storage or NULL on error.
* \ingroup	AlcAlloc
* \brief	Allocates space for the given number of bytes with each
*		each element set an undefined value.
* \param	byteCount 	 	Number of bytes.
*/
void		*AlcMalloc(size_t byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = malloc(byteCount);
  }
  return(data);
}

/*!
* \return	Allocated storage or NULL on error.
* \ingroup	AlcAlloc
* \brief	Re-allocates space for the given number of bytes with
*               the contents of given data being unchanged.
* \param	givenData 	 	Given storage.
* \param	byteCount 	 	Number of bytes required.
*/
void		*AlcRealloc(void *givenData, size_t byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = realloc(givenData, byteCount);
  }
  return(data);
}

/*!
* \return	void
* \ingroup	AlcAlloc
* \brief	Free's the given storage.
* \param	data 	 		Given storage.
*/
void		AlcFree(void *data)
{
  if(data)
  {
    free(data);
  }
}
