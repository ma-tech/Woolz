#pragma ident   "MRC HGU $Id$"
/*!
* \file         AlcAlloc.c
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
* \brief        Provides functions for basic storage allocation.
*               In their most basic form are simple wrappers for the
*               ANSI functions malloc(3), calloc(3), realloc(3) and
*               free(3) but they may be used to encapsulate more
*               complex allocation such as for persistant storage.
* \todo		-
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/*!
* \ingroup	Alc
* \defgroup	AlcAlloc
* @{
*/

/*!
* \return				Allocated storage or NULL on error.
* \brief	Allocates space for the given number of elements with
*		each element set to zero. At it's most basic this
*		function is a wrapper for calloc(3).
* \param	elCount 		Number of elements.
* \param	elSz 			Size of an element.
*/
void		*AlcCalloc(int elCount, int elSz)
{
  void		*data = NULL;

  if((elCount > 0) && (elSz > 0))
  {
    data = calloc(elCount, elSz);
  }
  return(data);
}

/*!
* \return		 		Allocated storage or NULL on error.
* \brief	Allocates space for the given number of bytes with each
*		each element set an undefined value.
* \param	byteCount 	 	Number of bytes.
*/
void		*AlcMalloc(int byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = malloc(byteCount);
  }
  return(data);
}

/*!
* \return		 		Allocated storage or NULL on error.
* \brief	Re-allocates space for the given number of bytes with
*               the contents of given data being unchanged.
* \param	givenData 	 	Given storage.
* \param	byteCount 	 	Number of bytes required.
*/
void		*AlcRealloc(void *givenData, int byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = realloc(givenData, byteCount);
  }
  return(data);
}

/*!
* \return	<void>
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

/*!
* @}
*/
