#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcString.c
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
* \brief        Provides functions for string duplication.
* \todo		-
* \bug          None found.
*/
#include <stdio.h>
#include <string.h>
#include <Alc.h>

/*!
* \ingroup      Alc
* \defgroup	AlcString
* @{
*/

/*!
* \return		 		Duplicated string or NULL on error.
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

/*!
* @}
*/
