#pragma ident   "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas					
* Title:        AlcAlloc.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:  	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Provides functions for basic storage allocation.
*		In their most basic form are simple wrappers for the
*		ANSI functions malloc(3), calloc(3), realloc(3) and
*		free(3) but they may be used to encapsulate more
*		complex allocation such as for persistant storage.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/************************************************************************
* Function:	AlcCalloc						*
* Returns:	void *:			Allocated storage or NULL on 	*
*					error.				*
* Purpose:	Allocates space for the given number of elements with	*
*		each element set to zero. At it's most basic this 	*
*		function is a wrapper for calloc(3).			*
* Global refs:	-							*
* Parameters:	int elCount:		Number of elements.		*
*		int elSz:		Size of an element.		*
************************************************************************/
void		*AlcCalloc(int elCount, int elSz)
{
  void		*data = NULL;

  if((elCount > 0) && (elSz > 0))
  {
    data = calloc(elCount, elSz);
  }
  return(data);
}

/************************************************************************
* Function:	AlcMalloc						*
* Returns:	void *:			Allocated storage or NULL on 	*
*					error.				*
* Purpose:	Allocates space for the given number of bytes with	*
*		each element set an undefined value.			*
* Global refs:	-							*
* Parameters:	int byteCount:		Number of bytes.		*
************************************************************************/
void		*AlcMalloc(int byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = malloc(byteCount);
  }
  return(data);
}

/************************************************************************
* Function:	AlcRealloc						*
* Returns:	void *:			Allocated storage or NULL on 	*
*					error.				*
* Purpose:	Re-allocates space for the given number of bytes with	*
*		the contents of given data being unchanged.		* 
* Global refs:	-							*
* Parameters:	void *givenData:	Given storage.			*
*		int byteCount:		Number of bytes required.	*
************************************************************************/
void		*AlcRealloc(void *givenData, int byteCount)
{
  void		*data = NULL;

  if(byteCount > 0)
  {
    data = realloc(givenData, byteCount);
  }
  return(data);
}

/************************************************************************
* Function:	AlcFree							*
* Returns:	void							*
* Purpose:	Free's the given storage.				*
* Global refs:	-							*
* Parameters:	void *data:		Given storage.			*
************************************************************************/
void		AlcFree(void *data)
{
  if(data)
  {
    free(data);
  }
}
