#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcString.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Provides functions for string duplication.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <string.h>
#include <Alc.h>

/************************************************************************
* Function:	AlcStrDup						*
* Returns:	char *:			Duplicated string or NULL on	*
*					error.				*
* Purpose:	Allocates space for and duplicates the given NULL	*
*		terminated character string.				*
* Global refs:	-							*
* Parameters:	const char *srcStr:	Given string.			*
************************************************************************/
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

