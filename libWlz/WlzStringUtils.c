#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzStringUtils.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Miscellaneous string handling functions.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <Wlz.h>

/************************************************************************
* Function:	WlzStringMatchValue					*
* Returns:	int:			Zero if no match found.		*
* Purpose:	Given a destination pointer, a string, and a null	*
*		terminated list of string, enum pairs. The strings are	*
*		matched and the	first match is used to fill in the	*
*		destination with the matched enum's value.		*
* Global refs:	-							*
* Parameters:	int *datum:		Required datum to be filled in.	*
*		const char *targetStr:	String to match.		*
*		const char *testStr:	First string of the...		*
*		...			Varargs list of string (char *)	*
*					enum (int) pairs which is null	*
*					terminated.			*
************************************************************************/
/*VARARGS*/
int		WlzStringMatchValue(int *datum,
				    const char *targetStr, const char *testStr,
				    ...)
{
  va_list	ap;
  int		matchedFlag = 0,
		testVal;

  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
	  ("WlzStringMatchValue FE 0x%lx 0x%lx 0x%lx...\n",
	   (unsigned long )datum, (unsigned long )targetStr,
	   (unsigned long)testStr));
  va_start(ap, testStr);
  while(testStr)
  {
    testVal = va_arg(ap, int);
    WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
	    ("WlzStringMatchValue 01 0x%lx>%s< %d\n",
	     (unsigned long )testStr, (testStr)?(testStr):"(null)", testVal));
    if(!matchedFlag)
    {
      matchedFlag = !strncmp(targetStr, testStr, strlen(testStr));
      if(matchedFlag)
      {
	*datum = testVal;
      }
    }
    testStr = va_arg(ap, char *);
  }
  va_end(ap);
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
	  ("WlzStringMatchValue FX %d\n",
	   matchedFlag));
  return(matchedFlag);
}
