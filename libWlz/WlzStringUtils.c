#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzStringUtils.c
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Miscellaneous string handling functions.
* \ingroup	WlzStrings
* \todo         -
* \bug          None known.
*/
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <Wlz.h>

/*!
* \return	Zero if no match found.
* \ingroup	WlzStrings
* \brief	Given a destination pointer, a string, and a null terminated
*		list of string, enum pairs. The strings are matched and the
*		first match is used to fill in the destination with the matched
*		enum's value.
* \param	datum			Pointer to datum for return.
* \param	targetStr		String to match.
* \param	testStr			First string of the list.
* \param	 ...			NULL terminated varargs list of
*					strings (char *) and enum (int)
*					pairs.
*/
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
  if(targetStr)
  {
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
  }
  WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
	  ("WlzStringMatchValue FX %d\n",
	   matchedFlag));
  return(matchedFlag);
}

/*!
* \return	Zero if no match found.
* \ingroup	WlzStrings
* \brief	Given a destination pointer, a value, and a null terminated
*		list of string, enum pairs. The values are matched and the
*		first match is used to fill in the destination with the matched
*		enum's string.
* \param	datum			Pointer to string for return.
* \param	targetVal		Value to match.
* \param	testStr			First string of the list.
* \param	 ...			NULL terminated varargs list of
*					strings (char *) and enum (int)
*					pairs.
*/
int		WlzValueMatchString(char **datum,
				    int targetVal, const char *testStr,
				    ...)
{
  va_list	ap;
  int		matchedFlag = 0,
		testVal;

  va_start(ap, testStr);
  while(testStr)
  {
    testVal = va_arg(ap, int);
    if(!matchedFlag)
    {
      if( targetVal == testVal ){
	matchedFlag = 1;
	*datum = (char *)testStr;
      }
    }
    testStr = va_arg(ap, char *);
  }
  va_end(ap);
  return(matchedFlag);
}
