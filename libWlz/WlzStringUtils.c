#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzStringUtils_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzStringUtils.c
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
* \brief	Miscellaneous string handling functions.
* \ingroup	WlzStrings
*/

#include <ctype.h>
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
	  ("WlzStringMatchValue FE %p %p %p...\n",
	   datum, targetStr, testStr));
  if(targetStr)
  {
    va_start(ap, testStr);
    while(testStr)
    {
      testVal = va_arg(ap, int);
      WLZ_DBG((WLZ_DBG_LVL_FN|WLZ_DBG_LVL_3),
	      ("WlzStringMatchValue 01 %p >%s< %d\n",
	       testStr, (testStr)? (testStr): "(null)", testVal));
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

/*!
* \return	The given string.
* \ingroup	WlzStrings
* \brief	Removes all white space characters (as determined by
* 		isspace(3)) from the given string.
* \param	str			Given string.
*/
char		*WlzStringWhiteSpSkip(char *str)
{
  char		*s0,
  		*s1;

  if(str)
  {
    s0 = s1 = str;
    while(*s1)
    {
      *s0 = *s1++;
      s0 += (isspace(*s0) == 0);
    }
    *s0 = '\0';
  }
  return(str);
}

/*!
* \return	The given string.
* \ingroup	WlzStrings
* \brief	Converts the string to all upper case using toupper(3).
* \param	str			Given string.
*/
char		*WlzStringToUpper(char *str)
{
  char		*s0;

  if(str)
  {
    s0 = str;
    while(*s0)
    {
      *s0 = toupper(*s0);
      ++s0;
    }
  }
  return(str);
}

/*!
* \return	The given string.
* \ingroup	WlzStrings
* \brief	Converts the string to all lower case using tolower(3).
* \param	str			Given string.
*/
char		*WlzStringToLower(char *str)
{
  char		*s0;

  if(str)
  {
    s0 = str;
    while(*s0)
    {
      *s0 = tolower(*s0);
      ++s0;
    }
  }
  return(str);
}
