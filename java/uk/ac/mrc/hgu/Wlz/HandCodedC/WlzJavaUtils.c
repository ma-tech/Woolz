#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzJavaUtils_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzJavaUtils.c
* \author       Bill Hill
* \date         January 1999
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
* \brief	Misc functions for the C side of Java Woolz.
* \ingroup	JWlz
*/
#include <WlzJava.h>

/*!
* \ingroup	JWlz
* \brief	Builds a fully qualified class name string which should
* 		be free'd (using free()) when no longer required.
* \param	name			Given class name.
*/
char		*WlzJavaBuildFQClassName(const char *name)
{
  int		fqLen;
  char		*fqName;
  const char 	pkg[] = "uk/ac/mrc/hgu/Wlz";

  fqLen = strlen(pkg) + strlen(name) + 2;
  if((fqName = (char *)malloc(sizeof(char) * fqLen)) != NULL)
  {
    (void )sprintf(fqName, "%s/%s", pkg, name);
  }
  return(fqName);
}

/*!
* \ingroup	JWlz
* \brief	Throws a WlzException with the Woolz error string as
* 		the message.
* \param	jEnv			Given JNI environment pointer.
* \param	errNum			Woolz error code.
*/
void		WlzJavaThrowWlzException(JNIEnv *jEnv, WlzErrorNum errNum)
{
  char		*eName;
  jclass	eCls;

  if((eName = WlzJavaBuildFQClassName("WlzException")) != NULL)
  {
    if((eCls = (*jEnv)->FindClass(jEnv, eName)) != NULL)
    {
      (void )((*jEnv)->ThrowNew(jEnv, eCls,
				WlzStringFromErrorNum(errNum, NULL)));
      (*jEnv)->DeleteLocalRef(jEnv, eCls);
    }
    free(eName);
  }
}
