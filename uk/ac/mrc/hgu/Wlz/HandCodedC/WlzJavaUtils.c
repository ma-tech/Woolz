/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaUtils.c
* Date:         January 1999
* Purpose:      Misc functions for the C side of Java Woolz.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/
#include <WlzJava.h>

/************************************************************************
* Function:	WlzJavaBuildFQClassName					*
* Returns:	void							*
* Purpose:	Builds a fully qualified class name string which should	*
*		be free'd (using free()) when no longer required.	*
* Global refs:	-							*
* Parameters:	const char *name:	Given class name.		*
************************************************************************/
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

/************************************************************************
* Function:	WlzJavaThrowWlzException				*
* Returns:	void							*
* Purpose:	Throws a WlzException with the Woolz error string as	*
*		the message.						*
* Global refs:	int errno:		System error number.		*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
************************************************************************/
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
