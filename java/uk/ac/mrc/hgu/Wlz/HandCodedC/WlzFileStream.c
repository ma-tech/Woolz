/************************************************************************
* Project:      Java Woolz
* Title:        WlzFileStream.c
* Date:         January 1999
* Purpose:      C binding for Java Woolz native file IO.
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

extern int	errno;

/************************************************************************
* Function:	ThrowFileIOException					*
* Returns:	void							*
* Purpose:	Throws an IOException with the system error string as	*
*		the message.						*
* Global refs:	int errno:		System error number.		*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
************************************************************************/
static void	ThrowFileIOException(JNIEnv *jEnv)
{
  jclass	eCls;

  if((eCls = (*jEnv)->FindClass(jEnv, "java/io/IOException")) != NULL)
  {
    (void )((*jEnv)->ThrowNew(jEnv, eCls, strerror(errno)));
    (*jEnv)->DeleteLocalRef(jEnv, eCls);
  }
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzClose		*
* Returns:	void							*
* Class:	WlzFileStream						*
* Method:	JWlzClose						*
* Purpose:	Native method implementation to close a WlzFileStream's	*
*		file.							*
* Global refs:	-							*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
*		jobject jObj:		Given java object.		*
*		jlong jValue:		Java long which holds the file	*
*					ptr.				*
************************************************************************/
JNIEXPORT void JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzClose(
				JNIEnv *jEnv, jobject jObj, jlong jValue)
{
  FILE		*fP = NULL;

  if(((fP = (FILE *)jValue) != NULL) &&
     (fP != stdin) && (fP != stdout) && (fP != stderr))
  {
    if(fclose(fP))
    {
      ThrowFileIOException(jEnv);
    }
  }
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzOpen		*
* Returns:	jlong:			Java long which holds the file	*
*					ptr.				*
* Class:	WlzFileStream						*
* Method:	JWlzClose						*
* Purpose:	Native method implementation to open the named file.	*
* Global refs:	-							*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
*		jobject jObj:		Given java object.		*
*		jstring jName:		Java file name string.		*
*		jstring jMode:		Java file mode string ("r" or	*
*					"w").				*
************************************************************************/
JNIEXPORT jlong JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzOpen(
				JNIEnv *jEnv, jobject jObj,
				jstring jName, jstring jMode)
{
  FILE		*fP = NULL;
  const char	*cName;
  jboolean      isCopyName;
  const char	*cMode;
  jboolean      isCopyMode;


/* #define JWLZ_DEBUG */
#ifdef JWLZ_DEBUG
  int hack = 1;
  (void )fprintf(stderr, "process id: %d", (int )getpid());
  while(hack)
  {
    sleep(2);
  }
#endif /* JWLZ_DEBUG */


  cMode = (*jEnv)->GetStringUTFChars(jEnv, jMode, &isCopyMode);
  cName = (*jEnv)->GetStringUTFChars(jEnv, jName, &isCopyName);

  if(strlen(cMode) == 0)
  {
    if(strcmp(cName, "stdin") == 0)
    {
      fP = stdin;
    }
    else if(strcmp(cName, "stdout") == 0)
    {
      fP = stdout;
    }
    else if(strcmp(cName, "stderr") == 0)
    {
      fP = stderr;
    }
  }
  else
  {
    if((fP = fopen(cName, cMode)) == NULL)
    {
      ThrowFileIOException(jEnv);
    }
  }

  if(isCopyName == JNI_TRUE) {
     (*jEnv)->ReleaseStringUTFChars(jEnv, jName, cName);
  }
  if(isCopyMode == JNI_TRUE) {
     (*jEnv)->ReleaseStringUTFChars(jEnv, jMode, cMode);
  }

  return((jlong )fP);
}
