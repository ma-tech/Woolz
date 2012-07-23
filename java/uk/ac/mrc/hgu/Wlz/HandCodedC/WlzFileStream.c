#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HandCodedC/WlzFileStream_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzFileStream.c
* \author       Bill Hill
* \date         July 1999
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
* \brief	C binding for Java Woolz native file IO.
* \ingroup	JWlz
*/
#include <WlzJava.h>

extern int	errno;

/*!
* \ingroup	JWlz
* \brief	Throws an IOException with the system error string as
* 		the message.
* \note		Global ref: errno.
* \param	jEnv			Given JNI environment pointer.
*/
static void	ThrowFileIOException(JNIEnv *jEnv)
{
  jclass	eCls;

  if((eCls = (*jEnv)->FindClass(jEnv, "java/io/IOException")) != NULL)
  {
    (void )((*jEnv)->ThrowNew(jEnv, eCls, strerror(errno)));
    (*jEnv)->DeleteLocalRef(jEnv, eCls);
  }
}

/*!
* \ingroup	JWlz
* \brief	Native method implementation to close a WlzFileStream's file.
* \param	jEnv			Given JNI environment pointer.
* \param	jObj			Given java object.
* \param	jValue			Java long which holds the file pointer.
*/
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

/*!
* \return	Java long which holds the file pointer.
* \ingroup	JWlz
* \brief	Native method implementation to open the named file.
* \param	jEnv			Given JNI environment pointer.
* \param	jObj			Given java object.
* \param	jName			Java file name string.
* \param	jMode			Java file mode string ("r" or "w").
*/
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
    sleep(1);
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
