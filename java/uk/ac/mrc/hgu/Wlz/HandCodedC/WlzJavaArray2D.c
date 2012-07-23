#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HandCodedC/WlzJavaArray2D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzJavaArray2D.c
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
* \brief	2D array functions for the C side of Java Woolz.
* \ingroup	JWlz
*/
#include <WlzJava.h>

/*!
* \return	The value held in the given Java object.
* \ingroup	JWlz
* \brief	Returns a 2D array built from the given 2D java array.
* \param	jEnv			Given JNI environment ptr.
* \param	cObjName		The Java woolz C object class string.	
* \param	jObjName		The Java woolz Java object class
* 					string.
* \param	jniObjName		The Java woolz JNI object class string.
* \param	idrCnt			Indirection count (ie 1 for *, 2 for
* 					**, ...).
* \param	pKey			Parameter key.
* \param	jWArray			The Java Woolz array.
* \param	wArraySz		The number of elements in the Java
* 					Woolz array.
* \param	isCpy			Destination pointer for JNI copy flag.
*/
jlong 		WlzJavaArray2DGet(JNIEnv *jEnv,
				       char *cObjName,
				       char *jObjName,
				       char *jniObjName,
				       int idrCnt, int pKey,
				       jarray jWArray,
				       WlzIVertex2 wArraySz,
				       jboolean *isCpy)
{
  int		idY,
  		ok = 0;
  void		*bufJ;
  jobject	aryJ1D;
  void		**aryW2D = NULL;
  jlong		rtnVal = 0;

  if(jWArray && (wArraySz.vtX > 0) && (wArraySz.vtY > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_BYTE_ARY2:
	ok = AlcChar2Malloc((char ***)&aryW2D,
			    wArraySz.vtY, wArraySz.vtX) == ALC_ER_NONE;
	if(ok)
	{
	  idY = 0;
	  while(ok && (idY < wArraySz.vtY))
	  {
	    ok = ((aryJ1D = (*jEnv)->GetObjectArrayElement(jEnv,
	    					(jobjectArray )jWArray,
						idY)) != NULL) &&
	         ((bufJ = (void *)(*jEnv)->GetByteArrayElements(jEnv,
				    		(jbyteArray )aryJ1D,
						isCpy)) == NULL);
	    if(ok)
	    {
	      (void )memcpy(*((WlzUByte **)aryW2D + idY), bufJ,
	      		    wArraySz.vtX * sizeof(jbyte));
	      if(*isCpy)
	      {
	        (*jEnv)->ReleaseByteArrayElements(jEnv, (jbyteArray )aryJ1D,
						  (jbyte *)bufJ, 0);
	        *isCpy = JNI_FALSE;
	      }
	    }
	    ++idY;
	  }
	}
	break;
      case WLZ_JPM_KEY_SHORT_ARY2:
	ok = AlcShort2Malloc((short ***)&aryW2D,
			     wArraySz.vtY, wArraySz.vtX) == ALC_ER_NONE;
	if(ok)
	{
	  idY = 0;
	  while(ok && (idY < wArraySz.vtY))
	  {
	    ok = ((aryJ1D = (*jEnv)->GetObjectArrayElement(jEnv,
	    					(jobjectArray )jWArray,
						idY)) != NULL) &&
	         ((bufJ = (void *)(*jEnv)->GetShortArrayElements(jEnv,
		 				(jshortArray )aryJ1D,
						isCpy)) != NULL);
	    if(ok)
	    {
	      (void )memcpy(*((WlzUByte **)aryW2D + idY), bufJ,
	      		    wArraySz.vtX * sizeof(short));
	      if(*isCpy)
	      {
	        (*jEnv)->ReleaseShortArrayElements(jEnv,
						(jshortArray )aryJ1D,
						(jshort *)bufJ, 0);
	        *isCpy = JNI_FALSE;
	      }
	    }
	    ++idY;
	  }
	}
	break;
      case WLZ_JPM_KEY_INT_ARY2:
	ok = AlcInt2Malloc((int ***)&aryW2D,
			   wArraySz.vtY, wArraySz.vtX) == ALC_ER_NONE;
	if(ok)
	{
	  idY = 0;
	  while(ok && (idY < wArraySz.vtY))
	  {
	    ok = ((aryJ1D = (*jEnv)->GetObjectArrayElement(jEnv,
	    					(jobjectArray )jWArray,
						idY)) != NULL) &&
	         ((bufJ = (void *)(*jEnv)->GetIntArrayElements(jEnv,
						(jintArray )aryJ1D,
						isCpy)) != NULL);
	    if(ok)
	    {
	      (void )memcpy(*((WlzUByte **)aryW2D + idY), bufJ,
	      		    wArraySz.vtX * sizeof(int));
	      if(*isCpy)
	      {
	        (*jEnv)->ReleaseIntArrayElements(jEnv,
						(jintArray )aryJ1D,
						(jint *)bufJ, 0);
	        *isCpy = JNI_FALSE;
	      }
	    }
	    ++idY;
	  }
	}
	break;
      case WLZ_JPM_KEY_FLOAT_ARY2:
	ok = AlcFloat2Malloc((float ***)&aryW2D,
			     wArraySz.vtY, wArraySz.vtX) == ALC_ER_NONE;
        if(ok)
	{
	  idY = 0;
	  while(ok && (idY < wArraySz.vtY))
	  {
	    ok = ((aryJ1D = (*jEnv)->GetObjectArrayElement(jEnv,
	    					(jobjectArray )jWArray,
						idY)) != NULL) &&
	       ((bufJ = (void *)(*jEnv)->GetFloatArrayElements(jEnv,
						(jfloatArray )aryJ1D,
						 isCpy)) != NULL);
	    if(ok)
	    {
	      (void )memcpy(*((WlzUByte **)aryW2D + idY), bufJ,
	      		    wArraySz.vtX * sizeof(float));
	      if(*isCpy)
	      {
	        (*jEnv)->ReleaseFloatArrayElements(jEnv,
						(jfloatArray )aryJ1D,
						(jfloat *)bufJ, 0);
	        *isCpy = JNI_FALSE;
	      }
	    }
	    ++idY;
	  }
	}
	break;
      case WLZ_JPM_KEY_DOUBLE_ARY2:
	ok = AlcDouble2Malloc((double ***)&aryW2D,
			      wArraySz.vtY, wArraySz.vtX) == ALC_ER_NONE;
	if(ok)
	{
	  idY = 0;
	  while(ok && (idY < wArraySz.vtY))
	  {
	    ok = ((aryJ1D = (*jEnv)->GetObjectArrayElement(jEnv,
	    					(jobjectArray )jWArray,
						idY)) != NULL) &&
	       ((bufJ = (void *)(*jEnv)->GetDoubleArrayElements(jEnv,
						(jdoubleArray )aryJ1D,
						isCpy)) != NULL);
	    if(ok)
	    {
	      (void )memcpy(*((WlzUByte **)aryW2D + idY), bufJ,
	      		    wArraySz.vtX * sizeof(double));
	      if(*isCpy)
	      {
	        (*jEnv)->ReleaseDoubleArrayElements(jEnv,
						(jdoubleArray )aryJ1D,
						(jdouble *)bufJ, 0);
	        *isCpy = JNI_FALSE;
	      }
	    }
	    ++idY;
	  }
	}
	break;
      default:
	break;
    }
    if(ok)
    {
      rtnVal = (jlong )aryW2D;
      *isCpy = JNI_TRUE;
    }
    else if(aryW2D)
    {
      Alc2Free(aryW2D);
    }
  }
  return(rtnVal);
}

/*!
* \ingroup	JWlz
* \brief	Sets the given value in a Java object for return to the Java
* 		side of th JNI.
* \param	jEnv			Given JNI environment pointer.
* \param	dstJObj			Destination Java array.
* \param	cObjName		The type of the value to set.
* \param	jObjName		The Java woolz Java object class
* 					string.
* \param	jniObjName		The Java woolz JNI object class string.
* \param	idrCnt			Indirection count (ie 1 for *, 2 for
* 					**, ...).
* \param	pKey			Parameter key.
* \param	aVal			C array to set value from.
* \param	aSz			Size of the 2D array.
*/
void		WlzJavaArray2DSet(JNIEnv *jEnv, jobjectArray dstJObj,
				char *cObjName,
				char *jObjName,
				char *jniObjName,
				int idrCnt, int pKey,
				void *aVal,
				WlzIVertex2 aSz)
{
  jobject	newJObj;

  if(aVal && (aSz.vtX > 0) && (aSz.vtY > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_INT_PTR1_ARY2:
	newJObj = WlzJavaArray2DWrap(jEnv, "int", "jint", "jintArray",
					  2, WLZ_JPM_KEY_INT_ARY2,
					  aVal, aSz);
	(*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
	break;
      case WLZ_JPM_KEY_SHORT_PTR1_ARY2:
	newJObj = WlzJavaArray2DWrap(jEnv, "short", "jshort", "jshortArray",
					  2, WLZ_JPM_KEY_SHORT_ARY2,
					  aVal, aSz);
	(*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
	break;
      case WLZ_JPM_KEY_BYTE_PTR1_ARY2:
	newJObj = WlzJavaArray2DWrap(jEnv, "char", "jbyte", "jbyteArray",
					  2, WLZ_JPM_KEY_BYTE_ARY2,
					  aVal, aSz);
	(*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
	break;
      case WLZ_JPM_KEY_FLOAT_PTR1_ARY2:
	newJObj = WlzJavaArray2DWrap(jEnv, "float", "jfloat", "jfloatArray",
					  2, WLZ_JPM_KEY_FLOAT_ARY2,
					  aVal, aSz);
	(*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
	break;
      case WLZ_JPM_KEY_DOUBLE_PTR1_ARY2:
	newJObj = WlzJavaArray2DWrap(jEnv, "double", "jdouble", "jdoubleArray",
					  2, WLZ_JPM_KEY_DOUBLE_ARY2,
					  aVal, aSz);
	(*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
	break;
      default:
	break;
    }
  }
}

/*!
* \return	New Java object.
* \ingroup	JWlz
* \brief	Wraps up the given array as a new Java array.
* \param	jEnv			Given JNI environment pointer.
* \param	cObjName		The Java woolz object class string.
* \param	jObjName		Java woolz Java object class string.
* \param	jniObjName		The Java woolz JNI object class string.
* \param	idrCnt			Indirection count (ie 1 for *, 2 for
* 					**, ...).
* \param	pKey			Parameter key.
* \param	aVal			C array to set values from.
* \param	aSz			Size of the 2D array.
*/
jobject		WlzJavaArray2DWrap(JNIEnv *jEnv,
				   char *cObjName,
				   char *jObjName,
				   char *jniObjName,
				   int idrCnt, int pKey,
				   void *aVal,
				   WlzIVertex2 aSz)
{
  int		idY;
  jclass	elmClass;
  jobject	jArray2D = NULL,
		jArray1D = NULL,
  		rtnJObj = NULL;

#ifdef JWLZ_DEBUG
  int hack = 1;
  while(hack)
  {
    sleep(1);
  }
#endif /* JWLZ_DEBUG */
  if(aVal && (aSz.vtX > 0) && (aSz.vtY > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_INT_ARY2:
        jArray1D = (jobject )((*jEnv)->NewIntArray(jEnv, aSz.vtX));
	if(jArray1D)
	{
	  elmClass = (*jEnv)->GetObjectClass(jEnv, jArray1D);
	  jArray2D = (*jEnv)->NewObjectArray(jEnv, aSz.vtY, elmClass,
	  				     (jobject )NULL);
	}
	if(jArray2D)
	{
	  idY = 0;
	  while(jArray1D && (idY < aSz.vtY))
	  {
	    (*jEnv)->SetIntArrayRegion(jEnv, (jintArray )jArray1D,
				      0, aSz.vtX,
				      (jint *)*((int **)aVal + idY));
	    (*jEnv)->SetObjectArrayElement(jEnv, jArray2D, idY, jArray1D);
	    jArray1D = (jobject )((*jEnv)->NewIntArray(jEnv, aSz.vtX));
	    ++idY;
	  }
	}
	if(jArray2D)
	{
	  rtnJObj = jArray2D;
	}
	break;
      case WLZ_JPM_KEY_SHORT_ARY2:
        jArray1D = (jobject )((*jEnv)->NewShortArray(jEnv, aSz.vtX));
	if(jArray1D)
	{
	  elmClass = (*jEnv)->GetObjectClass(jEnv, jArray1D);
	  jArray2D = (*jEnv)->NewObjectArray(jEnv, aSz.vtY, elmClass,
	  				     (jobject )NULL);
	}
	if(jArray2D)
	{
	  idY = 0;
	  while(jArray1D && (idY < aSz.vtY))
	  {
	    (*jEnv)->SetShortArrayRegion(jEnv, (jshortArray )jArray1D,
				    0, aSz.vtX,
				    (jshort *)*((short **)aVal + idY));
	    (*jEnv)->SetObjectArrayElement(jEnv, jArray2D, idY, jArray1D);
	    jArray1D = (jobject )((*jEnv)->NewShortArray(jEnv, aSz.vtX));
	    ++idY;
	  }
	}
	if(jArray2D)
	{
	  rtnJObj = jArray2D;
	}
	break;
      case WLZ_JPM_KEY_BYTE_ARY2:
        jArray1D = (jobject )((*jEnv)->NewByteArray(jEnv, aSz.vtX));
	if(jArray1D)
	{
	  elmClass = (*jEnv)->GetObjectClass(jEnv, jArray1D);
	  jArray2D = (*jEnv)->NewObjectArray(jEnv, aSz.vtY, elmClass,
	  				     (jobject )NULL);
	}
	if(jArray2D)
	{
	  idY = 0;
	  while(jArray1D && (idY < aSz.vtY))
	  {
	    (*jEnv)->SetByteArrayRegion(jEnv, (jbyteArray )jArray1D,
	    			     0, aSz.vtX,
				     (jbyte *)*((char **)aVal + idY));
	    (*jEnv)->SetObjectArrayElement(jEnv, jArray2D, idY, jArray1D);
	    jArray1D = (jobject )((*jEnv)->NewByteArray(jEnv, aSz.vtX));
	    ++idY;
	  }
	}
	if(jArray2D)
	{
	  rtnJObj = jArray2D;
	}
	break;
      case WLZ_JPM_KEY_FLOAT_ARY2:
        jArray1D = (jobject )((*jEnv)->NewFloatArray(jEnv, aSz.vtX));
	if(jArray1D)
	{
	  elmClass = (*jEnv)->GetObjectClass(jEnv, jArray1D);
	  jArray2D = (*jEnv)->NewObjectArray(jEnv, aSz.vtY, elmClass,
	  				     (jobject )NULL);
	}
	if(jArray2D)
	{
	  idY = 0;
	  while(jArray1D && (idY < aSz.vtY))
	  {
	    (*jEnv)->SetFloatArrayRegion(jEnv, (jfloatArray )jArray1D,
				    0, aSz.vtX,
				    (jfloat *)*((float **)aVal + idY));
	    (*jEnv)->SetObjectArrayElement(jEnv, jArray2D, idY, jArray1D);
	    jArray1D = (jobject )((*jEnv)->NewFloatArray(jEnv, aSz.vtX));
	    ++idY;
	  }
	}
	if(jArray2D)
	{
	  rtnJObj = jArray2D;
	}
	break;
      case WLZ_JPM_KEY_DOUBLE_ARY2:
        jArray1D = (jobject )((*jEnv)->NewDoubleArray(jEnv, aSz.vtX));
	if(jArray1D)
	{
	  elmClass = (*jEnv)->GetObjectClass(jEnv, jArray1D);
	  jArray2D = (*jEnv)->NewObjectArray(jEnv, aSz.vtY, elmClass,
	  				     (jobject )NULL);
	}
	if(jArray2D)
	{
	  idY = 0;
	  while(jArray1D && (idY < aSz.vtY))
	  {
	    (*jEnv)->SetDoubleArrayRegion(jEnv, (jdoubleArray )jArray1D,
	    			   0, aSz.vtX,
				   (jdouble *)*((double **)aVal + idY));
	    (*jEnv)->SetObjectArrayElement(jEnv, jArray2D, idY, jArray1D);
	    jArray1D = (jobject )((*jEnv)->NewDoubleArray(jEnv, aSz.vtX));
	    ++idY;
	  }
	}
	if(jArray2D)
	{
	  rtnJObj = jArray2D;
	}
	break;
      default:
	break;
    }
  }
  return(rtnJObj);
}

/*!
* \ingroup	JWlz
* \brief	Free's a temporary 2D array.
* \param	aDat			Array data structure.
* \param	aSz			Array size.
* \param	dSKey			Data structure identification.
* \param	isCpy			Copy flag for JNI functions.
*/
void		WlzJavaArray2DFree(void *aDat, WlzIVertex2 aSz,
			      int dSKey, jboolean isCpy)
{
  if(isCpy && aDat && (aSz.vtX > 0) && (aSz.vtY > 0))
  {
    switch(dSKey)
    {
      case WLZ_JPM_KEY_INT_ARY2:   /* FALLTHROUGH */
      case WLZ_JPM_KEY_SHORT_ARY2: /* FALLTHROUGH */
      case WLZ_JPM_KEY_BYTE_ARY2:  /* FALLTHROUGH */
      case WLZ_JPM_KEY_FLOAT_ARY2: /* FALLTHROUGH */
      case WLZ_JPM_KEY_DOUBLE_ARY2:
        Alc2Free((void **)aDat);
	break;
      default:
	break;
    }
  }
}
