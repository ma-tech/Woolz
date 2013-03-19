#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HandCodedC/WlzJavaArray1D_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzJavaArray1D.c
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
* \brief	Functions for 1D array access between java and native
* 		Woolz code. Many of these functions 'hardcode' Woolz
* 		field names and types, these will need to be maintained.
* \ingroup 	JWlz
*/

#include <WlzJava.h>

static jlong	WlzJavaArray1DGetPrim(
				JNIEnv *jEnv,
				int pKey,
				jarray jWArray,
				int wArraySz,
				jboolean *isCpy);
static jlong	WlzJavaArray1DGetNonPrim(
				JNIEnv *jEnv,
				char *cObjName,
			        char *jObjName,
			        char *jniObjName,
			        int idrCnt, int pKey,
			        jarray jWArray,
			        int wArraySz,
			        jboolean *isCpy);

/*!
* \return	The value held in the given Java object.
* \ingroup	JWlz
* \brief	Returns a native 1D array built from the given java 1D
* 		java array.
* \param	jEnv			Given JNI environment pointer.
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
jlong 		WlzJavaArray1DGet(JNIEnv *jEnv, char *cObjName,
				  char *jObjName, char *jniObjName,
				  int idrCnt, int pKey, jarray jWArray,
				  int wArraySz, jboolean *isCpy)
{
  jlong		rtnVal = 0;

  if(jWArray && (wArraySz > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_BYTE_ARY1:  /* FALLTHROUGH */
      case WLZ_JPM_KEY_SHORT_ARY1: /* FALLTHROUGH */
      case WLZ_JPM_KEY_INT_ARY1:   /* FALLTHROUGH */
      case WLZ_JPM_KEY_LONG_ARY1:  /* FALLTHROUGH */
      case WLZ_JPM_KEY_FLOAT_ARY1: /* FALLTHROUGH */
      case WLZ_JPM_KEY_DOUBLE_ARY1:
	rtnVal = WlzJavaArray1DGetPrim(jEnv, pKey, jWArray, wArraySz, isCpy);
	break;
      default:
	rtnVal = WlzJavaArray1DGetNonPrim(jEnv, cObjName, jObjName, jniObjName,
					  idrCnt, pKey, jWArray, wArraySz,
					  isCpy);
	break;
    }
  }
  return(rtnVal);
}

/*!
* \return	The value held in the given Java object.
* \ingroup	JWlz
* \brief	Returns a native 1D array built from the given java 1D java
* 		array of primative types (byte, short, int, long, float and
* 		double).
* \param	jEnv			Given JNI environment pointer.
* \param	pKey			Parameter key.
* \param	jWArray			The Java Woolz array.
* \param	wArraySz		The number of elements in the Java
* 					Woolz array.
* \param	isCpy			Destination pointer for JNI copy flag.
*/
static jlong	WlzJavaArray1DGetPrim(JNIEnv *jEnv, int pKey,
			   	      jarray jWArray, int wArraySz,
			   	      jboolean *isCpy)
{
  int		idN0;
  void		*bufJ = NULL,
  		*bufW = NULL;

  switch(pKey)
  {
    case WLZ_JPM_KEY_BYTE_ARY1:
      if((bufJ = (void *)(*jEnv)->GetByteArrayElements(jEnv,
      						       (jbyteArray )jWArray,
						       isCpy)) != NULL)
      {
	if((bufW = AlcMalloc(wArraySz * sizeof(jbyte))) != NULL)
	{
	  (void )memcpy(bufW, bufJ, wArraySz * sizeof(jbyte));
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseByteArrayElements(jEnv, (jbyteArray )jWArray,
					    (jbyte *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
    case WLZ_JPM_KEY_SHORT_ARY1:
      if((bufJ = (void *)(*jEnv)->GetShortArrayElements(jEnv,
      						       (jshortArray )jWArray,
						       isCpy)) != NULL)
      {
	if((bufW = AlcMalloc(wArraySz * sizeof(jshort))) != NULL)
	{
	  (void )memcpy(bufW, bufJ, wArraySz * sizeof(jshort));
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseShortArrayElements(jEnv, (jshortArray )jWArray,
					    (jshort *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
    case WLZ_JPM_KEY_INT_ARY1:
      if((bufJ = (void *)(*jEnv)->GetIntArrayElements(jEnv,
      						       (jintArray )jWArray,
						       isCpy)) != NULL)
      {
	if((bufW = AlcMalloc(wArraySz * sizeof(jint))) != NULL)
	{
	  (void )memcpy(bufW, bufJ, wArraySz * sizeof(jint));
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseIntArrayElements(jEnv, (jintArray )jWArray,
					    (jint *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
    case WLZ_JPM_KEY_LONG_ARY1:
      if((bufJ = (void *)(*jEnv)->GetLongArrayElements(jEnv,
      						       (jlongArray )jWArray,
						       isCpy)) != NULL)
      {
        if((bufW = AlcMalloc(wArraySz * sizeof(long))) != NULL)
	{
	  for(idN0 = 0; idN0 < wArraySz; ++idN0)
	  {
	    *((long *)bufW + idN0) = *((jlong *)bufJ + idN0);
	  }
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseLongArrayElements(jEnv, (jlongArray )jWArray,
	  				    (jlong *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
    case WLZ_JPM_KEY_FLOAT_ARY1:
      if((bufJ = (void *)(*jEnv)->GetFloatArrayElements(jEnv,
      						       (jfloatArray )jWArray,
						       isCpy)) != NULL)
      {
	if((bufW = AlcMalloc(wArraySz * sizeof(jfloat))) != NULL)
	{
	  (void )memcpy(bufW, bufJ, wArraySz * sizeof(jfloat));
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseFloatArrayElements(jEnv, (jfloatArray )jWArray,
					    (jfloat *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
    case WLZ_JPM_KEY_DOUBLE_ARY1:
      if((bufJ = (void *)(*jEnv)->GetDoubleArrayElements(jEnv,
      						       (jdoubleArray )jWArray,
						       isCpy)) != NULL)
      {
	if((bufW = AlcMalloc(wArraySz * sizeof(jdouble))) != NULL)
	{
	  (void )memcpy(bufW, bufJ, wArraySz * sizeof(jdouble));
	}
	if(*isCpy)
	{
	  (*jEnv)->ReleaseDoubleArrayElements(jEnv, (jdoubleArray )jWArray,
					    (jdouble *)bufJ, 0);
	}
	*isCpy = JNI_TRUE;
      }
      break;
  }
  return((jlong )bufW);
}

/*!
* \return	The value held in the given Java object.
* \ingroup	JWlz
* \brief	Returns a native 1D array built from the given java 1D java
* 		array of non-primative types.
* \param	jEnv			Given JNI environment pointer.
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
static jlong	WlzJavaArray1DGetNonPrim(JNIEnv *jEnv, char *cObjName,
				         char *jObjName, char *jniObjName,
				         int idrCnt, int pKey, jarray jWArray,
				         int wArraySz, jboolean *isCpy)
{
  int		idN0;
  jlong		rtnVal = 0;
  jfieldID	fldID[8];
  jvalue	fldVal[2];
  jclass	jWCls;
  jobject	jWObj;
  void		*bufW = NULL;

  switch(pKey)
  {
    case WLZ_JPM_KEY_WLZ_IVERTEX2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzIVertex2))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "I");
	((WlzIVertex2 *)bufW)->vtX = (*jEnv)->GetIntField(jEnv, jWObj,
							  fldID[0]);
	((WlzIVertex2 *)bufW)->vtY = (*jEnv)->GetIntField(jEnv, jWObj,
							  fldID[1]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzIVertex2 *)bufW + idN0)->vtX = (*jEnv)->GetIntField(jEnv, jWObj, 
	  							   fldID[0]);
	  ((WlzIVertex2 *)bufW + idN0)->vtY = (*jEnv)->GetIntField(jEnv, jWObj,
	  							   fldID[1]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzFVertex2))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "F");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "F");
	((WlzFVertex2 *)bufW)->vtX = (*jEnv)->GetFloatField(jEnv, jWObj,
							    fldID[0]);
	((WlzFVertex2 *)bufW)->vtY = (*jEnv)->GetFloatField(jEnv, jWObj,
							    fldID[1]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzFVertex2 *)bufW + idN0)->vtX = (*jEnv)->GetFloatField(jEnv,
							      jWObj, fldID[0]);
	  ((WlzFVertex2 *)bufW + idN0)->vtY = (*jEnv)->GetFloatField(jEnv,
	  						      jWObj, fldID[1]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzDVertex2))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "D");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "D");
	((WlzDVertex2 *)bufW)->vtX = (*jEnv)->GetDoubleField(jEnv, jWObj,
							     fldID[0]);
	((WlzDVertex2 *)bufW)->vtY = (*jEnv)->GetDoubleField(jEnv, jWObj,
							     fldID[1]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzDVertex2 *)bufW + idN0)->vtX = (*jEnv)->GetDoubleField(jEnv,
							      jWObj, fldID[0]);
	  ((WlzDVertex2 *)bufW + idN0)->vtY = (*jEnv)->GetDoubleField(jEnv,
							      jWObj, fldID[1]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IVERTEX3_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzIVertex3))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "I");
	((WlzIVertex3 *)bufW)->vtX = (*jEnv)->GetIntField(jEnv, jWObj,
							  fldID[0]);
	((WlzIVertex3 *)bufW)->vtY = (*jEnv)->GetIntField(jEnv, jWObj,
							  fldID[1]);
	((WlzIVertex3 *)bufW)->vtZ = (*jEnv)->GetIntField(jEnv, jWObj,
							  fldID[2]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzIVertex3 *)bufW + idN0)->vtX = (*jEnv)->GetIntField(jEnv, jWObj, 
	  							   fldID[0]);
	  ((WlzIVertex3 *)bufW + idN0)->vtY = (*jEnv)->GetIntField(jEnv, jWObj,
	  							   fldID[1]);
	  ((WlzIVertex3 *)bufW + idN0)->vtZ = (*jEnv)->GetIntField(jEnv, jWObj,
	  							   fldID[2]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX3_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzFVertex3))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "F");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "F");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "F");
	((WlzFVertex3 *)bufW)->vtX = (*jEnv)->GetFloatField(jEnv, jWObj,
							    fldID[0]);
	((WlzFVertex3 *)bufW)->vtY = (*jEnv)->GetFloatField(jEnv, jWObj,
							    fldID[1]);
	((WlzFVertex3 *)bufW)->vtZ = (*jEnv)->GetFloatField(jEnv, jWObj,
							    fldID[2]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzFVertex3 *)bufW + idN0)->vtX = (*jEnv)->GetFloatField(jEnv,
							      jWObj, fldID[0]);
	  ((WlzFVertex3 *)bufW + idN0)->vtY = (*jEnv)->GetFloatField(jEnv,
	  						      jWObj, fldID[1]);
	  ((WlzFVertex3 *)bufW + idN0)->vtZ = (*jEnv)->GetFloatField(jEnv,
	  						      jWObj, fldID[2]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX3_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzDVertex3))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "D");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "D");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "D");
	((WlzDVertex3 *)bufW)->vtX = (*jEnv)->GetDoubleField(jEnv, jWObj,
							     fldID[0]);
	((WlzDVertex3 *)bufW)->vtY = (*jEnv)->GetDoubleField(jEnv, jWObj,
							     fldID[1]);
	((WlzDVertex3 *)bufW)->vtZ = (*jEnv)->GetDoubleField(jEnv, jWObj,
							     fldID[2]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzDVertex3 *)bufW + idN0)->vtX = (*jEnv)->GetDoubleField(jEnv,
							      jWObj, fldID[0]);
	  ((WlzDVertex3 *)bufW + idN0)->vtY = (*jEnv)->GetDoubleField(jEnv,
							      jWObj, fldID[1]);
	  ((WlzDVertex3 *)bufW + idN0)->vtZ = (*jEnv)->GetDoubleField(jEnv,
							      jWObj, fldID[2]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IBOX2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzIBox2))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "I");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "I");
	fldID[3] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "I");
	((WlzIBox2 *)bufW)->xMin = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[0]);
	((WlzIBox2 *)bufW)->yMin = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[1]);
	((WlzIBox2 *)bufW)->xMax = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[2]);
	((WlzIBox2 *)bufW)->yMax = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[3]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzIBox2 *)bufW + idN0)->xMin = (*jEnv)->GetIntField(jEnv, jWObj,
							         fldID[0]);
	  ((WlzIBox2 *)bufW + idN0)->yMin = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[1]);
	  ((WlzIBox2 *)bufW + idN0)->xMax = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[2]);
	  ((WlzIBox2 *)bufW + idN0)->yMax = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[3]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IBOX3_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzIBox3))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "I");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "zMin", "I");
	fldID[3] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "I");
	fldID[4] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "I");
	fldID[5] = (*jEnv)->GetFieldID(jEnv, jWCls, "zMax", "I");
	((WlzIBox3 *)bufW)->xMin = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[0]);
	((WlzIBox3 *)bufW)->yMin = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[1]);
	((WlzIBox3 *)bufW)->zMin = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[2]);
	((WlzIBox3 *)bufW)->xMax = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[3]);
	((WlzIBox3 *)bufW)->yMax = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[4]);
	((WlzIBox3 *)bufW)->zMax = (*jEnv)->GetIntField(jEnv, jWObj,
						        fldID[5]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzIBox3 *)bufW + idN0)->xMin = (*jEnv)->GetIntField(jEnv, jWObj,
							         fldID[0]);
	  ((WlzIBox3 *)bufW + idN0)->yMin = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[1]);
	  ((WlzIBox3 *)bufW + idN0)->zMin = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[2]);
	  ((WlzIBox3 *)bufW + idN0)->xMax = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[3]);
	  ((WlzIBox3 *)bufW + idN0)->yMax = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[4]);
	  ((WlzIBox3 *)bufW + idN0)->zMax = (*jEnv)->GetIntField(jEnv, jWObj,
							  	 fldID[5]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DBOX2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzDBox2))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "D");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "D");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "D");
	fldID[3] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "D");
	((WlzDBox2 *)bufW)->xMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[0]);
	((WlzDBox2 *)bufW)->yMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[1]);
	((WlzDBox2 *)bufW)->xMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[2]);
	((WlzDBox2 *)bufW)->yMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[3]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzDBox2 *)bufW + idN0)->xMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
							         fldID[0]);
	  ((WlzDBox2 *)bufW + idN0)->yMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[1]);
	  ((WlzDBox2 *)bufW + idN0)->xMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[2]);
	  ((WlzDBox2 *)bufW + idN0)->yMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[3]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DBOX3_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzDBox3))) != NULL)
      {
	*isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "D");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "D");
	fldID[2] = (*jEnv)->GetFieldID(jEnv, jWCls, "zMin", "D");
	fldID[3] = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "D");
	fldID[4] = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "D");
	fldID[5] = (*jEnv)->GetFieldID(jEnv, jWCls, "zMax", "D");
	((WlzDBox3 *)bufW)->xMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[0]);
	((WlzDBox3 *)bufW)->yMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[1]);
	((WlzDBox3 *)bufW)->zMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[2]);
	((WlzDBox3 *)bufW)->xMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[3]);
	((WlzDBox3 *)bufW)->yMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[4]);
	((WlzDBox3 *)bufW)->zMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
						        fldID[5]);
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  ((WlzDBox3 *)bufW + idN0)->xMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
							         fldID[0]);
	  ((WlzDBox3 *)bufW + idN0)->yMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[1]);
	  ((WlzDBox3 *)bufW + idN0)->zMin = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[2]);
	  ((WlzDBox3 *)bufW + idN0)->xMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[3]);
	  ((WlzDBox3 *)bufW + idN0)->yMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[4]);
	  ((WlzDBox3 *)bufW + idN0)->zMax = (*jEnv)->GetDoubleField(jEnv, jWObj,
							  	 fldID[5]);
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_PIXELV_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(WlzPixelV))) != NULL)
      {
        *isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "type", "I");
	fldID[1] = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	fldVal[0].i = (*jEnv)->GetIntField(jEnv, jWObj, fldID[0]);
	fldVal[1].j = (*jEnv)->GetLongField(jEnv, jWObj, fldID[1]);
	((WlzPixelV *)bufW)->type = fldVal[0].i;
	switch(fldVal[0].i)
	{
          case WLZ_GREY_LONG:
            ((WlzPixelV *)bufW)->v.lnv = fldVal[1].j;
            break;
          case WLZ_GREY_INT:
            ((WlzPixelV *)bufW)->v.inv = fldVal[1].i;
            break;
          case WLZ_GREY_SHORT:
            ((WlzPixelV *)bufW)->v.shv = fldVal[1].s;
            break;
          case WLZ_GREY_UBYTE:
            ((WlzPixelV *)bufW)->v.ubv = fldVal[1].b;
            break;
          case WLZ_GREY_FLOAT:
            ((WlzPixelV *)bufW)->v.flv = fldVal[1].f;
            break;
          case WLZ_GREY_DOUBLE:
            ((WlzPixelV *)bufW)->v.dbv = fldVal[1].d;
            break;
	}
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  fldVal[0].i = (*jEnv)->GetIntField(jEnv, jWObj, fldID[0]);
	  fldVal[1].j = (*jEnv)->GetLongField(jEnv, jWObj, fldID[1]);
	  ((WlzPixelV *)bufW + idN0)->type = fldVal[0].i;
	  switch(fldVal[0].i)
	  { 
          case WLZ_GREY_LONG:
            ((WlzPixelV *)bufW + idN0)->v.lnv = fldVal[1].j;
            break;
          case WLZ_GREY_INT:
            ((WlzPixelV *)bufW + idN0)->v.inv = fldVal[1].i;
            break;
          case WLZ_GREY_SHORT:
            ((WlzPixelV *)bufW + idN0)->v.shv = fldVal[1].s;
            break;
          case WLZ_GREY_UBYTE:
            ((WlzPixelV *)bufW + idN0)->v.ubv = fldVal[1].b;
            break;
          case WLZ_GREY_FLOAT:
            ((WlzPixelV *)bufW + idN0)->v.flv = fldVal[1].f;
            break;
          case WLZ_GREY_DOUBLE:
            ((WlzPixelV *)bufW + idN0)->v.dbv = fldVal[1].d;
            break;
	  }
	}
	(*jEnv)->DeleteLocalRef(jEnv, jWObj);
      }
      break;
    case WLZ_JPM_KEY_WLZ_PTR1_ARY1: /* FALLTHROUGH */
    case WLZ_JPM_KEY_WLZ_PTR2_ARY1:
      if((bufW = AlcMalloc(wArraySz * sizeof(void *))) != NULL)
      {
        *isCpy = JNI_TRUE;
	jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, 0);
	jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	fldID[0] = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	*(void **)bufW = (void *)((*jEnv)->GetLongField(jEnv, jWObj,
							fldID[0]));
	for(idN0 = 1; idN0 < wArraySz; ++idN0)
	{
	  (*jEnv)->DeleteLocalRef(jEnv, jWObj);
	  jWObj = (*jEnv)->GetObjectArrayElement(jEnv, jWArray, idN0);
	  *((void **)bufW + idN0) = (void *)((*jEnv)->GetLongField(jEnv,
							     jWObj, fldID[0]));
	}
      }
      break;
    default:
      break;
  }
  rtnVal = (long )bufW;
  return(rtnVal);
}

/*!
* \ingroup	JWlz
* \brief	Sets the given value in a Java object for return to the Java
* 		side of the JNI.
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
* \param	aSz			Size of the 1D array.
*/
void		WlzJavaArray1DSet(JNIEnv *jEnv, jobjectArray dstJObj,
				char *cObjName,
				char *jObjName,
				char *jniObjName,
				int idrCnt, int pKey,
				void *aVal,
				int aSz)
{
  jobject	newJObj = NULL;

  if(aVal && (aSz > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_BYTE_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "char", "jbyte", "jbyteArray",
				     2, WLZ_JPM_KEY_BYTE_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_SHORT_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "short", "jshort", "jshortArray",
				     2, WLZ_JPM_KEY_LONG_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_INT_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "int", "jint", "jintArray",
				     2, WLZ_JPM_KEY_INT_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_LONG_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "long", "jlong", "jlongArray",
				     2, WLZ_JPM_KEY_LONG_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_FLOAT_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "float", "jfloat", "jfloatArray",
				     2, WLZ_JPM_KEY_FLOAT_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_DOUBLE_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, "double", "jdouble", "jdoubleArray",
				     2, WLZ_JPM_KEY_DOUBLE_ARY1,
				     aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_IVERTEX2_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzIVertex2", "WlzIVertex2", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_IVERTEX2_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_FVERTEX2_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzFVertex2", "WlzFVertex2", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_FVERTEX2_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_DVERTEX2_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzDVertex2", "WlzDVertex2", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_DVERTEX2_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_IVERTEX3_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzIVertex3", "WlzIVertex3", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_IVERTEX3_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_FVERTEX3_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzFVertex3", "WlzFVertex3", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_FVERTEX3_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_DVERTEX3_PTR1_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv,
				"WlzDVertex3", "WlzDVertex3", "jobjectArray",
			        2, WLZ_JPM_KEY_WLZ_DVERTEX3_ARY1,
			        aVal, aSz);
	break;
      case WLZ_JPM_KEY_WLZ_PTR2_ARY1:
	newJObj = WlzJavaArray1DWrap(jEnv, cObjName, jObjName, "jobjectArray",
				     2, WLZ_JPM_KEY_WLZ_PTR1_ARY1,
				     aVal, aSz);
	break;
      default:
	break;
    }
  }
  if(newJObj)
  {
    (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
  }
}

/*!
* \return	New Java object.
* \ingroup	JWlz
* \brief	Wraps up the given array as a new Java array.
* \param	jEnv			Given JNI environment pointer.
* \param	cObjName		The Java woolz object class string.
* \param	jObjName		The Java woolz Java object class
* 					string.
* \param	jniObjName		The Java woolz JNI object class string.
* \param	idrCnt			Indirection count (ie 1 for *, 2 for
* 					**, ...).
* \param	pKey			Parameter key.
* \param	aVal			C array to set values from.
* \param	aSz			Size of the 1D array.
*/
jobject		WlzJavaArray1DWrap(JNIEnv *jEnv,
				   char *cObjName,
				   char *jObjName,
				   char *jniObjName,
				   int idrCnt, int pKey,
				   void *aVal,
				   int aSz)
{
  int		isWlzObject,
  		idN0;
  char		*jWFqClassName;
  jclass	jWCls;
  jmethodID	mtdID;
  jfieldID	fldID;
  jobject	tObj0,
  		rtnJObj = NULL;

#ifdef JWLZ_DEBUG
  int hack = 1;
  while(hack)
  {
    sleep(1);
  }
#endif /* JWLZ_DEBUG */
  if(aVal && (aSz > 0))
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_BYTE_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewByteArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetByteArrayRegion(jEnv, (jbyteArray )rtnJObj, 0, aSz,
				      aVal);
	}
	break;
      case WLZ_JPM_KEY_SHORT_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewShortArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetShortArrayRegion(jEnv, (jshortArray )rtnJObj, 0, aSz,
				       aVal);
	}
	break;
      case WLZ_JPM_KEY_INT_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewIntArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetIntArrayRegion(jEnv, (jintArray )rtnJObj, 0, aSz,
				     aVal);
	}
	break;
      case WLZ_JPM_KEY_LONG_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewLongArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetLongArrayRegion(jEnv, (jlongArray )rtnJObj, 0, aSz,
				      aVal);
	}
	break;
      case WLZ_JPM_KEY_FLOAT_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewFloatArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetFloatArrayRegion(jEnv, (jfloatArray )rtnJObj, 0, aSz,
				       aVal);
	}
	break;
      case WLZ_JPM_KEY_DOUBLE_ARY1:
	if((rtnJObj = (jobject )((*jEnv)->NewDoubleArray(jEnv, aSz))) != NULL)
	{
	  (*jEnv)->SetDoubleArrayRegion(jEnv, (jdoubleArray )rtnJObj, 0, aSz,
					aVal);
	}
	break;
      case WLZ_JPM_KEY_WLZ_IVERTEX2_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(II)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzIVertex2 *)aVal + idN0)->vtX,
				       ((WlzIVertex2 *)aVal + idN0)->vtY);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_FVERTEX2_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(FF)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzFVertex2 *)aVal + idN0)->vtX,
				       ((WlzFVertex2 *)aVal + idN0)->vtY);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_DVERTEX2_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DD)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzDVertex2 *)aVal + idN0)->vtX,
				       ((WlzDVertex2 *)aVal + idN0)->vtY);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_IVERTEX3_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(III)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzIVertex3 *)aVal + idN0)->vtX,
				       ((WlzIVertex3 *)aVal + idN0)->vtY,
				       ((WlzIVertex3 *)aVal + idN0)->vtZ);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_FVERTEX3_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(FFF)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzFVertex3 *)aVal + idN0)->vtX,
				       ((WlzFVertex3 *)aVal + idN0)->vtY,
				       ((WlzFVertex3 *)aVal + idN0)->vtZ);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_DVERTEX3_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DDD)V");
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
	    			       ((WlzDVertex3 *)aVal + idN0)->vtX,
				       ((WlzDVertex3 *)aVal + idN0)->vtY,
				       ((WlzDVertex3 *)aVal + idN0)->vtZ);
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
	}
        break;
      case WLZ_JPM_KEY_WLZ_PTR1_ARY1:
	if((jWFqClassName = WlzJavaBuildFQClassName(cObjName)) != NULL)
	{
	  jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	  mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "()V");
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	  rtnJObj = (jobject )(*jEnv)->NewObjectArray(jEnv, aSz, jWCls, NULL);
	  AlcFree(jWFqClassName);
	}
	if(rtnJObj)
	{
	  isWlzObject = strcmp(cObjName, "WlzObject") == 0;
	  for(idN0 = 0; idN0 < aSz; ++idN0)
	  {
	    tObj0 = (*jEnv)->NewObject(jEnv, jWCls, mtdID);
	    if(isWlzObject)
	    {
	      (void )WlzAssignObject(*((WlzObject **)aVal + idN0), NULL);
	    }
	    (*jEnv)->SetLongField(jEnv, tObj0, fldID,
				  (jlong )(*((WlzObject **)aVal + idN0)));
	    (*jEnv)->SetObjectArrayElement(jEnv, (jobjectArray )rtnJObj,
					   idN0, tObj0);
	  }
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
* \brief	Free's a temporary native 1D array.
* \param	aDat			Array data structure.
* \param	aSz			Array size.
* \param	dSKey			Data structure identification.
* \param	isCpy			Copy flag for JNI functions.
*/
void		WlzJavaArray1DFree(void *aDat, int aSz,
			      int dSKey, jboolean isCpy)
{
  if(isCpy && aDat && (aSz > 0))
  {
    AlcFree(aDat);
  }
}
