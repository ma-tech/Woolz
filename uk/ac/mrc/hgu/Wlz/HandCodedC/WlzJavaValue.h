/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaValue.h
* Date:         January 1999
* Purpose:      Prototypes of value handling functions for the C side
*		of Java Woolz.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/

extern jlong 			WlzJavaValueGet(
				  JNIEnv *jEnv,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  jobject jWObj);
extern jobject 			WlzJavaValueWrap(
				  JNIEnv *jEnv,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  long val);
extern void			WlzJavaValueSet(
				  JNIEnv *jEnv,
				  jobjectArray dstJObj,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  long val);
extern void			WlzJavaValueFree(
				  void *dSP,
				  int dSKey,
				  jboolean isCopy);
