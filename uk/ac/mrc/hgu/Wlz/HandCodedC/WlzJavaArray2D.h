/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaArray2D.h
* Date:         January 1999
* Purpose:      Function prototypes for 2D array access from the C side
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

extern jlong 			WlzJavaArray2DGet(
				  JNIEnv *jEnv,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  jarray jWArray,
				  WlzIVertex2 wArraySz,
				  jboolean *isCpy);
extern void			WlzJavaArray2DSet(
				  JNIEnv *jEnv,
				  jobjectArray dstJObj,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  void *aVal,
				  WlzIVertex2 aSz);
extern jobject			WlzJavaArray2DWrap(JNIEnv *jEnv,
				   char *cObjName,
				   char *jObjName,
				   char *jniObjName,
				   int idrCnt, int pKey,
				   void *aVal,
				   WlzIVertex2 aSz);
extern void			WlzJavaArray2DFree(
				  void *aDat,
				  WlzIVertex2 aSz,
				  int dSKey,
				  jboolean isCopy);
