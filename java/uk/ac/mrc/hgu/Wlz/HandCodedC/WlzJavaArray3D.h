/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaArray3D.h
* Date:         January 1999
* Purpose:      Function prototypes for 3D array access from the C side
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

extern jlong 			WlzJavaArray3DGet(
				  JNIEnv *jEnv,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  jarray jWArray,
				  WlzIVertex3 wArraySz,
				  jboolean *isCpy);
extern void			WlzJavaArray3DSet(
				  JNIEnv *jEnv,
				  jobjectArray dstJObj,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  void *aVal,
				  WlzIVertex3 aSz);
extern jobject			WlzJavaArray3DWrap(JNIEnv *jEnv,
				   char *cObjName,
				   char *jObjName,
				   char *jniObjName,
				   int idrCnt, int pKey,
				   void *aVal,
				   WlzIVertex3 aSz);
extern void			WlzJavaArray3DFree(
				  void *aDat,
				  WlzIVertex3 aSz,
				  int dSKey,
				  jboolean isCopy);
