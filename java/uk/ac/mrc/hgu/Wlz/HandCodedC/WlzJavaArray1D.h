/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaArray1D.h
* Date:         January 1999
* Purpose:      Function prototypes for 1D array access from the C side
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

extern jlong 			WlzJavaArray1DGet(
				  JNIEnv *jEnv,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  jarray jWArray,
				  int wArraySz,
				  jboolean *isCpy);
extern void			WlzJavaArray1DSet(
				  JNIEnv *jEnv,
				  jobjectArray dstJObj,
				  char *cObjName,
				  char *jObjName,
				  char *jniObjName,
				  int idrCnt, int pKey,
				  void *aVal,
				  int aSz);
extern jobject			WlzJavaArray1DWrap(JNIEnv *jEnv,
				   char *cObjName,
				   char *jObjName,
				   char *jniObjName,
				   int idrCnt, int pKey,
				   void *aVal,
				   int aSz);
extern void			WlzJavaArray1DFree(
				  void *aDat,
				  int aSz,
				  int dSKey,
				  jboolean isCopy);
