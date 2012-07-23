#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HandCodedC/WlzJavaArray1D_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzJavaArray1D.h
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
* \brief	Function prototypes for 1D array access from the C side
* 		of Java Woolz.
* \ingroup	JWlz
*/

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
