#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzJavaArray3D_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzJavaArray3D.h
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
* \brief	Function prototypes for 3D array access from the C side
* 		of Java Woolz.
* \ingroup	JWlz
*/

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
