#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _HandCodedC/WlzFileStream_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzFileStream.h
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
* \brief 	Function prototypes C binding for Java Woolz native
*		file IO.
* \ingroup	JWlz
*/

JNIEXPORT void JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzClose(
				JNIEnv *jEnv, jobject jObj,
				jlong jValue);
JNIEXPORT jlong JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzOpen(
				JNIEnv *jEnv, jobject jObj,
				jstring jName, jstring jMode);
