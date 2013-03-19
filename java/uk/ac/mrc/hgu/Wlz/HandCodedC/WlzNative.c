#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzNative_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzNative.c
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
* \brief	C binding for Java Woolz Woolz native value access.
* \ingroup	JWlz
*/
#include <WlzJava.h>

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a long value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Java long value to set.
*/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromLong(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  return(jValL);
}

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a int value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValI			Java int value to set.
*/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromInt(
				JNIEnv *jEnv, jclass jCls,
				jint jValI)
{
  jvalue	jUnion;

  jUnion.i = jValI;
  return(jUnion.j);
}

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a short value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValS			Java short value to set.
*/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromShort(
				JNIEnv *jEnv, jclass jCls,
				jshort jValS)
{
  jvalue	jUnion;

  jUnion.s = jValS;
  return(jUnion.j);
}

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a byte value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValB			Java byte value to set.
*/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromByte(
				JNIEnv *jEnv, jclass jCls,
				jbyte jValB)
{
  jvalue	jUnion;

  jUnion.b = jValB;
  return(jUnion.j);
}

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a float value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValF			Java float value to set.
*/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromFloat(
				JNIEnv *jEnv, jclass jCls,
				jfloat jValF)
{
  jvalue	jUnion;

  jUnion.f = jValF;
  return(jUnion.j);
}

/*!
* \return	Java long.
* \ingroup	JWlz
* \brief	Native method implementation to set a double value within a
* 		union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValD			Java double value to set.
*/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromDouble(
				JNIEnv *jEnv, jclass jCls,
				jdouble jValD)
{
  jvalue	jUnion;

  jUnion.d = jValD;
  return(jUnion.j);
}

/*!
* \return	Long value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a long value from within
* 		a  union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToLong(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  return(jValL);
}

/*!
* \return	Int value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a int value from within
* 		a  union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jint JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToInt(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.i);
}

/*!
* \return	Short value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a short value from within
* 		a  union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jshort JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToShort(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.s);
}

/*!
* \return	Byte value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a byte value
*		within a union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jbyte JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToByte(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.b);
}

/*!
* \return	Float value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a float value
*		within a union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jfloat JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToFloat(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.f);
}

/*!
* \return	Double value held.
* \ingroup	JWlz
* \brief	Native method implementation to get a double value
*		within a union.
* \param	jEnv			Given JNI environment pointer.
* \param	jCls			The WlzNative class.
* \param	jValL			Given java long value which contains
* 					the union.
*/
JNIEXPORT jdouble JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToDouble(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.d);
}


