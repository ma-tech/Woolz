/************************************************************************
* Project:      Java Woolz
* Title:        WlzNative.c
* Date:         January 1999
* Purpose:      C binding for Java Woolz Woolz native value access.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/
#include <WlzJava.h>


/************************************************************************
* Function:    	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromLong
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromLong
* Purpose:	Native method implementation to set a long value within
*		a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValue:		Java long value to set.
************************************************************************/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromLong(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  return(jValL);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromInt
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromInt
* Purpose:	Native method implementation to set a int value within
*		a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jint jValue:		Java int value to set.
************************************************************************/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromInt(
				JNIEnv *jEnv, jclass jCls,
				jint jValI)
{
  jvalue	jUnion;

  jUnion.i = jValI;
  return(jUnion.j);
}

/************************************************************************
* Function:   	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromShort
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromShort
* Purpose:	Native method implementation to set a short value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jshort jValue:		Java short value to set.
************************************************************************/
JNIEXPORT jlong JNICALL	
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromShort(
				JNIEnv *jEnv, jclass jCls,
				jshort jValS)
{
  jvalue	jUnion;

  jUnion.s = jValS;
  return(jUnion.j);
}

/************************************************************************
* Function:    	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromByte
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromByte
* Purpose:	Native method implementation to set a byte value within
*		a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jbyte jValue:		Java byte value to set.
************************************************************************/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromByte(
				JNIEnv *jEnv, jclass jCls,
				jbyte jValB)
{
  jvalue	jUnion;

  jUnion.b = jValB;
  return(jUnion.j);
}

/************************************************************************
* Function:   	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromFloat
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromFloat
* Purpose:	Native method implementation to set a float value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jfloat jValue:		Java float value to set.
************************************************************************/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromFloat(
				JNIEnv *jEnv, jclass jCls,
				jfloat jValF)
{
  jvalue	jUnion;

  jUnion.f = jValF;
  return(jUnion.j);
}

/************************************************************************
* Function:  	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromDouble
* Returns:	jlong
* Class:	WlzNative
* Method:	JWlzNativeUnionFromDouble
* Purpose:	Native method implementation to set a double value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jdouble jValue:		Java double value to set.
************************************************************************/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionFromDouble(
				JNIEnv *jEnv, jclass jCls,
				jdouble jValD)
{
  jvalue	jUnion;

  jUnion.d = jValD;
  return(jUnion.j);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToLong
* Returns:	jlong:			Long value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToLong
* Purpose:	Native method implementation to get a long value
*		from within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jlong JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToLong(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  return(jValL);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToInt
* Returns:	jint:			Int value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToInt
* Purpose:	Native method implementation to get a int value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jint JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToInt(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.i);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToShort
* Returns:	jshort:			Short value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToShort
* Purpose:	Native method implementation to get a short value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jshort JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToShort(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.s);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToByte
* Returns:	jbyte:			Byte value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToByte
* Purpose:	Native method implementation to get a byte value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jbyte JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToByte(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.b);
}

/************************************************************************
* Function:	Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToFloat
* Returns:	jfloat:			Float value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToFloat
* Purpose:	Native method implementation to get a float value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jfloat JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToFloat(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.f);
}

/************************************************************************
* Function:    Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToDouble
* Returns:	jdouble:		Double value held.
* Class:	WlzNative
* Method:	JWlzNativeUnionToDouble
* Purpose:	Native method implementation to get a double value
*		within a union.
* Global refs:	-
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.
*		jclass jCls:		The WlzNative class.
*		jlong jValL:		Given java long value which
*					contains the union.
************************************************************************/
JNIEXPORT jdouble JNICALL
Java_uk_ac_mrc_hgu_Wlz_WlzNative_JWlzNativeUnionToDouble(
				JNIEnv *jEnv, jclass jCls,
				jlong jValL)
{
  jvalue	jUnion;

  jUnion.j = jValL;
  return(jUnion.d);
}


