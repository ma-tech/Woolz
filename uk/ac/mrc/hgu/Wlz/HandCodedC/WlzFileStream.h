/************************************************************************
* Project:      Java Woolz
* Title:        WlzFileStream.h
* Date:         January 1999
* Purpose:      Function prototypes C binding for Java Woolz native
*		file IO.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/
JNIEXPORT void JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzClose(
				JNIEnv *jEnv, jobject jObj,
				jlong jValue);
JNIEXPORT jlong JNICALL	Java_uk_ac_mrc_hgu_Wlz_WlzFileStream_JWlzOpen(
				JNIEnv *jEnv, jobject jObj,
				jstring jName, jstring jMode);
