/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaValue.c
* Date:         August 1999
* Purpose:      Value handling functions for the C side of Java Woolz.
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
* Function:	WlzJavaValueGet						*
* Returns:	jlong:			The value held in the given 	*
*					Java object.			*
* Purpose:	Returns the value held by the given Woolz Java object.	*
* Global refs:	-							*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
*		char *cObjName:		The Java woolz C object class	*
*					string.				*
*		char *jObjName:		The Java woolz Java object	*
*					class string.			*
*		char *jniObjName:	The Java woolz JNI object	*
*					class string.			*
*		int idrCnt:		Indirection count (ie 1 for *,	*
*					2 for **, ...).			*
*		int pKey:		Parameter key.			*
*		jobject jWObj:		The Java Woolz object.		*
************************************************************************/
jlong 		WlzJavaValueGet(JNIEnv *jEnv,
				char *cObjName,
				char *jObjName,
				char *jniObjName,
				int idrCnt, int pKey,
				jobject jWObj)
{
  jvalue	fldVal[16];
  jlong		rtnVal = 0;
  jfieldID	fldID;
  jclass	jWCls;

  if(jWObj)
  {
    switch(pKey)
    {
      case WLZ_JPM_KEY_WLZ_GREYV:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzGreyV))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	  fldVal[0].j = (*jEnv)->GetLongField(jEnv, jWObj, fldID);
	  (*(WlzGreyV *)rtnVal).lnv = fldVal[0].j;
	}
	break;
      case WLZ_JPM_KEY_WLZ_PIXELV:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzPixelV))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "type", "I");
	  fldVal[0].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	  fldVal[1].j = (*jEnv)->GetLongField(jEnv, jWObj, fldID);
	  ((WlzPixelV *)rtnVal)->type = fldVal[0].i;
	  switch(fldVal[0].i)
	  {
	    case WLZ_GREY_LONG:
	      ((WlzPixelV *)rtnVal)->v.lnv = fldVal[1].j;
	      break;
	    case WLZ_GREY_INT:
	      ((WlzPixelV *)rtnVal)->v.inv = fldVal[1].i;
	      break;
	    case WLZ_GREY_SHORT:
	      ((WlzPixelV *)rtnVal)->v.shv = fldVal[1].s;
	      break;
	    case WLZ_GREY_UBYTE:
	      ((WlzPixelV *)rtnVal)->v.ubv = fldVal[1].b;
	      break;
	    case WLZ_GREY_FLOAT:
	      ((WlzPixelV *)rtnVal)->v.flv = fldVal[1].f;
	      break;
	    case WLZ_GREY_DOUBLE:
	      ((WlzPixelV *)rtnVal)->v.dbv = fldVal[1].d;
	      break;
	  }
	}
	break;
      case WLZ_JPM_KEY_WLZ_IVERTEX2:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzIVertex2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "I");
	  fldVal[0].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "I");
	  fldVal[1].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  ((WlzIVertex2 *)rtnVal)->vtX = fldVal[0].i;
	  ((WlzIVertex2 *)rtnVal)->vtY = fldVal[1].i;
	}
	break;
      case WLZ_JPM_KEY_WLZ_FVERTEX2:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzFVertex2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "F");
	  fldVal[0].f = (float )((*jEnv)->GetFloatField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "F");
	  fldVal[1].f = (float )((*jEnv)->GetFloatField(jEnv, jWObj, fldID));
	  ((WlzFVertex2 *)rtnVal)->vtX = fldVal[0].f;
	  ((WlzFVertex2 *)rtnVal)->vtY = fldVal[1].f;
	}
	break;
      case WLZ_JPM_KEY_WLZ_DVERTEX2:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzDVertex2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "D");
	  fldVal[0].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "D");
	  fldVal[1].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  ((WlzDVertex2 *)rtnVal)->vtX = fldVal[0].i;
	  ((WlzDVertex2 *)rtnVal)->vtY = fldVal[1].i;
	}
	break;
      case WLZ_JPM_KEY_WLZ_IVERTEX3:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzIVertex3))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "I");
	  fldVal[0].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "I");
	  fldVal[1].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "I");
	  fldVal[2].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  ((WlzIVertex3 *)rtnVal)->vtX = fldVal[0].i;
	  ((WlzIVertex3 *)rtnVal)->vtY = fldVal[1].i;
	  ((WlzIVertex3 *)rtnVal)->vtZ = fldVal[2].i;
	}
	break;
      case WLZ_JPM_KEY_WLZ_FVERTEX3:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzFVertex3))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "F");
	  fldVal[0].f = (float )((*jEnv)->GetFloatField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "F");
	  fldVal[1].f = (float )((*jEnv)->GetFloatField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "F");
	  fldVal[2].f = (float )((*jEnv)->GetFloatField(jEnv, jWObj, fldID));
	  ((WlzFVertex3 *)rtnVal)->vtX = fldVal[0].f;
	  ((WlzFVertex3 *)rtnVal)->vtY = fldVal[1].f;
	  ((WlzFVertex3 *)rtnVal)->vtZ = fldVal[2].f;
	}
	break;
      case WLZ_JPM_KEY_WLZ_DVERTEX3:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzDVertex3))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtX", "D");
	  fldVal[0].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtY", "D");
	  fldVal[1].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "vtZ", "D");
	  fldVal[2].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  ((WlzDVertex3 *)rtnVal)->vtX = fldVal[0].d;
	  ((WlzDVertex3 *)rtnVal)->vtY = fldVal[1].d;
	  ((WlzDVertex3 *)rtnVal)->vtZ = fldVal[2].d;
	}
	break;
      case WLZ_JPM_KEY_WLZ_IBOX2:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzIBox2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "I");
	  fldVal[0].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "I");
	  fldVal[1].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "I");
	  fldVal[2].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "I");
	  fldVal[3].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  ((WlzIBox2 *)rtnVal)->xMin = fldVal[0].i;
	  ((WlzIBox2 *)rtnVal)->xMax = fldVal[1].i;
	  ((WlzIBox2 *)rtnVal)->yMin = fldVal[2].i;
	  ((WlzIBox2 *)rtnVal)->yMax = fldVal[3].i;
	}
	break;
      case WLZ_JPM_KEY_WLZ_DBOX2:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzDBox2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "D");
	  fldVal[0].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "D");
	  fldVal[1].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "D");
	  fldVal[2].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "D");
	  fldVal[3].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  ((WlzDBox2 *)rtnVal)->xMin = fldVal[0].d;
	  ((WlzDBox2 *)rtnVal)->xMax = fldVal[1].d;
	  ((WlzDBox2 *)rtnVal)->yMin = fldVal[2].d;
	  ((WlzDBox2 *)rtnVal)->yMax = fldVal[3].d;
	}
	break;
      case WLZ_JPM_KEY_WLZ_IBOX3:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzIBox3))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "I");
	  fldVal[0].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "I");
	  fldVal[1].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "I");
	  fldVal[2].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "I");
	  fldVal[3].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "zMin", "I");
	  fldVal[4].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "zMax", "I");
	  fldVal[5].i = (int )((*jEnv)->GetIntField(jEnv, jWObj, fldID));
	  ((WlzIBox3 *)rtnVal)->xMin = fldVal[0].i;
	  ((WlzIBox3 *)rtnVal)->xMax = fldVal[1].i;
	  ((WlzIBox3 *)rtnVal)->yMin = fldVal[2].i;
	  ((WlzIBox3 *)rtnVal)->yMax = fldVal[3].i;
	  ((WlzIBox3 *)rtnVal)->zMin = fldVal[4].i;
	  ((WlzIBox3 *)rtnVal)->zMax = fldVal[5].i;
	}
	break;
      case WLZ_JPM_KEY_WLZ_DBOX3:
	if((rtnVal = (long )AlcMalloc(sizeof(WlzDBox2))) != NULL)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMin", "D");
	  fldVal[0].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "xMax", "D");
	  fldVal[1].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMin", "D");
	  fldVal[2].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "yMax", "D");
	  fldVal[3].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "zMin", "D");
	  fldVal[5].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "zMax", "D");
	  fldVal[5].d = (double )((*jEnv)->GetDoubleField(jEnv, jWObj, fldID));
	  ((WlzDBox3 *)rtnVal)->xMin = fldVal[0].d;
	  ((WlzDBox3 *)rtnVal)->xMax = fldVal[1].d;
	  ((WlzDBox3 *)rtnVal)->yMin = fldVal[2].d;
	  ((WlzDBox3 *)rtnVal)->yMax = fldVal[3].d;
	  ((WlzDBox3 *)rtnVal)->zMin = fldVal[4].d;
	  ((WlzDBox3 *)rtnVal)->zMax = fldVal[5].d;
	}
	break;
      case WLZ_JPM_KEY_WLZ_PTR1:
	if(jWObj)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	  rtnVal = (long )((*jEnv)->GetLongField(jEnv, jWObj, fldID));
	}
	break;
      case WLZ_JPM_KEY_WLZ_FILE_PTR1:
	if(jWObj)
	{
	  jWCls = (*jEnv)->GetObjectClass(jEnv, jWObj);
	  fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	  rtnVal = (long )((*jEnv)->GetLongField(jEnv, jWObj, fldID));
	}
	break;
      default:
	break;
    }
  }
  return(rtnVal);
}

/************************************************************************
* Function:	WlzJavaValueWrap					*
* Returns:	jobject:		New Java object.		*
* Purpose:	Wraps up the given value to create a new Java Woolz	*
*		object of the required class.				*
* Global refs:	-							*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
*		char *cObjName:		The Java woolz object class	*
*					string.				*
*               char *jObjName:         The Java woolz Java object      *
*                                       class string.                   *
*               char *jniObjName:       The Java woolz JNI object       *
*                                       class string.                   *
*		int idrCnt:		Indirection count (ie 1 for *,	*
*					2 for **, ...).			*
*		int pKey:		Parameter key.			*
*		long val:		Woolz value to be wrapped up	*
*					into a Java object.		*
************************************************************************/
jobject 	WlzJavaValueWrap(JNIEnv *jEnv,
				 char *cObjName,
				 char *jObjName,
				 char *jniObjName,
				 int idrCnt, int pKey,
				 long val)
{
  char 		*jWFqClassName;
  jclass	jWCls;
  jmethodID	mtdID;
  jfieldID	fldID;
  jobject	rtnJObj = NULL;

  switch(pKey)
  {
    case WLZ_JPM_KEY_BYTE:
    case WLZ_JPM_KEY_SHORT:
      break;
    case WLZ_JPM_KEY_INT:
    case WLZ_JPM_KEY_INT_PTR1:
      jWCls = (*jEnv)->FindClass(jEnv, "java/lang/Integer");
      mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(I)V");
      rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID, *(int *)val);
      break;
    case WLZ_JPM_KEY_LONG:
    case WLZ_JPM_KEY_FLOAT:
      break;
    case WLZ_JPM_KEY_DOUBLE:
    case WLZ_JPM_KEY_DOUBLE_PTR1:
      jWCls = (*jEnv)->FindClass(jEnv, "java/lang/Double");
      mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(D)V");
      rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID, *(double *)val);
      break;
    case WLZ_JPM_KEY_WLZ_GREYV:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(J)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     (*(WlzGreyV *)val).lnv);
      }
      break;
    case WLZ_JPM_KEY_WLZ_PIXELV:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(IJ)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     (*(WlzPixelV *)val).type,
				     (*(WlzPixelV *)val).v.lnv);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IVERTEX2:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(II)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzIVertex2 *)val)->vtX,
				     ((WlzIVertex2 *)val)->vtY);
      }
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX2:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(FF)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzFVertex2 *)val)->vtX,
				     ((WlzFVertex2 *)val)->vtY);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX2:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DD)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzDVertex2 *)val)->vtX,
				     ((WlzDVertex2 *)val)->vtY);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IVERTEX3:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(III)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzIVertex3 *)val)->vtX,
				     ((WlzIVertex3 *)val)->vtY,
				     ((WlzIVertex3 *)val)->vtZ);
      }
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX3:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(FFF)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzFVertex3 *)val)->vtX,
				     ((WlzFVertex3 *)val)->vtY,
				     ((WlzFVertex3 *)val)->vtZ);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX3:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DDD)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzDVertex3 *)val)->vtX,
				     ((WlzDVertex3 *)val)->vtY,
				     ((WlzDVertex3 *)val)->vtZ);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IBOX2:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(IIII)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzIBox2 *)val)->xMin,
				     ((WlzIBox2 *)val)->yMin,
				     ((WlzIBox2 *)val)->xMax,
				     ((WlzIBox2 *)val)->yMax);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DBOX2:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DDDD)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzIBox2 *)val)->xMin,
				     ((WlzIBox2 *)val)->yMin,
				     ((WlzIBox2 *)val)->xMax,
				     ((WlzIBox2 *)val)->yMax);
      }
      break;
    case WLZ_JPM_KEY_WLZ_IBOX3:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(IIIIII)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzIBox3 *)val)->xMin,
				     ((WlzIBox3 *)val)->yMin,
				     ((WlzIBox3 *)val)->zMin,
				     ((WlzIBox3 *)val)->xMax,
				     ((WlzIBox3 *)val)->yMax,
				     ((WlzIBox3 *)val)->zMax);
      }
      break;
    case WLZ_JPM_KEY_WLZ_DBOX3:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
        jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "(DDDDDD)V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID,
				     ((WlzDBox3 *)val)->xMin,
				     ((WlzDBox3 *)val)->yMin,
				     ((WlzDBox3 *)val)->zMin,
				     ((WlzDBox3 *)val)->xMax,
				     ((WlzDBox3 *)val)->yMax,
				     ((WlzDBox3 *)val)->zMax);
      }
      break;
    case WLZ_JPM_KEY_STRING:
      rtnJObj = (jobject )
		((*jEnv)->NewStringUTF(jEnv,
				       (val)? (const char *)val: ""));
      break;
    case WLZ_JPM_KEY_WLZ_PTR1:
      if((jWFqClassName = WlzJavaBuildFQClassName(jObjName)) != NULL)
      {
	jWCls = (*jEnv)->FindClass(jEnv, jWFqClassName);
	mtdID = (*jEnv)->GetMethodID(jEnv, jWCls, "<init>", "()V");
	rtnJObj = (*jEnv)->NewObject(jEnv, jWCls, mtdID);
	fldID = (*jEnv)->GetFieldID(jEnv, jWCls, "value", "J");
	AlcFree(jWFqClassName);
	if(strcmp(cObjName, "WlzObject") == 0)
	{
	  (void )WlzAssignObject((WlzObject *)val, NULL);
	}
	(*jEnv)->SetLongField(jEnv, rtnJObj, fldID, val);
      }
      break;
    default:
      break;
  }
  return(rtnJObj);
}

/************************************************************************
* Function:	WlzJavaValueSet						*
* Returns:	void							*
* Purpose:	Sets the given value in a Java object for return to	*
*		the Java side of th JNI.				*
* Global refs:	-							*
* Parameters:	JNIEnv *jEnv:		Given JNI environment ptr.	*
*		jobjectArray dstJA:	Destination Java array.		*
*		char *cObjName:		The type of the value to set.	*
*               char *jObjName:         The Java woolz Java object      *
*                                       class string.                   *
*               char *jniObjName:       The Java woolz JNI object       *
*                                       class string.                   *
*		int idrCnt:		Indirection count (ie 1 for *,	*
*					2 for **, ...).			*
*		int pKey:		Parameter key.			*
*		long val:		Value to be set.		*
************************************************************************/
void		WlzJavaValueSet(JNIEnv *jEnv, jobjectArray dstJObj,
				char *cObjName,
				char *jObjName,
				char *jniObjName,
				int idrCnt, int pKey,
				long val)
{
  jboolean	isCpy;
  void		*aElmP;
  jobject	newJObj;

  switch(pKey)
  {
    case WLZ_JPM_KEY_LONG:
    case WLZ_JPM_KEY_INT:
    case WLZ_JPM_KEY_SHORT:
    case WLZ_JPM_KEY_BYTE:
    case WLZ_JPM_KEY_FLOAT:
    case WLZ_JPM_KEY_DOUBLE:
      break;
    case WLZ_JPM_KEY_LONG_PTR1:
      aElmP = (void *)(*jEnv)->GetLongArrayElements(jEnv, (jlongArray )dstJObj,
      					            &isCpy);
      *(jlong *)aElmP = *(long *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseLongArrayElements(jEnv, (jlongArray )dstJObj,
					  (jlong *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_INT_PTR1:
      aElmP = (void *)(*jEnv)->GetIntArrayElements(jEnv, (jintArray )dstJObj,
      					           &isCpy);
      *(jint *)aElmP = *(int *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseIntArrayElements(jEnv, (jintArray )dstJObj,
					 (jint *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_SHORT_PTR1:
      aElmP = (void *)(*jEnv)->GetShortArrayElements(jEnv,
      						     (jshortArray )dstJObj,
      					             &isCpy);
      *(jshort *)aElmP = *(short *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseShortArrayElements(jEnv, (jshortArray )dstJObj,
					   (jshort *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_BYTE_PTR1:
      aElmP = (void *)(*jEnv)->GetByteArrayElements(jEnv,
      						    (jbyteArray )dstJObj,
						    &isCpy);
      *(jbyte *)aElmP = *(char *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseByteArrayElements(jEnv, (jbyteArray )dstJObj,
					  (jbyte *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_FLOAT_PTR1:
      aElmP = (void *)(*jEnv)->GetFloatArrayElements(jEnv,
      						     (jfloatArray )dstJObj,
						     &isCpy);
      *(jfloat *)aElmP = *(float *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseFloatArrayElements(jEnv, (jfloatArray )dstJObj,
					   (jfloat *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_DOUBLE_PTR1:
      aElmP = (void *)(*jEnv)->GetDoubleArrayElements(jEnv,
      						      (jdoubleArray )dstJObj,
      					              &isCpy);
      *(jdouble *)aElmP = *(double *)val;
      if(isCpy)
      {
        (*jEnv)->ReleaseDoubleArrayElements(jEnv, (jdoubleArray )dstJObj,
					    (jdouble *)aElmP, 0);
      }
      break;
    case WLZ_JPM_KEY_STRING_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "char", "String", "jstring", 1,
				 WLZ_JPM_KEY_STRING, (long )*(char **)val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_GREYV:
      newJObj = WlzJavaValueWrap(jEnv, "WlzGreyV", "WlzGreyV", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_GREYV, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_PIXELV:
      newJObj = WlzJavaValueWrap(jEnv, "WlzPixelV", "WlzPixelV", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_PIXELV, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_IVERTEX2_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzIVertex2", "WlzIVertex2", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_IVERTEX2, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX2_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzFVertex2", "WlzFVertex2", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_FVERTEX2, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX2_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzDVertex2", "WlzDVertex2", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_IVERTEX2, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_IVERTEX3_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzIVertex3", "WlzIVertex3", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_IVERTEX3, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_FVERTEX3_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzFVertex3", "WlzFVertex3", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_FVERTEX3, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_DVERTEX3_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzDVertex3", "WlzDVertex3", "jobject",
      				 1, WLZ_JPM_KEY_WLZ_DVERTEX3, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_IBOX2_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzIBox2", "WlzIBox2", "jobject",
      			         1, WLZ_JPM_KEY_WLZ_IBOX2, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_DBOX2_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzDBox2", "WlzDBox2", "jobject",
      			         1, WLZ_JPM_KEY_WLZ_DBOX2, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_IBOX3_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzIBox3", "WlzIBox3", "jobject",
      			         1, WLZ_JPM_KEY_WLZ_IBOX3, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_DBOX3_PTR1:
      newJObj = WlzJavaValueWrap(jEnv, "WlzDBox3", "WlzDBox3", "jobject",
      			         1, WLZ_JPM_KEY_WLZ_DBOX3, val);
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
      break;
    case WLZ_JPM_KEY_WLZ_PTR2:
      newJObj = WlzJavaValueWrap(jEnv, "WlzObject", "WlzObject", "jobject",
      				 2, WLZ_JPM_KEY_WLZ_PTR1,
				 (long )(*(WlzObject **)val));
      (*jEnv)->SetObjectArrayElement(jEnv, dstJObj, 0, newJObj);
    default:
      break;
  }
}

/************************************************************************
* Function:	WlzJavaValueFree					*
* Returns:	void							*
* Purpose:	Free's a temporary data structure.			*
* Global refs:	-							*
* Parameters:	void *dSP:		Data structure pointer.		*
*		int dSKey:		Data structure identification.	*
*		jboolean isCopy:	Coppy flag for JNI functions.	*
************************************************************************/
void		WlzJavaValueFree(void *dSP, int dSKey, jboolean isCopy)
{
  if(isCopy && dSP)
  {
    switch(dSKey)
    {
      case WLZ_JPM_KEY_WLZ_GREYV:    /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_PIXELV:   /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_IVERTEX2: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_FVERTEX2: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_DVERTEX2: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_IVERTEX3: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_FVERTEX3: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_DVERTEX3: /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_IBOX2:    /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_DBOX2:    /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_IBOX3:    /* FALLTHROUGH */
      case WLZ_JPM_KEY_WLZ_DBOX3:
        AlcFree(dSP);
        break;
      case WLZ_JPM_KEY_WLZ_PTR1:
	break;
      case WLZ_JPM_KEY_WLZ_FILE_PTR1:
	break;
      default:
	break;
    }
  }
}
