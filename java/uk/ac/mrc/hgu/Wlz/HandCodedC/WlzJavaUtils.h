/************************************************************************
* Project:      Java Woolz
* Title:        WlzJavaUtils.h
* Date:         January 1999
* Purpose:      Prototypes of misc functions for the C side of Java
*		Woolz.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/

extern char			*WlzJavaBuildFQClassName(
				  const char *name);
extern void			WlzJavaThrowWlzException(
				  JNIEnv *jEnv,
				  WlzErrorNum errNum);
