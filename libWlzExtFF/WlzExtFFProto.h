#ifndef WLZEXTFFPROTO_H
#define WLZEXTFFPROTO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzExtFFProto.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file with function prototypes for external data
*		file format support for the MRC Human Genetics Unit
*		Woolz library.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */


/* From WlzExtFF.c */
WlzEffFormat	WlzEffStringToFormat(const char *fmtStr);
WlzEffFormat    WlzEffStringExtToFormat(const char *extStr);
const char	*WlzEffStringFromFormat(WlzEffFormat fileFmt,
				const char **dstExtStr);
WlzObject	*WlzEffReadObj(FILE *fP, const char *fName, WlzEffFormat fFmt,
			       WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObj(FILE *fP, const char *fName, WlzObject *obj,
			       WlzEffFormat fFmt);

/* From WlzExtFFBmp.c */
WlzObject 	*WlzEffReadObjBmp(const char *gvnFileName,
				WlzErrorNum *dstErr);
WlzErrorNum	WlzEffReadObjBmpData2D(FILE *fP, WlzIVertex2 *imgSz,
				unsigned char ***data),
		WlzEffWriteObjBmp2D(const char *fNameStr, WlzObject *obj,
				WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
				unsigned char *data, unsigned char bgd),
		WlzEffWriteObjBmp(const char *gvnFileName, 
				WlzObject *obj);
/* From WlzExtFFDen.c */
WlzObject	*WlzEffReadObjDen(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjDen(FILE *fP, WlzObject *obj);

/* From WlzExtFFIcs.c */
WlzErrorNum	WlzEffWriteObjIcs(const char *gvnFileName, WlzObject *obj);
WlzErrorNum	WlzEffIcsFileNames(char **fileBody,
				char **icsFileName,
				char **idsFileName,
				const char *gvnFileName);
WlzObject	*WlzEffReadObjIcs(const char *gvnFileName,
				WlzErrorNum *dstErr);

/* From WlzExtFFPic.c */
WlzObject 	*WlzEffReadObjPic(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjPic(FILE *fP, WlzObject *obj);

/* From WlzExtFFPnm.c */
WlzObject 	*WlzEffReadObjPnm(const char *gvnFileName,
				WlzErrorNum *dstErr);
WlzErrorNum	WlzEffReadObjPnmData2D(FILE *fP, WlzIVertex2 *imgSz,
				unsigned char ***data),
		WlzEffWriteObjPnm2D(const char *fNameStr, WlzObject *obj,
				WlzIVertex2 imgSz, WlzIVertex2 imgOrg,
				unsigned char *data, unsigned char bgd),
		WlzEffWriteObjPnm(const char *gvnFileName, 
				WlzObject *obj);

/* From WlzExtFFSlc.c */
WlzObject	*WlzEffReadObjSlc(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjSlc(FILE *fP, WlzObject *obj);

/* From WlzExtFFVff.c */
WlzObject	*WlzEffReadObjVff(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjVff(FILE *fP, WlzObject *obj);

/* From WlzExtFFVtk.c */
WlzObject	*WlzEffReadObjVtk(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjVtk(FILE *fP, WlzObject *obj);

/* From WlzExtFFSlc.c */
WlzObject	*WlzEffReadObjSlc(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjSlc(FILE *fP, WlzObject *obj);

/* From WlzExtFFStack.c */
WlzObject	*WlzEffReadObjStack(const char *gvnFileName, WlzEffFormat fFmt,
				WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjStack(const char *gvnFileName, WlzEffFormat fFmt,
				WlzObject *obj);

/* From WlzExtFFIPL.c */
WlzObject	*WlzEffReadObjIPL(FILE *fP, WlzErrorNum *dstErr);
WlzErrorNum	WlzEffWriteObjIPL(FILE *fP, WlzObject *obj);

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* ! WLZEXTFFPROTO_H */
