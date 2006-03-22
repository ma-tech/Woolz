#ifndef WLZEXTFF_PROTO_H
#define WLZEXTFF_PROTO_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzExtFFProto_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlzExtFF/WlzExtFFProto.h
* \author       Bill Hill
* \date         March 1999
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Header file with function prototypes for external data file
*		format support for the MRC Human Genetics Unit Woolz library.
* \ingroup	WlzExtFF
* \todo         -
* \bug          None known.
*/

#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */


/* From WlzExtFF.c */
extern WlzEffFormat 		WlzEffStringToFormat(
				  const char *fmtStr);
extern WlzEffFormat 		WlzEffStringExtToFormat(
				  const char *extStr);
extern const char 		*WlzEffStringFromFormat(
				  WlzEffFormat fileFmt,
				  const char **dstExtStr);
extern WlzObject 		*WlzEffReadObj(
				  FILE *fP,
				  const char *fName,
				  WlzEffFormat fFmt,
				  int split,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObj(
				  FILE *fP,
				  const char *fName,
				  WlzObject *obj,
				  WlzEffFormat fFmt);
extern int			WlzEffNumberOfFormats(void);

#ifndef WLZ_EXT_BIND
/* From WlzExtFFAm.c */
extern WlzObject		*WlzEffReadObjAm(
				  FILE *fP,
				  int split,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzEffWriteObjAm(
				  FILE *fP,
				  WlzObject *obj);
/* From WlzExtFFAnl.c */
extern WlzErrorNum		WlzEffAnlFileNames(
				  char **fileBody,
				  char **hdrFileName,
				  char **imgFileName,
				  const char *gvnFileName);
extern WlzObject		*WlzEffReadObjAnl(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjAnl(
				  const char *gvnFileName,
				  WlzObject *obj);
/* From WlzExtFFBmp.c */
extern WlzObject 		*WlzEffReadObjBmp(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffReadObjBmpData2D(
				  FILE *fP,
				  WlzIVertex2 *imgSz,
				  unsigned char ***data);
extern WlzErrorNum 		WlzEffWriteObjBmp2D(
				  const char *fNameStr,
				  WlzObject *obj,
				  WlzIVertex2 imgSz,
				  WlzIVertex2 imgOrg,
				  unsigned char *data,
				  unsigned char bgd);
extern WlzErrorNum 		WlzEffWriteObjBmp(
				  const char *gvnFileName,
				  WlzObject *obj);
/* From WlzExtFFDen.c */
extern WlzObject 		*WlzEffReadObjDen(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjDen(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFIcs.c */
extern WlzErrorNum 		WlzEffWriteObjIcs(
				  const char *gvnFileName,
				  WlzObject *obj);
extern WlzErrorNum 		WlzEffIcsFileNames(
				  char **fileBody,
				  char **icsFileName,
				  char **idsFileName,
				  const char *gvnFileName);
extern WlzObject 		*WlzEffReadObjIcs(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);

/* From WlzExtFFPic.c */
extern WlzObject 		*WlzEffReadObjPic(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjPic(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFPnm.c */
extern WlzObject 		*WlzEffReadObjPnm(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffReadObjPnmData2D(
				  FILE *fP,
				  WlzIVertex2 *imgSz,
				  unsigned char ***data);
extern WlzErrorNum 		WlzEffWriteObjPnm2D(
				  const char *fNameStr,
				  WlzObject *obj,
				  WlzIVertex2 imgSz,
				  WlzIVertex2 imgOrg,
				  unsigned char *data,
				  unsigned char bgd);
extern WlzErrorNum 		WlzEffWriteObjPnm(
				  const char *gvnFileName,
				  WlzObject *obj);

/* From WlzExtFFSlc.c */
extern WlzObject 		*WlzEffReadObjSlc(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjSlc(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFVff.c */
extern WlzObject 		*WlzEffReadObjVff(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjVff(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFVtk.c */
extern WlzObject 		*WlzEffReadObjVtk(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjVtk(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFSlc.c */
extern WlzObject 		*WlzEffReadObjSlc(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjSlc(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFStack.c */
extern WlzObject 		*WlzEffReadObjStack(
				  const char *gvnFileName,
				  WlzEffFormat fFmt,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjStack(
				  const char *gvnFileName,
				  WlzEffFormat fFmt,
				  WlzObject *obj);

/* From WlzExtFFIPL.c */
extern WlzObject 		*WlzEffReadObjIPL(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjIPL(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFTiff.c */
extern WlzErrorNum 		WlzEffWriteObjTiff(
				  const char *tiffFileName,
				  WlzObject *obj);
extern WlzObject 		*WlzEffReadObjTiff(
				  const char *tiffFileName,
				  WlzErrorNum *dstErr);

/* From WlzExtFFJpeg.c */
extern WlzObject		*WlzEffReadObjJpeg(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzEffWriteObjJpeg(
				  FILE *fP,
				  WlzObject *obj,
				  char	*params);

/* From	WlzExtFFBibUtils.c */
extern WlzErrorNum 		WlzEffBibWrite3DSectionViewParamsRecord(
				  FILE *fp,
				  char *recordName,
				  WlzThreeDViewStruct *wlzViewStr);
extern WlzErrorNum 		WlzEffBibParse3DSectionViewParamsRecord(
				  BibFileRecord *bibfileRecord,
				  WlzThreeDViewStruct *wlzViewStr);
extern WlzErrorNum 		WlzEffBibWriteWarpTransformParamsRecord(
				  FILE *fp,
				  char *recordName,
				  WlzFnType basisFnType,
				  WlzTransformType affineType,
				  WlzMeshGenMethod meshMthd,
				  int meshMinDst,
				  int meshMaxDst);
extern WlzErrorNum 		WlzEffBibParseWarpTransformParamsRecord(
				  BibFileRecord *bibfileRecord,
				  WlzFnType *basisFnType,
				  WlzTransformType *affineType,
				  WlzMeshGenMethod *meshMthd,
				  int *meshMinDst,
				  int *meshMaxDst);
extern WlzErrorNum 		WlzEffBibWriteTiePointVtxsRecord(
				  FILE *fp,
				  char *recordName,
				  int index,
				  WlzDVertex3 dstVtx,
				  WlzDVertex3 srcVtx,
				  int relativeFlg);
extern WlzErrorNum 		WlzEffBibParseTiePointVtxsRecord(
				  BibFileRecord *bibfileRecord,
				  int *index,
				  WlzDVertex3 *dstVtx,
				  WlzDVertex3 *srcVtx,
				  int *relativeFlg);
extern WlzErrorNum 		WlzEffBibWriteFileRecord(
				  FILE *fp,
				  char *recordName,
				  char *fileName,
				  WlzEffFormat fileType);
extern WlzErrorNum 		WlzEffBibParseFileRecord(
				  BibFileRecord *bibfileRecord,
				  int *dstIndex,
				  char **dstFileName,
				  WlzEffFormat *dstFileType);
extern WlzErrorNum 		WlzEffBibWriteWarpInputSegmentationParamsRecord(
				  FILE *fp,
				  char *recordName,
				  int normFlg,
				  int histoFlg,
				  int shadeFlg,
				  int gaussFlg,
				  double width);
extern WlzErrorNum 		WlzEffBibParseWarpInputSegmentationParamsRecord(
				  BibFileRecord *bibfileRecord,
				  int *normFlg,
				  int *histoFlg,
				  int *shadeFlg,
				  int *gaussFlg,
				  double *width,
				  int *threshLow,
				  int *threshHigh);
extern WlzErrorNum 		WlzEffBibWriteWarpInputThresholdParamsRecord(
				  FILE *fp,
				  char *recordName,
			WlzEffBibWarpInputThresholdParamsStruct *paramStruct);

extern WlzErrorNum		WlzEffBibParseWarpInputThresholdParamsRecord(
  				  BibFileRecord *bibfileRecord,
  			WlzEffBibWarpInputThresholdParamsStruct *paramStruct);
#endif /* WLZ_EXT_BIND */

#ifdef  __cplusplus
}
#endif /* __cplusplus */

#endif /* ! WLZEXTFF_PROTO_H */
