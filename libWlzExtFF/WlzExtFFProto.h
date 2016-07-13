#ifndef WLZEXTFF_PROTO_H
#define WLZEXTFF_PROTO_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFProto_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzExtFF/WlzExtFFProto.h
* \author       Bill Hill
* \date         March 1999
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
* \brief	Header file with function prototypes for external data file
*		format support for the MRC Human Genetics Unit Woolz library.
* \ingroup	WlzExtFF
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */


/* From WlzExtFF.c */
extern WlzEffFormat 		WlzEffStringToFormat(
				  const char *fmtStr);
extern WlzEffFormat 		WlzEffStringExtToFormat(
				  const char *extStr);
extern WlzEffFormat 		WlzEffStringFormatFromFileName(
				  const char *fNameStr);
extern const char 		*WlzEffStringFromFormat(
				  WlzEffFormat fileFmt,
				  const char **dstExtStr);
extern WlzObject 		*WlzEffReadObj(
				  FILE *fP,
				  const char *fName,
				  WlzEffFormat fFmt,
				  int split,
				  int sTrans,
				  int gTrans,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObj(
				  FILE *fP,
				  const char *fName,
				  WlzObject *obj,
				  WlzEffFormat fFmt);
extern int			WlzEffNumberOfFormats(void);
extern char			*WlzEffFormatTable(
				  WlzUInt indWth,
                                  WlzUInt desWth,
				  WlzUInt fmtWth,
				  WlzErrorNum *dstErr);

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
/* From WlzExtFFGif.c */
extern WlzObject		*WlzEffReadObjGif(
				  FILE *fP,
				  WlzErrorNum *dstErr);

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

/* From WlzExtFFNifti.c */
extern WlzObject 		*WlzEffReadObjNifti(
				  const char *gvnFileName,
				  int spatialTr,
				  int greySc,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjNifti(
				  const char *gvnFileName,
				  WlzObject *obj);

/* From WlzExtFFNodeEle.c */
extern WlzErrorNum 		WlzEffNodeEleFileNames(
				  char **fileBody,
				  char **nodeFileName,
				  char **eleFileName,
				  const char *gvnFileName);
extern WlzObject 		*WlzEffReadObjNodeEle(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjNodeEle(
				  const char *gvnFileName,
				  WlzObject *obj);

/* From WlzExtFFNrrd.c */
extern WlzObject		*WlzEffReadObjNrrd(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzEffWriteObjNrrd(
				  FILE *fP,
				  WlzObject *obj);
/* From WlzExtFFMesh.c */
extern WlzObject 		*WlzEffReadObjMesh(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjMesh(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFObj.c */
extern WlzObject 		*WlzEffReadObjObj(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjObj(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFEMT.c */
extern WlzObject 		*WlzEffReadObjEMT(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjEMT(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFPic.c */
extern WlzObject 		*WlzEffReadObjPic(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjPic(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFPly2.c */
extern WlzObject 		*WlzEffReadObjPly2(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjPly2(
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

/* From WlzExtFFPvl.c */
extern WlzObject		*WlzEffReadObjPvl(
				  const char *gvnFileName,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzEffWriteObjPvl(
				  const char *gvnFileName,
				  WlzObject *obj);

/* From WlzExtFFSlc.c */
extern WlzObject 		*WlzEffReadObjSlc(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjSlc(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFSMesh.c */
extern WlzObject 		*WlzEffReadObjSMesh(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjSMesh(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFVff.c */
extern WlzObject 		*WlzEffReadObjVff(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjVff(
				  FILE *fP,
				  WlzObject *obj);

/* From WlzExtFFVMesh.c */
extern WlzObject 		*WlzEffReadObjVMesh(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum 		WlzEffWriteObjVMesh(
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

/* From WlzExtFFStl.c */
extern WlzObject       		*WlzEffReadObjStl(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum     		WlzEffWriteObjStl(
				  FILE *fP,
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
				  int	split,
				  WlzErrorNum *dstErr);

/* From WlzExtFFJpeg.c */
extern WlzObject		*WlzEffReadObjJpeg(
				  FILE *fP,
				  WlzErrorNum *dstErr);
extern WlzErrorNum		WlzEffWriteObjJpeg(
				  FILE *fP,
				  WlzObject *obj,
				  char	*params);

/* From WlzExtFFTxt.c */
extern WlzErrorNum		WlzEffWriteObjTxt(
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
extern WlzThreeDViewStruct 	*WlzEffBibRead3DView(
				  FILE *fP,
				  char **dstMsg,
				  WlzErrorNum *dstErr);

#endif /* WLZ_EXT_BIND */

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif /* __cplusplus */
#endif /* WLZ_EXT_BIND */

#endif /* ! WLZEXTFF_PROTO_H */
