#ifndef RECONSTRUCT_PROTO_H
#define RECONSTRUCT_PROTO_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _ReconstructProto_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libReconstruct/ReconstructProto.h
* \author       Bill Hill
* \date         April 1999
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
* \brief	Header file with function prototypes for the Reconstruct
* 		library.
* \ingroup	Reconstruct
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif
#endif /* WLZ_EXT_BIND */

/* From ReconstructAuto.c */
extern RecError			RecAuto(
				  RecControl *rCtrl,
				  RecPPControl *ppCtrl,
				  HGUDlpList *secList,
				  int *cancelFlag,
				  RecSecUpdateFunction secFn,
				  void *secData,
				  RecWorkFunction workFn,
				  void *workData,
				  char **eMsg);

/* From ReconstructConstruct3D.c */
extern RecError			RecConstruct3DObj(
				  WlzObject **dstObj,
				  HGUDlpList *secList,
				  double confLimit,
				  int gaussSamFlg,
				  WlzInterpolationType interp,
				  int fastSamFlg,
				  int greedyFlg,
				  int intScaleFlg,
				  WlzDVertex3 scale,
				  int matchHistFlg,
				  int matchIdx,
				  WlzDBox3 *srcReg,
				  WlzDBox3 *dstReg,
				  char **eMsg);

/* From ReconstructCrossCor.c */
extern RecError			RecCrossCorrelate(
				  double **data0,
				  double **data1,
				  double *sSq0,
				  double *sSq1,
				  RecCcFlag ccFlags,
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzIVertex2 org,
				  WlzIVertex2 size,
				  RecPPControl *ppCtrl);
extern RecError			RecCrossCorrelateROI(
				  double *newCC,
				  double ***data0,
				  double ***data1,
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzIVertex2 roiCtr,
				  WlzIVertex2 roiSz,
				  RecPPControl *ppCtrl);

/* From ReconstructDebug.c */
extern RecDbgMask 		recDbgMask;
extern RecDbgMask		recDbgWlzMask;
extern void			*recDbgData;
extern void			*recDbgWlzData;

extern RecError			RecDbgWrite(
				  char *fmt,
				  ...);
extern RecError			RecDbgWlzWrite(
				  WlzObject *obj, int freeFlg);

/* From ReconstructExplode3D.c */
extern RecError			RecExplode3DObjToFile(
				  char *dstDirStr,
				  char *dstBodyStr,
				  WlzObject *srcObj,
				  char *srcFName,
				  WlzEffFormat srcFFormat,
				  char **eMsg);
extern RecError			RecExplode3DFileToFile(
				  char *dstDirStr,
				  char *dstBodyStr,
				  char *srcFile,
				  WlzEffFormat srcFmt,
				  char **eMsg);

/* From ReconstructFileIO.c */
extern RecError			RecFileSecListWrite(
				  FILE *fP,
				  RecSectionList *secList,
				  int numSec,
				  char **eMsg);
extern RecError			RecFileSecListRead(
				  RecSectionList *secList,
				  int *numSec,
				  FILE *fP,
				  char **eMsg);
extern RecError			RecFileSecObjRead(
				  RecSection *sec,
				  char **eMsg);
extern RecError			RecFileSecObjRead(
				  RecSection *sec,
				  char **eMsg);
extern RecError			RecFileSecObjsRead(
				  HGUDlpList *secList,
				  int minIdx,
				  int maxIdx,
				  int allSecFlag,
				  char **eMsg);
extern RecError			RecFileObjWlzRead(
				  FILE *fP,
				  WlzObject **obj);
extern RecError			RecFileObjWlzWrite(
				  FILE *fP,
				  WlzObject *obj);
extern void			RecFileSecObjFree(
				  RecSection *sec);
extern void			RecFileSecObjsFree(
				  HGUDlpList *secList,
				  int minIdx,
				  int maxIdx,
				  int allSecFlag);

/* From ReconstructMisc.c */
extern const char		*RecErrorToStr(
				  RecError errFlag);
extern const char		*RecMethodToStr(
				  RecMethod method);
extern int			RecPowerOfTwo(
				  int *ipLarger,
				  unsigned int given);
extern RecError			RecErrorFromWlz(
				  WlzErrorNum wlzErr);
extern WlzIVertex2		RecPowerOfTwoC2I(
				  WlzIVertex2 *ipLarger,
				  WlzIVertex2 given);

/* From ReconstructPreProc.c */
extern WlzObject		*RecPreProcObj(
				  WlzObject *obj,
				  RecPPControl *ppCtrl,
				  RecError *dstErr);

/* From ReconstructRegister.c */
extern RecError			RecRegisterPair(
				  WlzAffineTransform **dstTrans,
				  double *dstCrossC,
				  int *dstIter,
				  RecControl *rCtrl,
				  RecPPControl *ppCtrl,
				  WlzObject *obj0,
				  WlzObject *obj1,
				  RecWorkFunction workFn,
				  void *workData,
				  char **eMsg);
extern RecError			RecRegisterTiePoints(
				  WlzAffineTransform **dstTr,
				  double *dstED,
				  double *dstCC,
				  RecTiePointPair *tppVec,
				  int tppCount,
				  WlzObject *obj0,
				  WlzObject *obj1);

/* From ReconstructRotMatch.c */
extern RecError 		RecRotMatch(
				  double *angle,
				  double *value,
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzIVertex2 cRot,
				  double angleInc,
				  double distInc,
				  int maxRadiusFlag,
				  RecPPControl *ppCtrl);

/* From ReconstructSection.c */
extern int			RecSecIsEmpty(
				  RecSection *sec);
extern char			*RecSecToStr(
				  RecSection *sec,
				  unsigned int reqFields,
				  char **eMsg);
extern void			RecSecFree(
				  RecSection *sec);
extern void			RecSecListSort(
				  HGUDlpList *secList,
				  unsigned int sortMask);
extern void			RecSecListIndiciesSet(
				  HGUDlpList *secList,
				  int indexFirst,
				  int indexInc);
extern void			RecSecRecSetDefaults(
				  RecReconstruction *rec);
extern unsigned int	 	RecSecListInvalid(
				  HGUDlpList *secList,
				  HGUDlpListItem **invalidItem,
				  unsigned int mask);
extern RecSection		*RecSecAssign(
				  RecSection *sec);
extern RecSection		*RecSecDup(
				  RecSection *sec);
extern RecSection		*RecSecMake(
				  int index,
				  int iterations,
				  double correlation,
				  char *imageFile,
				  WlzAffineTransform *transform,
				  WlzObject *obj);
extern RecSection		*RecSecMakeEmpty(
				  int index);
extern RecSection		*RecSecNext(
				  HGUDlpList *secList,
				  HGUDlpListItem *cItem,
				  HGUDlpListItem **nItemP,
				  int skipEmpty);
extern RecSection		*RecSecPrev(
				  HGUDlpList *secList,
				  HGUDlpListItem *cItem,
				  HGUDlpListItem **pItemP,
				  int skipEmpty);
extern RecError			RecSecListToStrList(
				  char ***strList,
				  HGUDlpList *secList,
				  int numSec,
				  char **eMsg,
				  unsigned int reqFields);
extern RecError	 		RecSecAppendListFromFiles(
				  HGUDlpList *secList,
				  HGUDlpListItem *secMark,
				  char **fileVec,
				  int nFiles,
				  int indexFirst,
				  int indexInc);
extern RecError	 		RecSecCumTransfSet(
				  HGUDlpList *secList,
				  HGUDlpListItem *secItem);
extern RecError	 		RecSecCumTransfClear(
				  HGUDlpList *secList,
				  HGUDlpListItem *secItem);
extern HGUDlpListItem		*RecSecFindItemIndex(
				  HGUDlpList *list,
				  HGUDlpListItem *from,
				  int index,
				  HGUDlpListDirection dir);
extern RecSectionList		*RecSecNewSectionList(
				  RecError *dstErr);


/* From ReconstructTranMatch.c */
extern RecError			RecTranMatch(
				  WlzDVertex2 *shift,
				  double *value,
				  WlzObject *obj0,
				  WlzObject *obj1,
				  WlzIVertex2 maxShift,
				  RecPPControl *ppCtrl);

extern void			RecTranFindPeak(
				  WlzDVertex2 *shift,
				  double *pValue,
				  double **data,
				  WlzIVertex2 size,
				  WlzIVertex2 search);

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif
#endif /* WLZ_EXT_BIND */

#endif /* RECONSTRUCT_PROTO_H */
