#ifndef RECONSTRUCT_PROTO_H
#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructProto.h				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Header file with function prototypes for the MRC
*		Human Genetics Unit reconstruction library.	
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#define RECONSTRUCT_PROTO_H

#ifdef  __cplusplus
extern "C" {
#endif

/* From ReconstructAuto.c */
extern RecError RecAuto(RecControl *, RecPPControl *, HGUDlpList *, int *,
		        RecSecUpdateFunction, void *,
			RecWorkFunction, void *,
			char **);

/* From ReconstructConstruct3D.c */

extern RecError	RecConstruct3DObj(WlzObject **, HGUDlpList *, double,
			       int, WlzInterpolationType, int, int,
			       int, WlzDVertex3, int, int,
			       WlzDBox3 *, WlzDBox3 *, char **);

/* From ReconstructCrossCor.c */
extern RecError	RecCrossCorrelate(double **, double **, double *, double *,
				  RecCcFlag, WlzObject *, WlzObject *,
				  WlzIVertex2, WlzIVertex2, RecPPControl *),
		RecCrossCorrelateROI(double *, double ***, double ***,
				     WlzObject *, WlzObject *,
				     WlzIVertex2, WlzIVertex2,
				     RecPPControl *);
/* From ReconstructDebug.c */
extern RecDbgMask recDbgMask,
		recDbgWlzMask;
extern void	*recDbgData,
		*recDbgWlzData;

extern RecError	RecDbgWrite(char *, ...),
		RecDbgWlzWrite(WlzObject *, int);

/* From ReconstructExplode3D.c */
extern RecError	RecExplode3DObjToFile(char *, char *, WlzObject *,
				      char *, WlzEffFormat, char **),
		RecExplode3DFileToFile(char *, char *, char *,
				       WlzEffFormat, char **);
/* From ReconstructFileIO.c */
extern RecError	RecFileSecListWrite(FILE *, RecSectionList *, int, char **),
		RecFileSecListRead(RecSectionList *, int *,
				   FILE *, char **),
		RecFileSecObjRead(RecSection *, char **),
		RecFileSecObjRead(RecSection *, char **),
		RecFileSecObjsRead(HGUDlpList *, int, int, int, char **),
		RecFileObjWlzRead(FILE *, WlzObject **),
		RecFileObjWlzWrite(FILE *, WlzObject *);
extern void	RecFileSecObjFree(RecSection *),
		RecFileSecObjsFree(HGUDlpList *, int, int, int);

/* From ReconstructMisc.c */
extern const char *RecErrorToStr(RecError),
		*RecMethodToStr(RecMethod);
extern int	RecPowerOfTwo(int *ipLarger, unsigned int given);
extern RecError	RecErrorFromWlz(WlzErrorNum);
extern WlzIVertex2 RecPowerOfTwoC2I(WlzIVertex2 *ipLarger, WlzIVertex2 given);

/* From ReconstructPreProc.c */
extern WlzObject *RecPreProcObj(WlzObject *, RecPPControl *, RecError *);

/* From ReconstructRegister.c */
extern RecError	RecRegisterPair(WlzAffineTransform **,
				double *, int *,
				RecControl *, RecPPControl *,
			        WlzObject *, WlzObject *,
				RecWorkFunction, void *, char **),
		RecRegisterTiePoints(WlzAffineTransform **, double *, double *,
				     RecTiePointPair *, int,
				     WlzObject *, WlzObject *);

/* From ReconstructRotMatch.c */
extern RecError RecRotMatch(double *angle, double *value,
			    WlzObject *obj0, WlzObject *obj1,
			    WlzIVertex2 cRot, double angleInc, double distInc,
			    int maxRadiusFlag, RecPPControl *ppCtrl);

/* From ReconstructSection.c */
extern int	RecSecIsEmpty(RecSection *);
extern char	*RecSecToStr(RecSection *, unsigned int, char **);
extern void	RecSecFree(RecSection *),
		RecSecListSort(HGUDlpList *, unsigned int),
		RecSecListIndiciesSet(HGUDlpList *, int, int),
		RecSecRecSetDefaults(RecReconstruction *);
extern unsigned int RecSecListInvalid(HGUDlpList *, HGUDlpListItem **,
				    unsigned int);
extern RecSection *RecSecAssign(RecSection *),
		*RecSecDup(RecSection *),
		*RecSecMake(int, int, double, char *,
			    WlzAffineTransform *, WlzObject *),
		*RecSecMakeEmpty(int),
		*RecSecNext(HGUDlpList *, HGUDlpListItem *,
			    HGUDlpListItem **, int),
		*RecSecPrev(HGUDlpList *, HGUDlpListItem *,
			    HGUDlpListItem **, int);
extern RecError	RecSecListToStrList(char ***, HGUDlpList *, int, char **,
                                    unsigned int),
		RecSecAppendListFromFiles(HGUDlpList *, HGUDlpListItem *,
					  char **, int,
					  int, int),
		RecSecCumTransfSet(HGUDlpList *, HGUDlpListItem *),
		RecSecCumTransfClear(HGUDlpList *, HGUDlpListItem *);
HGUDlpListItem  *RecSecFindItemIndex(HGUDlpList *, HGUDlpListItem *,
				     int, HGUDlpListDirection);
RecSectionList  *RecSecNewSectionList(RecError *);

/* From ReconstructTranMatch.c */
extern RecError	RecTranMatch(WlzDVertex2 *shift, double *value,
			     WlzObject *obj0, WlzObject *obj1,
			     WlzIVertex2 maxShift, RecPPControl *ppCtrl);

extern void	RecTranFindPeak(WlzDVertex2 *shift, double *pValue,
                                double **data, WlzIVertex2 size,
				WlzIVertex2 search);

#ifdef  __cplusplus
}
#endif

#endif /* RECONSTRUCT_PROTO_H */
