#ifndef ALGPROTO_H
#define ALGPROTO_H
#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgProto.h
* Date:         September 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file with function prototypes for the MRC
*		Human Genetics Unit numerical algorithm library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
* 07-11-00 bill Add AlgShuffleIdx().
* 08-08-00 bill Add AlgBitNextSet().
* 12-07-00 bill	Add AlgHeapSort(), AlgHeapSortIdx() and AlgHeapElmSwap().
* 11-05-00 bill	Add AlgGammaLog(), AlgGammaP(), AlgLinearFit1D(),
*		AlgLinearFitIdx1D(), AlgRange1D(), AlgRangeIdx1D()
*		AlgBitSetCount() and AlgBitSetPositions().
* 31-03-00 bill Add AlgMixtureSyn() and tidy up prototypes.
* 26-01-00 bill Add AlgConvolve() and AlgMixture().
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif

/* From AlgBits.c */
extern int		AlgBitSetCount(
			  unsigned long gMsk);
extern int		AlgBitSetPositions(
			  int *posA,
			  unsigned long gMsk);
extern unsigned long	AlgBitNextNOfM(
			  unsigned long curMsk,
			  int n,
			  int m);
extern int		AlgBitNextSet(
			  unsigned long msk,
			  int idC);
/* From AlgComplexUtils.c */
extern double		AlgCMod(
			  ComplexD z);
extern double		AlgCArg(
			  ComplexD z);
extern double		AlgCRe(
			  ComplexD z);
extern double		AlgCIm(
			  ComplexD z);
extern ComplexD		AlgCConj(
			  ComplexD z);
extern ComplexD 	AlgCAdd(
			  ComplexD z1,
			  ComplexD z2);
extern ComplexD 	AlgCSub(
			  ComplexD z1,
			  ComplexD z2);
extern ComplexD 	AlgCMult(
			  ComplexD z1,
			  ComplexD z2);
extern ComplexD 	AlgCDiv(
			  ComplexD z1,
			  ComplexD z2);
extern ComplexD 	AlgCPow(
			  ComplexD z,
			  double pow);

/* From AlgConvolve.c */
extern AlgError 	AlgConvolve(
			  int sizeArrayCnv,
			  double *arrayCnv,
			  int sizeArrayKrn,
			  double *arrayKrn,
    	   	   	  int sizeArrayDat, 
			  double *arrayDat,
    			  AlgPadType pad);

/* From AlgDebug.c */
extern AlgError		AlgDbgWrite(
			  char *fmt, ...);

/* From AlgFourier.c */
extern void     	AlgFourHart1D(
			  double *data,
			  int num, 
			  int step,
			  int cThr);
extern void		AlgFour1D(
			  double *real,
			  double *imag,
			  int num,
			  int step,
			  int cThr);
extern void		AlgFourInv1D(
			  double *real,
			  double *imag,
			  int num,
			  int step,
			  int cThr);
extern void		AlgFourReal1D(
			  double *real,
			  int num,
			  int step,
			  int cThr);
extern void		AlgFourRealInv1D(
			  double *real,
			  int num,
			  int step,
			  int cThr);
extern void		AlgFour2D(
			  double **real,
			  double **imag,
			  double *reBuf,
			  double *imBuf,
			  int numX,
			  int numY,
			  int cThr);
extern void		AlgFourInv2D(
			  double **real,
			  double **imag,
			  double *reBuf,
			  double *imBuf,
			  int numX,
			  int numY,
			 int cThr);
extern void		AlgFourReal2D(
			  double **real,
			  double *reBuf,
			  double *imBuf,
			  int numX,
			  int numY,
			  int cThr);
extern void		AlgFourRealInv2D(
			  double **real,
			  double *reBuf,
			  double *imBuf,
			  int numX,
			  int numY,
			  int cThr);

/* From AlgGamma.c */
extern double		AlgGammaLog(
			  double x,
			  AlgError *dstErr);
extern double		AlgGammaP(
			  double a,
			  double x,
			  AlgError *dstErr);

/* From AlgGamma.c */
extern AlgError		AlgLinearFit1D(
			  int datSz,
			  double *datXA,
			  double *datYA,
			  double *dstA,
			  double *dstB,
			  double *dstSigA,
			  double *dstSigB,
			  double *dstQ);
extern AlgError		AlgLinearFitIdx1D(
			  double *datXA,
			  double *datYA,
			  int *idxXA,
			  int *idxYA,
			  int idxASz,
			  double *dstA,
			  double *dstB,
			  double *dstSigA,
			  double *dstSigB,
			  double *dstQ);

/* From AlgHeapSort.c */
extern AlgError		AlgHeapSort(
			  void *data,
			  unsigned nElm,
			  unsigned elmSz,
			  int (*cmpFn)(void *, void *));
extern AlgError		AlgHeapSortIdx(
			  void *data,
			  int *idx,
			  unsigned nElm,
			  int (*cmpFn)(void *, int *, int, int));
extern void		AlgHeapElmSwap(
			  void *elm0,
			  void *elm1,
			  int cnt);

/* From AlgMatrixGauss.c */
extern AlgError		AlgMatrixGaussSolve(
			  double **abMat,
			  int aSz,
			  double *xMat);

/* From AlgMatrixLU.c */
extern AlgError		AlgMatrixLUSolve(
			  double **aMat,
			  int aSz,
			  double *bMat,
			  int bSz);
extern AlgError		AlgMatrixLUInvert(
			  double **aMat,
			  int aSz);
extern AlgError		AlgMatrixLUDeterm(
			  double **aMat,
			  int aSz,
			  double *determ);
extern AlgError		AlgMatrixLUDecomp(
			  double **aMat,
			  int aSz,
			  int *idxVec,
			  double *evenOdd);
extern AlgError		AlgMatrixLUBackSub(
			  double **aMat,
			  int aSz,
			  int *idxVec,
			  double *bMat);

/* From AlgMatrixSV.c */
extern AlgError		AlgMatrixSVSolve(
			  double **aMat,
			  int nM,
			  int nN,
			  double *bMat,
			  double tol);
extern AlgError		AlgMatrixSVDecomp(double **aMat,
			  int nM,
			  int nN,
			  double *wMat,
			  double **vMat);
extern AlgError		AlgMatrixSVBackSub(double **uMat,
			  int nM,
			  int nN,
			  double *wMat,
			  double **vMat,
			  double *bMat);

/* From AlgMixture.c */
AlgError		AlgMixtureMLG(
			  int nDbn,
			  int nCls,
			  double samOrg,
			  double samItv,
			  double *freq,
			  double *alpha,
			  double *mean,
			  double *sd,
			  double tol,
			  double sumFreq,
			  double *dstLL,
			  int *dstNItn);
AlgError		AlgMixtureSyn(
			  int nCls,
			  int *synFreq,
			  int nObv,
			  double synOrg,
			  double synStep,
			  int nDbn,
			  double *alpha,
			  double *mu,
			  double *sigma);

/* From AlgPolyLSQ.c */
extern AlgError        	AlgPolynomialLSq(
			  double *xVec,
			  double *yVec,
			  int vecSz,
			  int polyDeg,
			  double *cVec);

/* From AlgRand.c */
extern void		AlgRandSeed(
			  long seed);
extern double		AlgRandUniform(
			  void);
extern double		AlgRandNormal(
			  double mu,
			  double sigma);
/* From AlgRange.c */
extern AlgError		AlgRange1D(
			  int datASz,
			  double *datA,
			  double *dstMin,
			  double *dstMax);
extern AlgError		AlgRangeIdx1D(
			  double *datA,
			  int idxASz,
			  int *idxA,
			  double *dstMin,
			  double *dstMax);

/* From AlgShuffle.c */
extern AlgError		AlgShuffleIdx(
			  int nShuffle,
			  int *shuffle,
			  int seed);
/* Debugging */
extern AlgDbgMask 	algDbgMask;
extern void		*algDbgData;
extern AlgDbgFn		algDbgOutFn;
extern AlgError		AlgDbgWrite(
			  char *,
			  ...);

#ifdef  __cplusplus
}
#endif

#endif /* ! ALGPROTO_H */
