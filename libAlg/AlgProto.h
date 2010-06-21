#ifndef ALGPROTO_H
#define ALGPROTO_H
#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgProto_h[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgProto.h
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
* \brief        Function prototypes for the Woolz numerical algorithm
*		library.
* \todo         -
* \bug          None known.
*/

#ifdef  __cplusplus
extern "C" {
#endif

/* From AlgAutoCorr.c */
extern AlgError			AlgAutoCorrelate2D(double **data,
				  int nX,
				  int nY);
/* From AlgBits.c */
extern int			AlgBitSetCount(
				  unsigned long gMsk);
extern int			AlgBitSetPositions(
				  int *posA,
				  unsigned long gMsk);
extern unsigned long		AlgBitNextNOfM(
				  unsigned long curMsk,
				  int n,
				  int m);
extern int			AlgBitNextSet(
				  unsigned long msk,
				  int idC);
extern int			AlgBitNextPowerOfTwo(
				  unsigned int *dstP2I,
				  unsigned int gI);
extern int			AlgBitIsPowerOfTwo(
				  unsigned int gI);

/* From AlgComplexUtils.c */
extern double			AlgCModSq(
				  ComplexD z);
extern double			AlgCMod(
				  ComplexD z);
extern double			AlgCArg(
			          ComplexD z);
extern double			AlgCRe(
			          ComplexD z);
extern double			AlgCIm(
			          ComplexD z);
extern ComplexD			AlgCConj(
			          ComplexD z);
extern ComplexD 		AlgCAdd(
			          ComplexD z1,
			          ComplexD z2);
extern ComplexD 		AlgCSub(
			          ComplexD z1,
			          ComplexD z2);
extern ComplexD 		AlgCMult(
			          ComplexD z1,
			          ComplexD z2);
extern ComplexD 		AlgCDiv(
			          ComplexD z1,
			          ComplexD z2);
extern ComplexD 		AlgCPow(
			          ComplexD z,
			          double pow);

/* From AlgConvolve.c */
extern AlgError 		AlgConvolve(
				  int sizeArrayCnv,
				  double *arrayCnv,
				  int sizeArrayKrn,
				  double *arrayKrn,
    	   	   		  int sizeArrayDat, 
				  double *arrayDat,
    				  AlgPadType pad);

/* From AlgCrossCorr.c */
extern AlgError        		AlgCrossCorrelate2D(
				  double **data0,
				  double **data1,
				  int nX,
				  int nY);
extern void            		AlgCrossCorrPeakXY(
				  int *dstMaxX,
				  int *dstMaxY,
				  double *dstMaxVal,
				  double **data,
				  int nX,
				  int nY,
				  int searchX,
				  int searchY);
extern void            		AlgCrossCorrPeakY(
				  int *dstMaxY,
				  double *dstMaxVal,
				  double **data,
				  int nY);

/* From AlgDebug.c */
extern AlgError			AlgDbgWrite(
				  char *fmt, ...);

/* From AlgFourier.c */
extern void     		AlgFourHart1D(
				  double *data,
				  int num, 
				  int step,
				  int cThr);
extern void			AlgFour1D(
				  double *real,
				  double *imag,
				  int num,
				  int step,
				  int cThr);
extern void			AlgFourInv1D(
				  double *real,
				  double *imag,
				  int num,
				  int step,
				  int cThr);
extern void			AlgFourReal1D(
				  double *real,
				  int num,
				  int step,
				  int cThr);
extern void			AlgFourRealInv1D(
				  double *real,
				  int num,
				  int step,
				  int cThr);
extern void			AlgFour2D(
				  double **real,
				  double **imag,
				  double *reBuf,
				  double *imBuf,
				  int numX,
				  int numY,
				  int cThr);
extern void			AlgFourInv2D(
				  double **real,
				  double **imag,
				  double *reBuf,
				  double *imBuf,
				  int numX,
				  int numY,
				 int cThr);
extern void			AlgFourReal2D(
				  double **real,
				  double *reBuf,
				  double *imBuf,
				  int numX,
				  int numY,
				  int cThr);
extern void			AlgFourRealInv2D(
				  double **real,
				  double *reBuf,
				  double *imBuf,
				  int numX,
				  int numY,
				  int cThr);

/* From AlgGamma.c */
extern double			AlgGammaLog(
				  double x,
				  AlgError *dstErr);
extern double			AlgGammaP(
				  double a,
				  double x,
				  AlgError *dstErr);

/* From AlgGamma.c */
extern AlgError			AlgLinearFit1D(
				  int datSz,
				  double *datXA,
				  double *datYA,
				  double *dstA,
				  double *dstB,
				  double *dstSigA,
				  double *dstSigB,
				  double *dstQ);
extern AlgError			AlgLinearFitIdx1D(
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
extern AlgError			AlgHeapSort(
				  void *data,
				  unsigned nElm,
				  unsigned elmSz,
				  int (*cmpFn)(void *, void *));
extern AlgError			AlgHeapSortIdx(
				  void *data,
				  int *idx,
				  unsigned nElm,
				  int (*cmpFn)(void *, int *, int, int));
extern void			AlgHeapElmSwap(
				  void *elm0,
				  void *elm1,
				  int cnt);
extern int             		AlgHeapSortCmpCFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpUFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpSFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpIFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpLFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpFFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpDFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpCFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpUFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpSFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpIFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpLFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpFFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortInvCmpDFn(
				  void *datum0,
				  void *datum1);
extern int             		AlgHeapSortCmpIdxCFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxUFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxSFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxIFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxLFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxFFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortCmpIdxDFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxCFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxUFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxSFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxIFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxLFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxFFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);
extern int             		AlgHeapSortInvCmpIdxDFn(
				  void *data,
				  int *idx,
				  int id0,
				  int id1);

/* From AlgMatrixGauss.c */
extern AlgError			AlgMatrixGaussSolve(
				  double **abMat,
				  int aSz,
				  double *xMat);

/* From AlgMatrixLSQR.c */
extern AlgError			AlgMatrixSolveLSQR(
				  AlgMatrixType aType,
				   double **aM,
				   size_t nR,
				   size_t nC,
				   double *bV,
				   double *xV,
				   double damping,
				   double relErrA,
				   double relErrB,
				   long maxItr,
				   long condLim,
				   int *dstTerm,
				   long *dstItr,
				   double *dstFNorm,
				   double *dstCondN,
				   double *dstResNorm,
				   double *dstResNormA,
				   double *dstNormX);
/* From AlgMatrixLU.c */
extern AlgError			AlgMatrixLUSolve(
				  double **aMat,
				  int aSz,
				  double *bMat,
				  int bSz);
extern AlgError			AlgMatrixLUInvert(
				  double **aMat,
				  int aSz);
extern AlgError			AlgMatrixLUDeterm(
				  double **aMat,
				  int aSz,
				  double *determ);
extern AlgError			AlgMatrixLUDecomp(
				  double **aMat,
				  int aSz,
				  int *idxVec,
				  double *evenOdd);
extern AlgError			AlgMatrixLUBackSub(
				  double **aMat,
				  int aSz,
				  int *idxVec,
				  double *bMat);
/* From AlgMatrixMath.c */
extern void			AlgMatrixAdd(
				  double **aM,
				  double **bM,
				  double **cM,
				  size_t nR,
				  size_t nC);
extern void            		AlgMatrixSub(
				  double **aM,
				  double **bM,
				  double **cM,
				  size_t nR,
				  size_t nC);
extern void            		AlgMatrixScale(
				  double **aM,
				  double **bM,
				  double sv,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixVectorScale(
				  double *aV,
				  double *bV,
				  double sv,
				  size_t nV);
extern void			AlgMatrixScaleAdd(
				  double **aM,
				  double **bM,
				  double **cM,
				  double sv,
				  size_t nR,
				  size_t nC);
extern void            		AlgMatrixMul(
				  double **aM,
				  double **bM,
				  double **cM,
				  size_t bR,
				  size_t bC,
				  size_t cC);
extern double          		AlgMatrixTrace(
				  double **aM,
				  size_t nRC);
extern void            		AlgMatrixTranspose(
				  double **aM,
				  double **bM,
				  size_t bR,
				  size_t bC);
extern void            		AlgMatrixCopy(
				  double **aM,
				  double **bM,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixVectorCopy(
				  double *aV,
				  double *bV,
				  size_t nV);
extern void			AlgMatrixScalar(
				  double **aM,
				  double sv,
				  size_t nRC);
extern void			AlgMatrixZero(
				  double **aM,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixVectorZero(
				  double *aV,
				  size_t nV);
extern void			AlgMatrixVectorMul(
				  double *aV,
				  AlgMatrixType bType,
				  double **bM,
				  double *cV,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixVectorMulAdd(
				  double *aV,
				  AlgMatrixType bType,
				  double **bM,
				  double *cV,
				  double *dV,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixVectorMulWAdd(
				  double *aV,
				  AlgMatrixType bType,
				  double **bM,
				  double *cV,
				  double *dV,
				  size_t nR,
				  size_t nC,
				  double s,
				  double t);
extern void			AlgMatrixTVectorMul(
				  double *aV,
                                  AlgMatrixType bType,
				  double **bM,
				  double *cV,
				  size_t nR,
				  size_t nC);
extern void			AlgMatrixTVectorMulAdd(
				  double *aV,
                                  AlgMatrixType bType,
				  double **bM,
				  double *cV,
				  double *dV,
				  size_t nR,
				  size_t nC);
extern double          		AlgMatrixVectorNorm(
				  double *aV,
				  size_t nV);

/* From AlgMatrixCG.c */
extern AlgError			AlgMatrixCGSolve(
				  AlgMatrixType aType,
				  double **aM,
				  double *xV,
				  double *bV,
				  double **wM,
				  size_t n,
				  void (*pFn)(void *,
				              double **, 
				  	      double *,
					      double *,
					      size_t),
				  void *pDat,
				  double tol,
				  int maxItr,
				  double *dstTol,
				  int *dstItr);

/* From AlgMatrixRSEigen.c */
extern AlgError			AlgMatrixRSEigen(
				  double **aM,
				  int aSz,
				  double *vM,
				  int reqEV);

/* From AlgMatrixRSTDiag.c */
extern AlgError			AlgMatrixRSTDiag(
				  double **aM,
				  int aSz,
				  double *dM,
				  double *oM);

/* From AlgMatrixSV.c */
extern AlgError			AlgMatrixSVSolve(
				  double **aMat,
				  int nM,
				  int nN,
				  double *bMat,
				  double tol,
				  int *dstIC);
extern AlgError			AlgMatrixSVDecomp(double **aMat,
				  int nM,
				  int nN,
				  double *wMat,
				  double **vMat);
extern AlgError			AlgMatrixSVBackSub(double **uMat,
				  int nM,
				  int nN,
				  double *wMat,
				  double **vMat,
				  double *bMat);

/* From AlgMatrixTDiagQLI.c */
extern AlgError			AlgMatrixTDiagQLI(
				  double *dM,
				  double *oM,
				  int aSz,
				  double **zM);

/* From AlgMixture.c */
AlgError			AlgMixtureMLG(
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
AlgError			AlgMixtureSyn(
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
extern AlgError	        	AlgPolynomialLSq(
				  double *xVec,
				  double *yVec,
				  int vecSz,
				  int polyDeg,
				  double *cVec);

/* From AlgQSort.c */
#ifndef WLZ_EXT_BIND
extern void			AlgQSort(
				  void *base,
				  size_t nElm,
				  size_t elmSz,
			          void *cData,
		                  int (*cmpFn)(const void *,
				               const void *,
					       const void *));
#endif /* WLZ_EXT_BIND */

/* From AlgRand.c */
extern void			AlgRandSeed(
				  long seed);
extern double			AlgRandUniform(
				  void);
extern double			AlgRandNormal(
				  double mu,
				  double sigma);

/* From AlgRandZig.c */
extern double			AlgRandZigNormal(
				  double mu,
				  double sigma);

/* From AlgRange.c */
extern AlgError			AlgRange1D(
				  int datASz,
				  double *datA,
				  double *dstMin,
				  double *dstMax);
extern AlgError			AlgRangeIdx1D(
				  double *datA,
				  int idxASz,
				  int *idxA,
				  double *dstMin,
				  double *dstMax);

/* From AlgRank.c */
extern void	 		AlgRankSelectI(
				  int *elm,
				  int nElm,
				  int rank);
extern void    	        	AlgRankSelectUB(
				  unsigned char *elm,
				  int nElm,
				  int rank);
extern void    	        	AlgRankSelectS(
				  short *elm,
				  int nElm,
				  int rank);
extern void    	        	AlgRankSelectF(
				  float *elm,
				  int nElm,
				  int rank);
extern void    	        	AlgRankSelectD(
				  double *elm,
				  int nElm,
				  int rank);
extern void    	        	AlgRankSelectV(
				  void *elm,
				  int nElm,
				  unsigned int elmSz,
				  int rank,
				  void *buf,
				  int (*compFn)(void *, void *));

/* From AlgShuffle.c */
extern AlgError			AlgShuffleIdx(
				  int nShuffle,
				  int *shuffle,
				  int seed);

/* From AlgVectorMath.c */
extern double			AlgVectorNorm(
				  double *aV,
				  size_t n);
extern double          		AlgVectorDot(
				  double *aV,
				  double *bV,
				  size_t n);
extern void            		AlgVectorAdd(
				  double *aV,
				  double *bV,
				  double *cV,
				  size_t n);
extern void            		AlgVectorSub(
				  double *aV,
				  double *bV,
				  double *cV,
				  size_t n);
extern void            		AlgVectorCopy(
				  double *aV,
				  double *bV,
				  size_t n);
extern void            		AlgVectorScaleAdd(
				  double *aV,
				  double *bV,
				  double *cV,
				  double s,
				  size_t n);

/* From AlgDPSearch.c */
extern int AlgDPSearch(int, int, double **, double **, int **,
		       double (*)(int, int, int, int **));

/* Debugging */
extern AlgDbgMask	 	algDbgMask;
extern void			*algDbgData;
extern AlgDbgFn			algDbgOutFn;
extern AlgError			AlgDbgWrite(
				  char *,
				  ...);

#ifdef  __cplusplus
}
#endif

#endif /* ! ALGPROTO_H */
