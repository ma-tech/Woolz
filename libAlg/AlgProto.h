#ifndef ALGPROTO_H
#define ALGPROTO_H
#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgProto_h[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgProto.h
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
* \brief        Function prototypes for the Woolz numerical algorithm
*		library.
* \todo         -
* \bug          None known.
*/

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
extern "C" {
#endif
#endif /* WLZ_EXT_BIND */

/* From AlgAutoCorr.c */
extern AlgError			AlgAutoCorrelate2D(double **data,
				  int nX,
				  int nY);
/* From AlgBits.c */
extern unsigned int		AlgBitSetCount(
				  unsigned long gMsk);
extern unsigned int		AlgBitSetPositions(
				  unsigned int *posA,
				  unsigned long gMsk);
extern unsigned long		AlgBitNextNOfM(
				  unsigned long curMsk,
				  int n,
				  int m);
extern int			AlgBitMostSigSet(
				  unsigned long gMsk);
extern int			AlgBitMostSigSetLL(
				  unsigned long long gMsk);
extern unsigned int		AlgBitRotateRight(
				  unsigned int g,
				  unsigned int n,
				  unsigned int d);
extern unsigned int		AlgBitRotateLeft(
				  unsigned int g,
				  unsigned int n,
				  unsigned int d);
extern int			AlgBitNextSet(
				  unsigned long msk,
				  int idC);
extern int			AlgBitNextPowerOfTwo(
				  unsigned int *dstP2I,
				  unsigned int gI);
extern int			AlgBitIsPowerOfTwo(
				  unsigned int gI);

/* From AlgBSpline.c */
extern void                     AlgBSplineBspl(
                                  double *t,
                                  int k,
                                  double x,
                                  int l,
                                  double *h);
extern AlgError			AlgBSplineNDFit(
                                  int iopt,
                                  int ipar,
                                  int idim,
                                  int m,
                                  double *u,
                                  int mx,
                                  double *x,
                                  double *w,
                                  double ub,
                                  double ue,
                                  int k,
                                  double s,
                                  int nest,
                                  int *n,
                                  double *t,
                                  int *nc,
                                  double *c,
                                  double *fp,
                                  double *wrk,
                                  int *iwrk);
extern AlgError			AlgBSplinePerFit(
				  int iopt,
				  int m,
				  double *x,
				  double *y,
				  double *w,
				  int k,
				  double s,
				  int nest,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *wrk,
				  int *iwrk);
extern AlgError			AlgBSplineFit(
				  int iopt,
				  int m,
				  double *x,
				  double *y,
				  double *w,
				  double xb,
				  double xe,
				  int k,
				  double s,
				  int nest,
				  int *n,
				  double *t,
				  double *c,
				  double *fp,
				  double *wrk,
				  int *iwrk);
extern AlgError			AlgBSplineEval(
                                  double *t,
                                  int n,
                                  double *c,
                                  int k,
                                  double * x,
                                  double *y,
                                  int m);
extern AlgError			AlgBSplineDer(
				  double *t,
                                  int n,
                                  double *c,
                                  int k,
                                  int nu,
                                  double *x,
                                  double *y,
                                  int m,
                                  double *wrk);
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
extern AlgError 		AlgConvolveD(
				  int sizeArrayCnv,
				  double *arrayCnv,
				  int sizeArrayKrn,
				  double *arrayKrn,
    	   	   		  int sizeArrayDat, 
				  double *arrayDat,
    				  AlgPadType pad,
				  double padVal);
extern AlgError 		AlgConvolveF(
				  int sizeArrayCnv,
				  float *arrayCnv,
				  int sizeArrayKrn,
				  float *arrayKrn,
    	   	   		  int sizeArrayDat, 
				  float *arrayDat,
    				  AlgPadType pad,
				  float padVal);

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
#ifndef CTYPESGEN
extern AlgError			AlgDbgWrite(
				  char *fmt, ...);
#endif

/* From AlgFourier.c */
extern void     		AlgFourHart1D(
				  double *data,
				  int num, 
				  int step);
extern void			AlgFour1D(
				  double *real,
				  double *imag,
				  int num,
				  int step);
extern void			AlgFourInv1D(
				  double *real,
				  double *imag,
				  int num,
				  int step);
extern void			AlgFourReal1D(
				  double *real,
				  int num,
				  int step);
extern void			AlgFourRealInv1D(
				  double *real,
				  int num,
				  int step);
extern AlgError			AlgFourHart2D(
				  double **data,
				  int useBuf,
				  int numX,
				  int numY);
extern AlgError			AlgFour2D(
				  double **real,
				  double **imag,
				  int useBuf,
				  int numX,
				  int numY);
extern AlgError			AlgFourInv2D(
				  double **real,
				  double **imag,
				  int useBuf,
				  int numX,
				  int numY);
extern AlgError			AlgFourReal2D(
				  double **real,
				  int useBuf,
				  int numX,
				  int numY);
extern AlgError			AlgFourRealInv2D(
				  double **real,
				  int useBuf,
				  int numX,
				  int numY);
extern AlgError			AlgFour3D(
				  double ***real,
				  double ***imag,
				  int useBuf,
				  int numX,
				  int numY,
				  int numZ);
extern AlgError			AlgFourInv3D(
				  double ***real,
				  double ***imag,
				  int useBuf,
				  int numX,
				  int numY,
				  int numZ);
extern AlgError        		AlgFourReal3D(
				  double ***real,
				  int useBuf,
				  int numX,
				  int numY,
				  int numZ);
extern AlgError        		AlgFourRealInv3D(
				  double ***real,
				  int useBuf,
				  int numX,
				  int numY,
				  int numZ);

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

/* From AlgGrayCode.c */
extern unsigned int		AlgGrayCode(
				  unsigned int g);
extern unsigned int    		AlgGrayCodeInv(
				  unsigned int g);
extern unsigned long long 	AlgGrayCodeLL(
				  unsigned long long g);
extern unsigned long long	AlgGrayCodeInvLL(
				  unsigned long long g);

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

/* From AlgHilbertIndex.c */
extern void			AlgHilbertIndex(
				  unsigned int *h,
				  unsigned int *p,
				  int n,
				  int o);
extern void			AlgHilbertIndexInv(
				  unsigned int *h,
				  unsigned int *p,
				  int n,
				  int o);

/* From AlgMatrix.c */
extern void			AlgMatrixFree(
				  AlgMatrix mat);
extern void			AlgMatrixRectFree(
				  AlgMatrixRect *mat);
extern void			AlgMatrixSymFree(
				  AlgMatrixSym *mat);
extern void			AlgMatrixLLRFree(
				  AlgMatrixLLR *mat);
extern AlgMatrixLLRE		*AlgMatrixLLRENew(
				  AlgMatrixLLR *mat);
extern void			AlgMatrixLLREFree(
				  AlgMatrixLLR *mat,
				  AlgMatrixLLRE *p);
extern AlgMatrix		AlgMatrixNew(
				  AlgMatrixType aType,
				  size_t nR,
				  size_t nC,
				  size_t nE,
				  double tol,
				  AlgError *dstErr);
extern AlgMatrixRect   		*AlgMatrixRectNew(
				  size_t nR,
				  size_t nC,
				  AlgError *dstErr);
extern AlgMatrixSym   		*AlgMatrixSymNew(
				  size_t nN,
				  AlgError *dstErr);
extern AlgMatrixLLR		*AlgMatrixLLRNew(
				  size_t nR,
				  size_t nC,
				  size_t nE,
				  double tol,
				  AlgError *dstErr);
extern AlgError			AlgMatrixLLRCopyInPlace(
				  AlgMatrixLLR *aM,
				  AlgMatrixLLR *bM);
extern void			AlgMatrixZero(
				  AlgMatrix mat);
extern void            		AlgMatrixRectZero(
				  AlgMatrixRect *mat);
extern void            		AlgMatrixSymZero(
				  AlgMatrixSym *mat);
extern void            		AlgMatrixLLRZero(
				  AlgMatrixLLR *mat);
extern AlgError			AlgMatrixSetAll(
  				  AlgMatrix mat,
				  double val);
extern void			AlgMatrixRectSetAll(
				  AlgMatrixRect *mat,
				  double val);
extern void			AlgMatrixRectSetAll(
				  AlgMatrixRect *mat,
				  double val);
extern void 			AlgMatrixSymSetAll(
				  AlgMatrixSym *mat,
				  double val);
extern AlgError			AlgMatrixLLRSetAll(
				  AlgMatrixLLR *mat,
				  double val);
extern void			AlgMatrixLLRERemove(
				  AlgMatrixLLR *mat,
				  size_t row,
				  size_t col);
extern AlgError			AlgMatrixSet(
				  AlgMatrix mat,
                                  size_t row,
				  size_t col,
				  double val);
extern AlgError			AlgMatrixLLRSet(
				  AlgMatrixLLR *mat,
                                  size_t row,
				  size_t col,
				  double val);
extern double			AlgMatrixValue(
				  AlgMatrix mat,
				  size_t row,
				  size_t col);
extern double			AlgMatrixLLRValue(
				  AlgMatrixLLR *mat,
				  size_t row,
				  size_t col);
extern AlgError			AlgMatrixLLRExpand(
				  AlgMatrixLLR *mat,
				  size_t nE);
extern AlgMatrix		AlgMatrixReadAscii(
				  AlgMatrixType mType,
				  double tol,
				  FILE *fP,
				  const char *fSep,
				  size_t recMax,
				  AlgError *dstErr);
extern AlgError        		AlgMatrixWriteAscii(
				  AlgMatrix mat,
				  FILE *fP);


/* From AlgMatrixGauss.c */
extern AlgError			AlgMatrixGaussSolve(
				  AlgMatrix aMat,
				  double *xMat);

/* From AlgMatrixLSQR.c */
extern AlgError        		AlgMatrixSolveLSQR(
				  AlgMatrix aM,
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
extern AlgError        		AlgMatrixLUSolveRaw3(
				  double **aM,
				  double *bV,
				  int bSz);
extern AlgError        		AlgMatrixLUSolveRaw4(
				  double **aM,
				  double *bV,
				  int bSz);
extern AlgError        		AlgMatrixLUSolve(
				  AlgMatrix aM,
				  double *bV,
				  int bSz);
extern AlgError        		AlgMatrixLUSolveRaw(
				  double **aM,
				  int aSz,
				  double *bV,
				  int bSz);
extern AlgError        		AlgMatrixLUInvertRaw3(
				  double **aM);
extern AlgError        		AlgMatrixLUInvertRaw4(
				  double **aM);
extern AlgError        		AlgMatrixLUInvert(
				  AlgMatrix aM);
extern AlgError        		AlgMatrixLUInvertRaw(
				  double **aM,
				  int aSz);
extern AlgError        		AlgMatrixLUDetermRaw3(
				  double **aM,
				  double *det);
extern AlgError        		AlgMatrixLUDetermRaw4(
				  double **aM,
				  double *det);
extern AlgError        		AlgMatrixLUDeterm(
				  AlgMatrix aM,
				  double *det);
extern AlgError        		AlgMatrixLUDetermRaw(
				  double **aM,
				  int aSz,
				  double *det);
extern AlgError        		AlgMatrixLUDecomp(
				  AlgMatrix aM,
				  int *iV,
				  double *evenOdd);
extern AlgError        		AlgMatrixLUDecompRaw(
				  double **aM,
				  int aSz,
				  int *iV,
				  double *evenOdd);
extern AlgError        		AlgMatrixLUBackSub(
				  AlgMatrix aM,
				  int *iV,
				  double *bV);
extern AlgError        		AlgMatrixLUBackSubRaw(
				  double **aM,
				  int aSz,
				  int *iV,
				  double *bV);

/* From AlgMatrixMath.c */
extern void			AlgMatrixAdd(
				  AlgMatrix aM,
				  AlgMatrix bM,
				  AlgMatrix cM);
extern void			AlgMatrixSub(
				  AlgMatrix aM,
				  AlgMatrix bM,
				  AlgMatrix cM);
extern void 			AlgMatrixMul(
				  AlgMatrix aM,
				  AlgMatrix bM,
				  AlgMatrix cM);
extern double          		AlgMatrixTrace(
				  AlgMatrix aM);
extern void 			AlgMatrixTranspose(
				  AlgMatrix aM,
				  AlgMatrix bM);
extern void			AlgMatrixCopy(
				  AlgMatrix aM,
				  AlgMatrix bM);
extern void 			AlgMatrixScale(
				  AlgMatrix aM,
				  AlgMatrix bM,
				  double sv);
extern void			AlgMatrixScaleAdd(
				  AlgMatrix aM,
				  AlgMatrix bM,
				  AlgMatrix cM,
                                  double sv);
extern void			AlgMatrixScalar(
				  AlgMatrix aM,
				  double sv);
extern void			AlgMatrixVectorMul(
				  double *aV,
				  AlgMatrix bM,
				  double *cV);
extern void			AlgMatrixVectorMulAdd(
				  double *aV,
				  AlgMatrix bM,
				  double *cV,
				  double *dV);
extern void 			AlgMatrixVectorMulWAdd(
				  double *aV,
				  AlgMatrix bM,
				  double *cV,
				  double *dV,
				  double s,
				  double t);
extern void 			AlgMatrixTVectorMul(
				  double *aV,
				  AlgMatrix bM,
				  double *cV);
extern void			AlgMatrixTVectorMulAdd(
				  double *aV,
				  AlgMatrix bM,
				  double *cV,
				  double *dV);
extern void			AlgMatrixZero(
				  AlgMatrix mat);
extern AlgError			AlgMatrixRawInv2x2(
				  double *a00,
				  double *a01,
				  double *a10,
				  double *a11);
extern AlgError			AlgMatrixRawInv3x3(
				  double *a00,
				  double *a01,
				  double *a02,
				  double *a10,
				  double *a11,
				  double *a12,
				  double *a20,
				  double *a21,
				  double *a22);

/* From AlgMatrixCG.c */
extern AlgError			AlgMatrixCGSolve(
				  AlgMatrix aM,
                                  double *xV, 
				  double *bV,
                                  AlgMatrix wM,
                                  void (*pFn)(void *,
				              AlgMatrix,
                                              double *,
					      double *),
                                  void *pDat,
				  double tol,
				  int maxItr,
                                  double *dstTol,
				  int *dstItr);

/* From AlgMatrixRSEigen.c */
extern AlgError        		AlgMatrixRSEigen(
				  AlgMatrix aM,
				  double *vM,
				  int reqEV);

/* From AlgMatrixRSTDiag.c */
extern AlgError			AlgMatrixRSTDiag(
				  AlgMatrix aMat,
				  double *dVec,
				  double *oVec);

/* From AlgMatrixSV.c */
extern AlgError			AlgMatrixSVSolve(
				  AlgMatrix aMat,
				  double *bVec,
				  double tol,
				  int *dstIC);
extern AlgError			AlgMatrixSVDecomp(
				  AlgMatrix aMat,
				  double *wVec,
				  AlgMatrix vMat);
extern AlgError			AlgMatrixSVBackSub(
				  AlgMatrix uMat,
				  double *wVec,
				  AlgMatrix vMat,
				  double *bVec);

/* From AlgMatrixTDiagQLI.c */
extern AlgError			AlgMatrixTDiagQLI(
				  double *dM,
				  double *oM,
				  int aSz,
				  AlgMatrix zM);

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

/* From AlgSort.c */
#ifndef WLZ_EXT_BIND
extern void			AlgSort(
				  void *base,
				  size_t nElm,
				  size_t elmSz,
		                  int (*cmpFn)(const void *,
					       const void *));
#endif /* WLZ_EXT_BIND */

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
extern void 			AlgVectorZero(
				  double *aV,
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
extern void            		AlgVectorScale(
				  double *aV,
				  double *bV,
				  double s,
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
#ifndef CTYPESGEN
extern AlgDbgMask	 	algDbgMask;
extern void			*algDbgData;
extern AlgDbgFn			algDbgOutFn;
extern AlgError			AlgDbgWrite(
				  char *,
				  ...);
#endif

#ifndef WLZ_EXT_BIND
#ifdef  __cplusplus
}
#endif
#endif /* WLZ_EXT_BIND */

#endif /* ! ALGPROTO_H */
