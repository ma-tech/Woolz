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
* 26-01-00 bill Add AlgConvolve() and AlgMixture().
************************************************************************/

#ifdef  __cplusplus
extern "C" {
#endif

/* From AlgComplexUtils.c */
extern double	AlgCMod(ComplexD z);
extern double	AlgCArg(ComplexD z);
extern double	AlgCRe(ComplexD z);
extern double	AlgCIm(ComplexD z);
extern ComplexD	AlgCConj(ComplexD z);
extern ComplexD AlgCAdd(ComplexD z1, ComplexD z2);
extern ComplexD AlgCSub(ComplexD z1, ComplexD z2);
extern ComplexD AlgCMult(ComplexD z1, ComplexD z2);
extern ComplexD AlgCDiv(ComplexD z1, ComplexD z2);
extern ComplexD AlgCPow(ComplexD z, double pow);

/* From AlgDebug.c */
extern AlgError	AlgDbgWrite(char *fmt, ...);

/* From AlgConvolve.c */
extern AlgError AlgConvolve(int sizeArrayCnv, double *arrayCnv,
    			    int sizeArrayKrn, double *arrayKrn,
    	   	   	    int sizeArrayDat, double *arrayDat,
    			    AlgPadType pad);

/* From AlgFourier.c */
extern void     AlgFourHart1D(double *data, int num, int step, int cThr),
		AlgFour1D(double *real, double *imag, int num, int step,
			  int cThr),
		AlgFourInv1D(double *real, double *imag, int num, int step,
			     int cThr),
		AlgFourReal1D(double *real, int num, int step, int cThr),
		AlgFourRealInv1D(double *real, int num, int step, int cThr),
		AlgFour2D(double **real, double **imag,
			  double *reBuf, double *imBuf,int numX, int numY,
			  int cThr),
		AlgFourInv2D(double **real, double **imag,
			     double *reBuf, double *imBuf, int numX, int numY,
			     int cThr),
		AlgFourReal2D(double **real, double *reBuf, double *imBuf,
			      int numX, int numY, int cThr),
		AlgFourRealInv2D(double **real, double *reBuf, double *imBuf,
				 int numX, int numY, int cThr);

/* From AlgMatrixGauss.c */
extern AlgError	AlgMatrixGaussSolve(double **abMat, int aSz,
                                    double *xMat);

/* From AlgMatrixLU.c */
extern AlgError	AlgMatrixLUSolve(double **aMat, int aSz,
                                 double *bMat, int bSz),
		AlgMatrixLUInvert(double **aMat, int aSz),
		AlgMatrixLUDeterm(double **aMat, int aSz, double *determ),
		AlgMatrixLUDecomp(double **aMat, int aSz,
                                  int *idxVec, double *evenOdd),
		AlgMatrixLUBackSub(double **aMat, int aSz,
                                  int *idxVec, double *bMat);

/* From AlgMatrixSV.c */
extern AlgError	AlgMatrixSVSolve(double **aMat, int nM, int nN,
				 double *bMat, double tol),
		AlgMatrixSVDecomp(double **aMat, int nM, int nN,
				  double *wMat, double **vMat),
		AlgMatrixSVBackSub(double **uMat, int nM, int nN,
				   double *wMat, double **vMat,
				   double *bMat);

/* From AlgMixture.c */
extern AlgError  AlgMixture(AlgDistribution dbnType, int nDbn, int nCls,
		    	    double *x, int *freq, double *alpha, double *mean,
		    	    double *sd, double tol, int nObv, double *dstLL,
		    	    int *nItn);

/* From AlgPolyLSQ.c */
extern AlgError        AlgPolynomialLSq(double *xVec, double *yVec,
                                        int vecSz, int polyDeg,
					double *cVec);
/* From AlgRand.c */
extern void	AlgRandSeed(long seed);
extern double	AlgRandUniform(void),
		AlgRandNormal(double mu, double sigma);
/* Debugging */
extern AlgDbgMask algDbgMask;
extern void	*algDbgData;
extern AlgDbgFn	algDbgOutFn;
extern AlgError	AlgDbgWrite(char *, ...);


#ifdef  __cplusplus
}
#endif

#endif /* ! ALGPROTO_H */
