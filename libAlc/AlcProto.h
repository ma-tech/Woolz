#ifndef ALCPROTO_H
#define ALCPROTO_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcProto.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Header file which contains function prototypes for
*		the MRC HGU memory allocation library.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

extern void	*AlcCalloc(int, int),
		*AlcMalloc(int),
		*AlcRealloc(void *, int);

extern void	AlcFree(void *);

extern char	*AlcStrDup(const char *srcStr);

extern AlcErrno	AlcBit1Calloc(unsigned char **, int),
		AlcChar1Calloc(char **, int),
		AlcUnchar1Calloc(unsigned char **, int),
		AlcShort1Calloc(short **, int),
		AlcInt1Calloc(int **, int),
		AlcInt1Calloc(int **, int),
		AlcFloat1Calloc(float **, int),
		AlcDouble1Calloc(double **, int),
		AlcBit1Malloc(unsigned char **, int),
		AlcChar1Malloc(char **, int),
		AlcUnchar1Malloc(unsigned char **, int),
		AlcShort1Malloc(short **, int),
		AlcInt1Malloc(int **, int),
		AlcFloat1Malloc(float **, int),
		AlcDouble1Malloc(double **, int),
		AlcBit2Calloc(unsigned char ***, int, int),
		AlcChar2Calloc(char ***, int, int),
		AlcUnchar2Calloc(unsigned char ***, int, int),
		AlcShort2Calloc(short ***, int, int),
		AlcInt2Calloc(int ***, int, int),
		AlcFloat2Calloc(float ***, int, int),
		AlcDouble2Calloc(double ***, int, int),
		AlcBit2Malloc(unsigned char ***, int, int),
		AlcChar2Malloc(char ***, int, int),
		AlcUnchar2Malloc(unsigned char ***, int, int),
		AlcShort2Malloc(short ***, int, int),
		AlcInt2Malloc(int ***, int, int),
		AlcFloat2Malloc(float ***, int, int),
		AlcDouble2Malloc(double ***, int, int),
		Alc2Free(void **),
		AlcBit2Free(unsigned char **),
		AlcChar2Free(char **),
		AlcUnchar2Free(unsigned char **),
		AlcShort2Free(short **),
		AlcInt2Free(int **),
		AlcFloat2Free(float **),
		AlcDouble2Free(double **),
		AlcBit3Calloc(unsigned char ****, int, int, int),
		AlcChar3Calloc(char ****, int, int, int),
		AlcUnchar3Calloc(unsigned char ****, int, int, int),
		AlcShort3Calloc(short ****, int, int, int),
		AlcInt3Calloc(int ****, int, int, int),
		AlcFloat3Calloc(float ****, int, int, int),
		AlcDouble3Calloc(double ****, int, int, int),
		AlcBit3Malloc(unsigned char ****, int, int, int),
		AlcChar3Malloc(char ****, int, int, int),
		AlcUnchar3Malloc(unsigned char ****, int, int, int),
		AlcShort3Malloc(short ****, int, int, int),
		AlcInt3Malloc(int ****, int, int, int),
		AlcFloat3Malloc(float ****, int, int, int),
		AlcDouble3Malloc(double ****, int, int, int),
		Alc3Free(void ***),
		AlcBit3Free(unsigned char ***),
		AlcChar3Free(char ***),
		AlcUnchar3Free(unsigned char ***),
		AlcShort3Free(short ***),
		AlcInt3Free(int ***),
		AlcFloat3Free(float ***),
		AlcDouble3Free(double ***);

#ifdef __cplusplus
}
#endif

#endif /* ALCPROTO_H */
