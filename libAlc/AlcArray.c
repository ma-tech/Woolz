#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcArray.c
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Provides functions for the allocation of 1, 2 and 3D
*		arrays of types char, short, int, float and double.
*		Extension to other types (including user defined types)
*		should be straight formward through templates defined
*		in AlcTemplates.h.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>
#include <AlcTemplates.h>

/************************************************************************
* Functions:	AlcChar1Calloc, AlcShort1Calloc, AlcInt1Calloc,		*
*		AlcFloat1Calloc, AlcDouble1Calloc and AlcUnchar1Calloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 1 dimensional zero'd	*
*		arrays of fundamental types.				*
*		All 1 dimensional arrays allocated by these functions	*
*		should be free'd using AlcFree().			*
* Global refs:	-							*
* Parameters:	<TYPE> **dest:		Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of elements in array.	*
************************************************************************/
AlcErrno	AlcChar1Calloc(char **dest, int mElem)
ALC_TEMPLATE_C1D(dest, char, mElem, "AlcChar1Calloc")

AlcErrno	AlcUnchar1Calloc(unsigned char **dest, int mElem)
{
  return(AlcChar1Calloc((char **)dest, mElem));
}

AlcErrno	AlcShort1Calloc(short **dest, int mElem)
ALC_TEMPLATE_C1D(dest, short, mElem, "AlcShort1Calloc")

AlcErrno	AlcInt1Calloc(int **dest, int mElem)
ALC_TEMPLATE_C1D(dest, int, mElem, "AlcInt1Calloc")

AlcErrno	AlcFloat1Calloc(float **dest, int mElem)
ALC_TEMPLATE_C1D(dest, float, mElem, "AlcFloat1Calloc")

AlcErrno	AlcDouble1Calloc(double **dest, int mElem)
ALC_TEMPLATE_C1D(dest, double, mElem, "AlcDouble1Calloc")

/************************************************************************
* Functions:	AlcChar1Malloc, AlcShort1Malloc, AlcInt1Malloc,		*
*		AlcFloat1Malloc, AlcDouble1Malloc and AlcUnchar1Malloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 1 dimensional non-	*
*		zero'd arrays of fundamental types.			*
*		All 1 dimensional arrays allocated by these functions	*
*		should be free'd using AlcFree().			*
* Global refs:	-							*
* Parameters:	<TYPE> **dest:		Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of elements in array.	*
************************************************************************/
AlcErrno	AlcChar1Malloc(char **dest, int mElem)
ALC_TEMPLATE_M1D(dest, char, mElem, "AlcChar1Malloc")

AlcErrno	AlcUnchar1Malloc(unsigned char **dest, int mElem)
{
  return(AlcChar1Malloc((char **)dest, mElem));
}

AlcErrno	AlcShort1Malloc(short **dest, int mElem)
ALC_TEMPLATE_M1D(dest, short, mElem, "AlcShort1Malloc")

AlcErrno	AlcInt1Malloc(int **dest, int mElem)
ALC_TEMPLATE_M1D(dest, int, mElem, "AlcInt1Malloc")

AlcErrno	AlcFloat1Malloc(float **dest, int mElem)
ALC_TEMPLATE_M1D(dest, float, mElem, "AlcFloat1Malloc")

AlcErrno	AlcDouble1Malloc(double **dest, int mElem)
ALC_TEMPLATE_M1D(dest, double, mElem, "AlcDouble1Malloc")

/************************************************************************
* Functions:	AlcChar2Calloc, AlcShort2Calloc, AlcInt2Calloc,		*
*		AlcFloat2Calloc, AlcDouble2Calloc and AlcUnchar2Calloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 2 dimensional zero'd	*
*		arrays of fundamental types.				*
*		All 2 dimensional arrays allocated by these functions	*
*		should be free'd using one of the associated free	*
*		functions, eg AlcChar2Calloc and AlcChar2Free.		*
* Global refs:	-							*
* Parameters:	<TYPE> ***dest:		Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of 1D arrays.		*
*		int nElem:		Number of elements in each 1D	*
*					array.				*
************************************************************************/
AlcErrno	AlcChar2Calloc(char ***dest, int mElem, int nElem)
ALC_TEMPLATE_C2D(dest, char, mElem, nElem, "AlcChar2Calloc")

AlcErrno	AlcUnchar2Calloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcChar2Calloc((char ***)dest, mElem, nElem));
}

AlcErrno	AlcShort2Calloc(short ***dest, int mElem, int nElem)
ALC_TEMPLATE_C2D(dest, short, mElem, nElem, "AlcShort2Calloc")

AlcErrno	AlcInt2Calloc(int ***dest, int mElem, int nElem)
ALC_TEMPLATE_C2D(dest, int, mElem, nElem, "AlcInt2Calloc")

AlcErrno	AlcFloat2Calloc(float ***dest, int mElem, int nElem)
ALC_TEMPLATE_C2D(dest, float, mElem, nElem, "AlcFloat2Calloc")

AlcErrno	AlcDouble2Calloc(double ***dest, int mElem, int nElem)
ALC_TEMPLATE_C2D(dest, double, mElem, nElem, "AlcDouble2Calloc")

/************************************************************************
* Functions:	AlcChar2Malloc, AlcShort2Malloc, AlcInt2Malloc,		*
*		AlcFloat2Malloc, AlcDouble2Malloc and AlcUnchar2Malloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 2 dimensional non-	*
*		zero'd arrays of fundamental types.			*
*		All 2 dimensional arrays allocated by these functions	*
*		should be free'd using one of the associated free	*
*		functions, eg AlcChar2Malloc and AlcChar2Free.		*
* Global refs:	-							*
* Parameters:	<TYPE> ***dest:		Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of 1D arrays.		*
*		int nElem:		Number of elements in each 1D	*
*					array.				*
************************************************************************/
AlcErrno	AlcChar2Malloc(char ***dest, int mElem, int nElem)
ALC_TEMPLATE_M2D(dest, char, mElem, nElem, "AlcChar2Malloc")

AlcErrno	AlcUnchar2Malloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcChar2Malloc((char ***)dest, mElem, nElem));
}

AlcErrno	AlcShort2Malloc(short ***dest, int mElem, int nElem)
ALC_TEMPLATE_M2D(dest, short, mElem, nElem, "AlcShort2Malloc")

AlcErrno	AlcInt2Malloc(int ***dest, int mElem, int nElem)
ALC_TEMPLATE_M2D(dest, int, mElem, nElem, "AlcInt2Malloc")

AlcErrno	AlcFloat2Malloc(float ***dest, int mElem, int nElem)
ALC_TEMPLATE_M2D(dest, float, mElem, nElem, "AlcFloat2Malloc")

AlcErrno	AlcDouble2Malloc(double ***dest, int mElem, int nElem)
ALC_TEMPLATE_M2D(dest, double, mElem, nElem, "AlcDouble2Malloc")

/************************************************************************
* Functions:	Alc2Free, AlcChar2Free, AlcShort2Free, AlcInt2Free,	*
*		AlcFloat2Free, AlcDouble2Free and AlcUnchar2Free.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for freeing 2 dimensional arrays	*
*		allocated by the 2 dimensional array allocation 	*
*		functions.						*
* Global refs:	-							*
* Parameters:	<TYPE> **dest:		Ptr with array to be free'd.	*
************************************************************************/
AlcErrno	Alc2Free(void **dest)
ALC_TEMPLATE_F2D(dest, "Alc2Free")

AlcErrno	AlcChar2Free(char **dest)
ALC_TEMPLATE_F2D(dest, "AlcChar2Free")

AlcErrno	AlcUnchar2Free(unsigned char **dest)
{
  return(AlcChar2Free((char **)dest));
}

AlcErrno	AlcShort2Free(short **dest)
ALC_TEMPLATE_F2D(dest, "AlcShort2Free")

AlcErrno	AlcInt2Free(int **dest)
ALC_TEMPLATE_F2D(dest, "AlcInt2Free")

AlcErrno	AlcFloat2Free(float **dest)
ALC_TEMPLATE_F2D(dest, "AlcFloat2Free")

AlcErrno	AlcDouble2Free(double **dest)
ALC_TEMPLATE_F2D(dest, "AlcDouble2Free")

/************************************************************************
* Functions:	AlcChar3Calloc, AlcShort3Calloc, AlcInt3Calloc,		*
*		AlcFloat3Calloc, AlcDouble3Calloc and AlcUnchar3Calloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 3 dimensional zero'd	*
*		arrays of fundamental types.				*
*		All 3 dimensional arrays allocated by these functions	*
*		should be free'd using one of the associated free	*
*		functions, eg AlcChar3Calloc and AlcChar3Free.		*
* Global refs:	-							*
* Parameters:	<TYPE> ****dest:	Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of 2D arrays.		*
*		int nElem:		Number of 1D arrays.		*
*		int oElem:		Number of elements in each 1D	*
*					array.				*
************************************************************************/
AlcErrno	AlcChar3Calloc(char ****dest, int mElem, int nElem,
			      int oElem)
ALC_TEMPLATE_C3D(dest, char, mElem, nElem, oElem, "AlcChar3Calloc")

AlcErrno	AlcUnchar3Calloc(unsigned char ****dest,
			   int mElem, int nElem, int oElem)
{
  return(AlcChar3Calloc((char ****)dest, mElem, nElem, oElem));
}

AlcErrno	AlcShort3Calloc(short ****dest, int mElem, int nElem,
			       int oElem)
ALC_TEMPLATE_C3D(dest, short, mElem, nElem, oElem, "AlcShort3Calloc")

AlcErrno	AlcInt3Calloc(int ****dest, int mElem, int nElem, int oElem)
ALC_TEMPLATE_C3D(dest, int, mElem, nElem, oElem, "AlcInt3Calloc")

AlcErrno	AlcFloat3Calloc(float ****dest, int mElem, int nElem,
			       int oElem)
ALC_TEMPLATE_C3D(dest, float, mElem, nElem, oElem, "AlcFloat3Calloc")

AlcErrno	AlcDouble3Calloc(double ****dest, int mElem, int nElem,
				int oElem)
ALC_TEMPLATE_C3D(dest, double, mElem, nElem, oElem, "AlcDouble3Calloc")

/************************************************************************
* Functions:	AlcChar3Malloc, AlcShort3Malloc, AlcInt3Malloc,		*
*		AlcFloat3Malloc, AlcDouble3Malloc and AlcUnchar3Malloc.	*
* Return:	AlcErrno:		Allocation error number.	*
* Purpose:	A set of functions for allocating 3 dimensional non -	*
*		zero'd arrays of fundamental types.			*
*		All 3 dimensional arrays allocated by these functions	*
*		should be free'd using one of the associated free	*
*		functions, eg AlcChar3Malloc and AlcChar3Free.		*
* Global refs:	-							*
* Parameters:	<TYPE> ****dest:	Destination for allocated array	*
*					pointer.			*
*		int mElem:		Number of 2D arrays.		*
*		int nElem:		Number of 1D arrays.		*
*		int oElem:		Number of elements in each 1D	*
*					array.				*
************************************************************************/
AlcErrno	AlcChar3Malloc(char ****dest, int mElem, int nElem,
			       int oElem)
ALC_TEMPLATE_M3D(dest, char, mElem, nElem, oElem, "AlcChar3Malloc")

AlcErrno	AlcUnchar3Malloc(unsigned char ****dest,
			      int mElem, int nElem, int oElem)
{
  return(AlcChar3Malloc((char ****)dest, mElem, nElem, oElem));
}

AlcErrno	AlcShort3Malloc(short ****dest, int mElem, int nElem,
				int oElem)
ALC_TEMPLATE_M3D(dest, short, mElem, nElem, oElem, "AlcShort3Malloc")

AlcErrno	AlcInt3Malloc(int ****dest, int mElem, int nElem,
			      int oElem)
ALC_TEMPLATE_M3D(dest, int, mElem, nElem, oElem, "AlcInt3Malloc")

AlcErrno	AlcFloat3Malloc(float ****dest, int mElem, int nElem,
				int oElem)
ALC_TEMPLATE_M3D(dest, float, mElem, nElem, oElem, "AlcFloat3Malloc")

AlcErrno	AlcDouble3Malloc(double ****dest, int mElem, int nElem,
				 int oElem)
ALC_TEMPLATE_M3D(dest, double, mElem, nElem, oElem, "AlcDouble3Malloc")

/************************************************************************
* Functions:	Alc3Free, AlcChar3Free, AlcShort3Free, AlcInt3Free,	*
*		AlcFloat3Free, AlcDouble3Free and AlcUnchar3Free.	*
* Return:	AlcErrno:		Type allocation error number.	*
* Purpose:	A set of functions for freeing 3 dimensional arrays	*
*		allocated by the 3 dimensional array allocation 	*
*		functions.						*
* Global refs:	-							*
* Parameters:	<TYPE> ***dest:		Ptr with array to be free'd.	*
************************************************************************/
AlcErrno	Alc3Free(void ***dest)
ALC_TEMPLATE_F3D(dest, "Alc3Free")

AlcErrno	AlcChar3Free(char ***dest)
ALC_TEMPLATE_F3D(dest, "AlcChar3Free")

AlcErrno	AlcUnchar3Free(unsigned char ***dest)
{
  return(AlcChar3Free((char ***)dest));
}

AlcErrno	AlcShort3Free(short ***dest)
ALC_TEMPLATE_F3D(dest, "AlcShort3Free")

AlcErrno	AlcInt3Free(int ***dest)
ALC_TEMPLATE_F3D(dest, "AlcInt3Free")

AlcErrno	AlcFloat3Free(float ***dest)
ALC_TEMPLATE_F3D(dest, "AlcFloat3Free")

AlcErrno	AlcDouble3Free(double ***dest)
ALC_TEMPLATE_F3D(dest, "AlcDouble3Free")
