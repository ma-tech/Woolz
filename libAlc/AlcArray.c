#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcArray.c
* \author       Bill Hill
* \date         April 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides functions for the allocation of 1, 2 and 3D
*		arrays of types char, short, int, float and double.
*               Extension to other types (including user defined types)
*               should be straight formward through templates defined
*               in AlcTemplates.h.
* \todo		-
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>
#include <AlcTemplates.h>

/*!
* \ingroup	Alc
* \defgroup	AlcArray
* @{
*/

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd bit array.
*		Should be free'd using AlcFree().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcBit1Calloc(unsigned char **dest, int mElem)
{
  return(AlcUnchar1Calloc(dest, (mElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd array of chars.
*		Should be free'd using AlcFree().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcChar1Calloc(char **dest, int mElem)
{
  ALC_TEMPLATE_C1D(dest, char, mElem, "AlcChar1Calloc")
}

/*!
* \return	 	 	 	Error code.
* \brief	Allocates a 1 dimensional zero'd array of unsigned chars.
*               Should be free'd using AlcFree().
* \param	dest 			 Destination for allocated array
*                                        pointer.
* \param        mElem 		         Number of elements in array.
*/
AlcErrno	AlcUnchar1Calloc(unsigned char **dest, int mElem)
{
  return(AlcChar1Calloc((char **)dest, mElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd array of shorts.
*               Should be free'd using AlcFree().
* \param	dest 			Destination for allocated array
*                                       pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcShort1Calloc(short **dest, int mElem)
{
  ALC_TEMPLATE_C1D(dest, short, mElem, "AlcShort1Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd array of ints.
*               Should be free'd using AlcFree().
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcInt1Calloc(int **dest, int mElem)
{
  ALC_TEMPLATE_C1D(dest, int, mElem, "AlcInt1Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd array of floats.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcFloat1Calloc(float **dest, int mElem)
{
  ALC_TEMPLATE_C1D(dest, float, mElem, "AlcFloat1Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional zero'd array of doubles.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcDouble1Calloc(double **dest, int mElem)
{
  ALC_TEMPLATE_C1D(dest, double, mElem, "AlcDouble1Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd bit array.
*               Should be free'd using AlcFree().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcBit1Malloc(unsigned char **dest, int mElem)
{
  return(AlcUnchar1Malloc(dest, (mElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of chars.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcChar1Malloc(char **dest, int mElem)
{
  ALC_TEMPLATE_M1D(dest, char, mElem, "AlcChar1Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of unsigned chars.
*               Should be free'd using AlcFree().
* \param	dest 			 Destination for allocated array
*					pointer.
* \param	mElem 			Number of elements in array.
*/
AlcErrno	AlcUnchar1Malloc(unsigned char **dest, int mElem)
{
  return(AlcChar1Malloc((char **)dest, mElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of shorts.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcShort1Malloc(short **dest, int mElem)
{
  ALC_TEMPLATE_M1D(dest, short, mElem, "AlcShort1Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of ints.
*               Should be free'd using AlcFree().
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcInt1Malloc(int **dest, int mElem)
{
  ALC_TEMPLATE_M1D(dest, int, mElem, "AlcInt1Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of floats.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcFloat1Malloc(float **dest, int mElem)
{
  ALC_TEMPLATE_M1D(dest, float, mElem, "AlcFloat1Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 1 dimensional non-zero'd array of doubles.
*               Should be free'd using AlcFree().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcDouble1Malloc(double **dest, int mElem)
{
  ALC_TEMPLATE_M1D(dest, double, mElem, "AlcDouble1Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd bit array.
*               Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit2Calloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcUnchar2Calloc(dest, mElem, (nElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of chars.
*               Should be free'd using Alc2Free().
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar2Calloc(char ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_C2D(dest, char, mElem, nElem, "AlcChar2Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of unsigned chars.
*               Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar2Calloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcChar2Calloc((char ***)dest, mElem, nElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of shorts.
*               Should be free'd using Alc2Free().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort2Calloc(short ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_C2D(dest, short, mElem, nElem, "AlcShort2Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of ints.
*               Should be free'd using Alc2Free().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt2Calloc(int ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_C2D(dest, int, mElem, nElem, "AlcInt2Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of floats.
*               Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat2Calloc(float ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_C2D(dest, float, mElem, nElem, "AlcFloat2Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional zero'd array of doubles.
*               Should be free'd using Alc2Free().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble2Calloc(double ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_C2D(dest, double, mElem, nElem, "AlcDouble2Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd bit array.
*		Should be free'd using Alc2Free().
* \param	dest 			 Destination for allocated array
*					  pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit2Malloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcUnchar2Malloc(dest, mElem, (nElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of chars.
*		Should be free'd using Alc2Free().
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar2Malloc(char ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_M2D(dest, char, mElem, nElem, "AlcChar2Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of unsigned chars.
*		Should be free'd using Alc2Free().
* \param	dest 			 Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar2Malloc(unsigned char ***dest, int mElem, int nElem)
{
  return(AlcChar2Malloc((char ***)dest, mElem, nElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of shorts.
*		Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort2Malloc(short ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_M2D(dest, short, mElem, nElem, "AlcShort2Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of ints.
*		Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt2Malloc(int ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_M2D(dest, int, mElem, nElem, "AlcInt2Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of floats.
*		Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat2Malloc(float ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_M2D(dest, float, mElem, nElem, "AlcFloat2Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 2 dimensional non-zero'd array of doubles.
*		Should be free'd using Alc2Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble2Malloc(double ***dest, int mElem, int nElem)
{
  ALC_TEMPLATE_M2D(dest, double, mElem, nElem, "AlcDouble2Malloc")
}

/*!
* \return		 		Error code.
* \brief	Free's a 2 dimensional array allocated by one of the
*		2 dimensional array allocation functions.
* \param	dest 			 Ptr with array to be free'd.
*/
AlcErrno	Alc2Free(void **dest)
{
  ALC_TEMPLATE_F2D(dest, "Alc2Free")
}

AlcErrno	AlcBit2Free(unsigned char **dest)
{
  return(AlcUnchar2Free(dest));
}

AlcErrno	AlcChar2Free(char **dest)
{
  ALC_TEMPLATE_F2D(dest, "AlcChar2Free")
}

AlcErrno	AlcUnchar2Free(unsigned char **dest)
{
  return(AlcChar2Free((char **)dest));
}

AlcErrno	AlcShort2Free(short **dest)
{
  ALC_TEMPLATE_F2D(dest, "AlcShort2Free")
}

AlcErrno	AlcInt2Free(int **dest)
{
  ALC_TEMPLATE_F2D(dest, "AlcInt2Free")
}

AlcErrno	AlcFloat2Free(float **dest)
{
  ALC_TEMPLATE_F2D(dest, "AlcFloat2Free")
}

AlcErrno	AlcDouble2Free(double **dest)
{
  ALC_TEMPLATE_F2D(dest, "AlcDouble2Free")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional zero'd bit array.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit3Calloc(unsigned char ****dest,
			      int mElem, int nElem, int oElem)
{
  return(AlcUnchar3Calloc(dest, mElem, nElem, (oElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional array of chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar3Calloc(char ****dest, int mElem, int nElem,
			      int oElem)
{
  ALC_TEMPLATE_C3D(dest, char, mElem, nElem, oElem, "AlcChar3Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional array of unsigned chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar3Calloc(unsigned char ****dest,
			   int mElem, int nElem, int oElem)
{
  return(AlcChar3Calloc((char ****)dest, mElem, nElem, oElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional array of chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort3Calloc(short ****dest, int mElem, int nElem,
			       int oElem)
{
  ALC_TEMPLATE_C3D(dest, short, mElem, nElem, oElem, "AlcShort3Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional array of chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt3Calloc(int ****dest, int mElem, int nElem, int oElem)
{
  ALC_TEMPLATE_C3D(dest, int, mElem, nElem, oElem, "AlcInt3Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional array of chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat3Calloc(float ****dest, int mElem, int nElem,
			       int oElem)
{
  ALC_TEMPLATE_C3D(dest, float, mElem, nElem, oElem, "AlcFloat3Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional zero'd array of doubles.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble3Calloc(double ****dest, int mElem, int nElem,
				int oElem)
{
  ALC_TEMPLATE_C3D(dest, double, mElem, nElem, oElem, "AlcDouble3Calloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd bit array.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit3Malloc(unsigned char ****dest,
			      int mElem, int nElem, int oElem)
{
  return(AlcUnchar3Malloc(dest, mElem, nElem, (oElem + 7) / 8));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of chars.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar3Malloc(char ****dest, int mElem, int nElem,
			       int oElem)
{
  ALC_TEMPLATE_M3D(dest, char, mElem, nElem, oElem, "AlcChar3Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of unsigned chars.
*		Should be free'd using Alc3Free().
* \param	dest 			 Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar3Malloc(unsigned char ****dest,
  			         int mElem, int nElem, int oElem)
{
  return(AlcChar3Malloc((char ****)dest, mElem, nElem, oElem));
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of shorts.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort3Malloc(short ****dest, int mElem, int nElem,
				int oElem)
{
  ALC_TEMPLATE_M3D(dest, short, mElem, nElem, oElem, "AlcShort3Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of ints.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt3Malloc(int ****dest, int mElem, int nElem,
			      int oElem)
{
  ALC_TEMPLATE_M3D(dest, int, mElem, nElem, oElem, "AlcInt3Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of floats.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat3Malloc(float ****dest, int mElem, int nElem,
				int oElem)
{
  ALC_TEMPLATE_M3D(dest, float, mElem, nElem, oElem, "AlcFloat3Malloc")
}

/*!
* \return		 		Error code.
* \brief	Allocates a 3 dimensional non-zero'd array of doubles.
*		Should be free'd using Alc3Free().
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble3Malloc(double ****dest, int mElem, int nElem,
				 int oElem)
{
  ALC_TEMPLATE_M3D(dest, double, mElem, nElem, oElem, "AlcDouble3Malloc")
}

/*!
* \return		 		Error code.
* \brief	Free's a 3 dimensional array allocated by one of the
*		3 dimensional array allocation functions.
* \param	dest 			 Ptr with array to be free'd.
*/
AlcErrno	Alc3Free(void ***dest)
{
  ALC_TEMPLATE_F3D(dest, "Alc3Free")
}

AlcErrno	AlcBit3Free(unsigned char ***dest)
{
  return(AlcUnchar3Free(dest));
}

AlcErrno	AlcChar3Free(char ***dest)
{
  ALC_TEMPLATE_F3D(dest, "AlcChar3Free")
}

AlcErrno	AlcUnchar3Free(unsigned char ***dest)
{
  return(AlcChar3Free((char ***)dest));
}

AlcErrno	AlcShort3Free(short ***dest)
{
  ALC_TEMPLATE_F3D(dest, "AlcShort3Free")
}

AlcErrno	AlcInt3Free(int ***dest)
{
  ALC_TEMPLATE_F3D(dest, "AlcInt3Free")
}

AlcErrno	AlcFloat3Free(float ***dest)
{
  ALC_TEMPLATE_F3D(dest, "AlcFloat3Free")
}

AlcErrno	AlcDouble3Free(double ***dest)
{
  ALC_TEMPLATE_F3D(dest, "AlcDouble3Free")
}

/*!
* @}
*/
