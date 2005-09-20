#pragma ident "MRC HGU $Id$"
/*!
* \file         libAlc/AlcArray.c
* \author       Bill Hill
* \date         April 2001
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
* \brief        Provides functions for the allocation of 1, 2 and 3D
*		arrays of types char, short, int, float and double.
*               Extension to other types (including user defined types)
*               should be straight formward through templates defined
*               in AlcTemplates.h.
* \ingroup	AlcArray
* \todo		-
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alc.h>
#include <AlcTemplates.h>

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd bit array.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcBit1Calloc(unsigned char **dest, size_t mElem)
{
  return(AlcUnchar1Calloc(dest, (mElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of pointers to void.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcPtr1Calloc(void ***dest, size_t mElem)
{
  AlcErrno	alcErrno = ALC_ER_NONE;
  
  /* Template doesn't work for pointer types. */
  if((dest) == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if(mElem < 1)
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if((*(dest) = (void **)AlcCalloc(mElem, sizeof(void *))) == NULL)
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno != ALC_ER_NONE)
  {
    if(dest)
    {
      *(dest) = NULL;
    }
  }
  return(alcErrno);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of chars.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcChar1Calloc(char **dest, size_t mElem)
{
  ALC_TEMPLATE_C1D(dest, char, mElem, "AlcChar1Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of unsigned chars.
* \note	        Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			 Destination for allocated array
*                                        pointer.
* \param        mElem 		         Number of elements in array.
*/
AlcErrno	AlcUnchar1Calloc(unsigned char **dest, size_t mElem)
{
  return(AlcChar1Calloc((char **)dest, mElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of shorts.
* \note	        Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*                                       pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcShort1Calloc(short **dest, size_t mElem)
{
  ALC_TEMPLATE_C1D(dest, short, mElem, "AlcShort1Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of ints.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcInt1Calloc(int **dest, size_t mElem)
{
  ALC_TEMPLATE_C1D(dest, int, mElem, "AlcInt1Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of floats.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcFloat1Calloc(float **dest, size_t mElem)
{
  ALC_TEMPLATE_C1D(dest, float, mElem, "AlcFloat1Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional zero'd array of doubles.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcDouble1Calloc(double **dest, size_t mElem)
{
  ALC_TEMPLATE_C1D(dest, double, mElem, "AlcDouble1Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd bit array.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcBit1Malloc(unsigned char **dest, size_t mElem)
{
  return(AlcUnchar1Malloc(dest, (mElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of pointers
*		to void.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcPtr1Malloc(void ***dest, size_t mElem)
{
  AlcErrno	alcErrno = ALC_ER_NONE;
  
  /* Template doesn't work for pointer types. */
  if(dest == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if(mElem < 1)
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if((*dest = (void **)AlcMalloc(mElem * sizeof(void *))) == NULL)
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno != ALC_ER_NONE)
  {
    if(dest)
    {
      *dest = NULL;
    }
  }
  return(alcErrno);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of chars.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcChar1Malloc(char **dest, size_t mElem)
{
  ALC_TEMPLATE_M1D(dest, char, mElem, "AlcChar1Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of unsigned chars.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 			 Destination for allocated array
*					pointer.
* \param	mElem 			Number of elements in array.
*/
AlcErrno	AlcUnchar1Malloc(unsigned char **dest, size_t mElem)
{
  return(AlcChar1Malloc((char **)dest, mElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of shorts.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcShort1Malloc(short **dest, size_t mElem)
{
  ALC_TEMPLATE_M1D(dest, short, mElem, "AlcShort1Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of ints.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcInt1Malloc(int **dest, size_t mElem)
{
  ALC_TEMPLATE_M1D(dest, int, mElem, "AlcInt1Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of floats.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcFloat1Malloc(float **dest, size_t mElem)
{
  ALC_TEMPLATE_M1D(dest, float, mElem, "AlcFloat1Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 1 dimensional non-zero'd array of doubles.
* \note		Should be free'd using AlcFree().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of elements in array.
*/
AlcErrno	AlcDouble1Malloc(double **dest, size_t mElem)
{
  ALC_TEMPLATE_M1D(dest, double, mElem, "AlcDouble1Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd bit array.
* \note         Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit2Calloc(unsigned char ***dest,
			      size_t mElem, size_t nElem)
{
  return(AlcUnchar2Calloc(dest, mElem, (nElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of pointers to void.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcPtr2Calloc(void ****dest, size_t mElem, size_t nElem)
{
  size_t	index;
  void		**dump0 = NULL;
  void		***dump1 = NULL;
  AlcErrno	alcErrno = ALC_ER_NONE;
 
  /* Template doesn't work for pointer types. */
  if(dest == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if((mElem < 1) || (nElem < 1))
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if(((dump0 = (void **)AlcCalloc(mElem * nElem,
  				       sizeof(void *))) == NULL) ||
          ((dump1 = (void ***)AlcMalloc(mElem * sizeof(void **))) == NULL))
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno == ALC_ER_NONE)
  {
    *dest = dump1;
    for(index = 0; index < mElem; ++index)
    {
      (*dest)[index] = dump0;
      dump0 += nElem;
    }
  }
  else
  {
    if(dest)
    {
      *dest = NULL;
    }
    if(dump0)
    {
      AlcFree(dump0);
    }
    if(dump1)
    {
      AlcFree(dump1);
    }
  }
  return(alcErrno);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of chars.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 	 		Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar2Calloc(char ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_C2D(dest, char, mElem, nElem, "AlcChar2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of unsigned chars.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar2Calloc(unsigned char ***dest,
				 size_t mElem, size_t nElem)
{
  return(AlcChar2Calloc((char ***)dest, mElem, nElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of shorts.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort2Calloc(short ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_C2D(dest, short, mElem, nElem, "AlcShort2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of ints.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt2Calloc(int ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_C2D(dest, int, mElem, nElem, "AlcInt2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of floats.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat2Calloc(float ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_C2D(dest, float, mElem, nElem, "AlcFloat2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd array of doubles.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble2Calloc(double ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_C2D(dest, double, mElem, nElem, "AlcDouble2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd bit array.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			 Destination for allocated array
*					  pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit2Malloc(unsigned char ***dest,
			      size_t mElem, size_t nElem)
{
  return(AlcUnchar2Malloc(dest, mElem, (nElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of pointers
*		to void.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcPtr2Malloc(void ****dest, size_t mElem, size_t nElem)
{
  size_t	index;
  void 		**dump0 = NULL;
  void   	***dump1 = NULL;
  AlcErrno	alcErrno = ALC_ER_NONE;
 
  /* Template doesn't work for pointer types. */
  if(dest == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if((mElem < 1) || (nElem < 1))
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if(((dump0 = (void **)AlcMalloc(mElem * nElem *
  				       sizeof(void *))) == NULL) ||
          ((dump1 = (void ***)AlcMalloc(mElem * sizeof(void **))) == NULL))
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno == ALC_ER_NONE)
  {
    *dest = dump1;
    for(index = 0; index < mElem; ++index)
    {
      (*dest)[index] = dump0;
      dump0 += nElem;
    }
  }
  else
  {
    if(dest)
    {
      *dest = NULL;
    }
    if(dump0)
    {
      AlcFree(dump0);
    }
    if(dump1)
    {
      AlcFree(dump1);
    }
  }
  return(alcErrno);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of chars.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 		 	Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar2Malloc(char ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_M2D(dest, char, mElem, nElem, "AlcChar2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of unsigned chars.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar2Malloc(unsigned char ***dest,
				 size_t mElem, size_t nElem)
{
  return(AlcChar2Malloc((char ***)dest, mElem, nElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of shorts.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort2Malloc(short ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_M2D(dest, short, mElem, nElem, "AlcShort2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of ints.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt2Malloc(int ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_M2D(dest, int, mElem, nElem, "AlcInt2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of floats.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat2Malloc(float ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_M2D(dest, float, mElem, nElem, "AlcFloat2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd array of doubles.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 1D arrays.
* \param	nElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble2Malloc(double ***dest, size_t mElem, size_t nElem)
{
  ALC_TEMPLATE_M2D(dest, double, mElem, nElem, "AlcDouble2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of chars
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymChar2Calloc(char ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_C2D(dest, char, nElem, "AlcSymChar2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of
*		unsigned chars in which only the lower triangle of
*		elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymUnchar2Calloc(unsigned char ***dest, size_t nElem)
{
  return(AlcSymChar2Calloc((char ***)dest, nElem));
}


/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of shorts
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymShort2Calloc(short ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_C2D(dest, short, nElem, "AlcSymShort2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of ints
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymInt2Calloc(int ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_C2D(dest, int, nElem, "AlcSymInt2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of floats
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymFloat2Calloc(float ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_C2D(dest, float, nElem, "AlcSymFloat2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional zero'd symetric array of doubles
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymDouble2Calloc(double ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_C2D(dest, double, nElem, "AlcSymDouble2Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of chars
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymChar2Malloc(char ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_M2D(dest, char, nElem, "AlcSymChar2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of
*		unsigned chars in which only the lower triangle of
*		elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymUnchar2Malloc(unsigned char ***dest, size_t nElem)
{
  return(AlcSymChar2Malloc((char ***)dest, nElem));
}


/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of shorts
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymShort2Malloc(short ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_M2D(dest, short, nElem, "AlcSymShort2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of ints
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymInt2Malloc(int ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_M2D(dest, int, nElem, "AlcSymInt2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of floats
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymFloat2Malloc(float ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_M2D(dest, float, nElem, "AlcSymFloat2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 2 dimensional non-zero'd symetric array of doubles
*		in which only the lower triangle of elements are stored.
* \note		Should be free'd using Alc2Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	nElem 	 		Number of row or column elements in
*					the symetric array.
*/
AlcErrno	AlcSymDouble2Malloc(double ***dest, size_t nElem)
{
  ALC_TEMPLATE_SYM_M2D(dest, double, nElem, "AlcSymDouble2Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Free's a 2 dimensional array allocated by one of the
*		2 dimensional array allocation functions.
* \param	dat 			 Ptr with array to be free'd.
*/
AlcErrno	Alc2Free(void **dat)
{
  ALC_TEMPLATE_F2D(dat, "Alc2Free")
}

AlcErrno	AlcBit2Free(unsigned char **dat)
{
  return(AlcUnchar2Free(dat));
}

AlcErrno	AlcChar2Free(char **dat)
{
  ALC_TEMPLATE_F2D(dat, "AlcChar2Free")
}

AlcErrno	AlcUnchar2Free(unsigned char **dat)
{
  return(AlcChar2Free((char **)dat));
}

AlcErrno	AlcShort2Free(short **dat)
{
  ALC_TEMPLATE_F2D(dat, "AlcShort2Free")
}

AlcErrno	AlcInt2Free(int **dat)
{
  ALC_TEMPLATE_F2D(dat, "AlcInt2Free")
}

AlcErrno	AlcFloat2Free(float **dat)
{
  ALC_TEMPLATE_F2D(dat, "AlcFloat2Free")
}

AlcErrno	AlcDouble2Free(double **dat)
{
  ALC_TEMPLATE_F2D(dat, "AlcDouble2Free")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional zero'd bit array.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit3Calloc(unsigned char ****dest,
			      size_t mElem, size_t nElem, size_t oElem)
{
  return(AlcUnchar3Calloc(dest, mElem, nElem, (oElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of pointers to void.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcPtr3Calloc(void *****dest, size_t mElem, size_t nElem,
			       size_t oElem)
{
  size_t	index0,
  		index1;
  void		**dump0 = NULL,
  		***dump1 = NULL,
		****dump2 = NULL;
  AlcErrno	alcErrno = ALC_ER_NONE;

  if((dest) == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if((mElem < 1) || (nElem < 1) || (oElem < 1))
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if(((dump0 = (void **)AlcCalloc(mElem * nElem * oElem,
  				       sizeof(void *))) == NULL) ||
          ((dump1 = (void ***)AlcMalloc(mElem * nElem *
	                                sizeof(void **))) == NULL) ||
          ((dump2 = (void ****)AlcMalloc(mElem * sizeof(void ***))) == NULL))
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno == ALC_ER_NONE)
  {
    *(dest) = dump2;
    for(index0 = 0; index0 < mElem; ++index0)
    {
      for(index1=0; index1 < nElem; ++index1)
      {
	dump1[index1] = dump0;
	dump0 += oElem;
      }
      (*(dest))[index0] = dump1;
      dump1 += nElem;
    }
  }
  else
  {
    if(dest)
    {
      *(dest) = NULL;
    }
    AlcFree(dump2);
    AlcFree(dump1);
    AlcFree(dump0);
  }
  return(alcErrno);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar3Calloc(char ****dest, size_t mElem, size_t nElem,
			       size_t oElem)
{
  ALC_TEMPLATE_C3D(dest, char, mElem, nElem, oElem, "AlcChar3Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of unsigned chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar3Calloc(unsigned char ****dest,
			   size_t mElem, size_t nElem, size_t oElem)
{
  return(AlcChar3Calloc((char ****)dest, mElem, nElem, oElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort3Calloc(short ****dest, size_t mElem, size_t nElem,
			       size_t oElem)
{
  ALC_TEMPLATE_C3D(dest, short, mElem, nElem, oElem, "AlcShort3Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt3Calloc(int ****dest, size_t mElem, size_t nElem,
			      size_t oElem)
{
  ALC_TEMPLATE_C3D(dest, int, mElem, nElem, oElem, "AlcInt3Calloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional array of chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat3Calloc(float ****dest, size_t mElem, size_t nElem,
			        size_t oElem)
{
  ALC_TEMPLATE_C3D(dest, float, mElem, nElem, oElem, "AlcFloat3Calloc")
}

/*!
* \return	Error code.
* \brief	Allocates a 3 dimensional zero'd array of doubles.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble3Calloc(double ****dest, size_t mElem, size_t nElem,
				 size_t oElem)
{
  ALC_TEMPLATE_C3D(dest, double, mElem, nElem, oElem, "AlcDouble3Calloc")
}

/*!
* \return	Error code.
* \brief	Allocates a 3 dimensional non-zero'd bit array.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcBit3Malloc(unsigned char ****dest,
			      size_t mElem, size_t nElem, size_t oElem)
{
  return(AlcUnchar3Malloc(dest, mElem, nElem, (oElem + 7) / 8));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of pointers
*		to void.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcPtr3Malloc(void *****dest, size_t mElem, size_t nElem,
			       size_t oElem)
{
  size_t	index0,
  		index1;
  void 		**dump0 = NULL,
  		***dump1 = NULL,
		****dump2 = NULL;
  AlcErrno	alcErrno = ALC_ER_NONE;

  if(dest == NULL)
  {
    alcErrno = ALC_ER_NULLPTR;
  }
  else if((mElem < 1) || (nElem < 1) || (oElem < 1))
  {
    alcErrno = ALC_ER_NUMELEM;
  }
  else if(((dump0 = (void **)AlcMalloc(mElem * nElem * oElem *
  				       sizeof(void *))) == NULL) ||
          ((dump1 = (void ***)AlcMalloc(mElem * nElem *
	                                sizeof(void **))) == NULL) ||
          ((dump2 = (void ****)AlcMalloc(mElem * sizeof(void ***))) == NULL))
  {
    alcErrno = ALC_ER_ALLOC;
  }
  if(alcErrno == ALC_ER_NONE)
  {
    *dest = dump2;
    for(index0 = 0; index0 < mElem; ++index0)
    {
      for(index1=0; index1 < nElem; ++index1)
      {
	dump1[index1] = dump0;
	dump0 += oElem;
      }
      (*dest)[index0] = dump1;
      dump1 += nElem;
    }
  }
  else
  {
    if(dest)
    {
      *dest = NULL;
    }
    AlcFree(dump2);
    AlcFree(dump1);
    AlcFree(dump0);
  }
  return(alcErrno); 
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcChar3Malloc(char ****dest, size_t mElem, size_t nElem,
			       size_t oElem)
{
  ALC_TEMPLATE_M3D(dest, char, mElem, nElem, oElem, "AlcChar3Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of unsigned chars.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			 Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcUnchar3Malloc(unsigned char ****dest,
  			         size_t mElem, size_t nElem, size_t oElem)
{
  return(AlcChar3Malloc((char ****)dest, mElem, nElem, oElem));
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of shorts.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcShort3Malloc(short ****dest, size_t mElem, size_t nElem,
				size_t oElem)
{
  ALC_TEMPLATE_M3D(dest, short, mElem, nElem, oElem, "AlcShort3Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of ints.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcInt3Malloc(int ****dest, size_t mElem, size_t nElem,
			      size_t oElem)
{
  ALC_TEMPLATE_M3D(dest, int, mElem, nElem, oElem, "AlcInt3Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of floats.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcFloat3Malloc(float ****dest, size_t mElem, size_t nElem,
				size_t oElem)
{
  ALC_TEMPLATE_M3D(dest, float, mElem, nElem, oElem, "AlcFloat3Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Allocates a 3 dimensional non-zero'd array of doubles.
* \note		Should be free'd using Alc3Free().
* \note		Array size is limited only by address space.
* \param	dest 			Destination for allocated array
*					pointer.
* \param	mElem 	 		Number of 2D arrays.
* \param	nElem 	 		Number of 1D arrays.
* \param	oElem 	 		Number of elements in each 1D
*					array.
*/
AlcErrno	AlcDouble3Malloc(double ****dest, size_t mElem, size_t nElem,
				 size_t oElem)
{
  ALC_TEMPLATE_M3D(dest, double, mElem, nElem, oElem, "AlcDouble3Malloc")
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Free's any 3 dimensional array allocated by one of the
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
* \return	Error code.
* \ingroup	AlcArray
* \brief	Reads a 1D double array from the given numeric ASCI file.
*		Each value should be on a seperate line.
* \param	fP:			File pointer.
* \param	dstA 			Destination pointer for the new
*					array.
* \param	dstNElem		Destination pointer for the number
*					of elements in the 1D array.
*/
AlcErrno	AlcDouble1ReadAsci(FILE *fP, double **dstA,
				   size_t *dstNElem)
{
  size_t	nR = 0;
  double	*dP0,
  		*aM = NULL;
  AlcVector	*vec = NULL;
  const size_t	vecCnt = 1024;		/* Initial number of elements in
  					 * the vector */
  const int	maxRecLen = 8192;	/* Maximum number of chars in an
  					 * input record */
  char		recS[8192 /* = maxRecLen */];
  AlcErrno	errNum = ALC_ER_NONE;

  vec = AlcVectorNew(vecCnt, sizeof(double), vecCnt, &errNum);
  while((errNum == ALC_ER_NONE) && (fgets(recS, maxRecLen, fP) != NULL))
  {
    if((dP0 = (double *)AlcVectorExtendAndGet(vec, nR)) == NULL)
    {
      errNum = ALC_ER_ALLOC;
    }
    else if(sscanf(recS, "%lg", dP0) != 1)
    {
      errNum = ALC_ER_READ;
    }
    ++nR;
  }
  if(errNum == ALC_ER_NONE)
  {
    aM = (double *)AlcVectorToArray1D(vec, 0, nR - 1, &errNum);
  }
  if(errNum == ALC_ER_NONE)
  {
    *dstA = aM;
    *dstNElem = nR;
  }
  (void )AlcVectorFree(vec);
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Reads a 2D double array from the given numeric ASCI file.
*		Fields in the file must be white space saperated and
*		records must be on separate lines. The number of fields
*		per record must be the same for all records.
* \param	fP:			File pointer.
* \param	dstA 			Destination pointer for the new
*					array.
* \param	dstMElem		Destination pointer for the number of
*					1D arrays (number of records).
* \param	dstNElem		Destination pointer for the number of
*					elements in each 1D array (number of
*					fields per record).
*/
AlcErrno	AlcDouble2ReadAsci(FILE *fP, double ***dstA,
				   size_t *dstMElem, size_t *dstNElem)
{
  size_t	iF,
  		nR = 0,
  		nF = 0,
		nV = 0;
  double	*dP0;
  double	**aM = NULL;
  char		*parseS,
  		*tokS;
  const size_t	vecCnt = 1024; 		/* Initial number of elements in
  					 * the vector */
  const int	maxRecLen = 8192; 	/* Maximum number of chars in an
  					 * input record */
  char		recS[8192 /* = maxRecLen */];
  AlcVector	*vec = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  vec = AlcVectorNew(vecCnt, sizeof(double), vecCnt, &errNum);
  while((errNum == ALC_ER_NONE) && (fgets(recS, maxRecLen, fP) != NULL))
  {
    iF = 0;
    parseS = recS;
    while((errNum == ALC_ER_NONE) &&
          ((tokS = (char *) strtok(parseS, " \t\n")) != NULL) && *tokS)
    {
      parseS = NULL;
      if((dP0 = (double *)AlcVectorExtendAndGet(vec, nV)) == NULL)
      {
        errNum = ALC_ER_ALLOC;
      }
      else if(sscanf(tokS, "%lg", dP0) != 1)
      {
        errNum = ALC_ER_READ;
      }
      else
      {
        ++iF;
	++nV;
      }
    }
    if((errNum == ALC_ER_NONE) && (iF > 0))
    {
      if(nR == 0)
      {
	 nF = iF;
      }
      else if(iF != nF)
      {
	errNum = ALC_ER_READ;
      }
      ++nR;
    }
  }
  if(errNum == ALC_ER_NONE)
  {
    aM = (double **)AlcVectorToArray2D(vec, 0, nV - 1, nR, nF, &errNum);
  }
  if(errNum == ALC_ER_NONE)
  {
    *dstA = aM;
    *dstNElem = nF;
    *dstMElem = nR;
  }
  (void )AlcVectorFree(vec);
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Writes a 1D double array in numeric ASCI format to the
*		given file file. Elements are on separate lines.
* \param	fP:			File pointer.
* \param	ar 			Given array.
* \param	nElem			Number of elements in the 1D array.
*/
AlcErrno	AlcDouble1WriteAsci(FILE *fP, double *ar,
				    size_t nElem)
{
  size_t	iR;
  AlcErrno	errNum = ALC_ER_NONE;

  iR = 0;
  while((iR < nElem) && (errNum == ALC_ER_NONE))
  {
    (void )fprintf(fP, "%lg ", ar[iR]);
    if(fprintf(fP, "\n") != 1)
    {
      errNum = ALC_ER_WRITE;
    }
    ++iR;
  }
  return(errNum);
}

/*!
* \return	Error code.
* \ingroup	AlcArray
* \brief	Writes a 2D double array in numeric ASCI format to the
*		given file file.
*		Fields in the file are white space seperated and
*		records are on separate lines.
* \param	fP:			File pointer.
* \param	ar 			Given array.
* \param	mElem			Number of 1D arrays (number of
*					records).
* \param	nElem			Number of elements in each 1D array
*					(number of fields per record).
*/
AlcErrno	AlcDouble2WriteAsci(FILE *fP, double **ar,
				    size_t mElem, size_t nElem)
{
  size_t	iC,
		iR;
  AlcErrno	errNum = ALC_ER_NONE;

  iR = 0;
  while((iR < mElem) && (errNum == ALC_ER_NONE))
  {
    for(iC = 0; iC < nElem; ++iC)
    {
      (void )fprintf(fP, "%lg ", ar[iR][iC]);
    }
    if(fprintf(fP, "\n") != 1)
    {
      errNum = ALC_ER_WRITE;
    }
    ++iR;
  }
  return(errNum);
}
