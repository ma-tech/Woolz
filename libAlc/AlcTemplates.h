#ifndef ALCTEMPLATES_H
#define ALCTEMPLATES_H
#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcTemplates.h
* Date:         March 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Templates used by the 'C' pre-processor to generate the
*		body of the MRC HGU allocation functions and the
*		associated freeing functions.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************
* Template:	ALC_TEMPLATE_C1D					*
* Purpose:	Allocates 1 dimensional zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T *).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of elements in array.	*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_C1D(D,T,M,F) \
{ \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if((M) < 1) \
    alcErrno = ALC_ER_NUMELEM; \
  else if((*(D) = (T *)AlcCalloc((M), sizeof(T))) == NULL) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno != ALC_ER_NONE) \
  { \
    if(D) \
      *(D) = NULL; \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_M1D					*
* Purpose:	Allocates 1 dimensional non-zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T *).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of elements in array.	*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_M1D(D,T,M,F) \
{ \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if((M) < 1) \
    alcErrno = ALC_ER_NUMELEM; \
  else if((*(D) = (T *)AlcMalloc((M) * sizeof(T))) == NULL) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno != ALC_ER_NONE) \
  { \
    if(D) \
      *(D) = NULL; \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_C2D					*
* Purpose:	Allocates 2 dimensional zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T **).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of 1D arrays.		*
*		N:			Number of elements in each 1D	*
*					array.				*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_C2D(D,T,M,N,F) \
{ \
  int		index; \
  T		*dump0 = NULL; \
  T		**dump1 = NULL; \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if(((M) < 1) || ((N) < 1)) \
    alcErrno = ALC_ER_NUMELEM; \
  else if(((dump0 = (T *)AlcCalloc((M) * (N), sizeof(T))) == NULL) || \
          ((dump1 = (T **)AlcMalloc((M) * sizeof(T *))) == NULL)) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno == ALC_ER_NONE) \
  { \
    *(D) = dump1; \
    for(index = 0; index < (M); ++index) \
    { \
      (*(D))[index] = dump0; \
      dump0 += (N); \
    } \
  } \
  else \
  { \
    if(D) \
      *(D) = NULL; \
    if(dump0) \
      AlcFree(dump0); \
    if(dump1) \
      AlcFree(dump1); \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_M2D					*
* Purpose:	Allocates 2 dimensional non-zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T **).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of 1D arrays.		*
*		N:			Number of elements in each 1D	*
*					array.				*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_M2D(D,T,M,N,F) \
{ \
  int		index; \
  T		*dump0 = NULL; \
  T  		**dump1 = NULL; \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if(((M) < 1) || ((N) < 1)) \
    alcErrno = ALC_ER_NUMELEM; \
  else if(((dump0 = (T *)AlcMalloc((M) * (N) * sizeof(T))) == NULL) || \
          ((dump1 = (T **)AlcMalloc((M) * sizeof(T *))) == NULL)) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno == ALC_ER_NONE) \
  { \
    *(D) = dump1; \
    for(index = 0; index < (M); ++index) \
    { \
      (*(D))[index] = dump0; \
      dump0 += (N); \
    } \
  } \
  else \
  { \
    if(D) \
      *(D) = NULL; \
    if(dump0) \
      AlcFree(dump0); \
    if(dump1) \
      AlcFree(dump1); \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_F2D					*
* Purpose:	Free's 2 dimensional arrays of any type, actualy no	*
*		type information is used in freeing the array.		*
* Global refs:	-							*
* Parameters:	D:			Pointer for array to be free'd.	*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_F2D(D,F) \
{ \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D == NULL) || (*(D) == NULL)) \
  { \
    alcErrno = ALC_ER_NULLPTR; \
  } \
  else \
  { \
    AlcFree(*(D)); \
    AlcFree(D); \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_C3D					*
* Purpose:	Allocates 3 dimensional zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T **).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of 2D arrays.		*
*		N:			Number of 1D arrays.		*
*		O:			Number of elements in each 1D	*
*					array.				*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_C3D(D,T,M,N,O,F) \
{ \
  int		index0, \
  		index1; \
  T		*dump0 = NULL, \
  		**dump1 = NULL, \
		***dump2 = NULL; \
  AlcErrno	alcErrno = ALC_ER_NONE; \
 \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if(((M) < 1) || ((N) < 1) || ((O) < 1)) \
    alcErrno = ALC_ER_NUMELEM; \
  else if(((dump0 = (T *)AlcCalloc((M) * (N) * (O), sizeof(T))) == NULL) || \
          ((dump1 = (T **)AlcMalloc((M) * (N) * sizeof(T *))) == NULL) || \
          ((dump2 = (T ***)AlcMalloc((M) * sizeof(T **))) == NULL)) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno == ALC_ER_NONE) \
  { \
    *(D) = dump2; \
    for(index0 = 0; index0 < (M); ++index0) \
    { \
      for(index1=0; index1 < (N); ++index1) \
      { \
	dump1[index1] = dump0; \
	dump0 += (O); \
      } \
      (*(D))[index0] = dump1; \
      dump1 += (N); \
    } \
  } \
  else \
  { \
    if(D) \
      *(D) = NULL; \
    if(dump2) \
      AlcFree(dump2); \
    if(dump1) \
      AlcFree(dump1); \
    if(dump0) \
      AlcFree(dump0); \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_M3D					*
* Purpose:	Allocates 3 dimensional non-zero'd arrays of any type.	*
* Global refs:	-							*
* Parameters:	D:			Destination pointer, of type	*
*					ie (T **).			*
*		T:			Type, eg char, short, int, ....	*
*		M:			Number of 2D arrays.		*
*		N:			Number of 1D arrays.		*
*		O:			Number of elements in each 1D	*
*					array.				*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_M3D(D,T,M,N,O,F) \
{ \
  int		index0, \
  		index1; \
  T		*dump0 = NULL, \
  		**dump1 = NULL, \
		***dump2 = NULL; \
  AlcErrno	alcErrno = ALC_ER_NONE; \
 \
  if((D) == NULL) \
    alcErrno = ALC_ER_NULLPTR; \
  else if(((M) < 1) || ((N) < 1) || ((O) < 1)) \
    alcErrno = ALC_ER_NUMELEM; \
  else if(((dump0 = (T *)AlcMalloc((M) * (N) * (O) * sizeof(T))) == NULL) || \
          ((dump1 = (T **)AlcMalloc((M) * (N) * sizeof(T *))) == NULL) || \
          ((dump2 = (T ***)AlcMalloc((M) * sizeof(T **))) == NULL)) \
    alcErrno = ALC_ER_ALLOC; \
  if(alcErrno == ALC_ER_NONE) \
  { \
    *(D) = dump2; \
    for(index0 = 0; index0 < (M); ++index0) \
    { \
      for(index1=0; index1 < (N); ++index1) \
      { \
	dump1[index1] = dump0; \
	dump0 += (O); \
      } \
      (*(D))[index0] = dump1; \
      dump1 += (N); \
    } \
  } \
  else \
  { \
    if(D) \
      *(D) = NULL; \
    if(dump2) \
      AlcFree(dump2); \
    if(dump1) \
      AlcFree(dump1); \
    if(dump0) \
      AlcFree(dump0); \
  } \
  return(alcErrno); \
}

/************************************************************************
* Template:	ALC_TEMPLATE_F3D					*
* Purpose:	Free's 3 dimensional arrays of any type, actualy no	*
*		type information is used in freeing the array.		*
* Global refs:	-							*
* Parameters:	D:			Pointer for array to be free'd.	*
*		F:			String with name of function.	*
************************************************************************/
#define ALC_TEMPLATE_F3D(D,F) \
{ \
  AlcErrno	alcErrno = ALC_ER_NONE; \
  \
  if((D == NULL) || (*(D) == NULL) || (**(D) == NULL)) \
  { \
    alcErrno = ALC_ER_NULLPTR; \
  } \
  else \
  { \
    AlcFree(**(D)); \
    AlcFree(*(D)); \
    AlcFree(D); \
  } \
  return(alcErrno); \
}

#ifdef __cplusplus
}  					       /* Close scope of 'extern "C" */
#endif
#endif /* ALCTEMPLATES_H */
