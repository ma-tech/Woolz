#pragma ident "MRC HGU $Id"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcBlockStack.c
* Date:         March 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A general purpose free stack which allows a single
*		pointer to be used to keep a list of data to be
*		free'd.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

typedef struct _AlcFreeStack
{
  void		*data;
  struct _AlcFreeStack *prev;
} AlcFreeStack;

/************************************************************************
* Function:	AlcFreeStackPush
* Returns:	void *:			New free stack pointer or
*					NULL on error.
* Purpose:	Push's the given pointer onto the free stack on top
*		of the previous free stack pointer.
* Global refs:	
* Parameters:	void *prev:		Previous free stack pointer.
*		void *data:		New pointer to push onto the
*					free stack.
*		AlcErrno *dstErr:	Destination error pointer,
*					may be NULL
************************************************************************/
void 		*AlcFreeStackPush(void *prev, void *data, AlcErrno *dstErr)
{
  AlcFreeStack *fPtr = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if((fPtr = (AlcFreeStack *)AlcMalloc(sizeof(AlcFreeStack))) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    fPtr->data = data;
    fPtr->prev = (AlcFreeStack *)prev;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return((void *)fPtr);
}

/************************************************************************
* Function:	AlcFreeStackPop
* Returns:	void *:			New free stack pointer or
*					NULL on error.
* Purpose:	Pop's the top entry from the free stack. Returns a
*		free stack pointer and set's the given destination
*		pointer to the entry's data. The entry's data is NOT
*		free'd
* Global refs:	-
* Parameters:	void *prev:		The free stack.
*		void **dstData:		Destination data pointer, may
*					be NULL.
*		AlcErrno *dstErr:	Destination error pointer,
*					may be NULL.
************************************************************************/
void 		*AlcFreeStackPop(void *prev, void **dstData, AlcErrno *dstErr)
{
  AlcFreeStack	*entry0,
  		*entry1 = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(prev == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    entry0 = (AlcFreeStack *)prev;
    entry1 = entry0->prev;
    if(dstData)
    {
      *dstData = entry0->data;
    }
    AlcFree(entry0);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return((void *)entry1);
}

/************************************************************************
* Function:	AlcFreeStackFree
* Returns:	AlcErrno:		Alc error code.
* Purpose:	Free's all entries on the given free stack.
* Global refs:	-
* Parameters:	void *stack:		The stack of pointers to be
*					free'd.
************************************************************************/
AlcErrno	AlcFreeStackFree(void *stack)
{
  AlcFreeStack	*entry0,
  		*entry1;
  AlcErrno	errNum = ALC_ER_NONE;

  if(stack == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else
  {
    entry0 = (AlcFreeStack *)stack;
    while(entry0)
    {
      entry1 = entry0;
      entry0 = entry1->prev;
      if(entry1->data)
      {
        AlcFree(entry1->data);
      }
      AlcFree(entry1);
    }
  }
  return(errNum);
}
