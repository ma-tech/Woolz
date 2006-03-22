#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcFreeStack_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlc/AlcFreeStack.c
* \author       Bill Hill
* \date         March 2000
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
* \brief        A general purpose free stack which allows a single
*               pointer to be used to keep a list of data to be free'd.
* \ingroup	AlcFreeStack
* \todo		-
* \bug          None known.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>


typedef struct _AlcFreeStack
{
  void		*data;
  struct _AlcFreeStack *prev;
} AlcFreeStack;

/*!
* \return	New free stack pointer or NULL on error.
* \ingroup	AlcFreeStack
* \brief	Push's the given pointer onto the free stack on top
*		of the previous free stack pointer.
* \param	prev 			Previous free stack pointer.
* \param	data 			New pointer to push onto the
*					free stack.
* \param	dstErr 			Destination error pointer,
*					may be NULL
*/
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

/*!
* \return	New free stack pointer or NULL on error.
* \ingroup	AlcFreeStack
* \brief	Pop's the top entry from the free stack. Returns a
*		free stack pointer and set's the given destination
*		pointer to the entry's data. The entry's data is NOT
*		free'd
* \param	prev 			The free stack.
* \param	dstData 		Destination data pointer, may
*					be NULL.
* \param	dstErr 			Destination error pointer,
*					may be NULL.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcFreeStack
* \brief	Free's all entries on the given free stack.
* \param	stack 			The stack of pointers to be free'd.
*/
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
