#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlcBlockStack_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlc/AlcBlockStack.c
* \author       Bill Hill
* \date         November 1999
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
* \brief        A general purpose memory block allocator. Blocks are
*		allocated and stored in on a stack.
* \ingroup	AlcBlockStack
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/*!
* \return	New block stack, or NULL on error.
* \ingroup	AlcBlockStack
* \brief	Creates a new memory block with the required number of
*               elements and adds it to the stack of existing blocks.
* \param	nElem 		 	Number of elements in block.
* \param	elmSz 			Size of elements.
* \param	tBlk 			Top block, may be NULL.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
AlcBlockStack	*AlcBlockStackNew(size_t nElem, size_t elmSz,
				  AlcBlockStack *tBlk, AlcErrno *dstErr)
{
  AlcBlockStack	*nBlk = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(((nBlk = (AlcBlockStack *)AlcCalloc(1, sizeof(AlcBlockStack))) == NULL) ||
     ((nBlk->elements = AlcCalloc(nElem, elmSz)) == NULL))
  {
    if(nBlk)
    {
      AlcFree(nBlk);
    }
    nBlk = NULL;
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    nBlk->elmCnt = 0;
    nBlk->maxElm = nElem;
    if((nBlk->next = tBlk) != NULL)
    {
      tBlk->prev = nBlk;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nBlk);
}

/*!
* \return	Error code.
* \ingroup	AlcBlockStack
* \brief	Free's the given block and it's elements together with
*               all those below it in the stack. It's not an error for
*               a block's elements pointer to be NULL.
* \param	blk 		 	Given top block.
*/
AlcErrno	AlcBlockStackFree(AlcBlockStack *blk)
{
  AlcBlockStack	*nBlk;
  AlcErrno	errNum = ALC_ER_NONE;

  if(blk == NULL)
  {
    errNum =  ALC_ER_NULLPTR;
  }
  else
  {
    while(blk)
    {
      nBlk = blk->next;
      if(blk->elements)
      {
        AlcFree(blk->elements);
      }
      AlcFree(blk);
      blk = nBlk;
    }
  }
  return(errNum);
}
