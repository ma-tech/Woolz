#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcBlockStack.c
* \author       Bill Hill
* \date         November 1999
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        A general purpose memory block allocator. Blocks are
*		allocated and stored in on a stack.
* \todo		-
* \bug          None found.
*/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/*!
* \ingroup      Alc
* \defgroup	AlcBlockStack
* @{
*/

/*!
* \return				New block stack, or NULL on
*					error.
* \brief	Creates a new memory block with the required number of
*               elements and adds it to the stack of existing blocks.
* \param	nElem 		 	Number of elements in block.
* \param	elmSz 			Size of elements.
* \param	tBlk 			Top block, may be NULL.
* \param	dstErr 			Destination pointer for error
*					code, may be NULL.
*/
AlcBlockStack	*AlcBlockStackNew(unsigned nElem, unsigned elmSz,
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
* \return		 		Error code.
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

/*!
* @}
*/
