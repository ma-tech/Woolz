#pragma ident "MRC HGU $Id"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcBlockStack.c
* Date:         November 1999
* Author:       Bill Hill
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A general purpose memory block allocator. Blocks are
*		allocated and stored in on a stack.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/************************************************************************
* Function:	AlcBlockStackNew
* Returns:	AlcBlockStack *:	New block stack, or NULL on
*					error.
* Purpose:	Creates a new memory block with the required number of
*		elements and adds it to the stack of existing blocks.
* Global refs:	-
* Parameters:	unsigned nElem:		Number of elements in block.
*		unsigned elmSz:		Size of elements.
*		AlcBlockStack *tBlk:	Top block, may be NULL.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
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

/************************************************************************
* Function:	AlcBlockStackFree
* Returns:	AlcErrno:		Error code.
* Purpose:	Free's the given block and it's elements together with
* 		all those below it in the stack. It's not an error for
*		a block's elements pointer to be NULL.
* Global refs:	-
* Parameters:	AlcBlockStack *blk:	Given top block.
************************************************************************/
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
