#pragma ident "MRC HGU $Id"
/***********************************************************************
* Project:      Mouse Atlas
* Title:        AlcVector.c
* Date:         March 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      A general purpose 1D vector (extensible array).
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <Alc.h>

/************************************************************************
* Function:	AlcVectorNew
* Returns:	AlcVector *:		New vector, or NULL on error.
* Purpose:	Creates a new 1D vector (extensible array) with
*		the required element size and initial number of
*		elements. Vector elements are initialised by setting
*		all bytes to zero.
* Global refs:	-
* Parameters:	unsigned int elmCnt:	Initial number of elements.
*		unsigned int elmSz:	Size of elements.
*		unsigned int int blkSz:	Number of elements per block,
*					also minimum initial number of
*					blocks allocated.
*					If blkSz <= 0 then a default
*					a defualt value is used.
*		AlcErrno *dstErr:	Destination pointer for error
*					code, may be NULL.
************************************************************************/
AlcVector	*AlcVectorNew(unsigned int elmCnt, unsigned int elmSz,
			      unsigned int blkSz, AlcErrno *dstErr)
{
  unsigned int	blkCnt,
		blkUse,
		blkIdx,
		blkInc;
  char		*data;
  AlcVector	*nVec = NULL;
  AlcErrno	errNum = ALC_ER_NONE;
  const unsigned int defaultBlkSz = 1024;

  if(elmSz <= 0)
  {
    errNum = ALC_ER_PARAM;
  }
  else
  {
    if(blkSz <= 0)
    {
      blkSz = defaultBlkSz;
    }
    blkUse = (elmCnt + blkSz - 1) / blkSz;
    blkCnt = ((blkUse + blkSz - 1) / blkSz) * blkSz;
    elmCnt = blkUse * blkSz;
    if(((nVec = (AlcVector *)AlcCalloc(1, sizeof(AlcVector))) == NULL) ||
       ((nVec->blocks = (void **)AlcCalloc(blkCnt, sizeof(void *))) == NULL) ||
       ((data = (char *)AlcCalloc(elmCnt, elmSz)) == NULL))
    {
      errNum = ALC_ER_ALLOC;
    }
  }
  if(errNum == ALC_ER_NONE)
  {
    /* Push block base pointer onto the free stack. */
    nVec->freeStack = AlcFreeStackPush(NULL, *(nVec->blocks), &errNum);
  }
  if(errNum == ALC_ER_NONE)
  {
    nVec->elmSz = elmSz;
    nVec->blkUse = blkUse;
    nVec->blkCnt = blkCnt;
    nVec->blkSz = blkSz;
    /* Set block pointers */
    blkIdx = 0;
    blkInc = elmSz * blkSz;
    while(blkIdx < blkUse)
    {
      *(nVec->blocks + blkIdx) = (void *)data;
      data += blkInc;
      ++blkIdx;
    }
  }
  else
  {
    /* Tidy up if hit memory allocation error. */
    if(nVec)
    {
      if(nVec->blocks)
      {
	if(*(nVec->blocks))
	{
	  AlcFree(*(nVec->blocks));
	}
        AlcFree(nVec->blocks);
      }
      AlcFree(nVec);
    }
    nVec = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(nVec);
}

/************************************************************************
* Function:	AlcVectorFree
* Returns:	AlcErrno:		Error code.
* Purpose:	Free's the given vector.
* Global refs:	-
* Parameters:	AlcVector *vec:		Vector to free.
************************************************************************/
AlcErrno	AlcVectorFree(AlcVector *vec)
{
  AlcErrno	errNum = ALC_ER_NONE;

  if(vec == NULL)
  {
    errNum =  ALC_ER_NULLPTR;
  }
  else
  {
    if(vec->freeStack)
    {
      AlcFreeStackFree(vec->freeStack);
    }
    if(vec->blocks)
    {
      AlcFree(vec->blocks);
    }
    AlcFree(vec);
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcVectorExtend
* Returns:	AlcErrno:		Error code.
* Purpose:	Extend the vector for at least the given number of 
*		elements.
* Global refs:	-
* Parameters:	AlcVector *vec:		Vector to extend.
*		unsigned int elmCnt:	Required number of elements.
************************************************************************/
AlcErrno	AlcVectorExtend(AlcVector *vec, unsigned int elmCnt)
{
  unsigned int	blkIdx,
  		nBlkUse,
  		nBlkCnt,
  		eBlkUse,
		eDataCnt,
		blkInc;
  char		*eData;
  void		**nBlocks,
  		**oBlocks;
  AlcErrno	errNum = ALC_ER_NONE;

  nBlkUse = (elmCnt + vec->blkSz - 1) / vec->blkSz;
  if((eBlkUse = nBlkUse - vec->blkUse) > 0)
  {
    /* Check the number of block pointers and allocate more if required */
    if(nBlkUse > vec->blkCnt)
    {
      nBlkCnt = ((nBlkUse + vec->blkSz - 1) / vec->blkSz) *
	        vec->blkSz;
      if((nBlocks = (void *)AlcCalloc(nBlkCnt, sizeof(void *))) == NULL)
      {
        errNum = ALC_ER_ALLOC;
      }
      else
      {
        blkIdx = 0;
	while(blkIdx < vec->blkUse)
	{
	  *(nBlocks + blkIdx) = *(vec->blocks + blkIdx);
	}
	while(blkIdx < nBlkCnt)
	{
	  *(nBlocks + blkIdx) = NULL;
	}
	oBlocks = vec->blocks;
	vec->blocks = nBlocks;
        vec->blkCnt = nBlkCnt;
	AlcFree(oBlocks);
      }
    }
    if(errNum == ALC_ER_NONE)
    {
      /* Allocate the required number of blocks */
      eDataCnt = eBlkUse * vec->blkSz;
      if((eData = (char *)AlcCalloc(eDataCnt, vec->elmSz)) == NULL)
      {
        errNum = ALC_ER_ALLOC;
      }
      else 
      {
	/* Push block base pointer onto the free stack. */
	vec->freeStack = AlcFreeStackPush(vec->freeStack, eData, &errNum);
	if(errNum != ALC_ER_NONE)
	{
	  AlcFree(eData);
	}
      }
    }
    if(errNum == ALC_ER_NONE)
    {
      /* Append the new data blocks */
      blkIdx = vec->blkUse;
      blkInc = vec->elmSz * vec->blkSz;
      while(blkIdx < nBlkUse)
      {
        *(vec->blocks + blkIdx) = (void *)eData;
	++blkIdx;
	eData += blkInc;
      }
      vec->blkUse = nBlkUse;
    }
  }
  return(errNum);
}

/************************************************************************
* Function:	AlcVectorItemGet
* Returns:	void *:			Vector item, or NULL on error.
* Purpose:	Gets a pointer to the vector item with the given index.
* Global refs:	-
* Parameters:	AlcVector *vec:		Vector to extend.
*		unsigned int idx:	Given item index.
************************************************************************/
void		*AlcVectorItemGet(AlcVector *vec, unsigned int idx)
{
  unsigned int	blkIdx;
  void		*data = NULL;

  if(vec && ((blkIdx = idx / vec->blkSz) < vec->blkUse))
  {
    data = (char *)(*(vec->blocks + blkIdx)) +
    	   ((idx % vec->blkSz) * vec->elmSz);
  }
  return(data);
}


/************************************************************************
* Function:	AlcVectorExtendAndGet
* Returns:	void *:			Vector item, or NULL on error.
* Purpose:	Extends the vector and gets the vector item with the 
* 		given index
* Global refs:	-
* Parameters:	AlcVector *vec:		Vector to extend.
*		unsigned int idx:	Given item index.
************************************************************************/
void		*AlcVectorExtendAndGet(AlcVector *vec, unsigned int idx)
{
  unsigned int	blkIdx;
  void		*data = NULL;

  if((AlcVectorExtend(vec, idx + 1) == ALC_ER_NONE) &&
     ((blkIdx = idx / vec->blkSz) < vec->blkUse))
  {
    
    data = (char *)(*(vec->blocks + blkIdx)) +
    	   ((idx % vec->blkSz) * vec->elmSz);
  }
  return(data);
}

/************************************************************************
* Function:	AlcVectorCount
* Returns:	unsigned int:		Number of elements in vector.
* Purpose:	Gets the number of elements that the vector can hold
*		before it needs to be extended.
* Global refs:	-
* Parameters:	AlcVector *vec:		Vector to extend.
************************************************************************/
unsigned int	AlcVectorCount(AlcVector *vec)
{
  unsigned int	cnt = 0;
  void		*data = NULL;

  if(vec)
  {
    cnt = vec->blkUse * vec->blkSz;
  }
  return(cnt);
}
