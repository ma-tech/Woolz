#pragma ident "MRC HGU $Id$"
/*!
* \file         AlcVector.c
* \author       Bill Hill
* \date         March 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        A general purpose 1D vector (extensible array).
* \ingroup	AlcVector
* \todo		-
* \bug          None found.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alc.h>

/*!
* \return	New vector, or NULL on error.
* \ingroup	AlcVector
* \brief	Creates a new 1D vector (extensible array) with
*		the required element size and initial number of
*		elements. Vector elements are initialised by setting
*		all bytes to zero.
* \param	elmCnt			Initial number of elements.
* \param	elmSz			Size of elements.
* \param	blkSz			Number of elements per block,
*					also minimum initial number of
*					blocks allocated.
*					If blkSz <= 0 then a default
*					a defualt value is used.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
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
    nVec->freeStack = AlcFreeStackPush(NULL, data, &errNum);
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

/*!
* \return	Error code.
* \ingroup	AlcVector
* \brief	Free's the given vector.
* \param	vec			Vector to free.
*/
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

/*!
* \return	Error code.
* \ingroup	AlcVector
* \brief	Extend the vector for at least the given number of 
*		elements.
* \param	vec			Vector to extend.
* \param	elmCnt			Required number of elements.
*/
AlcErrno	AlcVectorExtend(AlcVector *vec, unsigned int elmCnt)
{
  int		blkIdx,
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
	  ++blkIdx;
	}
	while(blkIdx < nBlkCnt)
	{
	  *(nBlocks + blkIdx) = NULL;
	  ++blkIdx;
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

/*!
* \return	Vector item, or NULL on error.
* \ingroup	AlcVector
* \brief	Gets a pointer to the vector item with the given index.
* \param	vec			Vector to extend.
* \param	idx			Given item index.
*/
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


/*!
* \return	Vector item, or NULL on error.
* \ingroup	AlcVector
* \brief	Extends the vector and gets the vector item with the 
* 		given index
* \param	vec			Vector to extend.
* \param	idx			Given item index.
*/
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

/*!
* \return	Number of elements in vector.
* \ingroup	AlcVector
* \brief	Gets the number of elements that the vector can hold
*		before it needs to be extended.
* \param	vec			Vector to extend.
*/
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

/*!
* \return	<void>
* \ingroup	AlcVector
* \brief	Copies elements from the vector into a 1 dimensional
*		array.
* \param	vec		 	Given vector.
* \param	fIdx	 		Index of the first element in the
*					vector to copy.
* \param	lIdx	 		Index of the last element in the
*					vector to copy.
* \param	aM			The 1 dimensional array.
*/
void		AlcVectorSetArray1D(AlcVector *vec, int fIdx, int lIdx,
				    void *aM)
{
  int		aIdx,
		fVIdx,
		lVIdx,
  		fVBlkIdx,
		lVBlkIdx;

  aIdx = 0;
  fVBlkIdx = fIdx / vec->blkSz;
  lVBlkIdx = lIdx / vec->blkSz;
  fVIdx = fIdx % vec->blkSz;
  while(fVBlkIdx <= lVBlkIdx)
  {
    lVIdx = (fVBlkIdx == lVBlkIdx)? lIdx % vec->blkSz: vec->blkSz - 1;
    (void )memcpy((void *)((char *)aM + (aIdx * vec->elmSz)),
		  (void *)((char *)(*(vec->blocks + fVBlkIdx)) + 
			   (fVIdx * vec->elmSz)),
		  (lVIdx - fVIdx + 1) * vec->elmSz);
    aIdx += lVIdx - fVIdx + 1;
    fVIdx = 0;
    ++fVBlkIdx;
  }
}

/*!
* \return	Array of copied elements.
* \ingroup	AlcVector
* \brief	Creates a 1 dimensional array which contains a copy
*		of the vectors elements.
* \param	vec		 	Given vector.
* \param	fIdx	 		Index of the first element in the
*					vector to copy, becomes the first
*					element of the array.
* \param	lIdx	 		Index of the last element in the
*					vector to copy, becomes the last
*					element of the array.
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
void		*AlcVectorToArray1D(AlcVector *vec, int fIdx, int lIdx,
				    AlcErrno *dstErr)
{
  int		aIdx,
		fVIdx,
		lVIdx,
  		fVBlkIdx,
		lVBlkIdx;
  void		*aM = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(vec == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if((fIdx < 0) || (lIdx < fIdx))
  {
    errNum = ALC_ER_PARAM;
  }
  else if((aM = AlcCalloc(lIdx - fIdx + 1, vec->elmSz)) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    AlcVectorSetArray1D(vec, fIdx, lIdx, aM);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(aM);
}

/*!
* \return	Array of copied elements.
* \ingroup	AlcVector
* \brief	Creates a 2 dimensional array which contains a copy
*		of the vectors elements.
* \param	vec		 	Given vector.
* \param	fIdx	 		Index of the first element in the
*					vector to copy, becomes the first
*					element of the array.
* \param	lIdx	 		Index of the last element in the
*					vector to copy, becomes the last
*					element of the array.
* \param	nR			Number of rows (1D arrays).
* \param	nC			Number of columns (elements in each
*					1D array).
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
void		**AlcVectorToArray2D(AlcVector *vec, int fIdx, int lIdx,
				     int nR, int nC,
				     AlcErrno *dstErr)
{
  int		iR,
  		nRC;
  void		**aM = NULL;
  
  AlcErrno	errNum = ALC_ER_NONE;

  if(vec == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if((fIdx < 0) || (lIdx < fIdx) || (nR < 1) || (nC < 1) ||
          ((nRC = nR * nC) != lIdx - fIdx + 1))
  {
    errNum = ALC_ER_PARAM;
  }
  else if(((aM = (void **)AlcCalloc(nR, sizeof(void *))) ==  NULL) ||
          ((*aM = (void *)AlcMalloc(sizeof(char) * nRC * vec->elmSz)) == NULL))
  {
    errNum = ALC_ER_PARAM;
  }
  else
  {
    for(iR = 0; iR < nR; ++iR)
    {
      *(aM + iR) = (void *)(*(char **)aM + (vec->elmSz * nC * iR));
      AlcVectorSetArray1D(vec, fIdx + (iR * nC), fIdx + ((iR + 1) * nC) - 1,
      			  *(aM + iR));
    }
  }
  if(errNum != ALC_ER_NONE)
  {
    if(aM)
    {
      AlcFree(*aM);
      AlcFree(aM);
      aM = NULL;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(aM);
}
