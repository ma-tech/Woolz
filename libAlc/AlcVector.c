#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlcVector_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlc/AlcVector.c
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
AlcVector	*AlcVectorNew(size_t elmCnt, size_t elmSz,
			      size_t blkSz, AlcErrno *dstErr)
{
  size_t	blkIdx,
		blkInc,
		blkCnt = 0,
		blkUse = 0;
  char		*data = NULL;
  AlcVector	*nVec = NULL;
  AlcErrno	errNum = ALC_ER_NONE;
  const size_t defaultBlkSz = 1024;

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
      (void )AlcFreeStackFree(vec->freeStack);
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
AlcErrno	AlcVectorExtend(AlcVector *vec, size_t elmCnt)
{
  size_t	blkIdx,
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
  if(nBlkUse > vec->blkUse)
  {
    /* Check the number of block pointers and allocate more if required */
    eBlkUse = nBlkUse - vec->blkUse;
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
void		*AlcVectorItemGet(AlcVector *vec, size_t idx)
{
  size_t	blkIdx;
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
void		*AlcVectorExtendAndGet(AlcVector *vec, size_t idx)
{
  size_t	blkIdx;
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
size_t		AlcVectorCount(AlcVector *vec)
{
  size_t	cnt = 0;

  if(vec)
  {
    cnt = vec->blkUse * vec->blkSz;
  }
  return(cnt);
}

/*!
* \return	void
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
void		AlcVectorSetArray1D(AlcVector *vec, size_t fIdx, size_t lIdx,
				    void *aM)
{
  size_t	aIdx,
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
void		*AlcVectorToArray1D(AlcVector *vec, size_t fIdx, size_t lIdx,
				    AlcErrno *dstErr)
{
  void		*aM = NULL;
  AlcErrno	errNum = ALC_ER_NONE;

  if(vec == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if(lIdx < fIdx)
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
* \param	nR			Number of rows (1D arrays).
* \param	nC			Number of columns (elements in each
*					1D array).
* \param	dstErr			Destination pointer for error
*					code, may be NULL.
*/
void		**AlcVectorToArray2D(AlcVector *vec, size_t nR, size_t nC,
				     AlcErrno *dstErr)
{
  size_t	iR;
  void		**aM = NULL;
  
  AlcErrno	errNum = ALC_ER_NONE;

  if(vec == NULL)
  {
    errNum = ALC_ER_NULLPTR;
  }
  else if((nR < 1) || (nC < 1))
  {
    errNum = ALC_ER_PARAM;
  }
  else if(((aM = (void **)AlcCalloc(nR, sizeof(void *))) ==  NULL) ||
          ((*aM = (void *)
	          AlcMalloc(sizeof(char) * nR * nC * vec->elmSz)) == NULL))
  {
    errNum = ALC_ER_PARAM;
  }
  else
  {
    for(iR = 0; iR < nR; ++iR)
    {
      *(aM + iR) = (void *)(*(char **)aM + (vec->elmSz * nC * iR));
      AlcVectorSetArray1D(vec, iR * nC, ((iR + 1) * nC) - 1, *(aM + iR));
    }
  }
  if(errNum != ALC_ER_NONE)
  {
    if(aM != NULL)
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


/*!
* \return	Extensible vector containing values for a 1D array or
* 		NULL on error.
* \ingroup	AlcVector
* \brief	Reads a 1D double array from the given numeric ASCI file.
*		Each value should be on a seperate line.
* \param	fP:			File pointer.
* \param	fs			Field seperator chracter string. If
* 					NULL whitespace assumed.
* \param	recMax			Maximum record length.
* \param	dstNV			Destination pointer for the number of
*					vector values. Must not be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlcVector	*AlcVecReadDouble1Asci(FILE *fP, const char *fSep,
				       size_t recMax, size_t *dstNV,
				       AlcErrno *dstErr)
{
  size_t	nV = 0;
  char		*recS = NULL;
  const char	*iFS;
  AlcVector	*vec = NULL;
  AlcErrno	errNum = ALC_ER_NONE;
  const size_t	vecCnt = 1024;		/* Initial number of elements in
  					 * the vector */

  iFS = (fSep == NULL)? " \t\n\r": fSep;
  if((recS = AlcMalloc(sizeof(char) * recMax)) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    vec = AlcVectorNew(vecCnt, sizeof(double), vecCnt, &errNum);
  }
  while((errNum == ALC_ER_NONE) && (fgets(recS, recMax, fP) != NULL))
  {
    char	*tokS,
    		*parseS;

    parseS = recS;
    while((errNum == ALC_ER_NONE) &&
          ((tokS = (char *)strtok(parseS, iFS)) != NULL) && *tokS)
    {
      double	v;

      if(sscanf(tokS, "%lg", &v) != 1)
      {
	errNum = ALC_ER_READ;
      }
      else
      {
        double	*dP0;

        if((dP0 = (double *)AlcVectorExtendAndGet(vec, nV)) == NULL)
	{
	  errNum = ALC_ER_ALLOC;
	}
	else
	{
	  *dP0 = v;
          ++nV;
	}
      }
    }
  }
  AlcFree(recS);
  if(errNum == ALC_ER_NONE)
  {
    *dstNV = nV;
  }
  else
  {
    (void )AlcVectorFree(vec);
    vec = NULL;
  }
  return(vec);
}
/*!
* \return	Extensible vector containing values for a 2D array or
* 		NULL on error.
* \ingroup	AlcVector
* \brief	Reads a 2D double array from the given numeric ASCI file
* 		into an extensible vector.
*		The number of fields per record must be the same for all
*		records.
* \param	fP:			File pointer.
* \param	fs			Field seperator chracter string. If
* 					NULL whitespace assumed.
* \param	recMax			Maximum record length.
* \param	dstNR			Destination pointer for the number of
*					rows (records) read. Must not be NULL.
* \param	dstNC			Destination pointer for the number of
*					columns, ie the number of elements in
*					each 1D array (number of fields per
*					record). Must not be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
AlcVector	*AlcVecReadDouble2Asci(FILE *fP, const char *fSep,
				       size_t recMax, 
                                       size_t *dstNR, size_t *dstNC,
				       AlcErrno *dstErr)
{
  size_t	iF,
  		nR = 0,
  		nF = 0,
		nV = 0;
  char		*recS = NULL;
  const char	*iFS;
  AlcVector	*vec = NULL;
  AlcErrno	errNum = ALC_ER_NONE;
  const size_t	vecCnt = 1024; 		/* Initial number of elements in
  					 * the vector */

  iFS = (fSep == NULL)? " \t\n\r": fSep;
  if((recS = AlcMalloc(sizeof(char) * recMax)) == NULL)
  {
    errNum = ALC_ER_ALLOC;
  }
  else
  {
    vec = AlcVectorNew(vecCnt, sizeof(double), vecCnt, &errNum);
  }
  while((errNum == ALC_ER_NONE) && (fgets(recS, recMax, fP) != NULL))
  {
    char	*tokS,
    		*parseS;

    iF = 0;
    parseS = recS;
    while((errNum == ALC_ER_NONE) &&
          ((tokS = (char *)strtok(parseS, iFS)) != NULL) && *tokS)
    {
      double	*dP0;

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
  AlcFree(recS);
  if(errNum == ALC_ER_NONE)
  {
    *dstNC = nF;
    *dstNR = nR;
  }
  else
  {
    (void )AlcVectorFree(vec);
    vec = NULL;
  }
  return(vec);
}
