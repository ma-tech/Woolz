#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:      Mouse Atlas
* Title:        AlgBits.c
* Date:         May 2000
* Author:       Bill Hill
* Copyright:	2000 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:	Provides bit fiddling functions for the MRC Human
*		Genetics Unit numerical algorithm library.
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.
* 31-01-01 bill	Add AlgBitNextPowerOfTwo().
* 08-08-00 bill Add AlgBitNextSet().
************************************************************************/
#include <Alg.h>
#include <math.h>
#include <float.h>

/************************************************************************
* Function:	AlgBitSetCount
* Returns:	int:			Number of bits set.
* Purpose:	Counts the number of bits set in the given mask.
* Global refs:	-
* Parameters:   unsigned long gMsk:	Given bit mask.
************************************************************************/
int		AlgBitSetCount(unsigned long gMsk)
{
  int		cnt = 0;

  while(gMsk)
  {
    cnt += gMsk & 0x01;
    gMsk >>= 1;
  }
  return(cnt);
}

/************************************************************************
* Function:	AlgBitSetPositions
* Returns:	int:			Number of bits set.
* Purpose:	Counts the number of bits set in the given mask
*		and sets the first elements of the given bit position
*		array.
* Global refs:	-
* Parameters:   int *posA:		Bit set position array, MUST
*					have a length >= the number of
*					bits in the unsigned long bit
*					mask.
* 		unsigned long gMsk:	Given bit mask.
************************************************************************/
int		AlgBitSetPositions(int *posA, unsigned long gMsk)
{
  int		bIdx = 0,
  		cnt = 0;

  while(gMsk)
  {
    *(posA + cnt) = bIdx++;
    cnt += gMsk & 0x01;
    gMsk >>= 1;
  }
  return(cnt);
}

/************************************************************************
* Function:	AlgBitNextNOfM
* Returns:	unsigned long:		Next bit mask for n bits, or
*					zero on error.
* Purpose:	Computes a bit mask which has n bits set and is <=
*		maxMsk, the new bit mask is the next after the given
*		current mask.
*		The maximum number of bits that can be stored in an
*		unsigned long is assumed to be 32. If the actual number
*		is less than this then this code won't work, if it's
*		higher everything should be fine.
*		This function can be used together with
*		AlgBitSetPositions() to select all ordered combinations
*		of N of M. Eg to select all ordered combinations of 3
*		out of 5 this function would return 00111, 01011,
*		01101, 01110, 10110, 11010, 11100.
* Global refs:	-
* Parameters:	unsigned long curMsk:	Current mask value. If the
*					current value is zero then the
*					first valid mask will be
*					computed, otherwise the current
*					mask is assumed to be valid.
*		int n:			Number of bits to be set.
*		int m:			Maximum number of bits to use.
************************************************************************/
unsigned long	AlgBitNextNOfM(unsigned long curMsk, int n, int m)
{
  int		idx,
		search,
  		ok = 0;
  unsigned long	tmpMsk,
  		maxMsk,
  		nxtMsk = 0;

  if((n > 0) && (m > n) && (m < 32))
  {
    ok = 1;
    if(curMsk == 0)
    {
      /* Compute initial value for mask. */
      nxtMsk = (1 << n) - 1;
    }
    else
    {
      maxMsk = ((1 << n) - 1) << (m - n);
      if(curMsk >= maxMsk)
      {
        nxtMsk = maxMsk;
      }
      else
      {
        /* Find most significant bit pattern of '01'. */
	idx = m - 2;
	search = 1;
	while(search && (idx >= 0))
	{
	  tmpMsk = curMsk >> idx;
	  switch(tmpMsk & 3)
	  {
	    case 1:
	      search = 0;			       /* Found '01' pattern */
	      break;
	    case 0: /* FALLTHROUGH */
	    case 2:
	      idx -= 1;
	      break;
	    case 3:
	      idx -= 2;
	      break;
	  }
	}
	/* Replace the '01' with '10'. */
	if(search)
	{
	  nxtMsk = maxMsk;
	}
	else
	{
	  nxtMsk = (curMsk | (2 << idx)) & (~(1 << idx));
	}
      }
    }
  }
  return(nxtMsk);
}

/************************************************************************
* Function:	AlgBitNextSet
* Returns:	int:			Index of the next bit set in
*					the given mask, will be -ve if
*					no next bit set.
* Purpose:	Computes the index of the next bit set in the given mask
*		using the index to the current bit as the start index.
*		maxMsk, the new bit mask is the next after the given
*		current mask.
* Global refs:	-
* Parameters:	unsigned long msk:	Given bit mask.
*		int idC:		Current bit index, -ve to find
*					first bit set.
************************************************************************/
int		AlgBitNextSet(unsigned long msk, int idC)
{
  int		idN = -1,
  		limit;

  limit = sizeof(unsigned long) * 8;
  if(++idC < 0)
  {
    idC = 0;
  }
  while((idC < limit) && ((msk & (1 << idC)) == 0))
  {
    ++idC;
  }
  if(idC < limit)
  {
    idN = idC;
  }
  return(idN);
}

/************************************************************************
* Function:	AlgBitNextPowerOfTwo
* Returns:	int:			Index of the single bit that
*					should be set for an unsigned
*					integer that's >= the given int.
* Purpose:	Computes the next integer that is >= the given integer
*		and has only a single bit set.
* Global refs:	-
* Parameters:	unsigned int *dstP2I:		Destination ptr for integer
*					that's >= the given integer
*					and is a power of two.
*		unsigned int gI:	Given integer.
************************************************************************/
int		AlgBitNextPowerOfTwo(unsigned int *dstP2I, unsigned int gI)
{
  unsigned int	pI = 0,
  		nI = 1;

  while((nI < gI) && (pI < 31))
  {
    nI <<= 1;
    ++pI;
  }
  if(dstP2I)
  {
    *dstP2I = nI;
  }
  return(pI);
}
