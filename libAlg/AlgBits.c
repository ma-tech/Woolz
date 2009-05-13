#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgBits_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libAlg/AlgBits.c
* \author       Bill Hill
* \date         May 2000
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
* \brief	Provides bit fiddling functions.
* \ingroup	AlgBits
* \todo         -
* \bug          None known.
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

/*!
* \return	Number of bits set.
* \ingroup	AlgBits
* \brief	Counts the number of bits set in the given mask.
* \param        gMsk			Given bit mask.
*/
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

/*!
* \return	Number of bits set.
* \ingroup	AlgBits
* \brief	Counts the number of bits set in the given mask
*		and sets the first elements of the given bit position
*		array.
* \param	posA 			Bit set position array, MUST
*					have a length >= the number of
*					bits in the unsigned long bit
*					mask.
* \param	gMsk		 	Given bit mask.
*/
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

/*!
* \return	Next bit mask for n bits, or zero on error.
* \ingroup	AlgBits
* \brief	Computes a bit mask which has n bits set and is <=
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
* \param	curMsk			Current mask value. If the
*					current value is zero then the
*					first valid mask will be
*					computed, otherwise the current
*					mask is assumed to be valid.
* \param	n			Number of bits to be set.
* \param	m			Maximum number of bits to use.
*/
unsigned long	AlgBitNextNOfM(unsigned long curMsk, int n, int m)
{
  int		idx,
		search;
  unsigned long	tmpMsk,
  		maxMsk,
  		nxtMsk = 0;

  if((n > 0) && (m > n) && (m < 32))
  {
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

/*!
* \return	Index of the next bit set in the given mask, will be -ve if
*		no next bit set.
* \ingroup	AlgBits
* \brief	Computes the index of the next bit set in the given mask
*		using the index to the current bit as the start index.
*		maxMsk, the new bit mask is the next after the given
*		current mask.
* \param	msk			Given bit mask.
* \param	idC			Current bit index, -ve to find
*					first bit set.
*/
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

/*!
* \return	Index of the single bit that should be set for an unsigned
*		integer that's >= the given int.
* \ingroup	AlgBits
* \brief	Computes the next integer that is >= the given integer
*		and has only a single bit set.
* \param	dstP2I			Destination ptr for integer
*					that's >= the given integer
*					and is a power of two.
* \param	gI			Given integer.
*/
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
