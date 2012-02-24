#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgShuffle_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgShuffle.c
* \author       Bill Hill
* \date         November 2000
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
* \brief        Functions for randomly permuting data.
* \ingroup	AlgRand
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Alg.h>

#if defined (CYGWIN) || defined (DARWIN) || defined (_WIN32)
#define drand48() (((double) rand()) / RAND_MAX)
#define srand48(X) (srand((unsigned int) X))
#define lrand48() ((long) ((((double) rand()) / RAND_MAX) * (1<<31)))
#endif /* CYGWIN || DARWIN */

/*!
* \return	Error code.
* \ingroup	AlgRand
* \brief	Inserts indicies into the given array which can be used
*		to shuffle data.
*		The permuted indicies are in the range [0-(nShuffle-1)]
*		with a random order.
*		This shuffling algorithm is based on a published
*		paper! Dursenfeld R. "Random Permutation" CACM July
*		1964 No 7 p420.
* \param	nShuffle		Number of permuted indicies
*					to insert into in shuffle.
* \param	shuffle			Array into which shuffled
*					indicies are placed.
* \param	seed			Seed fro random permutation.
*/
AlgError	AlgShuffleIdx(int nShuffle, int *shuffle, int seed)
{
  int		tI0,
		idx,
		rIdx;
  AlgError	algErr = ALG_ERR_NONE;

  if(shuffle == NULL)
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    srand48(seed);
    for(idx = 0; idx < nShuffle; ++idx)
    {
      *(shuffle + idx) = idx;
    }
    for(idx = 0; idx < nShuffle; ++idx)
    {
      rIdx = lrand48() % nShuffle;
      tI0 = *(shuffle + rIdx);
      *(shuffle + rIdx) = *(shuffle + idx);
      *(shuffle + idx) = tI0;
    }
  }
  return(algErr);
}

