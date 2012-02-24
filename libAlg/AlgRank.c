#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgRank_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgRank.c
* \author       Bill Hill
* \date         March 2002
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
* \brief	Rank selection algorithms which provide fast rank
*		selection from an array of values. This is the
*		general case of mimimum, median and maximum value
*		rank selection.
* \ingroup	AlgRank
*/

#include <stdio.h>
#include <Alg.h>

#ifndef _WIN32
   #include <unistd.h>
#endif

static void			AlgRankElmSwap(
				  void *elm0,
				  void *elm1,
				  unsigned int elmSz);
static void			AlgRankElmCopy(void *elm0,
				  void *elm1,
				  unsigned int elmSz);
/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of integers 
*		such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of elements.
* \param	nElm			Number of elements.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
*/
void		AlgRankSelectI(int *elm, int nElm, int rank)
{
  int		id0,
  		id1,
		id2,
		id3,
		tst,
		buf;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = nElm - 1;
  while(id2 < id3)
  {
    tst = *(elm + rank);
    id0 = id2;
    id1 = id3;
    do
    {
      while(*(elm + id0) < tst)
      {
	++id0;
      }
      while(tst < *(elm + id1))
      {
	--id1;
      }
      if(id0 <= id1)
      {
	buf = *(elm + id0);
	*(elm + id0++) = *(elm + id1);
	*(elm + id1--) = buf;
      }
    } while(id0 <= id1);
    if(id1 < rank)
    {
      id2 = id0;
    }
    if(rank < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of unsigned 
*		bytes such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of elements.
* \param	nElm			Number of elements.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
*/
void		AlgRankSelectUB(unsigned char *elm, int nElm, int rank)
{
  int		id0,
  		id1,
		id2,
		id3;
  unsigned char	tst,
		buf;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = nElm - 1;
  while(id2 < id3)
  {
    tst = *(elm + rank);
    id0 = id2;
    id1 = id3;
    do
    {
      while(*(elm + id0) < tst)
      {
	++id0;
      }
      while(tst < *(elm + id1))
      {
	--id1;
      }
      if(id0 <= id1)
      {
	buf = *(elm + id0);
	*(elm + id0++) = *(elm + id1);
	*(elm + id1--) = buf;
      }
    } while(id0 <= id1);
    if(id1 < rank)
    {
      id2 = id0;
    }
    if(rank < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of shorts 
*		such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of elements.
* \param	nElm			Number of elements.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
*/
void		AlgRankSelectS(short *elm, int nElm, int rank)
{
  int		id0,
  		id1,
		id2,
		id3;
  short		tst,
		buf;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = nElm - 1;
  while(id2 < id3)
  {
    tst = *(elm + rank);
    id0 = id2;
    id1 = id3;
    do
    {
      while(*(elm + id0) < tst)
      {
	++id0;
      }
      while(tst < *(elm + id1))
      {
	--id1;
      }
      if(id0 <= id1)
      {
	buf = *(elm + id0);
	*(elm + id0++) = *(elm + id1);
	*(elm + id1--) = buf;
      }
    } while(id0 <= id1);
    if(id1 < rank)
    {
      id2 = id0;
    }
    if(rank < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of floats 
*		such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of elements.
* \param	nElm			Number of elements.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
*/
void		AlgRankSelectF(float *elm, int nElm, int rank)
{
  int		id0,
  		id1,
		id2,
		id3;
  float		tst,
		buf;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = nElm - 1;
  while(id2 < id3)
  {
    tst = *(elm + rank);
    id0 = id2;
    id1 = id3;
    do
    {
      while(*(elm + id0) < tst)
      {
	++id0;
      }
      while(tst < *(elm + id1))
      {
	--id1;
      }
      if(id0 <= id1)
      {
	buf = *(elm + id0);
	*(elm + id0++) = *(elm + id1);
	*(elm + id1--) = buf;
      }
    } while(id0 <= id1);
    if(id1 < rank)
    {
      id2 = id0;
    }
    if(rank < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of doubles 
*		such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of elements.
* \param	nElm			Number of elements.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
*/
void		AlgRankSelectD(double *elm, int nElm, int rank)
{
  int		id0,
  		id1,
		id2,
		id3;
  double	tst,
		buf;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = nElm - 1;
  while(id2 < id3)
  {
    tst = *(elm + rank);
    id0 = id2;
    id1 = id3;
    do
    {
      while(*(elm + id0) < tst)
      {
	++id0;
      }
      while(tst < *(elm + id1))
      {
	--id1;
      }
      if(id0 <= id1)
      {
	buf = *(elm + id0);
	*(elm + id0++) = *(elm + id1);
	*(elm + id1--) = buf;
      }
    } while(id0 <= id1);
    if(id1 < rank)
    {
      id2 = id0;
    }
    if(rank < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Performs the minimum of sorting on an array of values
*		such that the n'th value in the array has rank n. That
*		is all values before the n'th are less than or equal to
*		it and all values after the n'th are greater than or equal
*		to it.
* \param	elm			Array of values.
* \param	nElm			Number of values.
* \param	elmSz			Size of each value.
* \param	rank			Required rank in the range
*					[0 - (nElm - 1)].
* \param	buf			A pointer to a single value to
*					be used as a workspace.
* \param	compFn			A value comparison function in
*					the same form as qsort(3).
*/
void		AlgRankSelectV(void *elm, int nElm, unsigned int elmSz,
			       int rank, void *buf,
			       int (*compFn)(void *, void *))
{
  int		id0,
  		id1,
		id2,
		id3,
		idR;
  char 		*elmP,
  		*tst;

  if(rank < 0)
  {
   rank = 0;
  }
  else if(rank >= nElm)
  {
    rank = nElm - 1;
  }
  id2 = 0;
  id3 = elmSz * (nElm - 1);
  idR = elmSz * rank;
  tst = (char *)buf;
  elmP = (char *)elm;
  while(id2 < id3)
  {
    AlgRankElmCopy(tst, elmP + idR, elmSz);
    id0 = id2;
    id1 = id3;
    do
    {
      while(compFn(elmP + id0, tst) < 0)
      {
	id0 += elmSz;
      }
      while(compFn(tst, elmP + id1) < 0)
      {
	id1 -= elmSz;
      }
      if(id0 <= id1)
      {
	AlgRankElmSwap(elmP + id0, elmP + id1, elmSz);
	id0 += elmSz;
	id1 -= elmSz;
      }
    } while(id0 <= id1);
    if(id1 < idR)
    {
      id2 = id0;
    }
    if(idR < id0)
    {
      id3 = id1;
    }
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Swaps to elements given their pointers and size.
* \param	elm0			Pointer to the first element.
* \param	elm1			Pointer to the second element.
* \param	elmSz			Element size.
*/
static void	AlgRankElmSwap(void *elm0, void *elm1, unsigned int elmSz)
{
  char		buf;
  char		*p0,
  		*p1;

  p0 = (char *)elm0;
  p1 = (char *)elm1;
  while(elmSz-- > 0)
  {
    buf = *p0;
    *p0++ = *p1;
    *p1++ = buf;
  }
}

/*!
* \return	void
* \ingroup	AlgRank
* \brief	Copies a single value.
* \param	elm0			Pointer to the element to be set.
* \param	elm1			Pointer to the element with value.
* \param	elmSz			Element size.
*/
static void	AlgRankElmCopy(void *elm0, void *elm1, unsigned int elmSz)
{
  char		*p0,
  		*p1;

  p0 = (char *)elm0;
  p1 = (char *)elm1;
  while(elmSz-- > 0)
  {
    *p0++ = *p1++;
  }
}

