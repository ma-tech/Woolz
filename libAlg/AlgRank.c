#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgRank.c
* \author       Bill Hill
* \date         March 2002
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Rank selection algorithms which provide fast rank
*		selection from an array of values. This is the
*		general case of mimimum, median and maximum value
*		rank selection.
* \ingroup	AlgRank
* \todo         -
* \bug          None known.
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
  char		buf;
  char		*p0,
  		*p1;

  p0 = (char *)elm0;
  p1 = (char *)elm1;
  while(elmSz-- > 0)
  {
    *p0++ = *p1++;
  }
}

/* #define ALG_RANK_TEST */
#ifdef ALG_RANK_TEST

int		AlgRankCmpI(void *v0, void *v1)
{
  int		i0,
  		i1;

  i0 = *(int *)v0;
  i1 = *(int *)v1;
  return(i0 - i1);
}


extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
    		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int  		id0,
		tI0,
		option,
		seed = 0,
  		nElm = 10,
  		rank = -1,
  		okFlg = 1,
		usageFlg = 0;
  int		*datI = NULL,
  		*datV = NULL;
  unsigned char	*datUB = NULL;
  short		*datS = NULL;
  float 	*datF = NULL;
  double 	*datD = NULL;
  static char	optList[] = "hn:r:s:";

  opterr = 0;
  while(okFlg && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
	if((sscanf(optarg, "%d", &nElm) != 1) || (nElm <= 0))
	{
	  okFlg = 0;
	  usageFlg = 1;
	}
	break;
      case 'r':
	if((sscanf(optarg, "%d", &rank) != 1) || (rank < 0))
	{
	  okFlg = 0;
	  usageFlg = 1;
	}
	break;
      case 's':
	if(sscanf(optarg, "%d", &seed) != 1)
	{
	  okFlg = 0;
	  usageFlg = 1;
	}
	break;
      case 'h':
	okFlg = 0;
	usageFlg = 1;
	break;
    }
  }
  if(okFlg)
  {
    if(rank < 0)
    {
      rank = nElm / 2;
    }
    if(((datI = (int *)AlcMalloc(nElm * sizeof(int))) == NULL) ||
       ((datS = (short *)AlcMalloc(nElm * sizeof(short))) == NULL) ||
       ((datUB = (unsigned char *)AlcMalloc(nElm *
       					    sizeof(unsigned char))) == NULL) ||
       ((datF = (float *)AlcMalloc(nElm * sizeof(float))) == NULL) ||
       ((datD = (double *)AlcMalloc(nElm * sizeof(double))) == NULL) ||
       ((datV = (int *)AlcMalloc(nElm * sizeof(int))) == NULL))
    {
      (void )fprintf(stderr, "%s: Failed to allocate memory.\n", argv[0]);
    }
  }
  if(okFlg)
  {
    srandom(seed);
    for(id0 = 0; id0 < nElm; ++id0)
    {
      tI0 = random() % 256;
      *(datI + id0) = tI0;
      *(datS + id0) = tI0;
      *(datUB + id0) = tI0;
      *(datF + id0) = tI0;
      *(datD + id0) = tI0;
      *(datV + id0) = tI0;
    }
    AlgRankSelectI(datI, nElm, rank);
    AlgRankSelectS(datS, nElm, rank);
    AlgRankSelectUB(datUB, nElm, rank);
    AlgRankSelectF(datF, nElm, rank);
    AlgRankSelectD(datD, nElm, rank);
    AlgRankSelectV(datV, nElm, sizeof(int), rank, &tI0, AlgRankCmpI);
    (void )printf("idx    I   S   UB  F   D   V\n");
    for(id0 = 0; id0 < nElm; ++id0)
    {
      (void )printf("%-6d %-3d %-3d %-3d %-3g %-3g %-3d\n",
                    id0, *(datI + id0), *(datS + id0), *(datUB + id0),
                    *(datF + id0), *(datD + id0), *(datV + id0));
    }
    (void )printf("\n       %-3d %-3d %-3d %-3g %-3g %-3d\n",
		  *(datI + rank), *(datS + rank), *(datUB + rank),
		  *(datF + rank), *(datD + rank), *(datV + rank));
  }
  AlcFree(datI);
  AlcFree(datS);
  AlcFree(datUB);
  AlcFree(datF);
  AlcFree(datD);
  AlcFree(datV);
  if(usageFlg)
  {
    fprintf(stderr, "Usage: %s [-h] [-n#] [-r#] [-s#]\n%s",
            argv[0],
	    "Test for rank select, selects the item of rank r from a\n"
	    "list of random numbers.\n"
	    "Options are:\n"
	    "  -h  Usage, prints this information.\n"
	    "  -n  The number of random elements to select from.\n"
	    "  -r  The rank to be selected, 0 smallest and n - 1 largest.\n"
	    "      The default is the median.\n"
	    "  -s  The random number seed.\n");
  }
  return(!okFlg);
}
#endif /* ALG_RANK_TEST */
