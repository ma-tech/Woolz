#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgShuffle.c
* \author       Bill Hill
* \date         November 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions for randomly permuting data.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgShuffle
* @{
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alg.h>

#if defined (CYGWIN) || defined (DARWIN)
#define drand48() (((double) rand()) / RAND_MAX)
#define srand48(X) (srand((unsigned int) X))
#define lrand48() ((long) ((((double) rand()) / RAND_MAX) * (1<<31)))
#endif /* CYGWIN || DARWIN */

/*!
* \return				Error code.
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

#ifdef ALG_SHUFFLE_TEST

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idx,
		nData = 10,
		seed = 0,
		option,
		ok = 1,
		usage = 0;
  int		*data = NULL;
  static char	optList[] = "n:s:";

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
	if((sscanf(optarg, "%d", &nData) != 1) || (nData < 1))
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 's':
	if(sscanf(optarg, "%d", &seed) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((data = (int *)AlcMalloc(nData * sizeof(int))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to allocate storage.\n",
		     *argv);
    }
  }
  if(ok)
  {
    (void )AlgShuffleIdx(nData, data, seed);
    for(idx = 0; idx < nData; ++idx)
    {
      (void )printf("%d\n", *(data + idx));
    }
  }
  if(data)
  {
    AlcFree(data);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    		   "Usage: %s [-n#] [-s#]\n",
		   "Options:\n"
		   "  -n  Number of indicies.\n"
		   "  -s  Permutation seed.\n"
		   "Test for AlgShuffleIdx().\n");
  }
  return(!ok);
}
#endif /* ALG_SHUFFLE_TEST */

/*!
* @}
*/
