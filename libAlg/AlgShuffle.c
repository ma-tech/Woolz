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
* \ingroup	AlgRand
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Alg.h>

#if defined (CYGWIN) || defined (DARWIN) || defined (WIN32)
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

#if (ALG_SHUFFLE_TEST == 1)

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
#endif /* ALG_SHUFFLE_TEST == 1 */

#if (ALG_SHUFFLE_TEST == 2)

#define MAX_BUF_LEN	(1000)
#define DATA_CHUNK_SZ	(10000)

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idx,
		maxData = 0,
		nData = 0,
		seed = 0,
		option,
		ok = 1,
		usage = 0;
  int		*iData = NULL;
  char		**sData = NULL;
  FILE		*fP;
  char		*inFileStr;
  char		rec[MAX_BUF_LEN];
  static char	optList[] = "hs:",
  		defFileStr[]= "-";

  opterr = 0;
  inFileStr = defFileStr;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 's':
        if(sscanf(optarg, "%d", &seed) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((fP = (strcmp(inFileStr, "-")?
    	     fopen(inFileStr, "r"): stdin)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: failed to open %s.\n", *argv, inFileStr);
      
    }
  }
  idx = 0;
  while(ok && (fgets(rec, MAX_BUF_LEN, fP) != NULL))
  {
    *(rec + MAX_BUF_LEN - 1) = '\0';
    *(rec + strlen(rec) - 1) = '\0';
    if(maxData <= idx)
    {
      maxData = (maxData)? maxData * 2: DATA_CHUNK_SZ;
      if((sData = (char **)AlcRealloc(sData,
				      maxData * sizeof(char *))) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to allocate memory.\n", *argv);
      }
    }
    if(ok)
    {
      if((*(sData + idx++) = AlcStrDup(rec)) == NULL)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to allocate memory.\n", *argv);
      }
    }
  }
  if(ok)
  {
    nData = idx;
    if((iData = (int *)AlcMalloc(nData * sizeof(int))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to allocate memory.\n", *argv);
    }
  }
  if(ok)
  {
    (void )AlgShuffleIdx(nData, iData, seed);
    for(idx = 0; idx < nData; ++idx)
    {
      (void )printf("%s\n", *(sData + *(iData + idx)));
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    		   "Usage: %s [-h] [-s#]\n",
		   "Randomize or shuffle the input data records.\n"
		   "Options:\n"
		   "  -h  Print this usage message.\n"
		   "  -s  Random number generator seed.\n",
		   *argv);
  }
  return(!ok);
}
#endif /* ALG_SHUFFLE_TEST == 2 */
