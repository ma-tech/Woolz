#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgRange.c
* \author       Bill Hill
* \date         May 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides functions for computing the range of values
*		within a given array.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgRange
* @{
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

/*!
* \return				Error code.
* \brief	Computes the range of the given data, ie it's minimum
*		and maximum values.
* \param	datASz			Number of elements in given data
*					array.
* \param	datA			Data array to examine for minimum
*					and maximum values.
* \param	dstMin			Destination ptr for minimum
*					value, may be NULL.
* \param	dstMax			Destination ptr for maximum
*					value, may be NULL.
*/
AlgError	AlgRange1D(int datASz, double *datA,
			   double *dstMin, double *dstMax)
{
  int		cnt;
  double	datV,
  		minV,
  		maxV;
  double	*datP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((datASz < 1) || (datA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = datASz;
    datP = datA;
    maxV = minV = *datP;
    while(--cnt > 0)
    {
      if((datV = *++datP) > maxV)
      {
        maxV = datV;
      }
      else if (datV < minV)
      {
        minV = datV;
      }
    }
    if(dstMin)
    {
      *dstMin = minV;
    }
    if(dstMax)
    {
      *dstMax = maxV;
    }
  }
  return(algErr);
}

/*!
* \return				Error code.
* \brief	Computes the range of the given indexed data, ie it's
*		minimum and maximum values.
* \param	datA			Data array to examine for minimum
*					and maximum values.
* \param	idxASz			Number of elements in given
*					index array.
* \param	idxA			Index array with indicies into
*					the data buffer for the values
*					to examine.
* \param	dstMin			Destination ptr for minimum
*					value, may be NULL.
* \param	dstMax			Destination ptr for maximum
*					value, may be NULL.
*/
AlgError	AlgRangeIdx1D(double *datA, int idxASz, int *idxA,
			      double *dstMin, double *dstMax)
{
  int		cnt;
  double	datV,
  		minV,
  		maxV;
  int		*idxP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((idxASz < 1) || (datA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = idxASz;
    idxP = idxA;
    maxV = minV = *(datA + *idxP);
    while(--cnt > 0)
    {
      if((datV = *(datA + *++idxP)) > maxV)
      {
        maxV = datV;
      }
      else if (datV < minV)
      {
        minV = datV;
      }
    }
    if(dstMin)
    {
      *dstMin = minV;
    }
    if(dstMax)
    {
      *dstMax = maxV;
    }
  }
  return(algErr);
}

#ifdef ALG_RANGE_TEST
#define ALG_RANGE_TEST_BUFSZ (256)

int		main(int argc, char **argv)
{
  int		tI0,
  		ok = 0,
		done = 0,
  		datCnt = 0,
  		idxCnt = 0;
  double	min,
		max;
  FILE		*fP = NULL;
  AlgError	algErr = ALG_ERR_NONE;
  static char	ioBuf[ALG_RANGE_TEST_BUFSZ];
  static int	idxA[ALG_RANGE_TEST_BUFSZ];
  static double	datA[ALG_RANGE_TEST_BUFSZ];

  switch(argc)
  {
    case 1:
      ok = 1;
      fP = stdin;
      break;
    case 2:
      if((fP = fopen(*(argv + 1), "r")) == NULL)
      {
        (void )fprintf(stderr, "%s: Failed to open input file %s.\n",
		       *argv, *(argv + 1));
      }
      else
      {
        ok = 1;
      }
      break;
    default:
      break;
  }
  if(ok == 0)
  {
    (void )fprintf(stderr, "Usage: %s [<file>]\n%s\n",
    *argv,
    "A test for AlgRange1D() and AlgRangeIdx1D() which compute the minimum\n"
    "and maximum values in the given data. Input data records should have\n"
    "the form:\n"
    "  <x> [<flag>]\n"
    "where the second field is optional and indicates an indexed record\n"
    "if non-zero. If any indexed records are input only records which\n"
    "are indexed will be used, this can be used to test the indexed\n"
    "function AlgRangeIdx1D().\n");
  }
  else
  {
    do
    {
      if(fgets(ioBuf, ALG_RANGE_TEST_BUFSZ, fP) == NULL)
      {
	done = 1;
      }
      else
      {
	ioBuf[ALG_RANGE_TEST_BUFSZ - 1] = '\0';
	if(sscanf(ioBuf, "%lg %d",
		  datA + datCnt, &tI0) == 2)
	{
	  if(tI0)
	  {
	    *(idxA + idxCnt) = datCnt;
	    ++idxCnt;
	  }
	  ++datCnt;
	}
	else if(sscanf(ioBuf, "%lg",
		       datA + datCnt) == 1)
	{
	  ++datCnt;
	}
	else
	{
	  ok = 0;
	}
      }
    } while(ok && (datCnt < ALG_RANGE_TEST_BUFSZ) && !done);
    if(!ok)
    {
      (void )fprintf(stderr, "%s: Failed to read input data!\n", *argv);
    }
    else if(datCnt >= ALG_RANGE_TEST_BUFSZ)
    {
      (void )fprintf(stderr,
      		     "%s: This test is restricted to %d data records.\n",
		     ALG_RANGE_TEST_BUFSZ);
    }
  }
  if(fP && (argc == 2))
  {
    (void )fclose(fP);
  }
  if(ok)
  {
    if(idxCnt > 0)
    {
      algErr = AlgRangeIdx1D(datA, idxCnt, idxA, &min, &max);
    }
    else
    {
      algErr = AlgRange1D(datCnt, datA, &min, &max);
    }
    if(algErr == ALG_ERR_NONE)
    {
      (void )printf("%8g %8g\n", min, max);
    }
  }
  return(!ok);
}

#endif /* ALG_RANGE_TEST */

/*!
* @}
*/
