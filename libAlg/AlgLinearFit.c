#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgLinearFit.c
* \author       Bill Hill
* \date         may 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Provides functions for fitting linear models to data,
*		ie linear regression.
* \todo         -
* \bug          None known.
*/

/*!
* \ingroup      Alg
* \defgroup     AlgLinearFit
* @{
*/

#include <Alg.h>
#include <math.h>
#include <float.h>

/*!
* \return				Error code.
* \brief	Computes the least squares best fit straight line
*		(y = a + bx) through the given data, ie linear
*		regression.
*		This function is based on the function fit():
*		Press W. H., Teukolsky S. A., Vetterling W. T.
*		and Flannery B. P, Numerical Recipies in C,
*		1992, CUP.
*  \param	datSz			Number of elements in given
*					data arrays array.
* \param	datXA			Data array with 'x' values.
* \param	datYA			Data array with 'y' values.
* \param	dstA			Destination ptr for intercept
*					'a', may be NULL.
* \param	dstB			Destination ptr for gradient
*					'b', may be NULL.
* \param	dstSigA			Destination ptr for std dev
*					of 'a', may be NULL.
* \param	dstSigB			Destination ptr for std dev
*					of 'b', may be NULL.
* \param	dstQ			Destination ptr for goodness
*					of fit, may be NULL.
*/
AlgError	AlgLinearFit1D(int datSz, double *datXA, double *datYA,
			       double *dstA, double *dstB,
			       double *dstSigA, double *dstSigB,
			       double *dstQ)
{
  int		cnt;
  double	tD0,
  		a = 0.0,
  		b = 0.0,
		t,
		sTSq = 0.0,
		chiSq = 0.0,
  		sX = 0.0,
		sY = 0.0,
		sigA,
		sigB,
		q;
  double	*datXP,
  		*datYP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((datSz < 3) || (datXA == NULL) || (datYA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      sX += *datXP++;
      sY += *datYP++;
    }
    tD0 = sX / datSz;
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      t = *datXP++ - tD0;
      sTSq += t * t;
      b += t * *datYP++;
    }
    b /= sTSq;
    a = (sY - (sX * b)) / datSz;
    sigA = sqrt((1.0 + ((sX * sX) / (datSz * sTSq))) / datSz);
    sigB = sqrt(1.0 / sTSq);
    cnt = datSz;
    datXP = datXA;
    datYP = datYA;
    while(cnt-- > 0)
    {
      tD0 = *datYP++ - a - (b * *datXP++);
      chiSq += tD0 * tD0;
    }
    tD0 = sqrt(chiSq / (datSz - 2));
    sigA *= tD0;
    sigB *= tD0;
    q = AlgGammaP(0.5 * (datSz - 2), 0.5 * (chiSq), &algErr);
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(dstA)
    {
      *dstA = a;
    }
    if(dstB)
    {
      *dstB = b;
    }
    if(dstSigA)
    {
      *dstSigA = sigA;
    }
    if(dstSigB)
    {
      *dstSigB = sigB;
    }
    if(dstQ)
    {
      *dstQ = q;
    }
  }
  return(algErr);
}

/*!
* \return				Error code.
* \brief	Computes the least squares best fit straight line
*		(y = a + bx) through the given data, ie linear
*		regression.
* \param	datXA			Data array with 'x' values.
* \param	datYA			Data array with 'y' values.
* \param	idxXA			Index array with indicies into
*					the 'x' data buffer.
*					for the values use examine.
* \param	idxYA			Index array with indicies into
*					the 'y' data buffer.
*					for the values use examine.
* \param	idxASz			Number of elements in each of
*					the given index arrays.
* \param	dstA			Destination ptr for intercept
*					'a', may be NULL.
* \param	dstB			Destination ptr for gradient
*					'b', may be NULL.
* \param	dstSigA			Destination ptr for std dev
*					of 'a', may be NULL.
* \param	dstSigB			Destination ptr for std dev
*					of 'b', may be NULL.
* \param	dstQ			Destination ptr for goodness
*					of fit, may be NULL.
*/
AlgError	AlgLinearFitIdx1D(double *datXA, double *datYA,
				  int *idxXA, int *idxYA, int idxASz,
				  double *dstA, double *dstB,
				  double *dstSigA, double *dstSigB,
				  double *dstQ)
{
  int		cnt;
  double	tD0,
  		a = 0.0,
  		b = 0.0,
		t,
		sTSq = 0.0,
		chiSq = 0.0,
  		sX = 0.0,
		sY = 0.0,
		sigA,
		sigB,
		q;
  int		*idxXP,
  		*idxYP;
  AlgError	algErr = ALG_ERR_NONE;

  /* Check parameters. */
  if((idxASz < 3) || (datXA == NULL) || (datYA == NULL) ||
     (idxXA == NULL) || (idxXA == NULL))
  {
    algErr = ALG_ERR_FUNC;
  }
  else
  {
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      sX += *(datXA + *idxXP++);
      sY += *(datYA + *idxYP++);
    }
    tD0 = sX / idxASz;
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      t = *(datXA + *idxXP++) - tD0;
      sTSq += t * t;
      b += t * *(datYA + *idxYP++);
    }
    b /= sTSq;
    a = (sY - (sX * b)) / idxASz;
    sigA = sqrt((1.0 + ((sX * sX) / (idxASz * sTSq))) / idxASz);
    sigB = sqrt(1.0 / sTSq);
    cnt = idxASz;
    idxXP = idxXA;
    idxYP = idxYA;
    while(cnt-- > 0)
    {
      tD0 = *(datYA + *idxYP++) - a - (b * *(datXA + *idxXP++));
      chiSq += tD0 * tD0;
    }
    tD0 = sqrt(chiSq / (idxASz - 2));
    sigA *= tD0;
    sigB *= tD0;
    q = AlgGammaP(0.5 * (idxASz - 2), 0.5 * (chiSq), &algErr);
  }
  if(algErr == ALG_ERR_NONE)
  {
    if(dstA)
    {
      *dstA = a;
    }
    if(dstB)
    {
      *dstB = b;
    }
    if(dstSigA)
    {
      *dstSigA = sigA;
    }
    if(dstSigB)
    {
      *dstSigB = sigB;
    }
    if(dstQ)
    {
      *dstQ = q;
    }
  }
  return(algErr);
}

#ifdef ALG_LINEARFIT_TEST
#define ALG_LINEARFIT_TEST_BUFSZ (256)

int		main(int argc, char **argv)
{
  int		tI0,
  		ok = 0,
		done = 0,
  		datCnt = 0,
  		idxCnt = 0;
  double	a,
  		b,
		sigA,
		sigB,
		q;
  FILE		*fP = NULL;
  AlgError	algErr = ALG_ERR_NONE;
  static char	ioBuf[ALG_LINEARFIT_TEST_BUFSZ];
  static int	idxA[ALG_LINEARFIT_TEST_BUFSZ];
  static double	datXA[ALG_LINEARFIT_TEST_BUFSZ],
		datYA[ALG_LINEARFIT_TEST_BUFSZ];

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
    "A test for AlgLinearFit1D() and AlgLinearFitIdx1D() which compute\n"
    "the least squares best fit linear model to the given data, ie linear\n"
    "regression. Input data records should have the form:\n"
    "  <x> <y> [<flag>]\n"
    "where the third field is optional and indicates an indexed record\n"
    "if non-zero. If any indexed records are input only records which\n"
    "are indexed will be used, this can be used to test the indexed\n"
    "function AlgLinearFitIdx1D().\n");
  }
  else
  {
    do
    {
      if(fgets(ioBuf, ALG_LINEARFIT_TEST_BUFSZ, fP) == NULL)
      {
	done = 1;
      }
      else
      {
	ioBuf[ALG_LINEARFIT_TEST_BUFSZ - 1] = '\0';
	if(sscanf(ioBuf, "%lg %lg %d",
		  datXA + datCnt, datYA + datCnt, &tI0) == 3)
	{
	  if(tI0)
	  {
	    *(idxA + idxCnt) = datCnt;
	    ++idxCnt;
	  }
	  ++datCnt;
	}
	else if(sscanf(ioBuf, "%lg %lg",
		       datXA + datCnt, datYA + datCnt) == 2)
	{
	  ++datCnt;
	}
	else
	{
	  ok = 0;
	}
      }
    } while(ok && (datCnt < ALG_LINEARFIT_TEST_BUFSZ) && !done);
    if(!ok)
    {
      (void )fprintf(stderr, "%s: Failed to read input data!\n", *argv);
    }
    else if(datCnt >= ALG_LINEARFIT_TEST_BUFSZ)
    {
      (void )fprintf(stderr,
      		     "%s: This test is restricted to %d data records.\n",
		     ALG_LINEARFIT_TEST_BUFSZ);
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
      algErr = AlgLinearFitIdx1D(datXA, datYA, idxA, idxA, idxCnt,
      			         &a, &b, &sigA, &sigB, &q);
    }
    else
    {
      algErr = AlgLinearFit1D(datCnt, datXA, datYA, &a, &b, &sigA, &sigB, &q);
    }
    if(algErr == ALG_ERR_NONE)
    {
      (void )printf("a = % 8g +/- % 8g, b = % 8g +/- % 8g, q = % 8g\n",
      		    a, sigA, b, sigB, q);
    }
  }
  return(!ok);
}

#endif /* ALG_LINEARFIT_TEST */

/*!
* @}
*/
