#pragma ident "MRC HGU $Id$"
/*!
* \file         AlgMatrixRSEigen.c
* \author       Bill Hill
* \date         May 2001
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        Functions to find the eigenvalues and eigenvectors of a
*		real symmetric matrix.
* \ingroup      AlgMatrix
* \bug          None known.
* \todo
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <Alg.h>

static void			AlgMatrixRSEigenSort(
				  double **vM,
				  double *xM,
				  int n,
				  int reqEV);

/*!
* \return       Error code.
* \ingroup      AlgMatrix
* \brief        Determines the eigenvalues and eigenvectors of a
*		real symmetric matrix by calling AlgMatrixRSTDiag()
*		to create a tridiagonal symetric matrix and then
*		AlgMatrixTDiagQLI() to compute its eigenvalues and
*		eigenvectors. The eigenvectors and eigenvalues 
*		are returned in descending eigenvalue order.
*		For efficiency, the eigenvectors should only be
*		computed if required.
* \param        aM 			Given real symmetric matrix
*					which contains the eigenvectors
*					in it's columns on return if
*					required.
* \param        aSz 		        Size of the array.
* \param	vM 			Given array for the return of the
* 					eigenvalues.
* \param	reqEV			Non zero if the eigenvectors are
*					required.
*/
AlgError	AlgMatrixRSEigen(double **aM, int aSz, double *vM, int reqEV)
{
  double	*oM = NULL;
  AlgError	errCode = ALG_ERR_NONE;


  if((aM == NULL) || (*aM == NULL) || (aSz < 2) || (vM == NULL))
  {
    errCode = ALG_ERR_FUNC;
  }
  else if((oM = (double *)AlcMalloc(sizeof(double) * aSz)) == NULL)
  {
    errCode = ALG_ERR_MALLOC;
  }
  if(errCode == ALG_ERR_NONE)
  {
    if((errCode = AlgMatrixRSTDiag(aM, aSz, vM, oM)) == ALG_ERR_NONE)
    {
      errCode = AlgMatrixTDiagQLI(vM, oM, aSz, reqEV? aM: NULL);
    }
  }
  if(errCode == ALG_ERR_NONE)
  {
    AlgMatrixRSEigenSort(aM, vM, aSz, reqEV);
  }
  if(oM)
  {
    AlcFree(oM);
  }
  return(errCode);
}

/*!
* \return       <void>
* \ingroup      AlgMatrix
* \brief	Sorts the eigenvectors and eigenvalues into descending
* 		eigenvalue order. Because AlgMatrixRSTDiag() runs in
*		O(N^3) an O(N^2) insertion sort algorithm is used
*		without it significantly impacting on the execution
*		time, but it shouldn't be used elsewhere!
*		This function is based on 'eigsrt' in Numerical Recipies
*		in C: The Art of ScientificComputing. Cambridge University
*		Press, 1992.
*		in 
* \param	vM			The eigenvectors.
* \param	xM			The eigenvalues.
* \param	n			The number of eigenvectors or
* 					eigenvalues.
* \param	reqEV			Non zero if the eigenvectors are
*					required.
*/
static void	AlgMatrixRSEigenSort(double **vM, double *xM, int n, int reqEV)
{
  int		id0,
  		id1,
		id2,
		n1;
  double 	keyVal;

  n1 = n - 1;
  for(id0 = 0; id0 < n1; ++id0)
  {
    keyVal = xM[id2 = id0];
    for(id1 = id0 + 1; id1 < n; ++id1)
    {
      if(xM[id1] >= keyVal)
      {
	keyVal = xM[id2 = id1];
      }
    }
    if(id2 != id0)
    {
      xM[id2] = xM[id0];
      xM[id0] = keyVal;
      if(reqEV)
      {
	for(id1 = 0; id1 < n; id1++)
	{
	  keyVal = vM[id1][id0];
	  vM[id1][id0] = vM[id1][id2];
	  vM[id1][id2] = keyVal;
	}
      }
    }
  }
}

#ifdef ALG_MATRIXRSEIGEN_TEST
extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idC,
  		idR,
		idV,
  		nC = 0,
		nR = 0,
		option,
		ok = 1,
		usage = 0,
		reqEV = 1;
  double	**aM = NULL;
  double	*dP0,
  		*vM = NULL;
  static char	optList[] = "hn";

  opterr = 0;
  /* Parse command line. */
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
        reqEV = 0;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  /* Read matrix from stdin. */
  if(AlcDouble2ReadAsci(stdin, &aM, &nR, &nC) != ALC_ER_NONE)
  {
    ok = 0;
    (void )fprintf(stderr, "%s: Failed to read matrix\n", *argv);
  }
  if(ok)
  {
    if(nR != nC)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Input matrix is not square\n", *argv);
    }
    else if(nR < 2)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Input matrix smaller than 2x2\n", *argv);
    }
  }
  /* Create space for the eigen values. */
  if(ok)
  {
    if((vM = (double *)AlcMalloc(sizeof(double) * nR)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to allocate storage.\n", *argv);
    }
  }
  /* Find the eigenvalues and eigenvectors. */
  if(ok)
  {
    if(AlgMatrixRSEigen(aM, nR, vM, reqEV) != ALG_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Failed to compute eigenvalues and eigenvectors.\n",
		     *argv);
    }
  }
  if(aM)
  {
    AlcDouble2Free(aM);
  }
  if(vM)
  {
    AlcFree(vM);
  }
  /* Output the eigenvalues and eigenvectors. */
  if(ok)
  {
    for(idR = 0, idV = 0; idR < nR; ++idR)
    {
      (void )printf("%g\n", vM[idR]);
    }
    if(reqEV)
    {
      for(idR = 0, idV = 0; idR < nR; ++idR)
      {
	for(idC = 0; idC < nC; ++idC)
	{
	  (void )printf("%g ", aM[idR][idC]);
	}
	(void )printf("\n");
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-n]\n%s",
    *argv,
    "Options:\n"
    "  -h  Show this usage message.\n"
    "  -n  Dont compute the eigen vectors.\n"
    "Test for AlgMatrixRSEigen(). Reads white space separated, ascii\n"
    "encoded, floating point matrix values from the standard input and\n"
    "then prints the matrix's eigen values followed by it's eigen vectors.\n");
  }
  return(!ok);
}
#endif /* ALG_MATRIXRSEIGEN_TEST */
