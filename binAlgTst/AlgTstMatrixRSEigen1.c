#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstMatrixRSEigen1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstMatrixRSEigen1.c
* \author       Bill Hill
* \date         November 2010
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
* \brief	Simple test for AlgMatrixRSEigen().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

int             main(int argc, char *argv[])
{
  int           iR,
                option,
                ok = 1,
		builtInEx = 0,
                usage = 0,
		verbose = 0,
                reqEV = 0;
  double        *vM = NULL;
  AlgMatrix	a;
  AlgMatrixType	aType = ALG_MATRIX_RECT;
  AlgError	errNum = ALG_ERR_NONE;
  static char   optList[] = "behvLRY";
  const double	tol = 0.000001;

  opterr = 0;
  a.core = NULL;
  /* Parse command line. */
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        builtInEx = 1;
        break;
      case 'e':
        reqEV = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'L':
        aType = ALG_MATRIX_LLR;
	break;
      case 'R':
        aType = ALG_MATRIX_RECT;
	break;
      case 'Y':
        aType = ALG_MATRIX_SYM;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
        break;
    }
  }
  ok = (usage == 0);
  /* Read matrix from stdin. */
  if(ok)
  {
    if(builtInEx)
    {
      a = AlgMatrixNew(aType, 3, 3, 9, tol, &errNum);
      AlgMatrixSet(a, 0, 0,  0.0);
      AlgMatrixSet(a, 0, 1,  1.0);
      AlgMatrixSet(a, 0, 2, -1.0);
      AlgMatrixSet(a, 1, 0,  1.0);
      AlgMatrixSet(a, 1, 1,  1.0);
      AlgMatrixSet(a, 1, 2,  0.0);
      AlgMatrixSet(a, 2, 0, -1.0);
      AlgMatrixSet(a, 2, 1,  0.0);
      AlgMatrixSet(a, 2, 2,  1.0);
    }
    else
    {
      a = AlgMatrixReadAscii(aType, tol, stdin, " \t\n\r", 4096, NULL);
      if(a.core == NULL)
      {
	ok = 0;
	(void )fprintf(stderr, "%s: Failed to read matrix\n", *argv);
      }
    }
  }
  if(ok)
  {
    if(a.core->nR != a.core->nC)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Input matrix is not square\n", *argv);
    }
    else if(a.core->nR < 2)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Input matrix smaller than 2x2\n", *argv);
    }
  }
  if(ok && verbose)
  {
    (void )printf("Matrix:\n");
    (void )AlgMatrixWriteAscii(a, stdout);
  }
  /* Create space for the eigen values. */
  if(ok)
  {
    if((vM = (double *)AlcMalloc(sizeof(double) * a.core->nR)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to allocate storage.\n", *argv);
    }
  }
  /* Find the eigenvalues and eigenvectors. */
  if(ok)
  {
    if(AlgMatrixRSEigen(a, vM, reqEV) != ALG_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to compute eigenvalues and eigenvectors.\n",
                     *argv);
    }
  }
  /* Output the eigenvalues and eigenvectors. */
  if(ok)
  {
    if(verbose)
    {
      if(builtInEx)
      {
	(void )printf("Eigen values should be:\n"
		      "  -1, 1, 2\n");
      }
      (void )printf("Eigen values:\n");
    }
    for(iR = 0; iR < a.core->nR; ++iR)
    {
      (void )printf("%g\n", vM[iR]);
    }
    if(reqEV)
    {
      if(verbose)
      {
	if(builtInEx)
	{
	(void )printf("Eigen vectors should be:\n"
                      "   8.01650  0.00000 -0.57735\n"
                      "  -0.40825  0.70711 -0.57735\n"
                      "   0.40825  0.70711  0.57735\n");
	}
	(void )printf("Eigen vectors:\n");
      }
      (void )AlgMatrixWriteAscii(a, stdout);
    }
  }
  AlgMatrixFree(a);
  if(vM)
  {
    AlcFree(vM);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-n]\n%s",
    *argv,
    "Options:\n"
    "  -b  Use built in example.\n"
    "  -h  Show this usage message.\n"
    "  -e  Compute the eigen vectors.\n"
    "  -v  Verbose output.\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -Y  Symetric matrices.\n"
    "Test for AlgMatrixRSEigen(). Reads white space separated, ascii\n"
    "encoded, floating point matrix values from the standard input and\n"
    "then prints the matrix's eigen values followed by it's eigen vectors.\n");
  }
  return(!ok);
}
