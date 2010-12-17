#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstMatrixCGSolve1_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstMatrixCGSolve1.c
* \author       Bill Hill
* \date         November 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Simple test for AlgMatrixCGSolve().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		id0,
		option,
		ok = 1,
		usage = 0,
  		itr = 1000;
  double	tol = 0.000001;
  AlgMatrix	a,
  		w;
  double	*b = NULL,
  		*x = NULL;
  AlgError	errCode = ALG_ERR_NONE;
  AlgMatrixType aType = ALG_MATRIX_RECT;
  const char	*optList = "hLRY";

  a.core = w.core = NULL;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
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
  if(ok)
  {
    switch(aType)
    {
      case ALG_MATRIX_LLR:
        a.llr = AlgMatrixLLRNew((size_t )3, (size_t )3, 9, tol, &errCode);
	break;
      case ALG_MATRIX_RECT:
        a.rect = AlgMatrixRectNew((size_t )3, (size_t )3, &errCode);
	break;
      case ALG_MATRIX_SYM:
        a.sym = AlgMatrixSymNew((size_t )3, &errCode);
	break;
      default:
        break;
    }
    w.rect = AlgMatrixRectNew((size_t )4, (size_t )3, &errCode);
    if((a.core != NULL) && (w.core != NULL))
    {
      b = (double *)AlcMalloc(3 * sizeof(double));
      x = (double *)AlcMalloc(3 * sizeof(double));
    }
    if((b == NULL) || (x == NULL))
    {
      (void )fprintf(stderr, "%s: Failed to allocate matrices\n", *argv);
      ok = 0;
    }
  }
  if(ok)
  {
    AlgMatrixSet(a, 0, 0, 6.0); 
    AlgMatrixSet(a, 0, 1, 1.0); 
    AlgMatrixSet(a, 0, 2, 0.0); 
    AlgMatrixSet(a, 1, 0, 1.0); 
    AlgMatrixSet(a, 1, 1, 5.0); 
    AlgMatrixSet(a, 1, 2, 1.0); 
    AlgMatrixSet(a, 2, 0, 0.0); 
    AlgMatrixSet(a, 2, 1, 1.0); 
    AlgMatrixSet(a, 2, 2, 8.0); 
    b[0] = 8.0;
    b[1] = 14.0;
    b[2] = 26.0;
    x[0] = 1.0;
    x[1] = 1.0;
    x[2] = 1.0;
    errCode = AlgMatrixCGSolve(a, x, b, w, NULL, NULL, tol, itr, &tol, &itr);
    if(errCode != NULL)
    {
      (void )fprintf(stderr, "%s: Failed to solve using CG, error code = %d\n",
                     *argv, (int )errCode);
      ok = 0;
    }
  }
  if(ok)
  {
    (void )printf("itr = %d\n", itr);
    (void )printf("tol = %g\n", tol);
    for(id0 = 0; id0 < 3; ++id0)
    {
      printf("%g\n", x[id0]);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-L] [-R] [-Y]\n%s",
    *argv,
    "Test for AlgMatrixCGSolve().\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -Y  Symetric matrices.\n");
  }
  exit(ok != 0);
}
