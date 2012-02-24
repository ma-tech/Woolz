#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstMatrixCGSolve2_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstMatrixCGSolve2.c
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
* \brief	Simple test for AlgMatrixCGSolve().
* \ingroup	binAlgTst
*/
#include <sys/time.h>
#include <stdio.h>
#include <float.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           id0,
                id1,
                rpt,
		option,
		ok = 0,
                nRpt = 1,
                itr = 1000,
		timer = 0,
		usage = 0,
		verbose = 0;
  long          seed = 0;
  size_t        sz = 8;
  double        del,
  		tol = 0.000001,
		zF = 0.9;
  AlgMatrix	w;
  AlgMatrix	a[2];
  double        *b0,
                *b1,
                *x0,
                *x1,
                *i0;
  AlgMatrixType	aType = ALG_MATRIX_RECT;
  AlgError      errCode = ALG_ERR_NONE;
  struct timeval times[3];
  const char	*optList = "d:hi:n:r:t:vz:LRTY";

  a[0].core = a[1].core = NULL;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        if(sscanf(optarg, "%ld", &seed) != 1)
	{
	  usage = 1;
	}
	break;
      case 'i':
        if((sscanf(optarg, "%d", &itr) != 1) || (itr < 0))
	{
	  usage = 1;
	}
	break;
      case 'n':
        if((sscanf(optarg, "%zd", &sz) != 1) || (sz < 1))
	{
	  usage = 1;
	}
	break;
      case 'r':
        if((sscanf(optarg, "%d", &nRpt) != 1) || (nRpt < 0))
	{
	  usage = 1;
	}
	break;
      case 't':
        if((sscanf(optarg, "%lg", &tol) != 1) || (tol < 0.0))
	{
	  usage = 1;
	}
	break;
      case 'v':
        verbose = 1;
	break;
      case 'z':
        if((sscanf(optarg, "%lg", &zF) != 1) || (zF < 0.0) || (zF > 1.0))
	{
	  usage = 1;
	}
	break;
      case 'L':
	aType = ALG_MATRIX_LLR;
	break;
      case 'R':
	aType = ALG_MATRIX_RECT;
	break;
      case 'T':
	timer = 1;
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
    size_t nInElm;

    switch(aType)
    {
      case ALG_MATRIX_LLR:
        nInElm = (size_t )
	         ((floor((double )sz * (double )sz * (1.0 - zF)) + 1.0) * 1.5);
	if(nInElm > sz * sz)
	{
	  nInElm = sz * sz;
	}
        a[0].llr = AlgMatrixLLRNew(sz, sz, nInElm, tol, &errCode);
	a[1].llr = AlgMatrixLLRNew(sz, sz, nInElm, tol, &errCode);
      case ALG_MATRIX_RECT:
        a[0].rect = AlgMatrixRectNew(sz, sz, &errCode);
        a[1].rect = AlgMatrixRectNew(sz, sz, &errCode);
        break;
      case ALG_MATRIX_SYM:
        a[0].sym = AlgMatrixSymNew(sz, &errCode);
        a[1].sym = AlgMatrixSymNew(sz, &errCode);
        break;
      default:
        break;
    }
    if(errCode == ALG_ERR_NONE)
    {
      w.rect = AlgMatrixRectNew(4, sz, &errCode);
    }
    if(errCode == ALG_ERR_NONE)
    {
      if(((b0 = (double *)AlcMalloc(sz * sizeof(double))) == NULL) ||
         ((b1 = (double *)AlcMalloc(sz * sizeof(double))) == NULL) ||
         ((x0 = (double *)AlcMalloc(sz * sizeof(double))) == NULL) ||
         ((x1 = (double *)AlcMalloc(sz * sizeof(double))) == NULL) ||
         ((i0 = (double *)AlcMalloc(sz * sizeof(double))) == NULL))
      {
        errCode = ALG_ERR_MALLOC;
      }
    }
    if(errCode != ALG_ERR_NONE)
    {
      (void )fprintf(stderr, "%s: Failed to allocate matrices\n", *argv);
      ok = 0;
    }
  }
  if(ok)
  {
    AlgRandSeed(seed);
    for(id1 = 0; id1 < sz; ++id1)
    {
      for(id0 = 0; id0 < id1; ++id0)
      {
	double 	v;

	v = AlgRandUniform();
	if((v < zF) || (zF > 0.999999))
	{
	  if(id0 == id1)
	  {
	    (void )AlgMatrixSet(a[0], id0, id0, 1.0);
	  }
	}
	else
	{
	  v = 0.5 * (v - zF) / (1.0 - zF);
	  (void )AlgMatrixSet(a[0], id0, id1, v);
	  (void )AlgMatrixSet(a[0], id1, id0, v);
	}
      }
      (void )AlgMatrixSet(a[0], id1, id1, 1.0);
      *(b0 + id1) = AlgRandUniform();
      *(i0 + id1) = 1.0;
    }
    if(verbose)
    {
      (void )printf("Ax = b, solution using AlgMatrixCGSolve():\n");
      (void )printf("A = {\n");
      AlgMatrixWriteAscii(a[0], stdout);
      (void )printf("}\n");
      (void )printf("b = {\n");
      for(id0 = 0; id0 < sz; ++id0)
      {
        (void )printf("%lg\n", b0[id0]);
      }
      (void )printf("}\n");
      (void )printf("x (initial) = {\n");
      for(id0 = 0; id0 < sz; ++id0)
      {
        (void )printf("%lg\n", i0[id0]);
      }
      (void )printf("}\n");
    }
    for(rpt = 0; rpt < nRpt; ++rpt)
    {
      AlgVectorCopy(x0, i0, sz);
      AlgVectorCopy(x1, b0, sz);
      AlgMatrixCopy(a[1], a[0]);
      del = tol;
      gettimeofday(times + 0, NULL);
      errCode = AlgMatrixCGSolve(a[1], x0, b0, w, NULL, NULL, tol, itr,
				 &del, &itr);
      gettimeofday(times + 1, NULL);
      if(errCode != ALG_ERR_NONE)
      {
	(void )fprintf(stderr, "%s: Failed to solve (%d).\n",
	               *argv, (int )errCode);
	ok = 0;
	break;
      }
      for(id1 = 0; id1 < sz; ++id1)
      {
	double v;

	v = *(x0 + id1) - *(x1 + id1);
	del += v * v;
      }
      del = (del > DBL_EPSILON)? sqrt(del) / sz: 0.0;
    }
  }
  if(ok)
  {
    if(verbose)
    {
      (void )printf("x (solution) = {\n");
      for(id0 = 0; id0 < sz; ++id0)
      {
        (void )printf("%lg\n", x0[id0]);
      }
      (void )printf("}\n");
    }
    (void )printf("%s: RMS = %g\n", *argv, del);
    (void )printf("%s: itr = %d\n", *argv, itr);
    if(timer)
    {
      timersub(times + 1, times + 0, times + 2);
      (void )printf("%s: Elapsed time = %g\n",
                    *argv, times[2].tv_sec + (0.000001 * times[2].tv_usec));
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-d#] [-i#] [-n#] [-r #] [-t#] [-z#]\n"
    "       [-L] [-R] [-T] [-Y]\n%s",
    *argv,
    "  -d  Seed for random matrix values.\n"
    "  -h  Output this help message.\n"
    "  -i  Output this help message.\n"
    "  -n  Matrix size.\n"
    "  -r  Number of repeats.\n"
    "  -t  Tollerance value.\n"
    "  -v  Verbose output.\n"
    "  -z  Fraction of zero entries.\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -T  Time function.\n"
    "  -Y  Symetric matrices.\n");
  }
  exit(ok != 0);
}
