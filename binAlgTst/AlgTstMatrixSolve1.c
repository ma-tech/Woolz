#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstMatrixSolve1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstMatrixSolve1.c
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
* \brief	Simple test for Ax = b matix solvers.
* \ingroup	binAlgTst
*/
#include <sys/time.h>
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

typedef enum _AlgMatrixTestMtd
{
  ALG_MATRIX_TST_MTD_CG,
  ALG_MATRIX_TST_MTD_LU,
  ALG_MATRIX_TST_MTD_LSQR,
  ALG_MATRIX_TST_MTD_SVD
} AlgMatrixTestMtd;

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
                ok = 1,
		usage = 0,
		special = 0,
                sz = 10,
                nRpt = 1,
                itr = 1000,
		timer = 0,
                verbose = 0;
  long		seed = 0;
  double        del,
  		tol = 0.000001,
		zF = 0.9;
  double        *b0 = NULL,
                *b1 = NULL,
                *x0 = NULL,
                *i0 = NULL,
                *r0 = NULL;
  AlgMatrix     a0,
                a1,
                w;
  AlgMatrixTestMtd mtd = ALG_MATRIX_TST_MTD_CG;
  AlgMatrixType aType = ALG_MATRIX_RECT;
  AlgError      errCode = ALG_ERR_NONE;
  struct timeval times[3];
  const char    *optList = "d:n:r:t:vz:CLPQSRTUY";

  a0.core = NULL;
  a1.core = NULL;
  w.core = NULL;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'd':
        if(sscanf(optarg, "%ld", &seed) != 1)
	{
	  usage = 1;
	}
      case 'n':
        if((sscanf(optarg, "%d", &sz) != 1) || (sz < 1))
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
        if((sscanf(optarg, "%lg", &tol) != 1) || (tol < 0.0) || (tol > 0.1))
        {
          usage = 1;
        }
        break;
      case 'v':
        verbose = 1;
        break;
      case 'z':
        if((sscanf(optarg, "%lg", &zF) != 1) || (zF < 0.0) || (tol > 1.0))
        {
          usage = 1;
        }
        break;
      case 'C':
        mtd = ALG_MATRIX_TST_MTD_CG;
        break;
      case 'L':
        aType = ALG_MATRIX_LLR;
        break;
      case 'P':
        special = 1;
        break;
      case 'Q':
        mtd = ALG_MATRIX_TST_MTD_LSQR;
        break;
      case 'S':
        mtd = ALG_MATRIX_TST_MTD_SVD;
        break;
      case 'R':
        aType = ALG_MATRIX_RECT;
        break;
      case 'T':
        timer = 1;
        break;
      case 'U':
        mtd = ALG_MATRIX_TST_MTD_LU;
        break;
      case 'Y':
        aType = ALG_MATRIX_SYM;
        break;
      default:
        usage = 1;
        break;
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    AlgRandSeed(seed);
    a0 = AlgMatrixNew(aType, sz, sz, sz * sz, tol, NULL);
    a1 = AlgMatrixNew(aType, sz, sz, sz * sz, tol, NULL);
    if(verbose)
    {
      r0 = (double *)AlcMalloc(sz * sizeof(double));
    }
    w.rect = AlgMatrixRectNew(4, sz, NULL);
    b0 = (double *)AlcCalloc(sz, sizeof(double));
    b1 = (double *)AlcCalloc(sz, sizeof(double));
    x0 = (double *)AlcCalloc(sz, sizeof(double));
    i0 = (double *)AlcCalloc(sz, sizeof(double));
    for(id1 = 0; id1 < sz; ++id1)
    {
      double v;

      for(id0 = 0; id0 < id1; ++id0)
      {
	v = AlgRandUniform();
	if((v > zF) && (zF < 0.999999))
	{
	  v = 0.5 * (v - zF) / (1.0 - zF);
	  AlgMatrixSet(a0, id1, id0, v);
	  AlgMatrixSet(a0, id0, id1, v);
	}
      }
      v = ((2 * id1) % 5) * 0.2 + 1.0;
      AlgMatrixSet(a0, id1, id1, v);
      *(b0 + id1) = ((id1) % 3) * AlgRandUniform();
      *(i0 + id1) = 1.0;
    }
    if(verbose != 0)
    {
      AlgVectorCopy(b1, b0, sz);
      (void )fprintf(stderr, "A = [\n");
      (void )AlgMatrixWriteAscii(a0, stderr);
      (void )fprintf(stderr, "]\nb = [\n");
      for(id0 = 0; id0 < sz; ++id0)
      {
        if(id0 < (sz - 1))
        {
          (void )fprintf(stderr, "%g;\n",  b0[id0]);
        }
        else
        {
          (void )fprintf(stderr, "%g]\n",  b0[id0]);
        }
      }
    }
    for(rpt = 0; rpt < nRpt; ++rpt)
    {
      long tL;

      switch(mtd)
      {
        case ALG_MATRIX_TST_MTD_CG:
          del = tol;
          AlgVectorCopy(x0, i0, sz);
	  gettimeofday(times + 0, NULL);
          errCode = AlgMatrixCGSolve(a0, x0, b0, w, NULL, NULL,
                                     tol, itr, &del, &itr);
	  gettimeofday(times + 1, NULL);
          break;
        case ALG_MATRIX_TST_MTD_LSQR:
	  gettimeofday(times + 0, NULL);
          errCode = AlgMatrixSolveLSQR(a0, b0, x0,
                                       1.0e-10, 1.0e-10, 1.0e-10, itr, 0,
                                       NULL, &tL, NULL, NULL, NULL,
                                       NULL, NULL);
	  gettimeofday(times + 1, NULL);
          itr = tL;
          break;
	case ALG_MATRIX_TST_MTD_LU:
          AlgVectorCopy(x0, b0, sz);
          (void )AlgMatrixCopy(a1, a0);
	  if((special != 0) && ((sz == 3) || (sz == 4)) &&
	     (aType == ALG_MATRIX_RECT))
	  {
	    if(sz == 3)
	    {
	      gettimeofday(times + 0, NULL);
	      errCode = AlgMatrixLUSolveRaw3(a1.rect->array, x0, 1);
	      gettimeofday(times + 1, NULL);
	    }
	    else /* sz == 4 */
	    {
	      gettimeofday(times + 0, NULL);
	      errCode = AlgMatrixLUSolveRaw4(a1.rect->array, x0, 1);
	      gettimeofday(times + 1, NULL);
	    }
	  }
	  else
	  {
	    gettimeofday(times + 0, NULL);
	    errCode = AlgMatrixLUSolve(a1, x0, 1);
	    gettimeofday(times + 1, NULL);
	  }
	  break;
        case ALG_MATRIX_TST_MTD_SVD:
          AlgVectorCopy(x0, b0, sz);
          (void )AlgMatrixCopy(a1, a0);
	  gettimeofday(times + 0, NULL);
          errCode = AlgMatrixSVSolve(a1, x0, tol, NULL);
	  gettimeofday(times + 1, NULL);
          break;
        default:
          break;
      }
    }
    if(verbose != 0)
    {
      (void )fprintf(stderr, "x = [\n");
      for(id1 = 0; id1 < sz; ++id1)
      {
        (void )fprintf(stderr, "%g", *(x0 + id1));
        if(id1 < (sz - 1))
        {
          (void )fprintf(stderr, ";\n");
        }
        else
        {
          (void )fprintf(stderr, "]\n");
        }
      }
      AlgMatrixVectorMul(r0, a0, x0);
      AlgVectorSub(r0, r0, b1, sz);
      (void )fprintf(stderr, "r = [\n");
      for(id1 = 0; id1 < sz; ++id1)
      {
        (void )fprintf(stderr, "%g", *(r0 + id1));
        if(id1 < (sz - 1))
        {
          (void )fprintf(stderr, ";\n");
        }
        else
        {
          (void )fprintf(stderr, "]\n");
        }
      }
    }
    else
    {
      for(id1 = 0; id1 < sz; ++id1)
      {
        (void )printf("%g\n", *(x0 + id1));
      }
    }
    if(timer)
    {
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
		     "%s: Elapsed time = %g\n",
		     *argv, times[2].tv_sec + (0.000001 * times[2].tv_usec));
    }
  }
  else
  {
    (void )fprintf(stderr,
                   "Usage: %s [-d#] [-n#] [-r#] [-t#] [-v] [-z#]\n"
		   "       [-C] [-S] [-R] [-T] [-Y]\n%s",
                   *argv,
                   "Test for timing solution of Ax = b using the Conjugate\n"
                   "Gradient and Singular Value Decomposition algorithms.\n"
		   "  -d  Seed value.\n"
                   "  -n  Matrix size.\n"
                   "  -r  Number of repeat solutions.\n"
		   "  -t  Tollerance value.\n"
                   "  -v  Verbose output (matching MATLAB) with residuals.\n"
		   "  -z  Fraction of zero entries.\n"
                   "  -C  Use the Conjugate Gradient algorithm.\n"
                   "  -L  Use a linked list row matrix.\n"
		   "  -P  Use special version code.\n"
                   "  -Q  Use the LSQR algorithm.\n"
                   "  -S  Use the Singular Value Decomposition algorithm.\n"
                   "  -T  Time solution.\n"
                   "  -R  Use a rectangular matrix.\n"
                   "  -Y  Use a symetric matrix.\n");
  }
  AlgMatrixFree(a0);
  AlgMatrixFree(a1);
  AlgMatrixFree(w);
  AlcFree(b0);
  AlcFree(b1);
  AlcFree(i0);
  AlcFree(r0);
  AlcFree(x0);
  exit(errCode);
}

