#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstMatrixArithmetic3_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstMatrixArithmetic3.c
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
* \brief	Basic matrix arithmetic tests.
* \ingroup	binAlgTst
*/
#include <sys/time.h>
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

typedef enum _AlgTstMatOp
{
  ALG_TST_MAT_OP_VMUL, 		/* -m
                                   \mathbf{a} = \mathbf{B} \mathbf{c}
				   AlgMatrixVectorMul() */
  ALG_TST_MAT_OP_VMULA,         /* -a
				   \mathbf{a} = \mathbf{B} \mathbf{c} +
				                \mathbf{d}
                                   a = B c + d
				   AlgMatrixVectorMulAdd() */
  ALG_TST_MAT_OP_VMULWA,	/* -w#,#
                                   \mathbf{a} = s \mathbf{B} \mathbf{c} +
				                t \mathbf{d}
                                   AlgMatrixVectorMulWAdd() */
  ALG_TST_MAT_OP_TVMUL,		/* -t
				   \mathbf{a} = \mathbf{B}^T \mathbf{c}
                                   AlgMatrixTVectorMul */
  ALG_TST_MAT_OP_TVMULA		/* -u
				   \mathbf{a} = \mathbf{B}^T \mathbf{c} +
				                \mathbf{d}
                                   AlgMatrixTVectorMulAdd */

} AlgTstMatOp;

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;
int             main(int argc, char *argv[])
{
  int           id0,
                id1,
                option,
		usage = 0,
                ok = 1,
                sz = 8,
		timer = 0,
                verbose = 0;
  long		seed = 0;
  double	zFrac = 0.2;
  double	param[2];
  const char	*opStr = NULL;
  AlgTstMatOp	op = ALG_TST_MAT_OP_VMUL;
  AlgMatrixType aType = ALG_MATRIX_RECT;
  double	*v[3];
  AlgMatrix	a;
  struct timeval times[3];
  AlgError      errCode = ALG_ERR_NONE;
  const double	tol = 1.0e-6;
  const char    *optList = "ad:hmn:tuvw:z:LRTY";

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        op = ALG_TST_MAT_OP_VMULA;
        break;
      case 'd':
        if(sscanf(optarg, "%ld", &seed) != 1)
        {
          usage = 1;
        }
        break;
      case 'm':
        op = ALG_TST_MAT_OP_VMUL;
        break;
      case 'n':
        if((sscanf(optarg, "%d", &sz) != 1) || (sz < 1))
        {
          usage = 1;
        }
        break;
      case 't':
        op = ALG_TST_MAT_OP_TVMUL;
        break;
      case 'u':
        op = ALG_TST_MAT_OP_TVMULA;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'w':
        op = ALG_TST_MAT_OP_VMULWA;
        if(sscanf(optarg, "%lg,%lg", param, param + 1) != 2)
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
      case 'z':
        if((sscanf(optarg, "%lg", &zFrac) != 1) ||
	   (zFrac < 0.0) || (zFrac > 1.0))
        {
          usage = 1;
        }
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
	a.llr = AlgMatrixLLRNew(sz, sz, sz * sz, tol, &errCode);
        break;
      case ALG_MATRIX_RECT:
	a.rect = AlgMatrixRectNew(sz, sz, &errCode);
        break;
      case ALG_MATRIX_SYM:
	a.sym = AlgMatrixSymNew(sz, &errCode);
        break;
      default:
        break;
    }
    if(errCode == ALG_ERR_NONE)
    {
      if(((v[0] = (double *)AlcMalloc(sizeof(double) * sz)) == NULL) ||
         ((v[1] = (double *)AlcMalloc(sizeof(double) * sz)) == NULL) ||
         ((v[2] = (double *)AlcMalloc(sizeof(double) * sz)) == NULL))
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
    for(id0 = 0; id0 < sz; ++id0)
    {
      for(id1 = 0; id1 <= id0; ++id1)
      {
	double v;

	v = AlgRandUniform();
	if((v < zFrac) || (zFrac > 0.999999))
	{
	  v = 0.0;
	}
	else
	{
	  v = (v - zFrac) / (1.0 - zFrac);
	  (void )AlgMatrixSet(a, id0, id1, v);
	  if(id0 != id1)
	  {
	    (void )AlgMatrixSet(a, id1, id0, v);
	  }
	}
      }
    }
    for(id0 = 0; id0 < sz; ++id0)
    {
      *(v[1] + id0) = AlgRandUniform();
      *(v[2] + id0) = AlgRandUniform();
    }
    switch(op)
    {
      case ALG_TST_MAT_OP_VMUL:
	gettimeofday(times + 0, NULL);
        AlgMatrixVectorMul(v[0], a, v[1]);
	gettimeofday(times + 1, NULL);
	opStr = "a = B c";
	break;
      case ALG_TST_MAT_OP_VMULA:
	gettimeofday(times + 0, NULL);
	AlgMatrixVectorMulAdd(v[0], a, v[1], v[2]);
	gettimeofday(times + 1, NULL);
	opStr = "a = B c + d";
	break;
      case ALG_TST_MAT_OP_VMULWA:
	gettimeofday(times + 0, NULL);
        AlgMatrixVectorMulWAdd(v[0], a, v[1], v[2], param[0], param[1]);
	gettimeofday(times + 1, NULL);
        opStr = "a = s B c + t d";
        break;
      case ALG_TST_MAT_OP_TVMUL:
	gettimeofday(times + 0, NULL);
        AlgMatrixTVectorMul(v[0], a, v[1]);
	gettimeofday(times + 1, NULL);
	opStr = "a = B^T c";
        break;
      case ALG_TST_MAT_OP_TVMULA:
	gettimeofday(times + 0, NULL);
        AlgMatrixTVectorMulAdd(v[0], a, v[1], v[2]);
	gettimeofday(times + 1, NULL);
	opStr = "a = B^T c + d";
        break;
      default:
        break;
    }
  }
  if(ok)
  {
    if(verbose)
    {
      AlgMatrixWriteAscii(a, stdout);
      (void )printf("%s\n", opStr);
      (void )printf("=\n");
    }
    for(id0 = 0; id0 < sz; ++id0)
    {
      (void )printf("%lg\n", *(v[0] + id0));
    }
    if(timer)
    {
      timersub(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
		     "%s: Elapsed time %g\n",
                     *argv, times[2].tv_sec + (0.000001 * times[2].tv_usec));
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-a] [-d#] [-m] [-n #] [-t] [-v]\n"
    "       [-w #,#] [-u] [-z#] [-L] [-R] [-Y]\n%s",
    *argv,
    "Test for basic matrix arithmatic\n"
    "  -a  a = B c + d.\n"
    "  -d  Seed for random matrix values.\n"
    "  -h  Output this help message.\n"
    "  -m  a = B c.\n"
    "  -n  Matrix size.\n"
    "  -t  a = B^T c.\n"
    "  -u  a = B^T c + d.\n"
    "  -v  Verbose output.\n"
    "  -w  a = s B c + t d.\n"
    "  -z  Fraction of zero entries.\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -T  Time function.\n"
    "  -Y  Symetric matrices.\n");
  }
  exit(errCode);
}
