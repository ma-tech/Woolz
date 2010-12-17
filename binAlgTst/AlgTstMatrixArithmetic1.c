#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstMatrixArithmetic1_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstMatrixArithmetic1.c
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
* \brief	Basic matrix arithmetic tests.
* \ingroup	binAlgTst
*/
#include <sys/time.h>
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

typedef enum _AlgTstMatOp
{
  ALG_TST_MAT_OP_COPY,
  ALG_TST_MAT_OP_SCALE,
  ALG_TST_MAT_OP_TRANSPOSE
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
		inPlace = 0,
		usage = 0,
                ok = 1,
                sz = 8,
		timer = 0,
                verbose = 0;
  long		seed = 0;
  double	param = 0.0,
  		zFrac = 0.2;
  const char	*opStr;
  AlgTstMatOp	op = ALG_TST_MAT_OP_TRANSPOSE;
  AlgMatrixType aType = ALG_MATRIX_RECT;
  AlgMatrix	a[2];
  AlgError      errCode = ALG_ERR_NONE;
  struct timeval times[3];
  const double	tol = 1.0e-6;
  const char    *optList = "cd:hn:s:tvz:ILRTY";

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
        op = ALG_TST_MAT_OP_COPY;
        break;
      case 'd':
        if(sscanf(optarg, "%ld", &seed) != 1)
        {
          usage = 1;
        }
        break;
      case 'n':
        if((sscanf(optarg, "%d", &sz) != 1) || (sz < 1))
        {
          usage = 1;
        }
        break;
      case 's':
        op = ALG_TST_MAT_OP_SCALE;
        if(sscanf(optarg, "%lg", &param) != 1)
        {
          usage = 1;
        }
        break;
      case 't':
        op = ALG_TST_MAT_OP_TRANSPOSE;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'I':
        inPlace = 1;
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
	a[0].llr = AlgMatrixLLRNew(sz, sz, sz * sz, tol, &errCode);
	a[1].llr = AlgMatrixLLRNew(sz, sz, sz * sz, tol, &errCode);
	a[0].llr->tol = a[1].llr->tol = tol;
        break;
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
	}
	(void )AlgMatrixSet(a[0], id0, id1, v);
	if(id0 != id1)
	{
	  (void )AlgMatrixSet(a[0], id1, id0, v);
	}
      }
    }
    switch(op)
    {
      case ALG_TST_MAT_OP_COPY:
	opStr = "copy";
	gettimeofday(times + 0, NULL);
        AlgMatrixCopy(a[1], a[0]);
	gettimeofday(times + 1, NULL);
	break;
      case ALG_TST_MAT_OP_SCALE:
	opStr = "transpose";
	if(inPlace)
	{
          AlgMatrixCopy(a[1], a[0]);
	  gettimeofday(times + 0, NULL);
          AlgMatrixScale(a[1], a[1], param);
	  gettimeofday(times + 1, NULL);
	}
	else
	{
	  gettimeofday(times + 0, NULL);
	  AlgMatrixScale(a[1], a[0], param);
	  gettimeofday(times + 1, NULL);
	}
	break;
      case ALG_TST_MAT_OP_TRANSPOSE:
	opStr = "transpose";
	gettimeofday(times + 0, NULL);
        AlgMatrixTranspose(a[1], a[0]);
	gettimeofday(times + 1, NULL);
	break;
      default:
        break;
    }
  }
  if(ok)
  {
    if(verbose)
    {
      (void )printf("%s\n(\n", opStr);
      AlgMatrixWriteAscii(a[0], stdout);
      (void )printf(")\n=\n", opStr);
    }
    AlgMatrixWriteAscii(a[1], stdout);
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
    "Usage: %s [-h] [-d#] [-i] [-n #] [-s#] [-t] [-v] [-z#]\n"
    "       [-L] [-R] [-Y]\n%s",
    *argv,
    "Test for basic matrix arithmatic\n"
    "  -d  Seed for random matrix values.\n"
    "  -h  Output this help message.\n"
    "  -n  Matrix size.\n"
    "  -s  Scale (multiply be scalar).\n"
    "  -t  Transpose matrix.\n"
    "  -v  Verbose output.\n"
    "  -z  Fraction of zero entries.\n"
    "  -I  In place (not always appropriate)\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -T  Time function.\n"
    "  -Y  Symetric matrices.\n");
  }
  exit(errCode);
}
