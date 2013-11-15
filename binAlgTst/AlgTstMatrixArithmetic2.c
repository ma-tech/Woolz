#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstMatrixArithmetic2_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstMatrixArithmetic2.c
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
  ALG_TST_MAT_OP_ADD,
  ALG_TST_MAT_OP_MUL,
  ALG_TST_MAT_OP_SUB
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
		id2,
                option,
		usage = 0,
                ok = 1,
                sz = 8,
		timer = 0,
                verbose = 0;
  long		seed = 0;
  double	zFrac = 0.2;
  const char	*opStr = NULL;
  AlgTstMatOp	op = ALG_TST_MAT_OP_ADD;
  AlgMatrixType aType = ALG_MATRIX_RECT;
  AlgMatrix	a[3];
  AlgError      errCode = ALG_ERR_NONE;
  struct timeval times[3];
  const double	tol = 1.0e-6;
  const char    *optList = "ad:hmn:svz:LRTY";

  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        op = ALG_TST_MAT_OP_ADD;
        break;
      case 'd':
        if(sscanf(optarg, "%ld", &seed) != 1)
        {
          usage = 1;
        }
        break;
      case 'm':
        op = ALG_TST_MAT_OP_MUL;
        break;
      case 'n':
        if((sscanf(optarg, "%d", &sz) != 1) || (sz < 1))
        {
          usage = 1;
        }
        break;
      case 's':
        op = ALG_TST_MAT_OP_SUB;
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
	a[2].llr = AlgMatrixLLRNew(sz, sz, sz * sz, tol, &errCode);
        break;
      case ALG_MATRIX_RECT:
	a[0].rect = AlgMatrixRectNew(sz, sz, &errCode);
	a[1].rect = AlgMatrixRectNew(sz, sz, &errCode);
	a[2].rect = AlgMatrixRectNew(sz, sz, &errCode);
        break;
      case ALG_MATRIX_SYM:
	a[0].sym = AlgMatrixSymNew(sz, &errCode);
	a[1].sym = AlgMatrixSymNew(sz, &errCode);
	switch(op)
	{
          case ALG_TST_MAT_OP_MUL:
	    a[2].rect = AlgMatrixRectNew(sz, sz, &errCode);
	    break;
	  default:
	    a[2].sym = AlgMatrixSymNew(sz, &errCode);
	    break;
	}
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
	for(id2 = 0; id2 < 2; ++id2)
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
	  (void )AlgMatrixSet(a[id2], id0, id1, v);
	  if(id0 != id1)
	  {
	    (void )AlgMatrixSet(a[id2], id1, id0, v);
	  }
	}
      }
    }
    switch(op)
    {
      case ALG_TST_MAT_OP_ADD:
	gettimeofday(times + 0, NULL);
        AlgMatrixAdd(a[2], a[0], a[1]);
	gettimeofday(times + 1, NULL);
	opStr = "+";
	break;
      case ALG_TST_MAT_OP_SUB:
	gettimeofday(times + 0, NULL);
        AlgMatrixSub(a[2], a[0], a[1]);
	gettimeofday(times + 1, NULL);
	opStr = "-";
	break;
      case ALG_TST_MAT_OP_MUL:
	gettimeofday(times + 0, NULL);
        AlgMatrixMul(a[2], a[0], a[1]);
	gettimeofday(times + 1, NULL);
	opStr = "*";
        break;
      default:
        break;
    }
  }
  if(ok)
  {
    if(verbose)
    {
      AlgMatrixWriteAscii(a[0], stdout);
      (void )printf("%s\n", opStr);
      AlgMatrixWriteAscii(a[1], stdout);
      (void )printf("=\n");
    }
    AlgMatrixWriteAscii(a[2], stdout);
    if(timer)
    {
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
		     "%s: Elapsed time %g\n",
                     *argv, times[2].tv_sec + (0.000001 * times[2].tv_usec));
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-a] [-d#] [-m] [-n #] [-s] [-v] [-z#]\n"
    "       [-L] [-R] [-Y]\n%s",
    *argv,
    "Test for basic matrix arithmatic\n"
    "  -a  Add matrices.\n"
    "  -d  Seed for random matrix values.\n"
    "  -h  Output this help message.\n"
    "  -m  Multiply matrices.\n"
    "  -n  Matrix size.\n"
    "  -s  Subtract matrices.\n"
    "  -v  Verbose output.\n"
    "  -z  Fraction of zero entries.\n"
    "  -L  Sparse (Linked List Row) matrices\n"
    "  -R  Rectangular matrices.\n"
    "  -T  Time function.\n"
    "  -Y  Symetric matrices.\n");
  }
  exit(errCode);
}
