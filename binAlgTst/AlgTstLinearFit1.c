#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstLinearFit1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstLinearFit1.c
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
* \brief	Simple test for AlgLinearFit1D().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

#define ALG_LINEARFIT_TEST_BUFSZ (256)

int             main(int argc, char **argv)
{
  int           tI0,
                ok = 0,
                done = 0,
                datCnt = 0,
                idxCnt = 0;
  double        a,
                b,
                sigA,
                sigB,
                q;
  FILE          *fP = NULL;
  AlgError      algErr = ALG_ERR_NONE;
  static char   ioBuf[ALG_LINEARFIT_TEST_BUFSZ];
  static int    idxA[ALG_LINEARFIT_TEST_BUFSZ];
  static double datXA[ALG_LINEARFIT_TEST_BUFSZ],
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

