#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstQSort1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstQSort1.c
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
* \brief	Simple test for AlgQSort().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

static int                      AlgQSortTestSortFn(
                                  const void *cData,
                                  const void *a,
                                  const void *b);
extern char     *optarg;
extern int      optind,
                optind,
                opterr;

int             main(int argc, char *argv[])
{
  int           idN,
                option,
                nElm = 20,
                elmSz,
                ok = 1,
                reverse = 0,
                usage = 0,
                verbose = 0;
  int           *ary = NULL;
  static char   optList[] = "hrvn:";

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
        if(sscanf(optarg, "%d", &nElm) != 1)
        {
          usage = 1;
        }
        break;
      case 'r':
        reverse = 1;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'h':
      default:
        usage = 1;
        break;
    }
  }
  ok = !usage;
  if(ok)
  {
    elmSz = sizeof(int);
    if((ary = AlcMalloc(elmSz * nElm)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Filed to allocate data array.\n",
                     *argv);
    }
  }
  if(ok)
  {
    for(idN = 0; idN < nElm; ++idN)
    {
      ary[idN] = (int )(10.0 * nElm * AlgRandUniform());
    }
    if(verbose)
    {
      for(idN = 0; idN < nElm; ++idN)
      {
        (void )printf("% 8d % 8d\n", idN, ary[idN]);
      }
    }
    (void )printf("\n");
    AlgQSort(ary, nElm, elmSz, (void *)reverse, AlgQSortTestSortFn);
    if(verbose)
    {
      for(idN = 0; idN < nElm; ++idN)
      {
        (void )printf("% 8d % 8d\n", idN, ary[idN]);
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
                   "Usage: %s [-ch [-v] [-n#]\n%s",
                   "Options:\n"
                   "  -h  Show this usage message.\n"
                   "  -r  Reverse order.\n"
                   "  -v  Be verbose.\n"
                   "  -n  Number of data to sort.\n"
                   "Test for AlgQSort().\n",
                   *argv);
  }
  return(!ok);
}

static int      AlgQSortTestSortFn(const void *cData,
                                   const void *a, const void *b)
{
  int           cmp;

  if(cData)
  {
    cmp = *(int *)b - *(int *)a;
  }
  else
  {
    cmp = *(int *)a - *(int *)b;
  }
  return(cmp);
}

