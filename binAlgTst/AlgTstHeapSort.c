#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstHeapSort1_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstHeapSort1.c
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
* \brief	Simple test for AlgHeapSort().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

#define ALG_HEAPSORT_TEST_DATA_SZ (17)

int             main(int argc, char *argv[])
{
  int           idI,
                option,
                ok = 1,
                usage = 0,
                elmSz,
                dataRange,
                nData = ALG_HEAPSORT_TEST_DATA_SZ,
                doChar = 0,
                doIdx = 0,
                doRand = 0;
  int           *dataI;
  char          *dataC;
  void          *data;
  int           *idx;
  static char   optList[] = "cirs:";

  opterr = 0;
  elmSz = sizeof(int);
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
        doChar = 1;
        elmSz = sizeof(char);
        break;
      case 'i':
        doIdx = 1;
        break;
      case 'r':
        doRand = 1;
        break;
      case 's':
        if((sscanf(optarg, "%d", &nData) != 1) || (nData < 0))
        {
          usage = 1;
          ok = 0;
	}
        break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
        ok = 0;
        break;
    }
  }
  if(ok)
  {
    dataRange = nData;
    if(doChar)
    {
      if(dataRange > 127)
      {
        dataRange = 127;
      }
      if(((dataC = (char *)AlcMalloc(nData * elmSz)) == NULL) ||
         (doIdx && ((idx = (int *)AlcMalloc(nData * sizeof(int))) == NULL)))
      {
        ok = 0;
      }
      data = dataC;
    }
    else
    {
      if(((dataI = (int *)AlcMalloc(nData * elmSz)) == NULL) ||
         (doIdx && ((idx = (int *)AlcMalloc(nData * sizeof(int))) == NULL)))
      {
        ok = 0;
      }
      data = dataI;
    }
    if(ok == 0)
    {
      (void )fprintf(stderr,
                     "%s: failed to allocate data storage.\n",
                     *argv);
    }
  }
  if(ok)
  {
    if(doIdx)
    {
      if(doChar)
      {
        for(idI = 0; idI < nData; ++idI)
        {
          idx[idI] = idI;
          if(doRand)
          {
            dataC[idI] = rand() % dataRange;
          }
          else
          {
            dataC[idI] = (nData - idI - 1) % dataRange;
          }
          printf("%3d %3d %3d\n", idI, idx[idI], dataC[idx[idI]]);
        }
      }
      else
      {
        for(idI = 0; idI < nData; ++idI)
        {
          idx[idI] = idI;
          if(doRand)
          {
            dataI[idI] = rand() % dataRange;
          }
          else
          {
            dataI[idI] = idI % dataRange;
            /* dataI[idI] = (nData - idI - 1) % dataRange; */
          }
          printf("%3d %3d %3d\n", idI, idx[idI], dataI[idx[idI]]);
        }
      }
    }
    else
    {
      if(doChar)
      {

        for(idI = 0; idI < nData; ++idI)
        {
          if(doRand)
          {
            dataC[idI] = rand() % nData;
          }
          else
          {
            dataC[idI] = (nData - idI - 1) % dataRange;
          }
          printf("%3d %3d\n", idI, dataC[idI]);
        }
      }
      else
      {
        for(idI = 0; idI < nData; ++idI)
        {
          if(doRand)
          {
            dataI[idI] = rand() % nData;
          }
          else
          {
            dataI[idI] = (nData - idI - 1) % dataRange;
          }
          printf("%3d %3d\n", idI, dataI[idI]);
        }
      }
    }
    if(doIdx)
    {
      if(doChar)
      {
        (void )AlgHeapSortIdx((char *)data, idx, nData, AlgHeapSortCmpIdxCFn);
      }
      else
      {
        (void )AlgHeapSortIdx((char *)data, idx, nData, AlgHeapSortCmpIdxIFn);
      }
    }
    else
    {
      if(doChar)
      {
        (void )AlgHeapSort((char *)data, nData, elmSz, AlgHeapSortCmpCFn);
      }
      else
      {
        (void )AlgHeapSort((char *)data, nData, elmSz, AlgHeapSortCmpIFn);
      }
    }
     printf("AlgHeapSort\n");
    if(doIdx)
    {
      if(doChar)
      {
        for(idI = 0; idI < nData; ++idI)
        {
          printf("%3d %3d %3d\n", idI, idx[idI], dataC[idx[idI]]);
        }
      }
      else
      {
        for(idI = 0; idI < nData; ++idI)
        {
          printf("%3d %3d %3d\n", idI, idx[idI], dataI[idx[idI]]);
        }
      }
    }
    else
    {
      if(doChar)
      {
        for(idI = 0; idI < nData; ++idI)
        {
          printf("%3d %3d\n", idI, dataC[idI]);
        }
      }
      else
      {
        for(idI = 0; idI < nData; ++idI)
        {
          printf("%3d %3d\n", idI, dataI[idI]);
        }
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
                   "Usage: %s [-c] [-i] [-r] [-s#]\n%s",
                   "Options:\n"
                   "  -h  Show this usage message.\n"
                   "  -c  Do a heap sort of chars instead of ints.\n"
                   "  -i  Do an indexed heap sort.\n"
                   "  -r  Use random instead of inverse values\n"
                   "  -s  Data size (number of data elements to sort).\n"
                   "Test for AlgHeapSort() and AlgHeapSortIdx().\n",
                   *argv);
  }
  return(!ok);
}


