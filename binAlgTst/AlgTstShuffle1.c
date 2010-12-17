#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstShuffle1_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstShuffle1.c
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
* \brief	Simple test for AlgShuffleIdx().
* \ingroup	binAlgTst
*/

#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char *argv[])
{
  int           idx,
                nData = 10,
                seed = 0,
                option,
                ok = 1,
                usage = 0;
  int           *data = NULL;
  static char   optList[] = "n:s:";

  opterr = 0;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
        if((sscanf(optarg, "%d", &nData) != 1) || (nData < 1))
        {
          usage = 1;
          ok = 0;
        }
        break;
      case 's':
        if(sscanf(optarg, "%d", &seed) != 1)
        {
          usage = 1;
          ok = 0;
        }
        break;
      default:
        usage = 1;
        ok = 0;
        break;
    }
  }
  if(ok)
  {
    if((data = (int *)AlcMalloc(nData * sizeof(int))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to allocate storage.\n",
                     *argv);
    }
  }
  if(ok)
  {
    (void )AlgShuffleIdx(nData, data, seed);
    for(idx = 0; idx < nData; ++idx)
    {
      (void )printf("%d\n", *(data + idx));
    }
  }
  if(data)
  {
    AlcFree(data);
  }
  if(usage)
  {
    (void )fprintf(stderr,
                   "Usage: %s [-n#] [-s#]\n",
                   "Options:\n"
                   "  -n  Number of indicies.\n"
                   "  -s  Permutation seed.\n"
                   "Test for AlgShuffleIdx().\n");
  }
  return(!ok);
}

