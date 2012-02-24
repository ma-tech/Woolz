#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstShuffle2_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstShuffle2.c
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
* \brief	Simple test for AlgShuffleIdx().
* \ingroup	binAlgTst
*/

#include <string.h>
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

#define MAX_BUF_LEN 4096
#define DATA_CHUNK_SZ 4096

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char *argv[])
{
  int           idx,
                maxData = 0,
                nData = 0,
                seed = 0,
                option,
                ok = 1,
                usage = 0;
  int           *iData = NULL;
  char          **sData = NULL;
  FILE          *fP;
  char          *inFileStr;
  char          rec[MAX_BUF_LEN];
  static char   optList[] = "hs:",
                defFileStr[]= "-";

  opterr = 0;
  inFileStr = defFileStr;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 's':
        if(sscanf(optarg, "%d", &seed) != 1)
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
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((fP = (strcmp(inFileStr, "-")?
             fopen(inFileStr, "r"): stdin)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: failed to open %s.\n", *argv, inFileStr);

    }
  }
  idx = 0;
  while(ok && (fgets(rec, MAX_BUF_LEN, fP) != NULL))
  {
    *(rec + MAX_BUF_LEN - 1) = '\0';
    *(rec + strlen(rec) - 1) = '\0';
    if(maxData <= idx)
    {
      maxData = (maxData)? maxData * 2: DATA_CHUNK_SZ;
      if((sData = (char **)AlcRealloc(sData,
                                      maxData * sizeof(char *))) == NULL)
      {
        ok = 0;
        (void )fprintf(stderr,
                       "%s: failed to allocate memory.\n", *argv);
      }
    }
    if(ok)
    {
      if((*(sData + idx++) = AlcStrDup(rec)) == NULL)
      {
        ok = 0;
        (void )fprintf(stderr,
                       "%s: failed to allocate memory.\n", *argv);
      }
    }
  }
  if(ok)
  {
    nData = idx;
    if((iData = (int *)AlcMalloc(nData * sizeof(int))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to allocate memory.\n", *argv);
    }
  }
  if(ok)
  {
    (void )AlgShuffleIdx(nData, iData, seed);
    for(idx = 0; idx < nData; ++idx)
    {
      (void )printf("%s\n", *(sData + *(iData + idx)));
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
                   "Usage: %s [-h] [-s#]\n"
                   "Randomize or shuffle the input data records.\n"
                   "Options:\n"
                   "  -h  Print this usage message.\n"
                   "  -s  Random number generator seed.\n",
                   *argv);
  }
  return(!ok);
}

