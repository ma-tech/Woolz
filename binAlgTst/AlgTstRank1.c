#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstRank1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstRank1.c
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
* \brief	Simple test for AlgRank().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

int             AlgRankCmpI(void *v0, void *v1)
{
  int           i0,
                i1;

  i0 = *(int *)v0;
  i1 = *(int *)v1;
  return(i0 - i1);
}


extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char *argv[])
{
  int           id0,
                tI0,
                option,
                seed = 0,
                nElm = 10,
                rank = -1,
                okFlg = 1,
                usageFlg = 0;
  int           *datI = NULL,
                *datV = NULL;
  unsigned char *datUB = NULL;
  short         *datS = NULL;
  float         *datF = NULL;
  double        *datD = NULL;
  static char   optList[] = "hn:r:s:";

  opterr = 0;
  while(okFlg && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
        if((sscanf(optarg, "%d", &nElm) != 1) || (nElm <= 0))
        {
          okFlg = 0;
          usageFlg = 1;
        }
        break;
      case 'r':
        if((sscanf(optarg, "%d", &rank) != 1) || (rank < 0))
        {
          okFlg = 0;
          usageFlg = 1;
        }
        break;
      case 's':
        if(sscanf(optarg, "%d", &seed) != 1)
        {
          okFlg = 0;
          usageFlg = 1;
        }
        break;
      case 'h':
        okFlg = 0;
        usageFlg = 1;
        break;
    }
  }
  if(okFlg)
  {
    if(rank < 0)
    {
      rank = nElm / 2;
    }
    if(((datI = (int *)AlcMalloc(nElm * sizeof(int))) == NULL) ||
       ((datS = (short *)AlcMalloc(nElm * sizeof(short))) == NULL) ||
       ((datUB = (unsigned char *)AlcMalloc(nElm *
                                            sizeof(unsigned char))) == NULL) ||
       ((datF = (float *)AlcMalloc(nElm * sizeof(float))) == NULL) ||
       ((datD = (double *)AlcMalloc(nElm * sizeof(double))) == NULL) ||
       ((datV = (int *)AlcMalloc(nElm * sizeof(int))) == NULL))
    {
      (void )fprintf(stderr, "%s: Failed to allocate memory.\n", argv[0]);
    }
  }
  if(okFlg)
  {
    srandom(seed);
    for(id0 = 0; id0 < nElm; ++id0)
    {
      tI0 = random() % 256;
      *(datI + id0) = tI0;
      *(datS + id0) = tI0;
      *(datUB + id0) = tI0;
      *(datF + id0) = tI0;
      *(datD + id0) = tI0;
      *(datV + id0) = tI0;
    }
    AlgRankSelectI(datI, nElm, rank);
    AlgRankSelectS(datS, nElm, rank);
    AlgRankSelectUB(datUB, nElm, rank);
    AlgRankSelectF(datF, nElm, rank);
    AlgRankSelectD(datD, nElm, rank);
    AlgRankSelectV(datV, nElm, sizeof(int), rank, &tI0, AlgRankCmpI);
    (void )printf("idx    I   S   UB  F   D   V\n");
    for(id0 = 0; id0 < nElm; ++id0)
    {
      (void )printf("%-6d %-3d %-3d %-3d %-3g %-3g %-3d\n",
                    id0, *(datI + id0), *(datS + id0), *(datUB + id0),
                    *(datF + id0), *(datD + id0), *(datV + id0));
    }
    (void )printf("\n       %-3d %-3d %-3d %-3g %-3g %-3d\n",
                  *(datI + rank), *(datS + rank), *(datUB + rank),
                  *(datF + rank), *(datD + rank), *(datV + rank));
  }
  AlcFree(datI);
  AlcFree(datS);
  AlcFree(datUB);
  AlcFree(datF);
  AlcFree(datD);
  AlcFree(datV);
  if(usageFlg)
  {
    fprintf(stderr, "Usage: %s [-h] [-n#] [-r#] [-s#]\n%s",
            argv[0],
            "Test for rank select, selects the item of rank r from a\n"
            "list of random numbers.\n"
            "Options are:\n"
            "  -h  Usage, prints this information.\n"
            "  -n  The number of random elements to select from.\n"
            "  -r  The rank to be selected, 0 smallest and n - 1 largest.\n"
            "      The default is the median.\n"
            "  -s  The random number seed.\n");
  }
  return(!okFlg);
}

