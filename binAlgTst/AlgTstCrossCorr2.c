#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstCrossCorr2_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstCrossCorr2.c
* \author       Bill Hill
* \date         November 2007
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2007 Medical research Council, UK.
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
* \brief	Simple test for the cross correlation code in libAlg.
* \ingroup	binAlgTst
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <Alc.h>
#include <Alg.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           iX,
                iY,
                oX,
                oY,
                nX = 16,
                nY = 16,
                rep,
                nRep = 1024,
                allRep = 0;
  double        sum;
  double        **data0 = NULL,
                **data1 = NULL;
  AlgError      errNum = ALG_ERR_NONE;

  AlgRandSeed(time(NULL));
  if((AlcDouble2Malloc(&data0, nY, nX) != ALC_ER_NONE) ||
     (AlcDouble2Malloc(&data1, nY, nX) != ALC_ER_NONE))
  {
    errNum = ALG_ERR_MALLOC;
  }
  else
  {
    for(rep = 0; rep < nRep; ++rep)
    {
      /* iX = (int )(8.0 * (0.5 - AlgRandUniform())); */
      /* iY = (int )(8.0 * (0.5 - AlgRandUniform())); */
      iX = 3;
      iY = -1;
      for(oY = 0; oY < nY; ++oY)
      {
        for(oX = 0; oX < nX; ++oX)
        {
          *(*(data0 + oY) + oX) = 0.0;
          *(*(data1 + oY) + oX) = 0.0;
        }
      }
      *(*(data0 + 8 + iY) + 8 + iX) = 1.0;
      *(*(data1 + 8) + 8) = 1.0;
      errNum = AlgCrossCorrelate2D(data0, data1, nX, nY);
      if(errNum == ALG_ERR_NONE)
      {
        sum = 0;
        for(oY = 0; oY < nY; ++oY)
        {
          for(oX = 0; oX < nX; ++oX)
          {
            sum += *(*(data0 + oY) + oX);
          }
        }
        if(sum < 1.0)
        {
          sum = 1.0;
        }
        AlgCrossCorrPeakXY(&oX, &oY, NULL, data0, nX, nY, 5, 5);
        if(allRep || (iX != oX) || (iY != oY))
        {
          (void )fprintf(stderr, "(%d %d) (%d %d)\n", iX, iY, oX, oY);
          for(oY = 0; oY < nY; ++oY)
          {
            for(oX = 0; oX < nX; ++oX)
            {
              if(*(*(data0 + oY) + oX) > sum / 2.0)
              {
                (void )fprintf(stderr, "#");
              }
              else
              {
                (void )fprintf(stderr, "+");
              }
            }
            (void )fprintf(stderr, "\n");
          }
          (void )fprintf(stderr, "\n");
        }
      }
    }
  }
  return(0);
}
