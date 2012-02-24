#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstAlgCrossCorr_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstAlgCrossCorr.c
* \author       Bill Hill
* \date         November 2007
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
* Boston, MA  02110-1301, USA.
* \brief	Simple test for the cross correlation code in libAlg.
* \ingroup	BinWlzTst
*/
#include <stdio.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		iX,
		iY,
  		oX,
  		oY,
		rep,
		option,
		ok = 1,
		usage = 0,
		nRep = 16,
             	dataSz = 256,
  		dataSz2,
		dataSz3,
		dataSz4;
  long		seed;
  double	sum;
  double	**data0 = NULL,
  		**data1 = NULL;
  AlgError	errNum = ALG_ERR_NONE;
  static char	optList[] = "hn:s:";

  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'n':
        if(sscanf(optarg, "%d", &nRep) != 1)
	{
	  usage = 1;
	}
	break;
      case 's':
        if(sscanf(optarg, "%ld", &seed) != 1)
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
    dataSz2 = (dataSz - 1) / 2;
    dataSz3 = (dataSz - 1) / 3;
    dataSz4 = (dataSz - 1) / 4;
    AlgRandSeed(seed);
    if((AlcDouble2Malloc(&data0, dataSz, dataSz) != ALC_ER_NONE) ||
       (AlcDouble2Malloc(&data1, dataSz, dataSz) != ALC_ER_NONE))
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to allocate data arrays.\n",* argv);
    }
  }
  if(ok)
  {
    for(rep = 0; rep < nRep; ++rep)
    {
      iX = (int )(dataSz4 * (0.5 - AlgRandUniform()));
      iY = (int )(dataSz4 * (0.5 - AlgRandUniform()));
      for(oY = 0; oY < dataSz; ++oY)
      {
	for(oX = 0; oX < dataSz; ++oX)
	{
	  *(*(data0 + oY) + oX) = 0.0;
	  *(*(data1 + oY) + oX) = 0.0;
        }
      }
      *(*(data0 + dataSz2 + iY) + dataSz2 + iX) = 1.0;
      *(*(data1 + dataSz2) + dataSz2) = 1.0;
      errNum = AlgCrossCorrelate2D(data0, data1, dataSz, dataSz);
      if(errNum == ALG_ERR_NONE)
      {
	sum = 0;
	for(oY = 0; oY < dataSz; ++oY)
	{
	  for(oX = 0; oX < dataSz; ++oX)
	  {
	    sum += *(*(data0 + oY) + oX);
	  }
	}
	if(sum < 1.0)
	{
	  sum = 1.0;
	}
	AlgCrossCorrPeakXY(&oX, &oY, NULL, data0, dataSz, dataSz,
			   dataSz3, dataSz3);
	(void )fprintf(stderr, "%d %d %d %d", iX, iY, oX, oY);
	if((iX != oX) || (iY != oY))
	{
	  (void )fprintf(stderr, " 0\n");
	}
	else
	{
	  (void )fprintf(stderr, " 1\n");
	}
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-h] [-n #] [-s #]\n"
    "Options:\n"
    "  -h   Prints this usage information.\n"
    "  -n  Number of repeats.\n"
    "  -s  Seed for random number generator.\n"
    "Tests AlgCrossCorrelate2D() and AlgCrossCorrPeakXY() by creating\n"
    "arrays with a single non-zero value, in one of which the non-zero\n"
    "value is offset. The output consists of the offset, the computed\n"
    "offset and finaly 1 if the two are equal or 0 if they are not.\n");
  }
  return(!ok);
}
