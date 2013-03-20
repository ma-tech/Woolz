#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomTetraAffineSolve_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstGeomTetraAffineSolve.c
* \author       Bill Hill
* \date         May 2007
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
* \brief	Test program for libWlz/WlzGeomTetraAffineSolve().
* \ingroup	BinWlzTst
*/

#include <stdio.h>
#include <float.h>
#include <sys/time.h>
#include <Wlz.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static double	MyRandom(void)
{
  double	rnd;

  rnd = AlgRandNormal(0.0, 1.0);
  return(rnd);
}

static void	SetVertices(WlzDVertex3 *sVx, WlzDVertex3 *dVx)
{
 int		idN;

 for(idN = 0; idN < 4; ++idN)
 {
   sVx[idN].vtX = MyRandom();
   sVx[idN].vtY = MyRandom();
   sVx[idN].vtZ = MyRandom();
   dVx[idN].vtX = (sVx[idN].vtX * MyRandom()) + MyRandom();
   dVx[idN].vtY = (sVx[idN].vtY * MyRandom()) + MyRandom();
   dVx[idN].vtZ = (sVx[idN].vtZ * MyRandom()) + MyRandom();
  }
  if(MyRandom() < 0.0)
  {
    if(MyRandom() < -0.3)
    {
      dVx[3].vtX = dVx[2].vtX = dVx[1].vtX = dVx[0].vtX;
    }
    if(MyRandom() < -0.3)
    {
      dVx[3].vtY = dVx[2].vtY = dVx[1].vtY = dVx[0].vtY;
    }
    if(MyRandom() < -0.3)
    {
      dVx[3].vtZ = dVx[2].vtZ = dVx[1].vtZ = dVx[0].vtZ;
    }
  }
}

int		main(int argc, char *argv[])
{
  int		idN,
		idR,
		option,
		squashed,
		timer = 0,
		useLU = 0,
		silent = 0,
  		repeats = 1,
		usage = 0;
  double	ss,
  		delta = 0.000001;
  WlzAffineTransform *tr = NULL;
  struct timeval times[3];
  WlzDVertex3	eVx;
  WlzDVertex3	sVx[4],
		tVx[4],
  		dVx[4];
  static char	optList[] = "hlstd:r:";

  opterr = 0;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	if(sscanf(optarg, "%lg", &delta) != 1)
	{
	  usage = 1;
	}
	break;
      case 'l':
	useLU = 1;
	break;
      case 't':
	timer = 1;
	break;
      case 'r':
	if(sscanf(optarg, "%d", &repeats) != 1)
	{
	  usage = 1;
	}
	break;
      case 's':
	silent = 1;
	break;
      case 'h':  /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    (void )AlgRandSeed(0);
    tr = WlzMakeAffineTransform(WLZ_TRANSFORM_3D_AFFINE, NULL);
    SetVertices(sVx, dVx);
    if(timer)
    {
      gettimeofday(&(times[0]), NULL);
    }
    for(idR = 0; idR < repeats; ++idR)
    {
      ss = 0.0;
      if(useLU)
      {
	squashed = WlzGeomTetraAffineSolveLU(*(tr->mat), sVx, dVx);
      }
      else
      {
	squashed = WlzGeomTetraAffineSolve(*(tr->mat), sVx, dVx, delta);
      }
      if(!silent)
      {
	for(idN = 0; idN < 4; ++idN)
	{
	  tVx[idN] = WlzAffineTransformVertexD3(tr, sVx[idN], NULL);
	  WLZ_VTX_3_SUB(eVx, tVx[idN], dVx[idN]);
	  ss += WLZ_VTX_3_SQRLEN(eVx);
	}
	if(ss > DBL_EPSILON)
	{
	  ss = sqrt(ss);
	}
	(void)
	printf("%d %+08e | %+08e %+08e %+08e | %+08e %+08e %+08e | "
	       "%+08e %+08e %+08e|\n",
	       squashed, ss,
	       sVx[0].vtX, sVx[0].vtY, sVx[0].vtZ, 
	       dVx[0].vtX, dVx[0].vtY, dVx[0].vtZ,
	       tVx[0].vtX - dVx[0].vtX,
	       tVx[0].vtY - dVx[0].vtY,
	       tVx[0].vtZ - dVx[0].vtZ);
	for(idN = 1; idN < 4; ++idN)
	{
	  (void)
	  printf("                | %+08e %+08e %+08e | %+08e %+08e %+08e | "
		 "%+08e %+08e %+08e|\n",
		 sVx[idN].vtX, sVx[idN].vtY, sVx[idN].vtZ, 
		 dVx[idN].vtX, dVx[idN].vtY, dVx[idN].vtZ,
		 tVx[idN].vtX - dVx[idN].vtX,
		 tVx[idN].vtY - dVx[idN].vtY,
		 tVx[idN].vtZ - dVx[idN].vtZ);
	}
	for(idN = 0; idN < 4; ++idN)
	{
	  (void )
	  printf("%+08e %+08e %+08e %+08e\n",
		 tr->mat[idN][0], tr->mat[idN][1],
		 tr->mat[idN][2], tr->mat[idN][3]);
	}
      }
      if(!silent)
      {
	SetVertices(sVx, dVx);
      }
    }
    if(timer)
    {
      double s = 0.0;

      gettimeofday(&(times[1]), NULL);
      timersub(&(times[1]), &(times[0]), &(times[2]));
      s = times[2].tv_sec + (0.000001 * times[2].tv_usec);
      (void )
      printf("Total elapsed time = %gs\n", s);
    }
    (void )WlzFreeAffineTransform(tr);
  }
  else
  {
    (void )fprintf(stderr,
      "Usage: %s [-h] [-l] [-s] [-t] [-r #] [-d #]\n"
      "Runs tests on functions WlzGeomTetraAffineSolve() and\n"
      "WlzGeomTetraAffineSolveLU().\n"
      "Version: %s\n"
      "Options are:\n"
      "  -h  Help, prints this usage message.\n"
      "  -l  Use WlzGeomTetraAffineSolveLU() instead of\n"
      "      WlzGeomTetraAffineSolve().\n"
      "  -s  Silent - not matrices etc output.\n"
      "  -t  Output timings.\n"
      "  -r  Number of repeat calls.\n"
      "  -d  Threshold volume (x 6) for WlzGeomTetraAffineSolve().\n",
      argv[0],
      WlzVersion());
  }
  return(0);
}
