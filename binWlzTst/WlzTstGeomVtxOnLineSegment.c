#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstGeomVtxOnLineSegment_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstGeomVtxOnLineSegment.c
* \author       Bill Hill
* \date         June 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Test for WlzGeomVtxOnLineSegment().
* \ingroup	binWlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int		main(int argc, char *argv[])
{
  int		option,
		onSeg = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  double	tol = WLZ_MESH_TOLERANCE;
  WlzDVertex2	seg0,
  		seg1,
		tst;
  static char   optList[] = "he:f:l:t:v";

  seg0.vtX = 0.0,
  seg0.vtY = 0.0;
  seg1.vtX = 1.0;
  seg1.vtY = 1.0;
  tst.vtX = 0.5;
  tst.vtY = 0.5;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'e':
	if(sscanf(optarg, "%lg,%lg", &(seg0.vtX), &(seg0.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'f':
	if(sscanf(optarg, "%lg,%lg", &(seg1.vtX), &(seg1.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'l':
	if(sscanf(optarg, "%lg", &tol) != 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	if(sscanf(optarg, "%lg,%lg", &(tst.vtX), &(tst.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'v':
        verbose = 1;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    usage = optind != argc;
  }
  ok = usage == 0;
  if(ok)
  {
    if(verbose)
    {
      (void )printf("Test using tollerance (%g) that vertex (%g, %g) lies\n"
                    "on line segment (%g, %g), (%g, %g):",
      		    tol, tst.vtX, tst.vtY,
		    seg0.vtX, seg0.vtY, seg1.vtX, seg1.vtY);
    }
    onSeg = WlzGeomVtxOnLineSegment(tst, seg0, seg1, tol);
    printf("%d\n", onSeg);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-e #,#] [-f #,#] [-l #] [-t #,#] [-v]\n"
    "Options are:\n"
    " -e  First vertex of line segment.\n"
    " -f  Second vertex of line segment.\n"
    " -l  Tollerance value.\n"
    " -t  Test vertex.\n"
    " -v  Verbose output.\n"
    "Tests whether the tst vertex lines on the given line segment, using the\n"
    "given tollerance value.\n",
    argv[0]);
  }
  exit(!ok);
}
