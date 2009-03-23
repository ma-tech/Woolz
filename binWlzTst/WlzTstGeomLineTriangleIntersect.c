#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstGeomLineTriangleIntersect_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstGeomLineTriangleIntersect.c
* \author       Bill Hill
* \date         March 2009
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2009 Medical research Council, UK.
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
* \brief	Command line test program for WlzGeomLineTriangleIntersect3D.
* \ingroup	binWlzTst
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
		dim,
		isn = 0,
		par = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  double	t,
  		u,
		v;
  WlzVertex	org,
  		dir;
  WlzVertex	tri[3];
  static char   optList[] = "23hd:g:t:v";

  dim = 3;
  WLZ_VTX_3_SET(org.d3,  0.0,  4.0,  0.0);
  WLZ_VTX_3_SET(dir.d3, -1.0,  2.0,  0.0);
  WLZ_VTX_3_SET(tri[0].d3,  2.0,  0.0,  0.0);
  WLZ_VTX_3_SET(tri[1].d3,  3.0,  1.0,  1.0);
  WLZ_VTX_3_SET(tri[2].d3,  3.0,  1.0,  0.0);
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case '2':
        dim = 3;
	break;
      case '3':
        dim = 3;
	break;
      case 'd':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg", &(dir.d2.vtX), &(dir.d2.vtY)) != 2)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg", &(dir.d3.vtX), &(dir.d3.vtY),
	            &(dir.d3.vtZ)) != 3)
	  {
	    usage = 1;
	  }
	}
        break;
      case 'g':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg", &(org.d2.vtX), &(org.d2.vtY)) != 2)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg", &(org.d3.vtX), &(org.d3.vtY),
	            &(org.d3.vtZ)) != 3)
	  {
	    usage = 1;
	  }
	}
        break;
      case 't':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg,%lg,%lg,%lg,%lg",
	            &(tri[0].d2.vtX), &(tri[0].d2.vtY),
	            &(tri[1].d2.vtX), &(tri[1].d2.vtY),
	            &(tri[2].d2.vtX), &(tri[2].d2.vtY)) != 6)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
	            &(tri[0].d3.vtX), &(tri[0].d3.vtY), &(tri[0].d3.vtZ),
	            &(tri[1].d3.vtX), &(tri[1].d3.vtY), &(tri[1].d3.vtZ),
	            &(tri[2].d3.vtX), &(tri[2].d3.vtY), &(tri[2].d3.vtZ)) != 9)
	  {
	    usage = 1;
	  }
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
    if(dim == 2)
    {
      (void )fprintf(stderr, "%s: 2D code not implimented yet.\n", argv[0]);
      ok = 0;
    }
    else /* dim == 3 */
    {
      if(verbose)
      {
	(void )printf(
	"Test line with origin (%g, %g, %g) and direction (%g, %g, %g)\n"
	"passes through the triangle (%g, %g, %g), (%g, %g, %g),\n"
	"(%g, %g, %g):",
	org.d3.vtX, org.d3.vtY, org.d3.vtZ,
	dir.d3.vtX, dir.d3.vtY, dir.d3.vtZ,
	tri[0].d3.vtX, tri[0].d3.vtY, tri[0].d3.vtZ,
	tri[1].d3.vtX, tri[1].d3.vtY, tri[1].d3.vtZ,
	tri[2].d3.vtX, tri[2].d3.vtY, tri[2].d3.vtZ);
      }
      isn = WlzGeomLineTriangleIntersect3D(org.d3, dir.d3,
      					   tri[0].d3, tri[1].d3, tri[2].d3,
					   &par, &t, &u, &v);
    }
  }
  if(ok)
  {
    printf("%d %d %lg %lg %lg\n", isn, par, t, u, v);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-e #,#] [-f #,#] [-l #] [-t #,#] [-v]\n"
    "Options are:\n"
    " -d  Line direction.\n"
    " -g  Line origin.\n"
    " -t  Triangle vertex coordinate list.\n"
    " -v  Verbose output.\n"
    "Tests whether a line passes through a triangle.\n",
    argv[0]);
  }
  exit(!ok);
}
