#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomLineTriangleIntersect_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstGeomLineTriangleIntersect.c
* \author       Bill Hill
* \date         March 2009
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
* \brief	Command line test program for WlzGeomLineTriangleIntersect3D.
* \ingroup	BinWlzTst
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
	"Test for line with origin (%g, %g, %g) and direction (%g, %g, %g)\n"
	"passing through triangle (%g, %g, %g), (%g, %g, %g), (%g, %g, %g):\n",
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
    if(verbose)
    {
      switch(isn)
      {
        case 0:
	  (void )printf("No intersection.\n");
	  break;
	case 1:
	  (void )printf("Intersection on edge of triangle.\n");
	  break;
	case 2:
	  (void )printf("Intersection in triangle.\n");
	  break;
	default:
	  (void )printf("Unknown intersection code %d!\n", isn);
	  break;
      }
      (void )printf("Line is%s parralel to the plane of the triangle.\n",
                    (par)? "": " not");
      if(isn > 1)
      {
	(void )printf("Distance from line origin to intersection is %g.\n",
		      t);
	(void )printf(
	       "Barycentric coordinates of the intersection are %g,%g.\n",
		      u, v);
      }
      (void )printf("\n");
    }
    else
    {
      printf("%d %d %lg %lg %lg\n", isn, par, t, u, v);
    }
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
