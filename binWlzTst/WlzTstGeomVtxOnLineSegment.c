#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstGeomVtxOnLineSegment_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstGeomVtxOnLineSegment.c
* \author       Bill Hill
* \date         June 2008
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
* \brief	Test for WlzGeomVtxOnLineSegment[23]D().
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
		dim = 2,
		onSeg = 0,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  double	tol = WLZ_MESH_TOLERANCE;
  WlzVertex	seg0,
  		seg1,
		tmp,
		tst;
  static char   optList[] = "3he:f:l:t:v";

  seg0.d3.vtX = 0.0,
  seg0.d3.vtY = 0.0;
  seg0.d3.vtZ = 0.0;
  seg1.d3.vtX = 1.0;
  seg1.d3.vtY = 1.0;
  seg1.d3.vtZ = 1.0;
  tst.d3.vtX = 0.5;
  tst.d3.vtY = 0.5;
  tst.d3.vtZ = 0.5;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case '3':
        dim = 3;
	break;
      case 'e':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg", &(seg0.d3.vtX), &(seg0.d3.vtY)) != 2)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg", &(seg0.d3.vtX), &(seg0.d3.vtY),
	            &(seg0.d3.vtZ)) != 3)
	  {
	    usage = 1;
	  }
	}
        break;
      case 'f':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg", &(seg1.d3.vtX), &(seg1.d3.vtY)) != 2)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg", &(seg1.d3.vtX), &(seg1.d3.vtY),
	            &(seg1.d3.vtZ)) != 3)
	  {
	    usage = 1;
	  }
	}
        break;
      case 'l':
	if(sscanf(optarg, "%lg", &tol) != 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	if(dim == 2)
	{
	  if(sscanf(optarg, "%lg,%lg", &(tst.d3.vtX), &(tst.d3.vtY)) != 2)
	  {
	    usage = 1;
	  }
	}
	else /* dim == 3 */
	{
	  if(sscanf(optarg, "%lg,%lg,%lg", &(tst.d3.vtX), &(tst.d3.vtY),
	            &(tst.d3.vtZ)) != 3)
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
      tmp.d2.vtX = tst.d3.vtX; tmp.d2.vtY = tst.d3.vtY; tst.d2 = tmp.d2;
      tmp.d2.vtX = seg0.d3.vtX; tmp.d2.vtY = seg0.d3.vtY; seg0.d2 = tmp.d2;
      tmp.d2.vtX = seg1.d3.vtX; tmp.d2.vtY = seg1.d3.vtY; seg1.d2 = tmp.d2;
      if(verbose)
      {
	(void )printf("Test using tollerance (%g) that vertex (%g, %g)\n"
	              "lies on line segment (%g, %g), (%g, %g):",
		      tol, tst.d2.vtX, tst.d2.vtY,
		      seg0.d2.vtX, seg0.d2.vtY,
		      seg1.d2.vtX, seg1.d2.vtY);
      }
      onSeg = WlzGeomVtxOnLineSegment2D(tst.d2, seg0.d2, seg1.d2, tol);
    }
    else /* dim == 3 */
    {
      if(verbose)
      {
	(void )printf("Test using tollerance (%g) that vertex (%g, %g, %g)\n"
	              "lies on line segment (%g, %g, %g), (%g, %g, %g):",
		      tol, tst.d3.vtX, tst.d3.vtY, tst.d3.vtZ,
		      seg0.d3.vtX, seg0.d3.vtY, seg0.d3.vtZ,
		      seg1.d3.vtX, seg1.d3.vtY, seg1.d3.vtZ);
      }
      onSeg = WlzGeomVtxOnLineSegment3D(tst.d3, seg0.d3, seg1.d3, NULL);
    }
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
