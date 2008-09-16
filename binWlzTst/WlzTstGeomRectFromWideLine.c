#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstGeomRectFromWideLine_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstGeomRectFromWideLine.c
* \author       Bill Hill
* \date         September 2008
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
* \brief	Test program for WlzGeomVxInTriangle() and
* 		WlzGeomVxInTetrahedron().
* \ingroup	WlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		gFlg = 0,
  		option,
  		ok = 1,
  		usage = 0;
  double	w;
  WlzDVertex2	s,
  		t;
  WlzDVertex2	rect[4];
  static char   optList[] = "ghs:t:w:";


  w = 1.0;
  s.vtX = 0.0; s.vtY = 0.0;
  t.vtX = 1.0; t.vtY = 1.0;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'g':
        gFlg = 1;
	break;
      case 's':
	if(sscanf(optarg, "%lg,%lg", &(s.vtX), &(s.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 't':
	if(sscanf(optarg, "%lg,%lg", &(t.vtX), &(t.vtY)) != 2)
	{
	  usage = 1;
	}
        break;
      case 'w':
	if(sscanf(optarg, "%lg", &w) != 1)
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
  ok = usage == 0;
  if(ok)
  {
    if(WlzGeomRectFromWideLine(s, t, w, rect) != 0)
    {
      (void )fprintf(stderr, "%s: Vertices are coincident.\n",
      		     argv[0]);
      ok = 0;
    }
    else
    {
      if(gFlg)
      {
	(void )printf("%lg %lg %lg,%lg\n"
	              "%lg %lg %lg %lg\n"
	              "%lg %lg %lg %lg\n"
	              "%lg %lg %lg %lg\n\n",
	              rect[0].vtX, rect[0].vtY, rect[1].vtX, rect[1].vtY,
	              rect[1].vtX, rect[1].vtY, rect[2].vtX, rect[2].vtY,
	              rect[2].vtX, rect[2].vtY, rect[3].vtX, rect[3].vtY,
	              rect[3].vtX, rect[3].vtY, rect[0].vtX, rect[0].vtY);
      }
      else
      {
	(void )printf("%lg,%lg\n"
	              "%lg,%lg\n"
		      "%lg,%lg\n"
		      "%lg,%lg\n",
	              rect[0].vtX, rect[0].vtY,
	              rect[1].vtX, rect[1].vtY,
	              rect[2].vtX, rect[2].vtY,
	              rect[3].vtX, rect[3].vtY);
      }
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-g] [-s <x>,<y>>] [-t <x>,<y>]\n"
    "Options are:\n"
    " -h  Prints this usage summary.\n"
    " -g  Print output as line segments suitable for ploting.\n"
    " -s  Coordinates of vertex S.\n"
    " -t  Coordinates of vertex T\n"
    " -w  Line width (W).\n"
    "Prints the coordinates of the vertices of a wide line between\n"
    "the two points S and T of width W.\n"
    "Default values of S and T are (0.0, 0.0) and (1.0, 1.0).\n"
    "Example:\n"
    "  %s -g -s 3,4 -t 7,9\n",
    argv[0], argv[0]);
  }
  exit(!ok);
}
