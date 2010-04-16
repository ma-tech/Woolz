#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzTstItrSpiral_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzTstItrSpiral.c
* \author       Bill Hill
* \date         July 2007
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
* \brief	Test program for WlzGeomItrSpiral2I() and WlzGeomItrSpiral3I().
* \ingroup	WlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		option,
		nStep = 100,
  		dim = 2,
		step,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  WlzIVertex3	v3;
  WlzIVertex2	v2;
  static char   optList[] = "23hn:";

  v2.vtX = 0;
  v2.vtY = 0;
  v3.vtX = 0;
  v3.vtY = 0;
  v3.vtZ = 0;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'n':
	if(sscanf(optarg, "%d", &nStep) != 1)
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
    (void )printf("# vtk DataFile Version 1.0\n"
                  "WlzTstItrSpiral output\n"
		  "ASCII\n"
		  "DATASET POLYDATA\n"
		  "POINTS %d float\n",
		  nStep);
    step = 0;
    while(step < nStep)
    {
      if(dim == 2)
      {
        step = WlzGeomItrSpiral2I(step, &(v3.vtX), &(v3.vtY));
	(void )printf("%d %d 0\n", v3.vtX, v3.vtY);
      }
      else /* dim == 3 */
      {
        step = WlzGeomItrSpiral3I(step, &(v3.vtX), &(v3.vtY), &(v3.vtZ));
	(void )printf("%d %d %d\n", v3.vtX, v3.vtY, v3.vtZ);
      }
    }
    (void )printf("LINES %d %d\n", nStep - 1, (nStep - 1) * 3);
    for(step = 0; step < nStep - 1; ++step)
    {
      (void )printf("2 %d %d\n", step, step + 1);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-2|3] [-n <number of steps>] [-v]\n"
    "Options are:\n"
    " -2  2D.\n"
    " -3  3D.\n"
    " -n  Number of spiral steps.\n"
    " -v  Verbose output.\n"
    "Draws a sprial in VTK polydata format, with the output going to the\n"
    "standard output.\n",
    argv[0]);
  }
  exit(!ok);
}
