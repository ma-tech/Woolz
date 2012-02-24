#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVtxDistance_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzVtxDistance.c
* \author       Richard Baldock
* \date         May 2006
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
* \brief        Calculate distance between two vertices.
* \ingroup      BinWlzApp
*               
* \par Binary
* \ref wlzvtxdistance "WlzVtxDistance"
*/
 
/*!
\ingroup      BinWlzApp
\defgroup     wlzvtxdistance WlzVtxDistance
\par Name
WlzVtxDistance - 
\par Synopsis
\verbatim
WlzVtxDistance

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-</b></td>
    <td> </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>

\par Description

\par Examples
\verbatim
\endverbatim

\par See Also
\par Bugs
None known
*/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s x1,y1,z1 x2,y2,z2\n"
	  "\tCalculate the distance between the two vertices\n"
	  "\twriting the distance to standard output\n"
	  "\tOptions are:\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  char 		optList[] = "h";
  int		option;
  double	x1, y1, z1, x2, y2, z2;
  double	dist;
    
  /* additional defaults */
  x1 = y1 = z1 = 0.0;
  x2 = y2 = z2 = 0.0;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    }
  }

  /* read input verticex values */
  if (sscanf(argv[1], "%lf,%lf,%lf", &x1, &y1, &z1) != 3) {
    fprintf(stderr, "%s: not enough input data\n", argv[0]);
    usage(argv[0]);
    return 1;
  }
  if (sscanf(argv[2], "%lf,%lf,%lf", &x2, &y2, &z2) != 3) {
    fprintf(stderr, "%s: not enough input data\n", argv[0]);
    usage(argv[0]);
    return 1;
  }

  dist = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
  dist = sqrt(dist);
  fprintf(stdout, "%f\n", dist);

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
