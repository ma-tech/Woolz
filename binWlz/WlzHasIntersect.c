#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzHasIntersect.c
* \author       Richard Baldock
* \date         December 2000
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Tests whether a pair of objects intersect.
* \ingroup	binWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzhasintersect "WlzHasIntersect"
*/

/*!
\ingroup BinWlz
\defgroup wlzhasintersect WlzHasIntersect
\par Name
WlzHasIntersect - tests for the intersection of a pair of objects.
\par Synopsis
\verbatim
WlzHasIntersect [-n] [-h] [-v] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>
    Use numerical output with 0 for no intersect and 1 for intersection.
    </td>
  </tr>
</table>
\par Description
Test if the input objects have a non-zero intersection.
\par Examples
WlzHasIntersect -n dom0.wlz dom1.wlz
\verbatim
Prints 1 if the domains of the objects read from dom0.wlz and dom1.wlz
intersect otherwise prints 0.
\endverbatim
\par File
\ref WlzHasIntersect.c "WlzHasIntersect.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzHasIntersection "WlzHasIntersection(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
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
	  "Usage:\t%s [-n] [-h] [-v] [<input file>]\n"
	  "\tTest if the input objects have a non-zero intersection\n"
	  "\tOptions are:\n"
	  "\t  -n        numerical output: 0 - no intersect, 1 - intersect\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj1, *obj2;
  FILE		*inFile;
  char 		optList[] = "hnv";
  int		option;
  int		verboseFlg=0;
  int		numOutputFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'n':
      numOutputFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* read objects */
  if( obj1 = WlzAssignObject(WlzReadObj(inFile, NULL), NULL) ){
    if( obj2 = WlzAssignObject(WlzReadObj(inFile, NULL), NULL) ){

      if( obj1->type != obj2->type ){
	fprintf(stderr, "%s: objects must be of the same type\n", argv[0]);
	return 1;
      }

      switch( obj1->type ){
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( WlzHasIntersection(obj1, obj2, &errNum) ){
	  if( numOutputFlg ){
	    fprintf(stdout, "1");
	  }
	  else {
	    fprintf(stdout, "Objects intersect\n");
	  }
	}
	else {
	  if( errNum == WLZ_ERR_NONE ){
	    if( numOutputFlg ){
	      fprintf(stdout, "0");
	    }
	    else {
	      fprintf(stdout, "Objects do not intersect\n");
	    }
	  }
	  else {
	    fprintf(stderr, "%s: some sort of error\n", argv[0]);
	  }
	}
	break;

      default:
	fprintf(stderr, "%s: wrong object types\n", argv[0]);
	break;

      }
    }
    else {
      fprintf(stderr, "%s: second object not found\n", argv[0]);
      return 1;
    }
  }
  else {
    fprintf(stderr, "%s: first object not found\n", argv[0]);
    return 1;
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
