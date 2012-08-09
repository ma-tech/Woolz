#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzArea_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzArea.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Calculates the area of 2D domain objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzarea "WlzArea"
*/

/*!
\ingroup BinWlz
\defgroup wlzarea WlzArea
\par Name
WlzArea - calculates the area of 2D domain objects.
\par Synopsis
\verbatim
WlzArea [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Calculates the area of 2D domain objects.
\par Examples
\verbatim
WlzArea fish4.wlz
\endverbatim
Computes the area of the 2D domain objects fish4.wlz.
\par File
\ref WlzArea.c "WlzArea.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzvolume "WlzVolume(1)"
\ref WlzArea "WlzArea(3)"
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
	  "Usage:\t%s [-h] [<input file>]\n"
	  "Version: %s"
	  "\n\tCalculate the area of the input 2D woolz objects\n"
	  "\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  int		count, area;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case '~': /* dummy to avoid compiler warning */
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  count = 0;
  while( (obj = WlzReadObj(inFile, NULL)) != NULL )
  {
    count++;

    switch( obj->type )
    {
     case WLZ_2D_DOMAINOBJ:
       if( (area = WlzArea( obj, NULL )) < 0 ){
	 fprintf(stderr, "Object %d: error in calculating the volume\n",
		 count);
       }
       else {
	 fprintf(stderr, "Object %d: number of pixels = %d\n", count, area);
       }
       break;

     default:
       fprintf(stderr, "Object %d: not 2D object type\n", count);
       break;

    }

    WlzFreeObj( obj );
  }

  return( WLZ_ERR_NONE );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
