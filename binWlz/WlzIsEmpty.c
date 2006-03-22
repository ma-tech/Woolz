#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzIsEmpty_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzIsEmpty.c
* \author       Richard Baldock
* \date         November 2000
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
* \brief	Determines whether objects are empty.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzisempty "WlzIsEmpty"
*/

/*!
\ingroup BinWlz
\defgroup wlzisempty WlzIsEmpty
\par Name
WlzIsEmpty - determines whether objects are empty.
\par Synopsis
\verbatim
WlzIsEmpty [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Determines whether objects are empty.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzIsEmpty.c "WlzIsEmpty.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
\ref WlzIsEmpty "WlzIsEmpty(3)"
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
	  "\tDetermines whether input object(s) are empty.\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  int		count, vol;
  const char    *errMsg;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
    
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

    if( WlzIsEmpty(obj, &errNum) ){
      fprintf(stderr, "Object %d: empty\n", count);
    }
    else {
      if( errNum != WLZ_ERR_NONE ){
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	fprintf(stderr,
		"%s: Object %d: error in determining if empty (%s).\n",
		*argv, count, errMsg);
      }
      else {
	fprintf(stderr, "Object %d: not empty\n", count);
      }
    }

    WlzFreeObj( obj );
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
