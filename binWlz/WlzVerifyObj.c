#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzVerifyObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzVerifyObj.c
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
* \brief	Checks objects and fixes them if possible.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzverifyobj "WlzVerifyObj"
*/


/*!
\ingroup BinWlz
\defgroup wlzverifyobj WlzVerifyObj
\par Name
WlzVerifyObj - checks objects and fixes them if possible.
\par Synopsis
\verbatim
WlzVerifyObj [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Check each input object, fix it if possible and write to stdout.
\par Examples
\verbatim
WlzVerifyObj in.wlz >in_fixed.wlz 
\endverbatim
Reads the object from in.wlz, attempts to fixes any errors and then
writes the fixed object to in_fixed.wlz.
Any error messages are written to the standard error output.
\par File
\ref WlzVerifyObj.c "WlzVerifyObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzVerifyObject "WlzVerifyObject(3)"
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
  (void )fprintf(stderr,
	  "Usage:\t%s [-h] [<input file>]\n"
	  "\tCheck each input object, fix if possible and write to stdout\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        Help - prints this usage message\n",
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
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
    
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

  /* loop reading objects correcting as required */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL )
  {
    if((errNum = WlzVerifyObject(obj, 1)) == WLZ_ERR_NONE) {
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: error when writing verified object (%s).\n",
		       argv[0], errMsg);
      }
    }
    else {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: error when verifying object (%s).\n",
      		     argv[0], errMsg);
    }
    WlzFreeObj(obj);
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
