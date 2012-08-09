#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSetBackground_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSetBackground.c
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
* \brief	Sets the background value of a domain object with
* 		grey values.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzsetbackground "WlzSetBackground"
*/

/*!
\ingroup BinWlz
\defgroup wlzsetbackground WlzSetBackground
\par Name
WlzSetBackground - sets the background value of a domain object with
                   grey values.
\par Synopsis
\verbatim
WlzSetBackground [-h] [-v] [-b#] [<input file>]
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
    <td><b>-b</b></td>
    <td>Background value, default 0.</td>
  </tr>
</table>
\par Description
Set the background value of a domain object with grey values to a
new value. The background value is converted to the grey type of
the input object.
The input object is read from the standard input unless a file is given on the
command line. The output object is always written to the standard output.
\par Examples
\verbatim
WlzSetBackground -b 255 in.wlz >out.wlz
\endverbatim
Reads an object from in.wlz,
set the objects background values to 255
and then writes the resulting object to out.wlz.
\par File
\ref WlzSetBackground.c "WlzSetBackground.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzfacts "WlzFacts(1)"
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
	  "Usage:\t%s [-b#] [-h] [-v] [<input file>]\n"
	  "\tSet the background of a grey-level woolz object\n"
	  "\twriting the new object to standard output.\n"
	  "\tThe background value is converted to the grey type\n"
	  "\tof the input object\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -b#       required background value, default 0\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "b:hv";
  int		option;
  WlzPixelV	bckgrnd;
  int		verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  bckgrnd.type = WLZ_GREY_DOUBLE;
  bckgrnd.v.dbv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      bckgrnd.v.dbv = atof(optarg);
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

  /* read objects and threshold if possible */
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      errNum = WlzSetBackground(obj, bckgrnd);
      if( errNum == WLZ_ERR_NONE ){
	if( verboseFlg ){
	  fprintf(stderr,
		  "%s:background set to %f\n",
		  argv[0], bckgrnd.v.dbv);
	}

	if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to write object (%s).\n",
			 argv[0], errMsg);
	  return(1);
	}
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(1);
      }
      break;
    }
    WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
