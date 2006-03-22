#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzLabel_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzLabel.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Labels (segments) the input objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzlabel "WlzLabel"
*/

/*!
\ingroup BinWlz
\defgroup wlzlabel WlzLabel
\par Name
WlzLabel - labels (segments) the input objects.
\par Synopsis
\verbatim
WlzLabel [i#] [-v] [-h] [<input file>]
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
    <td><b>-i</b></td>
    <td>Ignore objects with number of lines \f$<\f$ the given number.</td>
  </tr>
</table>
\par Description
Label (segment) the input objects and write the result
to stdout. Non-domain objects are ignored, the number
of segments found is written to stderr.
\par Examples
\verbatim
WlzThreshold -v 200 -L sec.wlz | WlzLabel | WlzExpolde -u -b sec -
\endverbatim
Thresholds the object read from the file sec.wlz
keeping all values less than 200,
labels the resulting domain and then expoldes the concatonated
objects into seperate files.

\verbatim
WlzThreshold -v 150 -L sec.wlz | WlzLabel | WlzArea
\endverbatim
Thresholds the object read from the file sec.wlz
keeping all values less than 200,
labels the resulting domain and then prints the area of each of the
resulting objects.
\par File
\ref WlzLabel.c "WlzLabel.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzthreshold "WlzThreshold(1)"
\ref wlzfacts "WlzFacts(1)"
\ref wlzarea "WlzArea(1)"
\ref WlzLabel "WlzLabel(3)"
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
	  "Usage:\t%s [i#] [-v] [-h] [<input file>]\n"
	  "\tLabel (segment) the input objects and write the result\n"
	  "\tto stdout. Non-domain objects are ignored, the number\n"
	  "\tof segments found is written to stderr\n"
	  "\tOptions are:\n"
	  "\t  -i#       Ignore objects with number of lines < #\n"
	  "\t  -v        Verbose flag\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str);
  return;
}
 
#define MAXOBJS 2000    /* Was 10000 */

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  WlzObject	**objlist = NULL;
  FILE		*inFile;
  char 		optList[] = "i:vh";
  int		option;
  int		count, numobj, i, verbose = 0;
  int		ignw = -1;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'i':
      ignw = atoi(optarg);
      if( ignw < 0 ){
        fprintf(stderr, "%s: ignw = %d is invalid\n", argv[0], ignw);
        usage(argv[0]);
        return( 1 );
      }
      break;
 
    case 'v':
      verbose = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return( 1 );

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "rb")) == NULL ){
      (void )fprintf(stderr,
      		     "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return( 1 );
    }
  }

  /* read objects and segment if possible */
  ignw = (ignw < 0) ? 1 : ignw;
  count = 0;
  while((obj = WlzReadObj(inFile, NULL)) != NULL) {
    count++;
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      errNum = WlzLabel(obj, &numobj, &objlist, MAXOBJS, ignw, WLZ_8_CONNECTED);
      if(errNum == WLZ_ERR_DOMAIN_TYPE) {
	errNum = WlzWriteObj(stdout, obj);
      }
      else {
	if(verbose) {
	  fprintf(stderr,"%s: writing %d objects from input object %d\n",
		  argv[0], numobj, count);
	}
	for(i=0; (i < numobj) && (errNum == WLZ_ERR_NONE); i++){
	  errNum = WlzWriteObj(stdout, *(objlist + i));
	  WlzFreeObj(*(objlist + i));
	}
      }
      break;

    default:
      errNum = WlzWriteObj(stdout, obj);
      break;

    }

    WlzFreeObj(obj);

  }
  if(errNum != WLZ_ERR_NONE) {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,"%s: failed to label object (%s).\n",
    		   argv[0],
		   errMsg);
    return( 1 );
  }

  return( 0 );
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
