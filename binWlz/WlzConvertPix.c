#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzConvertPix_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzConvertPix.c
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
* \brief	Converts the pixel type of a grey-level woolz object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzconvertpix "WlzConvertPix"
*/

/*!
\ingroup BinWlz
\defgroup wlzconvertpix WlzConvertPix
\par Name
WlzConvertPix - converts the pixel type of a grey-level woolz object.
\par Synopsis
\verbatim
WlzConvertPix [-t#] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td>
    <td>Output pixel type specified by an integer:
    <table width="500" border="0">
    <tr> <td>1</td> <td>integer</td> </tr>
    <tr> <td>2</td> <td>short</td> </tr>
    <tr> <td>3</td> <td>unsigned byte (default)</td> </tr>
    <tr> <td>4</td> <td>float</td> </tr>
    <tr> <td>5</td> <td>double</td> </tr>
    <tr> <td>7</td> <td>RGBA (colour)</td> </tr>
    </table>
    </tr>
    </td>
  </tr>
</table>
\par Description
Converts the pixel type of a grey-level woolz object writing the new
object to standard output.
\par Examples
\verbatim
WlzConvertPix lobster.wlz >lobster_ub.wlz
\endverbatim
Converts the values of the
domain object (with values) read from the file lobster.wlz
to unsigned byte and writes the resulting object to the file lobster_ub.wlz.
\par File
\ref WlzConvertPix.c "WlzConvertPix.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzConvertPix "WlzConvertPix(3)"
\ref WlzGreyType "WlzGreyType(3)"
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
	  "Usage:\t%s [-t#] [-h] [<input file>]\n"
	  "\tConvert the pixel type of a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -t#       Output pixel type:\n"
	  "\t            # = %d: integer\n"
	  "\t                %d: short\n"
	  "\t                %d: unsigned byte (default)\n"
	  "\t                %d: float\n"
	  "\t                %d: double\n"
	  "\t                %d: rgba\n"
	  "\t  -h        Help - prints this usage message\n"
	  "",
	  proc_str, WLZ_GREY_INT, WLZ_GREY_SHORT, WLZ_GREY_UBYTE,
	  WLZ_GREY_FLOAT, WLZ_GREY_DOUBLE, WLZ_GREY_RGBA);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *newobj;
  FILE		*inFile;
  char 		optList[] = "ht:";
  int		option;
  WlzGreyType	newpixtype = WLZ_GREY_UBYTE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 't':
      switch( newpixtype = (WlzGreyType) atoi(optarg) ){

      case WLZ_GREY_INT:
      case WLZ_GREY_SHORT:
      case WLZ_GREY_UBYTE:
      case WLZ_GREY_FLOAT:
      case WLZ_GREY_DOUBLE:
      case WLZ_GREY_RGBA:
	break;

      default:
        fprintf(stderr, "%s: grey type = %d is invalid\n",
		argv[0],newpixtype );
        usage(argv[0]);
        return 1;

      }
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
      return 1;
    }
  }

  /* read objects and convert if possible */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    if( obj->type == WLZ_2D_DOMAINOBJ || obj->type == WLZ_3D_DOMAINOBJ ){
      if( (newobj = WlzConvertPix(obj, newpixtype, NULL)) != NULL ){
	(void )WlzWriteObj(stdout, newobj);
	WlzFreeObj(newobj);
      }
    }
    else {
      (void )WlzWriteObj(stdout, obj);
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
		

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
