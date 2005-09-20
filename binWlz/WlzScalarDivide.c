#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzScalarDivide.c
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
* \brief	Divides the pixel values of a grey-level object by a
* 		scalar value.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzscalardivide "WlzScalarDivide"
*/

/*!
\ingroup BinWlz
\defgroup wlzscalardivide WlzScalarDivide
\par Name
WlzScalarDivide - Divide the pixel values of a grey-level object.
\par Synopsis
\verbatim
WlzScalarDivide -v#] [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Divisor value (0 is an error), default 1.0.</td>
  </tr>
</table>
\par Description
Divide the pixel values of a grey-level woolz object
by a scalar value.
\par Examples
\verbatim
WlzScalarDivide -v 0.5 in.wlz >out.wlz
\endverbatim
The grey values of the object read from in.wlz are halved in value and the
resulting object is written to out.wlz.
\par File
\ref WlzScalarDivide.c "WlzScalarDivide.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreysetrange "WlzGreySetRange(1)"
\ref WlzScalarDivide "WlzScalarDivide(3)"
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
	  "Usage:\t%s [-v#] [-h] [<input file>]\n"
	  "\tDivide the pixel values of a grey-level woolz object\n"
	  "\tby a scalar value\n"
	  "\tOptions are:\n"
	  "\t  -v#       divisor value (0 is an error), default 1.0:\n"
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
  char 		optList[] = "v:";
  int		option;
  WlzPixelV	div;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  div.type = WLZ_GREY_DOUBLE;
  div.v.dbv = 1.0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      div.v.dbv = atof(optarg);
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
  while((obj = WlzReadObj(inFile, NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
      if((obj = WlzScalarDivide(obj, div, NULL)) == NULL ){
	WlzWriteObj(stdout, obj);
      }
      break;

    default:
      WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
