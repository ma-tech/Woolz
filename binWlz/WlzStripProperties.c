#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzStripProperties.c
* \author       Richard Baldock
* \date         August 2002
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
* \brief	Strips the property list from an object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzstripproperties "WlzStripProperties"
*/

/*!
\ingroup BinWlz
\defgroup wlzstripproperties WlzStripProperties
\par Name
WlzStripProperties - strips the property list from an object.
\par Synopsis
\verbatim
WlzStripProperties [-h] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Strips the property list from an object.
This could be useful if you need to read an object
from a newer version of woolz with old code e.g. if
you have binaries for a different architecture than
supplied. Object written to stdout.
\par Examples
\verbatim
WlzStripProperties obj.wlz >obj_s.wlz
\endverbatim
Reads an object from obj.wlz and writes it back out,
without a property list,
to obj_s.wlz.
\par File
\ref WlzStripProperties.c "WlzStripProperties.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzchangeemapproperty "WlzChangeEMAPProperty(1)"
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
	  "\tStrip the property list from a woolz object.\n"
	  "\tThis could be useful if you need to read an object\n"
	  "\tfrom a newer version of woolz with old code e.g. if\n"
	  "\tyou have binaries for a different architecture than\n"
	  "\tsupplied. Object written to stdout\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'h':
      usage(argv[0]);
      return( 0 );

    case 'v':
      verboseFlg = 1;
      break;

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

  /* read objects and strip properties */
  while( (obj = WlzReadObj(inFile, NULL)) != NULL ){
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_3D_WARP_TRANS:
    case WLZ_PROPERTY_OBJ:
      if( obj->plist ){
	WlzFreePropertyList(obj->plist);
	obj->plist = NULL;
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      if(((WlzCompoundArray *)obj)->plist){
	WlzFreePropertyList(((WlzCompoundArray *)obj)->plist);
	((WlzCompoundArray *)obj)->plist = NULL;
      }
      break;

    default:
      break;
    }

    (void )WlzWriteObj(stdout, obj);
    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
		

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
