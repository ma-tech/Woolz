#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPolygonToObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzPolygonToObj.c
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
* \brief	Converts a polygon object to a domain object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzpolygontoobj "WlzPolygonToObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzpolygontoobj WlzPolygonToObj
\par Name
WlzPolygonToObj  -  converts a polygon object to a domain object.
\par Synopsis
\verbatim
WlzPolygonToObj  [-h] [-v] [-t#] [<input file>]
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
    <td><b>-t</b></td>
    <td>Polygon fill mode
    <table width="500" border="0">
    <tr><td>0</td> <td>Simple fill, default.</td> </tr>
    <tr><td>1</td> <td>Even-odd fill.</td> </tr>
    <tr><td>2</td> <td>Vertex fill.</td> </tr>
    </table></td>
  </tr>
</table>
\par Description
Converts a polygon object to a domain object.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzPolygonToObj.c "WlzPolygonToObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzboundarytoobj "WlzBoundaryToObj(1)"
\ref wlzrasterobj "WlzRasterObj(1)"
\ref WlzPolygonToObj "WlzPolygonToObj(3)"
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
	  "Usage:\t%s [-h] [-v] [-t#] [<input file>]\n"
	  "\tConvert a 2D or 3D polygon woolz object\n"
	  "\tto the corresponding domain object\n"
	  "\tOptions are:\n"
	  "\t  -t#       Polygon fill-mode:\n"
	  "\t            # = %d: simple fill (default)\n"
	  "\t                %d: even-odd fill\n"
	  "\t                %d: vertex-fill\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str, WLZ_SIMPLE_FILL, WLZ_EVEN_ODD_FILL,
	  WLZ_VERTEX_FILL);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "f:hv";
  int		option;
  int		verboseFlg=0;
  WlzPolyFillMode	fillMode=WLZ_SIMPLE_FILL;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'f':
      switch( fillMode = (WlzPolyFillMode) atoi(optarg) ){
      case WLZ_SIMPLE_FILL:
      case WLZ_EVEN_ODD_FILL:
      case WLZ_VERTEX_FILL:
	break;

      default:
	fprintf(stderr, "%s: grey typefill-mode = %d is invalid\n",
		argv[0], fillMode);
        usage(argv[0]);
        return 1;
      }
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
  while((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) 
  {
    switch( obj->type )
    {
    case WLZ_2D_POLYGON:
      if( (nobj = WlzPolygonToObj(obj, WLZ_SIMPLE_FILL,
      				NULL)) != NULL ){
	(void )WlzWriteObj(stdout, nobj);
	WlzFreeObj(nobj);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      if( obj->domain.p->type == WLZ_PLANEDOMAIN_POLYGON ){
	if( (nobj = WlzPolygonToObj(obj, WLZ_SIMPLE_FILL,
				    NULL)) != NULL ){
	  (void )WlzWriteObj(stdout, nobj);
	  WlzFreeObj(nobj);
	}
      }
      else {
	(void )WlzWriteObj(stdout, obj);
      }
      break;
	
    default:
      (void )WlzWriteObj(stdout, obj);
      break;
    }

    WlzFreeObj(obj);
  }

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
