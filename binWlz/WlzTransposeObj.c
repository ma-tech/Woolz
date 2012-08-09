#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransposeObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzTransposeObj.c
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
* \brief	Transposes a 2D object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlztransposeobj "WlzTransposeObj"
*/

/*!
\ingroup BinWlz
\defgroup wlztransposeobj WlzTransposeObj
\par Name
WlzTransposeObj - transposes a 2D object.
\par Synopsis
\verbatim
WlzTransposeObj [-h] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Transposes a 2D object.
\par Examples
\verbatim
WlzTransposeObj in.wlz >out.wlz
\endverbatim
Reads an object from in.wlz and then transposes the object,
such that \f$g'(x,y) = g(y,x)\f$.
The resulting object is then written to out.wlz.
\par File
\ref WlzTransposeObj.c "WlzTransposeObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref WlzTransposeObj "WlzTransposeObj(3)"
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
	  "\tTranspose a 2D woolz object\n"
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

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "h";
  int		option;
  const char	*errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  if( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

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
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
    case WLZ_TRANS_OBJ:
    case WLZ_2D_POLYGON:
    case WLZ_BOUNDLIST:
    case WLZ_EMPTY_OBJ:
      if( (nobj = WlzTransposeObj(obj, &errNum)) == NULL ){
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to transpose object (%s).\n",
		       argv[0], errMsg);
      }
      else {
	if((errNum = WlzWriteObj(stdout, nobj)) != WLZ_ERR_NONE) {
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr, "%s: failed to write object (%s).\n",
			 argv[0], errMsg);
	  return(1);
	}
	WlzFreeObj(nobj);
      }
      break;

    default:
      if((errNum = WlzWriteObj(stdout, obj)) != WLZ_ERR_NONE) {
        (void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr, "%s: failed to write object (%s).\n",
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
