#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzGreyNormalise_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzGreyNormalise.c
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
* \brief	Normalise the grey-range of a grey-level object.
* \ingroup	binWlz
*
* \par Binary
* \ref wlzgreynormalise "WlzGreyNormalise"
*/


/*!
\ingroup BinWlz
\defgroup wlzgreynormalise WlzGreyNormalise
\par Name
WlzGreyNormalise - normalise the grey-range of a grey-level woolz object.
\par Synopsis
\verbatim
WlzGreyNormalise [-h] [-d] [-v] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>dither values.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>verbatim operation.</td>
  </tr>
</table>
\par Description
tNormalise the grey-range of a grey-level woolz object
which resets the grey range to (0,255) and is a
convenience equivalent to WlzGreySetRange -U255 -L0.
\par Examples
\verbatim
WlzGreySetRange in.wlz >out.wlz
\endverbatim
Normalises the grey values of the object read from in.wlz to the range
0 - 255 and then writes the resulting object to out.wlz.
\par File
\ref WlzGreyNormalise.c "WlzGreyNormalise.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreysetrange "WlzGreySetRange(1)"
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
	  "Usage:\t%s [-h] [-v] [<input file>]\n"
	  "\tNormalise the grey-range of a grey-level woolz object\n"
	  "\twhich is resets the grey range to (0,255) and is a\n"
	  "\tconvenience routine equivalent to WlzGreySetRange -U255 -L0\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -d        dither values\n"
	  "\t  -u        convert to unsigne byte values\n"
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
  char 		optList[] = "dhu";
  int		option,
		dither = 0,
  		ubyteFlg = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'd':
      dither = 1;
      break;
    case 'u':
      ubyteFlg = 1;
      break;
    case 'h': /* FALLTHROUGH */
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
	if((errNum = WlzGreyNormalise(obj, dither)) == WLZ_ERR_NONE ){
	  if(ubyteFlg) {
	    WlzGreyType gType;
	    gType = WlzGreyTypeFromObj(obj, &errNum);
	    if((errNum == WLZ_ERR_NONE) && (gType != WLZ_GREY_UBYTE)){
	      WlzObject *obj1;
	      obj1 = WlzConvertPix(obj, WLZ_GREY_UBYTE, &errNum);
	      if(errNum == WLZ_ERR_NONE){
	        (void )WlzFreeObj(obj);
		obj = obj1;
	      }
	    }
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
    (void )WlzFreeObj(obj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
