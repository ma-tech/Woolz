#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTiffStackToShade_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzTiffStackToShade.c
* \author	Richard Baldock
* \date         June 2004
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
* \brief	Converts a TIFF stack to a shade object.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlztiffstacktoshade "WlzTiffStackToShade"
*/

/*!
\ingroup BinWlzApp
\defgroup wlztiffstacktoshade WlzTiffStackToShade
\par Name
WlzTiffStackToShade - converts a TIFF stack to a shade object.
\par Synopsis
\verbatim
WlzTiffStackToShade [-h] [-v] <input file>
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose output.</td>
  </tr>
</table>
\par Description
Converts a TIFF stack to a shade image by finding the
maximal pixel values in each channel.
The TIFF is read from the given file, output to stdout.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzTiffStackToShade.c "WlzTiffStackToShade.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>
#include <WlzExtFF.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s <input file>\n"
	  "\tConvert a TIFF stack to a shade image\n"
	  "\tby finding the maximal pixel values in each\n"
	  "\tchannel. The TIFF is read from the given file,\n"
	  "\toutput to stdout.\n"
	  "\tOptions are:\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "hv";
  int		option;
  int		verboseFlg=0;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*tiffFile=NULL;
  WlzObject	*inObj, *outObj, *obj, **objVec;
  int		i, objVecCount;
  WlzGreyType	gType;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  if( optind < argc ){
    tiffFile = *(argv+optind);
  }
  else {
    usage(argv[0]);
    return WLZ_ERR_UNSPECIFIED;
  }

  /* read the TIFF file */
  if( (inObj = WlzAssignObject(WlzEffReadObj(NULL, tiffFile,
  			                     WLZEFF_FORMAT_TIFF, 0, 0, 0,
					     &errNum), NULL)) == NULL ){
    usage(argv[0]);
    return errNum;
  }

  /* if 2D then that is the shade file, else take max of 3D stack */
  switch( inObj->type ){
  case WLZ_2D_DOMAINOBJ:
    outObj = inObj;
    break;

  case WLZ_3D_DOMAINOBJ:
    if( (errNum = WlzExplode3D(&objVecCount, &objVec, inObj)) != WLZ_ERR_NONE ){
      usage(argv[0]);
      return errNum;
    }
    gType = WlzGreyTypeFromObj(objVec[0], &errNum);
    outObj = WlzAssignObject(objVec[0], NULL);
    for(i=1; i < objVecCount; i++){
      if( gType == WLZ_GREY_RGBA ){
	obj = WlzRGBAImageArithmetic(outObj, objVec[i], WLZ_BO_MAX, 0, &errNum);
      }
      else {
	obj = WlzImageArithmetic(outObj, objVec[i], WLZ_BO_MAX, 0, &errNum);
      }
      if( obj ){
	WlzFreeObj(outObj);
	outObj = WlzAssignObject(obj, &errNum);
      }
      else {
	break;
      }
    }
    
    WlzFreeObj(inObj);
    break;

  default:
    WlzFreeObj(inObj);
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  /* write shade object */
  if( errNum == WLZ_ERR_NONE ){
    WlzWriteObj(stdout, outObj);
    WlzFreeObj(outObj);
  }
  else {
    usage(argv[0]);
  }
  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
