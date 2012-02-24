#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTiffStackToPatches_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzTiffStackToPatches.c
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
* \brief	Converts a TIFF stack to a compound object representing
* 		the set of patches or mosiac tiles.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlztiffstacktopatches "WlzTiffStackToPatches"
*/

/*!
\ingroup BinWlzApp
\defgroup wlztiffstacktopatches WlzTiffStackToPatches
\par Name
WlzTiffStackToPatches - converts a TIFF stack to a compound object representing
                        the set of patches or mosiac tiles.
\par Synopsis
\verbatim
WlzTiffStackToPatches  [-h] [-v] -n cols[,rows] -O overlap <input file>
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
  <tr> 
    <td><b>-n</b></td>
    <td>Number (columns, rows) estimated if omitted.</td>
  </tr>
  <tr> 
    <td><b>-O</b></td>
    <td>Estimated pixel overlap, default 50.</td>
  </tr>
</table>
\par Description
Converts a TIFF stack to a compound image
representing the set of patches or mosaic
tiles. It is assumed that the TIFF images
are arrannged in the stack row-wise starting
top-left with the numbers of rows and columns
input. Also required is the estimated pixel
overlap. The TIFF is read from the given file,
output to stdout.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzTiffStackToPatches.c "WlzTiffStackToPatches.c"
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
	  "Usage:\t%s -n cols[,rows] -O overlap -h -v <input file>\n"
	  "\tConvert a TIFF stack to a compound image\n"
	  "\trepresenting the set of patches or mosaic\n"
	  "\ttiles. It is assumed that the TIFF images\n"
	  "\tare arrannged in the stack row-wise starting\n"
	  "\ttop-left with the numbers of rows and columns\n"
	  "\tinput. Also required is the estimated pixel\n"
	  "\toverlap. The TIFF is read from the given file,\n"
	  "\toutput to stdout.\n"
	  "\tOptions are:\n" 
	  "\t  -n#,#     Number (columns, rows) estimated if omitted\n"
	  "\t  -O#       Estimated pixel overlap (default 50)\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n"
	  "",
	  proc_str);
  return;
}

static int	guessColsRows(
  int	numPatches,
  int	*cols,
  int	*rows)
{
  /* go to maximum 50 since the IP-lab system has a 7x7 max array */
  switch( numPatches ){
  case 1:
    *cols = 1; *rows = 1;
    return 1;
  case 2:
    *cols = 1; *rows =2;
    return 1;
  case 3:
    *cols = 1; *rows = 3;
    return 1;
  case 4:
    *cols = 2; *rows = 2;
    return 1;
  case 6:
    *cols = 2; *rows = 3;
    return 1;
  case 8:
    *cols = 2; *rows = 4;
    return 1;
  case 9:
    *cols = 3; *rows = 3;
    return 1;
  case 10:
    *cols = 2; *rows = 5;
    return 1;
  case 12:
    *cols = 3; *rows = 4;
    return 1;
  case 15:
    *cols = 3; *rows = 5;
    return 1;
  case 16:
    *cols = 4; *rows = 4;
    return 1;
  case 18:
    *cols = 3; *rows = 6;
    return 1;
  case 20:
    *cols = 4; *rows = 5;
    return 1;
  case 21:
    *cols = 3; *rows = 7;
    return 1;
  case 24:
    *cols = 4; *rows = 6;
    return 1;
  case 25:
    *cols = 5; *rows = 5;
    return 1;
  case 27:
    *cols = 3; *rows = 7;
    return 1;
  case 28:
    *cols = 4; *rows = 7;
    return 1;
  case 30:
    *cols = 5; *rows = 6;
    return 1;
  case 32:
    *cols = 4; *rows = 8;
    return 1;
  case 35:
    *cols = 5; *rows = 7;
    return 1;
  case 36:
    *cols = 6; *rows = 6;
    return 1;
  case 40:
    *cols = 5; *rows = 8;
    return 1;
  case 42:
    *cols = 6; *rows = 7;
    return 1;
  case 44:
    *cols = 4; *rows = 11;
    return 1;
  case 45:
    *cols = 5; *rows = 9;
    return 1;
  case 48:
    *cols = 6; *rows = 8;
    return 1;
  case 49:
    *cols = 7; *rows = 7;
    return 1;
  case 50:
    *cols = 5; *rows = 10;
    return 1;

    /* prime numbers */
  case 5:
  case 7:
  case 11:
  case 13:
  case 17:
  case 19:
  case 23:
  case 29:
  case 31:
  case 37:
  case 41:
  case 43:
  case 47:
    *cols = 1; *rows = numPatches;
    return 0;

    /* prime numbers * 2 */
  case 14:
  case 22:
  case 26:
  case 34:
  case 38:
  case 46:
    *cols = 2; *rows = numPatches/2;
    return 0;

    /* print numbers times 3 */
  case 33:
  case 39:
  case 51:
    *cols = 3; *rows = numPatches/3;
    return 0;

    /* the rest > 51 */
  default:
    *cols = 1; *rows = numPatches;
    return 0;
  }
}

int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "n:O:hv";
  int		option;
  int		verboseFlg=0;
  int		rows=0, cols=0;
  int		overlap=50;
  int		rowShift, colShift;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  char		*tiffFile=NULL;
  WlzObject	*inObj=NULL, *outObj=NULL, *obj, **objVec;
  WlzCompoundArray	*outCmpnd;
  int		i, j, objVecCount, objIndx;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'n':
      switch( sscanf(optarg, "%d,%d", &cols, &rows) ){
      default:
	usage(argv[0]);
	return 1;

      case 1:
	fprintf(stderr,
		"%s: number of rows not set and will be estimated\n",
		argv[0]);
	break;

      case 2:
	break;
      }
      break;

    case 'O':
      switch( sscanf(optarg, "%d", &overlap) ){
      default:
	fprintf(stderr,
		"%s: overlap not read", argv[0]);
	usage(argv[0]);
	return 1;

      case 1:
	break;
      }
      break;

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
    fprintf(stderr, "%s: input TIFF file required\n", argv[0]);
    usage(argv[0]);
    return WLZ_ERR_UNSPECIFIED;
  }

  /* read the TIFF file */
  if( (inObj = WlzAssignObject(
	 WlzEffReadObj(NULL, tiffFile, WLZEFF_FORMAT_TIFF, 0,
		       &errNum), NULL)) == NULL ){
    usage(argv[0]);
    return errNum;
  }

  /* if 2D then there is only one patch, simply output as a compound object */
  switch( inObj->type ){
  case WLZ_2D_DOMAINOBJ:
    outCmpnd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 3, 1, &inObj,
				    inObj->type, &errNum);
      outObj = (WlzObject *) outCmpnd;
    break;

  case WLZ_3D_DOMAINOBJ:
    if( (errNum = WlzExplode3D(&objVecCount, &objVec, inObj)) != WLZ_ERR_NONE ){
      usage(argv[0]);
      return errNum;
    }
    /* check if rows and columns plausible */
    if( cols && rows ){
      if( cols*rows != objVecCount ){
	fprintf(stderr,
		"%s: input columns and rows (%d,%d) does not match"
		"number of patches - %d\n Continuing anyway\n",
		argv[0], cols, rows, objVecCount);
      }
    }
    else if( cols ){
      rows = objVecCount / cols;
      if( cols*rows != objVecCount ){
	fprintf(stderr,
		"%s: input columns and estimated rows (%d,%d) does not match"
		"number of patches - %d\nContinuing anyway\n",
		argv[0], cols, rows, objVecCount);
      }
    }
    else {
      if( guessColsRows(objVecCount, &cols, &rows) ){
	fprintf(stderr,
		"%s: estimated cols and rows - (%d,%d)\nContinuing\n",
		argv[0], cols, rows);
      }
      else {
	fprintf(stderr, 
		"%s: cols,rows (%d, %d) estimated seems crazy,"
		"put values in explicitely\n",
		argv[0], cols, rows);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
    }

    /* shift patch objects, assume patches all the same size */
    rowShift = (objVec[0]->domain.i->lastln - objVec[0]->domain.i->line1
		+ 1 - overlap);
    colShift = (objVec[0]->domain.i->lastkl - objVec[0]->domain.i->kol1
		+ 1 - overlap);
    objIndx = 0;
    for(j=0; j < rows; j++){
      for(i=0; i < cols; i++, objIndx++){
	if( objIndx < objVecCount ){
	  obj = WlzShiftObject(objVec[objIndx], 
			       colShift*i, rowShift*j, 0, &errNum);
	  WlzFreeObj(objVec[objIndx]);
	  objVec[objIndx] = WlzAssignObject(obj, &errNum);
	}
      }
    }

    /* make array */
    outCmpnd = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 2,
				    objVecCount, objVec,
				    objVec[0]->type, &errNum);
    outObj = (WlzObject *) outCmpnd;
    break;

  default:
    errNum = WLZ_ERR_OBJECT_TYPE;
  }
  WlzFreeObj(inObj);

  /* write compound object */
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
