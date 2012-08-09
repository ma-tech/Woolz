#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzThreshold_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzThreshold.c
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
* \brief	Thresholds a grey-level object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzthreshold "WlzThreshold"
*/

/*!
\ingroup BinWlz
\defgroup wlzthreshold WlzThreshold
\par Name
WlzThreshold - thresholds a grey-level object.
\par Synopsis
\verbatim
WlzThreshold [-h] [-t#] [-v#] [-H] [-L] [-E] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-t</b></td> <td>Threshold pixel type:
      <table width="500" border="0">
	<tr><td>1</td> <td>integer, default</td></tr>
	<tr><td>2</td> <td>short</td></tr>
	<tr><td>3</td> <td>unsigned byte</td></tr>
	<tr><td>4</td> <td>float</td></tr>
	<tr><td>5</td> <td>double</td></tr>
      </table>
    </td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Threshold value, default 170.</td>
  </tr>
  <tr> 
    <td><b>-H</b></td>
    <td>Threshold high,
        keep pixels at or above threshold value,
	default.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Threshold low,
        keep pixels below threshold value.</td>
  </tr>
  <tr> 
    <td><b>-E</b></td>
    <td>Threshold equal,
        keep pixels equal to threshold value.</td>
  </tr>
</table>
\par Description
Thresholds a grey-level object writing the new object to the standard output.
\par Examples
\verbatim
WlzThreshold -v 100 -H in.wlz >out.wlz
\endverbatim
Reads an grey valued domain object from in.wlz,
thresholds the object so that it's domain consists of all locations at which
the grey value is at or above the threshold value.
The resulting object is written to out.wlz.
\par File
\ref WlzThreshold.c "WlzThreshold.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzThreshold "WlzThreshold(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

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
      "Usage:\t%s [-t#] [-v#] [-H] [-L] [-E] [-h] [<input file>]\n"
      "\tThreshold a grey-level woolz object\n"
      "\twriting the new object to standard output\n"
      "Version: %s\n"
      "Options:\n"
      "\t  -H        Threshold high, keep pixels above threshold value "
      "(default).\n"
      "\t  -L        Threshold low, keep pixels below threshold value.\n"
      "\t  -E        Threshold equal, keep pixels equal to threshold value.\n"
      "\t  -t#       Threshold pixel type:\n"
      "\t            # = %d: integer (default)\n"
      "\t                %d: short\n"
      "\t                %d: unsigned byte\n"
      "\t                %d: float\n"
      "\t                %d: double\n"
      "\t            Note -t option must precede -v\n"
      "\t  -v#       threshold value  - integer unless -t used\n"
      "\t  -h        Help - prints this usage message\n",
      proc_str,
      WlzVersion(),
      WLZ_GREY_INT, WLZ_GREY_SHORT, WLZ_GREY_UBYTE,
      WLZ_GREY_FLOAT, WLZ_GREY_DOUBLE);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFile;
  char 		optList[] = "HLEht:v:";
  int		option;
  WlzThresholdType highLow = WLZ_THRESH_HIGH;
  WlzGreyType	threshpixtype = WLZ_GREY_INT;
  WlzPixelV	thresh;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  thresh.type = threshpixtype;
  thresh.v.inv = 170;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 't':
      switch( threshpixtype = (WlzGreyType) atoi(optarg) ){

      case WLZ_GREY_INT: /* FALLTHROUGH */
      case WLZ_GREY_SHORT: /* FALLTHROUGH */
      case WLZ_GREY_UBYTE: /* FALLTHROUGH */
      case WLZ_GREY_FLOAT: /* FALLTHROUGH */
      case WLZ_GREY_DOUBLE:
	break;

      default:
        fprintf(stderr, "%s: grey type = %d is invalid\n",
		argv[0], threshpixtype );
        usage(argv[0]);
        return 1;

      }
      break;

    case 'v':
      switch( threshpixtype ){

      case WLZ_GREY_INT:
	thresh.v.inv = atoi(optarg);
	break;

      case WLZ_GREY_SHORT:
	thresh.v.shv = atoi(optarg);
	break;

      case WLZ_GREY_UBYTE:
	thresh.v.ubv = atoi(optarg);
	break;

      case WLZ_GREY_FLOAT:
	thresh.v.flv = atof(optarg);
	break;

      case WLZ_GREY_DOUBLE:
	thresh.v.dbv = atof(optarg);
	break;

      default:
        break;
      }
      break;

    case 'H':
      highLow = WLZ_THRESH_HIGH;
      break;

    case 'L':
      highLow = WLZ_THRESH_LOW;
      break;

    case 'E':
      highLow = WLZ_THRESH_EQUAL;
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
	if( (nobj = WlzThreshold(obj, thresh,
			         highLow, &errNum)) == NULL ){
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr, "%s: failed to threshold object (%s).\n",
	  		 argv[0], errMsg);
	  return(1);
	}
	else {
	  if((errNum = WlzWriteObj(stdout, nobj)) != WLZ_ERR_NONE) {
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
	    		   "%s: failed to write thresholded object (%s).\n",
			   argv[0], errMsg);
	    return(1);
	  }
	  WlzFreeObj(nobj);
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
  if((inFile != NULL) && (inFile != stdin)){
    (void )fclose(inFile);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
