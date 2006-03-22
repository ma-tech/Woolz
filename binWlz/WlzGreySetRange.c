#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreySetRange_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzGreySetRange.c
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
* \brief	Sets the grey-range of a grey-level woolz object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzgreysetrange "WlzGreySetRange"
*/

/*!
\ingroup BinWlz
\defgroup wlzgreysetrange WlzGreySetRange
\par Name
WlzGreySetRange  -  sets the grey-range of a grey-level woolz object.
\par Synopsis
\verbatim
WlzGreySetRange [-u#] [-l#] [-U#] [-L#] [-h] [-v] [<input file>]
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
    <td><b>-l</b></td>
    <td>Low grey value in source image, default min value in source.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Low grey value in destination image, default 0.</td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Upper grey value in source image, default max value in source.</td>
  </tr>
  <tr> 
    <td><b>-U</b></td>
    <td>Upper grey value in dest image, default 0.</td>
  </tr>
</table>
\par Description
Sets the grey-range of a grey-level woolz object
riting the new object to standard output.
The original grey-values are linearly transformed
from their original range. The transformation is
\f[
g_{out} = (\frac{U - L}{u - l})g_{in} + L
new_grey = ((U - L)/(u - l))*old_grey + L
\f]
where
\f$U\f$ and \f$L\f$ are the new upper and lower values and
\f$u\f$ and \f$l\f$ are the old upper and lower values.
If \f$u\f$ and \f$l\f$
are set it is up to the user to ensure that the actual
image values are within the range \f$[u,l]\f$.
Use WlzGreyRange to check.
\par Examples
\verbatim
\endverbatim
\par File
\ref WlzGreySetRange.c "WlzGreySetRange.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzgreyrange "WlzGreyRange(1)"
\ref WlzGreySetRange "WlzGreySetRange(3)"
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
	  "Usage:\t%s [-u#] [-l#] [-U#] [-L#] [-h] [-v] [<input file>]\n"
	  "\tReset the grey-range of a grey-level woolz object\n"
	  "\twriting the new object to standard output\n"
	  "\tThe original grey-values are linearly transformed\n"
	  "\tfrom their original range. The transformation is \n"
	  "\tnew_grey = ((U - L)/(u - l))*old_grey + L where\n"
	  "\t U and L are the new upper and lower values and \n"
	  "\t u and l are the old upper and lower values. If u and l\n"
	  "\t are set it is up to the user to ensure that the actual\n"
	  "\timage values are within the range [u,l]. Use WlzGreyRange\n"
	  "\tto check."
	  "\tOptions are:\n"
	  "\t  -l#       low grey value in source image, default min value"
	  "in source\n"
	  "\t  -L#       low grey value in dest image, default 0"
	  "\t  -u#       upper grey value in source image, default max value"
	  "in source\n"
	  "\t  -U#       upper grey value in dest image, default 0"
	  "in source\n"
	  "\t  -h        help - prints this usage message\n"
	  "\t  -v        verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj;
  FILE		*inFile;
  char 		optList[] = "l:L:u:U:hv";
  int		option;
  WlzPixelV	max, min, Max, Min;
  WlzPixelV	gmin, gmax;
  int		getminFlg=1, getmaxFlg=1, verboseFlg=0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  int		objCount=0;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  Max.type = Min.type = max.type = min.type = WLZ_GREY_DOUBLE;
  Max.v.dbv = 255.0;
  Min.v.dbv = 0.0;
  
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'l':
      min.v.dbv = atof(optarg);
      getminFlg = 0;
      break;

    case 'L':
      Min.v.dbv = atof(optarg);
      break;

    case 'u':
      max.v.dbv = atof(optarg);
      getmaxFlg = 0;
      break;

    case 'U':
      Max.v.dbv = atof(optarg);
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
  while(((obj = WlzAssignObject(WlzReadObj(inFile, NULL), NULL)) != NULL) &&
        (errNum == WLZ_ERR_NONE))
  {
    objCount++;
    switch( obj->type )
    {
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      /* get the existing min and max grey values */
      if( (errNum = WlzGreyRange(obj, &gmin, &gmax)) == WLZ_ERR_NONE ){
	if( getminFlg ){
	  WlzValueConvertPixel(&min, gmin, WLZ_GREY_DOUBLE);
	}
	if( getmaxFlg ){
	  WlzValueConvertPixel(&max, gmax, WLZ_GREY_DOUBLE);
	}

	if( verboseFlg ){
	  fprintf(stderr,
		  "%s:\nconverting object %d with parameters:\n"
		  "\tSource (l, u) = (%f, %f)\n"
		  "\tDestination (L, U) = (%f, %f)\n",
		  argv[0], objCount, min.v.dbv, max.v.dbv,
		  Min.v.dbv, Max.v.dbv);
	}

	errNum = WlzGreySetRange(obj, min, max, Min, Max);
	if( errNum == WLZ_ERR_NONE ){
	  errNum = WlzWriteObj(stdout, obj);
	}
      }
      if( errNum != WLZ_ERR_NONE ){
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: failed to write object (%s).\n",
		       argv[0], errMsg);
	return(errNum);
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

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
