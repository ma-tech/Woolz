#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzGreyStats_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzGreyStats.c
* \author       Bill Hill
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
* \brief	Calculates simple statistics for a domain object's grey values.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzgreystats "WlzGreyStats"
*/


/*!
\ingroup      BinWlz
\defgroup     wlzgreystats WlzGreyStats
\par Name
WlzGreyStats - calculates simple statistics for a domain object's grey values.
\par Synopsis
\verbatim
WlzGreyStats [-o<output file>] [-c] [-h] [-v] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-c</b></td>
    <td>Colour stats if RGB data else modulus data returned.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
By default the input object is read from the  standard  input  and  the
output is written to the standard output.

\par Description
WlzGreyStats  calculates  the  area,  minimum,  maximum,  sum,  sum  of
squares, mean and the standard deviation of the input Woolz  2D  or  3D
domain  object's  grey  values.   If the verbose output flag is not set
then the following are written to the  output  file  in  the  following
order:
\verbatim
 <area> <grey type> <min> <max> <sum> <sum of sq> <mean> <std dev>
\endverbatim

\par Examples
\verbatim
# WlzGreyStats reads an object from myobj.wlz, calculates
# its grey value statistics and then prints them to the
# standard output.

WlzGreyStats myobj.wlz
\endverbatim

\par File
\ref WlzGreyStats.c "WlzGreyStats.c"
\par See Also
\ref wlzgreyrange "WlzGreyRange(1)"
\ref WlzGreyStats "WlzGreyStats(3)"
\ref WlzRGBAGreyStats "WlzRGBAGreyStats(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		area,
  		option,
		ok = 1,
		verbose = 0,
                colFlg = 0,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL;
  double	min,
  		max,
		sum,
		sumSq,
		mean,
		stdDev;
  double	minA[4],
  		maxA[4],
		sumA[4],
		sumSqA[4],
		meanA[4],
		stdDevA[4];
  WlzGreyType	gType;
  char 		*outFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "co:vh",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
    case 'c':
      colFlg = 1;
      break;
    case 'o':
      outFileStr = optarg;
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h':
    default:
      usage = 1;
      ok = 0;
      break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzReadObj(fP, &errNum)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr,
		     errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if( colFlg ){
      area = WlzRGBAGreyStats(inObj, WLZ_RGBA_SPACE_RGB, &gType,	
			      minA, maxA, sumA, sumSqA, meanA,
			      stdDevA, &errNum);
    }
    else {
      area = WlzGreyStats(inObj, &gType, &min, &max, &sum, &sumSq, &mean,
			  &stdDev, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to compute object's grey value statistics"
		     " (%s)\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to write output file\n",
		     *argv);
    }
    else
    {
      if( colFlg ){
	if(verbose)
	{
	  (void )fprintf(fP,
			 "area        %d\n"
			 "grey type   %s\n"
			 "min         (%g, %g, %g, %g)\n"
			 "max         (%g, %g, %g, %g)\n"
			 "sum         (%g, %g, %g, %g)\n"
			 "sum sq      (%g, %g, %g, %g)\n"
			 "mean        (%g, %g, %g, %g)\n"
			 "std dev     (%g, %g, %g, %g)\n",
			 area, WlzStringFromGreyType(gType, NULL),
			 minA[0], minA[1], minA[2], minA[3],
			 maxA[0], maxA[1], maxA[2], maxA[3],
			 sumA[0], sumA[1], sumA[2], sumA[3],
			 sumSqA[0], sumSqA[1], sumSqA[2], sumSqA[3],
			 meanA[0], meanA[1], meanA[2], meanA[3],
			 stdDevA[0], stdDevA[1], stdDevA[2], stdDevA[3]);
	}
	else
	{
	  (void )fprintf(fP,
			 "%d %s (%g,%g,%g,%g) (%g %g %g %g)\n"
			 "(%g,%g,%g,%g) (%g %g %g %g)\n"
			 "(%g,%g,%g,%g) (%g %g %g %g)\n",
			 area, WlzStringFromGreyType(gType, NULL),
			 minA[0], minA[1], minA[2], minA[3],
			 maxA[0], maxA[1], maxA[2], maxA[3],
			 sumA[0], sumA[1], sumA[2], sumA[3],
			 sumSqA[0], sumSqA[1], sumSqA[2], sumSqA[3],
			 meanA[0], meanA[1], meanA[2], meanA[3],
			 stdDevA[0], stdDevA[1], stdDevA[2], stdDevA[3]);
	}
      }
      else {
	if(verbose)
	{
	  (void )fprintf(fP,
			 "area        % 14d\n"
			 "grey type    %14s\n"
			 "min         % 14g\n"
			 "max         % 14g\n"
			 "sum         % 14g\n"
			 "sum sq      % 14g\n"
			 "mean        % 14g\n"
			 "std dev     % 14g\n",
			 area, WlzStringFromGreyType(gType, NULL),
			 min, max, sum, sumSq, mean, stdDev);
	}
	else
	{
	  (void )fprintf(fP, "%d %s %g %g %g %g %g %g\n",
			 area, WlzStringFromGreyType(gType, NULL),
			 min, max, sum, sumSq, mean, stdDev);
	}
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out file>] [-c] [-v] [-h] [<in object>]\n"
    "Options:\n"
    "  -c  Colour stats if RGB, otherwise modulus stats calculated\n"
    "  -o  Output file name.\n"
    "  -v  Verbose output flag.\n"
    "  -h  Help, prints this usage message.\n"
    "Calculates the area, minimum, maximum, sum, sum of squares, mean\n"
    "and the standard deviation of the input 2D or 3D Woolz object's\n"
    "grey values.\n"
    "If the verbose output flag is not set then the following are written\n"
    "to the output file in the following order:\n"
    "  <area> <grey type> <min> <max> <sum> <sum of sq> <mean> <std dev>.\n"
    "The input object is read from stdin and values are written to stdout\n"
    "unless the filenames are given.\n",
    *argv,
    " -o stats.txt myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz. The statistics are\n"
    "calculated and  written to out.txt.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
