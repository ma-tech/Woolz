#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzHistogramEqualiseObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzHistogramEqualiseObj.c
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
* \brief	Histogram equalises a Woolz domain object.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzhistogramequaliseobj "WlzHistogramEqualiseObj"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzhistogramequaliseobj WlzHistogramEqualiseObj
\par Name
WlzHistogramEqualiseObj - Histogram equalises a Woolz domain object.
\par Synopsis
\verbatim
WlzHistogramEqualiseObj [-o<output object file>] [-h]
           [-D] [-s<smoothing>] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-D</b></td>
    <td>dither pixel value when mapping</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>histogram smoothing factor.</td>
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
By  default  the  input  object is read from the standard input and the
output object  is written to the standard output.
\par Description
Histogram equalises the given Woolz domain object so that the histogram
of  the modified object's grey values approximates a uniform histogram.
\par
The smoothing  factor  is  the  low-pass  gaussian  convolution  kernel
half-height  full-width in bins (default 0), used to smooth the objects
histogram before equalisation.

\par Examples
\verbatim
# An example which uses WlzHistogramEqualiseObj to histogram equalise
# an object.

WlzHistogramEqualiseObj myObj.wlz >myEqualisedObj.wlz

\endverbatim

\par File
\ref WlzHistogramEqualiseObj.c "WlzHistogramEqualiseObj.c"
\par See Also
\ref wlzhistogramobj "WlzHistogramObj(1)"
\ref WlzHistogramEqualiseObj "WlzHistogramEqualiseObj(3)"
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
  int		option,
  		dither = 0,
		smoothing = 0,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "Ds:o:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileStr = outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'D':
	dither = 1;
        break;
      case 's':
        if(sscanf(optarg, "%d", &smoothing) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if((inObj->type != WLZ_2D_DOMAINOBJ) && (inObj->type != WLZ_3D_DOMAINOBJ))
    {
      ok = 0;
      (void )fprintf(stderr,
		 "%s: input object read from file %s is not a domain object\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(inObj, NULL);
    if((errNum = WlzHistogramEqualiseObj(outObj, smoothing,
    				         dither)) != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to histogram equalise object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?  fopen(outObjFileStr, "w"):
	      				    stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(inObj)
  {
    WlzFreeObj(inObj);
  }
  if(outObj)
  {
    WlzFreeObj(outObj);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-D] [-s#] [-o<output object>] [<input object>]\n"
    "Options:\n"
    "  -D  Dither pixel values when mapping.\n"
    "  -s  Input object histogram smoothing.\n"
    "  -o  Output object/data file name.\n"
    "  -h  Help, prints this usage message.\n"
    "Histogram equalises the given Woolz domain object so that the\n"
    "histogram of the modified object's grey values approximates\n"
    "a uniform histogram.\n"
    "The smoothing parameter is a gaussian convolution half height full\n"
    "width specified in histogram bins.\n"
    "Objects/data are read from stdin and written to stdout unless the\n"
    "filenames are given.\n",
    *argv,
    " myobj.wlz -o equalised.wlz\n"
    "The input Woolz domain object is read from myobj.wlz, histogram\n"
    "equalised and written to the file equalised.wlz\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
