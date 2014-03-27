#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzSobel_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzSobel.c
* \author       Bill Hill
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
* \brief	A 3\f$\times\f$3 Sobel edge detection filter.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzsobel "WlzSobel"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzsobel WlzSobel
\par Name
WlzSobel - a 3\f$\times\f$3 Sobel edge detection filter.
\par Synopsis
\verbatim
WlzSobel [-o<output object file>] [-x] [-y] [-h]
           [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-x</b></td>
    <td>no horizontal pass. </td>
  </tr>
  <tr>
    <td><b>-y</b></td>
    <td>no vertical pass. </td>
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
output object is written to the standard output.
\par Description
Applies a 3X3 Sobel edge detection filter to the  input  Woolz  object,
which  must be a WLZ_2D_DOMAINOBJ with integral (int, short or unsigned
char) grey values.

\par Examples
\verbatim
# An example of using WlzSobel to detect only vertical edges.

WlzSobel -o outfile.wlz -x infile.wlz

# A simple example in which both horizontal and vertical edges
# are detected.

WlzSobel <infile.wlz >outfile.wlz
\endverbatim
\par File
\ref WlzSobel.c "WlzSobel.c"
\par See Also
\ref WlzSobel "WlzSobel(3)"
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
		hFlag = 1,
		vFlag = 1,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  static char	optList[] = "xyo:h",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'x':
        hFlag = 0;
	break;
      case 'y':
        vFlag = 0;
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
    if((outObj = WlzSobel(inObj, hFlag, vFlag, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to Sobel filter object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
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
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-x] [-y] [-h] [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -x  No horizontal pass.\n"
    "  -y  No vertical pass.\n"
    "  -h  Help, prints this usage message.\n"
    "Applies a Sobel edge detection filter to the input Woolz object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o edge.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, filtered and written\n"
    "to edge.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
