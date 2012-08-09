#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLaplacian_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzLaplacian.c
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
* \brief	Applies a Laplacian edge enhancement filter.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzlaplacian "WlzLaplacian"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzlaplacian WlzLaplacian
\par Name
WlzLaplacian - Applies a Laplacian edge enhancement filter.
\par Synopsis
\verbatim
WlzLaplacian [-o<output object file>] [-s kernel size]
           [-a] [-h] [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>absolute value of the convolution. </td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td> convolution kernel size (3, 5 or 7), default 3. </td>
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
output object is written to the standard output.  Unless the -a flag is
set  the  resulting pixel values will be signed (for int or short pixel
values) or centred around 128 for WlzUByte pixel values.

\par Description
Applies an edge enhancement filter to the  input  Woolz  object,  which
must  be a WLZ_2D_DOMAINOBJ with integral (int, short or unsigned char)
grey values.

\par Examples
\verbatim
# An example of using WlzLaplacian with the default 3X3 kernel:

WlzLaplacian -o outfile.wlz infile.wlz
\endverbatim

\par File
\ref WlzLaplacian.c "WlzLaplacian.c"
\par See Also
\ref wlzgauss "WlzGauss(1)"
\ref wlzrsvfilterobj "WlzRsvFilterObj(1)"
\ref wlzgreycrossing "WlzGreyCrossing(1)"
\ref WlzLaplacian "WlzLaplacian(3)"

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
		modulus = 0,
		kSize = 3,
		ok = 1,
		usage = 0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  static char	optList[] = "o:s:ah",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        modulus = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 's':
	if((sscanf(optarg, "%d", &kSize) != 1) ||
	   ((kSize != 3) && (kSize != 5) && (kSize != 7)))
	{
	  usage = 1;
	  ok = 0;
	}
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
    if((outObj = WlzLaplacian(inObj, kSize, 0, modulus, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to Laplacian filter object (%s).\n",
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
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%d%sExample: %s%s",
    *argv,
    " [-o<out object>] [-a] [-s #] [-h] [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -a  Take the absolute value of the convolution.\n"
    "  -s  Laplacian convolution kernel size (3, 5, or 7), default ",
    kSize,
    ".\n"
    "  -h  Help, prints this usage message.\n"
    "Applies a Laplacian edge enhancement filter to the input Woolz object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o enhanced.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, filtered and written\n"
    "to enhanced.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
