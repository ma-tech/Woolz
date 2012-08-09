#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzScalarFnObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzScalarFnObj.c
* \author       Bill Hill
* \date         August 2006
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
* \brief
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzscalarfnobj "WlzScalarFnObj"

*/

/*!
\ingroup BinWlz
\defgroup wlzscalarfnobj WlzScalarFnObj
\par Name
WlzScalarFnObj - applies a scalar function to the values of a Woolz
		 object.
\par Synopsis
\verbatim
WlzScalarFnObj [-o<out object>] [-h] 
		       [-e] [-m] [-l] [-s] [-S]
		       [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-e</b></td>
    <td>Exponential function (\f$g_{out} = e^{g_{in}}\f$).</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Modulus function (\f$g_{out} = |g_{in}|\f$) (default).</td>
  </tr>
  <tr> 
    <td><b>-l</b></td>
    <td>Log function (\f$g_{out} = \log(g_{in})\f$).</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Square root function (\f$g_{out} = \sqrt(g_{in})\f$).</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Inverse square root function (\f$g_{out} = \frac{1.0}{\sqrt(g_{in})}\f$).</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>File for the output object.</td>
  </tr>
</table>
\par Description
Applies  scalar function to the values of a Woolz
object.

\par Examples
\verbatim
WlzScalarFnObj -o abs.wlz -m in.wlz
\endverbatim
Reads a Woolz object from in.wlz and writes an output object to abs.wlz.
The output object has the same domain as the input object,
but the values of te output object are the modulus (absolute value)
of the values in the input object.
tie points read from the file points.tie. This transform is then
applied to the object is read from myobj.wlz. The resulting object
is then written to tied.wlz.
\par File
\ref WlzScalarFnObj.c "WlzScalarFnObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzScalarFn "WlzScalarFn(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		ok = 1,
		usage = 0;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  FILE		*fP = NULL;
  WlzFnType	fn = WLZ_FN_SCALAR_MOD;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*inObjFileStr,
  		*outObjFileStr;
  const char    *errMsg;
  static char	optList[] = "ehlmsSo:",
  		inObjFileStrDef[] = "-",
		outObjFileStrDef[] = "-";

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
      case 'e':
        fn = WLZ_FN_SCALAR_EXP;
	break;
      case 'm':
        fn = WLZ_FN_SCALAR_MOD;
	break;
      case 'l':
	fn = WLZ_FN_SCALAR_LOG;
	break;
      case 's':
	fn = WLZ_FN_SCALAR_SQRT;
	break;
      case 'S':
	fn = WLZ_FN_SCALAR_INVSQRT;
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
	((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP)
    {
      if(strcmp(inObjFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    outObj = WlzScalarFn(inObj, fn, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to apply scalar function to object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if(errNum == WLZ_ERR_NONE)
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
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%s",
    *argv,
    " [-o<out object>] [-h>]\n"
    "                  [-e] [-m] [-l]\n"
    "                  [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -e  Exponential function (g_out = exp(g_in)).\n"
    "  -m  Mudulus function (g_out = |g_in|).\n"
    "  -l  Log function (g_out = log(g_in)).\n"
    "  -s  Square root function (g_out = sqrt(g_in)).\n"
    "  -S  Inverse square root function (g_out = 1.0 / sqrt(g_in)).\n"
    "  -o  Output object file name.\n"
    "  -h  Help, prints this usage message.\n"
    "Applies a scalar function to the values of a Woolz object.\n");
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
