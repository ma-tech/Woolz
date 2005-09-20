#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlz/WlzIntRescaleObj.c
* \author       Bill Hill
* \date         March 2003
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
* \brief	Rescales an object using an integer scale factor.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzintrescaleobj "WlzIntRescaleObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzintrescaleobj WlzIntRescaleObj
\par Name
WlzIntRescaleObj - rescales an object using an integer scale factor.
\par Synopsis
\verbatim
WlzIntRescaleObj [-s#] [-c] [-e] [-h] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Scale factor, default 1.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Compress using 1/scale, default false.</td>
  </tr>
  <tr> 
    <td><b>-e</b></td>
    <td>Expand using scale, default true.</td>
  </tr>
</table>
\par Description
Rescales a woolz object using an integer scale factor.
\par Examples
\verbatim
WlzIntRescaleObj -o out.wlz -s 2 -c in.wlz
\endverbatim
Subsamples the object read from in.wlz using integer scaling so that 
the linear dimensions of the output object (written to out.wlz) 
are halved.
\par File
\ref WlzIntRescaleObj.c "WlzIntRescaleObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref wlzsampleobj "WlzSampleObj(1)"
\ref WlzIntRescaleObj "WlzIntRescaleObj(3)"
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
  int		expand = 1,
  		scale = 1,
  		option,
  		ok = 1,
		usage = 0;
  char		*outObjFileStr,
  		*inObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  const char	*errMsg;
  static char	optList[] = "o:s:ceh",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileStr = outObjFileStrDef;
  while(ok && (usage == 0) &&
  	((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'c':
	expand = 0;
	break;
      case 'e':
	expand = 1;
	break;
      case 's':
        if(sscanf(optarg, "%d", &scale) != 1)
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
    if((outObj = WlzIntRescaleObj(inObj, scale, expand, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to rescale object (%s).\n",
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
    "Usage: "
    "%s [-s#] [-c] [-e] [-h] [<in object>]\n"
    "Rescales a woolz object using an integer scale factor.\n"
    "Options are:\n"
    "  -s#   Scale factor (set to %d).\n"
    "  -c    Compress using 1/scale (%sset).\n"
    "  -e    Expand using scale (%sset).\n"
    "  -h    Display this usage information.\n",
    *argv,
    scale,
    (expand == 0) ? "" : "not ",
    (expand != 0) ? "" : "not ");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
