#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMinWidthRectAngle_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzMinWidthRectAngle.c
* \author       Bill Hill
* \date         June 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Computes the angle of the minimum width rectangle which
* 		encloses the given 2D object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzminwidthrectangle "WlzMinWidthRectAngle"
*/

/*!
\ingroup BinWlz
\defgroup wlzminwidthrectangle WlzMinWidthRectAngle
\par Name
WlzMinWidthRectAngle - computes the angle of the minimum width enclosing
                       rectangle.
\par Synopsis
\verbatim
WlzMinWidthRectAngle [-o<output file>] [-h] [-r] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Reports usage information.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Output angle using radians rather than the default (degrees).</td>
  </tr>
</table>
\par Description
Computes the angle of the minimum width rectangle which encloses
the given 2D object.
The input object is read from stdin and output data are written
to stdout unless filenames are given
\par Examples
\verbatim
WlzConvexHull in.wlz
\endverbatim
The input Woolz object is read from in.wlz, and the angle that the
longest edge of the minimum enclosing rectangle makes with the
x-axis is written to the standard output.
\par File
\ref WlzMwrAngle.c "WlzMwrAngle.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzObjToConvexHull "WlzObjToConvexHull(3)"
\ref WlzConvexHullToObj "WlzConvexHullToObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int           option,
  		ok = 1,
		radians = 0,
		usage = 0;
  double	angle = 0.0;
  FILE		*fP = NULL;
  char		*inFileStr,
  		*outFileStr;
  WlzObject     *cObj = NULL,
  		*inObj = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hro:";
  const char	outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  outFileStr = (char *)outFileStrDef;
  inFileStr = (char *)inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'r':
        radians = 1;
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
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
        inFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inFileStr == NULL) ||
       (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
	      fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s\n",
		     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      fclose(fP);
    }

  }
  if(ok)
  {
    if(inObj->type == WLZ_CONV_HULL)
    {
      cObj = WlzAssignObject(inObj, &errNum);
    }
    else
    {
      cObj = WlzAssignObject(
             WlzObjToConvexHull(inObj, &errNum), NULL);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		 "%s: Failed to compute or assign convex hull object (%s).\n",
		 argv[0], errMsgStr);
    }
  }
  if(ok)
  {
    angle = WlzMwrAngle(cObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		 "%s: Failed to compute angle (%s).\n",
		 argv[0], errMsgStr);
    }
  }
  if((fP = (strcmp(outFileStr, "-")?
	   fopen(outFileStr, "w"): stdout)) == NULL)
  {
    ok = 0;
    (void )fprintf(stderr,
		   "%s: Failed to open output file %s.\n",
		   argv[0], outFileStr);
  }
  if(ok)
  {
    if(radians == 0)
    {
      angle *= 180.0 / ALG_M_PI;
    }
    (void )fprintf(fP, "%g\n", angle);
    if(feof(fP))
    {
      ok = 0;
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to write output angle (%s).\n",
      		     argv[0], errMsgStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(cObj);
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-o<output object>] [-h] [-o] [-r] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -r  Output angle using radians rather than the default (degrees)\n"
    "Computes the angle of the minimum width rectangle which encloses\n"
    "the given 2D object.\n"
    "The input object is read from stdin and output data are written\n"
    "to stdout unless filenames are given.\n",
    *argv,
    " -r in.wlz\n"
    "The input Woolz object is read from in.wlz, and the angle that\n"
    "the longest edge of the minimum enclosing rectangle makes with the\n"
    "x-axis is written to the standard output.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
