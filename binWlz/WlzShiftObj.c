#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzShiftObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzShiftObj.c
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
* \brief	Shifts objects using an integer translation.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzshiftobj "WlzShiftObj"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzshiftobj WlzShiftObj
\par Name
WlzShiftObj - shifts objects using an integer translation.
\par Synopsis
\verbatim
WlzShiftObj [-h] [-o<output object file>]
            [-g] [-x<column shift>] [-y<line shift>] [-z<plane shift>]
            [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-g</b></td>
    <td>shift object to the origin. </td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>column shift. </td>
  </tr>
  <tr>
    <td><b>-y</b></td>
    <td>line shift. </td>
  </tr>
  <tr>
    <td><b>-z</b></td>
    <td>plane shift. </td>
  </tr>
</table>
\par Description
Shifts a Woolz object without copying it's grey values.
By  default  the  input  object is read from the standard input and the
output object is written to the standard output.
\par Examples
\verbatim
# The following c shell script uses awk(1), WlzBoundingBox(1) and
# WlzShiftObj(1) to shift infile.wlz so that its first
# column, line, plane are at the origin.

WlzShiftObj `WlzBoundingBox infile.wlz | \
       awk '{print "-x -" $1 " -y -" $2 " -z -" $3}'` \
       -o outfile.wlz infile.wlz




# The following commad also shifts the object to it's origin.
WlzShiftObj -g -o outfile.wlz infile.wlz

\endverbatim

\par File
\ref WlzShiftObj.c "WlzShiftObj.c"
\par See Also
\ref wlzboundingbox "WlzBoundingBox(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref WlzShiftObject "WlzShiftObject(3)"
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
		ok = 1,
		usage = 0,
		sftToOrg = 0,
		trX = 0,
		trY = 0,
		trZ = 0;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  WlzIBox3	bBox;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  char 		*outFileStr,
  		*inObjFileStr;
  static char	optList[] = "hgo:x:y:z:",
		outFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outFileStr = outFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'g':
        sftToOrg = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'x':
	if(sscanf(optarg, "%d", &trX) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'y':
	if(sscanf(optarg, "%d", &trY) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
        break;
      case 'z':
	if(sscanf(optarg, "%d", &trZ) != 1)
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
  if(ok)
  {
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
      if(fP && strcmp(inObjFileStr, "-"))
      {
	fclose(fP);
      }
    }
  }
  if(ok)
  {
    if(sftToOrg)
    {
      bBox = WlzBoundingBox3I(inObj, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        trX = -(bBox.xMin);
        trY = -(bBox.yMin);
        trZ = -(bBox.zMin);
      }
      else
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to get object origin\n",
		       *argv);
      }
    }
  }
  if(ok)
  {
    outObj = WlzAssignObject(WlzShiftObject(inObj, trX, trY, trZ,
					    &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to shift object\n",
		     *argv);
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"):
	      stdout)) == NULL) ||
       (WlzWriteObj(fP, outObj) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to write output object\n",
		     *argv);
    }
    if(fP && strcmp(outFileStr, "-"))
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
    " [-o<output object>] [-h] [-g] [-x#] [-y#] [-z#] [<input object>]\n"
    "Options:\n"
    "  -g  Shift object to the origin\n"
    "  -h  Prints this usage information\n"
    "  -o  Output object file name.\n"
    "  -x  Column (x) translation.\n"
    "  -y  Row (y) translation.\n"
    "  -z  Plane (z) translation.\n"
    "Shifts a Woolz object by applying an integer translation.\n"
    "The input object is read from stdin and the shifted object is\n"
    "written to stdout unless the filenames are given.\n",
    *argv,
    " -x100 -y200 -o shifted.wlz myobj.wlz\n"
    "The input Woolz object is read from myobj.wlz, shifted 100 columns\n"
    "and 200 lines and then written to shifted.wlz\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
