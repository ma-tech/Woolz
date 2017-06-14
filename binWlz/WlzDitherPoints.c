#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDitherPoints_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzDitherPoints.c
* \author       Bill Hill
* \date         April 2003
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
* \brief	Creates a domain with a marker located at the position
* 		of each point vertex using a points object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzditherpoints "WlzDitherPoints"
*/

/*!
\ingroup BinWlz
\defgroup wlzditherpoints  WlzDitherPoints
\par Name
WlzDitherPoints - dithers point positions.
\par Synopsis
\verbatim
WlzDitherPoints [-h] [-o<output file>] [-D #,#[#]] [-r <res obj>]
                [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr> <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr> <tr> 
    <td><b>-D</b></td>
    <td>Dither the points by applying a random offset. Supplied values
        are the maximum dither displacement. Default value is (1,1,1).</td>
  </tr> <tr> 
    <td><b>-r</b></td>
    <td>Restriction object.</td>
  </tr>
</table>
\par Description
Reads a 2D or 3D points object from either a given file or the
standard input (default). The vertices of the points are then
converted to floating point, dithered and output.
A restriction object may be given to restrict all dithered point
positions to within it's domain, but if given this must be a spatial
domain object of the appropriate dimension. Restriction assumes voxel
spacing of (1, 1, 1) and may not be appropriate.
When reading from the standard input the input object is read before
an optional restriction object.
\par Examples
\verbatim
WlzDitherPoints -o out.wlz -D 1,1,1 in.wlz
\endverbatim
Creates a new points file where the point vertices are double precision and
dithered by up to a voxel in each direction.
The input points are read from in.wlz and the output points written to out.wlz.
\par File
\ref binWlz/WlzDitherPoints.c "WlzDitherPoints.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzPointsDither "WlzPointsDither(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <sys/time.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		usage = 0;
  char		*inFileStr,
		*outFileStr,
		*resFileStr = NULL;
  FILE		*fP = NULL;
  WlzDVertex3	ditherSz = {1.0, 1.0, 1.0};
  WlzObject	*inObj = NULL,
		*outObj = NULL,
		*resObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "ho:D:r:";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'D':
	if(sscanf(optarg, "%lg,%lg,%lg",
	          &(ditherSz.vtX), &(ditherSz.vtY), &(ditherSz.vtZ)) < 2)
	{
	  usage = 1;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'r':
        resFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(ok)
  {
    ok = !usage;
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
      inFileStr = argv[optind];
    }
  }
  /* Read the input object. */
  if(ok)
  {
    if((*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read input object from file %s.\n",
                     argv[0], inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read restriction object if required. */
  if(ok && resFileStr)
  {
    if((*resFileStr == '\0') ||
       ((fP = (strcmp(resFileStr, "-")?
              fopen(resFileStr, "r"): stdin)) == NULL) ||
       ((resObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read restriction object from file %s.\n",
                     argv[0], resFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    WlzDomain dDom;

    dDom.pts = WlzPointsDither(inObj->domain.pts, ditherSz, resObj, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzMakeMain(inObj->type, dDom, inObj->values, inObj->plist,
      			   NULL, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to dither points %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Output the points object. */
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(resObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-D#,#[,#]] [-o <output file>] [-r <res obj>] [-h]\n"
    "\t\t[<input file>]\n"
    "Dithers point positions. A 2D or 3D points object is read either from\n"
    "a given file or the standard input (default). The vertices of the\n"
    "points are then converted to floating point, dithered and output\n"
    "either to the given file or the standard output (default).\n"
    "A restriction object may be given to restrict all dithered point\n"
    "positions to within it's domain, but if given this must be a spatial\n"
    "domain object of the appropriate dimension. Restriction assumes voxel\n"
    "spacing of (1, 1, 1) and may not be appropriate.\n"
    "When reading from the standard input the input object is read before\n"
    "an optional restriction object.\n"
    "Version: %s\n"
    "Options:\n"
    "  -D  Dither the points by applying a random offset. Supplied values\n"
    "      are the maximum dither displacement. Default value (1, 1, 1).\n"
    "  -o  Output object file.\n"
    "  -r  Restriction object.\n"
    "  -h  Help - prints this usage message\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o out.wlz -D 1,1,1 in.wlz\n"
    "The input points are read from in.wlz and points dithered using a\n"
    "maximum displacement of 1,1,1 are written to out.wlz.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
