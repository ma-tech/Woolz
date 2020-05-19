#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzLineSkeletonSegment_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzLineSkeletonSegment.c
* \author       Bill Hill
* \date         July 2016
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2016],
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
* \brief	Computes the line skeleton of the given object's domain.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzlineskeleton "WlzLineSkeletonSegment"
*/

#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzlineskeleton WlzLineSkeletonSegment
\par Name
WlzLineSkeletonSegment



    "Usage: %s [-h] [-p #,#,#] [-q #,#,#] [-v] [<input object>]\n"
    "Computes a line skeleton segment of the given object's domain.\n"
    "It is an error if both given points are not in the objects \n"
    "domain.\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -v  Verbose output messages.\n"
    "  -o  Output file.\n"
    "  -p  First point of line segment.\n"
    "  -q  Second point of line segment.\n"
    "By default the input object is read from stdin and the output\n"
    "object is written to stdout.\n"
    "Example:\n"
    "  %s -o path.wlz -p 249,325,37 -q 116,241,302 colon.wlz\n"
    "Computes a path from the point 249,325,37 to the point 116,241,302\n"
    "in the domain of the object read from the file colon.wlz. The output\n"
    "path is written to the file path.wlz",


\par Synopsis
\verbatim
WlzLineSkeletonSegment [-h] [-v] [-p #,#,#] [-q #,#,#] [-o<output object>]
                       [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>First point of line segment.</td>
  </tr>
  <tr>
    <td><b>-q</b></td>
    <td>Second point of line segment.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Verbose output messages.</td>
  </tr>
</table>
\par Description
Computes a line skeleton segment of the given object's domain.
It is an error if both given points are not in the objects domain.
By default the input object is read from stdin and the output object is written to stdout.
\par Examples
\verbatim
WlzLineSkeletonSegment -o path.wlz -p 249,325,37 -q 116,241,302 colon.wlz
\endverbatim
Computes a path from the point 249,325,37 to the point 116,241,302
in the domain of the object read from the file colon.wlz. The output
path is written to the file path.wlz
\par File
\ref WlzLineSkeletonSegment.c "WlzLineSkeletonSegment.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzLineSkeletonSegment "WlzLineSkeletonSegment(3)"
\ref WlzLineSkeleton "WlzLineSkeleton(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
  		ok = 1,
  		usage = 0,
		verbose = 0;
  FILE		*fP = NULL;
  char		*inObjStr,
  		*outObjStr;
  WlzIVertex3	p0,
  		p1;
  struct timeval times[3];
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "hvo:p:q:";
  const char    inObjStrDef[] = "-",
  	        outObjStrDef[] = "-";

  opterr = 0;
  WLZ_VTX_3_ZERO(p0);
  WLZ_VTX_3_ZERO(p1);
  inObjStr = (char *)inObjStrDef;
  outObjStr = (char *)outObjStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjStr = optarg;
	break;
      case 'p':
        if(sscanf(optarg, "%d,%d,%d", &(p0.vtX), &(p0.vtY), &(p0.vtZ)) < 2)
        {
	  usage = 1;
	}
	break;
      case 'q':
        if(sscanf(optarg, "%d,%d,%d", &(p1.vtX), &(p1.vtY), &(p1.vtZ)) < 2)
        {
	  usage = 1;
	}
	break;
      case 'v':
        verbose = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if((inObjStr == NULL) || (*inObjStr == '\0') ||
       (outObjStr == NULL) || (*outObjStr == '\0'))
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
        inObjStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjStr == NULL) ||
       (*inObjStr == '\0') ||
       ((fP = (strcmp(inObjStr, "-")?
              fopen(inObjStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjStr);
    }
    if(fP && strcmp(inObjStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    if(verbose)
    {
      gettimeofday(times + 0, NULL);
    }
    outObj = WlzLineSkeletonSegment(inObj, inObj->type, p0, p1, &errNum);
    if(verbose)
    {
      gettimeofday(times + 1, NULL);
      ALC_TIMERSUB(times + 1, times + 0, times + 2);
      (void )fprintf(stderr,
                     "%s: Elapsed time for WlzLineSkeletonSegment() %gus\n",
                     argv[0], (1.0e06 * times[2].tv_sec) + times[2].tv_usec);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to compute line skeleton, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outObjStr, "-")?  fopen(outObjStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outObjStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write object to file %s\n",
                     *argv, outObjStr);
    }
  }
  if(fP && strcmp(outObjStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-p #,#,#] [-q #,#,#] [-v] [<input object>]\n"
    "Computes a line skeleton segment of the given object's domain.\n"
    "It is an error if both given points are not in the objects \n"
    "domain.\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -v  Verbose output messages.\n"
    "  -o  Output file.\n"
    "  -p  First point of line segment.\n"
    "  -q  Second point of line segment.\n"
    "By default the input object is read from stdin and the output\n"
    "object is written to stdout.\n"
    "Example:\n"
    "  %s -o path.wlz -p 249,325,37 -q 116,241,302 colon.wlz\n"
    "Computes a path from the point 249,325,37 to the point 116,241,302\n"
    "in the domain of the object read from the file colon.wlz. The output\n"
    "path is written to the file path.wlz\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
