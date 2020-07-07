#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzFitBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzFitBSpline.c
* \author       Bill Hill
* \date         July 2020
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2020],
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
* \brief	Fits a B-spline to an object's domain.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzfitbspline "WlzFitBSpline"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>
/*!
\ingroup BinWlz
\defgroup wlzfitbspline WlzFitBSpline
\par Name
WlzFitBSpline
\par Synopsis
\verbatim
WlzFitBSpline [-h] [-c] [-k<order>] [-o<out object>] [-s<smoothing>]
                       [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Fit a closed (periodic) B-spline curve.</td>
  </tr>
  <tr>
    <td><b>-k</b></td>
    <td>Spline order.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Smoothing parameter (0.0 implies no smoothing).</td>
  </tr>
</table>
\par Description
Fits a B-spline to the domain of the given object and then outputs
a B-spline object.
By default the input object is read from stdin and the output object is
written to stdout.
\par Example
\verbatim
WlzFitBSpline -k3 -o spline.wlz -s 0.01 points.wlz
\endverbatim
Fits a B-spline to the ordered points read from points.wlz where the
spline has a smoothing factor of 0.01 and is of order 3. The output
B-spline object is written to the file spline.wlz.
\par File
\ref WlzFitBSpline.c "WlzFitBSpline.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineFromObj "WlzBSplineFromObj(3)"
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
		ok,
		order = 3,
		usage = 0,
		closed = 0;
  double	smoothing = 0.01;
  WlzDomain 	dom = {0};
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*ots;
  char		*iFile = NULL,
  		*oFile = NULL;
  static char	optList[] = "chk:o:s:",
		defFile[] = "-";
  iFile = oFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'c':
        closed = 1;
	break;
      case 'k':
         usage = (sscanf(optarg, "%d", &order) != 1) ||
	         (order < WLZ_BSPLINE_ORDER_MIN) ||
		 (order > WLZ_BSPLINE_ORDER_MAX);
        break;
      case 'o':
        oFile = optarg;
	break;
      case 's':
        usage = (sscanf(optarg, "%lg", &smoothing) != 1);
	break;
      case 'h':   /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(!usage && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      iFile = argv[optind];
    }
  }
  ok = !usage;
  if(ok)
  {
    if((fP = (strcmp(iFile, "-")?  fopen(iFile, "r"): stdin)) != NULL)
    {
      iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL);
    }
    if(fP && strcmp(iFile, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to read input object from  file %s (%s).\n",
	  *argv, iFile, ots);
    }
  }
  if(ok)
  {
    dom.bs = WlzBSplineFromObj(iObj, order, closed, smoothing, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to fit B-spline to input object (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    WlzValues nullVal = {0};

    oObj = WlzMakeMain(WLZ_SPLINE, dom, nullVal, NULL, NULL, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to create output object (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(oFile, "-")?  fopen(oFile, "w"): stdout)) != NULL)
    {
      errNum = WlzWriteObj(fP, oObj);
    }
    if(fP && strcmp(iFile, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
          "%s: Failed to write object to file %s (%s)\n",
	  *argv, oFile, ots);
    }
  }
  (void )WlzFreeObj(iObj);
  if(oObj)
  {
    (void )WlzFreeObj(oObj);
  }
  else
  {
    (void )WlzFreeDomain(dom);
  }
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%s%s",
    *argv,
    " [-h] [-c] [-k<order>] [-o<out object>] [-s<smoothing>] [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -c  Fit a closed (periodic) B-spline curve.\n"
    "  -k  Spline order.\n"
    "  -o  Output object file name.\n"
    "  -s  Smoothing parameter (0.0 implies no smoothing).\n"
    "Fits a B-spline to the domain of the given object and then outputs\n"
    "a B-spline object.\n"
    "By default the input object is read from stdin and the output object is\n"
    "written to stdout.\n"
    "Example:\n"
    "  ",
    *argv,
    " -k3 -o spline.wlz -s 0.01 points.wlz\n"
    "Fits a B-spline to the ordered points read from points.wlz where the\n"
    "spline has a smoothing factor of 0.01 and is of order 3. The output\n"
    "B-spline object is written to the file spline.wlz.\n");
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

