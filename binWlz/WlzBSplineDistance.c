#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBSplineDistance_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzBSplineDistance.c
* \author       Bill Hill
* \date         December 2021
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
* \brief	Computes the parametric coordinate at a given distance from a
* 		reference parametric coordinate in a B-spline curve.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzbsplinedistance "WlzBSplineDistance"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzbsplinedistance WlzBSplineDistance
\par Name
WlzBSplineDistance
\par Synopsis
\verbatim
WlzBSplineDistance [-h]  -d<dst> [-e<tol>] [-s <x>,<y>[,<z>]] [-t <ref>]
                   [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Distance from reference parametric coordinate (required).</td>
  </tr>
  <tr>
    <td><b>-e</b></td>
    <td>Tolerance for distances (default 0.01).</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Step size (default 1.0,1.0,1.0).</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Reference parametric coordinate along the B-spline curve
        (default 0.0).</td>
  </tr>
</table>
\par Description
Computes the parametric coordinate at a given distance from a reference
parametric coordinate.
The step size may be used to compute the spline distances for non-unit
pixel/poxel sizes.
The parametric coordinate is written as a floating point number to the
standard output.
By default the input object is read from stdin.
\par Example
\verbatim
WlzBSplineDistance -d 20.0 spline.wlz
\endverbatim
Prints out the parametric coordinate at the distance 20.0 from the start
of the B-spline read from the file spline.wlz to the standard output.
\par File
\ref WlzBSplineDistance.c "WlzBSplineDistance.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineDistance "WlzBSplineDistance(3)"
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
		usage = 0;
  double	t = 0.0,
  		d = NAN,
		tol = 0.01;
  WlzObject	*iObj = NULL;
  FILE		*fP = NULL;
  double	stepSz[3] = {1.0, 1.0, 1.0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char    *ots;
  char		*iFile = NULL;
  static char	optList[] = "hd:e:s:t:",
		defFile[] = "-";
  iFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
        if(sscanf(optarg, "%lg", &d) < 1)
	{
	  usage = 1;
	}
	break;
      case 'e':
        if(sscanf(optarg, "%lg", &tol) < 1)
	{
	  usage = 1;
	}
	break;
      case 's':
        if(sscanf(optarg, "%lg,%lg,%lg",
	          stepSz + 0, stepSz + 1, stepSz + 2) < 2)
        {
	  usage = 1;
	}
	break;
      case 't':
        if(sscanf(optarg, "%lg", &t) < 1)
	{
	  usage = 1;
	}
	break;
      case 'h':   /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(!usage && isnan(d))
  {
    usage = 1;
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
    if(iObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(iObj->type != WLZ_SPLINE)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(iObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
          "%s: Input object is inappropriate (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    t = WlzBSplineDistance(iObj->domain.bs, t, d,
                           stepSz[2], stepSz[1], stepSz[0], tol,  &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to compute B-spline parametric coordinate (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    (void )printf("%lg\n", t);
  }
  (void )WlzFreeObj(iObj);
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%s%s",
    *argv,
    " [-h] -d<dst> [-e<tol>] [-s<x>,<y>[,<z>]] [-t<ref>]\n"
    "\t\t[<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -d  Distance from reference parametric coordinate (required).\n"
    "  -e  Tolerance for distances (default 0.01).\n"
    "  -s  Step size (default 1.0,1.0,1.0).\n"
    "  -t  Reference parametric coordinate along the B-spline curve"
    "      (default 0.0).\n"
    "Computes the parametric coordinate at a given distance from a reference\n"
    "parametric coordinate.\n"
    "The step size may be used to compute the spline distances for non-unit\n"
    "pixel/poxel sizes.\n"
    "The parametric coordinate is written as a floating point number to the\n"
    "standard output.\n"
    "By default the input object is read from stdin.\n"
    "Example:\n"
    "  ",
    *argv,
    " -d 20.0 spline.wlz\n"
    "Prints out the parametric coordinate at the distance 20.0 from the start\n"
    "of the B-spline read from the file spline.wlz to the standard output.\n");
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

