#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzEvalBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzEvalBSpline.c
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
* \brief	Evaluates a B-spline to give a spatial domain or points
* 		object.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzevalbspline "WlzEvalBSpline"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzevalbspline WlzEvalBSpline
\par Name
WlzEvalBSpline
\par Synopsis
\verbatim
WlzEvalBSpline [-h] [-P] [-d<deriv>] [-n<evaluations>] [-o<out object>]
               [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-P</b></td>
    <td>Output a points object instead of a pixel / voxel spatial domain
        object.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Order of derivative, range [0-5]./td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Number of equally spaced evaluations of the B-spline within the
        parametrised curve domain. If zero then evaluations are made at
	the knots./td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
Evaluates a B-spline and then outputs either a pixel / voxel object or
a points object.
By default the input object is read from stdin and the output object is
written to stdout.
\par Example
\verbatim
WlzEvalBSpline -n 1000 -o out.wlz spline.wlz
\endverbatim
Evaluates the B-spline read from the file spline.wlz at 1000 points
spaced equally along the parametrised curve and then writes the
pixel / voxel spatial domain object to the file out.wlz.
\par File
\ref WlzEvalBSpline.c "WlzEvalBSpline.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineEval "WlzBSplineEval(3)"
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
		deriv = 0,
		nEval = 1000,
		usage = 0;
  WlzPoints	*pts = NULL;
  WlzObjectType outObjType = WLZ_NULL;
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*ots;
  char		*iFile = NULL,
  		*oFile = NULL;
  static char	optList[] = "hPd:n:o:",
		defFile[] = "-";
  iFile = oFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'P':
        outObjType = WLZ_POINTS;
	break;
      case 'd':
         usage = (sscanf(optarg, "%d", &deriv) != 1) ||
	         (deriv < 0) || (deriv > WLZ_BSPLINE_ORDER_MAX);
        break;
      case 'n':
         usage = (sscanf(optarg, "%d", &nEval) != 1) || (nEval < 0);
        break;
      case 'o':
        oFile = optarg;
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
    else
    {
      switch(iObj->domain.core->type)
      {
        case WLZ_BSPLINE_C2D:
          if(outObjType == WLZ_NULL)
	  {
	    outObjType = WLZ_2D_DOMAINOBJ;
	  }
	  break;
        case WLZ_BSPLINE_C3D:
          if(outObjType == WLZ_NULL)
	  {
	    outObjType = WLZ_3D_DOMAINOBJ;
	  }
	  break;
	default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  break;
      }
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
    pts = WlzBSplineEvalPoints(iObj->domain.bs, nEval, deriv, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to evaluate B-spline (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    if(outObjType == WLZ_POINTS)
    {
      WlzDomain	dom;
      WlzValues nullVal = {0};

      dom.pts = pts;
      oObj = WlzMakeMain(WLZ_POINTS, dom, nullVal, NULL, NULL, &errNum);
    }
    else
    {
      oObj = WlzPointsToDomObj(pts, 0.0, &errNum);
    }
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
    WlzDomain	dom;

    dom.pts = pts;
    (void )WlzFreeDomain(dom);
  }
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%s%s",
    *argv,
    " [-h] [-P] [-d<deriv>] [-n<evaluations>] [-o<out object>]\n"
    "\t\t[<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -P  Output a points object instead of a pixel / voxel spatial domain\n"
    "      object.\n"
    "  -d  Order of derivative, range [0-5].\n"
    "  -n  Number of equally spaced evaluations of the B-spline within the\n"
    "      parametrised curve domain. If zero then evaluations are made at\n"
    "      the knots.\n"
    "  -o  Output object file name.\n"
    "Evaluates a B-spline and then outputs either a pixel / voxel object or\n"
    "a points object.\n"
    "By default the input object is read from stdin and the output object is\n"
    "written to stdout.\n"
    "Example:\n"
    "  ",
    *argv,
    " -n 1000 -o out.wlz spline.wlz\n"
    "Evaluates the B-spline read from the file spline.wlz at 1000 points\n"
    "spaced equally along the parametrised curve and then writes the\n"
    "pixel / voxel spatial domain object to the file out.wlz.\n");
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

