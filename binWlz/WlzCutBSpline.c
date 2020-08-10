#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCutBSpline_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCutBSpline.c
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
* \brief	Cuts regions from a spatial domain using a B-spline define
* 		the region cut.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcutbspline "WlzCutBSpline"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcutbspline WlzCutBSpline
\par Name
WlzCutBSpline
\par Synopsis
\verbatim
WlzCutBSpline [-h] [-b<spline object>] [-B] [-N] [-P] [-o<out object>]
              [-r <rad>] [-t[<min>],[<max>]] [-v] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>The given B-spline object.</td>
  </tr>
  <tr>
    <td><b>-B</b></td>
    <td>Cut a region of the B-spline dilated by given radius (default).</td>
  </tr>
  <tr>
    <td><b>-N</b></td>
    <td>Don't fill grey values if they exist.</td>
  </tr>
  <tr>
    <td><b>-P</b></td>
    <td>Cut planes orthogonal to the B-spline. These may be limited
        to within a given radius of the B-spline.</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Radius of region cut, the orthogonal distance from the B-spline.
        If cutting the region of the dilated B-spline then this is it's
	dilation, otherwise if cutting orthogonal planes this is the
	radius of the sections. If zero then either just the pixels/voxels
	intersecting the B-spline or the entire (unclipped) orthogonal
	planes will be cut (default 0).</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Parametric range along the B-spline curve (default 0,1).</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose output, possibly useful for tests and debugging.</td>
  </tr>
</table>
\par Description
Cuts regions from a spatial domain using a B-spline define the region
cut. The region may be the domain of the B-spline or planes orthogonal
to the B-spline.
By default the input objects are read from stdin and the output object
is written to stdout. If both the input and spline objects are read
from stdin then the input object is read before the spline object.
\par Example
\verbatim
WlzCutBSpline -P -t 0.5, -o out.wlz -r -b spline.wlz in.wlz
\endverbatim
Outputs a object containing sections (of radius 50) through the object
read from in.wlz which are orthogonal to the B-spline (read from
spline.wlz) and are lie within the parametric coordinate range
[0.5-1.0]. The output object is written to the file out.wlz
\par File
\ref WlzCutBSpline.c "WlzCutBSpline.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineCut "WlzBSplineCut(3)"
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
		radius = 0,
		noGrey = 0,
		cutOrthog = 0,
		usage = 0,
		verbose = 0;
  WlzObject	*iObj = NULL,
  		*oObj = NULL,
		*sObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*ots;
  char		*iFile = NULL,
  		*oFile = NULL,
		*sFile = NULL;
  double	param[2] = {0.0, 1.0};
  static char	optList[] = "hBNPvb:o:r:t:",
		defFile[] = "-";

  sFile = iFile = oFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
        sFile = optarg;
	break;
      case 'B':
        cutOrthog = 0;
	break;
      case 'N':
        noGrey = 1;
	break;
      case 'P':
        cutOrthog = 1;
	break;
      case 'r':
         usage = (sscanf(optarg, "%d", &radius) != 1) || (radius < 0);
        break;
      case 't':
	{
	  int	i;
	  char	*str;
	  char 	*cut[2] = {0};

	  str = WlzStringWhiteSpSkipLeading(optarg);
	  if(*str == ',')
	  {
	    cut[1] = strtok(str, ",");
	  }
	  else
	  {
	    cut[0] = strtok(str, ",");
	    cut[1] = strtok(NULL, ",");
	  }
	  if((cut[0] == NULL) && (cut[1] == NULL))
	  {
	    usage = 1;
	  }
	  else
	  {
	    for(i = 0; i < 2; ++i)
	    {
	      usage |= (cut[i] && (sscanf(cut[i], "%lg", param +i) != 1));
	    }
	    usage |= (param[0] < 0.0) || (param[1] > 1.0) ||
	             (param[0] > param[1]);
	  }
	}
        break;
      case 'o':
        oFile = optarg;
	break;
      case 'v':
        verbose = 1;
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
  if(verbose)
  {
    (void )fprintf(stderr,
        "%s: iFile = %s, oFile = %s, sFile = %s,\n"
	"\t\tcutOrthog = %d, radius = %d, param = [%g, %g],\n"
	"\t\tnoGrey = %d, verbose = %d, usage = %d\n",
	*argv, iFile, oFile, sFile, cutOrthog, radius, param[0], param[1],
	noGrey, verbose, usage);
  }
  ok = !usage;
  if(ok)
  {
    int		i;
    FILE	*fP = NULL;
    char 	*files[2] = {iFile, sFile};
    WlzObject	**objs[2] = {&iObj, &sObj};
    char	*str[2] = {"input", "spline"};

    for(i = 0; i < 2; ++i)
    {
      errNum = WLZ_ERR_FILE_OPEN;
      if((fP = (strcmp(files[i], "-")?  fopen(files[i], "r"): stdin)) != NULL)
      {
	*(objs[i]) = WlzAssignObject(WlzReadObj(fP, &errNum), NULL);
      }
      if(fP && strcmp(files[i], "-"))
      {
	(void )fclose(fP);
	fP = NULL;
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	ots = WlzStringFromErrorNum(errNum, NULL);
	(void )fprintf(stderr,
	    "%s: Failed to read %s object from  file %s (%s).\n",
	    *argv, str[i], files[i], ots);
        break;
      }
    }
  }
  if(ok)
  {
    WlzBSpline *bs = NULL;

    if(sObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(sObj->type != WLZ_SPLINE)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((bs = sObj->domain.bs) == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      oObj = WlzBSplineCut(iObj, bs, cutOrthog, noGrey, radius,
          param[0], param[1], &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to cut domain object (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

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
  (void )WlzFreeObj(oObj);
  (void )WlzFreeObj(sObj);
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s [-h] [-b<spline object>] [-B] [-N] [-P]\n"
    "\t\t[-o<out object>] [-r <rad>] [-t[<min>],[<max>]] [-v]\n"
    "\t\t[<in object>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -b  The given B-spline object.\n"
    "  -B  Cut region of the B-spline dilated by given radius (default).\n"
    "  -N  Don't fill grey values if they exist.\n"
    "  -P  Cut planes orthogonal to the B-spline. These may be limited\n"
    "      to within a given radius of the B-spline.\n"
    "  -r  Radius of region cut, the orthogonal distance from the B-spline.\n"
    "      If cutting the region of the dilated B-spline then this is it's\n"
    "      dilation, otherwise if cutting orthogonal planes this is the\n"
    "      radius of the sections. If zero then either just the pixels/voxels\n"
    "      intersecting the B-spline or the entire (unclipped) orthogonal\n"
    "      planes will be cut (set to %d).\n"
    "  -t  Parametric range along the B-spline curve (set to %g,%g).\n"
    "  -o  Output object file name.\n"
    "  -v  Verbose output, possibly useful for tests and debugging.\n"
    "Cuts regions from a spatial domain using a B-spline define the region\n"
    "cut. The region may be the domain of the B-spline or planes orthogonal\n"
    "to the B-spline.\n"
    "By default the input objects are read from stdin and the output object\n"
    "is written to stdout. If both the input and spline objects are read\n"
    "from stdin then the input object is read before the spline object.\n"
    "Example:\n"
    "  %s -P -t 0.5, -o out.wlz -r -b spline.wlz in.wlz\n"
    "Outputs a object containing sections (of radius 50) through the object\n"
    "read from in.wlz which are orthogonal to the B-spline (read from\n"
    "spline.wlz) and are lie within the parametric coordinate range\n"
    "[0.5-1.0]. The output object is written to the file out.wlz\n",
    *argv,
    WlzVersion(),
    radius,
    param[0], param[1],
    *argv );
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

