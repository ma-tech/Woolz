#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBSplineToDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzBSplineToDomain.c
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
* \brief	Creates an interval or plane domain object corresponding to
* 		a B-spline domain.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzbsplinetodomain "WlzBSplineToDomain"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzbsplinetodomain WlzBSplineToDomain
\par Name
WlzBSplineToDomain
\par Synopsis
\verbatim
WlzBSplineToDomain [-h] [-V] [-o<out object>] [-t[<min>],[<max>]] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-V</b></td>
    <td>The output object will have values set to the B-spline parametric
        value.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Parametric range along the B-spline curve (default 0,1).</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
Creates an interval or plane domain object in which the pixels/voxels
intersect the given B-spline object. The output object may have
disconnected pixels/voxels unless dilated by 1.
By default the input object is from stdin and the output object is written
to stdout.
\par Example
\verbatim
WlzBSplineToDomain -t 0.5, -o out.wlz spline.wlz
\endverbatim
Outputs a domain object corresponding to the pixels/voxels intersected
by the given B-spline object in it's parametric range 0.5-1.0.
The B-spline object is read from spline.wlz and the output object is
written to the file out.wlz.
\par File
\ref WlzBSplineToDomain.c "WlzBSplineToDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineToDomain "WlzBSplineToDomain(3)"
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
		paramVal = 0,
		usage = 0;
  WlzObject	*iObj = NULL,
  		*oObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*ots;
  char		*iFile = NULL,
  		*oFile = NULL;
  double	param[2] = {0.0, 1.0};
  static char	optList[] = "hVo:t:",
		defFile[] = "-";

  iFile = oFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'V':
        paramVal = 1;
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
    if((iFile == NULL) || (*iFile == '\0') ||
       (oFile == NULL) || (*oFile == '\0'))
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
        iFile = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    FILE	*fP;

    if((iFile == NULL) ||
       (*iFile == '\0') ||
       ((fP = (strcmp(iFile, "-")?
              fopen(iFile, "r"): stdin)) == NULL) ||
       ((iObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, iFile);
    }
    if(fP && strcmp(iFile, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    WlzBSpline *bs = NULL;

    if(iObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(iObj->type != WLZ_SPLINE)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if((bs = iObj->domain.bs) == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      oObj = WlzBSplineToDomain(bs, param[0], param[1], paramVal, &errNum);
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
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s [-h] [-V] [-o<out object>] [-t[<min>],[<max>]] [<in object>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -V  The output object will have values set to the B-spline parametric\n"
    "      value.\n"
    "  -t  Parametric range along the B-spline curve (set to %g,%g).\n"
    "  -o  Output object file name.\n"
    "Creates an interval or plane domain object in which the pixels/voxels\n"
    "intersect the given B-spline object. The output object may have\n"
    "disconnected pixels/voxels unless dilated by 1.\n"
    "By default the input object is read from stdin and the output object\n"
    "is written to stdout.\n"
    "Example:\n"
    "  %s -t 0.5, -o out.wlz spline.wlz\n"
    "Outputs a domain object corresponding to the pixels/voxels intersected\n"
    "by the given B-spline object in it's parametric range [0.5-1.0].\n"
    "The B-spline object is read from spline.wlz and the output object is\n"
    "written to the file out.wlz\n",
    *argv,
    WlzVersion(),
    param[0], param[1],
    *argv );
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

