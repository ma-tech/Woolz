#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBSplineLength_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzBSplineLength.c
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
* \ref 		wlzbsplinelength "WlzBSplineLength"
*/


#include <string.h>
#include <stdio.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzbsplinelength WlzBSplineLength
\par Name
WlzBSplineLength
\par Synopsis
\verbatim
WlzBSplineLength [-h] [-t [<min>],[<max>]] [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Parametric range along the B-spline curve (default 0,1).</td>
  </tr>
</table>
\par Description
Computes the length of a path along a spline's curve between
a pair of parametric coordinates.
The length is written as a floating point number to the standard output.
By default the input object is read from stdin.
\par Example
\verbatim
WlzBSplineLength spline.wlz
\endverbatim
Prints out the length of the B-spline read from the file spline.wlz to the
standard output.
\par File
\ref WlzBSplineLength.c "WlzBSplineLength.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBSplineLength "WlzBSplineLength(3)"
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
  double	len = 0.0;
  WlzObject	*iObj = NULL;
  FILE		*fP = NULL;
  double        param[2] = {0.0, 1.0};
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*ots;
  char		*iFile = NULL;
  static char	optList[] = "ht:",
		defFile[] = "-";
  iFile = defFile;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 't':
        {
          int   i;
          char  *str;
          char  *cut[2] = {0};

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
    len = WlzBSplineLength(iObj->domain.bs, param[0], param[1], &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      ots = WlzStringFromErrorNum(errNum, NULL);
      (void )fprintf(stderr,
	  "%s: Failed to compute B-spline length (%s).\n",
	  *argv, ots);
    }
  }
  if(ok)
  {
    (void )printf("%lg\n", len);
  }
  (void )WlzFreeObj(iObj);
  if(usage)
  {
    errNum = WLZ_ERR_NONE;
    (void )fprintf(stderr,
    "Usage: %s%s%s%s%s%s",
    *argv,
    " [-h] [-t[<min>],[<max>]] [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -t  Parametric range along the B-spline curve (default 0,1).\n"
    "Computes the length of a path along a spline's curve between a pair\n"
    "of parametric coordinates.\n"
    "The length is written as a floating point number to the standard output.\n"
    "By default the input object is read from stdin.\n"
    "Example:\n"
    "  ",
    *argv,
    " spline.wlz\n"
    "Prints out the length of the B-spline read from the file spline.wlz to\n"
    "the standard output.\n");
  }
  exit(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

