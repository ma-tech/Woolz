#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPointsToText_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzPointsToText.c
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
* \brief	Writes out a Woolz points object as text.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzpointstotext	 "WlzPointsToText"
*/

/*!
\ingroup BinWlz
\defgroup wlzmarkerstodomain WlzPointsToText
\par Name
WlzPointsToText - writes out a Woolz points object as text.
\par Synopsis
\verbatim
WlzPointsToText [-h] [-n] [-o<output file>] [-f <fs>] [-r <rs>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Number of points as first record.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name, default standard output.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Field seperator (default single space).</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Record seperator (default newline).</td>
  </tr>
</table>
\par Description
Reads a 2D or 3D points object from either a given file or the
standard input (default) and then writes out the points as text.
\par Examples
\verbatim
WlzPointsToText -o out.num in.wlz
\endverbatim
Writes te point locations from the points object in file in.wlz to the
output file text out.num.
\par File
\ref WlzPointsToText.c "WlzPointsToText.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzPointsFromDomain "WlzPointsFromDomain(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <Wlz.h>

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		outNPts = 0,
  		usage = 0;
  FILE		*fP = NULL;
  char		*ptsFile,
  		*outFile,
		*fldSep,
		*recSep;
  const char	*errMsg;
  WlzObject	*ptsObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "hno:f:r:";
  const char    fileDef[] = "-",
		fldSepDef[] = " ",
  		recSepDef[] = "\n";

  opterr = 0;
  ptsFile = (char *)fileDef;
  outFile = (char *)fileDef;
  fldSep  = (char *)fldSepDef;
  recSep  = (char *)recSepDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outFile = optarg;
	break;
      case 'f':
	if((fldSep = WlzStringUnescape(optarg)) == NULL)
	{
	  usage = 1;
	}
        break;
      case 'n':
        outNPts = 1;
	break;
      case 'r':
	if((recSep = WlzStringUnescape(optarg)) == NULL)
	{
	  usage = 1;
	}
        break;
      case 'h': /* FALLTROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    if((outFile == NULL) || (*outFile == '\0') ||
       (ptsFile == NULL) || (*ptsFile == '\0'))
    {
      usage = 1;
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
      ptsFile = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    if((ptsFile == NULL) ||
       (*ptsFile == '\0') ||
       ((fP = (strcmp(ptsFile, "-")?
              fopen(ptsFile, "r"): stdin)) == NULL) ||
       ((ptsObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_READ_EOF;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read input file %s (%s).\n",
		     *argv, ptsFile, errMsg);
    }
  }
  if(fP)
  {
    if(strcmp(ptsFile, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(ok)
  {
    if(ptsObj->type != WLZ_POINTS)
    {
      ok = 0;
      errNum = WLZ_ERR_OBJECT_TYPE;
      (void )fprintf(stderr,
	             "%s: Object read from input file %s is not a points "
		     "object.\n",
		     *argv, ptsFile);
    }
    else if(ptsObj->domain.core == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_DOMAIN_NULL;
      (void )fprintf(stderr,
	             "%s: Object read from input file %s has NULL domain.",
		     *argv, ptsFile);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFile, "-")? fopen(outFile, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFile);
    }
  }
  if(ok)
  {
    int		i;
    WlzPoints	*pts;

    pts = ptsObj->domain.pts;
    if(outNPts)
    {
      (void )fprintf(fP, "%d%s", pts->nPoints, recSep);
    }
    switch(pts->type)
    {
      case WLZ_POINTS_2I:
	for(i = 0; i < pts->nPoints; ++i)
	{
	  WlzIVertex2 *p;

	  p = pts->points.i2 + i;
	  (void )fprintf(fP, "%d%s%d%s", p->vtX, fldSep, p->vtY, recSep);
	}
	break;
      case WLZ_POINTS_2D:
	for(i = 0; i < pts->nPoints; ++i)
	{
	  WlzDVertex2 *p;

	  p = pts->points.d2 + i;
	  (void )fprintf(fP, "%lg%s%lg%s", p->vtX, fldSep, p->vtY, recSep);
	}
	break;
      case WLZ_POINTS_3I:
	for(i = 0; i < pts->nPoints; ++i)
	{
	  WlzIVertex3 *p;

	  p = pts->points.i3 + i;
	  (void )fprintf(fP, "%d%s%d%s%d%s",
	                 p->vtX, fldSep, p->vtY, fldSep, p->vtZ, recSep);
	}
	break;
      case WLZ_POINTS_3D:
	for(i = 0; i < pts->nPoints; ++i)
	{
	  WlzDVertex3 *p;

	  p = pts->points.d3 + i;
	  (void )fprintf(fP, "%lg%s%lg%s%lg%s",
	                 p->vtX, fldSep, p->vtY, fldSep, p->vtZ, recSep);
	}
	break;
      default:
        break;
    }
    if(feof(fP))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to write to output file %s.\n",
		     argv[0], outFile);
    }
  }
  if(fP && strcmp(outFile, "-"))
  {
    (void )fclose(fP);
  }
  (void )WlzFreeObj(ptsObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-n] [-o<output file>] [-f <fs>] [-r <rs>]\n"
    "\t\t[<input file>]\n"
    "Version: %s\n"
    "Options:\n"
    "  -h  Output this usage message.\n"
    "  -n  Number of points as first record.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -f  Field seperator (default single space).\n"
    "  -r  Record seperator (default newline).\n"
    "Reads a 2D or 3D points object from either a given file or the\n"
    "standard input (default) and then writes out the points as text.\n"
    "standard input (default). The vertices of the points are then used\n"
    "Example:\n"
    "  %s -o out.num in.wlz\n"
    "Writes te point locations from the points object in file in.wlz to the\n"
    "output file text out.num.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
