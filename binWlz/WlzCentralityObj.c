#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCentralityObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCentralityObj.c
* \author       Bill Hill
* \date         April 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Computes the centrality of a feature object's domain
* 		with respect to a reference object's domain.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
* \par Binary
* \ref wlzcentralityobj "WlzCentralityObj"
*/

/*!
\ingroup BinWlz
\defgroup wlzcentralityobj WlzCentralityObj
\par Name
WlzCentralityObj - computes the centrality of a feature object's domain
		   with respect to a reference object's domain.
\par Synopsis
\verbatim
WlzCentralityObj [-b] [-h] [-n#] [-o <output file>] [-R]
                 [<feat obj>] [<ref obj>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-b</b></td>
    <td>Treat the feature object as a binary object with constant values.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help - prints a usage message.</td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Number of angle increments to use (default 360).</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Output maximum radius after centrality value.</td>
  </tr>
</table>
\par Description
Computes the centrality of the feature object's domain with respect to
the reference object's domain.
By default the feature and reference objects are read from the standard
input. This may be made explicit if only one of these is given by using
- to specify the standard input.
\par Examples
\verbatim
WlzCentralityObj -b spots.wlz fish.wlz
\endverbatim
Computes the centrality of the domain spots.wlz with respect to the reference
domain fish.wlz. The sports domain is treated as a binary domain (values are
ignored).
\par File
\ref WlzCentralityObj.c "WlzCentralityObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCentrality "WlzCentrality(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok,
		binFlg = 0,
  		nAng = 360,
		maxRFlg = 0,
  		option,
		usage = 0;
  double	cen = 0.0,
  		maxR = 0.0;
  char		*fFileStr,
  		*rFileStr,
		*oFileStr;
  FILE		*fP = NULL;
  WlzObject	*fObj = NULL,
		*rObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "bhn:o:R";
  const char    fileStrDef[] = "-";

  /* Parse the argument list and check for input files. */
  opterr = 0;
  fFileStr = (char *)fileStrDef;
  rFileStr = (char *)fileStrDef;
  oFileStr = (char *)fileStrDef;
  while((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'b':
        binFlg = 1;
	break;
      case 'n':
	if((sscanf(optarg, "%d", &nAng) != 1) || (nAng <= 0))
	{
	  usage = 1;
	}
	break;
      case 'o':
	oFileStr = optarg;
	break;
      case 'R':
        maxRFlg = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok && (optind < argc))
  {
    if((optind + 1) == argc)
    {
      fFileStr = argv[optind];
    }
    else if((optind + 2) == argc)
    {
      fFileStr = argv[optind];
      rFileStr = argv[optind + 1];
    }
    else
    {
      usage = 1;
    }
  }
  ok = !usage;
  /* Read the feature object. */
  if(ok)
  {
    if((fFileStr == NULL) ||
       (*fFileStr == '\0') ||
       ((fP = (strcmp(fFileStr, "-")?
              fopen(fFileStr, "r"): stdin)) == NULL) ||
       ((fObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read feature object from file %s.\n",
                     argv[0], fFileStr);
    }
    if(fP && strcmp(fFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the reference object. */
  if(ok)
  {
    if((rFileStr == NULL) ||
       (*rFileStr == '\0') ||
       ((fP = (strcmp(rFileStr, "-")?
              fopen(rFileStr, "r"): stdin)) == NULL) ||
       ((rObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read reference object from file %s.\n",
                     argv[0], rFileStr);
    }
    if(fP && strcmp(rFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Check object types and parse the promote the default distance function
   * if required. */
  if(ok)
  {
    if(rObj->type != fObj->type)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Feature and reference object types differ.\n",
		     argv[0]);
    }
  }
  /* Compute centrality feature. */
  if(ok)
  {
    cen = WlzCentrality(fObj, rObj, nAng, binFlg, &maxR, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to compute centrality feature value (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  /* Output the centrality feature value. */
  if(ok)
  {
    if((fP = (strcmp(oFileStr, "-")?
             fopen(oFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], oFileStr);
    }
  }
  if(ok)
  {
    if(((maxRFlg == 0) && (fprintf(fP, "%lg\n", cen) <= 0)) ||
       ((maxRFlg != 0) && (fprintf(fP, "%lg %lg\n", cen, maxR) <= 0)))
    {
      ok = 0;
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output feature value (%s).\n",
                     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(fObj);
  (void )WlzFreeObj(rObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-b] [-h] [-n#] [-o <output file>] [-R]\n"
    "Usage:                  [<feat obj>] [<ref obj>]\n"
    "Computes the centrality of the feature object's domain with respect to\n"
    "the reference object's domain.\n"
    "Options are:\n"
    "  -b  Treat the feature object as a binary object with constant values.\n"
    "  -h  Help - prints this usage message\n"
    "  -n  Number of angle increments to use (default 360).\n"
    "  -o  Output file for centrality feature value.\n"
    "  -R  Output maximum radius value after centrality value.\n"
    "By default the feature and reference objects are read from the standard\n"
    "input. This may be made explicit if only one of these is given by using\n"
    "- to specify the standard input.\n",
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
