#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshExtrapolate_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshExtrapolate.c
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
* \brief	Extrapolate values within a conforming mesh.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshextrapolate "WlzCMeshExtrapolate"
*/

#include <stdio.h>
#include <float.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshextrapolate WlzCMeshExtrapolate
\par Name
WlzCMeshExtrapolate - Extrapolates values within a conforming mesh.
\par Synopsis
\verbatim
WlzCMeshExtrapolate [-h] [-N] [-o<output object>] [-u<low>,<high>]
                    [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Use nearest neighbour rather than linear extrapolation.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Unknown value range (default is NaN).</td>
  </tr>
</table>
\par Description
Extrapolates values within a conforming mesh.
\par Examples
\verbatim
WlzCMeshExtrapolate -o allknownmesh.wlz partknownmesh.wlz
\endverbatim
Creates a new conforming mesh object that has all values known using the
given mesh and it's values to extrapolate unknown values.
\par File
\ref WlzCMeshExtrapolate.c "WlzCMeshExtrapolate.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
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
  		usage = 0;
  double	uvRange[2] = {NAN};
  WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;
  FILE		*fP = NULL;
  char		*o,
  		*inObjStr,
  		*outObjStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzUByte	*unk = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "hNo:u:";
  const char    inObjStrDef[] = "-",
  	        outObjStrDef[] = "-";

  opterr = 0;
  inObjStr = (char *)inObjStrDef;
  outObjStr = (char *)outObjStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'N':
        interp = WLZ_INTERPOLATION_NEAREST;
	break;
      case 'o':
        outObjStr = optarg;
	break;
      case 'u':
	for(*o = optarg; *o; ++o)
	{
	  *o = tolower(*o);
	}
        if(strstr(optarg, "nan") == NULL)
	{
	  if(sscanf(optarg, "%lg,%lg", &(uvRange[0]), &(uvRange[1])) != 2)
	  {
	    usage = 1;
	  }
	}
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
    unk = WlzCMeshIndexMaskFromValueRange(inObj, uvRange[0], uvRange[1],
    				          1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = WlzCMeshExpValues(inObj, unk, interp, &errNum);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to create conforming mesh, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  AlcFree(unk);
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
    fprintf(stderr,
            "Usage: %s [-h] [-N] [-o<output object>] [-u<low>,<high>]\n"
	    "\t\t[<input object>]\n"
            "Creates a new conforming mesh object that has all values known\n"
	    "using the given mesh and it's values to extrapolate unknown\n"
	    "values."
	    "Version: %s\n"
	    "Options:\n"
	    "  -h  Help, prints this usage message.\n"
            "  -N  Use nearest neighbour rather than linear extrapolation.\n"
	    "  -o  Output file.\n"
	    "  -u  Unknown value range (default is NaN).\n",
	    argv[0],
	    WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
