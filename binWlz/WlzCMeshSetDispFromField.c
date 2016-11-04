#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshSetDispFromField_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshSetDispFromField.c
* \author       Bill Hill
* \date         September 2016
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
* \brief	Creates a conforming mesh object with indexed
* 		values set using the given displacement field
* 		object.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcmeshsetdispfromfield "WlzCMeshSetDispFromField"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshsetdispfromfield WlzCMeshSetDispFromField
\par Name
WlzCMeshSetDispFromField - Sets mesh displacement values from field object.
\par Synopsis
\verbatim
WlzCMeshSetDispFromField [-h] [-A] [-L] [-N] [-R] [-o<output object>]
			 [-b #,#[,#]] [-m<mesh object>] [<field object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-A</b></td>
    <td>Displacement field values are absolute (default).</td>
  </tr> 
  <tr> 
    <td><b>-N</b></td>
    <td>Use nearest neighbour rather than linear extrapolation.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use linear interpolation.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-R</b></td>
    <td>Displacement field values are relative to position.</td>
  </tr> 
  <tr>
    <td><b>-b</b></td>
    <td>Background displacement value (default 0,0[,0]).</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Mesh object.</td>
  </tr>
</table>
\par Description
Creates a conforming mesh object that has values set from
the distance field object. Mesh values outside the distance
field object will be set to the background value.
The given field object must be a compound array object with
the appropriate number of components for the mesh dimension.
By default the mesh and field objects are read from the
standard input, with the mesh object being read before the
field object.
\par Examples
\verbatim
WlzCMeshSetDispFromField -o out.wlz -m mesh field.wlz
\endverbatim
Creates a conforming mesh object (out.wlz) that has all values set by
linear interpolation from the given field object (field.wlz) at the
nodes of the given mesh (mesh.wlz).
\par File
\ref WlzCMeshSetDispFromField.c "WlzCMeshSetDispFromField.c"
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
		abs = 1,
  		ok = 1,
  		usage = 0;
  WlzDVertex3	bgd = {0.0};
  WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;
  char		*inMeshStr,
		*inFieldStr,
  		*outObjStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inFieldObj = NULL,
		*inMeshObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "hANLRb:o:m:";
  const char    inObjStrDef[] = "-",
  	        outObjStrDef[] = "-";

  opterr = 0;
  inMeshStr = (char *)inObjStrDef;
  inFieldStr = (char *)inObjStrDef;
  outObjStr = (char *)outObjStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'A':
        abs = 1;
	break;
      case 'N':
        interp = WLZ_INTERPOLATION_NEAREST;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'R':
        abs = 0;
        break;
      case 'o':
        outObjStr = optarg;
	break;
      case 'b':
	if(sscanf(optarg, "%lg,%lg,%lg",
		  &(bgd.vtX), &(bgd.vtX), &(bgd.vtY)) < 2)
	{
	  usage = 1;
	}
	break;
      case 'm':
        inMeshStr = optarg;
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
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inFieldStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inFieldStr == NULL) || (*inFieldStr == '\0') ||
       (inMeshStr == NULL) || (*inMeshStr == '\0') ||
       (outObjStr == NULL) || (*outObjStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    int		idx;
    char	*fStr[2];
    WlzObject	**objs[2];

    fStr[0] = inMeshStr;
    fStr[1] = inFieldStr;
    objs[0] = &inMeshObj;
    objs[1] = &inFieldObj;
    for(idx = 0; idx < 2; ++idx)
    {
      FILE	*fP = NULL;

      if((fStr[idx] == NULL) ||
	 (*fStr[idx] == '\0') ||
	 ((fP = (strcmp(fStr[idx], "-")?
		fopen(fStr[idx], "r"): stdin)) == NULL) ||
	 ((*(objs[idx]) = WlzAssignObject(
	                  WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	 (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to read object from file %s\n",
		       *argv, fStr[idx]);
        break;
      }
      if(fP && strcmp(fStr[idx], "-"))
      {
	(void )fclose(fP); fP = NULL;
      }
    }
  }
  if(ok) {
    outObj = WlzCMeshSetDispFromField(inMeshObj, inFieldObj, bgd, interp,
                                      abs, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to set mesh displacements from field, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    FILE 	*fP = NULL;

    if((fP = (strcmp(outObjStr, "-")? fopen(outObjStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outObjStr);
    }
    else
    {
      errNum = WlzWriteObj(fP, outObj);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output object to file %s\n",
		       *argv, outObjStr);
      }
    }
    if(fP && strcmp(outObjStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(inMeshObj);
  (void )WlzFreeObj(inFieldObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-A] [-N] [-L] [-R] [-o<output object>]\n"
	    "\t\t[-b #,#[,#]] [-m<mesh object>] \t\t[<field object>]\n"
            "Creates a conforming mesh object that has values set from\n"
	    "the distance field object. Mesh values outside the distance\n"
	    "field object will be set to the background value.\n"
	    "The given field object must be a compound array object with\n"
	    "the appropriate number of components for the mesh dimension.\n"
	    "By default the mesh and field objects are read from the\n"
	    "standard input, with the mesh object being read before the\n"
	    "field object.\n"
            "Example:\n"
	    "%s -o out.wlz -m mesh field.wlz\n"
            "Creates a conforming mesh object (out.wlz) that has all values\n"
	    "set by linear interpolation from the given field object\n"
	    "(field.wlz) at the nodes of the given mesh (mesh.wlz).\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -A  Displacement field values are absolute (default).\n"
            "  -N  Use nearest neighbour rather than linear interpolation.\n"
            "  -L  Use linear interpolation.\n"
	    "  -R  Displacement field values are relative to position.\n"
	    "  -o  Output file.\n"
	    "  -b  Background displacement value (default 0,0[,0]).\n"
	    "  -m  Mesh object.\n",
	    argv[0],
	    argv[0],
	    WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
