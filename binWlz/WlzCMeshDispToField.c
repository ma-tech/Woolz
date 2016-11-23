#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshToDispField_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshDispToField.c
* \author       Bill Hill
* \date         October 2016
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
* \brief	Creates a displacement field object from the
* 		given conforming mesh object.
* 		object.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcmeshdisptofield "WlzCMeshDispToField"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshdisptofield WlzCMeshDispToField
\par Name
WlzCMeshDispToField - Creates a displacement field object from a conforming
		      mesh transform.
\par Synopsis
\verbatim
WlzCMeshDispToField [-h] [-A] [-L] [-N] [-R] [-i] [-b #,#[,#]
                    [-o<output object>] [<mesh transform object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-A</b></td>
    <td>Displacement field values are absolute.</td>
  </tr> 
  <tr> 
    <td><b>-N</b></td>
    <td>Use nearest neighbour rather than linear extrapolation.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use linear interpolation (default).</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-R</b></td>
    <td>Displacement field values are relative to position (default).</td>
  </tr> 
  <tr>
    <td><b>-b</b></td>
    <td>Background displacement value (default 0,0[,0]).</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Invert the transform.</td>
  </tr>
</table>
\par Description
Creates a displacement field object that has values set from
the conforming mesh transform. Field values outside the mesh
will be set to the background value.
The resulting field object will be a compound array object with
the appropriate number of components for the mesh dimension.
By default the mesh transform object is read from the
standard input and the field object is written to the standard output.
\par Examples
\verbatim
WlzCMeshDispToField -o field.wlz meshtr.wlz
\endverbatim
Creates a displacement field object (field.wlz) that has values set by
linear interpolation from the given mesh transform object (meshtr.wlz)
displacements.
\par File
\ref WlzCMeshDispToField.c "WlzCMeshDispToField.c"
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
		abs = 0,
		invert = 0,
  		ok = 1,
  		usage = 0;
  WlzDVertex3	bgd = {0.0};
  WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;
  char		*inMeshStr,
		*outFieldStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inMeshObj = NULL,
		*outFieldObj = NULL;
  static char   optList[] = "hANLRb:o:";
  const char    inObjStrDef[] = "-",
  	        outObjStrDef[] = "-";

  opterr = 0;
  inMeshStr = (char *)inObjStrDef;
  outFieldStr = (char *)outObjStrDef;
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
      case 'i':
        invert = 1;
	break;
      case 'o':
        outFieldStr = optarg;
	break;
      case 'b':
	if(sscanf(optarg, "%lg,%lg,%lg",
		  &(bgd.vtX), &(bgd.vtX), &(bgd.vtY)) < 2)
	{
	  usage = 1;
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
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inMeshStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inMeshStr == NULL) || (*inMeshStr == '\0') ||
       (outFieldStr == NULL) || (*outFieldStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    FILE	*fP = NULL;

    if(((fP = (strcmp(inMeshStr, "-"))?
              fopen(inMeshStr, "r"): stdin) == NULL) ||
       ((inMeshObj = WlzAssignObject(
		     WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inMeshStr);
    }
    if(fP && strcmp(inMeshStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok) {
    outFieldObj = WlzCMeshDispToField(inMeshObj, bgd, interp, invert,
                                      abs, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to set field from mesh displacements, %s.\n",
      		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    FILE 	*fP = NULL;

    if((fP = (strcmp(outFieldStr, "-")?
             fopen(outFieldStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outFieldStr);
    }
    else
    {
      errNum = WlzWriteObj(fP, outFieldObj);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to write output object to file %s\n",
		       *argv, outFieldStr);
      }
    }
    if(fP && strcmp(outFieldStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(inMeshObj);
  (void )WlzFreeObj(outFieldObj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-A] [-L] [-N] [-R] [-i] [-b #,#[,#]]\n"
	    "\t\t[-o<output object>] [<mesh transform object>]\n"
            "Creates a displacement field object that has values set from\n"
	    "the conforming mesh transform. Field values outside the mesh\n"
	    "will be set to the background value.\n"
	    "The resulting field object will be a compound array object with\n"
	    "the appropriate number of components for the mesh dimension.\n"
	    "By default the mesh transform object is read from the standard\n"
	    "input and the field object is written to the standard output.\n"
            "Example:\n"
	    "%s -o field meshtr.wlz\n"
            "Creates a displacement field object (field.wlz) that has values\n"
	    "set by linear interpolation from the given mesh transform object\n"
	    "(meshtr.wlz) displacements.\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -A  Displacement field values are absolute.\n"
            "  -L  Use linear interpolation.\n"
            "  -N  Use nearest neighbour rather than linear interpolation.\n"
	    "  -R  Displacement field values are relative to position\n"
	    "      (default).\n"
	    "  -i  Invert the transform.\n"
	    "  -b  Background displacement value (default 0,0[,0]).\n"
	    "  -o  Output file.\n",
	    argv[0],
	    argv[0],
	    WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
