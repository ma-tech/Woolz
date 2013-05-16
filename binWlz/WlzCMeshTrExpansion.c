#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshTrExpansion_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzCMeshTrExpansion.c
* \author       Bill Hill
* \date         February 2013
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2013],
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
* \brief	Creates a conforming mesh with a a normalised expansion
* 		scalar attached to each element from the given conforming
* 		mesh transform.
* \ingroup	BinWlz
*
* \par 		Binary
* \ref 		wlzcmeshtrexpansion "WlzCMeshTrExpansion"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshtrexpansion WlzCMeshTrExpansion
\par Name
WlzCMeshTrExpansion - computes strain tensors for the given conforming
                         mesh transform.
\par Synopsis
\verbatim
WlzCMeshTrExpansion [-h] [-i] [-I] [-L] [-o<out file>] [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-i</b></td>
    <td>Output spatial domain image not a conforming mesh.</td>
  </tr>
  <tr> 
    <td><b>-I</b></td>
    <td>Use transform inverse.</td>
  </tr>
  <tr> 
    <td><b>-E</b></td>
    <td>Use maximum tensor eigenvalue rather than trace.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use kriging rather than barycentric interpolation./td>
  </tr>
  <tr> 
    <td><b>-n</b></td>
    <td>Normalise expansion factors./td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
\par Examples
Computes the expansion of each element of the given mesh transform.
The expansion may be output as a mesh with scalar values attached to
the mesh elements or as a 2D spatial domain object with image values.
If requested the image values are normalised so that:
  [0:max(abs(min),abs(max))] -> [128:25]
By default all files are read from the standard input and are written
to the standard output.
\verbatim
WlzCMeshTrExpansion -o meshexp.wlz -i -n meshtr.wlz
\endverbatim
Creates a new spatial domain (image) object (written to meshexp.wlz)
with the mesh expansion factors normalised to fit within the unsigned
byte range and with 0 to the value 128.
\par File
\ref WlzCMeshTrExpansion.c "WlzCMeshTrExpansion.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzAssignObject "WlzAssignObject(WlzCMeshExpansion(3)"
\ref WlzCMeshValuesNormalise "WlzCMeshValuesNormalise(3)"
\ref WlzSetMeshInverse "WlzCMeshToDomObj(3)"
\ref WlzCMeshToDomObjValues "WlzCMeshToDomObjValues(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>


/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		option,
		image = 0,
		eigen = 0,
		interp = 0,
		inverse = 0,
		norm = 0,
  		ok = 1,
  		usage = 0;
  FILE		*fP = NULL;
  char		*inFileStr = NULL,
  		*outFileStr = NULL;
  WlzObject	*inObj = NULL,
  		*outObj = NULL,
		*meshExpObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char   optList[] = "hiIELno:";
  const char    fileStrDef[] = "-";

  opterr = 0;
  inFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'i':
        image = 1;
	break;
      case 'E':
        eigen = 1;
	break;
      case 'I':
        inverse = 1;
	break;
      case 'L':
        interp = 1;
	break;
      case 'n':
        norm = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if(optind < argc)
    {
      if((optind + 1) != argc)
      {
        usage = 1;
        ok = 0;
      }
      else
      {
        inFileStr = *(argv + optind);
      }
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if(((fP = (strcmp(inFileStr, "-")?
              fopen(inFileStr, "r"): stdin)) == NULL) ||
       ((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inFileStr);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    meshExpObj = WlzAssignObject(WlzCMeshExpansion(inObj, inverse, eigen,
                                                   &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to compute mesh expansion (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  WlzFreeObj(inObj);
  if(ok && norm)
  {
    errNum = WlzCMeshValuesNormalise(meshExpObj, 1, 128.0, 0.0, 255.0,
    				     0.001);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to normalise mesh expansion (%s)\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if(image)
    {
      WlzObject *objD;

      objD = WlzCMeshToDomObj(meshExpObj, 0, 1.0, &errNum);
      if(errNum == WLZ_ERR_NONE)
      {
        outObj = WlzAssignObject(
      	         WlzCMeshToDomObjValues(objD, meshExpObj,
					(interp)? WLZ_INTERPOLATION_LINEAR:
					          WLZ_INTERPOLATION_NEAREST,
					&errNum), NULL);
	         
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )fprintf(stderr,
		       "%s: failed to create spatial domain object (%s)\n",
		       *argv, outFileStr);
      }
    }
    else
    {
      outObj = WlzAssignObject(meshExpObj, NULL);
    }
  }
  (void )WlzFreeObj(meshExpObj);
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write object to file %s\n",
                     *argv, outFileStr);
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
      "Usage: %s [-h] [-i]  [-p] [-o<out file>] [-s<x>,<y>,<z>]\n"
      "\t\t[<input file>]\n"
      "Computes the expansion of each element of the given mesh transform.\n"
      "The expansion may be output as a mesh with scalar values attached to\n"
      "the mesh elements or as a 2D spatial domain object with image values.\n"
      "If requested the image values are normalised so that:\n"
      "  [0:max(abs(min),abs(max))] -> [128:255]\n"
      "By default all files are read from the standard input and are written\n"
      "to the standard output.\n"
      "Version: %s\n"
      "Options are:\n"
      "  -h  Help, prints this usage message.\n"
      "  -i  Output spatial domain image not a conforming mesh.\n"
      "  -E  Use maximum tensor eigenvalue rather than trace.\n"
      "  -I  Use inverse transform.\n"
      "  -L  Use kriging rather than barycentric interpolation.\n"
      "  -n  Normalise expansion factors.\n"
      "  -o  Output object.\n"
      "%s -o meshexp.wlz -i -n meshtr.wlz\n"
      "Creates a new spatial domain (image) object (written to meshexp.wlz)\n"
      "with the mesh expansion factors normalised to fit within the unsigned\n"
      "byte range and with 0 to the value 128.\n",
      argv[0],
      WlzVersion(),
      argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
