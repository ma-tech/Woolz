#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshCurvToImage_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshCurvToImage.c
* \author       Bill Hill
* \date         June 2010
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
* \brief	Creates a 2D domain object (image) in which the values are
*               the curvature of the 2D5 CMesh. The given CMesh must have
*               displacements that will flatten the mesh already computed.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshcurvtoimage "WlzCMeshCurvToImage"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshcurvtoimage WlzCMeshCurvToImage
\par Name
WlzCMeshCurvToImage - creates a grey scale image from a conforming mesh in
                      which the grey values are the interpolated mesh
		      Gaussian curvature values.
\par Synopsis
\verbatim
WlzCMeshCurvToImage [-h] [-m] [-N] [-L] [-Q] [-o<output object>] [-s]
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
    <td>Nearest neighbour interpolation from nearest element node value.</td>
  </tr>
  <tr>
    <td><b>-L</b></td>
    <td>Linear interpolation from element node values.</td>
  </tr>
  <tr>
    <td><b>-Q</b></td>
    <td>Interpolation from nodes surrounding element.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Computes the mean rather than the Gaussian curvature values..</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Scale factor from mesh to spatial domain.</td>
  </tr>
</table>
\par Description
Creates a 2D domain object with grey values (image) in which the values
that are interpolated from the Gaussian curvature of the mesh. The 2D
domain is created to cover the mesh following it's displacement. The
displaced mesh must have been computed so that the displacements
transform it to a plane. As for example by WlzCMeshSurfaceMap(1).
By default files are read from the standard input and written to the
standard output.
\par Examples
\verbatim
WlzCMeshCurvToImage -o out.wlz in.wlz
\endverbatim
Reads a 2D5 conforming mesh from the file in.wlz, applies the meshes
displacements to create a 2D domain object with the values set using
the Gaussian curvature of the mesh prior to he application of it's
displacements. The 2D domain object image is then writen to the file
out.wlz.
\par File
\ref WlzCMeshCurvToImage.c "WlzCMeshCurvToImage.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshCurvToImage "WlzCMeshCurvToImage(3)"
\ref WlzCMeshFromContour "WlzCMeshFromContour(1)"
\ref WlzCMeshSurfaceMap "WlzCMeshSurfaceMap(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
extern char     *optarg;
extern int      optind,
		opterr,
		optopt;

int             main(int argc, char **argv)
{
  int		option,
		meanCrv = 0,
  		ok = 1,
		usage = 0;
  double	scale = 1.0;
  WlzInterpolationType interp = WLZ_INTERPOLATION_LINEAR;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "hmNLQo:s:",
  		fileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = fileStrDef;
  outObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'N':
        interp = WLZ_INTERPOLATION_NEAREST;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'Q':
        interp = WLZ_INTERPOLATION_ORDER_2;
	break;
      case 'm':
        meanCrv = 1;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%lg", &scale) != 1)
	{
	  usage = 1;
	}
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inObjFileStr = *(argv + optind);
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    if((inObjFileStr == NULL) ||
	(*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
	       fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP)
    {
      if(strcmp(inObjFileStr, "-"))
      {
	(void )fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    outObj = WlzCMeshCurvToImage(inObj, scale, meanCrv, interp, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	     "%s: Failed to create contour from conforming mesh (%s).\n",
	     *argv, errMsg);

    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-m] [-N] [-L] [-Q] [-o<output object>] [-s#]\n"
    "      [<input object>]\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -N  Nearest neighbour interpolation from nearest element node value.\n"
    "  -L  Linear interpolation from element node values.\n"
    "  -Q  Interpolation from nodes surrounding element.\n"
    "  -m  Set image values to the mean rather than the Gaussian curvature.\n"
    "  -o  Output object file.\n"
    "  -s  Additional scale factor to be used in going from the mesh to the\n"
    "      spatial domain.\n"
    "Creates a 2D domain object with grey values (image) in which the values\n"
    "that are interpolated from the Gaussian curvature of the mesh. The 2D\n"
    "domain is created to cover the mesh following it's displacement. The\n"
    "displaced mesh must have been computed so that the displacements\n"
    "transform it to a plane (as by WlzCMeshSurfaceMap(1)).\n"
    "By default files are read from the standard input and written to the\n"
    "standard output.\n",
    *argv,
    " -o out.wlz in.wlz\n"
    "Reads a 2D5 conforming mesh from the file in.wlz, applies the meshes\n"
    "displacements to create a 2D domain object with the values set using\n"
    "the Gaussian curvature of the mesh prior to he application of it's\n"
    "displacements. The 2D domain object image is then writen to the file\n"
    "out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
