#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshCurvToImage_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshCurvToImage.c
* \author       Bill Hill
* \date         June 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Creates a contour from a conforming mesh.
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
WlzCMeshCurvToImage [-h] [-o<output object>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
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
  		ok = 1,
		usage = 0;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ho:",
  		fileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = fileStrDef;
  outObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
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
    outObj = WlzCMeshCurvToImage(inObj, &errNum);
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
    " [-h] [-o<output object>] [<input object>]\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -o  Output object file.\n"
    "Creates a 2D domain object with grey values (image) in which the values\n"
    "that are interpolated from the Gaussian curvature of the mesh. The 2D\n"
    "domain is created to cover the mesh following it's displacement. The\n"
    "displaced mesh must have been computed so that the displacements\n"
    "transform it to a plane. As for example by WlzCMeshSurfaceMap(1).\n"
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
