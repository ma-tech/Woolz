#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzCMeshIntersectDom_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzCMeshIntersectDom.c
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
* \brief	Creates a 2D domain object (image) in which the values are
*               the curvature of the 2D5 CMesh. The given CMesh must have
*               displacements that will flatten the mesh already computed.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshintersectdom "WlzCMeshIntersectDom"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshintersectdom WlzCMeshIntersectDom
\par Name
WlzCMeshIntersectDom - creates a grey scale image from a conforming mesh in
                      which the grey values are the interpolated mesh
		      Gaussian curvature values.
\par Synopsis
\verbatim
WlzCMeshIntersectDom [-h] [-m] [-o<output object>] [-s] [<input object>]
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
  <tr>
    <td><b>-s</b></td>
    <td>Scale factor from the mesh to the 2D spatial domain.</td>
  </tr>
</table>
\par Description
Creates a 2D domain object that is the (possibly scaled) intersection
of the given 2.5D conforming mesh and the given 3D domain, with the
intersection mapped from 3D to 2D using the displacements of the 2.5D mesh.
The mesh must have valid displacements that trfansform it to a plane
(as computed by WlzCMeshSurfaceMap(1)).
By default files are read from the standard input and written to the
standard output. If both the mesh and domain files are read from the
standard input, the mesh is read before the domain object.
\par Examples
\verbatim
WlzCMeshIntersectDom -s10 -o out.wlz mesh.wlz in.wlz
\endverbatim
Reads a 2D5 conforming mesh from the file mesh.wlz,
computes the domain of intersection and scales it by a factor of 10.
The 2D domain object is then writen to the file out.wlz.
\par File
\ref WlzCMeshIntersectDom.c "WlzCMeshIntersectDom.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshIntersectDom "WlzCMeshIntersectDom(3)"
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
  double	scale = 1.0;
  WlzObject     *inDomObj = NULL,
  		*inMeshObj = NULL,
  		*outObj = NULL;
  FILE		*fP = NULL;
  char		*inDomObjFileStr,
  		*inMeshObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ho:s:",
  		fileStrDef[] = "-";

  opterr = 0;
  inDomObjFileStr = fileStrDef;
  inMeshObjFileStr = fileStrDef;
  outObjFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
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
  if(usage == 0)
  {
    int		cnt;
    
    cnt = argc - optind;
    switch(cnt)
    {
      case 0:
        break;
      case 1:
        inMeshObjFileStr = *(argv + optind);
	break;
      case 2:
        inMeshObjFileStr = *(argv + optind);
        inDomObjFileStr  = *(argv + optind + 1);
	break;
      default:
        usage = 1;
	break;
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    if((inMeshObjFileStr == NULL) ||
	(*inMeshObjFileStr == '\0') ||
	((fP = (strcmp(inMeshObjFileStr, "-")?
	       fopen(inMeshObjFileStr, "r"): stdin)) == NULL) ||
	((inMeshObj = WlzAssignObject(
	              WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP)
    {
      if(strcmp(inMeshObjFileStr, "-"))
      {
	(void )fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    if((inDomObjFileStr == NULL) ||
	(*inDomObjFileStr == '\0') ||
	((fP = (strcmp(inDomObjFileStr, "-")?
	        fopen(inDomObjFileStr, "r"): stdin)) == NULL) ||
	((inDomObj = WlzAssignObject(
	             WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
    }
    if(fP)
    {
      if(strcmp(inDomObjFileStr, "-"))
      {
	(void )fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    outObj = WlzCMeshIntersectDom2D5(inDomObj, inMeshObj, scale, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	     "%s: Failed to compute domain of intersection (%s).\n",
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
  (void )WlzFreeObj(inDomObj);
  (void )WlzFreeObj(inMeshObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-h] [-s#] [<mesh object>] [<domain object>]\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -o  Output object file.\n"
    "  -s  Additional scale factor to be used in going from the mesh to the\n"
    "      2D spatial domain.\n"
    "Creates a 2D domain object that is the (possibly scaled) intersection\n"
    "of the given 2.5D conforming mesh and the given 3D domain, with the\n"
    "intersection mapped from 3D to 2D using the displacements of the 2.5D\n"
    "mesh.\n"
    "The mesh must have valid displacements that trfansform it to a plane\n"
    "(as computed by WlzCMeshSurfaceMap(1)).\n"
    "By default files are read from the standard input and written to the\n"
    "standard output. If both the mesh and domain files are read from the\n"
    "standard input, the mesh is read before the domain object.\n",
    *argv,
    " -s 10 -o out.wlz mesh.wlz in.wlz\n"
    "Reads a 2D5 conforming mesh from the file mesh.wlz, computes the\n"
    "domain of intersection and scales it by a factor of 10.\n"
    "The 2D domain object is then writen to the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
