#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshToContour_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshToContour.c
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
* \brief	Creates a contour from a conforming mesh.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshtocontour "WlzCMeshToContour"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshtocontour WlzCMeshToContour
\par Name
WlzCMeshToContour - creates a contour from a conforming mesh.
\par Synopsis
\verbatim
WlzCMeshToContour [-h] [-d #] [-o<output object>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Apply this scale factor to the displacements of mesh if they exist.
        By default no mesh displacements are applied.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
</table>
\par Description
Creates a contour object from a conforming mesh.
This is currently only possible for a 2D5 conforming mesh, which
results in a 3D contour (collection of surfaces).
By default the input file is read from the standard input and the output
file is written to the standard output.
\par Examples
\verbatim
WlzCMeshToContour -o out.wlz -d in.wlz
\endverbatim
Reads a 2D5 conforming mesh from the file in.wlz, applies the meshes
displacements to it's nodes and then outputs the corresponding contour
to the file out.wlz.
\par File
\ref WlzCMeshToContour.c "WlzCMeshToContour.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshToContour "WlzCMeshToContour(3)"
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
  double	disp = 0.0;
  WlzObject     *inObj = NULL,
  		*outObj = NULL;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outObjFileStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "hd:o:",
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
      case 'd':
	if(sscanf(optarg, "%lg", &disp) != 1)
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
    outObj = WlzCMeshToContour(inObj, disp, &errNum);
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
    " [-h] [-d #] [-o<output object>] [<input object>]\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -d  Apply given scale factor to mesh displacements if they exist.\n"
    "      The default is no mesh displacements are applied.\n"
    "  -o  Output object file.\n"
    "Creates a contour object from a conforming mesh. This is currently\n"
    "only possible for a 2D5 conforming mesh, which results in a 3D contour\n"
    "(collection of surfaces).  By default the input file is read from the\n"
    "standard input and the output file is written to the standard output.\n",
    *argv,
    "  -o out.wlz -d 1.0 in.wlz\n"
    "Reads a 2D5 conforming mesh from the file in.wlz, applies the meshes\n"
    "displacements to it's nodes and then outputs the corresponding contour\n"
    "to the file out.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
