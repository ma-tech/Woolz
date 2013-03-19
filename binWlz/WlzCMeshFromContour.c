#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCMeshFromContour_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCMeshFromContour.c
* \author       Bill Hill
* \date         April 2010
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
* \brief	Constructs a conforming simplical mesh from a contour
* 		with a geometric model.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcmeshfromcontour "WlzCMeshFromContour"
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

/*!
\ingroup BinWlz
\defgroup wlzcmeshfromcontour WlzCMeshFromContour
\par Name
WlzCMeshFromContour - constructs a conforming mesh that is equivalent to the
                      given contour.
\par Synopsis
\verbatim
WlzCMeshFromContour [-h] [-o<output object>] [<input object>]
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
Constructs an eqivalent conforming mesh to the given contour.
\par Examples
\verbatim
WlzCMeshFromContour -o mesh.wlz contour.wlz
\endverbatim
Creates a new conforming mesh object that has equivalent nodes and elements
to the vertices and faces of the given 3D contour.
\par File
\ref WlzCMeshFromContour.c "WlzCMeshFromContour.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshFromGM "WlzCMeshFromGM(3)"
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
  FILE		*fP = NULL;
  char		*inObjStr,
  		*outObjStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzDomain	dom;
  WlzValues	val;
  WlzObject	*inObj = NULL,
  		*outObj = NULL;
  static char   optList[] = "ho:";
  const char    inObjStrDef[] = "-",
  	        outObjStrDef[] = "-";

  opterr = 0;
  inObjStr = (char *)inObjStrDef;
  outObjStr = (char *)outObjStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjStr = optarg;
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
    if(inObj->type != WLZ_CONTOUR)
    {
      errNum = WLZ_ERR_OBJECT_TYPE;
    }
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      dom.cm2d5 = WlzCMeshFromGM(inObj->domain.ctr->model, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      val.core = NULL;
      outObj = WlzMakeMain(WLZ_CMESH_2D5, dom, val, NULL, NULL, &errNum);
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
            "Usage: %s [-h] [-o<output object>] [<input object>]\n"
            "Constructs a conforming mesh that is equivalent to the given\n"
	    "contour.\n"
	    "Version: %s\n"
	    "Options:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n",
	    argv[0],
	    WlzVersion());

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
