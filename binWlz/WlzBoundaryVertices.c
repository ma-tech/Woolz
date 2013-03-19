#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBoundaryVertices_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzBoundaryVertices.c
* \author       Bill Hill
* \date         May 2010
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
* \brief	Extracts boundary vertex positions from an object.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzboundaryvertices "WlzBoundaryVertices"
*/

/*!
\ingroup BinWlz
\defgroup wlzboundaryvertices WlzBoundaryVertices
\par Name
WlzBoundaryVertices - extracts boundary vertex positions from an object.
\par Synopsis
\verbatim
WlzBoundaryVertices [-h] [-o<output file>] [<input object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
</table>
\par Description
Given an object that is represented by connected vertices (eg a mesh)
this filter outputs the coordinates of the vertices that are on
the boundary of the object.
By default the input object if read from the standard input and the
vertices are written to the standard output.
\par Examples
\verbatim
WlzBoundaryVertices -o out.num in.wlz
\endverbatim
The boundary vertices of the object read from the file in.wlz are
written to the file out.num as white space seperated acsii text.
\par File
\ref WlzBoundaryVertices.c "WlzBoundaryVertices.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzBoundaryVertices "WlzBoundaryVertices(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

static WlzErrorNum 		WlzBoundaryVerticesCMesh(
				  FILE *fP,
				  WlzObject *obj);

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
  FILE		*fP = NULL;
  char 		*inObjFileStr,
  		*outFileStr;
  WlzObject     *inObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ho:",
  		fileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = fileStrDef;
  outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileStr = optarg;
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
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
    else
    {
      switch(inObj->type)
      {
        case WLZ_CMESH_2D:  /* FALLTHROUGH */
	case WLZ_CMESH_2D5: /* FALLTHROUGH */
	case WLZ_CMESH_3D:
	  errNum = WlzBoundaryVerticesCMesh(fP, inObj);
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to output boundary vertices (%s).\n",
		       *argv, errMsg);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-o<output file>] [<input object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -h  Help, prints usage message.\n"
    "  -o  Output file.\n"
    "Given an object that is represented by connected vertices (eg a mesh)\n"
    "this filter outputs the coordinates of the vertices that are on\n"
    "the boundary of the object.\n"
    "By default the input object if read from the standard input and the\n"
    "vertices are written to the standard output.\n",
    *argv,
    " -o out.num in.wlz\n"
    "The boundary vertices of the object read from the file in.wlz are\n"
    "written to the file out.num as white space seperated acsii text.\n");
  }
  return(!ok);
}

/*!
* \return	Woolz error code.
* \brief	Prints boundary node coordinates for a conforming mesh.
* \param	fP			File pointer opened for acsii write.
* \param	obj			Object containing the conforming mesh.
*/
static WlzErrorNum WlzBoundaryVerticesCMesh(FILE *fP, WlzObject *obj)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(obj == NULL)
  {
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else if(obj->domain.core == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else
  {
    int		idx,
    		mN;
    WlzCMeshP	mesh;
    AlcVector	*vec;

    switch(obj->type)
    {
      case WLZ_CMESH_2D:
        mesh.m2 = obj->domain.cm2;
	mN = mesh.m2->res.nod.maxEnt;
	vec = mesh.m2->res.nod.vec;
	for(idx = 0; idx < mN; ++idx)
	{
	  WlzCMeshNod2D *nod;

          nod = (WlzCMeshNod2D *)AlcVectorItemGet(vec, idx);
          if(nod->idx >= 0)
	  {
	    if(WlzCMeshNodIsBoundary2D(nod))
	    {
	      (void )fprintf(fP, "%lg %lg\n", nod->pos.vtX, nod->pos.vtY);
	    }
	  }
	}
	break;
      case WLZ_CMESH_2D5:
        mesh.m2d5 = obj->domain.cm2d5;
	mN = mesh.m2d5->res.nod.maxEnt;
	vec = mesh.m2d5->res.nod.vec;
	for(idx = 0; idx < mN; ++idx)
	{
	  WlzCMeshNod2D5 *nod;

          nod = (WlzCMeshNod2D5 *)AlcVectorItemGet(vec, idx);
          if(nod->idx >= 0)
	  {
	    if(WlzCMeshNodIsBoundary2D5(nod))
	    {
	      (void )fprintf(fP, "%lg %lg %lg\n",
	                     nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ);
	    }
	  }
	}
	break;
      case WLZ_CMESH_3D:
        mesh.m3 = obj->domain.cm3;
	mN = mesh.m3->res.nod.maxEnt;
	vec = mesh.m3->res.nod.vec;
	for(idx = 0; idx < mN; ++idx)
	{
	  WlzCMeshNod3D *nod;

          nod = (WlzCMeshNod3D *)AlcVectorItemGet(vec, idx);
          if(nod->idx >= 0)
	  {
	    if(WlzCMeshNodIsBoundary3D(nod))
	    {
	      (void )fprintf(fP, "%lg %lg %lg\n",
	                     nod->pos.vtX, nod->pos.vtY, nod->pos.vtZ);
	    }
	  }
	}
	break;
      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
    }
    if((errNum == WLZ_ERR_NONE) && feof(fP))
    {
      errNum = WLZ_ERR_WRITE_INCOMPLETE;
    }
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
