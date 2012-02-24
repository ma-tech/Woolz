#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTstCMeshVtxInMesh_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzTst/WlzTstCMeshVtxInMesh.c
* \author       Bill Hill
* \date         May 2009
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
* \brief	Test for vertex in conforming meshes.
* \ingroup	BinWlzTst
*/


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
  int		idN,
		eIdx = -1,
		nIdx = -1,
  		ok = 1,
		exhaustive = 0,
  		option,
		repeats = 1,
  		usage = 0;
  WlzDVertex3	vtx;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL;
  WlzCMeshP 	mesh;
  static char   optList[] = "eho:p:R:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  mesh.v = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'e':
        exhaustive = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'p':
	if(sscanf(optarg, "%lg,%lg,%lg",
		  &(vtx.vtX), &(vtx.vtY), &(vtx.vtZ)) != 3)
	{
	  usage = 1;
	}
	break;
      case 'R':
        if(sscanf(optarg, "%d", &repeats) != 1)
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
  if(usage == 0)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
    {
      usage = 1;
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
  }
  ok = usage == 0;
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
      (void )fprintf(stderr,
                     "%s: Failed to read object from file (%s)\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    switch(inObj->type)
    {
      case WLZ_CMESH_2D:
        mesh.m2 = inObj->domain.cm2;
	break;
      case WLZ_CMESH_3D:
        mesh.m3 = inObj->domain.cm3;
	break;
      default:
        ok = 0;
	errNum = WLZ_ERR_OBJECT_TYPE;
        (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		   "%s: Invalid object type, must be WLZ_CMESH_[23]D (%s),\n",
		       argv[0],
		       errMsgStr);
        break;
    }
  }
  if(ok)
  {
    for(idN = 0; idN < repeats; ++idN)
    {
      eIdx = WlzCMeshElmEnclosingPos(mesh, -1, vtx.vtX, vtx.vtY, vtx.vtZ,
                                     exhaustive, &nIdx);
    }
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    if(fprintf(fP, "%d %d\n", eIdx, nIdx) <= 0)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s.\n",
		     argv[0], outFileStr);
    }
  }
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>] [-p<position>] [-R<repeats>]\n"
    "       [<input cmesh object>]\n"
    "Reads a conforming mesh and then searches for the element which\n"
    "contains the given vertex.\n"
    "The output is either:\n"
    "  <element index> <node index>\n"
    "where the element index is the index of the mesh element enclosing the\n"
    "given vertex and the node index is the index of the node closest to\n"
    "the given vertex. If the vertex is not within a mesh element then the\n"
    "index will be a negative value. If the element index is not negative\n"
    "then the node index may not be the closest node.\n"
    "Options are:\n"
    "  -e  Exhaustive search (slow), every element in the mesh searched.\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file.\n"
    "  -p  Given vertex position given as <vx>,<vy>,<vz>.\n"
    "  -R  number of times to repeat the computation.\n",
    argv[0]);

  }
  return(!ok);
}
