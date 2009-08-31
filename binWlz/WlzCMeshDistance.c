#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzMeshGen_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzCMeshDistance.c
* \author       Bill Hill
* \date         January 2009
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief        Constructs a 2D or 3D domain object the values
* 		of which are the minimum distance from the seeds in
* 		a conforming mesh.
* \ingroup	Wlz
* \todo         -
* \bug          None known.
* \par 		Binary
* \ref 		wlzcmeshdistance "WlzCMeshDistance"
*/

/*!
\ingroup BinWlz
\defgroup wlzcmeshdistance WlzCMeshDistance
\par Name
WlzCMeshDistance - computes distances within conforming meshes.
\par Synopsis
\verbatim
WlzCMeshDistance [-b] [-h] [-o<out obj file>] [-r<ref obj file>] [-s #,#,#]
                 [<input mesh file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-b</b></td>
    <td>Use the boundary of the mesh for seed points.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Use the boundary of the reference object for seed points.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Use the single seed position.</td>
  </tr>
</table>
\par Description
Constructs a 2D or 3D domain object the values of which are the
minimum distance from the given seeds points in the given conforming
mesh. The domain of the output object covers the input mesh.
\par Examples
\verbatim
WlzCMeshDistance -s 100,200,0 -o out.wlz mesh.wlz
\endverbatim
Creates a new domain object with values, in which the domain covers the
given mesh and the values are distances from the seed (at column 100,
line 200, plane 0). The new domain object is written to the file out.wlz.
\par File
\ref WlzCMeshDistance.c "WlzDistanceTransform.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzCMeshDistance2D "WlzCMeshDistance2D(3)"
\ref WlzCMeshDistance3D "WlzCMeshDistance3D(3)"
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
  int		idN,
  		idM,
		nSeeds,
  		nBndSeeds,
		boundFlg = 0,
  		seedFlg = 0,
  		ok = 1,
  		option,
  		usage = 0;
  size_t	vtxSize = sizeof(WlzDVertex2);
  WlzVertexType vtxType,
  		bndVtxType;
  WlzVertex	seed;
  WlzVertexP	bndSeeds,
		seeds;
  FILE		*fP = NULL;
  char		*outObjFileStr,
  		*refObjFileStr = NULL,
		*meshFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*outObj = NULL,
		*meshObj = NULL,
  		*refObj = NULL;
  WlzCMeshP 	mesh;
  WlzCMeshNodP	nod;
  static char   optList[] = "bho:r:s:";
  const char    meshFileStrDef[] = "-",
  	        outObjFileStrDef[] = "-";

  mesh.v = NULL;
  opterr = 0;
  seeds.v = NULL;
  bndSeeds.v = NULL;
  outObjFileStr = (char *)outObjFileStrDef;
  meshFileStr = (char *)meshFileStrDef;
  seed.d3.vtX = seed.d3.vtY = seed.d3.vtZ = 0.0;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
	boundFlg = 1;
        break;
      case 's':
	seedFlg = 1;
        if(sscanf(optarg, "%lg,%lg,%lg",
	          &(seed.d3.vtX), &(seed.d3.vtY), &(seed.d3.vtZ)) < 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'r':
        refObjFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  ok = usage == 0;
  if(ok)
  {
    if((meshFileStr == NULL) || (*meshFileStr == '\0') ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
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
        meshFileStr = *(argv + optind);
      }
    }
  }
  if(ok && refObjFileStr)
  {
    if(((fP = (strcmp(refObjFileStr, "-")?
              fopen(refObjFileStr, "r"): stdin)) == NULL) ||
       ((refObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read reference object from file %s\n",
                     *argv, refObjFileStr);
    }
    if(fP && strcmp(refObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    if(((fP = (strcmp(meshFileStr, "-")?
              fopen(meshFileStr, "r"): stdin)) == NULL) ||
       ((meshObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read mesh object from file %s\n",
                     *argv, meshFileStr);
    }
    if(fP && strcmp(meshFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    switch(meshObj->type)
    {
      case WLZ_CMESH_2D:
	mesh.m2 = meshObj->domain.cm2;
	break;
      case WLZ_CMESH_3D:
	mesh.m3 = meshObj->domain.cm3;
	break;
      default:
        ok = 0;
	errNum = WLZ_ERR_OBJECT_TYPE;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		 "%s: Invalid mesh object, must be WLZ_CMESH_[23]D (%s),\n",
		 argv[0],
		 errMsgStr);
	break;
    }
  }
  if(ok && meshObj)
  {
    switch(meshObj->type)
    {
      case WLZ_CMESH_2D:
	vtxType = WLZ_VERTEX_D2;
	vtxSize = sizeof(WlzDVertex2);
        if((refObj != NULL) && (refObj->type != WLZ_2D_DOMAINOBJ))
	{
	  ok = 0;
	}
	break;
      case WLZ_CMESH_3D:
	vtxType = WLZ_VERTEX_D3;
	vtxSize = sizeof(WlzDVertex3);
        if((refObj != NULL) && (refObj->type != WLZ_3D_DOMAINOBJ))
	{
	  ok = 0;
	}
	break;
      default:
        ok = 0;
    }
    if(ok == 0)
    {
	errNum = WLZ_ERR_OBJECT_TYPE;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		 "%s: Mesh and reference objects must have the same\n"
		 "dimension (%s)\n",
		 argv[0],
		 errMsgStr);
    }
  }
  /* Create seed buffer. */
  if(ok)
  {
    nSeeds = (seedFlg != 0);
    if(boundFlg)
    {
      nSeeds += WlzCMeshSetBoundNodFlags(mesh);
    }
    if(refObj)
    {
      bndSeeds = WlzVerticesFromObjBnd(refObj, &nBndSeeds, &bndVtxType,
      				       &errNum);
      nSeeds += nBndSeeds;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
               "%s: Failed to extract seeds (%s)\n",
	       argv[0],
	       errMsgStr);
    }
    else if(nSeeds < 1)
    {
      ok = 0;
      (void )fprintf(stderr,
	       "%s: Number of seeds must be > 0.\n",
	       argv[0]);
    }
  }
  if(ok)
  {
    if((seeds.v = AlcMalloc(vtxSize * nSeeds)) == NULL)
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to allocate seed buffer (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  if(ok)
  {
    switch(meshObj->type)
    {
      case WLZ_CMESH_2D:
	idN = 0;
	if(seedFlg)
	{
	  idN = 1;
	  seeds.d2->vtX = seed.d2.vtX;
	  seeds.d2->vtY = seed.d2.vtY;
	}
	if(boundFlg)
	{
	  for(idM = 0; idM < mesh.m2->res.nod.maxEnt; ++idM)
	  {
	    nod.n2 = (WlzCMeshNod2D *)
	             AlcVectorItemGet(mesh.m2->res.nod.vec, idM);
	    if((nod.n2->idx >= 0) &&
	       ((nod.n2->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0))
	    {
	      seeds.d2[idN++] = nod.n2->pos;
	    }
	  }
	}
	if(refObj)
	{
	  switch(bndVtxType)
	  {
	    case WLZ_VERTEX_I2:
	      WlzValueCopyIVertexToDVertex(seeds.d2 + idN, bndSeeds.i2,
	      				   nBndSeeds);
	      break;
	    case WLZ_VERTEX_F2:
	      WlzValueCopyFVertexToDVertex(seeds.d2 + idN, bndSeeds.f2,
	      				   nBndSeeds);
	      break;
	    case WLZ_VERTEX_D2:
	      WlzValueCopyDVertexToDVertex(seeds.d2 + idN, bndSeeds.d2,
	      				   nBndSeeds);
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzCMeshDistance2D(mesh.m2, nSeeds, seeds.d2, &errNum);
	}
	break;
      case WLZ_CMESH_3D:
	idN = 0;
	if(seedFlg)
	{
	  idN = 1;
	  seeds.d3->vtX = seed.d3.vtX;
	  seeds.d3->vtY = seed.d3.vtY;
	  seeds.d3->vtZ = seed.d3.vtZ;
	}
	if(boundFlg)
	{
	  for(idM = 0; idM < mesh.m3->res.nod.maxEnt; ++idM)
	  {
	    nod.n3 = (WlzCMeshNod3D *)
	             AlcVectorItemGet(mesh.m3->res.nod.vec, idM);
	    if((nod.n3->idx >= 0) &&
	       ((nod.n3->flags & WLZ_CMESH_NOD_FLAG_BOUNDARY) != 0))
	    {
	      seeds.d3[idN++] = nod.n3->pos;
	    }
	  }
	}
	if(refObj)
	{
	  switch(bndVtxType)
	  {
	    case WLZ_VERTEX_I3:
	      WlzValueCopyIVertexToDVertex3(seeds.d3 + idN, bndSeeds.i3,
	      				    nBndSeeds);
	      break;
	    case WLZ_VERTEX_F3:
	      WlzValueCopyFVertexToDVertex3(seeds.d3 + idN, bndSeeds.f3,
	      				    nBndSeeds);
	      break;
	    case WLZ_VERTEX_D3:
	      WlzValueCopyDVertexToDVertex3(seeds.d3 + idN, bndSeeds.d3,
	      				    nBndSeeds);
	      break;
	    default:
	      errNum = WLZ_ERR_PARAM_TYPE;
	      break;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  outObj = WlzCMeshDistance3D(mesh.m3, nSeeds, seeds.d3, &errNum);
	}
	break;
	break;
      default:
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
		"%s: Failed to create domain object with distances (%s)\n",
		argv[0],
		errMsgStr);
    }
  }
  AlcFree(seeds.v);
  AlcFree(bndSeeds.v);
  WlzFreeObj(refObj);
  WlzFreeObj(meshObj);
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to write output object to file %s\n",
                     *argv, outObjFileStr);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-b] [-h] [-o<out obj file>] [-r<ref obj file>]\n"
	    "                        [-s #,#,#] [<input mesh file>]\n"
	    "Constructs a 2D or 3D domain object the values of which are\n"
	    "the minimum distance from the given seeds points in the given\n"
	    "conforming mesh. The domain of the output object covers the\n"
	    "input mesh. Options are:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output object.\n"
	    "  -b  Set seed points around the boundary of the mesh.\n"
	    "  -s  Single seed position.\n"
	    "  -r  Reference object with seed points around it's boundary.\n",
	    argv[0]);

  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
