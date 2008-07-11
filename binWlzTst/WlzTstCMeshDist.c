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
* \file         binWlzTst/WlzTstCMeshGen.c
* \author       Bill Hill
* \date         June 2003
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
* \brief	Test for 2D and 3D distance computation using fast marching
* 		methods within simplical conforming meshes.
*
* \ingroup	WlzTst
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <Wlz.h>

typedef enum _WlzTstParam
{
  WLZTST_SEED_BOUNDARY,
  WLZTST_SEED_SEEDS,
  WLZTST_OUT_IMAGE,
  WLZTST_OUT_TXT
} WlzTstParam;

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

int		main(int argc, char *argv[])
{
  int		idN,
  		nSeeds = 0,
  		ok = 1,
  		option,
		repeats = 1,
  		usage = 0,
		maxNod = 0,
		seedType = WLZTST_SEED_SEEDS,
		outType = WLZTST_OUT_TXT;
  double	*dist = NULL;
  double	**inSeeds = NULL;
  size_t	inRow = 0,
  		inCol = 0;
  WlzVertexP	seeds;
  WlzCMeshNod2D	*nod2;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  AlcErrno	alcErr = ALC_ER_NONE;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*inObj = NULL;
  WlzCMeshP 	mesh;
  static char   optList[] = "bhto:s:R:S:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  opterr = 0;
  mesh.v = NULL;
  seeds.v = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'b':
	seedType = WLZTST_SEED_BOUNDARY;
        break;
      case 't':
	outType = WLZTST_OUT_TXT;
        break;
      case 's':
	seedType = WLZTST_SEED_SEEDS;
	if(AlcDouble2Malloc(&inSeeds, 1, 3) != ALC_ER_NONE)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	  usage = 1;
	}
	else
	{
	  if(sscanf(optarg,
	  	    "%lg,%lg,%lg",
	            &(inSeeds[0][0]), &(inSeeds[0][1]), &(inSeeds[0][2])) != 3)
	  {
	    usage = 1;
	  }
	  else
	  {
	    inRow = 1;
	    inCol = 3;
	  }
	}
	break;
      case 'S':
	seedType = WLZTST_SEED_SEEDS;
	if((*optarg == '\0') ||
	   ((fP = (strcmp(optarg, "-")?
		  fopen(optarg, "r"): stdin)) == NULL) ||
	   ((alcErr = AlcDouble2ReadAsci(fP, &inSeeds,
	                                 &inRow, &inCol)) != ALC_ER_NONE))
	{
	  switch(alcErr)
	  {
	    case ALC_ER_READ:
	      errNum = WLZ_ERR_READ_EOF;
	    case ALC_ER_ALLOC:
	      errNum = WLZ_ERR_MEM_ALLOC;
	      break;
	    default:
	      break;
	  }
	  usage = 1;
	}
	if(fP && strcmp(optarg, "-"))
	{
	  (void )fclose(fP); fP = NULL;
	}
	break;
      case 'R':
        if(sscanf(optarg, "%d", &repeats) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
        outFileStr = optarg;
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
  if(ok && (seedType == WLZTST_SEED_SEEDS))
  {
    if((inRow < 1) || (inCol != 3))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s Failed to read input seeds (%s).\n",
		     argv[0],
		     errMsgStr);

    }
  }
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
	maxNod = mesh.m2->res.nod.maxEnt;
	break;
      case WLZ_CMESH_3D:
        mesh.m3 = inObj->domain.cm3;
	maxNod = mesh.m3->res.nod.maxEnt;
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
    if((dist = AlcCalloc(maxNod, sizeof(double))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to allocate distance array (%s),\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok && (seedType == WLZTST_SEED_SEEDS))
  {
    switch(inObj->type)
    {
      case WLZ_CMESH_2D:
	if((seeds.d2 = (WlzDVertex2 *)
	               AlcMalloc(sizeof(WlzDVertex2) * inRow)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idN = 0; idN < inRow; ++idN)
	  {
	    (seeds.d2 + idN)->vtX = *(*(inSeeds + idN) + 0);
	    (seeds.d2 + idN)->vtY = *(*(inSeeds + idN) + 1);
	  }
	  nSeeds = inRow;
	}
        break;
      case WLZ_CMESH_3D:
	if((seeds.d3 = (WlzDVertex3 *)
	               AlcMalloc(sizeof(WlzDVertex3) * inRow)) == NULL)
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  for(idN = 0; idN < inRow; ++idN)
	  {
	    (seeds.d3 + idN)->vtX = *(*(inSeeds + idN) + 0);
	    (seeds.d3 + idN)->vtY = *(*(inSeeds + idN) + 1);
	    (seeds.d3 + idN)->vtZ = *(*(inSeeds + idN) + 2);
	  }
	  nSeeds = inRow;
	}
        break;
      default:
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: failed to build seeds from input (%s),\n",
		     argv[0],
		     errMsgStr);
    }
  }
  (void )Alc2Free((void **)inSeeds);
  if(ok)
  {
    for(idN = 0; idN < repeats; ++idN)
    {
      errNum = WlzCMeshFMarNodes2D(mesh.m2, dist, NULL, nSeeds, seeds.d2);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		       "%s Failed to compute distances in mesh (%s).\n",
		       argv[0],
		       errMsgStr);
        break;
      }
    }
  }
  if(ok)
  {
    switch(outType)
    {
      case WLZTST_OUT_TXT:
	if((fP = (strcmp(outFileStr, "-")?
		 fopen(outFileStr, "w"): stdout)) == NULL)
	{
	  ok = 0;
	  (void )fprintf(stderr,
			 "%s: Failed to open output file %s.\n",
			 argv[0], outFileStr);
	}
	if(ok)
	{
	  switch(mesh.m2->type)
	  {
	    case WLZ_CMESH_TRI2D:
	      for(idN = 0; idN < maxNod; ++idN)
	      {
		nod2 = (WlzCMeshNod2D *)
		       AlcVectorItemGet(mesh.m2->res.nod.vec, idN);
		if(nod2->idx >= 0)
		{
		  if(fprintf(fP, "% 8d % 8lg % 8lg % 8lg % 8g\n",
			     nod2->idx,
			     nod2->pos.vtX,
			     nod2->pos.vtY,
			     0.0,
			     dist[nod2->idx]) <= 0)
		  {
		    errNum = WLZ_ERR_WRITE_INCOMPLETE;
		    break;
		  }
		}
	      }
	      break;
	    case WLZ_CMESH_TET3D:
	      errNum = WLZ_ERR_UNIMPLEMENTED;
	      break;
	    default:
	      break;
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsgStr);
	    (void )fprintf(stderr,
			   "%s Failed to output distance data (%s).\n",
			   argv[0],
			   errMsgStr);
	    
	  }
	}
	if(fP && strcmp(outFileStr, "-"))
	{
	  (void )fclose(fP); fP = NULL;
	}
	break;
      default:
	break;
    }
  }
  (void )WlzFreeObj(inObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>]\n"
    "       [-b] [-s<seed>] [-R<repeats>] [-S<seed file>]\n"
    "       [-t] [<input cmesh object>]\n"
    "Reads a conforming mesh and then computes distances from the given\n"
    "seeds or the boundary nodes.\n"
    "The distances are either printed to the output file as text or output\n"
    "as a Woolz domain object with double distance values.\n"
    "All seeds must be specified using 3D coordinates <x>,<y>,<z> on the\n"
    "command line or in a file:\n"
    "<x> <y> <z>\n"
    "ie one seed point per line with space seperated fields (this is so\n"
    "even for 2D objects).\n"
    "Text output has one node per line with the fields:\n"
    "<node index> <x> <y> <z> <distance>\n"
    "Image output is a double valued image with the same domain as the\n"
    "input object. Values within the output image are either 0.0 or the\n"
    "computed distance at the mesh nodes.\n"
    "Options are:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file.\n"
    "  -b  Compute distances from boundary nodes.\n"
    "  -s  Compute distances from given seed point, which must be within\n"
    "      the mesh.\n"
    "  -R  number of times to repeat the computation.\n"
    "  -S  File of seed points, each as of which must be within the mesh.\n"
    "  -t  Output text data.\n",
    argv[0]);

  }
  return(!ok);
}
