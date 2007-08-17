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
* \brief	Test for 2D and 3D conforming simplical mesh generation.
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


extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idE,
  		nElm = 0,
  		ok = 1,
  		option,
  		usage = 0,
		features = 0,
		laplacianItr = 0,
		lowPassItr = 0,
		smoothBnd = 0,
		verify = 0;
  double	laplacianAlpha = 0.1,
		lowPassLambda = 0.33,
		lowPassMu = -0.34,
  		minElmSz = 25.0,
  		maxElmSz = 100.0;
  int		*idx = NULL;
  double	*vol = NULL,
  		*minLen = NULL,
		*maxLen = NULL;
  FILE		*fP = NULL;
  char		*inObjFileStr,
  		*outFileStr;
  const char	*errMsgStr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*obj = NULL;
  WlzCMeshP 	mesh;
  static char   optList[] = "a:BFhl:L:m:M:o:u:VW:";
  const char    inObjFileStrDef[] = "-",
  	        outFileStrDef[] = "-";

  mesh.v = NULL;
  opterr = 0;
  inObjFileStr = (char *)inObjFileStrDef;
  outFileStr = (char *)outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'a':
	if(sscanf(optarg, "%lg", &laplacianAlpha) != 1)
	{
	  usage = 1;
	}
        break;
      case 'B':
        smoothBnd = 1;
	break;
      case 'F':
        features = 1;
	break;
      case 'l':
	if(sscanf(optarg, "%lg", &lowPassLambda) != 1)
	{
	  usage = 1;
	}
        break;
      case 'L':
	if(sscanf(optarg, "%d", &laplacianItr) != 1)
	{
	  usage = 1;
	}
        break;
      case 'm':
        if(sscanf(optarg, "%lg", &minElmSz) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
        if(sscanf(optarg, "%lg", &maxElmSz) != 1)
	{
	  usage = 1;
	}
        break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'u':
	if(sscanf(optarg, "%lg", &lowPassMu) != 1)
	{
	  usage = 1;
	}
        break;
      case 'W':
	if(sscanf(optarg, "%d", &lowPassItr) != 1)
	{
	  usage = 1;
	}
        break;
      case 'V':
        verify = 1;
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
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (outFileStr == NULL) || (*outFileStr == '\0'))
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
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
              fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((obj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
       (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: failed to read object from file %s\n",
                     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  if(ok)
  {
    (void )WlzAssignObject(obj, NULL);
    mesh = WlzCMeshFromObj(obj, minElmSz, maxElmSz, NULL, 1, &errNum);
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
    WlzCMeshSetBoundNodFlags(mesh);
    if(verify)
    {
      errNum = WlzCMeshVerify(mesh, NULL, 1, stderr);
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsgStr);
	(void )fprintf(stderr,
		       "%s Failed to verify mesh, %s.\n",
		       argv[0],
		       errMsgStr);
      }
    }
  }
  if(ok && (laplacianItr > 0))
  {
    errNum = WlzCMeshLaplacianSmooth(mesh, laplacianItr, laplacianAlpha,
    				     smoothBnd, 1);
    if((errNum == WLZ_ERR_NONE) && verify)
    {
      errNum = WlzCMeshVerify(mesh, NULL, 1, stderr);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to Laplacian smooth mesh, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok && (lowPassItr > 0))
  {
    errNum = WlzCMeshLPFilterLM(mesh, lowPassLambda, lowPassMu,
    				lowPassItr, smoothBnd, 1);
    if((errNum == WLZ_ERR_NONE) && verify)
    {
      errNum = WlzCMeshVerify(mesh, NULL, 1, stderr);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to low pass filter mesh, %s.\n",
		     argv[0],
		     errMsgStr);
    }
  }
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
	     fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to open output file %s.\n",
		     argv[0], outFileStr);
    }
  }
  if(ok && features)
  {
    errNum = WlzCMeshCmpElmFeat(mesh, &nElm, &idx, &vol, &minLen, &maxLen);
    if(errNum == WLZ_ERR_NONE)
    {
      for(idE = 0; idE < nElm; ++idE)
      {
        (void )printf("%d %lg %lg %lg\n",
	              *(idx + idE), *(vol + idE),
	              *(minLen + idE), *(maxLen + idE));
      }
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s Failed to compute mesh features, %s.\n",
		     argv[0],
		     errMsgStr);
    }
    AlcFree(idx);
    AlcFree(vol);
    AlcFree(minLen);
    AlcFree(maxLen);
  }
  if(ok)
  {
    WlzCMeshDbgOutVTK(fP, mesh);
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    fprintf(stderr,
            "Usage: %s [-h] [-o<output file>]\n"
	    "       [-L#] [-a#] [-W#] [-l#] [-u#] [-B]\n"
	    "       [-m#] [-M#] [-F] [-V] [<input object>]\n"
    	    "Computes a conforming mesh for the given input object.\n"
	    "Options are:\n"
	    "  -h  Help, prints this usage message.\n"
	    "  -o  Output file.\n"
	    "  -L  Number of Laplacian smoothing iterations.\n"
	    "  -a  Laplacian alpha parameter.\n"
	    "  -W  Number of low pass filter smoothing iterations.\n"
	    "  -l  Low pass filter lambda value.\n"
	    "  -u  Low pass filter mu value.\n"
	    "  -B  Smooth boundary (requires a smoothing method to be\n"
	    "      selected).\n"
	    "  -m  Minimum mesh element size.\n"
	    "  -M  Maximum mesh element size.\n"
	    "  -F  Prints features of the mesh elements to the standard\n"
	    "      output.\n"
	    "  -V  Verify mesh. This may take a long time and may give\n"
	    "      segmentation faults for buggy meshes.\n",
	    argv[0]);

  }
  return(!ok);
}
