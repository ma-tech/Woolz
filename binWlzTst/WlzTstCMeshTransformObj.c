#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBasisFnTransformObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file		binWlzTst/WlzTstCMeshTransformObj.c
* \author       Bill Hill
* \date         February 2005
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
* \brief	Applies a constrained mesh transform to a Woolz object.
* \ingroup	BinWlzTst
* \todo         -
* \bug          None known.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
static WlzErrorNum 		WlzTstCMeshSetDsp(
				  WlzCMeshTransform *cMesh);
static WlzErrorNum 		WlzTstCMeshSetDsp3D(
				  WlzCMeshTransform *cMesh);
static WlzDVertex3 		WlzTstCMeshCompDsp3D(
				  WlzDBox3 bBox,
				  WlzDVertex3 pos);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
  		deform = 1,
		ok = 1,
		usage = 0;
  double	minElmSz = 4.0,
  		maxElmSz = 8.0;
  WlzObject	*inObj = NULL,
		*outObj = NULL,
		*dilObj = NULL;
  WlzCMeshTransform *cMesh = NULL;
  WlzMeshGenMethod meshGenMth = WLZ_MESH_GENMETHOD_CONFORM;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*inObjFileStr,
  		*outObjFileStr,
		*outVTKFileStr = NULL;
  const char    *errMsg;
  static char	optList[] = "hIm:M:o:V:",
  		inObjFileStrDef[] = "-",
		outObjFileStrDef[] = "-";

  opterr = 0;
  inObjFileStr = inObjFileStrDef;
  outObjFileStr = outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
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
        outObjFileStr = optarg;
	break;
      case 'I':
        deform = 0;
	break;
      case 'V':
        outVTKFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
	(*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP)
    {
      if(strcmp(inObjFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    cMesh = WlzCMeshTransformFromObj(inObj,
			  meshGenMth, minElmSz, maxElmSz,
			  &dilObj, 1, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to compute transform from object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok && outVTKFileStr)
  {
    if(((fP = (strcmp(outVTKFileStr, "-")?
	      fopen(outVTKFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to write mesh to file.\n",
		     *argv);
    }
    if(ok)
    {
      WlzCMeshDbgOutVTK(fP, cMesh->mesh);
    }
    if(fP && strcmp(outVTKFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok && deform)
  {
    errNum = WlzTstCMeshSetDsp(cMesh);
  }
  if(ok)
  {
    outObj = WlzCMeshTransformObj(inObj, cMesh, interp, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to transform object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
	      fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to write out transformed object (%s).\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  (void )WlzFreeCMeshTransform(cMesh);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(dilObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s\n%s",
    *argv,
    " [-h] [-o<out object>] [-I] [-m #] [-M #] [-V <out VTK file>\n"
    "                  [<in object>]\n"
    "Options:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output object file name.\n"
    "  -I  Use identity transform, ie displacemets all zero.\n"
    "  -m  Minimum mesh element size.\n"
    "  -M  Maximum mesh element size.\n"
    "  -V  Output VTK mesh file name.\n"
    "Applies a constrained mesh transform to the given object.\n");
  }
  return(!ok);
}

static WlzErrorNum WlzTstCMeshSetDsp(WlzCMeshTransform *cMesh)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(cMesh->type)
  {
    case WLZ_TRANSFORM_2D_CMESH:
      errNum = WLZ_ERR_UNIMPLEMENTED;
      break;
    case WLZ_TRANSFORM_3D_CMESH:
      errNum = WlzTstCMeshSetDsp3D(cMesh);
      break;
    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
  }
  return(errNum);
}

static WlzErrorNum WlzTstCMeshSetDsp3D(WlzCMeshTransform *cMesh)
{
  int		idN,
  		nNod;
  WlzDBox3	bBox;
  AlcVector	*dVec,
  		*nVec;
  WlzCMeshNod3D	*nod;
  WlzDVertex3	*dsp;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  bBox = cMesh->mesh.m3->bBox;
  nNod = cMesh->mesh.m3->res.nod.maxEnt;
  dVec = cMesh->dspVec;
  nVec = cMesh->mesh.m3->res.nod.vec;
  for(idN = 0; idN < nNod; ++idN)
  {
    nod = (WlzCMeshNod3D *)AlcVectorItemGet(nVec, (size_t )idN);
    if(nod->idx >= 0)
    {
      dsp = (WlzDVertex3 *)AlcVectorItemGet(dVec, nod->idx);
      *dsp = WlzTstCMeshCompDsp3D(bBox, nod->pos);
    }
  }

  return(errNum);
}

static WlzDVertex3 WlzTstCMeshCompDsp3D(WlzDBox3 bBox, WlzDVertex3 pos)
{
  double	x,
  		y,
		z;
  WlzDVertex3	dsp;

  x = (0.5 * (bBox.xMin + bBox.xMax) - pos.vtX) / (bBox.xMax - bBox.xMin);
  y = (0.5 * (bBox.yMin + bBox.yMax) - pos.vtY) / (bBox.yMax - bBox.yMin);
  z = (0.5 * (bBox.zMin + bBox.zMax) - pos.vtZ) / (bBox.zMax - bBox.zMin);
  dsp.vtX = 0.3333 * (bBox.xMax - bBox.xMin) * sin(WLZ_M_PI * x * x);
  dsp.vtY = 0.3333 * (bBox.yMax - bBox.yMin) * sin(WLZ_M_PI * y * y);
  dsp.vtZ = 0.3333 * (bBox.zMax - bBox.zMin) * sin(WLZ_M_PI * z * z);
  return(dsp);
}
