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
* \file		binWlz/WlzBasisFnTransformObj.c
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
* \brief	Computes and applies Woolz basis function transforms.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzbasisfntransformobj "WlzBasisFnTransformObj"

*/

/*!
\ingroup BinWlz
\defgroup wlzbasisfntransformobj WlzBasisFnTransformObj
\par Name
WlzBasisFnTransformObj - computes and applies Woolz basis function transforms.
\par Synopsis
\verbatim
WlzBasisFnTransformObj [-o<out object>] [-p<tie points file>]
		       [-m<min mesh dist>] [-M<max mesh dist>]
		       [-b<basis fn transform>] [-Y<order of polynomial>]
		       [-D<flags>]
		       [-d] [-g] [-h] [-q] [-Q] [-s] [-t] [-y]
		       [-B] [-C] [-G] [-L] [-N] [-T]
		       [<in object>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-B</b></td>
    <td>Block mesh generation method (default).</td>
  </tr>
  <tr> 
    <td><b>-C</b></td>
    <td>Use conforming mesh.</td>
  </tr>
  <tr> 
    <td><b>-D</b></td>
    <td>Debug flags:
    <table width="500" border="0">
    <tr>
    <td>1</td>
    <td>Output the mesh as postscript to standard error output.</td>
    </tr>
    <tr>
    <td>2</td>
    <td>Output basis function distance maps to standard error
        output.</td>
    </tr>
    </table>
    </tr>
    </td>
  </tr>
  <tr> 
    <td><b>-G</b></td>
    <td>Gradient mesh generation method.</td>
  </tr>
  <tr> 
    <td><b>-L</b></td>
    <td>Use linear interpolation instead of nearest neighbour.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Minimum mesh node separation distance (default 10.0).</td>
  </tr>
  <tr> 
    <td><b>-M</b></td>
    <td>Maximum mesh node separation distance (default 100.0).</td>
  </tr>
  <tr> 
    <td><b>-N</b></td>
    <td>Don't transform the object (useful for testing).</td>
  </tr>
  <tr> 
    <td><b>-R</b></td>
    <td>Restrict the transformation to only a basis function transform
        instead of the default which uses a least squares affine and
	basis function transform.
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Basis function transform object.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Tie point file.</td>
  </tr>
  <tr> 
    <td><b>-c</b></td>
    <td>Use conformal polynomial basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Compute transform using conforming distance transforms, also
        sets the mesh generation to produce a conforming mesh and
	restricts the transformation to a basis function transform
	only (ie no least squares affine).</td>
  </tr>
  <tr> 
    <td><b>-g</b></td>
    <td>Use Gaussian basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-q</b></td>
    <td>Use multi-quadric basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-Q</b></td>
    <td>Use inverse-multi-quadric basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Use thin plate spline basis function (default) if tie points
        are given.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Target mesh. Only valid for conforming meshes. If given then the
        target mesh is used to compute a target-to-source radial basis
	function, with distances computed in the target mesh. This is
	then inverted as a mesh transform is created.
  </tr>
  <tr> 
    <td><b>-T</b></td>
    <td>Output a basis function transform instead of a transformed object.</td>
  </tr>
  <tr> 
    <td><b>-U</b></td>
    <td>Output a mesh transform instead of a transformed object.</td>
  </tr>
  <tr> 
    <td><b>-y</b></td>
    <td>Use polynomianl basis function if tie points are given.</td>
  </tr>
  <tr> 
    <td><b>-Y</b></td>
    <td>Polynomial order for polynomial basis function (default 3).</td>
  </tr>
</table>
\par Description
Computes and applies Woolz basis function transforms.
Tie points may be read from an ascii file with the format:
\verbatim
  <vertex x> <vertex y> <displacement x> <displacement y>
\endverbatim

\par Examples
\verbatim
WlzBasisFnTransformObj -o tied.wlz -p points.tie myobj.wlz
\endverbatim
A thin plate spline basis function transform is computed from the
tie points read from the file points.tie. This transform is then
applied to the object is read from myobj.wlz. The resulting object
is then written to tied.wlz.
\par File
\ref WlzBasisFnTransformObj.c "WlzBasisFnTransformObj.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref wlzaffinetransformobj "WlzAffineTransformObj(1)"
\ref wlzmeshtransformobj "WlzMeshTransformObj(1)"
\ref WlzBasisFnTransformObj "WlzBasisFnTransformObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>

#define	IN_RECORD_MAX   (1024)

static void			WlzMeshOutputPS(
				  WlzMeshTransform *meshTr);
static void			WlzCMeshOutputPS(
				  WlzObject *mObj,
				  WlzObject *sObj,
				  WlzObject *dilObj,
				  WlzObject *outObj);
static void			WlzBndOutputPS(
				  WlzBoundList *bnd,
			          double scale,
				  WlzDVertex2 offset);
static void			WlzPolyOutputPS(
				  WlzPolygonDomain *poly,
			          double scale,
				  WlzDVertex2 offset);
static WlzErrorNum 		WlzBFTOGetVertices2D(
				  int *dstNVx,
				  WlzDVertex2 **dstVx0,
				  WlzDVertex2 **dstVx1,
				  FILE *fP,
				  WlzCMesh2D *tMesh,
				  const char *prog,
				  const char *fStr);
static WlzErrorNum 		WlzBFTOGetVertices3D(
				  int *dstNVx,
				  WlzDVertex3 **dstVx0,
				  WlzDVertex3 **dstVx1,
				  FILE *fP,
				  WlzCMesh3D *tMesh,
				  const char *prog,
				  const char *fStr);

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
  		nTiePP = 0,
		dim = 0,                /* Dimension will be set by objects. */
		cdt = 0,
		cMesh = 0,
		basisFnPolyOrder = 3,
		dbgFlg = 0, 		  /* Bit wise debug flags see usage. */
		noTrObj = 0,
		outBasisTrFlag = 0,
		outMeshTrFlag = 0,
		restrictToBasisFn = 0,
		ok = 1,
		usage = 0;
  double	meshMinDist = 20.0,
  		meshMaxDist = 40.0;
  WlzVertexP	vxA0,
  		vxA1;
  WlzObject	*inObj = NULL,
		*outObj = NULL,
		*dilObj = NULL,
           	*tmpObj = NULL;
  WlzCMeshP	tarMesh;
  WlzTransform  meshTr;
  WlzMeshGenMethod meshGenMth = WLZ_MESH_GENMETHOD_GRADIENT;
  WlzBasisFnTransform *basisTr = NULL;
  WlzFnType basisFnType = WLZ_FN_BASIS_2DMQ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*inObjFileStr,
		*basisFnTrFileStr = NULL,
		*tarMeshFileStr = NULL,
		*tiePtFileStr = NULL,
  		*outObjFileStr;
  const int	delOut = 1;
  const char    *errMsg;
  static char	optList[] = "b:m:o:p:t:D:M:Y:cdghqsyBCGLNQRTU",
  		inObjFileStrDef[] = "-",
		outObjFileStrDef[] = "-";

  opterr = 0;
  vxA0.v = NULL;
  vxA1.v = NULL;
  tarMesh.v = NULL;
  meshTr.core = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'b':
        basisFnTrFileStr = optarg;
        break;
      case 'B':
        meshGenMth = WLZ_MESH_GENMETHOD_BLOCK;
	break;
      case 'C':
	meshGenMth = WLZ_MESH_GENMETHOD_CONFORM;
	break;
      case 'D':
	if(1 != sscanf(optarg, "%d", &dbgFlg))
	{
	  usage = 1;
	}
	break;
      case 'c':
        basisFnType = WLZ_FN_BASIS_2DCONF_POLY;
	break;
      case 'd':
        cdt = 1;
	break;
      case 'g':
        basisFnType = WLZ_FN_BASIS_2DGAUSS;
	break;
      case 'G':
        meshGenMth = WLZ_MESH_GENMETHOD_GRADIENT;
	break;
      case 'L':
        interp = WLZ_INTERPOLATION_LINEAR;
	break;
      case 'm':
	if(sscanf(optarg, "%lg", &meshMinDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'M':
	if(sscanf(optarg, "%lg", &meshMaxDist) != 1)
	{
	  usage = 1;
	}
	break;
      case 'N':
        noTrObj = 1;
	break;
      case 'p':
        tiePtFileStr = optarg;
	break;
      case 'q':
        basisFnType = WLZ_FN_BASIS_2DMQ;
	break;
      case 'Q':
        basisFnType = WLZ_FN_BASIS_2DIMQ;
	break;
      case 'R':
        restrictToBasisFn = 1;
	break;
      case 's':
        basisFnType = WLZ_FN_BASIS_2DTPS;
	break;
      case 't':
        tarMeshFileStr = optarg;
	break;
      case 'T':
        outBasisTrFlag = 1;
	break;
      case 'y':
        basisFnType = WLZ_FN_BASIS_2DPOLY;
	break;
      case 'U':
        outMeshTrFlag = 1;
	break;
      case 'Y':
        if(sscanf(optarg, "%d", &basisFnPolyOrder) != 1)
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
    if(cdt)
    {
      cMesh = 1;
      restrictToBasisFn = 1;
      meshGenMth = WLZ_MESH_GENMETHOD_CONFORM;
    }
    if((tarMeshFileStr != NULL) && (cMesh == 0))
    {
      usage = 1;
    }
  }
  if(usage == 0)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
       (((tiePtFileStr == NULL) || (*tiePtFileStr == '\0')) &&
	((basisFnTrFileStr == NULL) || (*basisFnTrFileStr == '\0'))) ||
       (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      usage = 1;
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
  /* Read input object that's to be transformed. */
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
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
	fclose(fP);
      }
      fP = NULL;
    }
    if(ok)
    {
      switch(inObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  dim = 2;
	  break;
	case WLZ_3D_DOMAINOBJ:
	  dim = 3;
	  break;
        default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
      }
    }
    if(ok == 0)
    {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: Failed to read object from file %s (%s)\n",
		     *argv, inObjFileStr, errMsg);
    }
  }
  /* Read target mesh if required. */
  if(ok)
  {
    if(tarMeshFileStr != NULL)
    {
      errNum = WLZ_ERR_READ_EOF;
      if((tarMeshFileStr == NULL) ||
	  (*tarMeshFileStr == '\0') ||
	  ((fP = (strcmp(tarMeshFileStr, "-")?
		  fopen(tarMeshFileStr, "r"): stdin)) == NULL) ||
	  ((tmpObj = WlzAssignObject(WlzReadObj(fP,
	                                        &errNum), NULL)) == NULL) ||
	  (errNum != WLZ_ERR_NONE))
      {
	ok = 0;
      }
      else
      {
        switch(tmpObj->type)
	{
	  case WLZ_CMESH_2D:
	    if(dim == 0)
	    {
	      dim = 2;
	    }
	    if(dim != 2)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    else
	    {
	      if(tmpObj->domain.core == NULL)
	      {
		errNum = WLZ_ERR_DOMAIN_NULL;
	      }
	      else if(tmpObj->domain.cm2->type != WLZ_CMESH_TRI2D)
	      {
		errNum = WLZ_ERR_DOMAIN_DATA;
	      }
	      else
	      {
		tarMesh.m2 = tmpObj->domain.cm2;
		tmpObj->domain.cm2 = NULL;
	      }
	    }
	    break;
	  case WLZ_CMESH_3D:
	    if(dim == 0)
	    {
	      dim = 3;
	    }
	    if(dim != 3)
	    {
	      errNum = WLZ_ERR_OBJECT_TYPE;
	    }
	    else
	    {
	      if(tmpObj->domain.core == NULL)
	      {
		errNum = WLZ_ERR_DOMAIN_NULL;
	      }
	      else if(tmpObj->domain.cm3->type != WLZ_CMESH_TET3D)
	      {
		errNum = WLZ_ERR_DOMAIN_DATA;
	      }
	      else
	      {
		tarMesh.m3 = tmpObj->domain.cm3;
		tmpObj->domain.cm3 = NULL;
	      }
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
        if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	}
      }
      (void )WlzFreeObj(tmpObj);
      tmpObj = NULL;
      if(ok == 0)
      {
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		 "%s: Failed to read target mesh object from file %s (%s).\n",
		       *argv, tarMeshFileStr, errMsg);
      }
    }
  }
  /* Either read the basis function transform from a file or compute it
   * from a set of tie points. */
  if(ok)
  {
    if(basisFnTrFileStr)
    {
      errNum = WLZ_ERR_UNSPECIFIED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Reading basis function transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else if(tiePtFileStr)
    {
      if((fP = (strcmp(tiePtFileStr, "-")?
                fopen(tiePtFileStr, "r"): stdin)) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_READ_EOF;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to open tie points file %s (%s).\n",
		       *argv, tiePtFileStr, errMsg);
      }
      else
      {
        if(dim == 2)
	{
	  errNum = WlzBFTOGetVertices2D(&nTiePP, &(vxA0.d2), &(vxA1.d2),
	  				fP, tarMesh.m2, argv[0], tiePtFileStr);
	}
	else /* dim == 3 */
	{
	  errNum = WlzBFTOGetVertices3D(&nTiePP, &(vxA0.d3), &(vxA1.d3),
	  				fP, tarMesh.m3, argv[0], tiePtFileStr);
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	}
      }
      if(fP && strcmp(tiePtFileStr, "-"))
      {
	fclose(fP);
      }
      fP = NULL;
    }
  }
  if(ok)
  {
    if(dim == 3)
    {
      /* Promote basis function type from 2D to 3D. */
      switch(basisFnType)
      {
        case WLZ_FN_BASIS_2DIMQ:
          basisFnType = WLZ_FN_BASIS_3DIMQ;
	  break;
        case WLZ_FN_BASIS_2DMQ:
          basisFnType = WLZ_FN_BASIS_3DMQ;
	  break;
        default:
	  errNum = WLZ_ERR_DOMAIN_TYPE;
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	                 "%s: invalid basis function type for 3D (%s).\n",
			 *argv, errMsg);
      }
    }
  }
  if(ok)
  {
    if(restrictToBasisFn)
    {
      if(cMesh)
      {
	if(tarMesh.v == NULL)
	{
	  meshTr.obj = WlzCMeshTransformFromObj(inObj,
				meshGenMth, meshMinDist, meshMaxDist,
				&dilObj, delOut, &errNum);
        }
      }
      else
      {
	meshTr.mesh = WlzMeshFromObj(inObj,
			      meshGenMth, meshMinDist, meshMaxDist,
			      &errNum);

      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to compute a mesh for the object (%s).\n",
		       *argv, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(tarMesh.v == NULL)
	{
	  if(dim == 2)
	  {
	    basisTr = WlzBasisFnTrFromCPts2D(basisFnType, basisFnPolyOrder,
					   nTiePP, vxA0.d2, nTiePP, vxA1.d2,
					   (cdt)? meshTr.obj->domain.cm2: NULL,
					   &errNum);
	  }
	  else /* dim == 3 */
	  {
	    basisTr = WlzBasisFnTrFromCPts3D(basisFnType, basisFnPolyOrder,
					   nTiePP, vxA0.d3, nTiePP, vxA1.d3,
					   (cdt)? meshTr.obj->domain.cm3: NULL,
					   &errNum);
	  }
	}
	else
	{
	  if(dim == 2)
	  {
	    basisTr = WlzBasisFnTrFromCPts2D(basisFnType, basisFnPolyOrder,
	    				     nTiePP, vxA1.d2, nTiePP, vxA0.d2,
					     tarMesh.m2, &errNum);
	  }
	  else /* dim == 3 */
	  {
	    basisTr = WlzBasisFnTrFromCPts3D(basisFnType, basisFnPolyOrder,
	    				     nTiePP, vxA1.d3, nTiePP, vxA0.d3,
					     tarMesh.m3, &errNum);
	  }
	}
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
		   "%s: Failed to compute basis function transform (%s).\n",
		       *argv, errMsg);
	}
      }
    }
    else /* Use affine and basis functions to define mesh. */
    {
      meshTr.mesh = WlzMeshTransformFromCPts(inObj,
      				basisFnType, basisFnPolyOrder,
      				nTiePP, vxA0.d2, nTiePP, vxA1.d2,
				meshGenMth, meshMinDist, meshMaxDist,
				&errNum);
    }
  }
  if(ok)
  {
    if(outBasisTrFlag)
    {
      errNum = WLZ_ERR_UNIMPLEMENTED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Writing basis function transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else /* Transform the object and then write it out. */
    {
      if(restrictToBasisFn)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  if(cMesh)
	  {
	    if(tarMesh.v == NULL)
	    {
	      errNum = WlzBasisFnSetCMesh(meshTr.obj, basisTr);
	    }
	    else
	    {
	       meshTr.obj = WlzBasisFnInvertMakeCMeshTr(basisTr, tarMesh,
	       				&errNum);
	    }
	  }
	  else
	  {
	    errNum = WlzBasisFnSetMesh(meshTr.mesh, basisTr);
	  }
	  if(errNum != WLZ_ERR_NONE)
	  {
	    ok = 0;
	    (void )WlzStringFromErrorNum(errNum, &errMsg);
	    (void )fprintf(stderr,
			   "%s: Failed to set the mesh displacements (%s).\n",
			   *argv, errMsg);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(outMeshTrFlag)
	{
	  if(cMesh)
	  {
	    errNum = WLZ_ERR_WRITE_EOF;
	    if(((fP = (strcmp(outObjFileStr, "-")?
		      fopen(outObjFileStr, "w"):
		      stdout)) == NULL) ||
	       ((errNum = WlzWriteObj(fP, meshTr.obj)) != WLZ_ERR_NONE))
	    {
	      ok = 0;
	      (void )WlzStringFromErrorNum(errNum, &errMsg);
	      (void )fprintf(stderr,
			     "%s: Failed to write CMesh object (%s).\n",
			     *argv, errMsg);
	    }
	    if(fP && strcmp(outObjFileStr, "-"))
	    {
	      fclose(fP);
	    }
	  }
	  else
	  {
	    errNum = WLZ_ERR_UNIMPLEMENTED;
	    ok = 0;
	    (void )fprintf(stderr,
			   "%s: Writing mesh transforms has not been\n"
			   "implemented yet\n",
			   *argv);
	  }
	}
	else
	{
	  if(noTrObj == 0)
	  {
	    if(cMesh)
	    {
	      outObj = WlzCMeshTransformObj(inObj, meshTr.obj, interp,
					    &errNum);
	    }
	    else
	    {
	      outObj = WlzMeshTransformObj(inObj, meshTr.mesh, interp,
					   &errNum);
	    }
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(dbgFlg & 1)
	{
	  if(cMesh)
	  {
	    WlzCMeshOutputPS(meshTr.obj, inObj, dilObj, outObj);
	  }
	  else
	  {
	    WlzMeshOutputPS(meshTr.mesh);
	  }
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
	ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
		       "%s: Failed to transform object (%s).\n",
		       *argv, errMsg);
      }
      if((errNum == WLZ_ERR_NONE) && (noTrObj == 0) &&
         (outBasisTrFlag == 0) && (outMeshTrFlag == 0))
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
			 "%s: Failed to write output object (%s).\n",
			 *argv, errMsg);
	}
	if(fP && strcmp(outObjFileStr, "-"))
	{
	  fclose(fP);
	}
      }
    }
  }
  AlcFree(vxA0.v);
  AlcFree(vxA1.v);
  (void )WlzCMeshFree(tarMesh);
  if(cMesh)
  {
    (void )WlzFreeObj(meshTr.obj);
  }
  else
  {
    (void )WlzMeshFreeTransform(meshTr.mesh);
  }
  (void )WlzBasisFnFreeTransform(basisTr);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(dilObj);
  if(usage != 0)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-p<tie points file>]\n"
    "                  [-m<min mesh dist>] [-M<max mesh dist>]\n"
    "                  [-b<basis fn transform>] [-Y<order of polynomial>]\n"
    "                  [-D<flags>]\n"
    "                  [-d] [-g] [-h] [-q] [-s] [-t] [-y]\n"
    "                  [-B] [-C] [-G] [-L] [-N] [-T]\n"
    "                  [<in object>]\n"
    "Options:\n"
    "  -b  Basis function transform object.\n"
    "  -B  Block mesh generation method (default).\n"
    "  -C  Use conforming mesh.\n"
    "  -D  Debug flags:\n"
    "        1  Output the mesh as postscript to standard error output.\n"
    "        2  Output basis function distance maps to standard error\n"
    "           output.\n"
    "      These debug flags are only intended for use when debuging and\n"
    "      they may be combined by an or operation (eg 11 = 1 | 2 | 8).\n"
    "  -G  Gradient mesh generation method.\n"
    "  -L  Use linear interpolation instead of nearest neighbour.\n"
    "  -m  Minimum mesh node separation distance (default 10.0)\n"
    "  -M  Maximum mesh node separation distance (default 100.0)\n"
    "  -N  Don't transform the object (useful for testing)\n"
    "  -R  Restrict the transformation to only a basis function transform\n"
    "      instead of the default which uses a least squares affine and\n"
    "      basis function transform.\n"
    "  -o  Output object file name.\n"
    "  -p  Tie point file.\n"
    "  -c  Use conformal polynomial basis function if tie points are given.\n"
    "  -d  Compute transform using conforming distance transforms, also\n"
    "      sets the mesh generation to produce a conforming mesh and\n"
    "      restricts the transformation to a basis function transform\n"
    "      only (ie no least squares affine).\n"
    "  -g  Use Gaussian basis function if tie points are given.\n"
    "  -h  Help, prints this usage message.\n"
    "  -q  Use multi-quadric basis function if tie points are given.\n"
    "  -Q  Use inverse-multi-quadric basis function if tie points are given.\n"
    "  -s  Use thin plate spline basis function (default) if tie points\n"
    "      are given.\n"
    "  -t  Target mesh. Only valid for conforming meshes. If given then the\n"
    "      target mesh is used to compute a target-to-source radial basis\n"
    "      function, with distances computed in the target mesh. This is\n"
    "      then inverted as a mesh transform is created.\n"
    "  -T  Output a basis function transform instead of a transformed\n"
    "      object.\n"
    "  -U  Output a mesh transform instead of a transformed object.\n"
    "  -y  Use polynomial basis function if tie points are given.\n"
    "  -Y  Polynomial order for polynomial basis function (default 3).\n"
    "Computes and applies Woolz basis function transforms.\n"
    "Tie points may be read from an ascii file with the 2D format:\n"
    "  <x> <y> <displacement x> <displacement y>\n"
    "and the 3D format:\n"
    "  <x> <y> <z> <displacement x> <displacement y> <displacement z>\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o tied.wlz -p points.tie myobj.wlz\n"
    "A thin plate spline basis function transform is computed from the\n"
    "tie points read from the file points.tie. This transform is then\n"
    "applied to the object is read from myobj.wlz. The resulting object\n"
    "is then written to tied.wlz.\n");
  }
  return(!ok);
}


/*!
* \return	Woolz error code.
* \brief	Reads 2D tie points from the file and checks that they
* 		are within the mesh if the mesh is given.
* \param	dstNVx			Destination pointer for number of
* 					tie point pairs, must not be NULL.
* \param	dstVx0			Destination pointer for source
* 					vertices, must not be NULL.
* \param	dstVx1			Destination pointer for displacement
* 					vertices, must not be NULL.
* \param	fP			Input file.
* \param	tMesh			Target mesh, may be NULL.
* \param	prog			Program name for error output, must
* 					not be NULL.
* \param	fStr			input file name for error output, must
* 					not be NULL.
*/
static WlzErrorNum WlzBFTOGetVertices2D(int *dstNVx,
				WlzDVertex2 **dstVx0, WlzDVertex2 **dstVx1,
				FILE *fP, WlzCMesh2D *tMesh,
				const char *prog, const char *fStr)
{
  int		vxCount = 0,
  		vxLimit = 0;
  char 		*rec;
  WlzDVertex2	*vx0,
  		*vx1,
		*vxA0 = NULL,
		*vxA1 = NULL;
  const char    *errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char  	inRecord[IN_RECORD_MAX];

  while((errNum == WLZ_ERR_NONE) &&
        (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
  {
    inRecord[IN_RECORD_MAX - 1] = '\0';
    rec = inRecord;
    while(*rec && isspace(*rec))
    {
      ++rec;
    }
    if(*rec && (*rec != '#'))
    {
      if(vxCount >= vxLimit)
      {
	vxLimit += 4096;
	if(((vxA0 = (WlzDVertex2 *)
		    AlcRealloc(vxA0,
			       vxLimit * sizeof(WlzDVertex2))) == NULL) ||
	    ((vxA1 = (WlzDVertex2 *)
		     AlcRealloc(vxA1,
				vxLimit * sizeof(WlzDVertex2))) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  vx0 = vxA0 + vxCount;
	  vx1 = vxA1 + vxCount;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(sscanf(rec, "%lg %lg %lg %lg", &(vx0->vtX), &(vx0->vtY),
	      &(vx1->vtX), &(vx1->vtY)) != 4)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  vx1->vtX += vx0->vtX;
	  vx1->vtY += vx0->vtY;
	}
      }
      if((errNum == WLZ_ERR_NONE) && (tMesh != NULL))
      {
	if(WlzCMeshElmEnclosingPos2D(tMesh, -1, vx1->vtX, vx1->vtY,
	                             0, NULL) < 0) 
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	      "%s: tie points line %d not in target mesh (%s).\n",
	      prog, vxCount, errMsg);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	++vx0;
	++vx1;
	++vxCount;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNVx = vxCount;
    *dstVx0 = vxA0;
    *dstVx1 = vxA1;
  }
  else
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
	"%s: Failed to read tie points file %s (%s).\n",
	prog, fStr, errMsg);
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \brief	Reads 3D tie points from the file and checks that they
* 		are within the mesh if the mesh is given.
* \param	dstNVx			Destination pointer for number of
* 					tie point pairs, must not be NULL.
* \param	dstVx0			Destination pointer for source
* 					vertices, must not be NULL.
* \param	dstVx1			Destination pointer for displacement
* 					vertices, must not be NULL.
* \param	fP			Input file.
* \param	tMesh			Target mesh, may be NULL.
* \param	prog			Program name for error output, must
* 					not be NULL.
* \param	fStr			input file name for error output, must
* 					not be NULL.
*/
static WlzErrorNum WlzBFTOGetVertices3D(int *dstNVx,
				WlzDVertex3 **dstVx0, WlzDVertex3 **dstVx1,
				FILE *fP, WlzCMesh3D *tMesh,
				const char *prog, const char *fStr)
{
  int		vxCount = 0,
  		vxLimit = 0;
  char 		*rec;
  WlzDVertex3	*vx0,
  		*vx1,
		*vxA0 = NULL,
		*vxA1 = NULL;
  const char    *errMsg;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char  	inRecord[IN_RECORD_MAX];

  while((errNum == WLZ_ERR_NONE) &&
        (fgets(inRecord, IN_RECORD_MAX - 1, fP) != NULL))
  {
    inRecord[IN_RECORD_MAX - 1] = '\0';
    rec = inRecord;
    while(*rec && isspace(*rec))
    {
      ++rec;
    }
    if(*rec && (*rec != '#'))
    {
      if(vxCount >= vxLimit)
      {
	vxLimit += 4096;
	if(((vxA0 = (WlzDVertex3 *)
		    AlcRealloc(vxA0,
			       vxLimit * sizeof(WlzDVertex3))) == NULL) ||
	    ((vxA1 = (WlzDVertex3 *)
		     AlcRealloc(vxA1,
				vxLimit * sizeof(WlzDVertex3))) == NULL))
	{
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
	else
	{
	  vx0 = vxA0 + vxCount;
	  vx1 = vxA1 + vxCount;
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(sscanf(rec, "%lg %lg %lg %lg %lg %lg",
	          &(vx0->vtX), &(vx0->vtY), &(vx0->vtZ),
	          &(vx1->vtX), &(vx1->vtY), &(vx1->vtZ)) != 6)
	{
	  errNum = WLZ_ERR_READ_INCOMPLETE;
	}
	else
	{
	  vx1->vtX += vx0->vtX;
	  vx1->vtY += vx0->vtY;
	  vx1->vtZ += vx0->vtZ;
	}
      }
      if((errNum == WLZ_ERR_NONE) && (tMesh != NULL))
      {
	if(WlzCMeshElmEnclosingPos3D(tMesh, -1, vx1->vtX, vx1->vtY, vx1->vtZ,
				     0, NULL) < 0) 
	{
	  errNum = WLZ_ERR_DOMAIN_DATA;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	      "%s: tie points line %d not in target mesh (%s).\n",
	      prog, vxCount, errMsg);
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	++vx0;
	++vx1;
	++vxCount;
      }
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    *dstNVx = vxCount;
    *dstVx0 = vxA0;
    *dstVx1 = vxA1;
  }
  else
  {
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
	"%s: Failed to read tie points file %s (%s).\n",
	prog, fStr, errMsg);
  }
  return(errNum);
}

/*!
* \return	void
* \brief	Outputs the mesh and displaced mesh of the given mesh
*		transform as Postscript to the standard error output.
* \param	mesh			Given mesh transform.
*/
static void	WlzMeshOutputPS(WlzMeshTransform *meshTr)
{
  int		elmCount,
  		nodCount;
  double	scale;
  WlzMeshNode	*nod;
  WlzMeshElem	*elm;
  WlzDVertex2	tVx,
  		offset;
  WlzDBox2	bBox;
  WlzDVertex2	elmVx[3];

  (void )fprintf(stderr, "%%!\n"
  			 "1 setlinecap\n"
			 "1 setlinejoin\n"
			 "0.01 setlinewidth\n");
  if(meshTr && (meshTr->type == WLZ_TRANSFORM_2D_MESH) &&
    (meshTr->nNodes > 0) && (meshTr->nElem > 0))
  {
    /* Compute the bounding box of the mesh transform. */
    nod = meshTr->nodes;
    bBox.xMin = bBox.xMax = nod->position.vtX;
    bBox.yMin = bBox.yMax = nod->position.vtY;
    nodCount = meshTr->nNodes;
    while(nodCount-- > 0)
    {
      tVx = nod->position;
      if(tVx.vtX < bBox.xMin)
      {
        bBox.xMin = tVx.vtX;
      }
      else if(tVx.vtX > bBox.xMax)
      {
        bBox.xMax = tVx.vtX;
      }
      if(tVx.vtY < bBox.yMin)
      {
        bBox.yMin = tVx.vtY;
      }
      else if(tVx.vtY > bBox.yMax)
      {
        bBox.yMax = tVx.vtY;
      }
      tVx.vtX += nod->displacement.vtX;
      tVx.vtY += nod->displacement.vtY;
      if(tVx.vtX < bBox.xMin)
      {
        bBox.xMin = tVx.vtX;
      }
      else if(tVx.vtX > bBox.xMax)
      {
        bBox.xMax = tVx.vtX;
      }
      if(tVx.vtY < bBox.yMin)
      {
        bBox.yMin = tVx.vtY;
      }
      else if(tVx.vtY > bBox.yMax)
      {
        bBox.yMax = tVx.vtY;
      }
      ++nod;
    }
    if(((tVx.vtX = fabs(bBox.xMax - bBox.xMin)) > 1.0) &&
       ((tVx.vtY = fabs(bBox.yMax - bBox.yMin)) > 1.0))
    {
      /* Compute scale and offset for an A4 page. */
      if((4.0 * tVx.vtX) > (3.0 * tVx.vtY))
      {
        scale = 500.0 / (bBox.xMax - bBox.xMin);
      }
      else
      {
        scale = 700.0 / (bBox.yMax - bBox.yMin);
      }
      offset.vtX = (bBox.xMin * scale) - 50;
      offset.vtY = (bBox.yMin * -scale) - 750;
      (void )fprintf(stderr, "255 0 0 setrgbcolor\n");
      elmCount = meshTr->nElem;
      elm = meshTr->elements;
      while(elmCount-- > 0)
      {
	nodCount = 3;
	while(--nodCount >= 0)
	{
	  nod = meshTr->nodes + elm->nodes[nodCount];
	  elmVx[nodCount] = nod->position;
	  elmVx[nodCount].vtX = (elmVx[nodCount].vtX * scale) - offset.vtX;
	  elmVx[nodCount].vtY = (elmVx[nodCount].vtY * -scale) - offset.vtY;
	}
	(void )fprintf(stderr,
		       "%g %g moveto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "stroke\n",
		       elmVx[0].vtX, elmVx[0].vtY,
		       elmVx[1].vtX, elmVx[1].vtY,
		       elmVx[2].vtX, elmVx[2].vtY,
		       elmVx[0].vtX, elmVx[0].vtY);
	++elm;
      }
      (void )fprintf(stderr, "0 0 255 setrgbcolor\n");
      elmCount = meshTr->nElem;
      elm = meshTr->elements;
      while(elmCount-- > 0)
      {
	nodCount = 3;
	while(--nodCount >= 0)
	{
	  nod = meshTr->nodes + elm->nodes[nodCount];
	  elmVx[nodCount].vtX = nod->position.vtX + nod->displacement.vtX;
	  elmVx[nodCount].vtY = nod->position.vtY + nod->displacement.vtY;
	  elmVx[nodCount].vtX = (elmVx[nodCount].vtX * scale) - offset.vtX;
	  elmVx[nodCount].vtY = (elmVx[nodCount].vtY * -scale) - offset.vtY;
	}
	(void )fprintf(stderr,
		       "%g %g moveto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "%g %g lineto "
		       "stroke\n",
		       elmVx[0].vtX, elmVx[0].vtY,
		       elmVx[1].vtX, elmVx[1].vtY,
		       elmVx[2].vtX, elmVx[2].vtY,
		       elmVx[0].vtX, elmVx[0].vtY);
	++elm;
      }
      (void )fprintf(stderr, "showpage\n");
    }
  }
}

/*!
* \return	void
* \brief	Outputs the mesh and displaced mesh of the given mesh
*		transform as Postscript to the standard error output.
* \param	mObj			Given mesh transform object.
* \param	sObj			Source object.
* \param	dilObj			Dilated object.
* \param	outObj			Transformed object.
*/
static void	WlzCMeshOutputPS(WlzObject *mObj,
				 WlzObject *sObj, WlzObject *dilObj,
				 WlzObject *outObj)
{
  int		idE,
  		idN,
		useElm;
  double	scale;
  double	*dsp;
  WlzDVertex2	pos,
  		offset;
  WlzDVertex2	elmVx[3];
  WlzDBox2	bBox;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D *elm;
  WlzCMesh2D	*mesh;
  WlzIndexedValues *ixv;
  WlzObject	*bObj = NULL;
  char		buf[100];

  (void )fprintf(stderr, "%%!\n"
			 "/Times-Courier findfont\n"
			 "1 scalefont\n"
			 "setfont\n"
  			 "1 setlinecap\n"
			 "1 setlinejoin\n"
			 "0.01 setlinewidth\n");
  if((mObj != NULL) && (mObj->type == WLZ_CMESH_2D) &&
     ((mesh = mObj->domain.cm2) != NULL) &&
     ((ixv = mObj->values.x) != NULL))
  {
    /* Compute the bounding box of the mesh transform. */
    bBox.xMin = bBox.yMin = DBL_MAX;
    bBox.xMax = bBox.yMax = DBL_MIN;
    for(idN = 0; idN < mesh->res.nod.maxEnt; ++idN)
    {
      nod = (WlzCMeshNod2D *)AlcVectorItemGet(mesh->res.nod.vec, idN);
      if(nod->idx >= 0)
      {
	if(nod->pos.vtX < bBox.xMin)
	{
	  bBox.xMin = nod->pos.vtX;
	}
	else if(nod->pos.vtX > bBox.xMax)
	{
	  bBox.xMax = nod->pos.vtX;
	}
	if(nod->pos.vtY < bBox.yMin)
	{
	  bBox.yMin = nod->pos.vtY;
	}
	else if(nod->pos.vtY > bBox.yMax)
	{
	  bBox.yMax = nod->pos.vtY;
	}
        dsp = (double *)WlzIndexedValueGet(ixv, idN);
	pos.vtX = nod->pos.vtX + dsp[0];
	pos.vtY = nod->pos.vtY + dsp[1];
	if(pos.vtX < bBox.xMin)
	{
	  bBox.xMin = pos.vtX;
	}
	else if(pos.vtX > bBox.xMax)
	{
	  bBox.xMax = pos.vtX;
	}
	if(pos.vtY < bBox.yMin)
	{
	  bBox.yMin = pos.vtY;
	}
	else if(pos.vtY > bBox.yMax)
	{
	  bBox.yMax = pos.vtY;
	}
      }
    }
    if(((pos.vtX = fabs(bBox.xMax - bBox.xMin)) > 1.0) &&
       ((pos.vtY = fabs(bBox.yMax - bBox.yMin)) > 1.0))
    {
      /* Compute scale and offset for an A4 page. */
      if((4.0 * pos.vtX) > (3.0 * pos.vtY))
      {
        scale = 500.0 / (bBox.xMax - bBox.xMin);
      }
      else
      {
        scale = 700.0 / (bBox.yMax - bBox.yMin);
      }
      offset.vtX = (bBox.xMin * scale) - 50;
      offset.vtY = (bBox.yMin * -scale) - 750;
      (void )fprintf(stderr, "255 0 0 setrgbcolor\n");
      for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
      {
        elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	if(elm->idx >= 0)
	{
	  for(idN = 0; idN < 3; ++idN)
	  {
	    pos = elm->edu[idN].nod->pos;
	    elmVx[idN].vtX = (pos.vtX * scale) - offset.vtX;
	    elmVx[idN].vtY = (pos.vtY * -scale) - offset.vtY;
	  }
	  (void )fprintf(stderr,
			 "%g %g moveto "
			 "%g %g lineto "
			 "%g %g lineto "
			 "%g %g lineto "
			 "stroke\n",
			 elmVx[0].vtX, elmVx[0].vtY,
			 elmVx[1].vtX, elmVx[1].vtY,
			 elmVx[2].vtX, elmVx[2].vtY,
			 elmVx[0].vtX, elmVx[0].vtY);
	  WLZ_VTX_2_ADD3(pos, elmVx[0], elmVx[1], elmVx[2]);
	  WLZ_VTX_2_SCALE(pos, pos, 0.333333);
	  (void )sprintf(buf, "%d", elm->idx);
	  (void )fprintf(stderr,
			 "%g %g moveto "
			 "(%s) show\n",
			 pos.vtX, pos.vtY,
			 buf);
	}
      }
      (void )fprintf(stderr, "0 0 255 setrgbcolor\n");
      for(idE = 0; idE < mesh->res.elm.maxEnt; ++idE)
      {
        elm = (WlzCMeshElm2D *)AlcVectorItemGet(mesh->res.elm.vec, idE);
	if(elm->idx >= 0)
	{
	  /* Display only those elements that have all nodes inside the
	   * domain.  
	  useElm = ((elm->edg[0].nod->flags | elm->edg[1].nod->flags |
	             elm->edg[2].nod->flags) &
		    WLZ_CMESH_NOD_FLAG_OUTSIDE) == 0;
	  */
	  useElm = 1;
	  if(useElm)
	  {
	    for(idN = 0; idN < 3; ++idN)
	    {
	      pos = elm->edu[idN].nod->pos;
	      dsp = (double *)WlzIndexedValueGet(ixv, elm->edu[idN].nod->idx);
	      pos.vtX += dsp[0];
	      pos.vtY += dsp[1];
	      elmVx[idN].vtX = (pos.vtX * scale) - offset.vtX;
	      elmVx[idN].vtY = (pos.vtY * -scale) - offset.vtY;
	    }
	    (void )fprintf(stderr,
			   "%g %g moveto "
			   "%g %g lineto "
			   "%g %g lineto "
			   "%g %g lineto "
			   "stroke\n",
			   elmVx[0].vtX, elmVx[0].vtY,
			   elmVx[1].vtX, elmVx[1].vtY,
			   elmVx[2].vtX, elmVx[2].vtY,
			   elmVx[0].vtX, elmVx[0].vtY);
	    WLZ_VTX_2_ADD3(pos, elmVx[0], elmVx[1], elmVx[2]);
	    WLZ_VTX_2_SCALE(pos, pos, 0.333333);
	    (void )sprintf(buf, "%d", elm->idx);
	    (void )fprintf(stderr,
	    		   "%g %g moveto "
			   "(%s) show\n",
			   pos.vtX, pos.vtY,
			   buf);
	  }
	}
      }
      (void )fprintf(stderr, "0 255 0 setrgbcolor\n");
      if((bObj = WlzObjToBoundary(sObj, 1, NULL)) != NULL)
      {
        WlzBndOutputPS(bObj->domain.b, scale, offset);
      }
      (void )fprintf(stderr, "0 255 255 setrgbcolor\n");
      if((bObj = WlzObjToBoundary(dilObj, 1, NULL)) != NULL)
      {
        WlzBndOutputPS(bObj->domain.b, scale, offset);
      }
      (void )fprintf(stderr, "127 0 127 setrgbcolor\n");
      if((bObj = WlzObjToBoundary(outObj, 1, NULL)) != NULL)
      {
        WlzBndOutputPS(bObj->domain.b, scale, offset);
      }
      (void )fprintf(stderr, "showpage\n");
    }
  }
}

/*!
* \return	void
* \brief	Outputs a boundary list as postscript.
* \param	bnd			Given boundary list.
* \param	scale			Scale.
* \param	offset			Offset.
*/
static void	WlzBndOutputPS(WlzBoundList *bnd,
			       double scale, WlzDVertex2 offset)
{
  WlzPolyOutputPS(bnd->poly, scale, offset);
  if(bnd->next)
  {
    WlzBndOutputPS(bnd->next, scale, offset);
  }
  if(bnd->down)
  {
    WlzBndOutputPS(bnd->down, scale, offset);
  }
}

/*!
* \return	void
* \brief	Outputs a poly line as postscript.
* \param	poly			Given poly line.
* \param	scale			Scale.
* \param	offset			Offset.
*/
static void	WlzPolyOutputPS(WlzPolygonDomain *poly,
			        double scale, WlzDVertex2 offset)
{
  int 		idN;

  if(poly && (poly->nvertices > 1))
  {
    for(idN = 1; idN < poly->nvertices; ++idN)
    {
      (void )fprintf(stderr,
		     "%g %g moveto "
                     "%g %g lineto stroke\n",
		     ((poly->vtx + idN - 1)->vtX * scale) - offset.vtX,
		     ((poly->vtx + idN - 1)->vtY * -scale) - offset.vtY,
		     ((poly->vtx + idN)->vtX * scale) - offset.vtX,
		     ((poly->vtx + idN)->vtY * -scale) - offset.vtY);
    }
    (void )fprintf(stderr,
		   "%g %g moveto "
                   "%g %g lineto stroke\n",
		   ((poly->vtx + poly->nvertices - 1)->vtX * scale) -
		   offset.vtX,
		   ((poly->vtx + poly->nvertices - 1)->vtY * scale) -
		   offset.vtY,
		   (poly->vtx->vtX * scale) - offset.vtX,
		   (poly->vtx->vtY * -scale) - offset.vtY);
  }
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
