#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPerturbLandmarks_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzPerturbLandmarks.c
* \author       Bill Hill
* \date         October 2015
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2015],
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
* \brief	Perturbs a set of landmark vertices using a random
* 		displacement, while constraining the landmarks to
* 		lie within a given domain.
* \ingroup	BinWlz
*/

/*!
\ingroup BinWlz
\defgroup wlzperturblandmarks WlzPerturbLandmarks
\par Name
WlzPerturbLandmarks - perturbs a of landmark vertices.
\par Synopsis
\verbatim
WlzPerturbLandmarks [-A] [-2] [-3] [-h] [-c <constraining domain>]
                    [-d <max displacement>] [-o <output landmark file>]
		    [-s <seed>] [<input landmark file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-A</b></td>
    <td>Landmarks are absolute (rather than relative). This flag is
        respected by both the input and output.</td>
  </tr>
  <tr>
    <td><b>-2</b></td>
    <td>Landmarks and optional constraining domain are 2D.</td>
  </tr>
  <tr>
    <td><b>-3</b></td>
    <td>Landmarks and optional constraining domain are 3D (default).</td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Constraining domain. If not supplied there is no constraint.
        The constraint domain may be either a mesh or spatial domain
	object. If a mesh is given then the landmarks will be constrained
	to lie within it and will remain valid for constrained mesh
	transforms.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Maximum displacement distance (default 1.0, one pixel/voxel).</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Seed for random displacements.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints a usage message.</td>
  </tr>
</table>
Perturbs a set of landmark vertices using a random displacement,
while constraining the landmarks to lie within a given domain.
Landmarks are read from a text file with the format:
\verbatim
  <vertex x> <vertex y> [<vertex z>] <disp x> <disp y> [<disp z>]
\endverbatim
unless the absolute flag is set in which case the displacement columns
are interpreted as target vertices.
By default files are read from the standard input and the perturbed landmarks
file is written to the standard output.
\par Examples
\verbatim
WlzPerturbLandmarks -o out.num -c mesh.wlz -d 10 landmarks.num
\endverbatim
3D landmarks are read from the white space seperated text file 'landmarks.num'
and perturbed by addding a random displacement which is within the maximum
distance 10 and the constraining domain (read from 'mesh.wlz').
The perturbed landmarks are written to the file 'out.num'.
\par File
\ref WlzPerturbLandmarks.c "WlzPerturbLandmarks.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzAffineTransformLSq  "WlzAffineTransformLSq(1)"
\ref WlzBasisFnTransformObj "WlzBasisFnTransformObj(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

#ifdef _OPENMP
#include <omp.h>
#endif

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern char     *optarg;
extern int      optind,
		opterr,
		optopt;


static double	WlzRRandUniform(unsigned int *seed)
{
  double	r;

#if defined(_OPENMP) && defined(HAVE_RAND_R)
  r = ((double )rand_r(seed)) / RAND_MAX;
#else
  r = AlgRandUniform();
#endif
  return(r);
}

/*!
* \return	Position in mesh.
* \brief	Ensures given vertex is within mesh object.
* \param	mshObj			Given mesh object, must be valid.
* \param	dim			Dimension (2 or 3).
* \param	gvnPos			Given position, must be of appropriate
* 					type: WlzDVertex2 for dim == 2
* 					or WlzDVertex3 for dim == 3.
* \param        maxDist			Maximum distance within which to
* 					search for position in the mesh,
* 					should always be greater than zero.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDVertex3 WlzEnsureVertexInMesh(WlzObject *mshObj, int dim,
                            	WlzVertex gvnPos, double maxDist,
				WlzErrorNum *dstErr)
{
  WlzDVertex3 pos;

  if(dim == 2)
  {
    WlzDVertex2 pos2;

    WLZ_VTX_2_SET(pos2, gvnPos.d2.vtX, gvnPos.d2.vtY);
    (void )WlzCMeshElmClosestPosIn2D(mshObj->domain.cm2, &pos2, pos2, maxDist);
    WLZ_VTX_3_SET(pos, pos2.vtX, pos2.vtY, 0.0);
  }
  else
  {
    pos = gvnPos.d3;
    (void )WlzCMeshElmClosestPosIn3D(mshObj->domain.cm3, &pos, pos, maxDist);
  }
  return(pos);
}

/*!
* \return	Random nearby vertex position.
* \brief	Computes a random vertex position which is nearby to the given
* 		position and within given spatial domain object.
* 		This isn't very efficient code, if this ever migrates into
* 		a library it should probably be rewritten.
* \param	cnsObj			Optional constraining spatial domain
* 					object.
* \param	mshObj			Optional given mesh object.
* \param	dim			Dimension (2 or 3).
* \param	gvnPos			Given position, must be of appropriate
* 					type: WlzDVertex2 for dim == 2
* 					or WlzDVertex3 for dim == 3.
* \param	maxDsp			Maximum displacement.
* \param	seedP			Seed pointer for rand_r();
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzDVertex3 WlzRandNearby(WlzObject *cnsObj, WlzObject *mshObj,
			         int dim, WlzVertex gvnPos, double maxDsp,
				 unsigned int *seedP, WlzErrorNum *dstErr)
{
  WlzDVertex3	pos;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const double	eps = 1.0e-04;

  if(dim == 2)
  {
    pos.vtX = gvnPos.d2.vtX;
    pos.vtY = gvnPos.d2.vtY;
    pos.vtZ = 0.0;
  }
  else /* dim == 3 */
  {
    pos = gvnPos.d3;
  }
  if(maxDsp > eps)
  {
    if(cnsObj == NULL)
    {
      double	md2;

      md2 = maxDsp * 2.0;
      pos.vtX += (WlzRRandUniform(seedP) - 0.5) * md2;
      pos.vtY += (WlzRRandUniform(seedP) - 0.5) * md2;
      if(dim == 3)
      {
	pos.vtZ += (WlzRRandUniform(seedP) - 0.5) * md2;
      }
    }
    else
    {
      WlzVertex  posV;
      WlzVertexP posP;
      WlzObject  *dObj = NULL;

      if(dim == 2)
      {
        WLZ_VTX_2_SET(posV.d2, pos.vtX, pos.vtY);
	posP.d2 = &(posV.d2);
      }
      else
      {
        posV.d3 = pos;
	posP.d3 = &(posV.d3);
      }
      dObj = WlzAssignObject(
	     WlzDomainNearby(cnsObj, 1, posP, WLZ_OCTAGONAL_DISTANCE,
	                     maxDsp, &errNum), NULL);
      if(errNum == WLZ_ERR_NONE)
      {
	if((dObj == NULL) || (dObj->type == WLZ_EMPTY_OBJ))
	{
	  pos = WlzRandNearby(NULL, mshObj, dim, gvnPos, maxDsp, seedP, NULL);
	  /* Reject any perturbed points with distance greater than the
	   * maximum, which can happen if the original point is outside
	   * the mesh. */
	  if(dim == 2)
	  {
	    WlzDVertex2 t;

	    WLZ_VTX_2_SUB(t, pos, gvnPos.d2);
	    if(WLZ_VTX_2_LENGTH(t) > maxDsp)
	    {
	      WLZ_VTX_3_SET(pos, 0.0, gvnPos.d2.vtY, gvnPos.d2.vtX);
	    }
	  }
	  else
	  {
	    WlzDVertex3 t;

	    WLZ_VTX_3_SUB(t, pos, gvnPos.d3);
	    if(WLZ_VTX_3_LENGTH(t) > maxDsp)
	    {
	      pos = gvnPos.d3;
	    }
	  }
	}
	else
	{
	  WlzPixelV  thrV;
	  WlzLong    vol = 0;
	  WlzObject  *tObj = NULL,
		     *vObj = NULL;
	  WlzIBox3   box = {0};

	  vol = WlzVolume(dObj, &errNum);
	  if(vol < 1)
	  {
	    errNum = WLZ_ERR_DOMAIN_DATA;
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    thrV.type = WLZ_GREY_INT;
	    thrV.v.inv = (int )floor((vol - 1) * WlzRRandUniform(seedP));
	    vObj = WlzAssignObject(
		   WlzGreyNewIncValues(dObj, &errNum), NULL);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tObj = WlzAssignObject(
		   WlzThreshold(vObj, thrV, WLZ_THRESH_EQUAL, &errNum), NULL);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    box = WlzBoundingBox3I(tObj, &errNum);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    pos.vtX = box.xMin;
	    pos.vtY = box.yMin;
	    pos.vtZ = box.zMin;
	  }
	  (void )WlzFreeObj(tObj);
	  (void )WlzFreeObj(vObj);
	}
      }
      (void )WlzFreeObj(dObj);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pos);
}

int             main(int argc, char **argv)
{
  int		option,
		rel = 1,
		dim = 3,
		nLmk = 0,
		seed = 0,
		ok = 1,
		usage = 0;
  double	maxDsp = 1.0;
  char		*inFile,
		*outFile,
		*cnsFile = NULL;
  WlzVertexP	srcV,
                dstV;
  WlzObject	*cnsObj = NULL,
		*mshObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  const double	eps = 1.0e-04;
  static char	optList[] = "23Ahc:d:o:s:",
  		fileStrDef[] = "-";

  opterr = 0;
  srcV.v = dstV.v = NULL;
  inFile = outFile = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case '2':
        dim = 2;
	break;
      case '3':
        dim = 3;
	break;
      case 'A':
        rel = 0;
	break;
      case 'c':
        cnsFile = optarg;
	break;
      case 'd':
        if((sscanf(optarg, "%lg", &maxDsp) != 1) || (maxDsp < 0.0))
	{
	  usage = 1;
	}
	break;
      case 'o':
        outFile = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%u", &seed) != 1)
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
      inFile = argv[optind];
    }
  }
  ok = !usage;
  /* Read the landmarks. */
  if(ok)
  {
    size_t	nR = 0,
    		nC = 0;
    double	**vtx = NULL;
    FILE	*fP = NULL;

    if((fP = (strcmp(inFile, "-")? fopen(inFile, "r"): stdin)) == NULL)
    {
      errNum = WLZ_ERR_READ_EOF;
    }
    else
    {
      if(AlcDouble2ReadAsci(fP, &vtx, &nR, &nC) != ALC_ER_NONE)
      {
	errNum = WLZ_ERR_READ_INCOMPLETE;
      }
      if(strcmp(inFile, "-"))
      {
        (void )fclose(fP);
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      if((nR < 1) || (dim * 2 != nC))
      {
        errNum = WLZ_ERR_PARAM_DATA;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
     size_t 	sz;

     nLmk = nR;
     sz = (dim == 2)? sizeof(WlzDVertex2): sizeof(WlzDVertex3);
     if(((srcV.v = AlcMalloc(sz * nLmk)) == NULL) ||
        ((dstV.v = AlcMalloc(sz * nLmk)) == NULL))
      {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
    if(errNum == WLZ_ERR_NONE)
    {
      int 	idx;

      for(idx = 0; idx < nLmk; ++idx)
      {
	double *v;

	v = vtx[idx];
	if(dim == 2)
	{
	  WLZ_VTX_2_SET(srcV.d2[idx], v[0], v[1]);
	  WLZ_VTX_2_SET(dstV.d2[idx], v[2], v[3]);
	  if(rel)
	  {
	    WLZ_VTX_2_ADD(dstV.d2[idx], dstV.d2[idx], srcV.d2[idx]);
	  }
	}
	else
	{
	  WLZ_VTX_3_SET(srcV.d3[idx], v[0], v[1], v[2]);
	  WLZ_VTX_3_SET(dstV.d3[idx], v[3], v[4], v[5]);
	  if(rel)
	  {
	    WLZ_VTX_3_ADD(dstV.d3[idx], dstV.d3[idx], srcV.d3[idx]);
	  }
	}
      }
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to read landmarks from the file %s (%s)\n",
	  *argv, inFile, errMsg);
    }
    (void )AlcDouble2Free(vtx);
  }
  /* Read the constraining domain if required. */
  if(ok && cnsFile)
  {
    FILE *fP = NULL;
    
    errNum = WLZ_ERR_FILE_OPEN;
    if((cnsFile == NULL) ||
        (*cnsFile == '\0') ||
        ((fP = (strcmp(cnsFile, "-")?
               fopen(cnsFile, "r"): stdin)) == NULL) ||
        ((cnsObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
        (errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to read constraining object from file %s (%s).\n",
	  *argv, cnsFile, errMsg);
    }
    else
    {
      switch(cnsObj->type)
      {
	case WLZ_2D_DOMAINOBJ: /* FALLTHROUGH */
	case WLZ_3D_DOMAINOBJ:
	  break;
	case WLZ_CMESH_2D: /* FALLTHROUGH */
	case WLZ_CMESH_3D:
	  {
	    WlzObject *tObj = NULL;

	    mshObj = cnsObj;
	    tObj = WlzAssignObject(
	           WlzCMeshToDomObj(mshObj, 0, 1.0, &errNum), NULL);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      cnsObj = WlzAssignObject(
	               WlzDilation(tObj, WLZ_26_CONNECTED, &errNum), NULL);
	    }
	    (void )WlzFreeObj(tObj);
	  }
	  break;
        default:
          errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        switch(cnsObj->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	   if(dim != 2)
	   {
	     errNum = WLZ_ERR_OBJECT_TYPE;
	   }
	   break;
	  case WLZ_3D_DOMAINOBJ:
	   if(dim != 3)
	   {
	     errNum = WLZ_ERR_OBJECT_TYPE;
	   }
	   break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	    "%s: Constraining object type is inappropriate (%s).\n",
	    *argv,  errMsg);
      }
    }
    if(fP && (strcmp(inFile, "-")))
    {
      (void )fclose(fP);
    }
  }
  /* Ensure landmarks are within the mesh (some may be very close, but
   * just outside). */
  if(ok && mshObj)
  {
    int		idx;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idx = 0; idx < nLmk; ++idx)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzVertex   posV;
	WlzDVertex3 pos3;
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

	if(dim == 2)
	{
	  posV.d2 = dstV.d2[idx];
	}
	else /* dim == 3 */
	{
	  posV.d3 = dstV.d3[idx];
	}
	pos3 = WlzEnsureVertexInMesh(mshObj, dim, posV, 10.0, &errNum2);
	if(errNum2 == WLZ_ERR_NONE)
	{
	  if(dim == 2)
	  {
	    WLZ_VTX_2_SET(dstV.d2[idx], pos3.vtX, pos3.vtY);
	  }
	  else /* dim == 3 */
	  {
	    dstV.d3[idx] = pos3;
	  }
	}
	else
	{
#ifdef _OPENMP
#pragma omp critical
	  {
#endif
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
#ifdef _OPENMP
	  }
#endif
	}
      }
    }
  }
  /* Perturb landmarks. */
  if(ok && (maxDsp > eps))  
  {
    int		idx,
    		nThr = 1;
    unsigned int *seeds = NULL;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        nThr = omp_get_num_threads();
      }
    }
#endif
    if((seeds = (unsigned int *)AlcCalloc(nThr, sizeof(unsigned int))) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      seeds[0] = seed;
      for(idx = 0; idx < nThr; ++idx)
      {
#if defined(_OPENMP) && defined(HAVE_RAND_R)
        seeds[idx] = rand_r(&(seeds[0]));
#else
        seeds[idx] = AlgRandUniform();
#endif
      }
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(idx = 0; idx < nLmk; ++idx)
      {
	if(errNum == WLZ_ERR_NONE)
	{
	  int	      thrId = 0;
	  WlzVertex   posV;
	  WlzErrorNum errNum2 = WLZ_ERR_NONE;

#ifdef _OPENMP
          thrId = omp_get_thread_num();
#endif
	  if(dim == 2)
	  {
	    posV.d2 = dstV.d2[idx];
	  }
	  else /* dim == 3 */
	  {
	    posV.d3 = dstV.d3[idx];
	  }
	  if(dim == 2)
	  {
	    WlzDVertex3 tmpV;

	    tmpV = WlzRandNearby(cnsObj, mshObj, dim, posV, maxDsp,
	                         &(seeds[thrId]), &errNum2);
	    if(errNum2 == WLZ_ERR_NONE)
	    {
	      WLZ_VTX_2_SET(dstV.d2[idx], tmpV.vtX, tmpV.vtY);
	    }
	  }
	  else
	  {
	    dstV.d3[idx] = WlzRandNearby(cnsObj, mshObj, dim, posV, maxDsp,
					 &(seeds[thrId]), &errNum2);
	  }
	  if(errNum2 != WLZ_ERR_NONE)
	  {
#ifdef _OPENMP
#pragma omp critical
	    {
#endif
	      if(errNum == WLZ_ERR_NONE)
	      {
		errNum = errNum2;
	      }
#ifdef _OPENMP
	    }
#endif
	  }
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
	  "%s: Failed to compute perturbed positions (%s).\n",
	  *argv,  errMsg);
    }
  }
  /* Again ensure landmarks are within the mesh. */
  if(ok && mshObj)
  {
    int		idx;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(idx = 0; idx < nLmk; ++idx)
    {
      if(errNum == WLZ_ERR_NONE)
      {
	WlzVertex   posV;
	WlzDVertex3 pos3;
	WlzErrorNum errNum2 = WLZ_ERR_NONE;

	if(dim == 2)
	{
	  posV.d2 = dstV.d2[idx];
	}
	else /* dim == 3 */
	{
	  posV.d3 = dstV.d3[idx];
	}
	pos3 = WlzEnsureVertexInMesh(mshObj, dim, posV, 10.0, &errNum2);
	if(errNum2 == WLZ_ERR_NONE)
	{
	  if(dim == 2)
	  {
	    WLZ_VTX_2_SET(dstV.d2[idx], pos3.vtX, pos3.vtY);
	  }
	  else /* dim == 3 */
	  {
	    dstV.d3[idx] = pos3;
	  }
	}
	else
	{
#ifdef _OPENMP
#pragma omp critical
	  {
#endif
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = errNum2;
	    }
#ifdef _OPENMP
	  }
#endif
	}
      }
    }
  }
  /* Write landmarks. */
  if(ok)  
  {
    int 	idx;
    FILE	*fP = NULL;

    if(rel)
    {
      if(dim == 2)
      {
	for(idx = 0; idx < nLmk; ++idx)
	{
	  WLZ_VTX_2_SUB(dstV.d2[idx], dstV.d2[idx], srcV.d2[idx]);
	}
      }
      else
      {
	for(idx = 0; idx < nLmk; ++idx)
	{
	  WLZ_VTX_3_SUB(dstV.d3[idx], dstV.d3[idx], srcV.d3[idx]);
	}
      }
    }
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFile, "-")? fopen(outFile, "w"): stdout)) != NULL)
    {
      if(dim == 2)
      {
	for(idx = 0; idx < nLmk; ++idx)
	{
	  if(fprintf(fP, "%lg %lg %lg %lg\n",
	             srcV.d2[idx].vtX, srcV.d2[idx].vtY,
	             dstV.d2[idx].vtX, dstV.d2[idx].vtY) < 3)
	  {
	    break;
	  }
	}
      }
      else
      {
	for(idx = 0; idx < nLmk; ++idx)
	{
	  if(fprintf(fP, "%lg %lg %lg %lg %lg %lg\n",
	             srcV.d3[idx].vtX, srcV.d3[idx].vtY, srcV.d3[idx].vtZ,
	             dstV.d3[idx].vtX, dstV.d3[idx].vtY, dstV.d3[idx].vtZ) < 6)
	  {
	    break;
	  }
	}
      }
      if(idx == nLmk)
      {
        errNum = WLZ_ERR_NONE;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to write landmarks to file %s (%s).\n",
	  *argv, outFile, errMsg);
    }
    if(fP && (strcmp(outFile, "-")))
    {
      (void )fclose(fP);
    }
  }
  AlcFree(srcV.v);
  AlcFree(dstV.v);
  WlzFreeObj(cnsObj);
  WlzFreeObj(mshObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-A] [-2] [-3] [-h] [-c <constraining domain>]\n"
    "\t\t[-d <max displacement>] [-o <output landmark file>]\n"
    "\t\t[-s <seed>] [<input landmark file>]\n"
    "Perturbs a set of landmark vertices using a random displacement,\n"
    "while constraining the landmarks to lie within a given domain.\n"
    "Landmarks are read from a text file with the format:\n"
    "  <vertex x> <vertex y> [<vertex z>] <disp x> <disp y> [<disp z>]\n"
    "unless the absolute flag is set in which case the displacement columns\n"
    "are interpreted as target vertices.\n"
    "By default files are read from the standard input and the perturbed\n"
    "landmarks file is written to the standard output.\n"
    "Version: %s\n"
    "Options are:\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file\n"
    "  -A  Landmarks are absolute (rather than relative). This flag is\n"
    "      respected by both the input and output.\n"
    "  -2  Landmarks and optional constraining domain are 2D.\n"
    "  -3  Landmarks and optional constraining domain are 3D (default).\n"
    "  -c  Constraining domain. If not supplied there is no constraint.\n"
    "      The constraint domain may be either a mesh or spatial domain\n"
    "      object. If a mesh is given then the landmarks will be constrained\n"
    "      to lie within it and will remain valid for constrained mesh\n"
    "      transforms.\n"
    "  -d  Maximum displacement distance (default 1.0, one pixel/voxel).\n"
    "  -s  Seed for random displacements.\n"
    "Example:\n"
    "  %s -o out.num -c mesh.wlz -d 10 landmarks.num\n"
    "3D landmarks are read from the white space seperated text file\n"
    "'landmarks.num' and perturbed by addding a random displacement\n"
    "to the target vertex which is within the maximum distance 10 and\n"
    "the constraining domain (read from 'mesh.wlz').  The perturbed\n"
    "landmarks are written to the file 'out.num'.\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  exit(ok != 0);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
