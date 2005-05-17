#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzBasisFnTransformObj.c
* \author       Bill Hill
* \date         February 2005
* \version      $Id$
* \note
*               Copyright
*               2005 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Woolz filter which applies a basis function to a given
		object.
* \todo         -
* \bug          None known.
*/
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
				  WlzCMeshTransform *meshTr,
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

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		idx,
  		nTiePP,
		option,
		cdt = 0,
		cMesh = 0,
		vxCount = 0,
		vxLimit = 0,
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
  WlzDVertex2	*vx0,
  		*vx1,
		*vxVec0 = NULL,
		*vxVec1 = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL,
		*dilObj = NULL;
  WlzTransform  meshTr;
  WlzMeshGenMethod meshGenMth = WLZ_MESH_GENMETHOD_GRADIENT;
  WlzBasisFnTransform *basisTr = NULL;
  WlzFnType basisFnType = WLZ_FN_BASIS_2DMQ;
  WlzInterpolationType interp = WLZ_INTERPOLATION_NEAREST;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  char 		*rec,
  		*inObjFileStr,
		*basisFnTrFileStr = NULL,
		*tiePtFileStr = NULL,
  		*outObjFileStr;
  const char    *errMsg;
  static char	optList[] = "m:o:p:t:D:M:Y:cdghqsyBCGLNRTU",
  		inObjFileStrDef[] = "-",
		outObjFileStrDef[] = "-",
  		inRecord[IN_RECORD_MAX];

  opterr = 0;
  meshTr.core = NULL;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'B':
        meshGenMth = WLZ_MESH_GENMETHOD_BLOCK;
	break;
      case 'C':
        cMesh = 1;
	break;
      case 'D':
	if(1 != sscanf(optarg, "%d", &dbgFlg))
	{
	  usage = 1;
	  ok = 0;
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
	  ok = 0;
	}
	break;
      case 'M':
	if(sscanf(optarg, "%lg", &meshMaxDist) != 1)
	{
	  usage = 1;
	  ok = 0;
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
      case 's':
        basisFnType = WLZ_FN_BASIS_2DTPS;
	break;
      case 'y':
        basisFnType = WLZ_FN_BASIS_2DPOLY;
	break;
      case 't':
        basisFnTrFileStr = optarg;
        break;
      case 'R':
        restrictToBasisFn = 1;
	break;
      case 'T':
        outBasisTrFlag = 1;
	break;
      case 'U':
        outMeshTrFlag = 1;
	break;
      case 'Y':
        if(sscanf(optarg, "%d", &basisFnPolyOrder) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(cdt)
  {
    cMesh = 1;
    restrictToBasisFn = 1;
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (((tiePtFileStr == NULL) || (*tiePtFileStr == '\0')) &&
      ((basisFnTrFileStr == NULL) || (*basisFnTrFileStr == '\0'))) ||
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
		       "%s: failed to open tie points file %s (%s).\n",
		       *argv, tiePtFileStr, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
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
	      vxLimit = (vxLimit + 1024) * 2;
	      if(((vxVec0 = (WlzDVertex2 *)AlcRealloc(vxVec0,
				     vxLimit * sizeof(WlzDVertex2))) == NULL) ||
		 ((vxVec1 = (WlzDVertex2 *)AlcRealloc(vxVec1,
				     vxLimit * sizeof(WlzDVertex2))) == NULL))
	      {
		errNum = WLZ_ERR_MEM_ALLOC;
	      }
	      else
	      {
		vx0 = vxVec0 + vxCount;
		vx1 = vxVec1 + vxCount;
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
		++vx0;
		++vx1;
		++vxCount;
	      }
	    }
	  }
	}
	nTiePP = vxCount;
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
			 "%s: failed to read tie points file %s (%s).\n",
			 *argv, tiePtFileStr, errMsg);
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
		     "%s: failed to read object from file %s\n",
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
    if(restrictToBasisFn)
    {
      if(cMesh)
      {
	meshTr.cMesh = WlzCMeshTransformFromObj(inObj,
			      meshGenMth, meshMinDist, meshMaxDist,
			      &dilObj, &errNum);
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
		       "%s: failed to compute a mesh for the object (%s).\n",
		       *argv, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	basisTr = WlzBasisFnTrFromCPts2D(basisFnType, basisFnPolyOrder,
					  nTiePP, vxVec0, nTiePP, vxVec1,
					  (cdt)? dilObj: NULL, &errNum);
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
		       "%s: failed to compute basis function transform (%s).\n",
		       *argv, errMsg);
	}
      }
    }
    else /* Use affine and basis functions to define mesh. */
    {
      meshTr.mesh = WlzMeshTransformFromCPts(inObj,
      				basisFnType, basisFnPolyOrder,
      				nTiePP, vxVec0, nTiePP, vxVec1,
				meshGenMth, meshMinDist, meshMaxDist,
				&errNum);
    }
  }
  if(ok)
  {
    if((dbgFlg & 2) && basisTr->basisFn->distMap)
    {
      for(idx = 0; idx < basisTr->basisFn->nVtx; ++idx)
      {
        (void )WlzWriteObj(stderr, basisTr->basisFn->distMap[idx]);
      }
    }
  }
  if(ok)
  {
    if(outMeshTrFlag)
    {
      errNum = WLZ_ERR_UNIMPLEMENTED;
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: Writing mesh transforms has not been\n"
		     "implemented yet\n",
		     *argv);
    }
    else if(outBasisTrFlag)
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
	    errNum = WlzBasisFnSetCMesh(meshTr.cMesh, basisTr);
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
			   "%s: failed to set the mesh displacements (%s).\n",
			   *argv, errMsg);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(noTrObj == 0)
	{
	  if(cMesh)
	  {
	    outObj = WlzCMeshTransformObj(inObj, meshTr.cMesh, interp,
	    			          &errNum);
	  }
	  else
	  {
	    outObj = WlzMeshTransformObj(inObj, meshTr.mesh, interp,
	    				 &errNum);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
	if(dbgFlg & 1)
	{
	  if(cMesh)
	  {
	    WlzCMeshOutputPS(meshTr.cMesh, inObj, dilObj, outObj);
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
		       "%s: failed to transform object (%s).\n",
		       *argv, errMsg);
      }
      if(errNum == WLZ_ERR_NONE)
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
			 "%s: failed to write output object (%s).\n",
			 *argv, errMsg);
	}
	if(fP && strcmp(outObjFileStr, "-"))
	{
	  fclose(fP);
	}
      }
    }
  }
  if(cMesh)
  {
    (void )WlzFreeCMeshTransform(meshTr.cMesh);
  }
  else
  {
    (void )WlzMeshFreeTransform(meshTr.mesh);
  }
  (void )WlzBasisFnFreeTransform(basisTr);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  (void )WlzFreeObj(dilObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-p<tie points file>]\n"
    "                  [-m<min mesh dist>] [-M<max mesh dist>]\n"
    "                  [-t<basis fn transform>] [-Y<order of polynomial>]\n"
    "                  [-D<flags>]\n"
    "                  [-d] [-g] [-h] [-q] [-s] [-y]\n"
    "                  [-B] [-C] [-G] [-L] [-N] [-T]\n"
    "                  [<in object>]\n"
    "Options:\n"
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
    "  -t  Basis function transform object.\n"
    "  -c  Use conformal polynomial basis function if tie points are given.\n"
    "  -d  Compute transform using conforming distance transforms, also\n"
    "      sets the mesh generation to produce a conforming mesh and\n"
    "      restricts the transformation to a basis function transform\n"
    "      only (ie no least squares affine).\n"
    "  -g  Use Gaussian basis function if tie points are given.\n"
    "  -h  Help, prints this usage message.\n"
    "  -q  Use multi-quadric basis function if tie points are given.\n"
    "  -s  Use thin plate spline basis function (default) if tie points\n"
    "      are given.\n"
    "  -T  Output a basis function transform instead of a transformed\n"
    "      object.\n"
    "  -U  Output a mesh transform instead of a transformed object.\n"
    "  -y  Use polynomianl basis function if tie points are given.\n"
    "  -Y  Polynomial order for polynomial basis function (default 3).\n"
    "Computes and applies Woolz basis function transforms.\n"
    "Tie points may be read from an ascii file with the format:\n"
    "  <vertex x> <vertex y> <displacement x> <displacement y>\n"
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
* \param	meshTr			Given mesh transform.
* \param	sObj			Source object.
* \param	dilObj			Dilated object.
* \param	outObj			Transformed object.
*/
static void	WlzCMeshOutputPS(WlzCMeshTransform *meshTr,
				 WlzObject *sObj, WlzObject *dilObj,
				 WlzObject *outObj)
{
  int		idE,
  		idN,
		useElm;
  double	scale;
  WlzDVertex2	pos,
  		offset;
  WlzDVertex2	*dsp;
  WlzDVertex2	elmVx[3];
  WlzDBox2	bBox;
  WlzCMeshNod2D	*nod;
  WlzCMeshElm2D *elm;
  WlzCMesh2D	*mesh;
  WlzObject	*bObj = NULL;
  char		buf[100];

  (void )fprintf(stderr, "%%!\n"
			 "/Times-Courier findfont\n"
			 "2 scalefont\n"
			 "setfont\n"
  			 "1 setlinecap\n"
			 "1 setlinejoin\n"
			 "0.01 setlinewidth\n");
  if((meshTr != NULL) &&
     (meshTr->type == WLZ_TRANSFORM_2D_CMESH) &&
     ((mesh = meshTr->mesh.m2) != NULL) &&
     (meshTr->dspVec != NULL))
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
        dsp = (WlzDVertex2 *)AlcVectorItemGet(meshTr->dspVec, idN);
	WLZ_VTX_2_ADD(pos, nod->pos, *dsp);
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
	    pos = elm->edg[idN].nod->pos;
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
	      pos = elm->edg[idN].nod->pos;
	      dsp = (WlzDVertex2 *)AlcVectorItemGet(meshTr->dspVec,
						    elm->edg[idN].nod->idx);
	      WLZ_VTX_2_ADD(pos, pos, *dsp);
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
